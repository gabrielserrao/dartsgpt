import numpy as np
from darts.reservoirs.struct_reservoir import StructReservoir
from darts.reservoirs.unstruct_reservoir import UnstructReservoir
from darts.models.darts_model import DartsModel

from darts.physics.super.physics import Compositional
from darts.physics.super.property_container import PropertyContainer

from darts.physics.properties.basic import PhaseRelPerm, ConstFunc
from darts.physics.properties.density import Garcia2001
from darts.physics.properties.viscosity import Fenghour1998, Islam2012
from darts.physics.properties.eos_properties import EoSDensity, EoSEnthalpy

from dartsflash.libflash import NegativeFlash
from dartsflash.libflash import CubicEoS, AQEoS, FlashParams, InitialGuess
from dartsflash.components import CompData
from darts.physics.properties.flash import ConstantK


class Model(DartsModel):
    def set_reservoir(self, nx=100):
        L = 1000
        H = 50
        ny = 1
        nz = 1
        nb = nx * ny * nz

        #x_axes = np.logspace(0, np.log10(L), nx+1)
        x_axes = np.linspace(0, L, nx+1)        
        dx_slice = x_axes[1:] - x_axes[:-1]
        self.dx = np.tile(dx_slice, nz)

        self.dz = H / nz
        depth = np.zeros(nb)
        n_layer = nx*ny
        for k in range(nz):
            depth[k*n_layer:(k+1)*n_layer] = 1000 + k * self.dz

        self.reservoir = StructReservoir(self.timer, nx, ny, nz, dx=self.dx, dy=10, dz=self.dz,
                                         permx=100, permy=100, permz=10, hcap=2200, rcond=100, poro=0.2, depth=depth)
        self.reservoir.boundary_volumes['yz_plus'] = 1e12

    def set_reservoir_core(self, nx=100, rad_grid=False):
        L = 0.1  # length of reservoir
        ny = 1
        nz = 1
        self.nb = nx * ny * nz

        x_axes = np.logspace(-4, np.log10(L), nx+1)
        dx_slice = x_axes[1:] - x_axes[:-1]
        self.dx = np.tile(dx_slice, nz)

        dy_slice = np.zeros(nx)
        if rad_grid:
            for i in range(nx):
                dy_slice[i] = 2 * np.pi * (x_axes[i+1])
        else:
            dy_slice[:] = 0.04

        dy = np.tile(dy_slice, nz)

        self.dz = 0.04
        depth = np.zeros(self.nb)
        n_layer = nx*ny
        for k in range(nz):
            depth[k*n_layer:(k+1)*n_layer] = 1.0 + k * self.dz

        self.reservoir = StructReservoir(self.timer, nx, ny, nz, dx=self.dx, dy=dy, dz=self.dz,
                                         permx=100, permy=100, permz=10, hcap=2200, rcond=100, poro=0.2, depth=depth)
        #self.reservoir.boundary_volumes['yz_plus'] = 1e8

    def set_reservoir_unstr(self):
        const_perm = 1000
        permx = const_perm
        permy = const_perm
        permz = const_perm
        poro = 0.3
        frac_aper = 0e-4
        mesh_file = 'original_wellnew.msh'
        self.reservoir = UnstructReservoir(self.timer, permx=permx, permy=permy, permz=permz, frac_aper=frac_aper,
                                           mesh_file=mesh_file, poro=poro)

    def set_wells(self):
        self.reservoir.add_well("I1")
        self.reservoir.add_perforation("I1", cell_index=(1, 1, 1), well_index=1e5, well_indexD=0)

        self.reservoir.add_well("P1")
        self.reservoir.add_perforation("P1", cell_index=(self.reservoir.nx, 1, 1), well_index=1e5, well_indexD=0)

    def set_wells_2D(self):
        self.reservoir.add_well("I1")
        for k in range(int(self.reservoir.nz/2), self.reservoir.nz):
            self.reservoir.add_perforation("I1", cell_index=(1, 1, k), well_radius=0.1, well_indexD=0)

        # self.reservoir.add_well("P1")
        # for k in range(self.reservoir.nz):
        #     self.reservoir.add_perforation("P1", cell_index=(self.reservoir.nx, self.reservoir.ny, k+1),
        #                                    well_index=100, well_indexD=0)

    def set_wells_unstr(self):
        self.reservoir.add_well("P1")
        self.reservoir.add_perforation("P1", cell_index=(1, 1, 1), well_index=1e5, well_indexD=0)

        # self.reservoir.add_well("P1")
        # for k in range(self.reservoir.nz):
        #     self.reservoir.add_perforation("P1", cell_index=(self.reservoir.nx, self.reservoir.ny, k+1),
        #                                    well_index=100, well_indexD=0)

    def set_physics(self,  zero, n_points, components, temperature=None, temp_inj=350.):
        """Physical properties"""
        # Fluid components, ions and solid
        phases = ["Aq", "V"]
        comp_data = CompData(components, setprops=True)

        pr = CubicEoS(comp_data, CubicEoS.PR)
        # aq = Jager2003(comp_data)
        #aq = AQEoS(comp_data, AQEoS.Ziabakhsh2012)
        aq_evaluators = {AQEoS.water: AQEoS.Jager2003,
                         AQEoS.solute: AQEoS.Ziabakhsh2012}
        aq = AQEoS(comp_data, aq_evaluators)

        flash_params = FlashParams(comp_data)

        # EoS-related parameters
        flash_params.add_eos("PR", pr)
        flash_params.add_eos("AQ", aq)
        eos_used = ["AQ", "PR"]
        split_initial_guesses = [InitialGuess.Henry_AV]

        # Flash-related parameters
        # flash_params.split_switch_tol = 1e-3

        if temperature is None:  # if None, then thermal=True
            thermal = True
        else:
            thermal = False

        """ properties correlations """
        property_container = PropertyContainer(phases_name=phases, components_name=components, Mw=comp_data.Mw,
                                               temperature=temperature, rock_comp=0, min_z=zero/10)
        diff_coef = -1e-1

        property_container.flash_ev = NegativeFlash(flash_params, eos_used, split_initial_guesses)
        # K_val = np.array([110, 0.016, 0.0015])
        # property_container.flash_ev = ConstantK(len(components), K_val, zero)
        property_container.density_ev = dict([('V', EoSDensity(pr, comp_data.Mw)),
                                              ('Aq', Garcia2001(components))])
        property_container.viscosity_ev = dict([('V', Fenghour1998()),
                                                ('Aq', Islam2012(components))])
        property_container.rel_perm_ev = dict([('V', PhaseRelPerm("gas")),
                                               ('Aq', PhaseRelPerm("oil", swc=0.4))])
        property_container.diffusion_ev = dict([('V', ConstFunc(np.ones(len(components)) * diff_coef)),
                                                ('Aq', ConstFunc(np.ones(len(components)) * diff_coef))])

        property_container.enthalpy_ev = dict([('V', EoSEnthalpy(pr)),
                                               ('Aq', EoSEnthalpy(aq))])
        property_container.conductivity_ev = dict([('V', ConstFunc(0.)),
                                                   ('Aq', ConstFunc(172)), ])

        property_container.output_props = {'satA': lambda: property_container.sat[0],
                                           'satV': lambda: property_container.sat[1],
                                           'densA': lambda: property_container.dens[0],
                                           'densV': lambda: property_container.dens[1],
                                           }

        state_spec = Compositional.StateSpecification.PT
        self.physics = Compositional(components, phases, self.timer, n_points, min_p=1, max_p=400, min_z=zero/10,
                                     max_z=1-zero/10, min_t=223.15, max_t=473.15, state_spec=state_spec, cache=False)
        self.physics.add_property_region(property_container)

    def set_initial_conditions(self):
        self.physics.set_initial_conditions_from_array(self.reservoir.mesh, self.initial_values)
        return

    def set_well_controls(self):
        from darts.engines import well_control_iface
        for i, w in enumerate(self.reservoir.wells):
            if 'I' in w.name:
                self.physics.set_well_controls(wctrl=w.control, is_inj=True, control_type=well_control_iface.BHP,
                                               target=self.p_inj, inj_composition=self.inj_comp, inj_temp=self.inj_temp)
            else:
                self.physics.set_well_controls(wctrl=w.control, is_inj=False, control_type=well_control_iface.BHP,
                                               target=self.p_prod)


