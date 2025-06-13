from darts.models.darts_model import DartsModel
from darts.engines import value_vector, index_vector, operator_set_evaluator_iface
import numpy as np

from darts.physics.super.physics import Compositional
from darts.physics.super.property_container import PropertyContainer

from darts.physics.properties.basic import ConstFunc, PhaseRelPerm
from darts.physics.properties.density import DensityBasic
from darts.physics.properties.enthalpy import EnthalpyBasic

from reservoir import UnstructReservoir

from darts.models.opt.opt_module_settings import OptModuleSettings

class Model(DartsModel, OptModuleSettings):
    def __init__(self, discr_type, mesh_file, T=2000, report_step=100.0, customize_new_operator=False):
        # call base class constructor
        super().__init__()
        OptModuleSettings.__init__(self)

        # measure time spend on reading/initialization
        self.timer.node["initialization"].start()

        self.T = T
        self.report_step = report_step
        self.customize_new_operator = customize_new_operator

        self.discr_type = discr_type
        self.physics_type = 'dead_oil'

        self.set_physics()
        self.set_reservoir(mesh_file)

        self.params.first_ts = 1.e-3
        self.params.mult_ts = 2
        self.params.max_ts = 100
        self.params.tolerance_newton = 1e-3
        self.params.tolerance_linear = 1e-6
        # self.params.newton_type = 2
        # self.params.newton_params = value_vector([0.2])

        self.timer.node["initialization"].stop()

    def set_reservoir(self, mesh_file):
        self.reservoir = UnstructReservoir(self.discr_type, mesh_file, n_vars=self.physics.n_vars)

        hcap = np.array(self.reservoir.mesh.heat_capacity, copy=False)
        rcond = np.array(self.reservoir.mesh.rock_cond, copy=False)
        hcap.fill(2200)
        rcond.fill(181.44)

        if self.discr_type == 'mpfa':
            self.reservoir.mesh.pz_bounds.resize(self.physics.n_vars * self.reservoir.n_bounds)
            pz_bounds = np.array(self.reservoir.mesh.pz_bounds, copy=False)
            pz_bounds[::self.reservoir.n_vars] = self.p_init
            for i in range(1, self.physics.nc):
                pz_bounds[i::self.reservoir.n_vars] = self.ini[i-1]

            pz_bounds[self.physics.n_vars-1::self.reservoir.n_vars] = self.init_temp

            self.reservoir.P_VAR = self.physics.engine.P_VAR

    def set_wells(self):
        # well model or boundary conditions
        if self.discr_type == 'mpfa':
            Lx = max([pt.values[0] for pt in self.reservoir.discr_mesh.nodes])
            Ly = max([pt.values[1] for pt in self.reservoir.discr_mesh.nodes])
            Lz = max([pt.values[2] for pt in self.reservoir.discr_mesh.nodes])
            dx = np.sqrt(np.mean(self.reservoir.volume_all_cells) / \
                    ( max([pt.values[2] for pt in self.reservoir.discr_mesh.nodes]) - min([pt.values[2] for pt in self.reservoir.discr_mesh.nodes])))

            n_cells = self.reservoir.discr_mesh.n_cells
            pt_x = np.array([c.values[0] for c in self.reservoir.discr_mesh.centroids])[:n_cells]
            pt_y = np.array([c.values[1] for c in self.reservoir.discr_mesh.centroids])[:n_cells]
            pt_z = np.array([c.values[2] for c in self.reservoir.discr_mesh.centroids])[:n_cells]

            x0 = 0.1 * Lx
            y0 = 0.1 * Ly
            self.id1 = ((pt_x - x0) ** 2 + (pt_y - y0) ** 2 + pt_z ** 2).argmin()

            x0 = 0.9 * Lx
            y0 = 0.9 * Ly
            self.id2 = ((pt_x - x0) ** 2 + (pt_y - y0) ** 2 + pt_z ** 2).argmin()
        else:
            pts = self.reservoir.unstr_discr.mesh_data.points
            Lx = np.max(pts[:,0])
            Ly = np.max(pts[:,1])
            Lz = np.max(pts[:,2])

            c = np.array([c.centroid for c in self.reservoir.unstr_discr.mat_cell_info_dict.values()])

            x0 = 0.1 * Lx
            y0 = 0.1 * Ly
            self.id1 = ((c[:,0] - x0) ** 2 + (c[:,1] - y0) ** 2 + c[:,2] ** 2).argmin()

            x0 = 0.9 * Lx
            y0 = 0.9 * Ly
            self.id2 = ((c[:,0] - x0) ** 2 + (c[:,1] - y0) ** 2 + c[:,2] ** 2).argmin()

        self.reservoir.add_well("P1", depth=0)
        self.reservoir.add_perforation(self.reservoir.wells[-1], int(self.id1), well_index=self.reservoir.well_index)

        self.reservoir.add_well("I1", depth=0)
        self.reservoir.add_perforation(self.reservoir.wells[-1], int(self.id2), well_index=self.reservoir.well_index)

    def set_physics(self):
        """Physical properties"""
        self.zero = 1e-13
        components = ['w', 'o']
        phases = ['wat', 'oil']
        property_container = ModelProperties(phases_name=phases, components_name=components, min_z=self.zero/10)

        self.cell_property = ['pressure'] + ['water']
        self.cell_property += ['temperature']

        # Define property evaluators based on custom properties
        property_container.density_ev = dict([('wat', DensityBasic(compr=1e-5, dens0=1014)),
                                              ('oil', DensityBasic(compr=5e-3, dens0=50))])
        property_container.viscosity_ev = dict([('wat', ConstFunc(0.3)),
                                                ('oil', ConstFunc(0.03))])
        property_container.rel_perm_ev = dict([('wat', PhaseRelPerm("gas", 0.1, 0.1)),
                                               ('oil', PhaseRelPerm("oil", 0.1, 0.1))])
        property_container.enthalpy_ev = dict([('wat', EnthalpyBasic(hcap=4.18)),
                                               ('oil', EnthalpyBasic(hcap=0.035))])
        property_container.conductivity_ev = dict([('wat', ConstFunc(1.)),
                                                   ('oil', ConstFunc(1.))])

        property_container.rock_energy_ev = EnthalpyBasic(hcap=1.0)

        # create physics
        thermal = True
        state_spec = Compositional.StateSpecification.PT if thermal else Compositional.StateSpecification.P
        self.physics = Compositional(components, phases, self.timer,
                                     n_points=400, min_p=0, max_p=1000, min_z=self.zero, max_z=1-self.zero,
                                     min_t=273.15 + 20, max_t=273.15 + 200, state_spec=state_spec)
        self.physics.add_property_region(property_container)

        self.runtime = 1000
        self.p_init = 200
        self.init_temp = 350
        self.inj_comp = value_vector([1 - self.zero])
        self.inj_temp = self.init_temp - 30
        self.ini = value_vector([self.zero])

    def set_initial_conditions(self):
        input_distribution = {'pressure': self.p_init, 'water': self.ini, 'temperature': self.init_temp}
        self.physics.set_initial_conditions_from_array(self.reservoir.mesh, input_distribution=input_distribution)

    def set_boundary_conditions(self):
        from darts.engines import well_control_iface
        for i, w in enumerate(self.reservoir.wells):
            if i == 0:
                self.physics.set_well_controls(well=w, is_control=True, is_inj=False, control_type=well_control_iface.BHP,
                                               target=self.p_init-10.)
            else:
                self.physics.set_well_controls(well=w, is_control=True, is_inj=True, control_type=well_control_iface.BHP,
                                               target=self.p_init+10., inj_composition=self.inj_comp, inj_temp=self.inj_temp)

    def init(self):
        """
        Function to initialize the model, which includes:
        - initialize well (perforation) position
        - initialize well rate parameters
        - initialize reservoir initial conditions
        - initialize well control settings
        - define list of operator interpolators for accumulation-flux regions and wells
        - initialize engine
        """
        self.set_boundary_conditions()
        self.reservoir.init_wells()
        self.physics.init_wells(self.reservoir.wells)
        self.set_initial_conditions()
        self.set_well_controls()
        self.set_op_list()
        self.reset()

    def run(self, export_to_vtk=False):
        import random
        # use seed to generate the same random values every run
        random.seed(3)
        if export_to_vtk:
            # Properties for writing to vtk format:
            # output_directory = 'trial_dir'  # Specify output directory here
            output_directory = 'sol_cpp_' + self.discr_type + '_{:s}'.format(self.physics_type)
            self.output_directory = output_directory
            # Write to vtk using class methods of unstructured discretizer (uses within meshio write to vtk function):
            if self.discr_type == 'mpfa':
                self.reservoir.write_to_vtk(output_directory, self.cell_property + ['perm'], 0, self.physics.engine)
            else:
                self.reservoir.write_to_vtk_old_discretizer(output_directory, self.cell_property, 0, self.physics.engine)

        ith_step = 0
        # now we start to run for the time report--------------------------------------------------------------
        time_step = self.report_step
        even_end = int(self.T / time_step) * time_step
        time_step_arr = np.ones(int(self.T / time_step)) * time_step
        if self.T - even_end > 0:
            time_step_arr = np.append(time_step_arr, self.T - even_end)

        for ts in time_step_arr:
            from darts.engines import well_control_iface
            for i, w in enumerate(self.reservoir.wells):
                if i == 0:
                    self.physics.set_well_controls(well=w, is_control=True, is_inj=False, control_type=well_control_iface.BHP,
                                                   target=self.p_init-10.)
                else:
                    self.physics.set_well_controls(well=w, is_control=True, is_inj=True, control_type=well_control_iface.BHP,
                                                   target=self.p_init+10., inj_composition=self.inj_comp, inj_temp=self.inj_temp)

            self.physics.engine.run(ts)
            self.physics.engine.report()
            if export_to_vtk:
                if self.discr_type == 'mpfa':
                    self.reservoir.write_to_vtk(output_directory, self.cell_property, ith_step + 1, self.physics.engine)
                else:
                    self.reservoir.write_to_vtk_old_discretizer(output_directory, self.cell_property, ith_step + 1, self.physics.engine)

                ith_step += 1

    def get_n_flux_multiplier(self):
        if self.discr_type == "mpfa":
            cell_m_mpfa = np.array(self.reservoir.mesh.block_m, copy=True)
            cell_p_mpfa = np.array(self.reservoir.mesh.block_p, copy=True)

            n_blocks = self.reservoir.mesh.n_blocks
            cell_m_one_way = []
            cell_p_one_way = []
            conn_index_to_one_way = []
            idx_one_way = 0
            for idx, cm in enumerate(cell_m_mpfa):
                cp = cell_p_mpfa[idx]
                # if cm < cp and cp < n_blocks:
                if cm < cp:
                    cell_m_one_way.append(cm)
                    cell_p_one_way.append(cp)
                    conn_index_to_one_way.append(idx_one_way)
                    idx_one_way += 1
                elif cm > cp:
                    cp_indices = np.where(np.array(cell_m_one_way) == cp)[0]  # note here we find cp from cm array
                    cm_indices = np.where(np.array(cell_p_one_way) == cm)[0]  # note here we find cm from cp array

                    # Find the intersection set of cp_indices and cm_indices
                    intersection_indices = np.intersect1d(cp_indices, cm_indices)

                    conn_index_to_one_way.append(intersection_indices[0])
                else:  # boundaries. This may be removed in the future
                    conn_index_to_one_way.append(-999)

            n_fm = len(cell_m_one_way) - len(self.reservoir.wells)  # minus the trans between well head and well body
            return n_fm


    def set_op_list(self):
        """
        Function to define list of operator interpolators for accumulation-flux regions and wells.

        Operator list is in order [acc_flux_itor[0], ..., acc_flux_itor[n-1], acc_flux_w_itor]
        """
        if type(self.physics.acc_flux_itor) == dict:
            self.op_list = [acc_flux_itor for acc_flux_itor in self.physics.acc_flux_itor.values()] + [self.physics.acc_flux_w_itor]
            self.op_num = np.array(self.reservoir.mesh.op_num, copy=False)
            # self.op_num[self.reservoir.nb:] = len(self.op_list) - 1
            self.op_num[self.reservoir.mesh.n_res_blocks:] = len(self.op_list) - 1

            if self.customize_new_operator:
                # customize your own operator, e.g. the Temperature
                temperature_etor = geothermal_customized_etor()

                temperature_itor = self.physics.create_interpolator(temperature_etor, timer_name="customized operator interpolation",
                                                                    n_ops=1, platform='cpu', algorithm='multilinear',
                                                                    mode='adaptive', precision='d')
                self.physics.create_itor_timers(temperature_itor, "customized operator interpolation")

                self.physics.engine.customize_operator = self.customize_new_operator

                self.op_list += [temperature_itor]

                # specify the index of blocks of customized operator
                op_num_new = np.array(self.reservoir.mesh.op_num, copy=True)
                op_num_new[:] = len(self.op_list) - 1  # set the last interpolator (i.e. "temperature_itor") from "self.op_list" to all blocks
                self.physics.engine.idx_customized_operator = len(self.op_list) - 1
                self.physics.engine.customize_op_num = index_vector(op_num_new)


class ModelProperties(PropertyContainer):
    def __init__(self, phases_name, components_name, min_z=1e-11):
        # Call base class constructor
        self.nph = len(phases_name)
        Mw = np.ones(self.nph)
        super().__init__(phases_name=phases_name, components_name=components_name, Mw=Mw, min_z=min_z, temperature=None)

    def evaluate(self, state):
        """
        Class methods which evaluates the state operators for the element based physics
        :param state: state variables [pres, comp_0, ..., comp_N-1]
        :param values: values of the operators (used for storing the operator values)
        :return: updated value for operators, stored in values
        """
        # Composition vector and pressure from state:
        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]

        zc = np.append(vec_state_as_np[1:self.nc], 1 - np.sum(vec_state_as_np[1:self.nc]))

        self.clean_arrays()
        # two-phase flash - assume water phase is always present and water component last
        for i in range(self.nph):
            self.x[i, i] = 1

        self.ph = [0, 1]

        for j in self.ph:
            # molar weight of mixture
            M = np.sum(self.x[j, :] * self.Mw)
            self.dens[j] = self.density_ev[self.phases_name[j]].evaluate(pressure)  # output in [kg/m3]
            self.dens_m[j] = self.dens[j] / M
            self.mu[j] = self.viscosity_ev[self.phases_name[j]].evaluate()  # output in [cp]

        self.nu = zc
        self.compute_saturation(self.ph)

        for j in self.ph:
            self.kr[j] = self.rel_perm_ev[self.phases_name[j]].evaluate(self.sat[j])
            self.pc[j] = 0

        mass_source = np.zeros(self.nc)

        return self.ph, self.sat, self.x, self.dens, self.dens_m, self.mu, self.kr, self.pc, mass_source

    def evaluate_at_cond(self, pressure, zc):

        self.sat[:] = 0

        ph = [0, 1]
        for j in ph:
            self.dens_m[j] = self.density_ev[self.phases_name[j]].evaluate(1, 0)

        self.dens_m = [1025, 0.77]  # to match DO based on PVT

        self.nu = zc
        self.compute_saturation(ph)

        return self.sat, self.dens_m


class geothermal_customized_etor(operator_set_evaluator_iface):
    def __init__(self):
        super().__init__()

    def evaluate(self, state, values):
        # state : [P, Z_1, ... Z_(NC-1), T]
        values[0] = np.array(state)[-1]

        return 0