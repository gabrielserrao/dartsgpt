import numpy as np
from darts.models.darts_model import DartsModel
from darts.reservoirs.struct_reservoir import StructReservoir
from physics import BrineVapour


class Model(DartsModel):
    def set_reservoir(self, nx, nz, dx=0.2, dz=0.2, poro=0.25, perm=100):
        """Reservoir"""
        nb = nx * nz

        depth = np.zeros(nb)
        for k in range(nz):
            depth[k * nx:(k + 1) * nx] = 0.5 * dz + k * dz

        self.reservoir = StructReservoir(self.timer, nx=nx, ny=1, nz=nz, dx=dx, dy=1, dz=dz,
                                         permx=perm, permy=perm, permz=perm, poro=poro, depth=depth)

    def set_physics(self, components: list, phases: list, salinity: float = 0., temperature: float = None,
                    n_points: int = 1001, zero: float = 1e-12, swc: float = 0.2, impurity: float = 0.):
        self.components = components
        self.salinity = salinity
        self.zero = zero
        self.swc = swc
        self.impurity = impurity

        self.physics = BrineVapour(components=components, phases=phases, ions=(salinity > 0.), swc=swc, timer=self.timer,
                                   n_points=n_points, min_p=50, max_p=300, min_t=273.15, max_t=373.15, zero=zero,
                                   temperature=temperature, cache=False)

    def set_initial_conditions(self):
        depths = np.asarray(self.reservoir.mesh.depth)
        min_depth = np.min(depths)
        max_depth = np.max(depths)
        nb = int(self.reservoir.nz)
        depths = np.linspace(min_depth, max_depth, nb)

        from darts.physics.super.initialize import Initialize
        init = Initialize(self.physics, aq_idx=0, h2o_idx=0)

        # Gas-water contact
        known_idx = int(nb / 2)
        mid_depth = depths[known_idx]

        nc = len(self.components)
        z = np.zeros((nb, nc))
        z[depths <= mid_depth, :] = [None for i in range(nc)]  # residual water is above GWC
        z[depths > mid_depth, :] = [1. - (nc-1) * self.zero] + [self.zero for i in range(1, nc-1)] + [None]  # liquid below GWC
        primary_specs = {var: z[:, i] for i, var in enumerate(self.components)}

        s = np.zeros(nb)
        s[depths <= mid_depth] = self.swc
        s[depths > mid_depth] = None
        secondary_specs = {'satAq': s}
        bound_specs = {'satAq': s[known_idx]}
        if nc > 2:
            y = np.zeros(nb)
            y[depths <= mid_depth] = self.impurity
            y[depths > mid_depth] = None
            secondary_specs.update({'x2V': y})
            bound_specs.update({'x2V': y[known_idx]})
        if self.salinity:
            mol = np.ones(nb) * self.salinity
            secondary_specs.update({'m' + str(nc): mol})
            bound_specs.update({'m' + str(nc): mol[known_idx]})
            primary_specs["H2O"][depths > mid_depth] = None
            primary_specs["CO2"][depths > mid_depth] = self.zero

        # Solve boundary state
        X0 = ([100., 0.5] +  # pressure, H2O
              ([0.5 * (1.-self.impurity)] if nc > 2 else []) +  # CO2
              ([1. - 0.5 - mol[known_idx] * 0.5 / 55.509] if self.salinity else []) +  # last component if ions
              ([300.] if init.thermal else []))  # temperature
        X0 = init.solve_state(X0, primary_specs={'pressure': 100., 'temperature': 300. if init.thermal else None},
                              secondary_specs=bound_specs)
        boundary_state = {v: X0[i] for i, v in enumerate(self.physics.vars)}

        # Solve vertical equilibrium
        X = init.solve(depth_bottom=max_depth, depth_top=min_depth, depth_known=mid_depth, nb=nb,
                       primary_specs=primary_specs, secondary_specs=secondary_specs,
                       boundary_state=boundary_state, dTdh=0.03).reshape((nb, self.physics.n_vars))

        self.physics.set_initial_conditions_from_depth_table(mesh=self.reservoir.mesh,
                                                             input_depth=init.depths,
                                                             input_distribution={v: X[:, i] for i, v in
                                                                                   enumerate(self.physics.vars)})
        return
