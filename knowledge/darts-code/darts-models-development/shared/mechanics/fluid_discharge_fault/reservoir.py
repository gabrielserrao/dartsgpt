from darts.engines import conn_mesh, ms_well, ms_well_vector, index_vector, value_vector, contact, contact_vector, vector_matrix
from darts.engines import matrix33, matrix, pm_discretizer, Face, vector_face_vector, face_vector, vector_matrix33, Stiffness, stf_vector
import numpy as np
from math import inf, pi
from darts.reservoirs.mesh.unstruct_discretizer import UnstructDiscretizer
from darts.reservoirs.mesh.geometrymodule import FType
from darts.engines import timer_node
from itertools import compress
import meshio
import os
import pandas as pd
from scipy.linalg import null_space
from matplotlib import pyplot as plt
from matplotlib import rcParams
from darts.reservoirs.mesh.transcalc import TransCalculations as TC
rcParams["text.usetex"]=False
# rcParams["font.sans-serif"] = ["Liberation Sans"]
# rcParams["font.serif"] = ["Liberation Serif"]
#
# # Definitions for the unstructured reservoir class:
class UnstructReservoir:
    def __init__(self, timer):
        self.timer = timer
        # Create mesh object (C++ object used by DARTS for all mesh related quantities):
        self.mesh = conn_mesh()

        # Specify elastic properties, mesh & boundaries
        self.mandel() # self.mandel() # self.porepy_cmp()

        # Discretization
        self.timer.node["discretization"] = timer_node()
        self.timer.node["discretization"].start()
        dt = 0.0
        self.pm.init(self.unstr_discr.mat_cells_tot, self.unstr_discr.frac_cells_tot, index_vector(self.ref_contact_cells))
        self.pm.reconstruct_gradients_per_cell(dt)
        self.pm.calc_all_fluxes_once(dt)

        # check sparsity of gradients
        # for cell_id in range(self.unstr_discr.mat_cells_tot):
        #     st, vals = self.pm.get_gradient(cell_id)
        #     st = np.array(st, dtype=np.intp)
        #     vals = np.array(vals).reshape(12, 4 * st.size)
        #     #assert((np.abs(vals[:9,3::4]) < 1.E-10).all())
        #     assert((np.abs(vals[9:12,0::4]) < 1.E-6).all())
        #     assert((np.abs(vals[9:12,1::4]) < 1.E-6).all())
        #     assert((np.abs(vals[9:12,2::4]) < 1.E-6).all())

        # check sparsity of coupled multi-point approximation
        # for i, cell_m in enumerate(self.pm.cell_m):
        #     cell_p = self.pm.cell_p[i]
        #     for k in range(self.pm.offset[i], self.pm.offset[i+1]):
        #         id = self.pm.stencil[k]
        #         tran = np.array(self.pm.tran[16*k:16*(k+1)]).reshape(4,4)
        #         tran_biot = np.array(self.pm.tran_biot[4 * k:4 * (k + 1)])
        #         assert((tran[3,:3] == 0.0).all()) # no displacements in flow
        #         #assert(tran_biot[3] == 0.0) # no pressures in volumetric strain
        #         #if id != cell_m and id != cell_p:
        #         #    assert(tran[3,3] == 0.0)    # TPFA

        self.mesh.init_pm(self.pm.cell_m, self.pm.cell_p, self.pm.stencil, self.pm.offset, self.pm.tran, self.pm.rhs,
                          self.pm.tran_biot, self.pm.rhs_biot,
                          self.unstr_discr.mat_cells_tot, self.unstr_discr.bound_cells_tot, self.unstr_discr.frac_cells_tot)
        #self.write_pm_conn_to_file(t_step=0)
        self.timer.node["discretization"].stop()

        self.unstr_discr.store_volume_all_cells()
        self.unstr_discr.store_depth_all_cells()

        # Create numpy arrays wrapped around mesh data (no copying, this will severely slow down the process!)
        self.poro = np.array(self.mesh.poro, copy=False)
        self.depth = np.array(self.mesh.depth, copy=False)
        self.volume = np.array(self.mesh.volume, copy=False)
        self.bc = np.array(self.mesh.bc, copy=False)
        self.bc_prev = np.array(self.mesh.bc_prev, copy=False)
        self.bc_ref = np.array(self.mesh.bc_ref, copy=False)
        self.mesh.f.resize(4 * (self.unstr_discr.frac_cells_tot + self.unstr_discr.mat_cells_tot))
        self.f = np.array(self.mesh.f, copy=False)
        self.biot_arr = np.array(self.mesh.biot, copy=False)
        self.kd = np.array(self.mesh.kd, copy=False)
        self.mesh.pz_bounds.resize(self.unstr_discr.bound_cells_tot)
        self.pz_bounds = np.array(self.mesh.pz_bounds, copy=False)
        self.p_ref = np.array(self.mesh.ref_pressure, copy=False)

        # Since we use copy==False above, we have to store the values by using the Python slicing option, if we don't
        # do this we will overwrite the variable, e.g. self.poro = poro --> overwrite self.poro with the variable poro
        # instead of storing the variable poro in self.mesh.poro (therefore "numpy array wrapped around mesh data!!!):
        self.poro[:] = self.porosity
        self.depth[:] = self.unstr_discr.depth_all_cells#[:self.unstr_discr.frac_cells_tot + self.unstr_discr.mat_cells_tot]
        self.volume[:self.unstr_discr.mat_cells_tot] = self.unstr_discr.volume_all_cells[self.unstr_discr.frac_cells_tot:]
        for i in range(self.unstr_discr.mat_cells_tot, self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot):
            self.volume[i] = self.unstr_discr.faces[i][4].area * self.frac_apers[i-self.unstr_discr.mat_cells_tot]
        self.bc_prev[:] = self.bc_rhs_prev
        self.bc[:] = self.bc_rhs
        self.bc_ref[:] = self.bc_rhs_ref
        self.biot_arr[:] = np.tile([self.biot,0,0,
                                    0,self.biot,0,
                                    0,0,self.biot], self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot)
        self.kd[:] = self.lam + 2 * self.mu / 3
        self.pz_bounds[:] = self.unstr_discr.pz_bounds
        self.p_ref[:] = self.unstr_discr.p_ref
        self.f[:] = self.unstr_discr.f

        # Calculate well_index (very primitive way....):
        rw = 0.1
        coords = self.unstr_discr.mat_cell_info_dict[0].coord_nodes_to_cell
        dx = np.max(coords[:,0]) - np.min(coords[:,0])
        dy = np.max(coords[:,1]) - np.min(coords[:,1])
        dz = np.max(coords[:,2]) - np.min(coords[:,2])
        # WIx
        wi_x = 0.0
        # WIy
        wi_y = 0.0
        # WIz
        hz = dz
        rp_z = 0.28 * np.sqrt((self.permy / self.permx) ** 0.5 * dx ** 2 +
                              (self.permx / self.permy) ** 0.5 * dy ** 2) / \
               ((self.permx / self.permy) ** 0.25 + (self.permy / self.permx) ** 0.25)
        wi_z = 2 * np.pi * np.sqrt(self.permx * self.permy) * hz / np.log(rp_z / rw)
        self.well_index = TC.darcy_constant * np.sqrt(wi_x ** 2 + wi_y ** 2 + wi_z ** 2)
        # Create empty list of wells:
        self.wells = []
    def set_equilibrium(self):
        # store original transmissibilities
        self.tran = np.array(self.mesh.tran, copy=True)
        self.rhs = np.array(self.mesh.rhs, copy=True)
        self.tran_biot = np.array(self.mesh.tran_biot, copy=True)
        self.rhs_biot = np.array(self.mesh.rhs_biot, copy=True)
        # turn off some terms for evaluation of momentum equilibrium
        tran = np.array(self.mesh.tran, copy=False)
        rhs = np.array(self.mesh.rhs, copy=False)
        tran_biot = np.array(self.mesh.tran_biot, copy=False)
        rhs_biot = np.array(self.mesh.rhs_biot, copy=False)

        tran[12::16] = 0.0
        tran[13::16] = 0.0
        tran[14::16] = 0.0
        tran[15::16] = 0.0
        rhs[3::4] = 0.0

        tran_biot[12::16] = 0.0
        tran_biot[13::16] = 0.0
        tran_biot[14::16] = 0.0
        tran_biot[15::16] = 0.0
        rhs_biot[3::4] = 0.0

        self.unstr_discr.f[3::4] = 0.0#self.p_init - self.unstr_discr.p_ref[:]
        self.f[:] = self.unstr_discr.f
    def turn_off_equilibrium(self):
        tran = np.array(self.mesh.tran, copy=False)
        rhs = np.array(self.mesh.rhs, copy=False)
        tran_biot = np.array(self.mesh.tran_biot, copy=False)
        rhs_biot = np.array(self.mesh.rhs_biot, copy=False)
        tran[12::16] = self.tran[12::16]
        tran[13::16] = self.tran[13::16]
        tran[14::16] = self.tran[14::16]
        tran[15::16] = self.tran[15::16]
        rhs[3::4] = self.rhs[3::4]

        tran_biot[12::16] = self.tran_biot[12::16]
        tran_biot[13::16] = self.tran_biot[13::16]
        tran_biot[14::16] = self.tran_biot[14::16]
        tran_biot[15::16] = self.tran_biot[15::16]
        rhs_biot[3::4] = self.rhs_biot[3::4]

        self.unstr_discr.f[:] = 0.0
        self.f[:] = 0.0
    def apply_geomehcanics_mode(self, physics):
        geom_mode = np.array(physics.engine.geomechanics_mode, copy=False)
        geom_mode[:] = 1

        cell_m = np.array(self.mesh.block_m, copy=False)
        cell_p = np.array(self.mesh.block_p, copy=False)
        offset = np.array(self.mesh.offset, copy=False)
        tran = np.array(self.mesh.tran, copy=False)
        tran_biot = np.array(self.mesh.tran_biot, copy=False)
        rhs = np.array(self.mesh.rhs, copy=False)
        rhs_biot = np.array(self.mesh.rhs_biot, copy=False)

        # tran[12::16] = 0.0
        # tran[13::16] = 0.0
        # tran[14::16] = 0.0
        # tran[15::16] = 0.0
        # rhs[3::4] = 0.0

        tran_biot[12::16] = 0.0
        tran_biot[13::16] = 0.0
        tran_biot[14::16] = 0.0
        tran_biot[15::16] = 0.0
        rhs_biot[3::4] = 0.0

        # for conn_id in range(len(cell_m)):
        #     cell_id1 = cell_m[conn_id]
        #     if cell_id1 < self.unstr_discr.mat_cells_tot:
        #         if self.unstr_discr.mat_cell_info_dict[cell_id1].prop_id == 99992:
        #             tran[12+16 * offset[conn_id]:12+16 * offset[conn_id + 1]:16] = 0.0
        #             tran[13+16 * offset[conn_id]:13+16 * offset[conn_id + 1]:16] = 0.0
        #             tran[14+16 * offset[conn_id]:14+16 * offset[conn_id + 1]:16] = 0.0
        #             tran[15+16 * offset[conn_id]:15+16 * offset[conn_id + 1]:16] = 0.0
        #             rhs[3 + 4 * cell_id1] = 0.0
        #             tran_biot[12+16 * offset[conn_id]:12+16 * offset[conn_id + 1]:16] = 0.0
        #             tran_biot[13+16 * offset[conn_id]:13+16 * offset[conn_id + 1]:16] = 0.0
        #             tran_biot[14+16 * offset[conn_id]:14+16 * offset[conn_id + 1]:16] = 0.0
        #             tran_biot[15+16 * offset[conn_id]:15+16 * offset[conn_id + 1]:16] = 0.0
        #             rhs_biot[3 + 4 * cell_id1] = 0.0
        #             continue
        #
        #     cell_id2 = cell_p[conn_id]
        #     if cell_id2 < self.unstr_discr.mat_cells_tot:
        #         if self.unstr_discr.mat_cell_info_dict[cell_id2].prop_id == 99992:
        #             tran[12+16 * offset[conn_id]:12+16 * offset[conn_id + 1]:16] = 0.0
        #             tran[13+16 * offset[conn_id]:13+16 * offset[conn_id + 1]:16] = 0.0
        #             tran[14+16 * offset[conn_id]:14+16 * offset[conn_id + 1]:16] = 0.0
        #             tran[15+16 * offset[conn_id]:15+16 * offset[conn_id + 1]:16] = 0.0
        #             rhs[3 + 4 * cell_id2] = 0.0
        #             tran_biot[12+16 * offset[conn_id]:12+16 * offset[conn_id + 1]:16] = 0.0
        #             tran_biot[13+16 * offset[conn_id]:13+16 * offset[conn_id + 1]:16] = 0.0
        #             tran_biot[14+16 * offset[conn_id]:14+16 * offset[conn_id + 1]:16] = 0.0
        #             tran_biot[15+16 * offset[conn_id]:15+16 * offset[conn_id + 1]:16] = 0.0
        #             rhs_biot[3 + 4 * cell_id2] = 0.0
        #             continue
    def update(self, dt, time):
        # update local array
        #if time > dt:
        self.bc_rhs_prev = np.copy(self.bc_rhs)
        self.pm.bc_prev = self.pm.bc
    def update_trans(self, dt, x):
        #self.pm.x_prev = value_vector(np.concatenate((x, self.bc_rhs_prev)))
        #self.pm.reconstruct_gradients_per_cell(dt)
        #self.pm.calc_all_fluxes(dt)
        #self.write_pm_conn_to_file(t_step=t_step)
        #self.mesh.init_pm(self.pm.cell_m, self.pm.cell_p, self.pm.stencil, self.pm.offset, self.pm.tran, self.pm.rhs,
        #                  self.unstr_discr.mat_cells_tot, self.unstr_discr.bound_cells_tot, 0)

        # update transient sources / sinks
        self.f[:] = self.unstr_discr.f
        # update boundaries at n+1 / n timesteps
        self.bc[:] = self.bc_rhs
        self.bc_prev[:] = self.bc_rhs_prev
        #self.init_wells()

    def porepy_cmp(self):
        n_dim = 3
        self.u_init = [0.0, 0.0, 0.0]
        self.p_init0 = 100.0
        self.porosity = 0.15
        self.permx = self.permy = self.permz = 1.0
        mesh_file = 'meshes/wedge_frac.msh'
        self.file_path = mesh_file
        self.unstr_discr = UnstructDiscretizer(permx=self.permx, permy=self.permy, permz=self.permz, frac_aper=0,
                                               mesh_file=mesh_file)
        self.unstr_discr.eps_t = 1.E+0
        self.unstr_discr.eps_n = 1.E+0
        self.unstr_discr.mu = 3.2
        self.unstr_discr.P12 = 0
        self.unstr_discr.Prol = 1
        self.unstr_discr.n_dim = 3
        self.unstr_discr.bcf_num = 3
        self.unstr_discr.bcm_num = self.unstr_discr.n_dim + 3
        #lam = 1.0 * 10000  # in bar
        #mu = 1.0 * 10000
        #nu = lam / 2 / (lam + mu)
        #E = lam * (1 + nu) * (1 - 2 * nu) / nu

        #self.mu = 65000 # in bars
        #nu = 0.15
        self.lam = 1.0#E * nu / (1 + nu) / (1 - 2 * nu)
        self.mu = 1.0#E / 2 / (1 + nu)
        E = self.mu * (3 * self.lam + 2 * self.mu) / (self.lam + self.mu)
        nu = self.lam / (self.lam + self.mu) / 2

        self.unstr_discr.init_matrix_stiffness({99991: {'E': E, 'nu': nu}})
        self.unstr_discr.physical_tags['matrix'] = [99991]
        self.unstr_discr.physical_tags['fracture_shape'] = [1, 2]
        self.unstr_discr.physical_tags['fracture'] = [91]
        self.unstr_discr.physical_tags['output'] = []
        self.unstr_discr.physical_tags['boundary'] = [991, 992, 993, 994, 995, 996]
        # General representation of BC: a*p + b*f = r (a=1,b=0 - Dirichlet, a=0,b=1 - Neumann)

        NO_FLOW = {'a': 0.0, 'b': 1.0, 'r': 0.0}
        AQUIFER = lambda p: {'a': 1.0, 'b': 0.0, 'r': p}
        ROLLER =    {'an': 1.0, 'bn': 0.0, 'rn': 0.0, 'at': 0.0, 'bt': 1.0, 'rt': np.array([0, 0, 0])}
        FREE =      {'an': 0.0, 'bn': 1.0, 'rn': 0.0, 'at': 0.0, 'bt': 1.0, 'rt': np.array([0, 0, 0])}
        STUCK = lambda un, ut: {'an': 1.0, 'bn': 0.0, 'rn': un, 'at': 1.0, 'bt': 0.0, 'rt': np.array(ut)}
        LOAD = lambda Fn, Ft: {'an': 0.0, 'bn': 1.0, 'rn': Fn, 'at': 0.0, 'bt': 1.0, 'rt': np.array(Ft)}
        STUCK_ROLLER = lambda un: {'an': 1.0, 'bn': 0.0, 'rn': un, 'at': 0.0, 'bt': 1.0, 'rt': np.array([0.0, 0.0, 0.0])}
        STUCK_T_LOAD_N = lambda Fn, ut: {'an': 0.0, 'bn': 1.0, 'rn': Fn, 'at': 1.0, 'bt': 0.0, 'rt': np.array(ut)}

        self.unstr_discr.boundary_conditions[991] = {'flow': NO_FLOW,           'mech': FREE, 'cells': []}
        self.unstr_discr.boundary_conditions[992] = {'flow': NO_FLOW,           'mech': FREE, 'cells': []}
        self.unstr_discr.boundary_conditions[993] = {'flow': NO_FLOW,           'mech': STUCK(0.0, [0.0, 0.0, 0.0]), 'cells': []}
        self.unstr_discr.boundary_conditions[994] = {'flow': NO_FLOW,           'mech': STUCK(-0.001, [0.01, 0.0, 0.0]), 'cells': []}
        self.unstr_discr.boundary_conditions[995] = {'flow': NO_FLOW,           'mech': ROLLER, 'cells': []}
        self.unstr_discr.boundary_conditions[996] = {'flow': NO_FLOW,           'mech': ROLLER, 'cells': []}
        self.unstr_discr.boundary_conditions[1] =   {'flow': {'a': 0.0, 'b': 1.0, 'r': 0.0},
                                                        'mech': {'an': 1.0, 'bn': 0.0, 'rn': 0.0, 'at': 1.0, 'bt': 0.0, 'rt': np.array([0.0, 0.0, 0.0])}, 'cells': [] }
        self.unstr_discr.boundary_conditions[2] =   {'flow': {'a': 0.0, 'b': 1.0, 'r': 0.0},
                                                        'mech': {'an': 1.0, 'bn': 0.0, 'rn': 0.0, 'at': 0.0, 'bt': 1.0, 'rt': np.array([0.0, 0.0, 0.0])}, 'cells': [] }
        self.unstr_discr.fracture_aperture = 1
        self.unstr_discr.load_mesh_with_bounds()
        self.unstr_discr.calc_cell_neighbours()

        # init poromechanics discretizer
        self.pm = pm_discretizer()
        self.pm.grav = matrix([0.0, 0.0, 0.0], 1, 3)
        self.pm.visc = 1#9.81e-2
        self.biot = 1
        self.unstr_discr.f = np.zeros(4 * (self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot))
        self.p_init = np.zeros(self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot)
        self.unstr_discr.p_ref = np.zeros(self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot)
        for cell_id in range(self.unstr_discr.mat_cells_tot):
            faces = self.unstr_discr.faces[cell_id]
            fs = face_vector()
            for face_id in range(len(faces)):
                face = faces[face_id]
                fs.append(Face(face.type.value, face.cell_id1, face.cell_id2,
                                        face.face_id1, face.face_id2,
                                        face.area, list(face.n), list(face.centroid)))
            self.pm.faces.append(fs)

            cell = self.unstr_discr.mat_cell_info_dict[cell_id]

            self.pm.cell_centers.append(matrix(list(cell.centroid), cell.centroid.size, 1))
            self.pm.perms.append(matrix33(self.permx, self.permy, self.permz))
            self.pm.biots.append(matrix33(self.biot))

            E = self.unstr_discr.E[cell.prop_id]
            nu = self.unstr_discr.nu[cell.prop_id]
            lam = E * nu / (1 + nu) / (1 - 2 * nu)
            mu = E / (1 + nu) / 2
            self.pm.stfs.append(Stiffness(lam, mu))
        # fracture
        self.frac_aper = 1.E-5
        self.frac_apers = self.frac_aper * np.ones(self.unstr_discr.frac_cells_tot)
        for frac_id in range(self.unstr_discr.mat_cells_tot,
                             self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot):
            frac = self.unstr_discr.frac_cell_info_dict[frac_id]
            faces = self.unstr_discr.faces[frac_id]
            self.pm.cell_centers.append(matrix(list(frac.centroid), frac.centroid.size, 1))
            self.pm.frac_apers.append(self.frac_apers[frac_id - self.unstr_discr.mat_cells_tot])
            fs = face_vector()
            for face_id, face in faces.items():
                face = faces[face_id]
                fs.append(Face(face.type.value, face.cell_id1, face.cell_id2,
                               face.face_id1, face.face_id2,
                               face.area, list(face.n), list(face.centroid)))
            self.pm.faces.append(fs)

            face1 = faces[4]
            face2 = faces[5]
            self.mesh.fault_normals.append(face1.n[0])
            self.mesh.fault_normals.append(face1.n[1])
            self.mesh.fault_normals.append(face1.n[2])
            # Local basis
            S = np.zeros((n_dim, n_dim))
            S[:n_dim - 1] = null_space(np.array([face1.n])).T
            S[n_dim - 1] = face1.n
            Sinv = np.linalg.inv(S)
            K = np.zeros((n_dim, n_dim))
            K[0, 0] = K[1, 1] = self.permx
            K[n_dim - 1, n_dim - 1] = self.permx
            K = Sinv.dot(K).dot(S)
            self.pm.perms.append(matrix33(list(K.flatten())))
            self.pm.biots.append(matrix33(self.biot))
        # contact
        self.contacts = contact_vector()
        for tag in self.unstr_discr.physical_tags['fracture']:
            con = contact()
            con.f_scale = 1.E+4
            cell_ids = [cell_id for cell_id, cell in self.unstr_discr.frac_cell_info_dict.items() if
                        cell.prop_id == tag]
            fric_coef = 0.7 * np.ones(len(cell_ids))
            con.init_geometry(tag, self.pm, self.mesh, index_vector(cell_ids))
            con.init_friction(value_vector(fric_coef))
            self.contacts.append(con)

        self.unstr_discr.p_ref[:] = self.p_init0
        self.p_init[:] = self.p_init0
        self.bc_rhs_ref = np.zeros(4 * len(self.unstr_discr.bound_cell_info_dict))
        self.bc_rhs = np.zeros(4 * len(self.unstr_discr.bound_cell_info_dict))
        self.bc_rhs_prev = np.zeros(4 * len(self.unstr_discr.bound_cell_info_dict))
        self.unstr_discr.pz_bounds = np.zeros(self.unstr_discr.bound_cells_tot)
        self.unstr_discr.pz_bounds[:] = self.p_init0
        for bound_id in range(len(self.unstr_discr.bound_cell_info_dict)):
            b_cell = self.unstr_discr.bound_cell_info_dict[bound_id]
            prop_id = b_cell.prop_id
            n = self.get_normal_to_bound_face(bound_id)
            P = np.identity(3) - np.outer(n, n)
            mech = self.unstr_discr.boundary_conditions[prop_id]['mech']
            flow = self.unstr_discr.boundary_conditions[prop_id]['flow']
            bc = [mech['an'], mech['bn'], mech['at'], mech['bt'], flow['a'], flow['b']]
            self.pm.bc.append(matrix(bc, len(bc), 1))
            self.bc_rhs[4 * bound_id:4 * bound_id + 3] = mech['rn'] * n + mech['rt']
            self.bc_rhs[4 * bound_id + 3] = flow['r']
            self.bc_rhs_prev[4 * bound_id:4 * bound_id + 3] = np.array([0, 0, 0])
            self.bc_rhs_prev[4 * bound_id + 3] = flow['r']
            self.bc_rhs_ref[4 * bound_id:4 * bound_id + 3] = np.array([0, 0, 0])
            self.bc_rhs_ref[4 * bound_id + 3] = flow['r']
        for bound_id in range(len(self.unstr_discr.bound_cell_info_dict), self.unstr_discr.bound_cell_count):
            mech = self.unstr_discr.boundary_conditions[self.unstr_discr.frac_bound_cell_info_dict[bound_id].prop_id]['mech']
            flow = self.unstr_discr.boundary_conditions[self.unstr_discr.frac_bound_cell_info_dict[bound_id].prop_id]['flow']
            bc = [mech['an'], mech['bn'], mech['at'], mech['bt'], flow['a'], flow['b']]
            self.pm.bc.append(matrix(bc, len(bc), 1))
        self.bc_rhs_prev = np.copy(self.bc_rhs)
        self.pm.bc_prev = self.pm.bc
    def mandel(self):
        n_dim = 3
        self.u_init = [0.0, 0.0, 0.0]
        self.p_init0 = 0.0
        self.porosity = 0.15
        self.permx = self.permy = self.permz = 1.0
        mesh_file = 'meshes/wedge_frac.msh'
        self.file_path = mesh_file
        self.unstr_discr = UnstructDiscretizer(permx=self.permx, permy=self.permy, permz=self.permz, frac_aper=0,
                                               mesh_file=mesh_file)
        self.unstr_discr.eps_t = 1.E+0
        self.unstr_discr.eps_n = 1.E+0
        self.unstr_discr.mu = 3.2
        self.unstr_discr.P12 = 0
        self.unstr_discr.Prol = 1
        self.unstr_discr.n_dim = 3
        self.unstr_discr.bcf_num = 3
        self.unstr_discr.bcm_num = self.unstr_discr.n_dim + 3

        E = 10000 # in bars
        nu = 0.25
        self.lam = E * nu / (1 + nu) / (1 - 2 * nu)
        self.mu = E / 2 / (1 + nu)
        #E = self.mu * (3 * self.lam + 2 * self.mu) / (self.lam + self.mu)
        #nu = self.lam / (self.lam + self.mu) / 2

        self.unstr_discr.init_matrix_stiffness({99991: {'E': E, 'nu': nu}})
        self.unstr_discr.physical_tags['matrix'] = [99991]
        self.unstr_discr.physical_tags['fracture_shape'] = [1, 2]
        self.unstr_discr.physical_tags['fracture'] = [91]
        self.unstr_discr.physical_tags['output'] = []
        self.unstr_discr.physical_tags['boundary'] = [991, 992, 993, 994, 995, 996]
        # General representation of BC: a*p + b*f = r (a=1,b=0 - Dirichlet, a=0,b=1 - Neumann)

        NO_FLOW = {'a': 0.0, 'b': 1.0, 'r': 0.0}
        AQUIFER = lambda p: {'a': 1.0, 'b': 0.0, 'r': p}
        ROLLER =    {'an': 1.0, 'bn': 0.0, 'rn': 0.0, 'at': 0.0, 'bt': 1.0, 'rt': np.array([0, 0, 0])}
        FREE =      {'an': 0.0, 'bn': 1.0, 'rn': 0.0, 'at': 0.0, 'bt': 1.0, 'rt': np.array([0, 0, 0])}
        STUCK = lambda un, ut: {'an': 1.0, 'bn': 0.0, 'rn': un, 'at': 1.0, 'bt': 0.0, 'rt': np.array(ut)}
        LOAD = lambda Fn, Ft: {'an': 0.0, 'bn': 1.0, 'rn': Fn, 'at': 0.0, 'bt': 1.0, 'rt': np.array(Ft)}
        STUCK_ROLLER = lambda un: {'an': 1.0, 'bn': 0.0, 'rn': un, 'at': 0.0, 'bt': 1.0, 'rt': np.array([0.0, 0.0, 0.0])}
        STUCK_T_LOAD_N = lambda Fn, ut: {'an': 0.0, 'bn': 1.0, 'rn': Fn, 'at': 1.0, 'bt': 0.0, 'rt': np.array(ut)}

        self.unstr_discr.boundary_conditions[991] = {'flow': NO_FLOW,               'mech': ROLLER, 'cells': []}
        self.unstr_discr.boundary_conditions[992] = {'flow': AQUIFER(self.p_init0),  'mech': FREE, 'cells': []}
        self.unstr_discr.boundary_conditions[993] = {'flow': NO_FLOW,               'mech': ROLLER, 'cells': []}
        self.unstr_discr.boundary_conditions[994] = {'flow': NO_FLOW,               'mech': LOAD(-10, [0.0, 0.0, 0.0]), 'cells': []}
        self.unstr_discr.boundary_conditions[995] = {'flow': NO_FLOW,               'mech': ROLLER, 'cells': []}
        self.unstr_discr.boundary_conditions[996] = {'flow': NO_FLOW,               'mech': ROLLER, 'cells': []}
        self.unstr_discr.boundary_conditions[1] =   {'flow': {'a': 0.0, 'b': 1.0, 'r': 0.0},
                                                        'mech': {'an': 1.0, 'bn': 0.0, 'rn': 0.0, 'at': 1.0, 'bt': 0.0, 'rt': np.array([0.0, 0.0, 0.0])}, 'cells': [] }
        self.unstr_discr.boundary_conditions[2] =   {'flow': {'a': 0.0, 'b': 1.0, 'r': 0.0},
                                                        'mech': {'an': 1.0, 'bn': 0.0, 'rn': 0.0, 'at': 0.0, 'bt': 1.0, 'rt': np.array([0.0, 0.0, 0.0])}, 'cells': [] }
        self.unstr_discr.fracture_aperture = 1
        self.unstr_discr.load_mesh_with_bounds()
        self.unstr_discr.calc_cell_neighbours()

        # init poromechanics discretizer
        self.pm = pm_discretizer()
        self.pm.grav = matrix([0.0, 0.0, 0.0], 1, 3)
        self.fluid_compressibility = 1.e-5
        self.fluid_viscosity = 1.0
        self.pm.visc = 1#9.81e-2
        self.biot = 1
        self.unstr_discr.f = np.zeros(4 * (self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot))
        self.p_init = np.zeros(self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot)
        self.unstr_discr.p_ref = np.zeros(self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot)
        for cell_id in range(self.unstr_discr.mat_cells_tot):
            faces = self.unstr_discr.faces[cell_id]
            fs = face_vector()
            for face_id in range(len(faces)):
                face = faces[face_id]
                fs.append(Face(face.type.value, face.cell_id1, face.cell_id2,
                                        face.face_id1, face.face_id2,
                                        face.area, list(face.n), list(face.centroid)))
            self.pm.faces.append(fs)

            cell = self.unstr_discr.mat_cell_info_dict[cell_id]

            self.pm.cell_centers.append(matrix(list(cell.centroid), cell.centroid.size, 1))
            self.pm.perms.append(matrix33(self.permx, self.permy, self.permz))
            self.pm.biots.append(matrix33(self.biot))

            E = self.unstr_discr.E[cell.prop_id]
            nu = self.unstr_discr.nu[cell.prop_id]
            lam = E * nu / (1 + nu) / (1 - 2 * nu)
            mu = E / (1 + nu) / 2
            self.pm.stfs.append(Stiffness(lam, mu))
        # fracture
        self.frac_aper = 1.E-5
        self.frac_apers = self.frac_aper * np.ones(self.unstr_discr.frac_cells_tot)
        for frac_id in range(self.unstr_discr.mat_cells_tot,
                             self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot):
            frac = self.unstr_discr.frac_cell_info_dict[frac_id]
            faces = self.unstr_discr.faces[frac_id]
            self.pm.cell_centers.append(matrix(list(frac.centroid), frac.centroid.size, 1))
            self.pm.frac_apers.append(self.frac_apers[frac_id - self.unstr_discr.mat_cells_tot])
            fs = face_vector()
            for face_id, face in faces.items():
                face = faces[face_id]
                fs.append(Face(face.type.value, face.cell_id1, face.cell_id2,
                               face.face_id1, face.face_id2,
                               face.area, list(face.n), list(face.centroid)))
            self.pm.faces.append(fs)

            face1 = faces[4]
            face2 = faces[5]
            self.mesh.fault_normals.append(face1.n[0])
            self.mesh.fault_normals.append(face1.n[1])
            self.mesh.fault_normals.append(face1.n[2])
            # Local basis
            S = np.zeros((n_dim, n_dim))
            S[:n_dim - 1] = null_space(np.array([face1.n])).T
            S[n_dim - 1] = face1.n
            Sinv = np.linalg.inv(S)
            K = np.zeros((n_dim, n_dim))
            K[0, 0] = K[1, 1] = self.permx
            K[n_dim - 1, n_dim - 1] = self.permx
            K = Sinv.dot(K).dot(S)
            self.pm.perms.append(matrix33(list(K.flatten())))
            self.pm.biots.append(matrix33(self.biot))
        # contact
        self.ref_contact_cells = np.zeros(self.unstr_discr.frac_cells_tot, dtype=np.intc)
        self.contacts = contact_vector()
        for tag in self.unstr_discr.physical_tags['fracture']:
            con = contact()
            con.f_scale = 1.E+4
            cell_ids = [cell_id for cell_id, cell in self.unstr_discr.frac_cell_info_dict.items() if
                        cell.prop_id == tag]
            fric_coef = 0.7 * np.ones(len(cell_ids))
            self.ref_contact_cells[np.array(cell_ids, dtype=np.intp) - self.unstr_discr.mat_cells_tot] = cell_ids[0]
            con.mu0 = value_vector(fric_coef)
            con.mu = con.mu0
            con.cell_ids = index_vector(cell_ids)
            con.fault_tag = tag
            self.contacts.append(con)

        self.unstr_discr.p_ref[:] = self.p_init0
        self.p_init[:] = self.p_init0
        self.bc_rhs_ref = np.zeros(4 * len(self.unstr_discr.bound_cell_info_dict))
        self.bc_rhs = np.zeros(4 * len(self.unstr_discr.bound_cell_info_dict))
        self.bc_rhs_prev = np.zeros(4 * len(self.unstr_discr.bound_cell_info_dict))
        self.unstr_discr.pz_bounds = np.zeros(self.unstr_discr.bound_cells_tot)
        self.unstr_discr.pz_bounds[:] = self.p_init0
        for bound_id in range(len(self.unstr_discr.bound_cell_info_dict)):
            b_cell = self.unstr_discr.bound_cell_info_dict[bound_id]
            prop_id = b_cell.prop_id
            n = self.get_normal_to_bound_face(bound_id)
            P = np.identity(3) - np.outer(n, n)
            mech = self.unstr_discr.boundary_conditions[prop_id]['mech']
            flow = self.unstr_discr.boundary_conditions[prop_id]['flow']
            bc = [mech['an'], mech['bn'], mech['at'], mech['bt'], flow['a'], flow['b']]
            self.pm.bc.append(matrix(bc, len(bc), 1))
            self.bc_rhs[4 * bound_id:4 * bound_id + 3] = mech['rn'] * n + mech['rt']
            self.bc_rhs[4 * bound_id + 3] = flow['r']
            self.bc_rhs_prev[4 * bound_id:4 * bound_id + 3] = np.array([0, 0, 0])
            self.bc_rhs_prev[4 * bound_id + 3] = flow['r']
            self.bc_rhs_ref[4 * bound_id:4 * bound_id + 3] = np.array([0, 0, 0])
            self.bc_rhs_ref[4 * bound_id + 3] = flow['r']
        for bound_id in range(len(self.unstr_discr.bound_cell_info_dict), self.unstr_discr.bound_cell_count):
            mech = self.unstr_discr.boundary_conditions[self.unstr_discr.frac_bound_cell_info_dict[bound_id].prop_id]['mech']
            flow = self.unstr_discr.boundary_conditions[self.unstr_discr.frac_bound_cell_info_dict[bound_id].prop_id]['flow']
            bc = [mech['an'], mech['bn'], mech['at'], mech['bt'], flow['a'], flow['b']]
            self.pm.bc.append(matrix(bc, len(bc), 1))
        self.bc_rhs_prev = np.copy(self.bc_rhs)
        self.pm.bc_prev = self.pm.bc

    def add_well(self, name, depth):
        """
        Class method which adds wells heads to the reservoir (Note: well head is not equal to a perforation!)
        :param name:
        :param depth:
        :return:
        """
        well = ms_well()
        well.name = name
        well.segment_volume = 0.0785 * 40  # 2.5 * pi * 0.15**2 / 4
        well.well_head_depth = depth
        well.well_body_depth = depth
        well.segment_transmissibility = 1e5
        well.segment_depth_increment = 1
        self.wells.append(well)
        return 0
    def add_perforation(self, well, res_block, well_index):
        """
        Class method which ads perforation to each (existing!) well
        :param well: data object which contains data of the particular well
        :param res_block: reservoir block in which the well has a perforation
        :param well_index: well index (productivity index)
        :return:
        """
        well_block = 0
        well.perforations = well.perforations + [(well_block, res_block, well_index)]
        return 0
    def init_wells(self):
        """
        Class method which initializes the wells (adding wells and their perforations to the reservoir)
        :return:
        """
        # Add injection well:
        # self.add_well("I1", 0.5)
        # # Perforate all boundary cells:
        # for nth_perf in range(len(self.left_boundary_cells)):
        #     well_index = self.mesh.volume[self.left_boundary_cells[nth_perf]] / self.max_well_vol * self.well_index
        #     self.add_perforation(well=self.wells[-1], res_block=self.left_boundary_cells[nth_perf],
        #                          well_index=well_index)
        #
        # # Add production well:
        # self.add_well("P1", 0.5)
        # # Perforate all boundary cells:
        # for nth_perf in range(len(self.right_boundary_cells)):
        #     well_index = self.mesh.volume[self.right_boundary_cells[nth_perf]] / self.max_well_vol * self.well_index
        #     self.add_perforation(self.wells[-1], res_block=self.right_boundary_cells[nth_perf],
        #                          well_index=well_index)
        #
        # # Add wells to the DARTS mesh object and sort connection (DARTS related):
        self.mesh.add_wells_mpfa(ms_well_vector(self.wells), self.P_VAR)
        self.mesh.reverse_and_sort_pm()
        #self.mesh.init_grav_coef()
        return 0
    def get_normal_to_bound_face(self, b_id):
        cell = self.unstr_discr.bound_cell_info_dict[b_id]
        cells = [self.unstr_discr.mat_cells_to_node[pt] for pt in cell.nodes_to_cell]
        cell_id = next(iter(set(cells[0]).intersection(*cells)))
        for face in self.unstr_discr.faces[cell_id].values():
            if face.cell_id1 == face.cell_id2 and face.face_id2 == b_id:
                t_face = cell.centroid - self.unstr_discr.mat_cell_info_dict[cell_id].centroid
                n = face.n
                if np.inner(t_face, n) < 0: n = -n
                return n
    def get_parametrized_fault_props(self):
        ref_id = next(iter(self.unstr_discr.frac_cell_info_dict))
        tags = np.array([cell.prop_id for cell in self.unstr_discr.frac_cell_info_dict.values()])
        tag_ids = {}
        t0 = {}
        coords = {}
        z_coords = {}
        for tag in self.unstr_discr.physical_tags['fracture']:
            ids = np.argwhere(tags == tag)[:,0]
            tag_ids[tag] = ref_id + ids
            n0 = self.unstr_discr.faces[ref_id + ids[0]][4].n[:2]
            t0[tag] = np.identity(2) - np.outer(n0, n0)
            coords[tag] = np.array([self.unstr_discr.frac_cell_info_dict[i].centroid for i in tag_ids[tag]])
            z_coords[tag] = np.unique(coords[tag][:,2])
        def dist_sort_key(id):
            c_ref = self.unstr_discr.frac_cell_info_dict[list(self.unstr_discr.frac_cell_info_dict.keys())[0]].centroid
            c = self.unstr_discr.frac_cell_info_dict[id].centroid
            return (c[0] - c_ref[0]) ** 2 + (c[1] - c_ref[1]) ** 2
        def eval_frac_proj(tag, coords):
            c_ref = self.unstr_discr.frac_cell_info_dict[list(self.unstr_discr.frac_cell_info_dict.keys())[0]].centroid
            coords1 = np.copy(coords)
            coords1[0] -= c_ref[0]
            coords1[1] -= c_ref[1]
            return np.linalg.norm(t0[tag].dot(coords1), axis=0)

        output_layers = 1
        output_var_num = {tag: int(output_layers * inds.size / z_coords[tag].size) for tag, inds in tag_ids.items()}
        faults_num = len(self.unstr_discr.physical_tags['fracture'])

        s = {tag: np.zeros(num) for tag, num in output_var_num.items()}
        #gap = np.zeros( (output_var_num, 3) )
        #Ftan = np.zeros( (output_var_num, 3) )
        #Fnorm = np.zeros( output_var_num )
        inds = {tag: np.zeros((output_layers, num), dtype=np.int64) for tag, num in output_var_num.items()}
        s_ref_prev = 0
        for tag, ids in tag_ids.items():
            counter = 0
            for l, z in enumerate(z_coords[tag][:output_layers]):
                z_inds = list(ids[np.argwhere( np.logical_and(coords[tag][:,2] > z-1.E-5, coords[tag][:,2] < z+1.E-5) )[:,0]])
                z_inds.sort(key=dist_sort_key)
                pts = self.unstr_discr.frac_cell_info_dict[z_inds[0]].coord_nodes_to_cell
                s_ref = np.min(eval_frac_proj(tag, pts[:,:2].T))
                inds[tag][l] = np.array(z_inds) - ref_id
                for id in z_inds:
                    c = self.unstr_discr.frac_cell_info_dict[id].centroid[:2]
                    s[tag][counter] = eval_frac_proj(tag, c) - s_ref + s_ref_prev
                    #gap[counter] = g[id - ref_id]
                    #Ftan[counter] = Ft[id - ref_id]
                    #Fnorm[counter] = Fn[id - ref_id]
                    counter += 1
                pts = self.unstr_discr.frac_cell_info_dict[z_inds[-1]].coord_nodes_to_cell
                s_ref_prev += np.max(eval_frac_proj(tag, pts[:,:2].T)) - s_ref
        z_output = {tag: z[:output_layers] for tag, z in z_coords.items()}
        return s, z_output, inds#gap, Ftan, Fnorm
    def write_data_field(self, filename, u, s = None):
        r = np.array([cell.centroid for cell in self.unstr_discr.mat_cell_info_dict.values()])
        inds = list(np.arange(len(r)))
        inds.sort(key=lambda id: r[id][0] + 1000 * r[id][1] + 1.E+6 * r[id][2])
        if s is self.write_data_field.__defaults__[0]:
            np.savetxt(filename, np.c_[r[inds,0], r[inds,1], r[inds, 2], u[inds,0], u[inds,1], u[inds,2]])
        else:
            np.savetxt(filename, np.c_[r[inds, 0], r[inds, 1], r[inds, 2], u[inds, 0], u[inds, 1], u[inds,2],
                s[inds, 0], s[inds, 1], s[inds, 2], s[inds, 3], s[inds, 4], s[inds, 5]])
    def check_positive_negative_sides(self):
        block_size = 4
        cell_m = np.array(self.pm.cell_m,dtype=np.intp)
        cell_p = np.array(self.pm.cell_p,dtype=np.intp)
        for i, cell_id1 in enumerate(cell_m):
            cell_id2 = cell_p[i]
            if cell_id2 < self.unstr_discr.mat_cells_tot:
                # find other one
                st1 = np.array(self.pm.stencil[self.pm.offset[i]:self.pm.offset[i + 1]], dtype=np.intp)
                all_trans1 = np.array(self.pm.tran)[(self.pm.offset[i] * block_size) * block_size: (self.pm.offset[i + 1] * block_size) * block_size].reshape(self.pm.offset[i + 1] - self.pm.offset[i], block_size, block_size)
                j = np.where(np.logical_and(cell_m == cell_id2, cell_p == cell_id1))[0][0]
                st2 = np.array(self.pm.stencil[self.pm.offset[j]:self.pm.offset[j + 1]], dtype=np.intp)
                all_trans2 = np.array(self.pm.tran)[(self.pm.offset[j] * block_size) * block_size: (self.pm.offset[j + 1] * block_size) * block_size].reshape(self.pm.offset[j + 1] - self.pm.offset[j], block_size, block_size)
                assert(set(st1) == set(st2) or cell_id1 >= self.unstr_discr.mat_cells_tot or cell_id2 >= self.unstr_discr.mat_cells_tot)
                ids1 = np.argsort(st1)
                ids2 = np.argsort(st2)
                diff = all_trans1[ids1] + all_trans2[ids2]
                assert((np.abs(diff) < 1.E-6).all())
    def write_pm_conn_to_file(self, t_step, path='pm_conn.dat'):
        #self.check_positive_negative_sides()
        path = 'pm_conn' + str(t_step) + '.dat'
        block_size = 4
        f = open(path, 'w')
        f.write(str(len(self.pm.cell_m)) + '\n')
        for i, cell_id1 in enumerate(self.pm.cell_m):
            cell_id2 = self.pm.cell_p[i]
            f.write(str(cell_id1) + '\t' + str(cell_id2) + '\n')

            for k in range(block_size):
                row = 'F' + str(k) + '\t{:.5e}'.format(self.pm.rhs[i * block_size + k])
                for j in range(self.pm.offset[i], self.pm.offset[i + 1]):
                    #if (self.pm.stencil[j] < self.unstr_discr.mat_cells_tot):
                    row += '\t' + str(self.pm.stencil[j]) + '\t[' + ', '.join(
                            ['{:.5e}'.format(self.pm.tran[n]) for n in range((j * block_size + k) * block_size, (j * block_size + k + 1) * block_size)]) + str(']')
                f.write(row + '\n')
            # Biot
            for k in range(block_size):
                row = 'b' + str(k) + '\t{:.5e}'.format(self.pm.rhs_biot[i * block_size + k])
                for j in range(self.pm.offset[i], self.pm.offset[i + 1]):
                    #if (self.pm.stencil[j] < self.unstr_discr.mat_cells_tot):
                    row += '\t' + str(self.pm.stencil[j]) + '\t[' + ', '.join(
                            ['{:.5e}'.format(self.pm.tran_biot[n]) for n in range((j * block_size + k) * block_size, (j * block_size + k + 1) * block_size)]) + str(']')
                f.write(row + '\n')

            #if self.pm.cell_p[i] < self.unstr_discr.mat_cells_tot:
                #st = np.array(self.pm.stencil[self.pm.offset[i]:self.pm.offset[i+1]],dtype=np.intp)
                #all_trans = np.array(self.pm.tran)[(self.pm.offset[i] * block_size) * block_size: (self.pm.offset[i + 1] * block_size) * block_size].reshape(self.pm.offset[i + 1] - self.pm.offset[i], block_size, block_size)
                #sum = np.sum(all_trans, axis=0)
                #sum_no_bound = np.sum(all_trans[st < self.unstr_discr.mat_cells_tot], axis=0)
                #assert((abs(sum[:3,:3]) < 1.E-10).all())
        f.close()

    # def write_to_vtk(self, output_directory, ith_step, physics):
    #     """
    #     Class method which writes output of unstructured grid to VTK format
    #     :param output_directory: directory of output files
    #     :param property_array: np.array containing all cell properties (N_cells x N_prop)
    #     :param cell_property: list with property names (visible in ParaView (format strings)
    #     :param ith_step: integer containing the output step
    #     :return:
    #     """
    #     # First check if output directory already exists:
    #     if not os.path.exists(output_directory):
    #         os.makedirs(output_directory)
    #
    #     # Temporarily store mesh_data in copy:
    #     Mesh = meshio.read(self.unstr_discr.mesh_file)
    #     Mesh.cell_data.pop('quad', None)
    #     Mesh.cell_data.pop('triangle', None)
    #     Mesh.cells.pop('quad', None)
    #     Mesh.cells.pop('triangle', None)
    #     # Allocate empty new cell_data dictionary:
    #     cell_property = ['u_x', 'u_y', 'u_z', 'p']
    #     props_num = len(cell_property)
    #     property_array = np.array(physics.engine.X, copy=False)
    #     cell_data = {}
    #
    #     # if ith_step != 0:
    #     fluxes = np.array(physics.engine.fluxes, copy=False)
    #     fluxes_biot = np.array(physics.engine.fluxes_biot, copy=False)
    #     # vels = self.reconstruct_velocities(fluxes[physics.engine.P_VAR::physics.engine.N_VARS],
    #     #                                   fluxes_biot[physics.engine.P_VAR::physics.engine.N_VARS])
    #     self.mech_operators.eval_porosities(physics.engine.X, self.mesh.bc)
    #     self.mech_operators.eval_stresses(physics.engine.fluxes, physics.engine.fluxes_biot, physics.engine.X,
    #                                       self.mesh.bc, physics.engine.op_vals_arr)
    #     # else:
    #     #    self.mech_operators.eval_porosities(physics.engine.X, self.mesh.bc_prev)
    #     #    self.mech_operators.eval_stresses(physics.engine.X, self.mesh.bc_prev, physics.engine.op_vals_arr)
    #
    #     for ith_geometry in self.unstr_discr.mesh_data.cells:
    #         if ith_geometry == 'hexahedron' or ith_geometry == 'wedge' or ith_geometry == 'tetra':
    #             # Add matrix data to dictionary:
    #             cell_data[ith_geometry] = {}
    #             for i in range(props_num):
    #                 cell_data[ith_geometry][cell_property[i]] = property_array[
    #                                                             i:props_num * self.unstr_discr.mat_cells_tot:props_num]
    #
    #             # cell_data[ith_geometry]['velocity'] = vels
    #             # if hasattr(self.unstr_discr, 'E') and hasattr(self.unstr_discr, 'nu'):
    #             #     cell_data[ith_geometry]['E'] = np.zeros(self.unstr_discr.mat_cells_tot, dtype=np.float64)
    #             #     cell_data[ith_geometry]['nu'] = np.zeros(self.unstr_discr.mat_cells_tot, dtype=np.float64)
    #             #     for id, cell in enumerate(self.unstr_discr.mat_cell_info_dict.values()):
    #             #         cell_data[ith_geometry]['E'][id] = self.unstr_discr.E[cell.prop_id]
    #             #         cell_data[ith_geometry]['nu'][id] = self.unstr_discr.nu[cell.prop_id]
    #
    #             cell_data[ith_geometry]['eps_vol'] = np.array(self.mech_operators.eps_vol, copy=False)
    #             cell_data[ith_geometry]['porosity'] = np.array(self.mech_operators.porosities, copy=False)
    #             cell_data[ith_geometry]['stress'] = np.zeros((self.unstr_discr.mat_cells_tot, 6), dtype=np.float64)
    #             cell_data[ith_geometry]['tot_stress'] = np.zeros((self.unstr_discr.mat_cells_tot, 6), dtype=np.float64)
    #
    #             stress = np.array(self.mech_operators.stresses, copy=False)
    #             total_stress = np.array(self.mech_operators.total_stresses, copy=False)
    #             for i in range(6):
    #                 cell_data[ith_geometry]['stress'][:, i] = stress[i::6]
    #                 cell_data[ith_geometry]['tot_stress'][:, i] = total_stress[i::6]
    #
    #     if self.unstr_discr.frac_cells_tot > 0:
    #         self.write_fault_props(output_directory, property_array, ith_step, fluxes)
    #     # Store solution for each time-step:
    #     Mesh.cell_data = cell_data
    #     print('Writing data to VTK file for {:d}-th reporting step'.format(ith_step))
    #     meshio.write("{:s}/solution{:d}.vtk".format(output_directory, ith_step), Mesh)
    #     return 0
    def write_to_vtk(self, output_directory, ith_step, physics):
        """
        Class method which writes output of unstructured grid to VTK format
        :param output_directory: directory of output files
        :param property_array: np.array containing all cell properties (N_cells x N_prop)
        :param cell_property: list with property names (visible in ParaView (format strings)
        :param ith_step: integer containing the output step
        :return:
        """
        # First check if output directory already exists:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        # Temporarily store mesh_data in copy:
        Mesh = meshio.read(self.unstr_discr.mesh_file)

        # Allocate empty new cell_data dictionary:
        cell_property = ['u_x', 'u_y', 'u_z', 'p']
        props_num = len(cell_property)
        property_array = np.array(physics.engine.X, copy=False)
        available_matrix_geometries = ['hexahedron', 'wedge', 'tetra']
        available_fracture_geometries = ['quad', 'triangle']


        # vels = self.reconstruct_velocities(fluxes[physics.engine.P_VAR::physics.engine.N_VARS],
        #                                   fluxes_biot[physics.engine.P_VAR::physics.engine.N_VARS])
        self.mech_operators.eval_porosities(physics.engine.X, self.mesh.bc)
        self.mech_operators.eval_stresses(physics.engine.fluxes, physics.engine.fluxes_biot, physics.engine.X,
                                          self.mesh.bc, physics.engine.op_vals_arr)
        #self.mech_operators.eval_stresses(physics.engine.fluxes, physics.engine.fluxes_biot, physics.engine.X,
        #                                  self.mesh.bc, physics.engine.op_vals_arr)

        # Matrix
        Mesh.cells = []
        cell_data = {}
        start_geom_cell_id = 0
        for ith_geometry in self.unstr_discr.geom_order:
            if ith_geometry in available_matrix_geometries:
                Mesh.cells.append(next(cell for cell in self.unstr_discr.mesh_data.cells if cell.type == ith_geometry))
                # Add matrix data to dictionary:
                cell_size = self.unstr_discr.mesh_data.cells_dict[ith_geometry].shape[0]
                for i in range(props_num):
                    if cell_property[i] not in cell_data: cell_data[cell_property[i]] = []
                    cell_data[cell_property[i]].append(property_array[props_num * start_geom_cell_id + i:props_num * (cell_size + start_geom_cell_id):props_num])

                # if 'porosity' not in cell_data: cell_data['porosity'] = []
                # if 'eps_vol' not in cell_data: cell_data['eps_vol'] = []
                if 'stress' not in cell_data: cell_data['stress'] = []
                if 'tot_stress' not in cell_data: cell_data['tot_stress'] = []

                # if ith_step == 0:
                #     cell_data['eps_vol'].append(np.array(self.mech_operators.eps_vol, copy=False))
                #     cell_data['porosity'].append(np.array(self.mech_operators.porosities[start_geom_cell_id:cell_size + start_geom_cell_id], copy=False))
                # else:
                #     eps_vol = np.array(physics.engine.eps_vol, copy=False)
                #     cell_data['eps_vol'].append(eps_vol)
                #     biot = self.biot_arr[::9]
                #     comp_mult = np.where(biot != 0.0, (biot - self.poro) * (1 - biot) / self.kd, 1.0 / self.kd)
                #     eps_vol_ref = np.array(self.mesh.ref_eps_vol, copy=False)
                #     poro = self.porosity + (comp_mult * (property_array[self.P_VAR:props_num * (self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot):props_num] - self.p_ref))[:self.unstr_discr.mat_cells_tot] + (eps_vol - eps_vol_ref)
                #     poro[poro < 0.0] = self.poro[:self.unstr_discr.mat_cells_tot][poro < 0.0]
                #     cell_data['porosity'].append(poro[start_geom_cell_id:cell_size + start_geom_cell_id])
                cell_data['stress'].append(np.zeros((cell_size, 6), dtype=np.float64))
                cell_data['tot_stress'].append(np.zeros((cell_size, 6), dtype=np.float64))

                stress = np.array(self.mech_operators.stresses, copy=False)
                total_stress = np.array(self.mech_operators.total_stresses, copy=False)
                for i in range(6):
                    cell_data['stress'][-1][:, i] = stress[6 * start_geom_cell_id + i:6 * (start_geom_cell_id + cell_size):6]
                    cell_data['tot_stress'][-1][:, i] = total_stress[6 * start_geom_cell_id + i:6 * (start_geom_cell_id + cell_size):6]

                # if 'cell_id' not in cell_data: cell_data['cell_id'] = []
                # cell_data['cell_id'].append(np.array([cell_id for cell_id, cell in self.unstr_discr.mat_cell_info_dict.items() if cell.geometry_type == ith_geometry], dtype=np.intc))
                # if ith_step == 0:
                #     cell_data[ith_geometry]['permx'] = self.permx[:]
                #     cell_data[ith_geometry]['permy'] = self.permy[:]
                #     cell_data[ith_geometry]['permz'] = self.permz[:]

                start_geom_cell_id += cell_size

        # Store solution for each time-step:
        mesh = meshio.Mesh(
            Mesh.points,
            Mesh.cells,
            cell_data=cell_data)
        meshio.write("{:s}/solution{:d}.vtk".format(output_directory, ith_step), mesh)

        # Fractures
        geom_id = 0
        Mesh.cells = []
        cell_data = {}
        for ith_geometry in self.unstr_discr.mesh_data.cells_dict.keys():
            if ith_geometry in available_fracture_geometries:
                # fracture geometry
                frac_ids = np.argwhere(np.in1d(self.unstr_discr.mesh_data.cell_data['gmsh:physical'][geom_id],
                                                                                              self.unstr_discr.physical_tags['fracture']))[:, 0]
                if len(frac_ids):
                    Mesh.cells.append(meshio.CellBlock(ith_geometry, data=self.unstr_discr.mesh_data.cells[geom_id].data[frac_ids]))
                    data = property_array[4 * self.unstr_discr.mat_cells_tot:4 * (
                            self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot)].reshape(
                        self.unstr_discr.frac_cells_tot, 4)
                    for i in range(props_num):
                        if cell_property[i] not in cell_data: cell_data[cell_property[i]] = []
                        cell_data[cell_property[i]].append(data[:, i])

                # output geometry
                out_ids = np.argwhere(np.in1d(self.unstr_discr.mesh_data.cell_data['gmsh:physical'][geom_id],
                                              self.unstr_discr.physical_tags['output']))[:, 0]
                if len(out_ids):
                    Mesh.cells.append(meshio.CellBlock(type=ith_geometry, data=self.unstr_discr.mesh_data.cells[geom_id].data[out_ids]))

            geom_id += 1
        # fracture output
        if self.unstr_discr.frac_cells_tot > 0:
            self.write_fault_props(output_directory, property_array, ith_step, physics)
            frac_data = self.get_fault_props(property_array, ith_step, physics)
            for key, val in frac_data.items():
               if key not in cell_data: cell_data[key] = []
               cell_data[key].append(val)
        # just output
        if self.unstr_discr.output_face_tot > 0:
            out_data = self.get_props_over_output(property_array, ith_step, physics)
            for key, val in out_data.items():
               if key not in cell_data: cell_data[key] = []
               cell_data[key].append(val)

        # Store solution for each time-step:
        mesh = meshio.Mesh(
            Mesh.points,
            Mesh.cells,
            cell_data=cell_data)
        meshio.write("{:s}/solution_fault{:d}.vtk".format(output_directory, ith_step), mesh)

        print('Writing data to VTK file for {:d}-th reporting step'.format(ith_step))
        return 0

    def get_props_over_output(self, property_array, ith_step, physics):
        n_vars = 4
        n_dim = 3
        fluxes = np.array(physics.engine.fluxes, copy=False)
        fluxes_biot = np.array(physics.engine.fluxes_biot, copy=False)
        cell_m = np.array(self.mesh.block_m, copy=False)
        cell_p = np.array(self.mesh.block_p, copy=False)

        #S_eng = vector_matrix(physics.engine.contacts[0].S)
        #frac_prop = property_array[n_vars * self.unstr_discr.mat_cells_tot:n_vars * (self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot)].reshape(self.unstr_discr.frac_cells_tot, n_vars)
        #fstress = np.array(physics.engine.contacts[0].fault_stress, copy=False)

        S = np.zeros((n_dim, n_dim))
        frac_data = {}
        #frac_data['tag'] = np.zeros(self.unstr_discr.frac_cells_tot, dtype=np.intp)
        frac_data['f_local'] = np.zeros((self.unstr_discr.output_face_tot, n_dim))

        ref_id = self.unstr_discr.output_face_to_face[0]
        ref_face = self.unstr_discr.faces[ref_id[0]][ref_id[1]]
        S[1:n_dim] = null_space(np.array([ref_face.n])).T
        S[0] = ref_face.n
        for face_id, ids in self.unstr_discr.output_face_to_face.items():
            #cell_id -= self.unstr_discr.mat_cells_tot
            #frac_data['tag'][cell_id] = int(cell.prop_id)
            face = self.unstr_discr.faces[ids[0]][ids[1]]
            sign = np.sign((self.unstr_discr.mat_cell_info_dict[face.cell_id2].centroid - \
                            self.unstr_discr.mat_cell_info_dict[face.cell_id1].centroid).dot(ref_face.n))
            flux_ids = np.argwhere(np.logical_and(cell_m == face.cell_id1, cell_p == face.cell_id2))[0]
            flux = sign * fluxes[n_vars * flux_ids[0]:n_vars * flux_ids[0] + n_dim] / face.area

            if len(flux_ids):
                frac_data['f_local'][face_id] = S.dot(flux)
            else:
                return 0

        for face_id in range(self.unstr_discr.output_face_tot):
            print(str(face_id + self.unstr_discr.mat_cells_tot) + ' ' + str(frac_data['f_local'][face_id][0] * 1.E+5) + ' ' +
                  str(frac_data['f_local'][face_id][1] * 1.E+5) + ' ' + str(frac_data['f_local'][face_id][2] * 1.E+5))

        return frac_data
    def get_fault_props(self, property_array, ith_step, physics):
        n_vars = 4
        n_dim = 3
        fluxes = np.array(physics.engine.fluxes, copy=False)
        fluxes_biot = np.array(physics.engine.fluxes_biot, copy=False)
        S_eng = vector_matrix(physics.engine.contacts[0].S)
        frac_prop = property_array[n_vars * self.unstr_discr.mat_cells_tot:n_vars * (self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot)].reshape(self.unstr_discr.frac_cells_tot, n_vars)
        fstress = np.array(physics.engine.contacts[0].fault_stress, copy=False)

        frac_data = {}
        #frac_data['tag'] = np.zeros(self.unstr_discr.frac_cells_tot, dtype=np.intp)
        frac_data['g_local'] = np.zeros((self.unstr_discr.frac_cells_tot, n_dim))
        frac_data['f_local'] = np.zeros((self.unstr_discr.frac_cells_tot, n_dim))
        frac_data['mu'] = np.array(physics.engine.contacts[0].mu, copy=False)

        for cell_id, cell in self.unstr_discr.frac_cell_info_dict.items():
            cell_id -= self.unstr_discr.mat_cells_tot
            #frac_data['tag'][cell_id] = int(cell.prop_id)
            face = self.unstr_discr.faces[cell_id + self.unstr_discr.mat_cells_tot][4]
            S = np.array(S_eng[cell_id].values).reshape((n_dim, n_dim))
            f = fstress[n_dim * cell_id:n_dim * (cell_id + 1)] / face.area
            frac_data['f_local'][cell_id] = S.dot(f)
            frac_data['g_local'][cell_id] = S.dot(frac_prop[cell_id,:n_dim])

        phi = np.array(physics.engine.contacts[0].phi, copy=False)
        #states = phi > 0
        frac_data['phi'] = phi

        return frac_data
    def write_fault_props(self, output_directory, property_array, ith_step, physics):
        n_vars = 4
        n_dim = 3
        fluxes = np.array(physics.engine.fluxes, copy=False)
        fluxes_biot = np.array(physics.engine.fluxes_biot, copy=False)
        s, z_coords, inds = self.get_parametrized_fault_props()
        x = property_array.reshape(int(property_array.size / n_vars), n_vars)[self.unstr_discr.mat_cells_tot:self.unstr_discr.mat_cells_tot + self.unstr_discr.frac_cells_tot]
        g = {}
        glocal = {}
        flocal = {}
        mu = {}
        for tag, ids in inds.items():
            if ids.size * z_coords[tag].size < self.unstr_discr.frac_cells_tot: return 0

            g[tag] = np.array(x[ids[0],:n_dim])
            glocal[tag] = np.zeros((len(ids[0]), 3))
            flocal[tag] = np.zeros((len(ids[0]), 3))
            S_eng = vector_matrix(physics.engine.contacts[0].S)
            mu[tag] = np.array(physics.engine.contacts[0].mu, copy=False)
            fstress = np.array(physics.engine.contacts[0].fault_stress, copy=False)
            for i, id in enumerate(ids[0]):
                face = self.unstr_discr.faces[id + self.unstr_discr.mat_cells_tot][4]
                S = np.array(S_eng[id].values).reshape((n_dim, n_dim))
                flocal[tag][i] = S.dot(fstress[n_dim * id:n_dim * (id + 1)] / face.area)
                #n = self.unstr_discr.faces[self.unstr_discr.mat_cells_tot][max(self.unstr_discr.faces[self.unstr_discr.mat_cells_tot].keys())].n
                #S = np.zeros((n_dim, n_dim))
                #S[:n_dim - 1] = null_space(np.array([-n])).T
                #S[n_dim - 1] = -n
                glocal[tag][i] = S.dot(g[tag][i])

            #if ith_step == 0:
            self.fig, self.ax = plt.subplots(nrows=2, sharex=True, figsize=(12, 10))
            self.ax0 = self.ax[0].twinx()
            self.ax1 = self.ax[1].twinx()
            self.ax11 = self.ax[1].twinx()
            #self.ax[0].set_ylabel('normal gap, $g_N$')
            self.ax[0].set_ylabel('friction coefficient, $\mu$')
            self.ax0.set_ylabel('slip, $g_T$')
            self.ax[1].set_ylabel('normal traction, $F_N$')
            self.ax1.set_ylabel('tangential traction, $F_T$')
            self.ax[1].set_xlabel('distance')
                #self.ax1.set_ylabel('distance along fault')

            phi = np.array(physics.engine.contacts[0].phi, copy=False)
            states = phi > 0
            for tag, s_cur in s.items():
                Fn = flocal[tag][:,0]
                Ft = flocal[tag][:,1]
                #self.ax[0].plot(s_cur, glocal[tag][:,0], color='b', linestyle='-', marker='o', label=str(tag) + r': $g_N$')
                if (mu[tag] != 0.0).all() and (Fn != 0.0).all():
                    self.ax[0].plot(s_cur, mu[tag], color='b', linestyle='-', marker='o', label=str(tag) + r': $\mu$')
                    self.ax[0].plot(s_cur, Ft / Fn, color='r', linestyle='--', marker='o', label=str(tag) + r': $\mu * SCU$')
                self.ax0.plot(s_cur, -glocal[tag][:,1], color='r', linestyle='-', marker='o', label=str(tag) + r': $g_T$')
                self.ax[1].plot(s_cur, Fn, color='b', linestyle='-', marker='o', label=str(tag) + r': $F_N$')
                self.ax1.plot(s_cur, Ft, color='r', linestyle='-', marker='o', label=str(tag) + r': $F_T$')
                self.ax11.plot(s_cur, states[ids[0]], color='g', linestyle=':', marker='x')
                if states[ids[0]][0] == 0:
                    self.ax11.text(0, 0, 'STUCK', fontsize=15)
                elif states[ids[0]][0] == 1:
                    self.ax11.text(0, 1, 'SLIP', fontsize=15)

                np.savetxt(output_directory + '/fault_step_' + str(ith_step) + '_tag_' + str(tag) + ".txt",
                           np.c_[s_cur, glocal[tag][:, 0], glocal[tag][:, 1], glocal[tag][:, 2], Fn, Ft, mu[tag]])

            self.ax[0].grid(axis='x')
            self.ax[1].grid(axis='x')
            self.ax0.grid(axis='y')
            self.ax1.grid(axis='y')
            self.ax[0].legend(loc='upper left')
            self.ax0.legend(loc='upper right')
            self.ax[1].legend(loc='upper left')
            self.ax1.legend(loc='upper right')
            if mu[tag].min() != mu[tag].max():
                self.ax[0].set_ylim([0.98 * mu[tag].min(), 1.02 * mu[tag].max()])
            else:
                self.ax[0].set_ylim([0.0, 1.0])
            #self.ax0.set_ylim([-0.0002, 0.0006])
            #self.ax[1].set_ylim([19, 21])
            #self.ax1.set_ylim([-0.1, 0.1])
            self.ax11.get_yaxis().set_visible(False)

            self.fig.tight_layout()
            self.fig.savefig(output_directory + '/fig_' + str(ith_step) + '.png')
            plt.close(self.fig)