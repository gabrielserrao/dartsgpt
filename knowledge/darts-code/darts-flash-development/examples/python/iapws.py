import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from dartsflash.libflash import EoS, IAPWS95
from dartsflash.libflash import PXFlash

from dartsflash.pyflash import PyFlash
from dartsflash.mixtures import Mixture
from dartsflash.plot import *


vmax = None
components = ["H2O"]
z = [1.]
n = z
pressure = np.arange(1., 251.1, 25)
# pressure = np.array([60.])
temperature = np.array([270., 640., 647.2, 800.])
# temperature = np.array([203.15])
temperaturePT = np.arange(403.15, 434.15, 10)
vmax = 6e-3
nvol = 1001
prange = [0., 1000]

mix = Mixture(components)
f = PyFlash(mixture=mix)

f.add_eos("IAPWS95", IAPWS95(mix.comp_data, iapws_ideal=True))
iapws = f.eos["IAPWS95"]

if 0:
    """ Compressibility factor and P-V diagram of Cubic EoS """
    """ Evaluate properties at PT """
    state_spec = {"pressure": pressure,
                  "temperature": temperature
                  }
    compositions = {comp: n[i] for i, comp in enumerate(components)}
    properties = {"V": iapws.V,
                  "Z": iapws.Z,
                  }

    # pprops = f.evaluate_phase_properties_1p(state_spec=state_spec, compositions=compositions,
    #                                   properties_to_evaluate=properties, mole_fractions=True)
    pprops = None

    """ Evaluate properties at VT """
    vmin = iapws.V(p=pressure[-1], T=temperature[0], n=n)
    vmin = 3e-5
    vmax = iapws.V(p=pressure[0], T=temperature[-1], n=n) if vmax is None else vmax

    state_spec = {"temperature": temperature,
                  # "volume": np.linspace(vmin, vmax, nvol)
                  "volume": np.linspace(1e-8, 3., 100)
                  }
    properties = {"Z": iapws.Zd,
                  "P": iapws.Pd
                  }

    vprops = f.evaluate_phase_properties_1p(state_spec=state_spec, compositions=compositions,
                                      properties_to_evaluate=properties, mole_fractions=True)
    # vprops = None

    """ Plot PV and PZ diagrams """
    pv = PlotEoS.pressure_volume(f, temperatures=temperature, compositions=n,
                                 p_props=pprops, v_props=vprops,
                                 v_range=[0, 2],
                                 # p_range=prange,
                                 p_range=[8e1, 1e2],
                                 logy=True,
                                 )

    pz = PlotEoS.compressibility(f, temperatures=temperature, compositions=n,
                                 p_props=None, v_props=vprops, z_range=[-0.1, 1.1], p_range=prange)
    plt.savefig("P-Z-" + mix.filename + ".pdf")

if 1:
    """ PT-diagram of properties """
    """ Evaluate properties at PT """
    state_spec = {"pressure": np.logspace(-2, 2.5, 20),
                  "temperature": np.arange(270, 700, 1e0)
                  }
    compositions = {comp: n[i] for i, comp in enumerate(components)}
    properties = {#"Volume": iapws.V,
                  # "Density": iapws.rho,
                  # "VolumeIterations": iapws.volume_iterations,
                  "Root": iapws.is_root_type,
                  # "CriticalT": iapws.is_critical,
                  }

    # pprops = f.evaluate_phase_properties_1p(state_spec=state_spec, compositions=compositions,
    #                                         properties_to_evaluate=properties, mole_fractions=True)
    props = {
        # "H": EoS.Property.ENTHALPY,
        "S": EoS.Property.ENTROPY,
        # "G": EoS.Property.GIBBS
    }
    h = f.evaluate_properties_1p(state_spec=state_spec, compositions=compositions, mole_fractions=True,
                                 properties_to_evaluate=props, )
    # cp = iapws.critical_point(n)

    # plot = PlotEoS.surf(f, props=pprops, x_var='temperature', y_var='pressure', prop_names=list(properties.keys()),
    #                     composition=n, logy=True)
    plot2 = PlotEoS.surf(f, props=h, x_var='temperature', y_var='pressure', prop_names=list(props.keys()),
                         composition=n, logy=True)
    # for i in range(len(properties.keys())):
    #     plot.subplot_idx = i
    #     plot.draw_point(X=cp.Tc, Y=cp.Pc, colours='red')
    # plt.savefig("IAPWS-enthalpy.pdf")

if 0:
    f.flash_params.eos_order = ["IAPWS95"]
    f.flash_params.eos_params["IAPWS95"].root_order = [EoS.MAX, EoS.MIN]
    f.flash_params.T_min = 250.
    f.flash_params.T_max = 700.
    f.flash_params.phflash_Htol = 1e-3
    f.flash_params.phflash_Ttol = 1e-8
    # f.flash_params.split_switch_tol = 1e-4

    f.np_max = 2

    Trange = [270., 700.]
    Prange = [1., 250.]
    # Hrange = [-40000., 10000.]
    logy = False

    if 0:
        Xspec = PXFlash.ENTHALPY
        Xrange = [f.eos["IAPWS95"].H(Prange[1], Trange[i], n, 0, pt=True) * R for i in range(2)]
    else:
        Xspec = PXFlash.ENTROPY
        Xrange = [f.eos["IAPWS95"].S(Prange[1], Trange[i], n, 0, pt=True) * R for i in range(2)]

    f.flash_params.verbose = 0
    f.f = PXFlash(f.flash_params, Xspec)

    state_spec = {"pressure": np.linspace(Prange[0], Prange[1], 100) if not logy
    else np.logspace(np.log10(Prange[0]), np.log10(Prange[1]), 100),
                  "enthalpy" if Xspec == PXFlash.ENTHALPY else "entropy":
                      np.linspace(Xrange[0], Xrange[1], 100),
                  }
    compositions = {comp: np.arange(min_z[i], max_z[i] + 0.1 * dz, dz) for i, comp in enumerate(components[:-1])}
    compositions[components[-1]] = 1.

    if len(components) == 1:
        flash_results = f.evaluate_flash_1c(state_spec=state_spec, print_state="Flash")
    else:
        flash_results = f.evaluate_flash(state_spec=state_spec, compositions=compositions, mole_fractions=True,
                                         print_state="Flash")

    plot_pt = False
    if plot_pt:
        state_pt = {"pressure": state_spec["pressure"],
                    "temperature": np.linspace(Trange[0], Trange[1], 100), }
        pt_props = f.evaluate_properties_1p(state_spec=state_pt, compositions=compositions, mole_fractions=True,
                                            properties_to_evaluate={"H": f.eos["CEOS"].H} if Xspec == PXFlash.ENTHALPY
                                            else {"S": f.eos["CEOS"].S}
                                            )
    else:
        pt_props = None

    plot_method = PlotFlash.ph if Xspec == PXFlash.ENTHALPY else PlotFlash.ps
    plot = plot_method(f, flash_results, composition=z, min_temp=Trange[0], max_temp=Trange[1], min_val=0., max_val=1.,
                       plot_phase_fractions=True, pt_props=pt_props)

    plt.savefig(mix.filename + "-" + "-".join(str(int(zi * 100)) for zi in z) + "-ph.pdf")

plt.show()
