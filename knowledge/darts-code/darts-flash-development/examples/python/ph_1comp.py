import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from dartsflash.libflash import FlashParams, EoS, InitialGuess
from dartsflash.libflash import PXFlash
from dartsflash.libflash import CubicEoS, AQEoS

from dartsflash.mixtures import Mixture
from dartsflash.pyflash import PyFlash, R
from dartsflash.plot import *


if 1:
    if 0:
        components = ["CO2"]
        mix = Mixture(components=components, setprops=True)
        Trange = [260., 350.]
        Prange = [10., 100.]
        np_max = 2
    elif 1:
        components = ["CO2", "C1"]
        mix = Mixture(components=components, setprops=True)

        dz = 0.01
        min_z = [0.]
        max_z = [1.]

        Trange = [260., 350.]
        Prange = [10., 100.]
        np_max = 3
    else:
        components = ["C1", "nC4"]
        mix = Mixture(components=components, setprops=True)
        mix.comp_data.Pc = [46.0, 38.0]
        mix.comp_data.Tc = [190.60, 425.20]
        mix.comp_data.ac = [0.008, 0.193]
        mix.comp_data.kij = np.zeros(4)
        mix.comp_data.T0 = 273.15
        z = [0.99, 0.01]
        Trange = [100., 300.]
        Prange = [10., 100.]
        np_max = 3

    f = PyFlash(mixture=mix)

    f.add_eos("CEOS", CubicEoS(mix.comp_data, CubicEoS.PR),
              initial_guesses=[InitialGuess.Wilson, InitialGuess.Wilson13],
              switch_tol=1e-3)
    f.flash_params.eos_params["CEOS"].root_order = [EoS.MAX, EoS.MIN]
    f.flash_params.phflash_Htol = 1e-3
    f.flash_params.phflash_Ttol = 1e-8
    # f.flash_params.split_switch_tol = 1e-5
    f.np_max = np_max

    logy = False

elif 0:
    if 0:
        components = ["H2O"]
        z = [1.]
    else:
        components = ["H2O", "CO2"]
        z = [0.5, 0.5]
    mix = Mixture(components=components, setprops=True)

    f = PyFlash(mixture=mix)

    f.add_eos("CEOS", CubicEoS(mix.comp_data, CubicEoS.PR),
              initial_guesses=[0, 1],
              switch_tol=1e-3, preferred_roots=[(0, 0.75, EoS.MAX)])

    f.add_eos("AQ", AQEoS(mix.comp_data, {AQEoS.water: AQEoS.Jager2003,
                                          AQEoS.solute: AQEoS.Ziabakhsh2012}),
              initial_guesses=[0], eos_range={0: [0.6, 1.]}, max_iter=10, use_gmix=True)

    f.flash_params.eos_order = ["AQ", "CEOS"]
    f.flash_params.eos_params["CEOS"].root_order = [EoS.MAX, EoS.MIN]
    f.flash_params.T_min = 250.
    f.flash_params.T_max = 575.
    f.flash_params.phflash_Htol = 1e-3
    f.flash_params.phflash_Ttol = 1e-8
    # f.flash_params.split_switch_tol = 1e-4

    f.np_max = 3

    Trange = [270., 575.]
    Prange = [10., 100.]
    Hrange = [-40000., 10000.]
    logy = False
else:
    components = ["H2O"]
    mix = Mixture(components=components, setprops=True)
    z = [1.]

    f = PyFlash(mixture=mix)

    from dartsflash.libflash import IdealGas, IAPWS_Ice, PureSolid, CubicEoS
    # f.add_eos("V", IdealGas(mix.comp_data))
    f.add_eos("V", CubicEoS(mix.comp_data, CubicEoS.PR),
              preferred_roots=[(0, 0.75, EoS.MAX)])
    f.add_eos("Ice", PureSolid(mix.comp_data, "Ice"))
    # f.add_eos("Ice", IAPWS_Ice(mix.comp_data))

    f.flash_params.eos_order = ["Ice", "V"]
    f.flash_params.T_min = 100.
    f.flash_params.T_max = 300.
    f.flash_params.phflash_Htol = 1e-3
    f.flash_params.phflash_Ttol = 1e-8
    # f.flash_params.split_switch_tol = 1e-4

    f.np_max = 2

    Trange = [100., 200.]
    Prange = [1e-5, 1e-2]
    Hrange = [f.eos["Ice"].H(Prange[1], Trange[0], z) * 8.314472,
              f.eos["V"].H(Prange[1], Trange[1], z) * 8.314472]

    logy = True
    if 1:
        Xspec = PXFlash.ENTHALPY
        z = [0.5, 0.5]
        Xrange = [f.eos["CEOS"].H(Prange[1], Trange[i], z, 0, pt=True) * R for i in range(2)]
    else:
        Xspec = PXFlash.ENTROPY
        z = [0.5, 0.5]
        Xrange = [f.eos["CEOS"].S(Prange[1], Trange[i], z, 0, pt=True) * R for i in range(2)]

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
                "temperature": np.linspace(Trange[0], Trange[1], 100),}
    pt_props = f.evaluate_properties_1p(state_spec=state_pt, compositions=compositions, mole_fractions=True,
                                        properties_to_evaluate={"H": f.eos["CEOS"].H} if Xspec == PXFlash.ENTHALPY
                                                          else {"S": f.eos["CEOS"].S}
                                        )
else:
    pt_props = None

if 1:
    plot_method = PlotFlash.ph if Xspec == PXFlash.ENTHALPY else PlotFlash.ps
    plot = plot_method(f, flash_results, composition=z, min_temp=250., max_temp=350., min_val=0., max_val=1.,
                       plot_phase_fractions=True, pt_props=pt_props)

    plt.savefig(mix.filename + "-" + "-".join(str(int(zi*100)) for zi in z) + "-ph.pdf")

if 1 and len(components) > 1:
    plot = PlotFlash.binary_xz(f, flash_results, variable_comp_idx=0, dz=dz, state=state_spec,

                               )
    plt.savefig(mix.filename + "-" + "-".join(str(int(zi*100)) for zi in z) + "-hx.pdf")


# plt.tight_layout()
plt.show()
