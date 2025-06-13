import numpy as np
import xarray as xr
import pytest

from dartsflash.libflash import InitialGuess, CubicEoS, AQEoS
from dartsflash.pyflash import PyFlash
from .conftest import mixture_brine_co2


# Test of full compositional space for H2O-CO2 mixture
@pytest.fixture()
def flash_obj(mixture_brine_co2):
    mix = mixture_brine_co2
    f = PyFlash(mixture=mix)

    f.add_eos("PR", CubicEoS(mix.comp_data, CubicEoS.PR),
              initial_guesses=[InitialGuess.Yi.Wilson])
    f.add_eos("AQ", AQEoS(mix.comp_data, {AQEoS.water: AQEoS.Jager2003,
                                          AQEoS.solute: AQEoS.Ziabakhsh2012}), initial_guesses=[int(mix.comp_data.components.index("H2O"))])
    f.init_flash(flash_type=PyFlash.FlashType.NegativeFlash, eos_order=["AQ", "PR"],
                 initial_guess=[InitialGuess.Henry_AV])
    return f


@pytest.fixture()
def ref_flash_brine_co2():
    return "./tests/python/data/ref_flash_brine_co2.h5"


# def test_stability_flash_brine_co2(flash_obj, ref_flash_brine_co2):
#     flash_obj.init_flash(stabilityflash=True)
#
#     state_spec = ["pressure", "temperature"]
#     dimensions = {"pressure": np.array([1., 2., 5., 10., 50., 100., 200., 400.]),
#                   "temperature": np.arange(273.15, 373.15, 10),
#                   "H2O": np.concatenate([np.array([1e-14, 1e-12, 1e-10, 1e-8]), np.linspace(1e-6, 1. - 1e-6, 10),
#                                          np.array([1. - 1e-8, 1. - 1e-10, 1. - 1e-12, 1. - 1e-14])]),
#                   }
#     constants = {"CO2": 1.}
#
#     results = flash_obj.evaluate_flash(state_spec=state_spec, dimensions=dimensions, constants=constants,
#                                        mole_fractions=True)
#     results = results[["pressure", "temperature", "H2O"]]
#
#     if 1:
#         # Use assertions from xarray
#         ref = xr.open_dataset(ref_flash_brine_co2, engine='h5netcdf')
#         ref = ref[["pressure", "temperature", "H2O"]]
#
#         xr.testing.assert_allclose(results, ref, rtol=1e-5)
#     else:
#         # Update reference file
#         results.to_netcdf(ref_flash_brine_co2, engine='h5netcdf')


def test_negative_flash_brine_co2(flash_obj, ref_flash_brine_co2):
    state_spec = {"pressure": np.array([1., 2., 5., 10., 50., 100., 200., 400.]),
                  "temperature": np.arange(273.15, 373.15, 10),
                  }
    compositions = {"H2O": np.concatenate([np.array([1e-14, 1e-12, 1e-10, 1e-8]), np.linspace(1e-6, 1. - 1e-6, 10),
                                           np.array([1. - 1e-8, 1. - 1e-10, 1. - 1e-12, 1.-1e-14])]),
                    "CO2": 1.,
                    }

    results = flash_obj.evaluate_flash(state_spec=state_spec, compositions=compositions, mole_fractions=True)
    results = results[["pressure", "temperature", "H2O"]]

    if 1:
        # Use assertions from xarray
        ref = xr.open_dataset(ref_flash_brine_co2, engine='h5netcdf')
        ref = ref[["pressure", "temperature", "H2O"]]

        xr.testing.assert_allclose(results, ref, rtol=1e-5)
    else:
        # Update reference file
        results.to_netcdf(ref_flash_brine_co2, engine='h5netcdf')
