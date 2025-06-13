import numpy as np
import xarray as xr
import pytest

from dartsflash.libflash import InitialGuess, CubicEoS, AQEoS
from dartsflash.pyflash import PyFlash
from .conftest import mixture_brine_vapour


# Test of full compositional space for H2O-CO2-C1 mixture
@pytest.fixture()
def flash_obj(mixture_brine_vapour):
    mix = mixture_brine_vapour
    f = PyFlash(mixture=mix)

    f.add_eos("PR", CubicEoS(mix.comp_data, CubicEoS.PR),
              initial_guesses=[InitialGuess.Yi.Wilson])
    f.add_eos("AQ", AQEoS(mix.comp_data, {AQEoS.water: AQEoS.Jager2003,
                                          AQEoS.solute: AQEoS.Ziabakhsh2012}))
    f.flash_params.eos_order = ["AQ", "PR"]
    f.init_flash(flash_type=PyFlash.FlashType.NegativeFlash, eos_order=["AQ", "PR"],
                 initial_guess=[InitialGuess.Henry_AV])
    return f


@pytest.fixture()
def ref_brine_vapour2():
    return "./tests/python/data/ref_flash_brine_vapour2.h5"


def test_brine_vapour2(flash_obj, ref_brine_vapour2):
    state_spec = {"pressure": np.array([1., 2., 5., 10., 50., 100., 200., 400.]),
                  "temperature": np.arange(273.15, 373.15, 10),
                  }
    zrange = np.concatenate([np.array([1e-12, 1e-10, 1e-8]),
                             np.linspace(1e-6, 1. - 1e-6, 10),
                             np.array([1. - 1e-8, 1. - 1e-10, 1. - 1e-12])])
    compositions = {"H2O": zrange,
                    "CO2": zrange,
                    "C1": 1.,
                    }

    results = flash_obj.evaluate_flash(state_spec=state_spec, compositions=compositions, mole_fractions=True)
    results = results[["pressure", "temperature", "H2O", "CO2"]]

    if 1:
        # Use assertions from xarray
        ref = xr.open_dataset(ref_brine_vapour2, engine='h5netcdf')
        ref = ref[["pressure", "temperature", "H2O", "CO2"]]

        xr.testing.assert_allclose(results, ref, rtol=1e-5)
    else:
        # Update reference file
        results.to_netcdf(ref_brine_vapour2, engine='h5netcdf')
