import numpy as np
import xarray as xr
import pytest

from dartsflash.libflash import PXFlash, EoS, CubicEoS, AQEoS
from dartsflash.pyflash import PyFlash
from .conftest import Mixture


# Test of PH diagram for CO2 vapour and liquid
@pytest.fixture()
def flash_obj_co2():
    mix = Mixture(components=["CO2"], setprops=True)
    f = PyFlash(mixture=mix)

    f.add_eos("CEOS", CubicEoS(mix.comp_data, CubicEoS.PR))

    f.init_flash(eos_order=["AQ", "CEOS"], flash_type=PyFlash.FlashType.PHFlash, np_max=2)

    return f


@pytest.fixture()
def ref_phflash_co2():
    return "./tests/python/data/ref_phflash_co2.h5"


# Test of PH diagram for H2O vapour and liquid
@pytest.fixture()
def flash_obj_h2o():
    mix = Mixture(components=["H2O"], setprops=True)
    f = PyFlash(mixture=mix)
    f.flash_params.T_min = 250.
    f.flash_params.T_max = 575.

    f.add_eos("CEOS", CubicEoS(mix.comp_data, CubicEoS.PR),
              preferred_roots=[(0, 0.75, EoS.MAX)])
    f.add_eos("AQ", AQEoS(mix.comp_data, {AQEoS.water: AQEoS.Jager2003,
                                          AQEoS.solute: AQEoS.Ziabakhsh2012}),
              eos_range={0: [0.6, 1.]}
              )

    f.flash_params.eos_order = ["AQ", "CEOS"]
    f.flash_params.T_min = 250.
    f.flash_params.T_max = 575.
    f.flash_params.phflash_Htol = 1e-3
    f.flash_params.phflash_Ttol = 1e-8
    f.flash_params.split_switch_tol = 1e-4
    f.flash_params.verbose = 0

    f.init_flash(eos_order=["CEOS"], flash_type=PyFlash.FlashType.PHFlash, np_max=2)

    return f


@pytest.fixture()
def ref_phflash_h2o():
    return "./tests/python/data/ref_phflash_h2o.h5"



def test_phflash_co2(flash_obj_co2, ref_phflash_co2):
    Trange = [260., 350.]
    Hrange = [flash_obj_co2.eos["CEOS"].H(100., Ti, [1.], 0, pt=True) * 8.314472 for Ti in Trange]
    state_spec = {"pressure": np.linspace(10., 100., 20),
                  "enthalpy": np.linspace(Hrange[0], Hrange[1], 20)
                  }

    results = flash_obj_co2.evaluate_flash_1c(state_spec=state_spec)
    results = results[["pressure", "enthalpy"]]

    if 1:
        # Use assertions from xarray
        ref = xr.open_dataset(ref_phflash_co2, engine='h5netcdf')
        ref = ref[["pressure", "enthalpy"]]

        xr.testing.assert_allclose(results, ref, rtol=1e-5)

    else:
        # Update reference file
        results.to_netcdf(ref_phflash_co2, engine='h5netcdf')


def test_phflash_h2o(flash_obj_h2o, ref_phflash_h2o):
    Prange = [10., 100.]
    Hrange = [-40000., 10000.]

    state_spec = {"pressure": np.linspace(Prange[0], Prange[1], 20),
                  "enthalpy": np.linspace(Hrange[0], Hrange[1], 10),
                  }

    results = flash_obj_h2o.evaluate_flash_1c(state_spec=state_spec)
    results = results[["pressure", "enthalpy"]]

    if 1:
        # Use assertions from xarray
        ref = xr.open_dataset(ref_phflash_h2o, engine='h5netcdf')
        ref = ref[["pressure", "enthalpy"]]

        xr.testing.assert_allclose(results, ref, rtol=1e-5)

    else:
        # Update reference file
        results.to_netcdf(ref_phflash_h2o, engine='h5netcdf')
