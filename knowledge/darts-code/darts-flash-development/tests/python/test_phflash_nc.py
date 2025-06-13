import numpy as np
import xarray as xr
import pytest

from dartsflash.libflash import PXFlash, EoS, CubicEoS, AQEoS, InitialGuess
from dartsflash.pyflash import PyFlash
from .conftest import Mixture


# Test of PH diagram for CO2 vapour and liquid
@pytest.fixture()
def flash_obj_co2c1():
    mix = Mixture(components=["CO2", "C1"], setprops=True)
    f = PyFlash(mixture=mix)

    f.add_eos("CEOS", CubicEoS(mix.comp_data, CubicEoS.PR),
              initial_guesses=[InitialGuess.Wilson, InitialGuess.Wilson13],
              switch_tol=1e-3)

    f.flash_params.phflash_Htol = 1e-3
    f.flash_params.phflash_Ttol = 1e-8
    # f.flash_params.split_switch_tol = 1e-5

    f.init_flash(eos_order=["CEOS"], flash_type=PyFlash.FlashType.PHFlash, np_max=2)

    return f


@pytest.fixture()
def ref_phflash_co2c1():
    return "./tests/python/data/ref_phflash_co2c1.h5"



def test_phflash_co2c1(flash_obj_co2c1, ref_phflash_co2c1):
    z = [0.9, 0.1]
    Trange = [260., 350.]
    Hrange = [flash_obj_co2c1.eos["CEOS"].H(100., Ti, z, 0, pt=True) * 8.314472 for Ti in Trange]
    state_spec = {"pressure": np.linspace(10., 100., 20),
                  "enthalpy": np.linspace(Hrange[0], Hrange[1], 20)
                  }
    compositions = {comp: z[i] for i, comp in enumerate(flash_obj_co2c1.components)}

    results = flash_obj_co2c1.evaluate_flash(state_spec=state_spec, compositions=compositions, mole_fractions=True)
    results = results[["pressure", "enthalpy"]]

    if 1:
        # Use assertions from xarray
        ref = xr.open_dataset(ref_phflash_co2c1, engine='h5netcdf')
        ref = ref[["pressure", "enthalpy"]]

        xr.testing.assert_allclose(results, ref, rtol=1e-5)

    else:
        # Update reference file
        results.to_netcdf(ref_phflash_co2c1, engine='h5netcdf')
