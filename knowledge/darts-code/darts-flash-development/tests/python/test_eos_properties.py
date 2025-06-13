import numpy as np
import xarray as xr
import pytest

from dartsflash.libflash import CubicEoS, AQEoS, Ballard, InitialGuess, EoS
from dartsflash.pyflash import PyFlash
from .conftest import mixture_brine_vapour


# Test of full compositional space for H2O-CO2-C1 mixture
@pytest.fixture()
def flash_obj_cubic(mixture_brine_vapour):
    mix = mixture_brine_vapour
    f = PyFlash(mixture=mix)

    f.add_eos("PR", CubicEoS(mix.comp_data, CubicEoS.PR),
              initial_guesses=[InitialGuess.Yi.Wilson])
    return f


@pytest.fixture()
def flash_obj_hybrid(mixture_brine_vapour):
    mix = mixture_brine_vapour
    f = PyFlash(mixture=mix)

    f.add_eos("PR", CubicEoS(mix.comp_data, CubicEoS.PR),
              initial_guesses=[InitialGuess.Yi.Wilson])
    f.add_eos("AQ", AQEoS(mix.comp_data, {AQEoS.water: AQEoS.Jager2003,
                                          AQEoS.solute: AQEoS.Ziabakhsh2012}))
    f.add_eos("sI", Ballard(mix.comp_data, "sI"))
    return f


@pytest.fixture()
def ref_properties_pt_cubic():
    return "./tests/python/data/ref_properties_pt_cubic.h5"


@pytest.fixture()
def ref_properties_pt_hybrid():
    return "./tests/python/data/ref_properties_pt_hybrid.h5"


@pytest.fixture()
def ref_properties_vt():
    return "./tests/python/data/ref_properties_vt.h5"


def test_properties_pt_cubic(flash_obj_cubic, ref_properties_pt_cubic):
    state_spec = {"pressure": np.array([1., 2., 5., 10., 50., 100., 200., 400.]),
                  "temperature": np.arange(273.15, 373.15, 10)
                  }
    zrange = np.concatenate([np.array([1e-12, 1e-10, 1e-8]),
                             np.linspace(1e-6, 1. - 1e-6, 10),
                             np.array([1. - 1e-8, 1. - 1e-10, 1. - 1e-12])])
    compositions = {"H2O": zrange,
                    "CO2": zrange,
                    "C1": 1.
                    }
    eos = flash_obj_cubic.eos["PR"]
    properties_cubic = {"V": eos.V,
                        "Cv": eos.Cv,
                        "Cp": eos.Cp,
                        "JT": eos.JT,
                        "vs": eos.vs,
                        }

    results = flash_obj_cubic.evaluate_phase_properties_1p(state_spec=state_spec, compositions=compositions,
                                                           properties_to_evaluate=properties_cubic, mole_fractions=True)
    results = results[["pressure", "temperature", "H2O", "CO2"]]

    if 1:
        # Use assertions from xarray
        ref = xr.open_dataset(ref_properties_pt_cubic, engine='h5netcdf')
        ref = ref[["pressure", "temperature", "H2O", "CO2"]]

        xr.testing.assert_allclose(results, ref, rtol=1e-5)

    else:
        # Update reference file
        results.to_netcdf(ref_properties_pt_cubic, engine='h5netcdf')


def test_properties_pt_hybrid(flash_obj_hybrid, ref_properties_pt_hybrid):
    state_spec = {"pressure": np.array([1., 2., 5., 10., 50., 100., 200., 400.]),
                  "temperature": np.arange(273.15, 373.15, 10)
                  }
    zrange = np.concatenate([np.array([1e-12, 1e-10, 1e-8]),
                             np.linspace(1e-6, 1. - 1e-6, 10),
                             np.array([1. - 1e-8, 1. - 1e-10, 1. - 1e-12])])
    compositions = {"H2O": zrange,
                    "CO2": zrange,
                    "C1": 1.
                    }
    properties_hybrid = {"S": EoS.Property.ENTROPY,
                         "G": EoS.Property.GIBBS,
                         "H": EoS.Property.ENTHALPY,
                         }

    results = flash_obj_hybrid.evaluate_properties_1p(state_spec=state_spec, compositions=compositions,
                                                      properties_to_evaluate=properties_hybrid,
                                                      mix_properties_to_evaluate=properties_hybrid,
                                                      mole_fractions=True)
    results = results[["pressure", "temperature", "H2O", "CO2"]]

    if 1:
        # Use assertions from xarray
        ref = xr.open_dataset(ref_properties_pt_hybrid, engine='h5netcdf')
        ref = ref[["pressure", "temperature", "H2O", "CO2"]]

        xr.testing.assert_allclose(results, ref, rtol=1e-5)

    else:
        # Update reference file
        results.to_netcdf(ref_properties_pt_hybrid, engine='h5netcdf')


def test_properties_vt(flash_obj_cubic, ref_properties_vt):
    pr = flash_obj_cubic.eos["PR"]
    vmax = pr.V(p=1., T=273.15, n=[1e-12, 1e-12, 1.-2e-12])
    vmin = pr.V(p=100., T=273.15, n=[1.-2e-12, 1e-12, 1e-12])
    state_spec = {"volume": np.linspace(vmin, vmax, 100),
                  "temperature": np.arange(273.15, 373.15, 10)
                  }
    zrange = np.concatenate([np.array([1e-12, 1e-10, 1e-8]),
                             np.linspace(1e-6, 1. - 1e-6, 10),
                             np.array([1. - 1e-8, 1. - 1e-10, 1. - 1e-12])])
    compositions = {"H2O": zrange,
                    "CO2": zrange,
                    "C1": 1.
                    }

    eos = flash_obj_cubic.eos["PR"]
    properties = {"P": eos.P}

    results = flash_obj_cubic.evaluate_phase_properties_1p(state_spec=state_spec, compositions=compositions,
                                                           properties_to_evaluate=properties, mole_fractions=True)
    results = results[["volume", "temperature", "H2O", "CO2"]]

    if 1:
        # Use assertions from xarray
        ref = xr.open_dataset(ref_properties_vt, engine='h5netcdf')
        ref = ref[["volume", "temperature", "H2O", "CO2"]]

        xr.testing.assert_allclose(results, ref, rtol=1e-5)

    else:
        # Update reference file
        results.to_netcdf(ref_properties_vt, engine='h5netcdf')
