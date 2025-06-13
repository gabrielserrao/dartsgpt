import numpy as np
import pickle
import pytest

import dartsflash.libflash
from .conftest import compdata, compdata_ions


# Tests of class CompData
def test_compdata(compdata):
    assert isinstance(compdata, dartsflash.libflash.CompData)
    assert compdata.nc == 3
    assert compdata.ni == 0
    assert compdata.ns == 3


def test_compdata_ions(compdata_ions):
    assert isinstance(compdata_ions, dartsflash.libflash.CompData)
    assert compdata_ions.nc == 3
    assert compdata_ions.ni == 2
    assert compdata_ions.ns == 5


# def test_eos_properties():
#
#
# def test_enthalpy_evaluate(enthalpy_ideal_object):
#     Hi = enthalpy_ideal_object.evaluate(100, [1e-6, 2e-6])  # @Michiel check these values
#     assert Hi == pytest.approx(-1.0491059264390992)  # It also works in the other way
#
#
# # Tests of class enthalpy CompData
#
#
# def test_CompData_twu_dij(comp_data_object, ref_comp_data_twu_dij):
#     dij = comp_data_object.twu_dij()  # put parameters necessary
#     np.testing.assert_allclose(ref_comp_data_twu_dij, dij, rtol=1e-4)