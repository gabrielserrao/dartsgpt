import pytest

from dartsflash.components import CompData
from dartsflash.mixtures import Mixture


# Define pytest.fixtures: Mocking objects only necessary in the test of this files.
# If one of this is needed in other files, consider defining it in conftest.py

# If you have a function that generate this data and is used in multiple files, you can also put it in conftest.py.
# In that case (function): you may have to import the file as from tests.python.conftest import function_in_conftest

@pytest.fixture()
def compdata() -> CompData:
    return CompData(components=["H2O", "CO2", "C1"], setprops=True)


@pytest.fixture()
def compdata_ions() -> CompData:
    return CompData(components=["H2O", "CO2", "C1"], ions=["Na+", "Cl-"], setprops=True)


@pytest.fixture()
def mixture_brine_co2() -> Mixture:
    return Mixture(components=["H2O", "CO2"], setprops=True)


@pytest.fixture()
def mixture_brine_co2_ions() -> Mixture:
    return Mixture(components=["H2O", "CO2"], ions=["Na+", "Cl-"], setprops=True)


@pytest.fixture()
def mixture_brine_vapour() -> Mixture:
    return Mixture(components=["H2O", "CO2", "C1"], setprops=True)


@pytest.fixture()
def mixture_brine_vapour_ions() -> Mixture:
    return Mixture(components=["H2O", "CO2", "C1"], ions=["Na+", "Cl-"], setprops=True)
