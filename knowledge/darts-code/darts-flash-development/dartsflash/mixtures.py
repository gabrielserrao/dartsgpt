import numpy as np
from dartsflash.components import CompData
from enum import Enum


class ConcentrationUnits(Enum):
    MOLALITY = 0
    WEIGHT = 1


class Mixture:
    def __init__(self, components: list, ions: list = None, setprops: bool = True,
                 name: str = None, composition: list = None):
        self.comp_data = CompData(components, ions, setprops=setprops)
        self.composition = composition
        self.nc = self.comp_data.nc

        self.filename = "-".join(comp for comp in components)
        self.name = "-".join(label for label in self.comp_data.comp_labels[:self.nc]) if name is None else name

        self.H2O_idx = components.index("H2O") if "H2O" in components else None
        self.salt_stoich = {"NaCl": {0: 1, 1: 1}, "CaCl2": {0: 1, 1: 2}, "KCl": {0: 1, 1: 1}}
        self.salt_mass = {"NaCl": 53.99, "CaCl2": 110.98, "KCl": 74.5513}

    def set_binary_coefficients(self, i: int, kij: list):
        self.comp_data.set_binary_coefficients(i, kij)

    def calculate_concentrations(self, ni: np.ndarray, mole_fractions: bool = False, concentrations: dict = None,
                                 concentration_unit: ConcentrationUnits = ConcentrationUnits.MOLALITY):
        # Translate concentration into composition of dissolved components
        nH2O = ni[self.H2O_idx]
        ni = np.append(ni, np.zeros(self.comp_data.ni))
        if concentration_unit == ConcentrationUnits.WEIGHT:
            M_H2O = nH2O * self.comp_data.Mw[self.H2O_idx]
            for comp, ci in concentrations.items():
                Mw_comp = self.salt_mass[comp]
                Ni = 1. / Mw_comp * M_H2O * (1./(1.-ci) - 1.)  # weight fraction to mole number conversion
                for i, stoich in self.salt_stoich[comp].items():
                    ni[self.nc + i] = Ni * stoich

        elif concentration_unit == ConcentrationUnits.MOLALITY:
            for comp, ci in concentrations.items():
                Ni = ci * nH2O / 55.509
                for i, stoich in self.salt_stoich[comp].items():
                    ni[self.nc + i] = Ni * stoich

        if mole_fractions:
            ni = ni / np.sum(ni)
        return ni


class M7(Mixture):
    """
    Seven-component gas mixture (Michelsen, 1982a fig. 2)
    """
    def __init__(self):
        super().__init__(name="M7", components=["C1", "C2", "C3", "nC4", "nC5", "nC6", "N2"], setprops=False,
                         composition=[0.9430, 0.0270, 0.0074, 0.0049, 0.0027, 0.0010, 0.0140])
        self.comp_data.Pc = [45.99, 48.72, 42.48, 33.70, 27.40, 21.10, 34.00]
        self.comp_data.Tc = [190.56, 305.32, 369.83, 469.70, 540.20, 617.70, 126.20]
        self.comp_data.ac = [0.011, 0.099, 0.152, 0.252, 0.350, 0.490, 0.0377]
        self.comp_data.Mw = [16.043, 30.07, 44.097, 58.124, 72.151, 86.178, 28.013]
        self.comp_data.kij = np.zeros(self.nc * self.nc)


class M3(Mixture):
    """
    Ternary mixture (Michelsen, 1982a fig. 4)
    """
    def __init__(self):
        super().__init__(components=["C1", "CO2", "H2S"], setprops=False,
                         composition=[0.50, 0.10, 0.40])
        self.comp_data.Pc = [46.04, 73.75, 89.63]
        self.comp_data.Tc = [190.58, 304.10, 373.53]
        self.comp_data.ac = [0.012, 0.239, 0.0942]
        self.comp_data.Mw = [16.043, 44.01, 34.1]
        self.comp_data.kij = np.zeros(self.nc * self.nc)


class Y8(Mixture):
    """
    Y8 mixture
    """
    def __init__(self):
        super().__init__(name="Y8", components=["C1", "C2", "C3", "nC5", "nC7", "nC10"], setprops=False,
                         composition=[0.8097, 0.0566, 0.0306, 0.0457, 0.0330, 0.0244])
        self.comp_data.Pc = [45.99, 48.72, 42.48, 33.70, 27.40, 21.10]
        self.comp_data.Tc = [190.56, 305.32, 369.83, 469.70, 540.20, 617.70]
        self.comp_data.ac = [0.011, 0.099, 0.152, 0.252, 0.350, 0.490]
        self.comp_data.Mw = [16.043, 30.07, 44.097, 72.151, 100.205, 142.2848]
        self.comp_data.kij = np.zeros(self.nc * self.nc)


class MY10(Mixture):
    """
    MY10 mixture
    """
    def __init__(self):
        super().__init__(name="MY10", components=["C1", "C2", "C3", "nC4", "nC5", "nC6", "nC7",
                                                  "nC8", "nC10", "nC14"], setprops=False,
                         composition=[0.35, 0.03, 0.04, 0.06, 0.04, 0.03, 0.05, 0.05, 0.30, 0.05])
        self.comp_data.Pc = [45.99, 48.72, 42.48, 37.96, 33.70, 30.25, 27.40, 24.9, 21.10, 15.7]
        self.comp_data.Tc = [190.56, 305.32, 369.83, 425.12, 469.70, 507.6, 540.20, 568.7, 617.70, 693.0]
        self.comp_data.ac = [0.011, 0.099, 0.152, 0.2, 0.252, 0.3, 0.350, 0.399, 0.490, 0.644]
        self.comp_data.Mw = [16.043, 30.07, 44.097, 58.124, 72.151, 86.178, 100.205, 114.231, 142.2848, 0.0]
        self.comp_data.kij = np.zeros(self.nc * self.nc)
        self.set_binary_coefficients(0, [0., 0., 0., 0.02, 0.02, 0.025, 0.025, 0.035, 0.045, 0.045])


class DepletedOilMixture(Mixture):
    """
    Depleted oil field with injected CO2 mixture
    """
    def __init__(self):
        super().__init__(components=["C1", "nC10", "CO2"], setprops=False)
        self.comp_data.Pc = [46.0, 21.2, 73.75]
        self.comp_data.Tc = [190.58, 617.7, 304.10]
        self.comp_data.ac = [0.012, 0.489, 0.239]
        self.comp_data.Mw = [16.043, 142.2848, 44.01]
        self.comp_data.kij = np.zeros(self.nc * self.nc)
        self.set_binary_coefficients(0, [0., 0.048388, 0.0936])
        self.set_binary_coefficients(1, [0.048388, 0., 0.1])


class OilB(Mixture):
    """
    Oil B mixture (Shelton and Yarborough, 1977), data from (Li, 2012)
    """
    def __init__(self):
        super().__init__(name="OilB", components=["CO2", "N2", "C1", "C2", "C3", "iC4", "nC4", "iC5", "nC5", "C6",
                                                  "PC1", "PC2", "PC3", "PC4", "PC5", "PC6"], setprops=False,
                         composition=[0.0011, 0.0048, 0.1630, 0.0403, 0.0297, 0.0036, 0.0329, 0.0158, 0.0215, 0.0332,
                                      0.181326, 0.161389, 0.125314, 0.095409, 0.057910, 0.022752])
        self.comp_data.Pc = [73.819, 33.5, 45.4, 48.2, 41.9, 36., 37.5, 33.4, 33.3, 33.9, 25.3, 19.1, 14.2, 10.5, 7.5, 4.76]
        self.comp_data.Tc = [304.211, 126.2, 190.6, 305.4, 369.8, 408.1, 425.2, 460.4, 469.6, 506.35, 566.55, 647.06, 719.44, 784.93, 846.33, 919.39]
        self.comp_data.ac = [0.225, 0.04, 0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.299, 0.3884, 0.5289, 0.6911, 0.8782, 1.1009, 1.4478]
        self.comp_data.Mw = [44.01, 28.01, 16.04, 30.07, 44.1, 58.12, 58.12, 72.15, 72.15, 84., 112.8, 161.2, 223.2, 304.4, 417.5, 636.8]
        self.comp_data.kij = np.zeros(self.nc * self.nc)
        self.set_binary_coefficients(0, [0., -0.02, 0.075, 0.08, 0.08, 0.085, 0.085, 0.085, 0.085, 0.095, 0.095, 0.095, 0.095, 0.095, 0.095, 0.095])
        self.set_binary_coefficients(1, [-0.02, 0., 0.08, 0.07, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.1, 0.12, 0.12, 0.12, 0.12, 0.12])
        self.set_binary_coefficients(2, [0.075, 0.08, 0., 0.003, 0.01, 0.018, 0.018, 0.025, 0.026, 0.036, 0.049, 0.073, 0.098, 0.124, 0.149, 0.181])


class MaljamarRes(Mixture):
    """
    Maljamar reservoir mixture (Orr, 1981), data from (Li, 2012)
    """
    def __init__(self):
        super().__init__(name="MaljamarRes", components=["CO2", "C1", "C2", "C3", "nC4", "C5-7", "C8-10", "C11-14",
                                                         "C15-20", "C21-28", "C29+"], setprops=False,
                         composition=[0., 0.2939, 0.1019, 0.0835, 0.0331, 0.1204, 0.1581, 0.0823,
                                      0.0528, 0.0276, 0.0464])
        self.comp_data.Pc = [73.819, 45.4, 48.2, 41.9, 37.5, 28.82, 23.743, 18.589, 14.8, 11.954, 8.523]
        self.comp_data.Tc = [304.211, 190.6, 305.4, 369.8, 425.2, 516.667, 590., 668.611, 745.778, 812.667, 914.889]
        self.comp_data.ac = [0.225, 0.008, 0.098, 0.152, 0.193, 0.2651, 0.3644, 0.4987, 0.6606, 0.8771, 1.2789]
        self.comp_data.Mw = [44.01, 16.043, 30.1, 44.1, 58.1, 89.9, 125.7, 174.4, 240.3, 336.1, 536.7]
        self.comp_data.kij = np.zeros(self.nc * self.nc)
        self.set_binary_coefficients(0, [0., 0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115])
        self.set_binary_coefficients(1, [0.115, 0., 0., 0., 0., 0.045, 0.055, 0.055, 0.06, 0.08, 0.28])


class MaljamarSep(Mixture):
    """
    Maljamar separator mixture (Orr, 1981), data from (Li, 2012)
    """
    def __init__(self):
        super().__init__(name="MaljamarSep", components=["CO2", "C5-7", "C8-10", "C11-14", "C15-20",
                                                         "C21-28", "C29+"], setprops=False,
                         composition=[0.0, 0.2354, 0.3295, 0.1713, 0.1099, 0.0574, 0.0965])
        self.comp_data.Pc = [73.9, 28.8, 23.7, 18.6, 14.8, 12.0, 8.5]
        self.comp_data.Tc = [304.2, 516.7, 590.0, 668.6, 745.8, 812.7, 914.9]
        self.comp_data.ac = [0.225, 0.265, 0.364, 0.499, 0.661, 0.877, 1.279]
        self.comp_data.Mw = [44.01, 89.9, 125.7, 174.4, 240.3, 336.1, 536.7]
        self.comp_data.kij = np.zeros(self.nc * self.nc)
        self.set_binary_coefficients(0, [0.0, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115])


class SourGas(Mixture):
    """
    Sour gas mixture, data from (Li, 2012)
    """
    def __init__(self):
        super().__init__(name="SourGas", components=["CO2", "N2", "H2S", "C1", "C2", "C3"], setprops=False,
                         composition=[0.70592, 0.07026, 0.01966, 0.06860, 0.10559, 0.02967])
        self.comp_data.Pc = [73.819, 33.9, 89.4, 45.992, 48.718, 42.462]
        self.comp_data.Tc = [304.211, 126.2, 373.2, 190.564, 305.322, 369.825]
        self.comp_data.ac = [0.225, 0.039, 0.081, 0.01141, 0.10574, 0.15813]
        self.comp_data.Mw = [44.01, 28.013, 34.1, 16.043, 30.07, 44.097]
        self.comp_data.kij = np.zeros(self.nc * self.nc)
        self.set_binary_coefficients(0, [0., -0.02, 0.12, 0.125, 0.135, 0.150])
        self.set_binary_coefficients(1, [-0.02, 0., 0.2, 0.031, 0.042, 0.091])
        self.set_binary_coefficients(2, [0.12, 0.2, 0., 0.1, 0.08, 0.08])


class BobSlaughterBlock(Mixture):
    """
    Bob Slaughter Block mixture, data from (Li, 2012)
    """
    def __init__(self):
        super().__init__(name="BSB", components=["CO2", "C1", "PC1", "PC2"], setprops=False,
                         composition=[0.0337, 0.0861, 0.6478, 0.2324])
        self.comp_data.Pc = [73.77, 46., 27.32, 17.31]
        self.comp_data.Tc = [304.2, 160., 529.03, 795.33]
        self.comp_data.ac = [0.225, 0.008, 0.481, 1.042]
        self.comp_data.Mw = [44.01, 16.043, 98.45, 354.2]
        self.comp_data.kij = np.zeros(self.nc * self.nc)
        self.set_binary_coefficients(0, [0., 0.055, 0.081, 0.105])


class NorthWardEstes(Mixture):
    """
    NorthWardEstes mixture, data from (Li, 2012)
    """
    def __init__(self):
        super().__init__(name="NWE", components=["CO2", "C1", "PC1", "PC2", "PC3", "PC4", "PC5"], setprops=False,
                         composition=[0.0077, 0.2025, 0.118, 0.1484, 0.2863, 0.149, 0.0881])
        self.comp_data.Pc = [73.77, 46., 45.05, 33.51, 24.24, 18.03, 17.26]
        self.comp_data.Tc = [304.2, 190.6, 343.64, 466.41, 603.07, 733.79, 923.2]
        self.comp_data.ac = [0.225, 0.008, 0.13, 0.244, 0.6, 0.903, 1.229]
        self.comp_data.Mw = [44.01, 16.04, 38.4, 72.82, 135.82, 257.75, 479.95]
        self.comp_data.kij = np.zeros(self.nc * self.nc)
        self.set_binary_coefficients(0, [0., 0.12, 0.12, 0.12, 0.09, 0.09, 0.09])
