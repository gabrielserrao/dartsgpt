import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from dartsflash.libflash import CubicEoS

from dartsflash.pyflash import PyFlash
from dartsflash.mixtures import Mixture
from dartsflash.plot import *


vmax = None
if 1:
    components = ["CO2"]
    n = [1.]
    pressure = np.arange(47.68, 151, 20)
    # pressure = np.array([60.])
    temperature = np.array([273.15, 293.15, 304.10, 313.15, 418, 500])
    # temperature = np.array([203.15])
    temperaturePT = np.arange(403.15, 434.15, 10)
    vmax = 6e-3
    nvol = 1001
    prange = [-25, 200]
elif 0:
    components = ["H2O"]
    n = [1.]
    pressure = np.arange(1, 300, 0.1)
    temperature = np.array([373, 473, 573, 673])
    temperaturePT = np.arange(273.15, 743.15, 10)
    vmax = 6e-3
    nvol = 3001
    prange = [-25, 300]
elif 1:
    components = ["H2O", "CO2"]
    n = [0.9, 0.1]
    pressure = np.arange(1, 300, 1)
    temperature = np.array([423, 473, 523, 573, 623, 673])
    temperaturePT = np.arange(273.15, 643.15, 10)
    vmax = 6e-3
    nvol = 3001
    prange = [-25, 300]
elif 1:
    components = ["H2O", "H2S"]
    n = [0.9, 0.1]
    pressure = np.arange(1, 300, 1)
    temperature = np.array([423, 473, 523, 573, 623, 673])
    vmax = 6e-3
    nvol = 3001
    prange = [-25, 300]
elif 0:
    components = ["H2O", "CO2"]
    n = [np.linspace(0., 1., 11), 1.]
    pressure = np.arange(1, 300, 1)
    temperature = np.array([304.1, 320.5, 340.5, 365.0, 392.8, 425.6, 462.6, 503.1, 547.5, 595.2, 647.5])
    vmax = 6e-3
    nvol = 3001
    prange = [-25, 300]
else:
    components = ["H2O", "H2S"]
    n = [np.linspace(0., 1., 11), 1.]
    pressure = np.arange(1, 300, 1)
    temperature = np.array([373.5, 391.2, 411.2, 433.5, 457.9, 485.4, 513.2, 543.9, 576.9, 611.3, 647.5])
    vmax = 6e-3
    nvol = 3001
    prange = [-25, 300]

mix = Mixture(components)
f = PyFlash(mixture=mix)

f.add_eos("CEOS", CubicEoS(mix.comp_data, CubicEoS.PR))
ceos = f.eos["CEOS"]

if 0:
    """ Compressibility factor and P-V diagram of Cubic EoS """
    """ Evaluate properties at PT """
    state_spec = {"pressure": pressure,
                  "temperature": temperature
                  }
    compositions = {comp: n[i] for i, comp in enumerate(components)}
    properties = {"V": ceos.V,
                  "Z": ceos.Zp,
                  }

    pprops = f.evaluate_properties_1p(state_spec=state_spec, compositions=compositions,
                                      properties_to_evaluate=properties, mole_fractions=True)

    """ Evaluate properties at VT """
    vmin = ceos.V(p=pressure[-1], T=temperature[0], n=n)
    vmin = 1e-4
    vmax = ceos.V(p=pressure[0], T=temperature[-1], n=n) if vmax is None else vmax

    state_spec = {"temperature": temperature,
                  "volume": np.linspace(vmin, vmax, nvol)
                  }
    properties = {"Z": ceos.Zv,
                  "P": ceos.P
                  }

    vprops = f.evaluate_properties_1p(state_spec=state_spec, compositions=compositions,
                                      properties_to_evaluate=properties, mole_fractions=True)

    """ Plot PV and PZ diagrams """
    pv = PlotEoS.pressure_volume(f, temperatures=temperature, compositions=n,
                                 pt_props=pprops, vt_props=vprops, vrange=[0, vmax], prange=[0, prange[1]])

    pz = PlotEoS.compressibility(f, temperatures=temperature, compositions=n,
                                 pt_props=None, vt_props=vprops, zrange=[-0.1, 1.1],) #prange=prange)
    plt.savefig("P-Z-" + mix.filename + ".pdf")

if 0:
    """ Cubic polynomial and phase identification """
    state_spec = {"pressure": np.array([150, 250, 350, 500]),
                  "temperature": 399.5916961
                  }
    compositions = {comp: n[i] for i, comp in enumerate(components)}
    properties = {"coeff": ceos.calc_coefficients,
                  }

    pprops = f.evaluate_properties_1p(state_spec=state_spec, compositions=compositions,
                                      properties_to_evaluate=properties, mole_fractions=True)

    """ Plot f-Z diagrams """
    def F(Z):
        # Slice Dataset at composition
        comps = {comp: compositions[i] for i, comp in enumerate(f.components[:-1]) if comp in pprops.dims}

        coeff = pprops.sel(comps, method='nearest').coeff.transpose('temperature', ...).values
        return np.transpose(np.array([z ** 3 + coeff[..., 0] * z ** 2 + coeff[..., 1] * z + coeff[..., 2] for z in Z]), axes=(1, 2, 0))

    # Initialize Plot object
    from dartsflash.diagram import Plot
    plot = Plot(figsize=(8, 4), nrows=1, ncols=len(temperature))
    plot.add_attributes(suptitle="Cubic polynomial of " + f.mixture.name +
                                 " at T = {} K".format(state_spec['temperature']), ax_labels=["Z", "f"])

    Z = np.linspace(0.01, 1., 100)
    y = F(Z)
    for i, temp in enumerate(temperature):
        plot.subplot_idx = i
        plot.draw_plot(xdata=Z, ydata=y[i], style="solid", xlim=[0, 1],# ylim=[-0.1, 0.1],
                       datalabels=['{} bar'.format(pres) for pres in state_spec['pressure']])
        plot.add_attributes(legend=True)
        plot.draw_line(x=[0., 1.], y=[0., 0.], colours='k', linestyle='dashed')

    plt.savefig("f-Z-" + mix.filename + ".pdf")
    plt.show()

if 1:
    """ PT-diagram of properties """
    """ Evaluate properties at PT """
    state_spec = {"pressure": np.arange(10, 701, 1e0),
                  "temperature": np.arange(100, 700, 1e0)
                  }
    compositions = {comp: n[i] for i, comp in enumerate(components)}
    properties = {"Volume": ceos.V,
                  "Root": ceos.is_root_type,
                  "CriticalT": ceos.is_critical,
                  }

    pprops = f.evaluate_properties_1p(state_spec=state_spec, compositions=compositions,
                                      properties_to_evaluate=properties, mole_fractions=True)
    cp = ceos.critical_point(n)

    plot = PlotEoS.surf(f, props=pprops, x_var='temperature', y_var='pressure', prop_names=list(properties.keys()),
                        composition=n)
    for i in range(len(properties.keys())):
        plot.subplot_idx = i
        plot.draw_point(X=cp.Tc, Y=cp.Pc, colours='red')

plt.show()
