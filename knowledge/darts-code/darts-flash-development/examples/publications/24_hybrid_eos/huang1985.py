import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from dartsflash.libflash import FlashParams, EoSParams, EoS, InitialGuess
from dartsflash.libflash import Flash
from dartsflash.libflash import CubicEoS, AQEoS

from dartsflash.mixtures import Mixture
from dartsflash.pyflash import PyFlash
from dartsflash.plot import *


# Huang et al. (1985) H2O-C1-CO2-H2S mixtures
components = ["H2O", "C1", "CO2", "H2S"]

mix = Mixture(components, setprops=True)
f = PyFlash(mix)

f.add_eos("CEOS", CubicEoS(mix.comp_data, CubicEoS.PR),
          initial_guesses=[0, 1, 2, 3], preferred_roots=[(0, 0.75, EoS.MAX)],
          switch_tol=1e-1, stability_tol=1e-20, max_iter=50, use_gmix=0)
f.add_eos("AQ", AQEoS(mix.comp_data, {AQEoS.water: AQEoS.Jager2003,
                                      AQEoS.solute: AQEoS.Ziabakhsh2012}),
          initial_guesses=[0], eos_range={0: [0.6, 1.]}, max_iter=10, use_gmix=True)

f.flash_params.stability_variables = FlashParams.alpha
f.flash_params.split_variables = FlashParams.nik
f.flash_params.split_tol = 1e-20
# f.flash_params.split_switch_tol = 1e-2
f.flash_params.tpd_close_to_boundary = 1e-2
f.flash_params.split_negative_flash_iter = 3
f.flash_params.tpd_tol = 1e-11
f.flash_params.verbose = 0

eos_order = ["AQ", "CEOS"]
np_max = 3

# Ref data mixture 1
if 1:
    composition = [0.5, 0.15, 0.30, 0.05]  # mixture 1

    ref_p = [[48.2, 76.0, 125.2, 169.3], [83.6, 129.3, 171.7], [118.0, 173.1]]
    ref_T = np.array([37.8, 107.2, 176.7]) + 273.15
    colours = [Diagram.colours[i] for i in [0, 2, 3]]

    ref_h2o_V = [[1.91e-3, 1.71e-3, 1.87e-3, 1.99e-3], [0.0225, 0.0196, 0.0179], [0.0950, 0.0848]]
    ref_c1_V = [[0.3040, 0.3031, 0.3029, 0.3021], [0.2907, 0.2929, 0.2935], [0.2641, 0.2762]]
    ref_co2_V = [[0.5945, 0.5970, 0.5967, 0.5963], [0.5919, 0.5920, 0.5916], [0.5575, 0.5520]]
    ref_h2s_V = [[0.0998, 0.0982, 0.0985, 0.0996], [0.0967, 0.0963, 0.0970], [0.0848, 0.0870]]

    ref_h2o_A = [[0.9854, 0.9816, 0.9781, 0.9777], [0.9894, 0.9854, 0.9834], [0.9877, 0.9827]]
    ref_c1_A = [[2.76e-4, 4.66e-4, 7.96e-4, 9.90e-4], [3.79e-4, 5.78e-4, 7.79e-4], [7.17e-4, 1.10e-3]]
    ref_co2_A = [[9.30e-3, 0.0121, 0.0151, 0.0154], [6.98e-3, 9.59e-3, 0.0113], [7.95e-3, 0.0114]]
    ref_h2s_A = [[5.03e-3, 5.40e-3, 5.95e-3, 6.08e-3], [3.42e-3, 4.47e-3, 4.73e-3], [3.74e-3, 4.78e-3]]

    ref_h2o_L = [[None, None, None, None], [None, None, None], [None, None]]
    ref_c1_L = [[None, None, None, None], [None, None, None], [None, None]]
    ref_co2_L = [[None, None, None, None], [None, None, None], [None, None]]
    ref_h2s_L = [[None, None, None, None], [None, None, None], [None, None]]

    ref_dew_p = [[32.5, 27.0, 19.4, 16.7, 13.9, 9.5, 6.4, 4.2]]
    ref_dew_T = np.array([[198.8, 190.6, 176.7, 169.8, 163.1, 148.8, 135.0, 121.1]]) + 273.15
    ref_labels = ["Aq-V dew point"]
    mixname = "1"
else:
    composition = [0.5, 0.05, 0.05, 0.40]  # mixture 2

    ref_p = [[130.0, 164.6, 62.6], [84.3], [75.6, 122.7, 169.2], [110.0, 181.7]]
    ref_T = np.array([37.8, 65.6, 107.2, 176.7]) + 273.15
    colours = [Diagram.colours[i] for i in [0, 1, 2, 3]]

    ref_h2o_V = [[None, None, 2.14e-3], [8.66e-3], [0.0253, 0.0264, 0.0295], [0.093, 0.113]]
    ref_c1_V = [[None, None, 0.3213], [0.1872], [0.1182, 0.1060, 0.1207], [0.1092, 0.0928]]
    ref_co2_V = [[None, None, 0.1739], [0.1484], [0.1112, 0.1148, 0.1176], [0.1078, 0.0914]]
    ref_h2s_V = [[None, None, 0.5028], [0.6557], [0.7485, 0.7528, 0.7322], [0.6896, 0.704]]

    ref_h2o_A = [[0.9666, 0.9672, 0.9677], [0.9684], [0.9682, 0.9613, 0.9568], [0.9694, 0.9454]]
    ref_c1_A = [[8.59e-4, 8.82e-4, 4.90e-4], [3.85e-4], [1.55e-4, 3.32e-4, 6.06e-4], [3.50e-4, 7.15e-4]]
    ref_co2_A = [[3.62e-3, 3.81e-3, 3.50e-3], [2.72e-3], [1.25e-3, 2.26e-3, 3.34e-3], [1.64e-3, 2.92e-3]]
    ref_h2s_A = [[0.0291, 0.0281, 0.0284], [0.0321], [0.0304, 0.0361, 0.0392], [0.0286, 0.0517]]

    ref_h2o_L = [[9.32e-3, 9.05e-3, 0.0101], [0.0212], [None, None, None], [None, None]]
    ref_c1_L = [[0.0891, 0.0891, 0.0653], [0.0580], [None, None, None], [None, None]]
    ref_co2_L = [[0.0994, 0.1061, 0.1049], [0.0904], [None, None, None], [None, None]]
    ref_h2s_L = [[0.8016, 0.7958, 0.8197], [0.8287], [None, None, None], [None, None]]

    ref_dew_p = [[4.2, 9.7, 19.5, 27.0, 37.0],  # A-V dew point
                 [35.9, 45.0, 57.2, 72.3, 94.1],  # V-L dew point
                 [78.1, 86.5, 94.7, 102.1, 105.2]]  # V-L bubble point
    ref_dew_T = np.array([[120.8, 148.6, 176.7, 190.7, 204.7],
                          [38.7, 49.0, 59.9, 71.0, 82.1],
                          [37.7, 48.8, 59.9, 71.1, 76.0]]) + 273.15
    ref_labels = ["Aq-V dew point", "V-L dew point", "V-L bubble point"]
    mixname = "2"

f.flash_params.verbose = 0
f.init_flash(eos_order=eos_order, stabilityflash=True, np_max=np_max)

compositions = {comp: composition[i] for i, comp in enumerate(components)}

if 1:
    """ PT DIAGRAM """
    N = 2
    state_spec = {"pressure": np.linspace(1, 200, int(N * 200)),
                  "temperature": np.linspace(273.15, 473.15, int(N * 200 + 1)),
                  }

    flash_results = f.evaluate_flash(state_spec=state_spec, compositions=compositions, mole_fractions=True, print_state="Flash")

    plot = PlotFlash.pt(f, flash_results, composition=composition, cmap='RdBu')
    plot.subplot_idx = 0
    plot.draw_point(X=ref_dew_T, Y=ref_dew_p, colours=plot.colours[1:])
    plot.add_attributes(title="Mixture " + mixname)
    plt.savefig("huang_pt" + mixname + ".pdf")

if 1:
    """ SOLUBILITY CURVES """
    state_spec = {"pressure": np.linspace(1, 200, 200),
                  "temperature": ref_T}

    flash_results = f.evaluate_flash(state_spec=state_spec, compositions=compositions, mole_fractions=True, print_state="Flash")

    plotV = PlotFlash.solubility(f, flash_results, dissolved_comp_idx=0, phase_idx=1, x_var='pressure', state=state_spec,
                                 labels=[r"{:.1f} $\degree$C".format(t-273.15) for t in ref_T], plot_1p=False, logy=True)
    plotV.draw_refdata(xref=ref_p, yref=ref_h2o_V)
    plotV.draw_refdata(xref=ref_p, yref=ref_h2o_L, style="^")

    plot_c1_A = PlotFlash.solubility(f, flash_results, dissolved_comp_idx=1, phase_idx=0, x_var='pressure', state=state_spec,
                                     labels=[r"{:.1f} $\degree$C".format(t-273.15) for t in ref_T], plot_1p=False, logy=True,
                                     colours=colours, xlim=[0., 200.], ylim=[1e-5, 2e-3], legend_loc='lower right')
    plot_c1_A.draw_refdata(xref=ref_p, yref=ref_c1_A, colour=colours)
    plt.savefig("huang_sol_c1_" + mixname + ".pdf")

    plot_co2_A = PlotFlash.solubility(f, flash_results, dissolved_comp_idx=2, phase_idx=0, x_var='pressure', state=state_spec,
                                      labels=[r"{:.1f} $\degree$C".format(t - 273.15) for t in ref_T], plot_1p=False, logy=True,
                                      colours=colours, xlim=[0., 200.], ylim=[1e-4, 5e-2], legend_loc='lower right')
    plot_co2_A.draw_refdata(xref=ref_p, yref=ref_co2_A, style="^", colour=colours)
    plt.savefig("huang_sol_co2_" + mixname + ".pdf")

    plot_h2s_A = PlotFlash.solubility(f, flash_results, dissolved_comp_idx=3, phase_idx=0, x_var='pressure', state=state_spec,
                                      labels=[r"{:.1f} $\degree$C".format(t - 273.15) for t in ref_T], plot_1p=False, logy=True,
                                      colours=colours, xlim=[0., 200.], ylim=[1e-4, 5e-2], legend_loc='lower right')
    plot_h2s_A.draw_refdata(xref=ref_p, yref=ref_h2s_A, style="*", colour=colours)
    plt.savefig("huang_sol_h2s_" + mixname + ".pdf")

plt.show()