import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from dartsflash.libflash import FlashParams, EoSParams, EoS, InitialGuess
from dartsflash.libflash import Flash

from dartsflash.mixtures import Mixture
from dartsflash.pyflash import PyFlash


if 1:
    mix = Mixture(components=["H2O"])

    from dartsflash.libflash import CubicEoS, AQEoS, IAPWS_Ice, PureSolid
    pr = CubicEoS(mix.comp_data, CubicEoS.PR)
    pr.set_root_flag(EoS.MAX)
    prL = CubicEoS(mix.comp_data, CubicEoS.PR)
    prL.set_root_flag(EoS.MIN)
    aq = AQEoS(mix.comp_data, {AQEoS.water: AQEoS.Jager2003,
                               AQEoS.solute: AQEoS.Ziabakhsh2012})
    # ice = IAPWS_Ice(mix.comp_data)
    ice = PureSolid(mix.comp_data, "Ice")

    # Calculate triple point (V-Aq-I)
    triplePT = np.array([6e-3, 273.15])
    tol = 1e-12
    max_iter = 500
    it = 0
    while True:
        aq.solve_PT(triplePT[0], triplePT[1], [1.])
        gA = aq.G(triplePT[0], triplePT[1], [1.], 0, pt=True)
        dgA_dT = aq.dlnphi_dT()
        dgA_dP = aq.dlnphi_dP()

        pr.solve_PT(triplePT[0], triplePT[1], [1.])
        gV = pr.G(triplePT[0], triplePT[1], [1.], 0, pt=True)
        dgV_dT = pr.dlnphi_dT()
        dgV_dP = pr.dlnphi_dP()

        ice.solve_PT(triplePT[0], triplePT[1], [1.])
        gS = ice.G(triplePT[0], triplePT[1], [1.], 0, pt=True)
        dgS_dT = ice.dlnphi_dT()
        dgS_dP = ice.dlnphi_dP()

        res = np.array([gA - gV, gV - gS])
        Jac = np.array([[dgA_dP[0] - dgV_dP[0], dgA_dT[0] - dgV_dT[0]],
                        [dgV_dP[0] - dgS_dP[0], dgV_dT[0] - dgS_dT[0]]])
        triplePT -= np.linalg.solve(Jac, res)

        if np.linalg.norm(res) < tol:
            print(gA, gV, gS)
            break
        elif it > max_iter:
            print("max iter")
            break
    print("triple point", triplePT)

    pres_SV = np.logspace(-14, np.log10(triplePT[0]), 100)
    pres_SV = np.array(list(reversed(pres_SV)))
    T_sub = np.empty(np.shape(pres_SV))
    T_sub0 = triplePT[1]

    # pres_SV = pres_SV[65:]
    for i, p in enumerate(pres_SV):
        # Calculate sublimation curve
        it = 0

        T_sub[i] = T_sub0
        while True:
            pr.solve_PT(p, T_sub[i], [1.])
            gV = pr.G(p, T_sub[i], [1.], 0, pt=True)
            dgV_dT = pr.dlnphi_dT()

            ice.solve_PT(p, T_sub[i], [1.])
            gS = ice.G(p, T_sub[i], [1.], 0, pt=True)
            dgS_dT = ice.dlnphi_dT()

            res = gV - gS
            dres_dT = dgV_dT[0] - dgS_dT[0]

            T_sub[i] -= res / dres_dT

            if np.isnan(T_sub[i]):
                print("S-V NANs", it)
                break
            if np.abs(res) < tol:
                # T_sub0 = T_sub[i]
                break
            elif it > max_iter:
                print("S-V it >", max_iter, T_sub0, T_sub[i])
                T_sub0 = T_sub[i]
                break
            it += 1

    pres = np.logspace(np.log10(triplePT[0]), np.log10(140.), 100)
    arr_shape = np.shape(pres)
    T_sol, T_vap, T_vapPR = np.empty(arr_shape), np.empty(arr_shape), np.empty(arr_shape)
    T_sol0, T_vap0, T_vapPR0 = 273.15, 273.15, 273.15

    for i, p in enumerate(pres):
        # Calculate evaporation curve Aq-V
        it = 0

        T_vap[i] = T_vap0
        while True:
            gA = aq.G(p, T_vap[i], [1.], 0, pt=True)
            dgA_dT = aq.dlnphi_dT()

            gV = pr.G(p, T_vap[i], [1.], 0, pt=True)
            dgV_dT = pr.dlnphi_dT()

            res = gA - gV
            dres_dT = dgA_dT[0] - dgV_dT[0]

            T_vap[i] -= 0.1 * res/dres_dT

            if np.abs(res) < tol:
                T_vap0 = T_vap[i]
                break
            elif it > max_iter:
                print("A-V it >", max_iter)
                T_vap0 = T_vap[i]
                break
            it += 1

        # Calculate evaporation curve L-V full PR
        it = 0

        T_vapPR[i] = np.nan
        # T_vapPR[i] = T_vapPR0
        # while True:
        #     gA = prL.G(p, T_vapPR[i], [1.], 0, pt=True)
        #     dgA_dT = prL.dlnphi_dT()
        #
        #     gV = pr.G(p, T_vapPR[i], [1.], 0, pt=True)
        #     dgV_dT = pr.dlnphi_dT()
        #
        #     res = gA - gV
        #     dres_dT = dgA_dT[0] - dgV_dT[0]
        #
        #     T_vapPR[i] -= res / dres_dT
        #
        #     if np.abs(res) < tol:
        #         T_vapPR0 = T_vapPR[i]
        #         break
        #     elif it > max_iter:
        #         print("PR-VL it >", max_iter)
        #         T_vapPR0 = T_vapPR[i]
        #         break
        #     it += 1

        # Calculate solid curve
        it = 0

        T_sol[i] = T_sol0
        while True:
            gA = aq.G(p, T_sol[i], [1.], 0, pt=True)
            dgA_dT = aq.dlnphi_dT()

            gS = ice.G(p, T_sol[i], [1.], 0, pt=True)
            dgS_dT = ice.dlnphi_dT()

            res = gA - gS
            dres_dT = dgA_dT[0] - dgS_dT[0]

            T_sol[i] -= 0.1 * res / dres_dT

            if np.abs(res) < tol:
                T_sol0 = T_sol[i]
                break
            elif it > max_iter:
                print("V-S it >", max_iter)
                T_sol0 = T_sol[i]
                break
            it += 1

    from dartsflash.diagram import Plot
    plot = Plot()
    plot.add_attributes(suptitle=r"Phase diagram of H$_2$O", ax_labels=["Temperature [K]", "Pressure [bar]"], grid=True)

    plot.draw_plot(xdata=[T_sub, T_sol, T_vap], ydata=[pres_SV, pres, pres], logy=True)
    plot.draw_plot(xdata=T_vapPR, ydata=pres, style='dotted', colour=plot.colours[2])
    plot.draw_plot(xdata=[triplePT[1]], ydata=[triplePT[0]], plot_type="scatter")

    plt.savefig("H2O-PT.pdf")

plt.show()
