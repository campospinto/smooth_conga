# -*- coding: UTF-8 -*-

import numpy as np
from scipy.integrate import quad

# from spl.core.interface import make_open_knots
from local_projections import SplineProjectionOperators

# import os.path
# import spl.core as core
# print ('writing path : ')
# print (os.path.abspath(core.__file__))
# exit()

import matplotlib
matplotlib.use('Agg')   # backend configuration: high quality PNG images using the Anti-Grain Geometry engine
import matplotlib.pyplot as plt



####################################################################################
if __name__ == '__main__':


    p = 2
    m = 1
    M = p+m+1
    ncs = 2
    ncs_max = max(ncs, 3*p+m+1)
    ncs_good = int(np.ceil(ncs_max/M))*M

    print("ncs = ", ncs)
    print("ncs_max = ", ncs_max)
    print("ncs_good = ", ncs_good)

    fk = 'pol1'
    # fk = 'pol'
    proj_kinds = ['db', 'ps', 'ms']
    use_macro_elem_duals = ('ms' in proj_kinds)

    # proj_types = ['P','D', 'M'] #,'M']
    # use_macro_elem_duals = ('M' in proj_types)

    spo = SplineProjectionOperators(p, m, N_cells_sub=ncs_good, N_subdomains=2, use_macro_elem_duals=use_macro_elem_duals)

    print()

    if fk == 'runge':
        f = lambda x: 1./((5*( 2*(x-0.5)))**2+1)
    elif fk == 'pol1':
        f = lambda x: x
    elif fk == 'pol2':
        f = lambda x: x*(x-0.5)
    elif fk == 'abs':
        f = lambda x: abs(x-0.5)
    elif fk == 'pws':
        f = lambda x: int(x > 0.5)* np.sin(x*np.pi)
    elif fk == 'pwh':
        f = lambda x: int(x > 0.5)* np.sin(5*x*np.pi)
    elif fk == 'step':
        f = lambda x: int(x > 0.5)
    else:
        raise ValueError("unkown value for f_kind: ", fk)

    for pk in proj_kinds:

        print(" - - - - - - - - - - - - - - - - - - - -")
        print("## - tests with proj = ", pk)

        # A - L2 proj on tilde Vh
        print("## A - L2 proj of f="+fk+" on tilde V_h...")
        spo.l2_proj_tilde_V(f)

        error = spo.plot_spline(
            filename="tfh_L2_"+pk+".png",
            f_ref= f,
            N_points=200,
            # spline_kind="discontinuous",
            iqnorm=0.5,
            # save_plot=True,
        )
        print("tfh - f, error = ", error)

        # B - try the smooth projection on tfh
        print("## B - smooth projection in Vh...")
        tc = spo.tilde_coefs
        P = spo.proj_matrix_from_tV_to_V(kind=pk)
        # P = spo.get_smooth_proj_on_tilde_V(kind=pk)
        c = P.dot(tc)
        spo.set_coefs(c)

        error = spo.plot_spline(
            filename="fh_"+pk+".png",
            f_ref= f,
            N_points=200,
            # spline_kind="continuous",
            iqnorm=0.5,
            # save_plot=True,
        )
        print("fh - f, error = ", error)

        # C - P* approx on tilde Vh
        if 1:
            print("## C - approx f="+fk+" with P* on tilde V_h, using proj="+pk+"...")
            spo.P_star_tilde_V(f, kind=pk)
            error = spo.plot_spline(
                filename="tfh_Pstar_"+pk+".png",
                f_ref= f,
                N_points=200,
                # spline_kind="discontinuous",
                iqnorm=0.5,
                # save_plot=True,
            )
            print("tfh - f, error = ", error)



        print(" - - - - - - - - - - - - - - - - - - - -")