# -*- coding: UTF-8 -*-

# import pickle

import numpy as np
from scipy.integrate import quad

# from spl.core.interface import make_open_knots
from local_projections import SplineProjectionOperators
from local_projections import is_valid

# import os.path
# import spl.core as core
# print ('writing path : ')
# print (os.path.abspath(core.__file__))
# exit()

import matplotlib
matplotlib.use('Agg')   # backend configuration: high quality PNG images using the Anti-Grain Geometry engine
import matplotlib.pyplot as plt

# pickle_dir = 'pickled_data/'

def file_str(f_kind, p, space_kind='V', proj_kind=None, N_cells=None, iqn=None):
    str = "f="+f_kind+"_p="+repr(p)+"_sp="+space_kind
    if proj_kind is not None:
        str += "_proj="+proj_kind
    if N_cells is not None:
        str += "_ncells="+repr(N_cells)
    if iqn is not None:
        if abs(iqn-0.5) < 1e-5:
            str += "_iqn=0_5"
        else:
            str += "_iqn="+repr(iqn)
    return str


def title_str(f_kind, p, space_kind='V', proj_kind=None, N_cells=None, iqn=None):
    str = "f="+f_kind+" p="+repr(p)+" space="+space_kind
    if proj_kind is not None:
        str += " proj="+proj_kind
    if N_cells is not None:
        str += " N_cells="+repr(N_cells)
    if iqn is not None:
        if abs(iqn) < 1e-5:
            str += " norm = L_inf"
        if abs(iqn-0.5) < 1e-5:
            str += " norm = L_2"
    return str


def run():

    plot_dir = 'plots/'
    import os
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    iqn = 0.5  # Lq norm to measure approximation errors
    conv_dir = 'conv_curves_iqn='+repr(iqn)+'/'
    if not os.path.exists(conv_dir):
        os.makedirs(conv_dir)

    save_plot = True
    test_set = 2

    N_subdomains = 3

    if test_set == 1:
        f_kinds = ['abs', 'runge']  # ['pol1']  # ['runge']
        ncs_values = [2, 10, 20, 40, 80]  # approx nb of cells on each subdomain
        p_values = [1]
        proj_kinds = ['db', 'ps']  # ''ms', 'ps', 'db']
        space_kinds = ['V']

    elif test_set == 2:
        f_kinds = ['abs', 'runge']  #['step', 'pws', 'abs', 'runge']
        ncs_values = [10, 20, 40, 80]   # approx nb of cells on each subdomain
        p_values = [1, 2, 3]  # [1, 2, 3]
        proj_kinds = ['L2', 'db', 'ps', 'ms']  # , 'tL2+db', 'tL2+ps', 'tL2+ms']
        space_kinds = ['V', 'tV']

    else:
        raise ValueError("Error -- Bad test_set value:", test_set)

    color = dict([
        ('L2', 'r'),
        ('db', 'b'),
        ('ps', 'b'),
        ('ms', 'b'),
        ('tL2+db', 'g'),
        ('tL2+ps', 'g'),
        ('tL2+ms', 'g'),
    ])
    linetype = dict([
        ('L2', '-+'),
        ('db', ':'),
        ('ps', '-.'),
        ('ms', '--'),
        ('tL2+db', ':'),
        ('tL2+ps', '-.'),
        ('tL2+ms', '--'),
    ])


    with open('product_errors.txt', 'w') as the_file:
        the_file.write(
            "-- product errors |<f-f_h, x**r>| -- \n\n"
        )

    for p in p_values:

        m = p   ## we could also study m = p-1
        M = p + m + 1

        slope = dict([
            ('step', 0.5),
            ('pws', 0.5),
            ('abs', 1.5),
            ('runge', p+1),
            ('pol3', p+1),
        ])

        cst = dict([
            ('step', 1),
            ('pws', 1),
            ('abs', 0.1),
            ('runge', 100),
            ('pol3', 100),
        ])

        h_list = []
        errors_list = dict(
            [(sk, dict(
                [(pk, dict(
                    [(fk, []) for fk in f_kinds]
                )) for pk in proj_kinds]
            )) for sk in space_kinds]
        )
        slopes_list = dict([(fk, []) for fk in f_kinds])

        for ncs in ncs_values:
            ncs_aux = max(ncs, 3*p+m+1)
            N_cells_sub = int(np.ceil(ncs_aux/M))*M

            spo = SplineProjectionOperators(p, m, N_cells_sub=N_cells_sub, N_subdomains=N_subdomains,
                        use_macro_elem_duals=('ms' in proj_kinds))

            print("ncs = ",ncs," -- N_cells_sub = ", N_cells_sub, " -- N_cells = ", spo.N_cells)
            h = 1./spo.N_cells
            h_list.append(h)

            for fk in f_kinds:

                x_sep = 1/np.pi  # not a regular spline node...

                if fk == 'runge':
                    f = lambda x: 1./((5*( 2*(x-0.5)))**2+1)
                elif fk == 'abs':
                    f = lambda x: abs(x-x_sep)
                elif fk == 'pol1':
                    f = lambda x: x
                elif fk == 'pol3':
                    f = lambda x: (x-0.5)*x**2
                elif fk == 'pws':
                    f = lambda x: int(x > x_sep)* np.sin(3*x*np.pi)
                elif fk == 'step':
                    f = lambda x: int(x > x_sep)
                else:
                    raise ValueError("unkown value for f_kind: ", fk)

                slopes_list[fk].append(cst[fk] * h**slope[fk])

                for sk in space_kinds:

                    for pk in proj_kinds:

                        if is_valid(space_kind=sk, proj_kind=pk):
                            
                            print("## approx f="+fk+" on space "+sk+" with proj="+pk+"...")
                            spo.proj(f, kind=pk, space_kind = sk, localize_quadratures=True)

                            approx_name = file_str(fk, p, space_kind = sk, proj_kind=pk, N_cells=spo.N_cells)
                            if save_plot:
                                filename = plot_dir+approx_name+".png"
                            else:
                                filename = None
                            error = spo.plot_spline(
                                filename=filename,
                                f_ref= f,
                                N_points=400,
                                space_kind=sk,
                                iqnorm=iqn,
                            )
                            errors_list[sk][pk][fk].append(error)

                            if sk == 'V':
                                # print("product errors |<f-f_h, x**r>|: ")
                                with open('product_errors.txt', 'a') as the_file:
                                    the_file.write(
                                        "\n"
                                        +"run: "+title_str(fk,p,space_kind=sk, proj_kind=pk,N_cells=spo.N_cells)
                                        +"\n"
                                    )
                                for r in range(p+1):
                                    f_moment_r = quad(
                                        lambda x: f(x) * x**r, 0, 1,
                                        points=[0.5]
                                    )[0]
                                    fh_moment_r = quad(
                                        lambda x: spo.eval_continuous_spline_splev(x) * x**r, 0, 1,
                                    )[0]
                                    with open('product_errors.txt', 'a') as the_file:
                                        the_file.write(
                                            " (r = "+repr(r)+") -- error = "+repr(abs(f_moment_r-fh_moment_r))+"\n"
                                        )

        print("----- CONVERGENCE CURVES ---- ")
        for sk in space_kinds:
            for fk in f_kinds:
                fig = plt.figure()
                plt.clf()
                #plt.ylim( ymin=1e-7, ymax=1e0 )
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('h')
                plt.title("approx errors, "+ title_str(fk, p, space_kind=sk))
                slope_curve, = plt.plot(h_list, slopes_list[fk], '-', color='m', label='C h^'+repr(slope[fk]))
                #handles=[slope_curve]
                for pk in proj_kinds:
                    if is_valid(space_kind=sk, proj_kind=pk):
                        print("plotting errors for " + title_str(fk, p, space_kind=sk))
                        error_curve, = plt.plot(
                            h_list,
                            errors_list[sk][pk][fk],
                            linetype[pk],
                            color=color[pk],
                            label=pk,
                        )  #  label=r'$||f-Pf||_{L^2}$')
                    #handles.append(error_curve)
                # plt.legend(handles=handles)
                plt.legend(loc='lower right')
                figname = conv_dir+"conv-curve_"+file_str(fk, p, space_kind=sk, iqn=iqn)+".png"
                print("writing convergence curve in "+figname)
                fig.savefig(figname)

                plt.clf()
                plt.close(fig)


####################################################################################
if __name__ == '__main__':

    run()
