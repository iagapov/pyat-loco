import numpy as np
from at import *
from at.load import load_mat
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
import at.plot
import numpy as np
from pylab import *
import pandas as pd
import csv
from random import random

def plot_closedOrbit_all(ring, refpts):

    elements_indexes = get_refpts(ring, refpts)

    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=elements_indexes)
    closed_orbitx = lindata['closed_orbit'][:, 0]
    closed_orbity = lindata['closed_orbit'][:, 2]
    s_pos = lindata['s_pos']
    closed_orbit = lindata['closed_orbit']
    beta_x= lindata['beta'][:, 0]
    beta_y= lindata['beta'][:, 1]
    dx = lindata['dispersion'][:, 0]
    dy = lindata['dispersion'][:, 2]
    plt.rc('font', size=13)

    plt.plot(s_pos, closed_orbitx/1.e-06)

    # Label for x-axis
    plt.xlabel("s_pos[m]", )

    # Label for y-axis
    plt.ylabel("closed_orbit x [m]*1.e-06")

    # for display

    i = 0
    S_pos2 = []
    plt.title("Closed orbit x")
    plt.show()


    plt.plot(s_pos, closed_orbity/1.e-06)

    # Label for x-axis
    plt.xlabel("s_pos[m]")

    # Label for y-axis
    plt.ylabel("closed_orbit y [m]*1.e-06")

    # for display

    i = 0
    S_pos2 = []

    plt.title("Closed orbit y")
    plt.show()

    n = len(closed_orbitx)
    rmsx =np.sqrt(np.mean(closed_orbitx ** 2))
    rmsy =np.sqrt(np.mean(closed_orbity ** 2))

    return rmsx, rmsy

def ResponseMatrix_x(dkick, ring):

    cxx = []
    correctors_indexes = get_refpts(ring, elements.Corrector)
    bpm_indexes = get_refpts(ring, elements.Monitor)

    for j in range(len(correctors_indexes)):
        ring[correctors_indexes[j]].KickAngle[0] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cxx.append(closed_orbitx)

        ring[correctors_indexes[j]].KickAngle[0] = 0.00

    R = np.squeeze(cxx) / dkick

    return R



def ResponseMatrix_y(dkick, ring):

    cyy = []
    correctors_indexes = get_refpts(ring, elements.Corrector)
    bpm_indexes = get_refpts(ring, elements.Monitor)

    for j in range(len(correctors_indexes)):
        ring[correctors_indexes[j]].KickAngle[1] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cyy.append(closed_orbity)

        ring[correctors_indexes[j]].KickAngle[1] = 0.00

    R = np.squeeze(cyy) / dkick

    return R


def simulateError(lattice, errorQ, tiltQ, shiftQx, shiftQy, debug=False):
    ## Insert the errors

    print('simulating perturbed machine...')

    Quad_strength = []
    Quad_strength_err = []
    elements_name_q = []

    quad_indexes = get_refpts(lattice, elements.Quadrupole)

    for i in quad_indexes:
        quad_strength = lattice[i].K
        Quad_strength.append(quad_strength)
        element_name = lattice[i].FamName
        elements_name_q.append(element_name)

    i = 0
    while (i < len(quad_indexes)):

        if debug: print(
            f"sumulateError: qaudrupole name {lattice[quad_indexes[i]].FamName} , #{i} out of {len(quad_indexes)}")
        lattice[quad_indexes[i]].K *= (1 + errorQ * np.random.randn())
        at.tilt_elem(lattice[quad_indexes[i]], tiltQ * np.random.randn(), relative=False)
        at.shift_elem(lattice[quad_indexes[i]], shiftQx * np.random.randn(), shiftQy * np.random.randn(),
                      relative=False)

        i += 1

    if debug: print(f"done")

    for j in quad_indexes:
        quad_strength_err = lattice[j].K
        Quad_strength_err.append(quad_strength_err)

    opt = at.linopt(lattice, refpts=quad_indexes, get_chrom=True)
    s_pos_q = opt[3].s_pos

    output1 = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 1)]

    quad = {'s_pos': s_pos_q, 'Quad_strength': Quad_strength_err,
            'elements_name': elements_name_q, 'occ': output1}
    quads = pd.DataFrame(quad)

    print('Done...')

    return quads



def show_lindata(lat, refpts):
    '''
    calculate and plot linear optical functions
    '''
    lindata0, tune, chrom, lindata = lat.linopt(get_chrom=True, refpts=refpts)
    print(lindata.dtype.names)
    s_pos = lindata['s_pos']
    closed_orbit = lindata['closed_orbit']
    closed_orbitx = lindata['closed_orbit'][:, 0]
    closed_orbity = lindata['closed_orbit'][:, 1]
    betax = lindata['beta'][:, 0]
    betay = lindata['beta'][:, 1]
    dispersion = lindata['dispersion'][:, 0]
    f = plt.figure(figsize=(10, 3))
    ax = f.add_subplot(121)
    ax2 = f.add_subplot(122)
    ax.plot(s_pos, closed_orbitx)
    ax2.plot(s_pos, closed_orbity)
    ax.title.set_text("Closed_orbit_x [m]")
    ax2.title.set_text("Closed_orbit_y [m]")
    # plt.title("Closed orbit")

    f = plt.figure(figsize=(10, 3))
    ax = f.add_subplot(121)
    ax2 = f.add_subplot(122)
    ax.plot(s_pos, betax)
    ax2.plot(s_pos, betay)
    ax.title.set_text("Beta_x [m]")
    ax2.title.set_text("Beta_y [m]")
    # plt.title("Beta")

    plt.show()

def getOptics(ring, refpts):


    elements_indexes = get_refpts(ring, refpts)

    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=elements_indexes)
    closed_orbitx = lindata['closed_orbit'][:, 0]
    closed_orbity = lindata['closed_orbit'][:, 2]
    s_pos = lindata['s_pos']
    closed_orbit = lindata['closed_orbit']
    beta_x= lindata['beta'][:, 0]
    beta_y= lindata['beta'][:, 1]
    dx = lindata['dispersion'][:, 0]
    dy = lindata['dispersion'][:, 2]

    print("preparing twiss ..")
    print(f"Tunes={ring.get_tune()}")
    print(f"Chrom={ring.get_chrom()}")


    elements_name = []
    elements_strength = []
    elements_type =[]
    i = 0
    while (i < len(elements_indexes)):
        element_name = ring[i].FamName
        elements_name.append(element_name)

        i +=1
    output = [elements_name[:e].count(v) for e, v in enumerate(elements_name,1)]
    twiss1 = {'s_pos': s_pos,'betax': beta_x,
              'betay': beta_y, 'dx': dx, 'dy': dy}

    twiss = pd.DataFrame(twiss1)



    return twiss

def make_plot(twiss, plot_name):
    from mpl_toolkits.axes_grid1 import host_subplot
    import matplotlib.pyplot as plt
    plt.rc('font', size=15)

    host = host_subplot(111)

    par = host.twinx()

    host.set_xlabel("s_pos")
    host.set_ylabel(r'$\beta_x$')
    host.set_ylabel(r'$\beta_y$')
    par.set_ylabel("dx")

    p1, = host.plot(twiss.s_pos, twiss.betax, label=r'$\beta_x$')
    p2, = host.plot(twiss.s_pos, twiss.betay, label=r'$\beta_y$')
    p3, = par.plot(twiss.s_pos, twiss.dx, label=r'$\eta_x$')
    p4, = par.plot(twiss.s_pos, twiss.dy, label=r'$\eta_y$')
    leg = plt.legend()

    host.yaxis.get_label().set_color(p1.get_color())
    leg.texts[0].set_color(p1.get_color())

    host.yaxis.get_label().set_color(p2.get_color())
    leg.texts[1].set_color(p2.get_color())

    par.yaxis.get_label().set_color(p3.get_color())
    leg.texts[2].set_color(p3.get_color())

    plt.title(plot_name)

    plt.show()

