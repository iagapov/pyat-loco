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


def getBetaBeat(twiss, twiss_error):
    print("getBetaBeat bx and by: ")

    bxi =[]
    for i in range(len(twiss.betax)):
        bxx = (twiss_error.betax[i] - twiss.betax[i]) / twiss.betax[i]
        bxi.append(bxx)

    byi =[]
    for i in range(len(twiss.betay)):
        byy = (twiss_error.betay[i] - twiss.betay[i]) / twiss.betay[i]
        byi.append(byy)

    bx = np.std((twiss_error.betax - twiss.betax) / twiss.betax)
    by = np.std((twiss_error.betay - twiss.betay) / twiss.betay)
    print("Simulated beta beat, x:" + str(bx * 100) + "%   y: " + str(by* 100) + "%")

def make_plot(twiss, plot_name):
    from mpl_toolkits.axes_grid1 import host_subplot
    import matplotlib.pyplot as plt

    host = host_subplot(111)

    par = host.twinx()

    host.set_xlabel("s_pos")
    host.set_ylabel(r'$\beta_x$')
    host.set_ylabel(r'$\beta_y$')
    par.set_ylabel("dx")

    p1, = host.plot(twiss.s_pos, twiss.betax, label=r'$\beta_x$')
    p2, = host.plot(twiss.s_pos, twiss.betay, label=r'$\beta_y$')
    p3, = par.plot(twiss.s_pos, twiss.dx, label=r'$\eta_x$')

    leg = plt.legend()

    host.yaxis.get_label().set_color(p1.get_color())
    leg.texts[0].set_color(p1.get_color())

    host.yaxis.get_label().set_color(p2.get_color())
    leg.texts[1].set_color(p2.get_color())

    par.yaxis.get_label().set_color(p3.get_color())
    leg.texts[2].set_color(p3.get_color())

    plt.title(plot_name)

    plt.show()





def getOptics(ring, refpts):


    elements_indexes = get_refpts(ring, refpts)

    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=elements_indexes)
    closed_orbitx = lindata['closed_orbit'][:, 0]
    closed_orbity = lindata['closed_orbit'][:, 1]
    s_pos = lindata['s_pos']
    closed_orbit = lindata['closed_orbit']
    beta_x= lindata['beta'][:, 0]
    beta_y= lindata['beta'][:, 1]
    dx = lindata['dispersion'][:, 0]
    dy = lindata['dispersion'][:, 1]

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
        if (ring[i].FamName == 'QF' or ring[i].FamName == 'QD'):

            elements_Strength = ring[i].K
            elements_strength.append(elements_Strength)
            element_type = "elements.Quadrupole"
            elements_type.append(element_type)


            i += 1

        else:
            elements_Strength = 0
            elements_strength.append(elements_Strength)
            element_type = 0
            elements_type.append(element_type)
            i += 1

    output = [elements_name[:e].count(v) for e, v in enumerate(elements_name,1)]
    twiss1 = {'s_pos': s_pos,'betax': beta_x,
              'betay': beta_y, 'dx': dx, 'dy': dy}

    twiss = pd.DataFrame(twiss1)



    return twiss


def ORM_x(dkick, ring):
    cxx = []
    cxy = []
    correctors_indexes = get_refpts(ring, elements.Corrector)
    bpm_indexes = get_refpts(ring, elements.Monitor)

    for j in range(len(correctors_indexes)):
        ring[correctors_indexes[j]].KickAngle[0] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cxx.append(closed_orbitx)
        cxy.append(closed_orbity)

        ring[correctors_indexes[j]].KickAngle[0] = 0.00

    Cxx = np.squeeze(cxx) / dkick
    Cxy = np.squeeze(cxy) / dkick

    return Cxx, Cxy


def ORM_y(dkick, ring):
    cyy = []
    cyx = []
    correctors_indexes = get_refpts(ring, elements.Corrector)
    bpm_indexes = get_refpts(ring, elements.Monitor)

    for j in range(len(correctors_indexes)):
        ring[correctors_indexes[j]].KickAngle[1] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cyy.append(closed_orbity)
        cyx.append(closed_orbitx)

        ring[correctors_indexes[j]].KickAngle[1] = 0.00

    Cyy = np.squeeze(cyy) / dkick
    Cyx = np.squeeze(cyx) / dkick

    return Cyy, Cyx








def ORM_x_q(dkick, ring):
    cxx = []
    cxy = []
    correctors_indexes = get_refpts(ring, elements.Quadrupole)
    bpm_indexes = get_refpts(ring, elements.Monitor)

    for j in range(len(correctors_indexes)):
        ring[correctors_indexes[j]].PolynomB[0] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cxx.append(closed_orbitx)
        cxy.append(closed_orbity)

        ring[correctors_indexes[j]].PolynomB[0] = 0.00

    Cxx = np.squeeze(cxx) / dkick
    Cxy = np.squeeze(cxy) / dkick

    return Cxx, Cxy


def ORM_y_q(dkick, ring):
    cyy = []
    cyx = []
    correctors_indexes = get_refpts(ring, elements.Quadrupole)
    bpm_indexes = get_refpts(ring, elements.Monitor)

    for j in range(len(correctors_indexes)):
        ring[correctors_indexes[j]].PolynomB[1] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cyy.append(closed_orbity)
        cyx.append(closed_orbitx)

        ring[correctors_indexes[j]].PolynomB[1] = 0.00

    Cyy = np.squeeze(cyy) / dkick
    Cyx = np.squeeze(cyx) / dkick

    return Cyy, Cyx



















def ORM_x_error(dkick, ring):
    cxx = []
    cxy = []


    correctors_indexes = get_refpts(ring, elements.Corrector)
    bpm_indexes = get_refpts(ring, elements.Monitor)


    for j in range(len(correctors_indexes)):
        ring[correctors_indexes[j]].KickAngle[0] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        s_pos = lindata['s_pos']
        closed_orbit = lindata['closed_orbit']
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cxx.append(closed_orbitx)
        cxy.append(closed_orbity)


        ring[correctors_indexes[j]].KickAngle[0] = 0.00

    Cxx = np.squeeze(cxx) / dkick
    Cxy = np.squeeze(cxy) / dkick


    return Cxx, Cxy





def ORM_y_error(dkick, ring):
    cyy = []
    cyx = []


    correctors_indexes = get_refpts(ring, elements.Corrector)
    bpm_indexes = get_refpts(ring, elements.Monitor)


    for j in range(len(correctors_indexes)):
        ring[correctors_indexes[j]].KickAngle[1] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        s_pos = lindata['s_pos']
        closed_orbit = lindata['closed_orbit']
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cyy.append(closed_orbity)
        cyx.append(closed_orbitx)


        ring[correctors_indexes[j]].KickAngle[1] = 0.00

    Cyy = np.squeeze(cyy) / dkick
    Cyx = np.squeeze(cyx) / dkick


    return Cyy, Cyx


def quad_info(lattice):

    Quad_strength = []
    Quad_strength_err = []
    elements_name_q = []

    quad_indexes = get_refpts(lattice, elements.Quadrupole)

    for i in quad_indexes:
        quad_strength = lattice[i].K
        Quad_strength.append(quad_strength)
        element_name = lattice[i].FamName
        elements_name_q.append(element_name)

    elements_indexes = get_refpts(lattice, elements.Quadrupole)
    opt = at.linopt(lattice, refpts=elements_indexes, get_chrom=True)
    s_pos_q = opt[3].s_pos

    output = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 0)]
    output1 = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 1)]


    quad = {'s_pos': s_pos_q, 'Quad_strength': Quad_strength,
            'elements_name': elements_name_q, 'occ': output, 'occ1':output1}
    quads = pd.DataFrame(quad)
    quads.to_csv("C:/Users/musa/pyat-loco-1/fodo_loco/mydata/quad_info.csv")

    print('Done...')

    return quads


def simulateError(lattice, errorQF, errorQD):

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

        if (lattice[quad_indexes[i]].FamName == 'QF'):
            lattice[quad_indexes[i]].K *= (1 + errorQF * random())

            i += 1

        elif (lattice[quad_indexes[i]].FamName == 'QD'):
            lattice[quad_indexes[i]].K *= (1 + errorQD * random())

        i += 1

        #elif (lattice[quad_indexes[i]].FamName == 'QS'):
        #    lattice[quad_indexes[i]].K = 0.0

        #    i += 1

    for j in quad_indexes:
        quad_strength_err = lattice[j].K
        Quad_strength_err.append(quad_strength_err)

    elements_indexes = get_refpts(lattice, elements.Quadrupole)
    opt = at.linopt(lattice, refpts=elements_indexes, get_chrom=True)
    s_pos_q = opt[3].s_pos

    output1 = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 1)]


    quad = {'s_pos': s_pos_q, 'Quad_strength': Quad_strength_err,
            'elements_name': elements_name_q, 'occ': output1}
    quads = pd.DataFrame(quad)
    quads.to_csv("C:/Users/musa/pyat-loco-1/fodo_loco/mydata/quad_info_error.csv")


    print('Done...')

    return quads



def simulateError_tilt(lattice, errorQF, errorQD, tiltQF, tiltQD):

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

        if (lattice[quad_indexes[i]].FamName == 'QF'):
            lattice[quad_indexes[i]].K *= (1 + errorQF * random())
            at.tilt_elem(lattice[quad_indexes[i]], tiltQF, relative=False)

            i += 1

        elif (lattice[quad_indexes[i]].FamName == 'QD'):
            lattice[quad_indexes[i]].K *= (1 + errorQD * random())
            at.tilt_elem(lattice[quad_indexes[i]], tiltQD, relative=False)

        i += 1

       # elif (lattice[quad_indexes[i]].FamName == 'QS'):
       #     lattice[quad_indexes[i]].K *= (1 + errorQS * random())
       #     at.tilt_elem(lattice[quad_indexes[i]], tiltQS, relative=False)


        #    i += 1

    for j in quad_indexes:
        quad_strength_err = lattice[j].K
        Quad_strength_err.append(quad_strength_err)

    elements_indexes = get_refpts(lattice, elements.Quadrupole)
    opt = at.linopt(lattice, refpts=elements_indexes, get_chrom=True)
    s_pos_q = opt[3].s_pos

    output1 = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 1)]


    quad = {'s_pos': s_pos_q, 'Quad_strength': Quad_strength_err,
            'elements_name': elements_name_q, 'occ': output1}
    quads = pd.DataFrame(quad)
    quads.to_csv("C:/Users/musa/pyat-loco-1/fodo_loco/mydata/quad_info_error.csv")


    print('Done...')

    return quads


def getQuadFamilies(Quads):

    n_list = len(Quads.s_pos)
    eocc_a = {}
    vals = {}

    for idx in range(n_list):
        ename = Quads.elements_name[idx]
        par = Quads.Quad_strength[idx]
        eocc =  Quads.occ[idx]
        eocc_a[ename] = int(eocc)
        vals[ename, eocc] = float(par)

    return eocc_a, vals

def computeOpticsD(ring, qname, i, dk, quad_vals):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    quad_indexes = get_refpts(ring, qname)

    ring[quad_indexes[i]].K = quad_vals[qname,i] + dk

    qxx, qxy = ORM_x(dk, ring)
    qyy, qyx = ORM_y(dk, ring)

    ring[quad_indexes[i]].K = quad_vals[qname,i]


    return  qxx, qxy, qyy, qyx