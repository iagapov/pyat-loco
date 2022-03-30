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

    plt.figure()
    plt.plot(twiss.s_pos, bxi)
    plt.xlabel('s(m)')
    plt.ylabel(r'$\Delta \beta_x / \beta_x$')
    plt.show()

    plt.figure()
    plt.plot(twiss.s_pos, byi)
    plt.xlabel('s(m)')
    plt.ylabel(r'$\Delta \beta_y / \beta_y$')
    plt.show()


    print("Simulated beta beat, x:" + str(bx * 100) + "%   y: " + str(by* 100) + "%")








def getOptics(ring, refpts, optics, error):


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

    if (optics == 'closed_orbit'):
        p1, = plt.plot(s_pos, closed_orbitx)
        p2, = plt.plot(s_pos, closed_orbity)
        plt.legend([p1, p2], [r'$x$', r'$y$'])

    if (optics == 'beta'):
        # plotting the bet + the dispersion
        fig, ax = plt.subplots()
        ax.plot(s_pos, beta_x, color="blue")
        ax.plot(s_pos, beta_y, color="red")
        ax.set_xlabel(r's(m)', color="blue", fontsize=15)
        ax.set_ylabel(r'$\beta_x$', color="blue", fontsize=15)
        ax.set_ylim([3, 11])
        ax.set_ylabel(r'$\beta_x \beta_y$', color="red", fontsize=15)

        # twin object for two different y-axis on the sample plot
        ax2 = ax.twinx()
        # make a plot with different y-axis using second axis object
        ax2.plot(s_pos, dx, color="green")
        ax2.set_ylabel(r'$\eta_x$', color="green", fontsize=14)
        ax2.set_ylim([3.4, 4.6])
        plt.show()

    if (optics == 'dispersion'):
        plt.figure()
        p1, = plt.plot(s_pos, dx)
        p2, = plt.plot(s_pos, dy)
        plt.legend([p1, p2], ['Dx', 'Dy'])

    elements_indexes = get_refpts(ring, '*')


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

    twiss1 = {'s_pos': s_pos, 'closed_orbitx': closed_orbitx, 'closed_orbity': closed_orbity, 'betax': beta_x, 'betay': beta_y,'dx':dx, 'dy':dy
              , 'elements_strength': elements_strength,
              'elements_name': elements_name, 'occ':output, 'elements_type':elements_type}
    twiss2 = {'s_pos': s_pos,'betax': beta_x,
              'betay': beta_y, 'dx': dx, 'dy': dy}

    twiss = pd.DataFrame(twiss1)
    if error == 'True':
       twiss_error = pd.DataFrame(twiss1)
       twiss_error.to_csv('C:/Users/musa/pyat-loco-1/fodo_loco/mydata/twiss_error.csv')
    twiss = pd.DataFrame(twiss2)
    twiss.to_csv('C:/Users/musa/pyat-loco-1/fodo_loco/mydata/twiss.csv')


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
            #at.shift_elem(lattice[quad_indexes[i]], deltax=shiftQF, deltaz=0, relative=False)
            #at.tilt_elem(lattice[quad_indexes[i]], tiltQF  , relative=False)

            i += 1

        elif (lattice[quad_indexes[i]].FamName == 'QD'):
            lattice[quad_indexes[i]].K *= (1 + errorQD * random())
            #at.shift_elem(lattice[quad_indexes[i]], deltax=shiftQD, deltaz=0, relative=False)
            #at.tilt_elem(lattice[quad_indexes[i]],tiltQD  , relative=False)

            i += 1

        elif (lattice[quad_indexes[i]].FamName == 'QS'):
            lattice[quad_indexes[i]].K *= 0.0
            # at.shift_elem(lattice[quad_indexes[i]], deltax=shiftQD, deltaz=0, relative=False)
            # at.tilt_elem(lattice[quad_indexes[i]], tiltQD, relative=False)

            i += 1

    for j in quad_indexes:
        quad_strength_err = lattice[j].K
        Quad_strength_err.append(quad_strength_err)

    elements_indexes = get_refpts(lattice, elements.Quadrupole)
    opt = at.linopt(lattice, refpts=elements_indexes, get_chrom=True)
    s_pos_q = opt[3].s_pos

    #output = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 0)]
    output1 = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 1)]
    #occurrences = {}
    #elements_occurence = []
    #for i in elements_name_q:
    #    if i in occurrences:
    #        occurrences[i] += 1
    #    else:
    #        occurrences[i] = 1
    #    elements_occurence.append(occurrences[i])


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
        #etype = Quads_info.elements_type[idx]
        ename = Quads.elements_name[idx]
        par = Quads.Quad_strength[idx]
        eocc =  Quads.occ[idx]
        eocc_a[ename] = int(eocc)
        vals[ename, eocc] = float(par)

    return eocc_a, vals
