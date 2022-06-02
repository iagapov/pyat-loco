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



def correctionType(alpha1,alpha2,alpha3):

    if alpha1 == 1:
       correction_type = "optics correction"
    if alpha2 == 1:
       correction_type = "dispersion correction"
    else:
       correction_type = "optics and dispersion correction"
    print("This code performs: ", correction_type)

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

    opt = at.linopt(lattice, refpts=quad_indexes, get_chrom=True)
    s_pos_q = opt[3].s_pos

    # calculate the first occurrence if elements
    output = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 0)]
    output1 =[elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 1)]


    quad = {'s_pos': s_pos_q, 'Quad_strength': Quad_strength,
            'elements_name': elements_name_q, 'occ': output, 'occ1':output1}
    quads = pd.DataFrame(quad)


    return quads


def getQuadFamilies(quads_info):

    n_list = len(quads_info.s_pos)
    eocc_a = {}
    vals = {}

    for idx in range(n_list):
        ename = quads_info.elements_name[idx]
        par = quads_info.Quad_strength[idx]
        eocc =  quads_info.occ[idx]
        eocc_a[ename] = int(eocc)
        vals[ename, eocc] = float(par)

    return eocc_a, vals

def used_quads_f(ring, used_quads_list, quad_dict):

    #elements_name =
    elements_indexes = []
    quad_dict_ = []
    elements_name =[]
    quads = pd.DataFrame()
    for i in used_quads_list:

        quad_dict_= int(quad_dict[i])
        elements_numbers = quad_dict_
        elements_index = get_refpts(ring, i)
        element_name = ring[elements_index[0]].FamName
        elements_indexes.append(np.squeeze(elements_index))

        df1 = {
                    str(i) + str("=") + str(" ")+ str( quad_dict_)+ str(" ")+ str('quads'): elements_index,
                }
        df2 = pd.concat([pd.DataFrame(v, columns=[k]) for k, v in df1.items()], axis=1)


        elements_name.append(element_name)
        quads = pd.concat([quads, df2], axis=1)


    quads_indexes =[]
    s_pos1 =[]
    for i in used_quads_list:
        element_index = get_refpts(ring, i)
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=element_index)
        s_pos_ = lindata['s_pos']
        s_pos1.append(np.squeeze(s_pos_))
    #s_pos = sorted(numpy.concatenate(s_pos1))
    s_pos= sorted(np.concatenate([x.flatten() for x in s_pos1]))
    return quads, elements_indexes, s_pos



def used_elements_plot(lattice, id_Q):

    elements_indexes = get_refpts(lattice, '*')
    lindata0, tune, chrom, lindata = lattice.linopt(get_chrom=True, refpts=elements_indexes)
    closed_orbitx = lindata['closed_orbit'][:, 0]
    closed_orbity = lindata['closed_orbit'][:, 2]
    s_pos = lindata['s_pos']
    closed_orbit = lindata['closed_orbit']
    beta_x = lindata['beta'][:, 0]
    beta_y = lindata['beta'][:, 1]
    dx = lindata['dispersion'][:, 0]
    dy = lindata['dispersion'][:, 2]

    plt.rc('font', size=15)
    plt.plot(s_pos, closed_orbitx)

    # Label for x-axis
    plt.xlabel("s[m]")
    # Label for y-axis
    plt.ylabel("closed_orbit_x[m]")

    i = 0
    for i in id_Q:
        scatter(i, 0)
    plt.title("used quadrupoles indices")
    plt.show()

    #second_plot

    plt.rc('font', size=15)
    plt.plot(s_pos, beta_x)

    # Label for x-axis
    plt.xlabel("s[m]")
    # Label for y-axis
    plt.ylabel("beta_x[m]")

    i = 0
    for i in id_Q:
        scatter(i, 0)
    plt.title("used quadrupoles indices")
    plt.show()

def func(j, mylist):
    deduped = list(dict.fromkeys(mylist))
    # Slice off all but the part you care about:
    return deduped[::j]


def used_cor(lattice, used_correctors_list):

    elements_name = used_correctors_list
    correctors_indexes = []
    quad_dict = []
    for i in used_correctors_list:
        # print(i)
        corrector_index = get_refpts(lattice, i)
        correctors_indexes.append(np.squeeze(corrector_index))


        j = 0
        s_pos1=[]
        Ind= []
        while (j < len( corrector_index)):
            corrector_indexx = corrector_index[j]
            lindata0, tune, chrom, lindata = lattice.linopt(get_chrom=True, refpts=correctors_indexes)
            s_poss = lindata['s_pos']
            #s_pos.append(s_poss)
            s_pos1.append(np.squeeze(s_poss))
            Ind.append(corrector_indexx)
            df1 = {'Used elements names': elements_name, 'S_pos': correctors_indexes}
            j += 1

    return df1, s_pos1, correctors_indexes

def correctors_phi(ring, ind_c):

    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=ind_c)
    s_pos = lindata['s_pos']
    mux = lindata['mu'][:, 0]
    muy = lindata['mu'][:, 1]

    print("phase advance of the used correctors")

    plt.axes(projection='polar')

    r = 1

    # plotting the circle
    for i in mux:
        plt.polar(i, r, 'g.')

    # display the Polar plot
    plt.show()

    plt.axes(projection='polar')

    r = 1

    # plotting the circle
    for i in muy:
        plt.polar(i, r, 'g.')

    # display the Polar plot
    plt.show()



def getOptics(ring, refpts_):

    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=refpts_)
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
    id = refpts_

    for i in id:
        element_name = ring[i].FamName
        elements_name.append(element_name)

    output = [elements_name[:e].count(v) for e, v in enumerate(elements_name,1)]
    twiss1 = {'s_pos': s_pos,'betax': beta_x,
              'betay': beta_y, 'dx': dx, 'dy': dy}

    twiss = pd.DataFrame(twiss1)



    return twiss

def make_plot(twiss, plot_name):

    from mpl_toolkits.axes_grid1 import host_subplot
    import matplotlib.pyplot as plt
    plt.rc('font', size=15)
    plt.figure(figsize=(8, 3))
    host = host_subplot(111)

    par = host.twinx()

    host.set_xlabel(r'$s$ (m) ', fontsize=14)
    host.set_ylabel(r'$B_x B_y$ (m) ', fontsize=16)
    par.set_ylabel(r'$\eta_x \eta_y$', fontsize=16)

    p1, = host.plot(twiss.s_pos, twiss.betax, 'b', label=r'$\beta_x[m]$')
    p2, = host.plot(twiss.s_pos, twiss.betay, 'g', label=r'$\beta_y[m]$')
    p3, = par.plot(twiss.s_pos, twiss.dx, 'r', label=r'$\eta_x[m]$')
    p4, = par.plot(twiss.s_pos, twiss.dy, label=r'$\eta_y[m]$')
    leg = plt.legend()

    host.yaxis.get_label().set_color(p1.get_color())
    leg.texts[0].set_color(p1.get_color())

    host.yaxis.get_label().set_color(p2.get_color())
    leg.texts[1].set_color(p2.get_color())

    par.yaxis.get_label().set_color(p3.get_color())
    leg.texts[2].set_color(p3.get_color())

    plt.xlim(1000, 1500)

    plt.title(plot_name)

    plt.show()



def eta(dkick, ring):
    bpm_indexes = get_refpts(ring, elements.Monitor)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    Eta_xx = lindata['dispersion'][:, 0]
    Eta_yy = lindata['dispersion'][:, 2]

    return Eta_xx,Eta_yy

def ORM_x_eta(dkick, ring, ind):
    cxx = []
    cxy = []
    c_indexes =[]
    bpm_indexes = get_refpts(ring, elements.Monitor)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    closed_orbitx0 = lindata['closed_orbit'][:, 0]
    closed_orbity0 = lindata['closed_orbit'][:, 2]
    for j in ind:
        ring[j].KickAngle[0] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cxx.append(closed_orbitx-closed_orbitx0)
        cxy.append(closed_orbity-closed_orbity0)

        ring[j].KickAngle[0] = 0.00
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    Eta_xx = lindata['dispersion'][:, 0]
    Eta_yy = lindata['dispersion'][:, 2]
    cxx.append(Eta_xx)
    cxx.append(Eta_yy)
    cxy.append(Eta_xx)
    cxy.append(Eta_yy)

    Cxx = np.squeeze(cxx) / dkick
    Cxy = np.squeeze(cxy) / dkick

    return Cxx, Cxy


def ORM_y_eta(dkick, ring, ind):
    cyy = []
    cyx = []
    c_indexes = []
    #correctors_indexes = get_refpts(ring, used_correctors)
    bpm_indexes = get_refpts(ring, elements.Monitor)
    bpm_indexes = get_refpts(ring, elements.Monitor)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    closed_orbitx0 = lindata['closed_orbit'][:, 0]
    closed_orbity0 = lindata['closed_orbit'][:, 2]

    for j in ind:
        ring[j].KickAngle[1] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cyy.append(closed_orbity-closed_orbity0)
        cyx.append(closed_orbitx-closed_orbitx0)

        ring[j].KickAngle[1] = 0.00

    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    Eta_xx = lindata['dispersion'][:, 0]
    Eta_yy = lindata['dispersion'][:, 2]
    cyy.append(Eta_xx)
    cyy.append(Eta_yy)
    cyx.append(Eta_xx)
    cyx.append(Eta_yy)

    Cyy = np.squeeze(cyy) / dkick
    Cyx = np.squeeze(cyx) / dkick

    return Cyy, Cyx


def plotORM(orm):
    plt.rc('font', size=15)
    plt.figure()
    imshow(orm)
    plt.show()



def generatingQuadsResponse(ring, Cxx, Cyy,Cxy, Cyx , used_quads,indc):
    # %%time

    quads_info = quad_info(ring)
    quad_dict, quad_vals = getQuadFamilies(quads_info)
    quads = [k for k in quad_dict.keys()]
    dk = 0.0001
    qxx = []
    qxy = []
    qyy = []
    qyx = []
    quad_names = quads

    for qname in quad_names:
        if qname in used_quads:
            print('generating response to {}'.format(qname))
            t0 = time.time()
            nq = quad_dict[qname] + 1
            #nqq=nq/2
            #for i in range(0, int(nq)):
            Qxx, Qxy, Qyy, Qyx = computeOpticsD(ring, qname, dk, quad_vals, indc)
            qxx.append(Qxx)
            qxy.append(Qxy)
            qyy.append(Qyy)
            qyx.append(Qyx)
            t1 = time.time()
            print(f"Execution time: {t1 - t0} sec")

    C0x = Cxx
    C0y = Cyy
    C0xy = Cxy
    C0yx = Cyx

    dCx = []
    dCy = []
    dCxy = []
    dCyx = []
    quad_names = quads

    i = 0
    while (i < len(qxx)):
        C1x = qxx[i]
        C1y = qyy[i]
        C1xy = qxy[i]
        C1yx = qyx[i]
        dcxx = ((C1x - C0x) / dk)
        dcyy = ((C1y - C0y) / dk)

        dCxy.append((C1xy - C0xy) / dk)
        dCyx.append((C1yx - C0yx) / dk)

        dCx.append(dcxx)
        dCy.append(dcyy)

        i += 1

    return C0x, C0y, C0xy, C0yx, dCx, dCy, dCxy,dCyx



def computeOpticsD(ring, qname, dk, quad_vals,indc):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    quad_indexes = get_refpts(ring, qname)
    for i in quad_indexes:
        a = ring[i].K
        #print('a',a)
        ring[i].K = a + dk
        b = ring[i].K
        #print('b',b)
    qxx, qxy = ORM_x_eta(dk, ring, indc)
    qyy, qyx = ORM_y_eta(dk, ring, indc)

    for i in quad_indexes:
        ring[i].K = a
        c = ring[i].K
        #print('c',c)
    return  qxx, qxy, qyy, qyx


def simulateErrorQFQD(lattice, errorQ, tiltQ, shiftQx,shiftQy, debug=False):

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


        if lattice[quad_indexes[i]].FamName == 'QF':
            lattice[quad_indexes[i]].K *= (1 + errorQ * np.random.randn())
            at.tilt_elem(lattice[quad_indexes[i]], tiltQ * np.random.randn(), relative=False)
            at.shift_elem(lattice[quad_indexes[i]], shiftQx * np.random.randn(), shiftQy * np.random.randn(),
                          relative=False)

        if lattice[quad_indexes[i]].FamName == 'QD':
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


def simulateError(lattice, errorQ, tiltQ, shiftQx,shiftQy, debug=False):

    ## Insert the errors
    np.random.seed(1)

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

        if debug: print(f"sumulateError: qaudrupole name {lattice[quad_indexes[i]].FamName} , #{i} out of {len(quad_indexes)}")
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


def getBetaBeat(twiss, twiss_error):

    print("getBetaBeat bx and by: ")

    bxx = np.array((twiss_error.betax - twiss.betax) / twiss.betax)
    byy = np.array((twiss_error.betay - twiss.betay) / twiss.betay)

    bx = np.sqrt(np.mean(bxx ** 2))
    by = np.sqrt(np.mean(byy ** 2))

    print("Simulated beta beat, x:" + str(bx * 100) + "%   y: " + str(by* 100) + "%")


def defineMatrices_w_eta(W, alpha1, alpha2, alpha3, C0x, C0y, C0xy, C0yx, Cxx_err, Cyy_err, Cxy_err, Cyx_err, dCx, dCy,
                         dCxy, dCyx):
    Nk = len(dCx)  # number of free parameters
    Nm = len(dCx)  # number of measurements
    print('NK:', Nk)
    print('Nm:', Nm)

    Ax = np.zeros([Nk, Nk])
    Ay = np.zeros([Nk, Nk])
    Axy = np.zeros([Nk, Nk])
    Ayx = np.zeros([Nk, Nk])

    A = np.zeros([4 * Nk, Nk])

    ##

    Bx = np.zeros([Nk, 1])
    By = np.zeros([Nk, 1])
    Bxy = np.zeros([Nk, 1])
    Byx = np.zeros([Nk, 1])

    B = np.zeros([4 * Nk, 1])
    ##

    Dx = (Cxx_err[:, :] - C0x[:, :])  # - error_variance)  ### dk ?
    Dy = (Cyy_err[:, :] - C0y[:, :])
    Dxy = (Cxy_err[:, :] - C0xy[:, :])
    Dyx = (Cyx_err[:, :] - C0yx[:, :])

    ##

    for i in range(Nk):  ## i represents each quad
        # print('done A:', 100.* i ,'%')
        for j in range(Nk):
            Ax[i, j] = np.sum(np.dot(np.dot(dCx[i][0: -2, :], W * alpha1), dCx[j][0: -2, :].T)) + np.sum(
                np.dot(np.dot(dCx[i][-2::, :], W * alpha2), dCx[j][-2::, :].T)) + np.sum(
                np.dot(np.dot(dCx[i], W * alpha3), dCx[j].T))
            Ay[i, j] = np.sum(np.dot(np.dot(dCy[i][0: -2, :], W * alpha1), dCy[j][0: -2, :].T)) + np.sum(
                np.dot(np.dot(dCy[i][-2::, :], W * alpha2), dCy[j][-2::, :].T)) + np.sum(
                np.dot(np.dot(dCy[i], W * alpha3), dCy[j].T))
            Axy[i, j] = np.sum(np.dot(np.dot(dCxy[i][0: -2, :], W * alpha1), dCxy[j][0: -2, :].T)) + np.sum(
                np.dot(np.dot(dCxy[i][-2::, :], W * alpha2), dCxy[j][-2::, :].T)) + np.sum(
                np.dot(np.dot(dCxy[i], W * alpha3), dCxy[j].T))
            Ayx[i, j] = np.sum(np.dot(np.dot(dCyx[i][0: -2, :], W * alpha1), dCyx[j][0: -2, :].T)) + np.sum(
                np.dot(np.dot(dCyx[i][-2::, :], W * alpha2), dCyx[j][-2::, :].T)) + np.sum(
                np.dot(np.dot(dCyx[i], W * alpha3), dCyx[j].T))
        A[i, :] = Ax[i, :]
        A[i + Nk, :] = Ay[i, :]
        A[i + 2 * Nk, :] = Axy[i, :]
        A[i + 3 * Nk, :] = Ayx[i, :]

    ##

    for i in range(Nk):
        Bx[i] = np.sum(np.dot(np.dot(dCx[i][0: -2, :], W * alpha1), Dx[0: -2, :].T)) + np.sum(
            np.dot(np.dot(dCx[i][-2::, :], W * alpha2), Dx[-2::, :].T)) + np.sum(
            np.dot(np.dot(dCx[i], W * alpha3), Dx.T))
        By[i] = np.sum(np.dot(np.dot(dCy[i][0: -2, :], W * alpha1), Dy[0: -2, :].T)) + np.sum(
            np.dot(np.dot(dCy[i][-2::, :], W * alpha2), Dy[-2::, :].T)) + np.sum(
            np.dot(np.dot(dCy[i], W * alpha3), Dy.T))
        Bxy[i] = np.sum(np.dot(np.dot(dCxy[i][0: -2, :], W * alpha1), Dxy[0: -2, :].T)) + np.sum(
            np.dot(np.dot(dCxy[i][-2::, :], W * alpha2), Dxy[-2::, :].T)) + np.sum(
            np.dot(np.dot(dCxy[i], W * alpha3), Dxy.T))
        Byx[i] = np.sum(np.dot(np.dot(dCyx[i][0: -2, :], W * alpha1), Dyx[0: -2, :].T)) + np.sum(
            np.dot(np.dot(dCyx[i][-2::, :], W * alpha2), Dyx[-2::, :].T)) + np.sum(
            np.dot(np.dot(dCyx[i], W * alpha3), Dyx.T))

        B[i] = Bx[i]
        B[i + Nk] = By[i]
        B[i + 2 * Nk] = Bxy[i]
        B[i + 3 * Nk] = Byx[i]

    return A, B,

def getInverse(A, B,Nk, sCut):

    u, s, v = np.linalg.svd(A, full_matrices=True)

    smat = 0.0 * A
    si = s ** -1
    n_sv = sCut
    si[n_sv:] *= 0.0

    print("number of singular values {}".format(len(si)))
    smat[:Nk, :Nk] = np.diag(si)

    print('A' + str(A.shape), 'B' + str(B.shape), 'U' + str(u.shape), 'smat' + str(smat.shape), 'v' + str(v.shape))
    plt.rc('font', size=15)
    plt.plot(np.log(s), 'd--')
    plt.title('singular value')
    plt.show()

    plt.plot(si, 'd--')
    plt.title('singular value inverse')
    plt.show()

    Ai = np.dot(v.transpose(), np.dot(smat.transpose(), u.transpose()))

    ###
    #r = (np.dot(Ai, B)).reshape(-1)
    r = (np.dot(Ai, B)).reshape(-1)
    plot(r, 'd')
    plt.show()

    # error
    e = np.dot(A, r).reshape(-1) - B.reshape(-1)
    plt.plot(e)
    plt.show()
    plt.plot(B)
    plt.show()

    return Ai, r, e

def setCorrection(ring, quads_info_error,quad_names, r , quads_info,n_list, used_quads):

    quad_dict, quad_vals = getQuadFamilies(quads_info_error)
    n_list = len(quads_info_error.s_pos)
    # print(n_list)

    quad_names = quad_names
    iq = 0
    frac = 1.0
    cor_dict = {}
    DK = []

    for qname in quad_names:
         if qname in used_quads:
            cor_dict[qname] = -r[iq] * frac
            iq += 1
    print("define correction : Done")


    quads_indexes = get_refpts(ring, elements.Quadrupole)
    for qname in quads_indexes:
        if ring[qname].FamName in used_quads:
            dk1 = cor_dict[ring[qname].FamName]
            DK.append(dk1)

        else:
            DK.append(0)



    quads_indexes = get_refpts(ring, elements.Quadrupole)
    i = 0
    while (i < len(quads_indexes)):


        ring[quads_indexes[i]].K += DK[i]


        i += 1

    print("set correction : Done")

def compare_orm(Cxy, Cxy_err, Cxy_corr, no):
    # plot the 3 sets
    plt.rc('font', size=15)
    plt.plot(Cxy[no], label='C')
    plt.plot(Cxy_err[no], label='C_err')
    plt.plot(Cxy_corr[no], label='C_corr')

    # call with no parameters
    plt.legend()

    plt.show()

def compare_drm(Cxy, Cxy_err, Cxy_corr):
    plt.rc('font', size=15)
    # plot the 3 sets
    plt.plot(Cxy, label='$\eta$')
    plt.plot(Cxy_err, label='$\eta_{err}$')
    plt.plot(Cxy_corr, label='$\eta_{corr}$')

    # call with no parameters
    plt.legend()

    plt.show()