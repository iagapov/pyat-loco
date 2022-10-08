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

def getOptics(ring, refpts, bpm_random_noise):


    elements_indexes = get_refpts(ring, refpts)

    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=elements_indexes)
    closed_orbitx, closed_orbity = closed_orbit_bpm(ring, bpm_random_noise)
    s_pos = lindata['s_pos']
    closed_orbit = lindata['closed_orbit']
    beta_x= lindata['beta'][:, 0]
    beta_y= lindata['beta'][:, 1]
    dx = lindata['dispersion'][:, 0]
    dy = lindata['dispersion'][:, 2]

    twiss1 = {'s_pos': s_pos,'betax': beta_x,
              'betay': beta_y, 'dx': dx, 'dy': dy}

    twiss = pd.DataFrame(twiss1)

    return twiss

def closed_orbit_bpm(ring, bpm_random_noise):

    bpm_indexes = get_refpts(ring, elements.Monitor)

    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    closed_orbitx1 = lindata['closed_orbit'][:, 0]
    closed_orbity1 = lindata['closed_orbit'][:, 2]

    Eta_xx1 = lindata['dispersion'][:, 0]
    Eta_yy1 = lindata['dispersion'][:, 2]

    bpm_noise = bpm_random_noise * np.random.randn(len(bpm_indexes))

    closed_orbitx = closed_orbitx1 + bpm_noise
    closed_orbity = closed_orbity1 + bpm_noise


    return closed_orbitx, closed_orbity

def make_plot(twiss, plot_name,file_name, min, max, ymin, ymax):
    import matplotlib.pyplot as plt
    plt.plot(twiss.s_pos, twiss.betax, label=r'$\beta_x$')
    plt.plot(twiss.s_pos, twiss.betay, label=r'$\beta_y$')


    plt.xlabel("s[m]")
    plt.ylabel(r'$\beta_x \beta_y$')

    plt.xlim(min, max)
    plt.ylim(ymin, ymax)

    plt.title(plot_name)
    plt.show()

def cor_info(ring, CfamilyNames=None):

    cor_names_ = []
    elements_indexes = get_refpts(ring, elements.Corrector)
    for i in elements_indexes:
        Cname = ring[i].FamName
        cor_names_.append(Cname)

    if CfamilyNames == None:
        c_names = cor_names_
    else:
        c_names = CfamilyNames

    Cor_strength = []
    elements_name = []
    s_pos =[]
    cor_indexes = get_refpts(ring, elements.Corrector)
    for i in cor_indexes:
        #print(ring[i].FamName)
        #print(CfamilyNames)
        if ring[i].FamName in c_names:

            cor_strength = ring[i].KickAngle
            Cor_strength.append(cor_strength)
            element_name_ = ring[i].FamName
            elements_name.append(element_name_)


    output = [elements_name[:e].count(v) for e, v in enumerate(elements_name, 0)]
    output1 = [elements_name[:e].count(v) for e, v in enumerate(elements_name, 1)]


    cors = {'Cor_strength': Cor_strength,
            'elements_name': elements_name, 'occ': output, 'occ1':output1}
    cor = pd.DataFrame(cors)



    return cor

def used_cor(j, mylist):
    deduped = list(dict.fromkeys(mylist))
    # Slice off all but the part you care about:
    return deduped[::j]



def quad_info(ring, familyNames=None):

    quad_names_ = []
    elements_indexes = get_refpts(ring, elements.Quadrupole)
    for i in elements_indexes:
        Qname = ring[i].FamName
        quad_names_.append(Qname)

    if familyNames == None:
        quad_names = quad_names_
    else:
        quad_names = familyNames

    Quad_strength = []
    elements_name = []
    s_pos =[]
    quad_indexes = get_refpts(ring, elements.Quadrupole)
    for i in quad_indexes:
        if ring[i].FamName in quad_names:

            quad_strength = ring[i].K
            Quad_strength.append(quad_strength)
            element_name_ = ring[i].FamName
            elements_name.append(element_name_)


    output = [elements_name[:e].count(v) for e, v in enumerate(elements_name, 0)]
    output1 = [elements_name[:e].count(v) for e, v in enumerate(elements_name, 1)]


    quads = {'Quad_strength': Quad_strength,
            'elements_name': elements_name, 'occ': output, 'occ1':output1}
    quad = pd.DataFrame(quads)



    return quad


def getQuadFamilies(Quads):

    n_list = len(Quads.elements_name)
    eocc_a = {}
    vals = {}

    for idx in range(n_list):
        ename = Quads.elements_name[idx]
        par = Quads.Quad_strength[idx]
        eocc =  Quads.occ[idx]
        eocc_a[ename] = int(eocc)
        vals[ename, eocc] = float(par)

    return eocc_a, vals


def ORM_x(dkick, ring, BPMs_random_noise, used_correctors_Names=None):
    cxx = []
    cxy = []

    for i in range(len(used_correctors_Names)):

        cor_index = get_refpts(ring, used_correctors_Names[i])
        cor_index = cor_index[0]

        ring[cor_index].KickAngle[0] = dkick

        closed_orbitx, closed_orbity = closed_orbit_bpm(ring, BPMs_random_noise)
        cxx.append(closed_orbitx)
        cxy.append(closed_orbity)

        ring[cor_index].KickAngle[0] = 0.00

    Cxx = np.squeeze(cxx) / dkick
    Cxy = np.squeeze(cxy) / dkick

    return Cxx, Cxy


def ORM_y(dkick, ring, bpm_random_noise,  CfamilyNames=None):
    cyy = []
    cyx = []

    for i in range(len(CfamilyNames)):
        cor_index  = get_refpts(ring, CfamilyNames[i])
        cor_index = cor_index[0]

        ring[cor_index].KickAngle[1] = dkick

        closed_orbitx, closed_orbity = closed_orbit_bpm(ring, bpm_random_noise)

        cyy.append(closed_orbity)
        cyx.append(closed_orbitx)

        ring[cor_index].KickAngle[1] = 0.00
            #print(j)

    Cyy = np.squeeze(cyy) / dkick
    Cyx = np.squeeze(cyx) / dkick

    return Cyy, Cyx







def getCorFamilies(Cor):

    n_list = len(Cor.elements_name)
    eocc_a = {}
    vals = {}

    for idx in range(n_list):
        ename = Cor.elements_name[idx]
        par = Cor.Cor_strength[idx]
        eocc =  Cor.occ[idx]
        eocc_a[ename] = int(eocc)
        #vals[ename, eocc] = float(par)

    return eocc_a
def used_elements(j, mylist):
    # dedup, preserving order (dict is insertion-ordered as a language guarantee as of 3.7):
    deduped = list(dict.fromkeys(mylist))
    # Slice off all but the part you care about:
    return deduped[::j]













def make_plot_all(twiss, plot_name,file_name, min, max, ymin, ymax):
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
    p4, = par.plot(twiss.s_pos, twiss.dy, label=r'$\eta_y$')
    plt.xlim(min, max)
    plt.ylim(ymin, ymax)
    leg = plt.legend()

    host.yaxis.get_label().set_color(p1.get_color())
    leg.texts[0].set_color(p1.get_color())

    host.yaxis.get_label().set_color(p2.get_color())
    leg.texts[1].set_color(p2.get_color())

    par.yaxis.get_label().set_color(p3.get_color())
    leg.texts[2].set_color(p3.get_color())

    plt.title(plot_name)
    plt.savefig(file_name)

    plt.show()



def  defineMatrices_eta(C0x, C0y, Cxx_err, Cyy_err, dCx, dCy):
    Nk = len(dCx)  # number of free parameters
    Nm = len(dCx)   # number of measurements
    print('NK:', Nk)
    print('Nm:', Nm)

    Ax = np.zeros([Nk, Nk])
    Ay = np.zeros([Nk, Nk])

    A = np.zeros([2 * Nk, Nk])

    ##

    Bx = np.zeros([Nk, 1])
    By = np.zeros([Nk, 1])

    B = np.zeros([2 * Nk, 1])

    ##

    Dx = (Cxx_err - C0x)  ### dk ?
    Dy = (Cyy_err - C0y)


    ##

    for i in range(Nk):  ## i represents each quad
        # print('done A:', 100.* i ,'%')
        for j in range(Nk):
            Ax[i, j] = np.sum(np.dot(dCx[i], dCx[j].T))
            Ay[i, j] = np.sum(np.dot(dCy[i], dCy[j].T))

        A[i, :] = Ax[i, :]
        A[i + Nk, :] = Ay[i, :]


    ##

    for i in range(Nk):
        Bx[i] = np.sum(np.dot(dCx[i], Dx.T))
        By[i] = np.sum(np.dot(dCy[i], Dy.T))


        B[i] = Bx[i]
        B[i + Nk] = By[i]


    return A, B


def svd_plot_(A, B,Nk):

    u, s, v = np.linalg.svd(A, full_matrices=True)

    smat = 0.0 * A
    si = s ** -1


    print("number of singular values {}".format(len(si)))
    smat[:Nk, :Nk] = np.diag(si)

    print('A' + str(A.shape), 'B' + str(B.shape), 'U' + str(u.shape), 'smat' + str(smat.shape), 'v' + str(v.shape))

    plt.plot(np.log(s), 'd--')
    plt.title('singular value')
    plt.show()


    return si , smat, u, v


def compare_orm(Cxy, Cxy_err, Cxy_corr, no):
    # plot the 3 sets
    plt.plot(Cxy[no], label='C')
    plt.plot(Cxy_err[no], label='C_err')
    plt.plot(Cxy_corr[no], label='C_corr')

    # call with no parameters
    plt.legend()

    plt.show()






def getEta(etax, etay, eta_errorx, eta_errory, ring):
    print("getBetaBeat bx and by: ")

    elements_indexes = get_refpts(ring, elements.Monitor)

    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=elements_indexes)
    s_pos = lindata['s_pos']

    bxi =[]
    for i in range(len(etax)):
        bxx = (eta_errorx[i] - etax[i]) / etax[i]
        bxi.append(bxx)
    fig = plt.figure()
    plt.plot(s_pos, bxi)
    plt.xlabel('s[m]')
    plt.ylabel(r'$d\eta_x$')
    #plt.ioff()
    #plt.savefig(file_namex)
    plt.show()


    byi =[]
    for i in range(len(etay)):
        byy = (eta_errory[i] - etay[i]) / etay[i]
        byi.append(byy)
    fig = plt.figure()
    plt.plot(s_pos, byi)
    plt.xlabel('s[m]')
    plt.ylabel(r'$d\eta_y$')

    #plt.ioff()

    #plt.savefig(file_namey)

    plt.show()


    bx = np.std((eta_errorx- etax) / etax)
    by = np.std((eta_errory - etay) / etay)

    bxrms = np.array((eta_errorx- etax) / etax)
    byrms = np.array((eta_errory - etay) / etay)

    bx_rms = np.sqrt(np.mean(bxrms ** 2))
    by_rms = np.sqrt(np.mean(byrms ** 2))

    print("RMS Dispersion, x:" + str(bx_rms * 100) + "%   y: " + str(by_rms * 100) + "%")
    print("STD Dispersion, x:" + str(bx * 100) + "%   y: " + str(by* 100) + "%")

    return bx_rms*100, by_rms*100



def plotORM(orm, file_name):
    plt.figure()
    imshow(orm)
    #plt.savefig(file_name)
    plt.show()


def getBetaBeat(twiss, twiss_error, file_namex, file_namey):
    #print("getBetaBeat bx and by: ")

    bxi =[]
    for i in range(len(twiss.betax)):
        bxx = (twiss_error.betax[i] - twiss.betax[i]) / twiss.betax[i]
        bxi.append(bxx)
    fig = plt.figure()
    plt.plot(twiss.s_pos, bxi)
    plt.xlabel('s[m]')
    plt.ylabel('Beta_Beating_x')
    plt.ioff()
    plt.savefig(file_namex)
    plt.show()


    byi =[]
    for i in range(len(twiss.betay)):
        byy = (twiss_error.betay[i] - twiss.betay[i]) / twiss.betay[i]
        byi.append(byy)
    fig = plt.figure()
    plt.plot(twiss.s_pos, byi)
    plt.xlabel('s[m]')
    plt.ylabel('Beta_Beating_y')

    plt.ioff()
    plt.savefig(file_namey)

    plt.show()


    bx = np.std((twiss_error.betax - twiss.betax) / twiss.betax)
    by = np.std((twiss_error.betay - twiss.betay) / twiss.betay)

    bxrms = np.array((twiss_error.betax - twiss.betax) / twiss.betax)
    byrms = np.array((twiss_error.betay - twiss.betay) / twiss.betay)

    bx_rms = np.sqrt(np.mean(bxrms ** 2))
    by_rms = np.sqrt(np.mean(byrms ** 2))

    print("RMS beta beat, x:" + str(bx_rms * 100) + "%   y: " + str(by_rms * 100) + "%")
    print("STD beta beat, x:" + str(bx * 100) + "%   y: " + str(by* 100) + "%")

    return bx_rms*100, by_rms*100


def plot_closedOrbit(ring, refpts):

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

    plt.plot(s_pos, closed_orbitx/1.e-6)

    # Label for x-axis
    plt.xlabel("s_pos[m]", )

    # Label for y-axis
    plt.ylabel("closed_orbit x [m]")

    # for display

    i = 0
    S_pos2 = []
    plt.xlim(0, 4000)
    plt.title("Closed orbit x")
    plt.show()


    plt.plot(s_pos, closed_orbity)

    # Label for x-axis
    plt.xlabel("s_pos[m]")

    # Label for y-axis
    plt.ylabel("closed_orbit y[m]")

    # for display

    i = 0
    S_pos2 = []
    plt.xlim(0, 50)
    plt.title("Closed orbit y")
    plt.show()




def eta(ring):
    bpm_indexes = get_refpts(ring, elements.Monitor)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    Eta_xx = lindata['dispersion'][:, 0]
    Eta_yy = lindata['dispersion'][:, 2]

    return Eta_xx,Eta_yy


def compare_drm(Cxy, Cxy_err, Cxy_corr):
    # plot the 3 sets
    plt.plot(Cxy, label='$\eta$')
    plt.plot(Cxy_err, label='$\eta_{err}$')
    plt.plot(Cxy_corr, label='$\eta_{corr}$')

    # call with no parameters
    plt.legend()

    plt.show()



def simulateGradientErrors(lattice, gradErr, familyNames=None, debug=False):

    #for i in range(len(familyNames)):
    #    qnames = familyNames[i]
    #    if qnames[0].startswith('qf'):
    #        quad_index = at.get_refpts(lattice, 'qf*')
    #        quads_index.append(np.squeeze(quad_index))
    #    else:
    #        quad_index = at.get_refpts(lattice, 'qd*')
    #        quads_index.append(np.squeeze(quad_index))

    for i in range(len(familyNames)):
        ind = familyNames[i]

        for i in ind:

            a = (1 + gradErr * 0.1)
            lattice[i].K *= a


    if debug: print(f"done")

def simulateAlignmentErrors(lattice, tilt, familyNames=None, debug=False):


    # for i in range(len(familyNames)):
    #    qnames = familyNames[i]
    #    if qnames[0].startswith('qf'):
    #        quad_index = at.get_refpts(lattice, 'qf*')
    #        quads_index.append(np.squeeze(quad_index))
    #    else:
    #        quad_index = at.get_refpts(lattice, 'qd*')
    #        quads_index.append(np.squeeze(quad_index))

    for i in range(len(familyNames)):
        ind = familyNames[i]

        for i in ind:

            at.tilt_elem(lattice[i], tilt, relative=False)

    if debug: print(f"done")



def compare_orm(Cxy, Cxy_err, Cxy_corr, no):
    # plot the 3 sets
    plt.plot(Cxy[no], label='C')
    plt.plot(Cxy_err[no], label='C_err')
    plt.plot(Cxy_corr[no], label='C_corr')

    # call with no parameters
    plt.legend()

    plt.show()

########################################################



def used_quads_f(ring, used_quads_list):

    quads_indexes =[]
    s_pos1 =[]
    for i in used_quads_list:
        element_index = get_refpts(ring, i)
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=element_index)
        s_pos_ = lindata['s_pos']
        s_pos1.append(np.squeeze(s_pos_))
    s_pos= sorted(np.concatenate([x.flatten() for x in s_pos1]))
    return s_pos


def used_elements_plot(lattice, id_Q, min, max):

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
    plt.plot(s_pos, beta_x)
    plt.xlim(min, max)

    # Label for x-axis
    plt.xlabel("s[m]")
    # Label for y-axis
    plt.ylabel("beta_x[m]")

    i = 0
    for i in id_Q:
        for j in i:
            scatter(j, 0)
    plt.title("used quadrupoles indices")
    plt.show()

    #second_plot

    plt.rc('font', size=15)
    plt.plot(s_pos, beta_x)

    plt.xlim(min, max)

    # Label for x-axis
    plt.xlabel("s[m]")
    # Label for y-axis
    plt.ylabel("beta_x[m]")

    i = 0
    for i in id_Q:
        scatter(i, 0)
    plt.title("used quadrupoles indices")
    plt.show()



def coupling_parameters(ring,refpt, qs_family_name, noise):
    elements_indexes = get_refpts(ring, refpt)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=elements_indexes)
    closed_orbitx, closed_orbity = closed_orbit_bpm_r(ring, noise)
    s_pos = lindata['s_pos']
    closed_orbit = lindata['closed_orbit']
    beta_x = lindata['beta'][:, 0]
    beta_y = lindata['beta'][:, 1]

    alls = np.append(np.array(s_pos), ring.circumference)
    ds = diff(alls)

    phase_x = np.dot(1 / beta_x, ds)
    phase_y = np.dot(1 / beta_y, ds)

    #print('phase_x', phase_x)
    #print('phase_y', phase_y)

    #print("Tunex", phase_x / 2 * (pi))
    #print("Tuney", phase_y / 2 * (pi))

    # ring.get_tune()

    i = 0
    k_qs = []
    while (i < len(elements_indexes)):

        if (ring[i].FamName).startswith('qs'):

            skew_quad_strengths = ring[i].K


        else:
            skew_quad_strengths = 0

        k_qs.append(skew_quad_strengths)

        i += 1

    diff_resonances = k_qs * np.sqrt(beta_x * beta_y) * np.exp(
        1j * (phase_x - phase_y - (tune[0] - tune[1])) * 2 * pi * s_pos / ring.circumference)
    sum_resonances = k_qs * np.sqrt(beta_x * beta_y) * np.exp(
        1j * (phase_x + phase_y - (tune[0] + tune[1])) * 2 * pi * s_pos / ring.circumference)

    I1 = np.dot(sum_resonances, ds)
    I2 = np.dot(diff_resonances, ds)

    k1 = I1 / 2 * pi
    k2 = I2 / 2 * pi
    print('Coupling Coefficients (driving term for sum_resonances)  k : ', k1)
    print('Coupling Coefficients (driving term for diff_resonances)  k : ', k2)

def generatingQuadsResponse(ring,dk, Cxx, Cyy,Cxy, Cyx, bpm_random_noise, familyIND=None,CfamilyNames= None, useUniqueDevices=True):

    quad_names_ = []
    qxx = []
    qxy = []
    qyy = []
    qyx = []

    if familyIND == None:
        familyIND = quad_names_
    else:
        familyIND = familyIND


        i=0
        while(i < len(familyIND)):
            print('generating response to family ',i+1)
            t0 = time.time()

            Qxx, Qxy, Qyy, Qyx = QsensitivityMatrices(ring, familyIND[i], dk, bpm_random_noise, CfamilyNames)

            qxx.append(Qxx)
            qxy.append(Qxy)
            qyy.append(Qyy)
            qyx.append(Qyx)
            t1 = time.time()
            print(f"Execution time: {t1 - t0} sec")

            i+=1

    C0x = Cxx
    C0y = Cyy
    C0xy = Cxy
    C0yx = Cyx

    dCx = []
    dCy = []
    dCxy = []
    dCyx = []
    # quad_names = quads
    # for qname in quad_names:
    #     # nquad = quad_dict[qname]
    #     print('loading response to:', qname)
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


def simulateFixedGradientErrors(lattice, gradErr, familyNames=None, debug=False):

    #for i in range(len(familyNames)):
    #    qnames = familyNames[i]
    #    if qnames[0].startswith('qf'):
    #        quad_index = at.get_refpts(lattice, 'qf*')
    #        quads_index.append(np.squeeze(quad_index))
    #    else:
    #        quad_index = at.get_refpts(lattice, 'qd*')
    #        quads_index.append(np.squeeze(quad_index))

    for i in range(len(familyNames)):
        ind = familyNames[i]

        a = (1 + gradErr * np.random.randn())

        for j in ind:
            lattice[j].K *= a


    if debug: print(f"done")


def simulateRandomGradientErrors(lattice, gradErr, familyNames=None, debug=False):


    # for i in range(len(familyNames)):
    #    qnames = familyNames[i]
    #    if qnames[0].startswith('qf'):
    #        quad_index = at.get_refpts(lattice, 'qf*')
    #        quads_index.append(np.squeeze(quad_index))
    #    else:
    #        quad_index = at.get_refpts(lattice, 'qd*')
    #        quads_index.append(np.squeeze(quad_index))

    for i in range(len(familyNames)):
        ind = familyNames[i]

        for i in ind:

            lattice[i].K *= (1 + gradErr * np.random.randn())

    if debug: print(f"done")


def getBetaBeat(twiss, twiss_error, file_namex, file_namey):
    #print("getBetaBeat bx and by: ")

    bxi =[]
    for i in range(len(twiss.betax)):
        bxx = (twiss_error.betax[i] - twiss.betax[i]) / twiss.betax[i]
        bxi.append(bxx)
    fig = plt.figure()
    plt.plot(twiss.s_pos, bxi)
    plt.xlabel('s[m]')
    plt.ylabel('Beta_Beating_x')
    plt.ioff()
    plt.savefig(file_namex)
    plt.show()


    byi =[]
    for i in range(len(twiss.betay)):
        byy = (twiss_error.betay[i] - twiss.betay[i]) / twiss.betay[i]
        byi.append(byy)
    fig = plt.figure()
    plt.plot(twiss.s_pos, byi)
    plt.xlabel('s[m]')
    plt.ylabel('Beta_Beating_y')

    plt.ioff()
    plt.savefig(file_namey)

    plt.show()


    bx = np.std((twiss_error.betax - twiss.betax) / twiss.betax)
    by = np.std((twiss_error.betay - twiss.betay) / twiss.betay)

    bxrms = np.array((twiss_error.betax - twiss.betax) / twiss.betax)
    byrms = np.array((twiss_error.betay - twiss.betay) / twiss.betay)

    bx_rms = np.sqrt(np.mean(bxrms ** 2))
    by_rms = np.sqrt(np.mean(byrms ** 2))

    print("RMS beta beat, x:" + str(bx_rms * 100) + "%   y: " + str(by_rms * 100) + "%")
    print("STD beta beat, x:" + str(bx * 100) + "%   y: " + str(by* 100) + "%")

    return bx_rms*100, by_rms*100



def defineMatrices(C0x, C0y, C0xy, C0yx, Cxx_err, Cyy_err, Cxy_err, Cyx_err, dCx, dCy, dCxy,dCyx):
    Nk = len(dCx)  # number of free parameters
    Nm = len(C0x)  # number of measurements
    #print('NK:', Nk)
    #print('Nm:', Nm)

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

    Dx = (Cxx_err[0:Nm, :] - C0x[0:Nm, :])  ### dk ?
    Dy = (Cyy_err[0:Nm, :] - C0y[0:Nm, :])
    Dxy = (Cxy_err[0:Nm, :] - C0xy[0:Nm, :])
    Dyx = (Cyx_err[0:Nm, :] - C0yx[0:Nm, :])

    ##

    for i in range(Nk):  ## i represents each quad
        # print('done A:', 100.* i ,'%')
        for j in range(Nk):
            Ax[i, j] = np.sum(np.dot(dCx[i], dCx[j].T))
            Ay[i, j] = np.sum(np.dot(dCy[i], dCy[j].T))
            Axy[i, j] = np.sum(np.dot(dCxy[i], dCxy[j].T))
            Ayx[i, j] = np.sum(np.dot(dCyx[i], dCyx[j].T))
        A[i, :] = Ax[i, :]
        A[i + Nk, :] = Ay[i, :]
        A[i + 2 * Nk, :] = Axy[i, :]
        A[i + 3 * Nk, :] = Ayx[i, :]

    ##

    for i in range(Nk):
        Bx[i] = np.sum(np.dot(dCx[i], Dx.T))
        By[i] = np.sum(np.dot(dCy[i], Dy.T))
        Bxy[i] = np.sum(np.dot(dCxy[i], Dxy.T))
        Byx[i] = np.sum(np.dot(dCyx[i], Dyx.T))
        B[i] = Bx[i]
        B[i + Nk] = By[i]
        B[i + 2 * Nk] = Bxy[i]
        B[i + 3 * Nk] = Byx[i]

    return A, B



def getInverse(A, B,Nk, sCut):

    u, s, v = np.linalg.svd(A, full_matrices=True)
    #plt.plot(np.log(s), 'd--')
    #plt.title('singular value')
    #plt.show()

    smat = 0.0 * A
    si = s ** -1
    n_sv = sCut
    si[n_sv:] *= 0.0

    #print("number of singular values {}".format(len(si)))
    smat[:Nk, :Nk] = np.diag(si)

    #print('A' + str(A.shape), 'B' + str(B.shape), 'U' + str(u.shape), 'smat' + str(smat.shape), 'v' + str(v.shape))



    Ai = np.dot(v.transpose(), np.dot(smat.transpose(), u.transpose()))

    ###
    r = (np.dot(Ai, B)).reshape(-1)
    #plot(r, 'd')
    #plt.show()

    # error
    e = np.dot(A, r).reshape(-1) - B.reshape(-1)
    #plt.plot(e)
    #plt.show()
    #plt.plot(B)
    #plt.show()

    return Ai, r, e


def setCorrection(ring, r ,familyInd, familyNames = None):

    iq = 0
    frac = 1.0
    cor_dict = {}

    quads_index_ = []
    qnames_all = []
    for i in range(len(familyNames)):
        qnames = familyNames[i]

        for qname in qnames:
            cor_dict[qname] = -r[i] * frac

    quads_index = [item for sublist in familyInd for item in sublist]
    qnames_all = [item for sublist in familyNames for item in sublist]

    DK = []
    for idx in qnames_all:

        dk = cor_dict[idx]
        DK.append(dk)

    #print("DK------------------", DK)

    for j in range(len(quads_index)):
        ring[quads_index[j]].K += DK[j]

    return qnames_all


def QsensitivityMatrices(ring, qfamilyIndices, dk, BPMs_random_noise, used_correctors_Names):

    strength_before = []
    for i in qfamilyIndices:

        strength_before.append(ring[i].K)
        a = ring[i].K
        ring[i].K = a + dk

    qxx, qxy = ORM_x(dk, ring, BPMs_random_noise, used_correctors_Names)
    qyy, qyx = ORM_y(dk, ring, BPMs_random_noise, used_correctors_Names)

    for i in range(len(qfamilyIndices)):

        ring[qfamilyIndices[i]].K = strength_before[i]


    return  qxx, qxy, qyy, qyx



