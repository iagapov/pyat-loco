import numpy as np
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
    [elemdata0, beamdata, elemdata] = at.get_optics(ring, elements_indexes)
    twiss = elemdata

    return twiss

def make_plot_all2(twiss, plot_name):

    from mpl_toolkits.axes_grid1 import host_subplot
    import matplotlib.pyplot as plt

    plt.plot(twiss.s_pos, twiss.beta)
    plt.xlabel('s [m]')
    plt.ylabel(r'$\beta$ [m]')

    host = host_subplot(111)

    par = host.twinx()

    host.set_xlabel("s_pos")
    host.set_ylabel(r'$\beta_x$')
    host.set_ylabel(r'$\beta_y$')
    par.set_ylabel("dx")

    p1, = host.plot(twiss.s_pos, twiss.beta[:, 0], label=r'$\beta_x$')
    p2, = host.plot(twiss.s_pos, twiss.beta[:, 1], label=r'$\beta_y$')
    p3, = par.plot(twiss.s_pos, twiss.dispersion[:, 0], label=r'$\eta_x$')
    p4, = par.plot(twiss.s_pos, twiss.dispersion[:, 2], label=r'$\eta_y$')

    leg = plt.legend()

    host.yaxis.get_label().set_color(p1.get_color())
    leg.texts[0].set_color(p1.get_color())

    host.yaxis.get_label().set_color(p2.get_color())
    leg.texts[1].set_color(p2.get_color())

    par.yaxis.get_label().set_color(p3.get_color())
    leg.texts[2].set_color(p3.get_color())

    plt.title(plot_name)
    plt.show()


def used_cor(j, mylist):
    #deduped = list(dict.fromkeys(mylist))
    # Slice off all but the part you care about:
    return mylist[::j]


def used_quadrpoles(ring, steps):

    quad_indexes = get_refpts(ring, elements.Quadrupole)
    # all_quads_names
    all_quads_names = [ring[ind].FamName for ind in quad_indexes]


    qf_indexes = get_refpts(ring, 'QF*')
    qd_indexes = get_refpts(ring, 'QD*')

    print("# of QF:", len(qf_indexes), '# of QD:', len(qd_indexes), )



    # choosing the used Quads families
    qd_names = []
    qf_names = []
    for i in all_quads_names:
        if i.startswith('QD'):
            qd_names.append(i)
        if i.startswith('QF'):
            qf_names.append(i)

    used_qf = [qf_names[x:x + steps] for x in range(0, len(qf_names), steps)]  # 720
    used_qd = [qd_names[x:x + steps] for x in range(0, len(qd_names), steps)]

    used_qf_ind = [qf_indexes[x:x + steps] for x in range(0, len(qf_indexes), steps)]
    used_qd_ind = [qd_indexes[x:x + steps] for x in range(0, len(qd_indexes), steps)]


    used_quadrpoles_families = used_qd + used_qf
    used_quadrpoles_families_ind = used_qd_ind + used_qf_ind
    print("used_quadrpoles_families_ind: ",len(used_quadrpoles_families_ind))

    return used_quadrpoles_families, used_quadrpoles_families_ind


def used_skew(ring, steps):
    quad_indexes = get_refpts(ring, elements.Quadrupole)
    # all_quads_names
    all_quads_names = [ring[ind].FamName for ind in quad_indexes]

    qs_indexes = get_refpts(ring, 'QS*')
    #qd_indexes = get_refpts(ring, 'QD*')

    print("# of QS:", len(qs_indexes) )

    # choosing the used Quads families

    qs_names = np.array([ring[ind].FamName for ind in qs_indexes])

    used_qs = [qs_names[x:x + steps] for x in range(0, len(qs_names), steps)]  # 720

    used_qs_ind = [qs_indexes[x:x + steps] for x in range(0, len(qs_indexes), steps)]

    used_qs_families = used_qs
    used_qs_families_ind = used_qs_ind
    print("used_qs_families_ind: ", len(used_qs_families_ind))

    return used_qs_families, used_qs_families_ind


def used_quadrpoles_fcc(ring, step):

    quad_indexes = get_refpts(ring, elements.Quadrupole)
    # all_quads_names
    all_quads_names = [ring[ind].FamName for ind in quad_indexes]


    qfg_indexes = get_refpts(ring, 'QFG*')
    qdg_indexes = get_refpts(ring, 'QDG*')
    qf2_indexes = get_refpts(ring, 'QF2*')
    qd1_indexes = get_refpts(ring, 'QD1*')
    qf4_indexes = get_refpts(ring, 'QF4*')
    qd3_indexes = get_refpts(ring, 'QD3*')
    print("# of QDG1:", len(qdg_indexes), '# of QFG2:', len(qfg_indexes), '# of QD1:', len(qd1_indexes), '# of QF2:',
          len(qf2_indexes), '# of QD3:', len(qd3_indexes), '# of QF4:', len(qf4_indexes))

    # choosing the used Quads families
    qdg1_names = []
    qfg2_names = []
    qd1_names = []
    qf2_names = []
    qd3_names = []
    qf4_names = []
    used_quadrpoles_families = ['qd1', 'qf2', 'qd3', 'qf4']  # did not choose from 'qdg', 'qfg',
    for i in all_quads_names:
        if i.startswith('qdg'):
            qdg1_names.append(i)
        if i.startswith('qfg'):
            qfg2_names.append(i)
        if i.startswith('qd1'):
            qd1_names.append(i)
        if i.startswith('qf2'):
            qf2_names.append(i)
        if i.startswith('qd3'):
            qd3_names.append(i)
        if i.startswith('qf4'):
            qf4_names.append(i)
    used_qf4 = [qf4_names[x:x + step] for x in range(0, len(qf4_names), step)]  # 720
    used_qd1 = [qd1_names[x:x + step] for x in range(0, len(qd1_names), step)]
    used_qf2 = [qf2_names[x:x + step] for x in range(0, len(qf2_names), step)]
    used_qd3 = [qd3_names[x:x + step] for x in range(0, len(qd3_names), step)]

    used_qf4_ind = [qf4_indexes[x:x + step] for x in range(0, len(qf4_indexes), step)]
    used_qd1_ind = [qd1_indexes[x:x + step] for x in range(0, len(qd1_indexes), step)]
    used_qf2_ind = [qf2_indexes[x:x + step] for x in range(0, len(qf2_indexes), step)]
    used_qd3_ind = [qd3_indexes[x:x + step] for x in range(0, len(qd3_indexes), step)]

    #print(len(used_qf4))
    used_quadrpoles_families = used_qd1 + used_qf2 + used_qd3 + used_qf4
    used_quadrpoles_families_ind = used_qd1_ind + used_qf2_ind + used_qd3_ind + used_qf4_ind
    print(len(used_quadrpoles_families_ind))

    return used_quadrpoles_families, used_quadrpoles_families_ind



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


def closed_orbit_bpm(ring, bpm_random_noise):

    bpm_indexes = get_refpts(ring, elements.Monitor)

    [elemdata0, beamdata, elemdata] = at.get_optics(ring, bpm_indexes)
    closed_orbitx1 = elemdata.closed_orbit[:, 0]
    closed_orbity1 = elemdata.closed_orbit[:, 2]

    Eta_xx1 = elemdata.dispersion[:, 0]
    Eta_yy1 = elemdata.dispersion[:, 2]

    bpm_noise = bpm_random_noise * np.random.randn(len(bpm_indexes))

    closed_orbitx = closed_orbitx1 + bpm_noise
    closed_orbity = closed_orbity1 + bpm_noise


    return closed_orbitx, closed_orbity


def generatingQuadsResponse(ring, dk, Cxx, Cyy, Cxy, Cyx, bpm_random_noise, familyIND=None, CfamilyNames=None,
                            useUniqueDevices=True):

    quad_names_ = []
    C0x = Cxx
    C0y = Cyy
    C0xy = Cxy
    C0yx = Cyx

    dCx = []
    dCy = []
    dCxy = []
    dCyx = []

    if familyIND == None:
        familyIND = quad_names_
    else:
        familyIND = familyIND


    i = 0
    while (i < len(familyIND)):
        print('generating response to family ', i + 1)
        t0 = time.time()

        Qxx, Qxy, Qyy, Qyx = QsensitivityMatrices(ring, familyIND[i], dk, bpm_random_noise, CfamilyNames)

        C1x = Qxx
        C1y = Qyy
        C1xy = Qxy
        C1yx = Qyx

        dCx.append((C1x - C0x) / dk)
        dCy.append((C1y - C0y) / dk)

        dCxy.append((C1xy - C0xy) / dk)
        dCyx.append((C1yx - C0yx) / dk)

        t1 = time.time()
        print(f"Execution time: {t1 - t0} sec")

        i += 1

    return C0x, C0y, C0xy, C0yx, dCx, dCy, dCxy, dCyx


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



def simulateFixedGradientErrors(lattice, gradErr, familyIND=None, debug=False):


    for i in range(len(familyIND)):
        ind = familyIND[i]

        a = (1 + gradErr * np.random.randn())

        for j in ind:
            lattice[j].K *= a


    if debug: print(f"done")


def simulateRandomGradientErrors(lattice, gradErr, familyIND=None, debug=False):


    for i in range(len(familyIND)):
        ind = familyIND[i]

        for i in ind:

            lattice[i].K *= (1 + gradErr * np.random.randn())

    if debug: print(f"done")

def simulateAlignmentErrors(lattice, tilt, familyIND=None, debug=False):

    for i in range(len(familyIND)):
        ind = familyIND[i]

        for i in ind:

            at.tilt_elem(lattice[i], tilt, relative=False)

    if debug: print(f"done")

def getBetaBeat(twiss, twiss_error):
    #print("getBetaBeat bx and by: ")

    bxrms = np.array((twiss_error.beta[:, 0] - twiss.beta[:, 0]) / twiss.beta[:, 0])
    byrms = np.array((twiss_error.beta[:, 1] - twiss.beta[:, 1]) / twiss.beta[:, 1])

    bx_std = np.std(bxrms)
    by_std = np.std(byrms)

    bx_rms = np.sqrt(np.mean(bxrms ** 2))
    by_rms = np.sqrt(np.mean(byrms ** 2))

    #fig = plt.figure()
    #plt.plot(twiss.s_pos, bxrms)
    #plt.xlabel('s[m]')
    #plt.ylabel('Beta_Beating_x')

    #plt.show()

    #fig = plt.figure()
    #plt.plot(twiss.s_pos, byrms)
    #plt.xlabel('s[m]')
    #plt.ylabel('Beta_Beating_y')

    #plt.show()

    print("RMS beta beat, x:" + str(bx_rms * 100) + "%   y: " + str(by_rms * 100) + "%")
    print("STD beta beat, x:" + str(bx_std * 100) + "%   y: " + str(by_std* 100) + "%")

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
    #e = np.dot(A, r).reshape(-1) - B.reshape(-1)
    #plt.plot(e)
    #plt.show()
    #plt.plot(B)
    #plt.show()

    return  r


def setCorrection(ring, r ,familyInd):

    iq = 0
    frac = 1.0
    cor_dict = {}
    quads_index = [item for sublist in familyInd for item in sublist]


    for i in range(len(familyInd)):
        qInds = familyInd[i]

        for qInd in qInds:
            ring[qInd].K += -r[i]


def eta(ring):
    elements_indexes = get_refpts(ring, '*')
    [elemdata0, beamdata, elemdata] = at.get_optics(ring, elements_indexes)
    Eta_xx = elemdata.dispersion[:, 0]
    Eta_yy = elemdata.dispersion[:, 2]

    return Eta_xx,Eta_yy

def generatingQuadsResponseEta(ring, dk, Cxx, Cyy, familyIND=None,
                            useUniqueDevices=True):

    quad_names_ = []
    C0x = Cxx
    C0y = Cyy

    dCx = []
    dCy = []
    dCxy = []
    dCyx = []

    if familyIND == None:
        familyIND = quad_names_
    else:
        familyIND = familyIND


    i = 0
    while (i < len(familyIND)):
        print('generating response to family ', i + 1)
        t0 = time.time()

        Qxx, Qyy = QsensitivityMatricesEta(ring, familyIND[i], dk)

        C1x = Qxx
        C1y = Qyy


        dCx.append((C1x - C0x) / dk)
        dCy.append((C1y - C0y) / dk)



        t1 = time.time()
        print(f"Execution time: {t1 - t0} sec")

        i += 1

    return C0x, C0y, dCx, dCy

def QsensitivityMatricesEta(ring, qfamilyIndices, dk):

    strength_before = []
    for i in qfamilyIndices:

        strength_before.append(ring[i].K)
        a = ring[i].K
        ring[i].K = a + dk

    qxx, qyy = eta(ring)
    qxx, qyy = eta(ring)

    for i in range(len(qfamilyIndices)):

        ring[qfamilyIndices[i]].K = strength_before[i]


    return  qxx,  qyy


def coupling_parameters(ring, IND):

    elements_indexes = get_refpts(ring, IND)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=elements_indexes)
    s_pos = lindata['s_pos']
    beta_x = lindata['beta'][:, 0]
    beta_y = lindata['beta'][:, 1]

    alls = np.append(np.array(s_pos), ring.circumference)
    ds = diff(alls)

    phase_x = np.dot(1 / beta_x, ds)
    phase_y = np.dot(1 / beta_y, ds)


    i = 0
    k_qs = []
    while (i < len(elements_indexes)):

        if (ring[i].FamName).startswith('QS'):

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
    k1 =abs(k1)
    k2 =abs(k2)


    return k1, k2


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


def getEta(etax, etay, eta_errorx, eta_errory, ring):
    print("get Eeta dx and dy: ")

    bxx = (eta_errorx- etax) / etax
    byy = (eta_errory - etay)


    bx_std = np.std(bxx)
    by_std = np.std(byy)

    bx_rms = np.sqrt(np.mean(bxx ** 2))
    by_rms = np.sqrt(np.mean(byy ** 2))

    print("RMS Dispersion, x:" + str(bx_rms * 100) + "%   y: " + str(by_rms * 100) + "%")
    print("STD Dispersion, x:" + str(bx_std * 100) + "%   y: " + str(by_std* 100) + "%")

    return bx_rms*100, by_rms*100

def compare_drm(Cxy, Cxy_err, Cxy_corr):
    # plot the 3 sets
    plt.plot(Cxy, label='$\eta$')
    plt.plot(Cxy_err, label='$\eta_{err}$')
    plt.plot(Cxy_corr, label='$\eta_{corr}$')

    # call with no parameters
    plt.legend()

    plt.show()

def compare_orm(Cxy, Cxy_err, Cxy_corr, no):
    # plot the 3 sets
    plt.plot(Cxy[no], label='C')
    plt.plot(Cxy_err[no], label='C_err')
    plt.plot(Cxy_corr[no], label='C_corr')

    # call with no parameters
    plt.legend()

    plt.show()

def ORM_x_eta(dkick, ring, BPMs_random_noise, used_correctors_Names=None):
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

    bpm_indexes = get_refpts(ring, elements.Monitor)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    Eta_xx = lindata['dispersion'][:, 0]
    Eta_yy = lindata['dispersion'][:, 2]
    cxx.append(Eta_xx)
    cxx.append(Eta_yy)

    Cxx = np.squeeze(cxx) / dkick
    Cxy = np.squeeze(cxy) / dkick

    return Cxx, Cxy

def ORM_y_eta(dkick, ring, BPMs_random_noise, used_correctors_Names=None):
    cyy = []
    cyx = []

    for i in range(len(used_correctors_Names)):
        cor_index = get_refpts(ring, used_correctors_Names[i])
        cor_index = cor_index[0]

        ring[cor_index].KickAngle[1] = dkick

        closed_orbitx, closed_orbity = closed_orbit_bpm(ring, BPMs_random_noise)
        cyy.append(closed_orbitx)
        cyx.append(closed_orbity)

        ring[cor_index].KickAngle[1] = 0.00

    bpm_indexes = get_refpts(ring, elements.Monitor)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    Eta_xx = lindata['dispersion'][:, 0]
    Eta_yy = lindata['dispersion'][:, 2]
    cyy.append(Eta_xx)
    cyy.append(Eta_yy)

    Cyy = np.squeeze(cyy) / dkick
    Cyx = np.squeeze(cyx) / dkick

    return Cyy, Cyx

def generatingQuadsResponse_eta(ring, dk, Cxx, Cyy, Cxy, Cyx, bpm_random_noise, familyIND=None, CfamilyNames=None,
                            useUniqueDevices=True):

    quad_names_ = []
    C0x = Cxx
    C0y = Cyy
    C0xy = Cxy
    C0yx = Cyx

    dCx = []
    dCy = []
    dCxy = []
    dCyx = []

    if familyIND == None:
        familyIND = quad_names_
    else:
        familyIND = familyIND


    i = 0
    while (i < len(familyIND)):
        print('generating response to family ', i + 1)
        t0 = time.time()

        Qxx, Qxy, Qyy, Qyx = QsensitivityMatrices_eta(ring, familyIND[i], dk, bpm_random_noise, CfamilyNames)

        C1x = Qxx
        C1y = Qyy
        C1xy = Qxy
        C1yx = Qyx

        dCx.append((C1x - C0x) / dk)
        dCy.append((C1y - C0y) / dk)

        dCxy.append((C1xy - C0xy) / dk)
        dCyx.append((C1yx - C0yx) / dk)

        t1 = time.time()
        print(f"Execution time: {t1 - t0} sec")

        i += 1

    return C0x, C0y, C0xy, C0yx, dCx, dCy, dCxy, dCyx


def QsensitivityMatrices_eta(ring, qfamilyIndices, dk, BPMs_random_noise, used_correctors_Names):

    strength_before = []
    for i in qfamilyIndices:

        strength_before.append(ring[i].K)
        a = ring[i].K
        ring[i].K = a + dk

    qxx, qxy = ORM_x_eta(dk, ring, BPMs_random_noise, used_correctors_Names)
    qyy, qyx = ORM_y_eta(dk, ring, BPMs_random_noise, used_correctors_Names)

    for i in range(len(qfamilyIndices)):

        ring[qfamilyIndices[i]].K = strength_before[i]


    return  qxx, qxy, qyy, qyx
