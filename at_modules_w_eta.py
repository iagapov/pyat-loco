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

    plt.plot(s_pos, closed_orbitx)

    # Label for x-axis
    plt.xlabel("s_pos")

    # Label for y-axis
    plt.ylabel("closed_orbit x")

    # for display

    i = 0
    S_pos2 = []
    plt.title("Closed orbit x")
    plt.show()


    plt.plot(s_pos, closed_orbity)

    # Label for x-axis
    plt.xlabel("s_pos")

    # Label for y-axis
    plt.ylabel("closed_orbit y")

    # for display

    i = 0
    S_pos2 = []
    plt.title("Closed orbit y")
    plt.show()


def correctionType(alpha1,alpha2, alpha3):

    if alpha1 == 1:
       type = "optics correction"

    if alpha2 == 1:
       type = "dispersion correction"


    if alpha3 == 1:
       type = "optics and dispersion correction"
    print("This code performs: ", type)

    #return type


def func(j, mylist):
    # dedup, preserving order (dict is insertion-ordered as a language guarantee as of 3.7):
    deduped = list(dict.fromkeys(mylist))
    # Slice off all but the part you care about:
    return deduped[::j]


def defineMatrices_w_eta(W, alpha1, alpha2,alpha3, C0x, C0y, C0xy, C0yx, Cxx_err, Cyy_err, Cxy_err, Cyx_err, dCx, dCy, dCxy,dCyx):
    Nk = len(dCx)  # number of free parameters
    Nm = len(dCx) # number of measurements
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

    Dx = (Cxx_err[:, :] - C0x[:, :] )#- error_variance)  ### dk ?
    Dy = (Cyy_err[:, :] - C0y[:, :] )
    Dxy = (Cxy_err[:, :] - C0xy[:, :])
    Dyx = (Cyx_err[:, :] - C0yx[:, :] )


    ##

    for i in range(Nk):  ## i represents each quad
        # print('done A:', 100.* i ,'%')
        for j in range(Nk):
            Ax[i, j] = np.sum(np.dot(np.dot(dCx[i][0: -2, :],W*alpha1), dCx[j][0: -2, :].T)) + np.sum(np.dot(np.dot(dCx[i][ -2 ::, :],W*alpha2), dCx[j][ -2 ::, :].T)) + np.sum(np.dot(np.dot(dCx[i],W*alpha3), dCx[j].T))
            Ay[i, j] = np.sum(np.dot(np.dot(dCy[i][0: -2, :],W*alpha1), dCy[j][0: -2, :].T)) +  np.sum(np.dot(np.dot(dCy[i][ -2 ::, :],W*alpha2), dCy[j][ -2 ::, :].T))+ np.sum(np.dot(np.dot(dCy[i],W*alpha3), dCy[j].T))
            Axy[i, j] = np.sum(np.dot(np.dot(dCxy[i][0: -2, :],W*alpha1), dCxy[j][0: -2, :].T)) + np.sum(np.dot(np.dot(dCxy[i][ -2 ::, :],W*alpha2), dCxy[j][ -2 ::, :].T))+ np.sum(np.dot(np.dot(dCxy[i],W*alpha3), dCxy[j].T))
            Ayx[i, j] = np.sum(np.dot(np.dot(dCyx[i][0: -2, :],W*alpha1), dCyx[j][0: -2, :].T)) + np.sum(np.dot(np.dot(dCyx[i][ -2 ::, :],W*alpha2), dCyx[j][ -2 ::, :].T))+ np.sum(np.dot(np.dot(dCyx[i],W*alpha3), dCyx[j].T))
        A[i, :] = Ax[i, :]
        A[i + Nk, :] = Ay[i, :]
        A[i + 2 * Nk, :] = Axy[i, :]
        A[i + 3 * Nk, :] = Ayx[i, :]


    ##

    for i in range(Nk):

        Bx[i] = np.sum(np.dot(np.dot(dCx[i][0: -2, :],W*alpha1), Dx[0: -2, :].T))+  np.sum(np.dot(np.dot(dCx[i][ -2 ::, :],W*alpha2), Dx[ -2 ::, :].T)) + np.sum(np.dot(np.dot(dCx[i],W*alpha3), Dx.T))
        By[i] = np.sum(np.dot(np.dot(dCy[i][0: -2, :],W*alpha1), Dy[0: -2, :].T)) + np.sum(np.dot(np.dot(dCy[i][ -2 ::, :],W*alpha2), Dy[ -2 ::, :].T))+np.sum(np.dot(np.dot(dCy[i],W*alpha3), Dy.T))
        Bxy[i] = np.sum(np.dot(np.dot(dCxy[i][0: -2, :],W*alpha1), Dxy[0: -2, :].T))+ np.sum(np.dot(np.dot(dCxy[i][ -2 ::, :],W*alpha2), Dxy[ -2 ::, :].T))+np.sum(np.dot(np.dot(dCxy[i],W*alpha3), Dxy.T))
        Byx[i] = np.sum(np.dot(np.dot(dCyx[i][0: -2, :],W*alpha1), Dyx[0: -2, :].T))+ np.sum(np.dot(np.dot(dCyx[i][ -2 ::, :],W*alpha2), Dyx[ -2 ::, :].T))+np.sum(np.dot(np.dot(dCyx[i],W*alpha3), Dyx.T))


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

    plt.plot(np.log(s), 'd--')
    plt.title('singular value')
    plt.show()

    plt.plot(si, 'd--')
    plt.title('singular value inverse')
    plt.show()

    Ai = np.dot(v.transpose(), np.dot(smat.transpose(), u.transpose()))

    ###
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



def compare_orm(Cxy, Cxy_err, Cxy_corr, no):
    # plot the 3 sets
    plt.plot(Cxy[no], label='C')
    plt.plot(Cxy_err[no], label='C_err')
    plt.plot(Cxy_corr[no], label='C_corr')

    # call with no parameters
    plt.legend()

    plt.show()


def compare_drm(Cxy, Cxy_err, Cxy_corr):
    # plot the 3 sets
    plt.plot(Cxy, label='$\eta$')
    plt.plot(Cxy_err, label='$\eta_{err}$')
    plt.plot(Cxy_corr, label='$\eta_{corr}$')

    # call with no parameters
    plt.legend()

    plt.show()





def generatingQuadsResponse1(ring, Cxx, Cyy,Cxy, Cyx , used_correctors):
    # %%time

    quads_info = quad_info(ring)
    quad_dict, quad_vals = getQuadFamilies(quads_info)
    quads = [k for k in quad_dict.keys()]
    quad_names = quads
    dk = 0.0001
    qxx = []
    qxy = []
    qyy = []
    qyx = []
    quad_names = quads
    for qname in quad_names:
        print('generating response to {}, n={}'.format(qname, quad_dict[qname]))
        t0 = time.time()
        nq = quad_dict[qname] + 1
        for i in range(0, nq):
            Qxx, Qxy, Qyy, Qyx = computeOpticsD1(ring, qname, i, dk, quad_vals, used_correctors)
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
    for qname in quad_names:
        # nquad = quad_dict[qname]
        print('loading response to:', qname)
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

def setCorrection1(ring, quads_info_error,quad_names, r , quads_info,n_list):

    quad_dict, quad_vals = getQuadFamilies(quads_info_error)
    n_list = len(quads_info_error.s_pos)
    # print(n_list)

    quad_names = quad_names
    iq = 0
    frac = 1.0
    cor_dict = {}
    for qname in quad_names:
        nquad = quad_dict[qname]
        # print(qname, quad_dict[qname])
        for i in range(0, nquad):
            cor_dict[qname, i + 1] = -r[iq] * frac
            iq += 1
    print("define correction : Done")

    DK = []
    for idx in range(n_list):
        qname_ = quads_info.elements_name[idx]  # ElementName

        occ = quads_info_error.occ[idx]
        dk = cor_dict[qname_, occ]
        DK.append(dk)

    quads_indexes = get_refpts(ring, elements.Quadrupole)
    i = 0
    while (i < len(quads_indexes)):
        ring[quads_indexes[i]].K += DK[i]
        i += 1

    print("set correction : Done")


def plotORM(orm):
    plt.figure()
    imshow(orm)
    plt.show()


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

    bxx = np.array((twiss_error.betax - twiss.betax) / twiss.betax)
    byy = np.array((twiss_error.betay - twiss.betay) / twiss.betay)

    bx = np.sqrt(np.mean(bxx ** 2))
    by = np.sqrt(np.mean(byy ** 2))
    #bx = np.std((twiss_error.betax - twiss.betax) / twiss.betax)
    #by = np.std((twiss_error.betay - twiss.betay) / twiss.betay)
    print("Simulated beta beat, x:" + str(bx * 100) + "%   y: " + str(by* 100) + "%")
def used_elements_plot(lattice, elements_indexes, used_quad):
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



    plt.plot(s_pos, closed_orbitx)
    #plt.plot(s_pos, closed_orbitx)

    # Label for x-axis
    plt.xlabel("elements_indexes")

    # Label for y-axis
    plt.ylabel("closed_orbit_x")

      # for display


    i = 0
    S_pos2 = []
    while (i < used_quad.shape[1]):
        S_pos1 = used_quad.iloc[:, i]
        S_pos_ = df = pd.concat([S_pos1])
        S_pos2.append(S_pos_)
        i += 1
        for i in S_pos_:
            scatter(i, 0)
    plt.title("used quadrupoles indices")
    plt.show()



    plt.plot(s_pos, beta_x)
    #plt.plot(s_pos, beta_x)
    # Label for x-axis
    plt.xlabel("elements_indexes")

    # Label for y-axis
    plt.ylabel("beta_x")

     # for display

    S_pos2 = []
    i = 0
    S_pos2 = []
    while (i < used_quad.shape[1]):
        S_pos1 = used_quad.iloc[:, i]
        S_pos_ = df = pd.concat([S_pos1])
        S_pos2.append(S_pos_)
        i += 1
        for i in S_pos_:
            scatter(i, 0)
    plt.title("used quadrupoles indices")
    plt.show()


def used_elements_plot1(lattice, s_poss, used_quad):
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



    plt.plot(s_pos, closed_orbitx)
    #plt.plot(s_pos, closed_orbitx)

    # Label for x-axis
    plt.xlabel("elements_indexes")

    # Label for y-axis
    plt.ylabel("closed_orbit_x")

      # for display


    i = 0

    for i in s_poss:
        scatter(i, 0)
    plt.title("used quadrupoles indices")
    plt.show()



    plt.plot(s_pos, beta_x)
    #plt.plot(s_pos, beta_x)
    # Label for x-axis
    plt.xlabel("elements_indexes")

    # Label for y-axis
    plt.ylabel("beta_x")

     # for display

    S_pos2 = []
    i = 0
    S_pos2 = []


    for i in s_poss:
        scatter(i, 0)
    plt.title("used quadrupoles indices")
    plt.show()



def getDispersion(twiss, twiss_error, twiss_corrected):
    plt.plot(twiss.dx, label='$\eta_x$')
    plt.plot(twiss_error.dx, label='$\eta_x_err$')
    plt.plot(twiss_corrected.dx, label='$\eta_x_corr$')

    plt.legend()

    plt.show()

    plt.plot(twiss.dy, label='$\eta_y$')
    plt.plot(twiss_error.dy, label='$\eta_y_err$')
    plt.plot(twiss_corrected.dy, label='$\eta_y_corr$')

    plt.legend()

    plt.show()


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
def used_quads_f1(ring, used_correctors_list, quad_dict):
    # elements_name = used_correctors_list
    correctors_indexes = []
    quad_dict_ = []
    elements_name = []
    quads = pd.DataFrame()
    s_pos =[]
    for i in used_correctors_list:

        # quad_dict_.append(int(quad_dict[i]))
        quad_dict_ = int(quad_dict[i])
        elements_numbers = quad_dict_
        corrector_indexx = get_refpts(ring, i)
        # print(corrector_index)
        element_name = ring[corrector_indexx[0]].FamName
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=corrector_indexx)
        s_poss = lindata['s_pos']
        s_pos.append(s_poss)

        df1 = {
            str(i) + str("=") + str(" ") + str(quad_dict_) + str(" ") + str('quads'): corrector_indexx,
        }
        df2 = pd.concat([pd.DataFrame(v, columns=[k]) for k, v in df1.items()], axis=1)


        quads = pd.concat([quads, df2], axis=1)
        for j in range(len(s_pos)):
            array1 = numpy.append(s_pos[0], s_pos[j])


    return quads, s_pos

def used_quads_f(ring, used_correctors_list, quad_dict):

    #elements_name = used_correctors_list
    correctors_indexes = []
    quad_dict_ = []
    elements_name =[]
    quads = pd.DataFrame()
    for i in used_correctors_list:

        #quad_dict_.append(int(quad_dict[i]))
        quad_dict_= int(quad_dict[i])
        elements_numbers = quad_dict_
        corrector_index = get_refpts(ring, i)
        #print(corrector_index)
        element_name = ring[corrector_index[0]].FamName
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=corrector_index)
        s_poss = lindata['s_pos']
        #print(element_name)


        df1 = {
                    str(i) + str("=") + str(" ")+ str( quad_dict_)+ str(" ")+ str('quads'): s_poss,
                }
        df2 = pd.concat([pd.DataFrame(v, columns=[k]) for k, v in df1.items()], axis=1)


        correctors_indexes.append(np.squeeze(corrector_index))
        elements_name.append(element_name)
        quads = pd.concat([quads, df2], axis=1)


    return quads


def used_elements_plot(lattice, used_quad):
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



    plt.plot(s_pos, closed_orbitx)
    #plt.plot(s_pos, closed_orbitx)

    # Label for x-axis
    plt.xlabel("elements_indexes")

    # Label for y-axis
    plt.ylabel("closed_orbit_x")

      # for display


    i = 0
    S_pos2 = []
    while (i < used_quad.shape[1]):
        S_pos1 = used_quad.iloc[:, i]
        S_pos_ = df = pd.concat([S_pos1])
        S_pos2.append(S_pos_)
        i += 1
        for i in S_pos_:
            scatter(i, 0)
    plt.title("used quadrupoles indices")
    plt.show()



    plt.plot(s_pos, beta_x)
    #plt.plot(s_pos, beta_x)
    # Label for x-axis
    plt.xlabel("elements_indexes")

    # Label for y-axis
    plt.ylabel("beta_x")

     # for display

    S_pos2 = []
    i = 0
    S_pos2 = []
    while (i < used_quad.shape[1]):
        S_pos1 = used_quad.iloc[:, i]
        S_pos_ = df = pd.concat([S_pos1])
        S_pos2.append(S_pos_)
        i += 1
        for i in S_pos_:
            scatter(i, 0)
    plt.title("used quadrupoles indices")
    plt.show()



def used_correctors_f(lattice, used_correctors_list):

    elements_name = used_correctors_list
    correctors_indexes = []
    quad_dict = []
    for i in used_correctors_list:
        # print(i)
        corrector_index = get_refpts(lattice, i)
        correctors_indexes.append(np.squeeze(corrector_index))


        j = 0
        s_pos=[]
        Ind= []
        while (j < len( corrector_index)):
            corrector_indexx = corrector_index[j]

            lindata0, tune, chrom, lindata = lattice.linopt(get_chrom=True, refpts=correctors_indexes)
            s_poss = lindata['s_pos']
            #s_pos.append(s_poss)
            s_pos.append(np.squeeze(s_poss))
            Ind.append(corrector_indexx)
            df1 = {'Used elements names': elements_name, 'S_pos': correctors_indexes}
            j += 1


    return df1, s_pos, correctors_indexes



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




def ORM_x_eta_w(dkick, ring, ind):
    cxx = []
    cxy = []
    c_indexes =[]

    #correctors_indexes = get_refpts(ring, used_correctors)
    bpm_indexes = get_refpts(ring, elements.Monitor)

    for j in ind:
        ring[j].KickAngle[0] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cxx.append(closed_orbitx)
        cxy.append(closed_orbity)

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


def ORM_y_eta_w(dkick, ring, ind):
    cyy = []
    cyx = []
    c_indexes = []
    #correctors_indexes = get_refpts(ring, used_correctors)
    bpm_indexes = get_refpts(ring, elements.Monitor)

    for j in ind:
        ring[j].KickAngle[1] = dkick
        lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 2]

        cyy.append(closed_orbity)
        cyx.append(closed_orbitx)

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


### Dispersion_ORM


def eta(dkick, ring):
    bpm_indexes = get_refpts(ring, elements.Monitor)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    Eta_xx = lindata['dispersion'][:, 0]
    Eta_yy = lindata['dispersion'][:, 2]

    return Eta_xx,Eta_yy



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
    #quads.to_csv("C:/Users/musa/pyat-loco-1/fodo_loco/mydata/quad_info.csv")

    #print('Done...')

    return quads


def simulateErrorQFQD(lattice, errorQ, tiltQ, shiftQx,shiftQy, debug=False):

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


        if lattice[quad_indexes[i]].FamName == 'QF':
            lattice[quad_indexes[i]].K *= (1 + errorQ * np.random.randn())
            at.tilt_elem(lattice[quad_indexes[i]], tiltQ * np.random.randn(), relative=False)
            at.shift_elem(lattice[quad_indexes[i]], shiftQx * np.random.randn(), shiftQy * np.random.randn(),
                          relative=False)
        #if debug: print(
         #   f"sumulateError: qaudrupole name {lattice[quad_indexes[i]].FamName} , #{i} out of {len(quad_indexes)}")
        if lattice[quad_indexes[i]].FamName == 'QD':
            lattice[quad_indexes[i]].K *= (1 + errorQ * np.random.randn())
            at.tilt_elem(lattice[quad_indexes[i]], tiltQ * np.random.randn(), relative=False)
            at.shift_elem(lattice[quad_indexes[i]], shiftQx * np.random.randn(), shiftQy * np.random.randn(),
                          relative=False)
            #at.tilt_elem(lattice[quad_indexes[i]], tiltQ , relative=False)
            #at.shift_elem(lattice[quad_indexes[i]], shiftQx,shiftQx, relative=False)



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

def simulateErrorTest(qname, lattice, errorQ, tiltQ, shiftQx,shiftQy, debug=False):

    quad_indexes = get_refpts(lattice, qname)
    np.random.seed(1)

    Quad_strength = []
    Quad_strength_err = []
    elements_name_q = []

    quad_indexes = get_refpts(lattice, qname)
    for i in quad_indexes:
        quad_strength = lattice[i].K
        Quad_strength.append(quad_strength)
        element_name = lattice[i].FamName
        elements_name_q.append(element_name)


        lattice[i].K *= (1 + errorQ * np.random.randn())
        np.random.seed(1)
        at.tilt_elem(lattice[i], tiltQ * np.random.randn(), relative=False)
        np.random.seed(1)
        at.shift_elem(lattice[i], shiftQx * np.random.randn(), shiftQy * np.random.randn(),relative=False)
    Quad_strength_err =[]


    for j in quad_indexes:
        quad_strength_err = lattice[j].K
        Quad_strength_err.append(quad_strength_err)



    opt = at.linopt(lattice, refpts=quad_indexes, get_chrom=True)
    s_pos_q = opt[3].s_pos

    output1 = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 1)]


    quad = {'s_pos': s_pos_q, 'Quad_strength': Quad_strength_err,
            'elements_name': elements_name_q, 'occ': output1}
    quads = pd.DataFrame(quad)




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
        #at.tilt_elem(lattice[quad_indexes[i]], tiltQ , relative=False)
        #at.shift_elem(lattice[quad_indexes[i]], shiftQx,shiftQx, relative=False)



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


def computeOpticsD1(ring, qname, i, dk, quad_vals,used_correctors):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    quad_indexes = get_refpts(ring, qname)

    ring[quad_indexes[i]].K = quad_vals[qname,i] + dk

    qxx, qxy = ORM_x_eta_w(dk, ring, used_correctors)
    qyy, qyx = ORM_y_eta_w(dk, ring, used_correctors)

    ring[quad_indexes[i]].K = quad_vals[qname,i]


    return  qxx, qxy, qyy, qyx



def generatingQuadsResponse(ring, Cxx, Cyy,Cxy, Cyx , used_quads,used_correctors ):
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
            #for i in range(0, nq):
            Qxx, Qxy, Qyy, Qyx = computeOpticsD(ring, qname, dk, quad_vals, used_correctors)
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
    #for qname in quad_names:
        # nquad = quad_dict[qname]
        #print('loading response to:', qname)
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





def computeOpticsD(ring, qname, dk, quad_vals,used_correctors):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    quad_indexes = get_refpts(ring, qname)
    for i in quad_indexes:
        a = ring[i].K
        #print('a',a)
        ring[i].K = a + dk
        b = ring[i].K
        #print('b',b)
    qxx, qxy = ORM_x_eta_w(dk, ring, used_correctors)
    qyy, qyx = ORM_y_eta_w(dk, ring, used_correctors)

    for i in quad_indexes:
        ring[i].K = a
        c = ring[i].K
        #print('c',c)
    return  qxx, qxy, qyy, qyx

def computeOptics_eta(ring, qname, i, dk, quad_vals):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    quad_indexes = get_refpts(ring, qname)

    ring[quad_indexes[i]].K = quad_vals[qname,i] + dk

    qx, qy = eta(dk, ring)
    #qyy, qyx = eta(dk, ring)

    ring[quad_indexes[i]].K = quad_vals[qname,i]


    return  qx, qy



def getBetaBeat_plot(betax, betax_error,betay, betay_error, s_pos):
    print("getBetaBeat bx and by: ")

    bxi =[]
    for i in range(len(betax)):
        bxx = (betax_error[i] - betax[i]) / betax[i]
        bxi.append(bxx)

    byi =[]
    for i in range(len(betay)):
        byy = (betay_error[i] - betay[i]) / betay[i]
        byi.append(byy)

    bx = np.std((betax_error - betax) / betax)
    by = np.std((betay_error - betay) / betay)
    print("Simulated beta beat, x:" + str(bx * 100) + "%   y: " + str(by* 100) + "%")
    plt.figure()
    plt.plot(s_pos, bxi)
    plt.xlabel('s(m)')
    plt.ylabel(r'$\Delta \beta_x / \beta_x$')
    plt.show()

    plt.figure()
    plt.plot(s_pos, byi)
    plt.xlabel('s(m)')
    plt.ylabel(r'$\Delta \beta_y / \beta_y$')
    plt.show()