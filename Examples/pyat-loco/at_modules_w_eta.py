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


def defineMatrices_w_eta_old(W, alpha1, alpha2, C0x, C0y, C0xy, C0yx, Cxx_err, Cyy_err, Cxy_err, Cyx_err, dCx, dCy, dCxy,dCyx, eta_0x, eta_0y, eta_xx_err, eta_yy_err, d_eta_x, d_eta_y ):
    Nk = len(dCx)  # number of free parameters
    Nm = len(dCx) # number of measurements
    print('NK:', Nk)
    print('Nm:', Nm)

    Ax = np.zeros([Nk, Nk])
    Ay = np.zeros([Nk, Nk])
    Axy = np.zeros([Nk, Nk])
    Ayx = np.zeros([Nk, Nk])
    Ax_eta = np.zeros([Nk, Nk])
    Ay_eta = np.zeros([Nk, Nk])
    Axy_eta = np.zeros([Nk, Nk])
    Ayx_eta = np.zeros([Nk, Nk])

    A = np.zeros([4 * Nk, Nk])
    A_eta = np.zeros([4 * Nk, Nk])

    ##

    Bx = np.zeros([Nk, 1])
    By = np.zeros([Nk, 1])
    Bxy = np.zeros([Nk, 1])
    Byx = np.zeros([Nk, 1])
    Bx_eta = np.zeros([Nk, 1])
    By_eta= np.zeros([Nk, 1])
    Bxy_eta = np.zeros([Nk, 1])
    Byx_eta = np.zeros([Nk, 1])

    B = np.zeros([4 * Nk, 1])
    B_eta = np.zeros([4 * Nk, 1])
    ##

    Dx = (Cxx_err[:, :] - C0x[:, :] )#- error_variance)  ### dk ?
    Dy = (Cyy_err[:, :] - C0y[:, :] )
    Dxy = (Cxy_err[:, :] - C0xy[:, :])
    Dyx = (Cyx_err[:, :] - C0yx[:, :] )
    Dx_eta = (eta_xx_err - eta_0x)  # - error_variance)  ### dk ?
    Dy_eta = (eta_yy_err - eta_0y)

    ##

    for i in range(Nk):  ## i represents each quad
        # print('done A:', 100.* i ,'%')
        for j in range(Nk):
            Ax[i, j] = np.sum(np.dot(np.dot(dCx[i],W*alpha1), dCx[j].T))
            Ay[i, j] = np.sum(np.dot(np.dot(dCy[i],W*alpha1), dCy[j].T))
            Axy[i, j] = np.sum(np.dot(np.dot(dCxy[i],W*alpha1), dCxy[j].T))
            Ayx[i, j] = np.sum(np.dot(np.dot(dCyx[i],W*alpha1), dCyx[j].T))
        A[i, :] = Ax[i, :]
        A[i + Nk, :] = Ay[i, :]
        A[i + 2 * Nk, :] = Axy[i, :]
        A[i + 3 * Nk, :] = Ayx[i, :]

    for i in range(Nk):  ## i represents each quad
        # print('done A:', 100.* i ,'%')
        for j in range(Nk):
            Ax_eta[i, j] = np.sum(np.dot(np.dot(d_eta_x[i], W*alpha2), d_eta_x[j].T))
            Ay_eta[i, j] = np.sum(np.dot(np.dot(d_eta_y[i], W*alpha2), d_eta_y[j].T))

        A_eta[i, :] = Ax_eta[i, :]
        A_eta[i + Nk, :] = Ay_eta[i, :]


    ##

    for i in range(Nk):
        Bx[i] = np.sum(np.dot(np.dot(dCx[i],W*alpha1), Dx.T))
        By[i] = np.sum(np.dot(np.dot(dCy[i],W*alpha1), Dy.T))
        Bxy[i] = np.sum(np.dot(np.dot(dCxy[i],W*alpha1), Dxy.T))
        Byx[i] = np.sum(np.dot(np.dot(dCyx[i],W*alpha1), Dyx.T))
        B[i] = Bx[i]
        B[i + Nk] = By[i]
        B[i + 2 * Nk] = Bxy[i]
        B[i + 3 * Nk] = Byx[i]
    for i in range(Nk):
        Bx_eta[i] = np.sum(np.dot(np.dot(d_eta_x[i], W * alpha2), Dx_eta.T))
        By_eta[i] = np.sum(np.dot(np.dot(d_eta_y[i], W * alpha2), Dy_eta.T))

        B_eta[i] = Bx_eta[i]
        B_eta[i + Nk] = By_eta[i]

    return A, B, A_eta, B_eta


def defineMatrices_w_eta(W, alpha1, alpha2, C0x, C0y, C0xy, C0yx, Cxx_err, Cyy_err, Cxy_err, Cyx_err, dCx, dCy, dCxy,dCyx):
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
            Ax[i, j] = np.sum(np.dot(np.dot(dCx[i],W*alpha1), dCx[j].T))
            Ay[i, j] = np.sum(np.dot(np.dot(dCy[i],W*alpha1), dCy[j].T))
            Axy[i, j] = np.sum(np.dot(np.dot(dCxy[i],W*alpha1), dCxy[j].T))
            Ayx[i, j] = np.sum(np.dot(np.dot(dCyx[i],W*alpha1), dCyx[j].T))
        A[i, :] = Ax[i, :]
        A[i + Nk, :] = Ay[i, :]
        A[i + 2 * Nk, :] = Axy[i, :]
        A[i + 3 * Nk, :] = Ayx[i, :]


    ##

    for i in range(Nk):
        Bx[i] = np.sum(np.dot(np.dot(dCx[i],W*alpha1), Dx.T))
        By[i] = np.sum(np.dot(np.dot(dCy[i],W*alpha1), Dy.T))
        Bxy[i] = np.sum(np.dot(np.dot(dCxy[i],W*alpha1), Dxy.T))
        Byx[i] = np.sum(np.dot(np.dot(dCyx[i],W*alpha1), Dyx.T))
        B[i] = Bx[i]
        B[i + Nk] = By[i]
        B[i + 2 * Nk] = Bxy[i]
        B[i + 3 * Nk] = Byx[i]


    return A, B,





def defineMatrices_w_eta_t(W, alpha1, alpha2, C0x, C0y, C0xy, C0yx, Cxx_err, Cyy_err, Cxy_err, Cyx_err, dCx, dCy, dCxy,dCyx, eta_0x, eta_0y,eta_0xy, eta_0yx, eta_xx_err, eta_yy_err,eta_xy_err,eta_yx_err, d_eta_x, d_eta_y,d_eta_xy ,d_eta_yx ):
    Nk = len(dCx)  # number of free parameters
    Nm = len(dCx) # number of measurements
    print('NK:', Nk)
    print('Nm:', Nm)

    Ax = np.zeros([Nk, Nk])
    Ay = np.zeros([Nk, Nk])
    Axy = np.zeros([Nk, Nk])
    Ayx = np.zeros([Nk, Nk])
    Ax_eta = np.zeros([Nk, Nk])
    Ay_eta = np.zeros([Nk, Nk])
    Axy_eta = np.zeros([Nk, Nk])
    Ayx_eta = np.zeros([Nk, Nk])

    A = np.zeros([4 * Nk, Nk])
    A_eta = np.zeros([4 * Nk, Nk])

    ##

    Bx = np.zeros([Nk, 1])
    By = np.zeros([Nk, 1])
    Bxy = np.zeros([Nk, 1])
    Byx = np.zeros([Nk, 1])
    Bx_eta = np.zeros([Nk, 1])
    By_eta= np.zeros([Nk, 1])
    Bxy_eta = np.zeros([Nk, 1])
    Byx_eta = np.zeros([Nk, 1])

    B = np.zeros([4 * Nk, 1])
    B_eta = np.zeros([4 * Nk, 1])
    ##

    Dx = (Cxx_err[:, :] - C0x[:, :] )#- error_variance)  ### dk ?
    Dy = (Cyy_err[:, :] - C0y[:, :] )
    Dxy = (Cxy_err[:, :] - C0xy[:, :])
    Dyx = (Cyx_err[:, :] - C0yx[:, :] )
    Dx_eta = (eta_xx_err - eta_0x)  # - error_variance)  ### dk ?
    Dy_eta = (eta_yy_err - eta_0y)
    Dxy_eta = (eta_xy_err - eta_0xy)  # - error_variance)  ### dk ?
    Dyx_eta = (eta_yx_err - eta_0yx)
    ##

    for i in range(Nk):  ## i represents each quad
        # print('done A:', 100.* i ,'%')
        for j in range(Nk):
            Ax[i, j] = np.sum(np.dot(np.dot(dCx[i],W*alpha1), dCx[j].T))
            Ay[i, j] = np.sum(np.dot(np.dot(dCy[i],W*alpha1), dCy[j].T))
            Axy[i, j] = np.sum(np.dot(np.dot(dCxy[i],W*alpha1), dCxy[j].T))
            Ayx[i, j] = np.sum(np.dot(np.dot(dCyx[i],W*alpha1), dCyx[j].T))
        A[i, :] = Ax[i, :]
        A[i + Nk, :] = Ay[i, :]
        A[i + 2 * Nk, :] = Axy[i, :]
        A[i + 3 * Nk, :] = Ayx[i, :]

    for i in range(Nk):  ## i represents each quad
        # print('done A:', 100.* i ,'%')
        for j in range(Nk):
            Ax_eta[i, j] = np.sum(np.dot(np.dot(d_eta_x[i], W*alpha2), d_eta_x[j].T))
            Ay_eta[i, j] = np.sum(np.dot(np.dot(d_eta_y[i], W*alpha2), d_eta_y[j].T))
            Axy_eta[i, j] = np.sum(np.dot(np.dot(d_eta_xy[i], W * alpha2), d_eta_xy[j].T))
            Ayx_eta[i, j] = np.sum(np.dot(np.dot(d_eta_yx[i], W * alpha2), d_eta_yx[j].T))
        A_eta[i, :] = Ax_eta[i, :]
        A_eta[i + Nk, :] = Ay_eta[i, :]
        A_eta[i + 2 * Nk, :] = Axy_eta[i, :]
        A_eta[i + 3 * Nk, :] = Ayx_eta[i, :]

    ##

    for i in range(Nk):
        Bx[i] = np.sum(np.dot(np.dot(dCx[i],W*alpha1), Dx.T))
        By[i] = np.sum(np.dot(np.dot(dCy[i],W*alpha1), Dy.T))
        Bxy[i] = np.sum(np.dot(np.dot(dCxy[i],W*alpha1), Dxy.T))
        Byx[i] = np.sum(np.dot(np.dot(dCyx[i],W*alpha1), Dyx.T))
        B[i] = Bx[i]
        B[i + Nk] = By[i]
        B[i + 2 * Nk] = Bxy[i]
        B[i + 3 * Nk] = Byx[i]
    for i in range(Nk):
        Bx_eta[i] = np.sum(np.dot(np.dot(d_eta_x[i], W * alpha2), Dx_eta.T))
        By_eta[i] = np.sum(np.dot(np.dot(d_eta_y[i], W * alpha2), Dy_eta.T))
        Bxy_eta[i] = np.sum(np.dot(np.dot(d_eta_xy[i], W * alpha2), Dxy_eta.T))
        Byx_eta[i] = np.sum(np.dot(np.dot(d_eta_yx[i], W * alpha2), Dyx_eta.T))


        B_eta[i] = Bx_eta[i]
        B_eta[i + Nk] = By_eta[i]
        B_eta[i + 2 * Nk] = Bxy_eta[i]
        B_eta[i + 3 * Nk] = Byx_eta[i]
    return A, B, A_eta, B_eta













def defineMatrices(C0x, C0y, C0xy, C0yx, Cxx_err, Cyy_err, Cxy_err, Cyx_err, dCx, dCy, dCxy,dCyx ):
    Nk = len(dCx)  # number of free parameters
    Nm = len(dCx)   # number of measurements
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


def  defineMatrices_eta_w(W, C0x, C0y, Cxx_err, Cyy_err, dCx, dCy):
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
            #Ax[i, j] = np.sum(np.dot(dCx[i], dCx[j].T))
            #Ay[i, j] = np.sum(np.dot(dCy[i], dCy[j].T))
            Ax[i, j] = np.sum(np.dot(np.dot(dCx[i], W), dCx[j].T))
            Ay[i, j] = np.sum(np.dot(np.dot(dCy[i], W), dCy[j].T))

        A[i, :] = Ax[i, :]
        A[i + Nk, :] = Ay[i, :]


    ##

    for i in range(Nk):
        #Bx[i] = np.sum(np.dot(dCx[i], Dx.T))
        #By[i] = np.sum(np.dot(dCy[i], Dy.T))
        Bx[i] = np.sum(np.dot(np.dot(dCx[i], W), Dx.T))
        By[i] = np.sum(np.dot(np.dot(dCy[i], W), Dy.T))

        B[i] = Bx[i]
        B[i + Nk] = By[i]


    return A, B





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




def generatingQuadsResponse_weta(ring, Cxx, Cyy,Cxy, Cyx ):
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
        nq = quad_dict[qname] + 1
        for i in range(0, nq):
            Qxx, Qxy, Qyy, Qyx = computeOpticsD_weta(ring, qname, i, dk, quad_vals)
            qxx.append(Qxx)
            qxy.append(Qxy)
            qyy.append(Qyy)
            qyx.append(Qyx)

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



def generatingQuadsResponse(ring, Cxx, Cyy,Cxy, Cyx ):
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
        nq = quad_dict[qname] + 1
        for i in range(0, nq):
            Qxx, Qxy, Qyy, Qyx = computeOpticsD(ring, qname, i, dk, quad_vals)
            qxx.append(Qxx)
            qxy.append(Qxy)
            qyy.append(Qyy)
            qyx.append(Qyx)

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





def generatingQuadsResponse_eta(ring, Cxx, Cyy):
    # %%time

    quads_info = quad_info(ring)
    quad_dict, quad_vals = getQuadFamilies(quads_info)
    quads = [k for k in quad_dict.keys()]
    quad_names = quads
    dk = 0.0001
    qxx = []
    qyy = []

    quad_names = quads
    for qname in quad_names:
        print('generating response to {}, n={}'.format(qname, quad_dict[qname]))
        nq = quad_dict[qname] + 1
        for i in range(0, nq):
            Qxx, Qyy= computeOptics_eta(ring, qname, i, dk, quad_vals)
            qxx.append(Qxx)
            qyy.append(Qyy)

    C0x = Cxx
    C0y = Cyy


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
        dcxx = ((C1x - C0x) / dk)
        dcyy = ((C1y - C0y) / dk)



        dCx.append(dcxx)
        dCy.append(dcyy)

        i += 1

    return C0x, C0y, dCx, dCy

def generatingQuadsResponse_eta_t(ring, Cxx, Cyy, Cxy,Cyx):
    # %%time

    quads_info = quad_info(ring)
    quad_dict, quad_vals = getQuadFamilies(quads_info)
    quads = [k for k in quad_dict.keys()]
    quad_names = quads
    dk = 0.0001
    qxx = []
    qyy = []
    qxy = []
    qyx = []
    quad_names = quads
    for qname in quad_names:
        print('generating response to {}, n={}'.format(qname, quad_dict[qname]))
        nq = quad_dict[qname] + 1
        for i in range(0, nq):
            Qxx, Qxy, Qyy,Qyx= computeOptics_eta_t(ring, qname, i, dk, quad_vals)
            qxx.append(Qxx)
            qyy.append(Qyy)
            qxy.append(Qxy)
            qyx.append(Qyx)

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
        dcxy = ((C1xy - C0xy) / dk)
        dcyx = ((C1yx - C0yx) / dk)


        dCx.append(dcxx)
        dCy.append(dcyy)
        dCxy.append(dcxy)
        dCyx.append(dcyx)

        i += 1

    return C0x, C0y, C0xy,C0yx, dCx, dCy,dCxy,dCyx






def setCorrection(ring, quads_info_error,quad_names, r , quads_info,n_list):

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

        if qname_ == 'QF':
            occ = quads_info_error.occ[idx]
            dk = cor_dict['QF', occ]
            DK.append(dk)

        if qname_ == 'QD':
            occ = quads_info_error.occ[idx]
            dk = cor_dict['QD', occ]
            DK.append(dk)

        if qname_ == 'QS':
            occ = quads_info_error.occ[idx]
            dk = cor_dict['QS', occ]
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

    bx = np.std((twiss_error.betax - twiss.betax) / twiss.betax)
    by = np.std((twiss_error.betay - twiss.betay) / twiss.betay)
    print("Simulated beta beat, x:" + str(bx * 100) + "%   y: " + str(by* 100) + "%")






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




def ORM_x_eta_w(dkick, ring):
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


def ORM_y_eta_w(dkick, ring):
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

def eta_t(dkick, ring):
    bpm_indexes = get_refpts(ring, elements.Monitor)
    lindata0, tune, chrom, lindata = ring.linopt(get_chrom=True, refpts=bpm_indexes)
    Eta_xx = lindata['dispersion'][:, 0]
    Eta_xy = lindata['dispersion'][:, 1]
    Eta_yy = lindata['dispersion'][:, 2]
    Eta_yx = lindata['dispersion'][:, 1]
    return Eta_xx, Eta_xy, Eta_xy, Eta_yx


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





def quad_info_qs(lattice):

    Quad_strength = []
    Quad_strength_err = []
    elements_name_q = []

    quad_indexes = get_refpts(lattice, 'QS')

    for i in quad_indexes:
        quad_strength = lattice[i].K
        Quad_strength.append(quad_strength)
        element_name = lattice[i].FamName
        elements_name_q.append(element_name)

    elements_indexes = get_refpts(lattice, 'QS')
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






def simulateError1(lattice, errorQF, errorQD, tiltQF, tiltQD,shiftQF, shiftQD):

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
            at.shift_elem(lattice[quad_indexes[i]], tiltQF, relative=False)

            i += 1

        elif (lattice[quad_indexes[i]].FamName == 'QD'):
            lattice[quad_indexes[i]].K *= (1 + errorQD * random())
            at.tilt_elem(lattice[quad_indexes[i]], tiltQD, relative=False)
            at.shift_elem(lattice[quad_indexes[i]], shiftQD, relative=False)

        i += 1

       # elif (lattice[quad_indexes[i]].FamName == 'QS'):
       #     lattice[quad_indexes[i]].K *= (1 + errorQS * random())
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
    #quads.to_csv("C:/Users/musa/pyat-loco-1/fodo_loco/mydata/quad_info_error.csv")


    #print('Done...')

    return quads

def simulateError(lattice, errorQF, errorQD, errorQS, tiltQF, tiltQD, shiftQF, shiftQD):

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
            at.shift_elem(lattice[quad_indexes[i]], shiftQF, relative=False)

            i += 1

        elif (lattice[quad_indexes[i]].FamName == 'QD'):
            lattice[quad_indexes[i]].K *= (1 + errorQD * random())
            at.tilt_elem(lattice[quad_indexes[i]], tiltQD, relative=False)
            at.shift_elem(lattice[quad_indexes[i]], shiftQD, relative=False)

            i += 1

        elif (lattice[quad_indexes[i]].FamName == 'QS'):
            lattice[quad_indexes[i]].K += (errorQS * random())
            #at.tilt_elem(lattice[quad_indexes[i]], tiltQS, relative=False)

            i += 1




    for j in quad_indexes:
        quad_strength_err = lattice[j].K
        Quad_strength_err.append(quad_strength_err)

   # elements_indexes = get_refpts(lattice, 'QS')
    opt = at.linopt(lattice, refpts=quad_indexes, get_chrom=True)
    s_pos_q = opt[3].s_pos

    output1 = [elements_name_q[:e].count(v) for e, v in enumerate(elements_name_q, 1)]


    quad = {'s_pos': s_pos_q, 'Quad_strength': Quad_strength_err,
            'elements_name': elements_name_q, 'occ': output1}
    quads = pd.DataFrame(quad)
    #quads.to_csv("C:/Users/musa/pyat-loco-1/fodo_loco/mydata/quad_info_error.csv")


    #print('Done...')

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

def computeOpticsD_weta(ring, qname, i, dk, quad_vals):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    quad_indexes = get_refpts(ring, qname)

    ring[quad_indexes[i]].K = quad_vals[qname,i] + dk

    qxx, qxy = ORM_x_eta_w(dk, ring)
    qyy, qyx = ORM_y_eta_w(dk, ring)

    ring[quad_indexes[i]].K = quad_vals[qname,i]


    return  qxx, qxy, qyy, qyx



def computeOpticsD(ring, qname, i, dk, quad_vals):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    quad_indexes = get_refpts(ring, qname)

    ring[quad_indexes[i]].K = quad_vals[qname,i] + dk

    qxx, qxy = ORM_x(dk, ring)
    qyy, qyx = ORM_y(dk, ring)

    ring[quad_indexes[i]].K = quad_vals[qname,i]


    return  qxx, qxy, qyy, qyx


def computeOptics_eta(ring, qname, i, dk, quad_vals):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    quad_indexes = get_refpts(ring, qname)

    ring[quad_indexes[i]].K = quad_vals[qname,i] + dk

    qx, qy = eta(dk, ring)
    #qyy, qyx = eta(dk, ring)

    ring[quad_indexes[i]].K = quad_vals[qname,i]


    return  qx, qy


def computeOptics_eta_t(ring, qname, i, dk, quad_vals):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    quad_indexes = get_refpts(ring, qname)

    ring[quad_indexes[i]].K = quad_vals[qname,i] + dk

    qx, qxy,qy,qyx = eta_t(dk, ring)
    #qyy, qyx = eta(dk, ring)

    ring[quad_indexes[i]].K = quad_vals[qname,i]


    return  qx, qxy,qy,qyx
def computeOpticsQS(ring, qname, i, dk, qs_indexes):

    bpm_indexes = get_refpts(ring, elements.Monitor)
    #print('before '+str(ring[qs_indexes[i]].K))
    Qs_old_strength = ring[qs_indexes[i]].K
    ring[qs_indexes[i]].K = ring[qs_indexes[i]].K + dk
    #print('after '+str(ring[qs_indexes[i]].K))
    qxx, qxy = ORM_x(dk, ring)
    qyy, qyx = ORM_y(dk, ring)

    ring[qs_indexes[i]].K = Qs_old_strength
    #print('final qs ' +str(ring[qs_indexes[i]].K))


    return  qxx, qxy, qyy, qyx


def getBetaBeat_t(betax, betax_error,betay, betay_error, s_pos):
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