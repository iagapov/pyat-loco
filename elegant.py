"""
elegant post-/preprocessing
"""

import sys, time, os
import sdds
import subprocess
from pylab import *

elegant="/Users/ia/Products/elegant/darwin-x86/elegant"


def getTwiss(idx=0,fname='twiss'):
     x = sdds.SDDS(0)
     x.load(fname + '.twi')
     s = x.columnData[0][idx]
     betaX = np.array(x.columnData[1][idx])
     alphaX = np.array(x.columnData[2][idx])
     betaY = np.array(x.columnData[7][idx])
     alphaY = np.array(x.columnData[8][idx])
     etaX = np.array(x.columnData[4][idx])

     res = {'s':s, 'beta_x': betaX, 'alpha_x': alphaX, 'beta_y':betaY, 'alpha_y':alphaY, 'eta_x' : etaX }
     return res


def plotOptics(idx=0, fname='twiss', limits=None, showAlpha=False, plotTitle=""):
     x = sdds.SDDS(0)
     x.load(fname + '.twi')
     #print(x.columnName)
     #print(x.columnName[0], x.columnName[1], x.columnName[6], x.columnName[7])
     #print(x.columnData[0])
     ax = plt.figure().add_subplot(111)
     ax.set_title(plotTitle)
     s = x.columnData[0][idx]
     betaX = np.array(x.columnData[1][idx])
     alpha_x = np.array(x.columnData[2][idx])
     betaY = np.array(x.columnData[7][idx])
     alpha_y = np.array(x.columnData[8][idx])
     etaX = np.array(x.columnData[4][idx])
     if limits == None:
         p1,=ax.plot(s, betaX, 'b-', lw=3)
         p2,=ax.plot(s, betaY, 'r-', lw=3)
         ax2 = ax.twinx()
         if showAlpha:
             p3,=ax2.plot(s, alpha_x, 'b--', lw=3)
             p4,=ax2.plot(s, alpha_y, 'r--', lw=3)
         else:
             p3,=ax2.plot(s, etaX, 'g-', lw=3)
     else:
         i1 = 0
         i2 = len(s)-1
         for i in range(1,len(s)):
             if s[i-1] < limits[0] and s[i] >= limits[0]: i1 = i
             if s[i-1] < limits[1] and s[i] >= limits[1]: i2 = i

         p1,=ax.plot(s[i1:i2], betaX[i1:i2], 'b-', lw=3)
         p2,=ax.plot(s[i1:i2], betaY[i1:i2], 'r-', lw=3)
         ax2 = ax.twinx()
         if showAlpha:
            p3,=ax2.plot(s[i1:i2], alpha_x[i1:i2], 'b--', lw=3)
            p4,=ax2.plot(s[i1:i2], alpha_y[i1:i2], 'r--', lw=3)
         else:
             p3,=ax2.plot(s[i1:i2], etaX[i1:i2], 'g-', lw=3)
     if showAlpha:
         ax.legend([p1,p2,p3,p4],['$\\beta_x$','$\\beta_y$','$\\alpha_x$','$\\alpha_y$'])
     else:
         ax.legend([p1,p2,p3],['$\\beta_x$','$\\beta_y$','$\\eta_x$'])
     plt.show()

def plotCoupledOptics(idx=0,fname='twiss', limits=None, plotTitle=""):
     f = sdds.SDDS(0)
     f.load(fname + '.ctwi')
     ax = plt.figure().add_subplot(111)
     ax.set_title(plotTitle)
     s = np.array(x.columnData[1][idx])
     betaX1 = np.array(f.columnData[8][idx]) # TODO: use some index finding function
     betaX2 = np.array(f.columnData[9][idx])
     betaY1 = np.array(f.columnData[10][idx])
     betaY2 = np.array(f.columnData[11][idx])
     etaX = np.array(f.columnData[12][idx])

     if limits == None:
         p1,=ax.plot(s, betaX1, 'b-', lw=3)
         p2,=ax.plot(s, betaY1, 'r--', lw=3)
         ax.plot(s, betaX2, 'b--', lw=2)
         ax.plot(s, betaY2, 'r-', lw=2)
         ax2 = ax.twinx()
         p3,=ax2.plot(s, etaX, 'g-', lw=3)
     else:
         i1 = 0
         i2 = len(s)-1
         print(i2)
         for i in range(1,len(s)):
             if s[i-1] < limits[0] and s[i] >= limits[0]: i1 = i
             if s[i-1] < limits[1] and s[i] >= limits[1]: i2 = i

         p1,=ax.plot(s[i1:i2], betaX1[i1:i2], 'b-', lw=3)
         p2,=ax.plot(s[i1:i2], betaY1[i1:i2], 'r-', lw=3)
         ax.plot(s[i1:i2], betaX1[i1:i2], 'b--', lw=2)
         ax.plot(s[i1:i2], betaY1[i1:i2], 'r--', lw=2)
         ax2 = ax.twinx()
         p3,=ax2.plot(s[i1:i2], etaX[i1:i2], 'g-', lw=3)
     ax.legend([p1,p2,p3],['$\\beta_x$','$\\beta_y$','$\\eta_x$'])
     plt.show()

def getSdict(twissFile):
    #print(x.columnName[0],x.columnName[14],x.columnName[15],x.columnName[16])

    n_list = len(twissFile.columnData[0][0])

    s_dict = {}

    for idx in range(n_list):
        s = twissFile.columnData[0][0][idx]
        ename = twissFile.columnData[14][0][idx]
        eocc = twissFile.columnData[15][0][idx]
        etype = twissFile.columnData[16][0][idx]

        s_dict[ename, eocc] = s
    return s_dict

def showCorrectors(ipage=0,fnameCalc='twiss', fname='twiss'):
    '''
    fnameCalc contains element s-positions
    fname contains parameter table of quadrupole settings
    '''
    twissFile = sdds.SDDS(0)
    twissFile.load(fnameCalc + '.twi')

    s_dict = getSdict(twissFile)

    x2 = sdds.SDDS(0)
    x2.load(fname+'.param')

    #['ElementName', 'ElementParameter', 'ParameterValue', 'ElementType', 'ElementOccurence', 'ElementGroup']
    #print(x2.columnName[0], x2.columnName[1], x2.columnName[2], x2.columnName[3], x2.columnName[4])
    n_list = len(x2.columnData[0][ipage])

    s_a = []
    hc_a = []
    vc_a = []

    for idx in range(n_list):
        ename = x2.columnData[0][ipage][idx] # ElementName
        par = x2.columnData[1][ipage][idx] #ElementParameter
        eval = x2.columnData[2][ipage][idx] #ParameterValue
        etype = x2.columnData[4][ipage][idx] #ElementType
        eocc = x2.columnData[5][ipage][idx] #ElementOcuurance
        if etype in ['KICKER']:
            if par == "HKICK":
                hc_a.append(eval)
                s_a.append(s_dict[ename, eocc])
            if par == "VKICK":
                vc_a.append(eval)


    ax = plt.figure().add_subplot(111)
    ax.plot(s_a, 1.e3*np.array(hc_a), 'bd-')
    ax.plot(s_a, 1.e3*np.array(vc_a), 'gd--')
    ax.set_xlabel('$s \, (m)$')
    ax.set_ylabel('$Corrector kick$ ($m rad$)')
    plt.show()
    return 1.e3*np.array(hc_a), 1.e3*np.array(vc_a)

def showQuads(ipage=0,fnameCalc='twiss', fname='twiss'):
    '''
    fnameCalc contains element s-positions
    fname contains parameter table of quadrupole settings
    '''
    twissFile = sdds.SDDS(0)
    twissFile.load(fnameCalc + '.twi')

    s_dict = getSdict(twissFile)

    x2 = sdds.SDDS(0)
    x2.load(fname+'.param')

    #['ElementName', 'ElementParameter', 'ParameterValue', 'ElementType', 'ElementOccurence', 'ElementGroup']
    #print(x2.columnName[0], x2.columnName[1], x2.columnName[2], x2.columnName[3], x2.columnName[4])
    n_list = len(x2.columnData[0][ipage])

    s_a = []
    k1_a = []

    for idx in range(n_list):
        ename = x2.columnData[0][ipage][idx] # ElementName
        par = x2.columnData[1][ipage][idx] #ElementParameter
        eval = x2.columnData[2][ipage][idx] #ParameterValue
        etype = x2.columnData[4][ipage][idx] #ElementType
        eocc = x2.columnData[5][ipage][idx] #ElementOcuurance
        if etype in ['KQUAD']:
            if par == "K1":
                k1_a.append(eval)
                s_a.append(s_dict[ename, eocc])


    ax = plt.figure().add_subplot(111)
    ax.plot(s_a, np.array(k1_a), 'rd')
    ax.set_xlabel('$s \, (m)$')
    ax.set_ylabel('k1')
    plt.show()

def showQuadsDiff(ipage=0,fnameCalc='twiss', fnameRef='twiss', fname='twiss'):
    '''
    fnameCalc contains element s-positions
    fname contains parameter table of quadrupole settings
    '''
    twissFile = sdds.SDDS(0)
    twissFile.load(fnameCalc + '.twi')

    s_dict = getSdict(twissFile)

    f1 = sdds.SDDS(0)
    f1.load(fnameRef+'.param')

    f2 = sdds.SDDS(0)
    f2.load(fname+'.param')

    #['ElementName', 'ElementParameter', 'ParameterValue', 'ElementType', 'ElementOccurence', 'ElementGroup']
    #print(x2.columnName[0], x2.columnName[1], x2.columnName[2], x2.columnName[3], x2.columnName[4])
    n_list = len(f1.columnData[0][ipage])

    s_a = []
    k1_a = []
    k1_a_2 = []

    for idx in range(n_list):
        ename = f1.columnData[0][ipage][idx] # ElementName
        par = f1.columnData[1][ipage][idx] #ElementParameter
        eval1 = f1.columnData[2][ipage][idx] #ParameterValue
        eval2 = f2.columnData[2][ipage][idx] #ParameterValue
        etype = f1.columnData[4][ipage][idx] #ElementType
        eocc = f1.columnData[5][ipage][idx] #ElementOcuurance
        if etype in ['KQUAD']:
            if par == "K1":
                k1_a.append(eval1)
                k1_a_2.append(eval2)
                s_a.append(s_dict[ename, eocc])


    v = np.array(k1_a_2) / np.array(k1_a)

    ax = plt.figure().add_subplot(111)
    ax.plot(s_a, v, 'rd')
    ax.set_xlabel('$s \, (m)$')
    ax.set_ylabel('k1')
    plt.show()

def getQuadFamilies(fname = 'twiss'):
    '''
    return:
        dictionary of quadrupole names->number of occurances
        dictionary of (quadrupole name, occurance) -> strength
    '''
    #twissFile = sdds.SDDS(0)
    #twissFile.load(fname + '.twi')
    x2 = sdds.SDDS(0)
    x2.load(fname+'.param')
    ipage = 0
    #['ElementName', 'ElementParameter', 'ParameterValue', 'ElementType', 'ElementOccurence', 'ElementGroup']
    #print(x2.columnName[0], x2.columnName[1], x2.columnName[2], x2.columnName[3], x2.columnName[4])
    n_list = len(x2.columnData[0][ipage])

    eocc_a = {}
    vals = {}

    for idx in range(n_list):
        ename = x2.columnData[0][ipage][idx] # ElementName
        par = x2.columnData[1][ipage][idx] #ElementParameter
        eval = x2.columnData[2][ipage][idx] #ParameterValue
        etype = x2.columnData[4][ipage][idx] #ElementType
        eocc = x2.columnData[5][ipage][idx] #ElementOcuurance
        if etype in ['KQUAD']:
            if par == "K1":
                eocc_a[ename] = int(eocc)
                vals[ename, eocc] = float(eval)

    return eocc_a, vals

def plotOrb(idx=0, fname='twiss'):
     x = sdds.SDDS(0)
     x.load(fname + '.clo')

     ax = plt.figure().add_subplot(111)
     p1, = ax.plot(x.columnData[0][idx], 1.e6*np.array(x.columnData[1][idx]), 'b-',lw=3)
     p2, = ax.plot(x.columnData[0][idx], 1.e6*np.array(x.columnData[3][idx]), 'r-',lw=3)
     ax.set_xlabel('$s \, (m)$')
     ax.set_ylabel('Orbit ($\mu m$)')
     ax.legend([p1,p2],['X','Y'])

     print('rms x [mu m]: ', 1.e6*np.std(np.array(x.columnData[1][idx])))
     print('rms y [mu m]: ', 1.e6*np.std(np.array(x.columnData[3][idx])))
     print('max x [mu m]: ', 1.e6*np.max(np.array(x.columnData[1][idx])))
     print('max y [mu m]: ', 1.e6*np.max(np.array(x.columnData[3][idx])))


     plt.show()

def getBpms(sddsFile, ipage):
    n_list = len(sddsFile.columnData[0][ipage]) # length
    #print(sddsFile.columnName)
    #print('orbit data length length:', n_list)
    xs = []
    ys = []
    for idx in range(n_list):
        s = sddsFile.columnData[0][ipage][idx]
        x = sddsFile.columnData[1][ipage][idx]
        y = sddsFile.columnData[3][ipage][idx]
        etype = sddsFile.columnData[7][ipage][idx]
        #print(etype)
        if etype in ['MONI']:
            #print('BPM!!!')
            xs.append(x)
            ys.append(y)

    return np.array(xs), np.array(ys)

def getBetaBeat(fname, fnameRef, idx=0):
    bx, by = 0.0, 0.0
    f1 = sdds.SDDS(0)
    f1.load(fname + '.twi')

    f2 = sdds.SDDS(0)
    f2.load(fnameRef + '.twi')

    #print(x.columnName[0], x.columnName[1], x.columnName[6], x.columnName[7])
    #print(x.columnData[0])
    s = f1.columnData[0][idx]
    betaX_1 = np.array(f1.columnData[1][idx])
    betaY_1 = np.array(f1.columnData[7][idx])

    betaX_2 = np.array(f2.columnData[1][idx])
    betaY_2 = np.array(f2.columnData[7][idx])

    #bx = np.sqrt( np.sum([ (betaX_1[i] - betaX_2[i])**2/betaX_2[i]**2 for i in range(len(s))]) ) / len(s)
    #by = np.sqrt( np.sum([ (betaY_1[i] - betaY_2[i])**2/betaY_2[i]**2 for i in range(len(s))]) ) / len(s)

    #bx = np.mean( (betaX_1 - betaX_2)/betaX_2 )
    #by = np.mean( (betaY_1 - betaY_2)/betaY_2 )

    bx = np.std( (betaX_1 - betaX_2)/betaX_2 )
    by = np.std( (betaY_1 - betaY_2)/betaY_2 )


    return bx, by

def getParameter(f,pname):
    id1 = f.parameterName.index(pname)
    return f.parameterData[id1]

def plotDA(fname='da.aper'):
     f = sdds.SDDS(0)
     f.load(fname)
     print(f.columnName[0], f.columnName[1])
     plt.plot(f.columnData[0][0], f.columnData[1][0], 'bd-')
     plt.show()

def evaluateDa(fname, nlines):
    f =  sdds.SDDS(0)
    f.load(fname)
    x = f.columnData[0][0]
    y = f.columnData[1][0]
    i1 = 0
    i2 = int(nlines/2)
    i3 = nlines-1
    area = 0
    for i in range(nlines-1):
        area += 0.5*(y[i]+ y[i+1])*(x[i+1]-x[i])

    return( x[i1]*1.e3, y[i2]*1.e3, x[i3]*1.e3, area * 1.e6)

# TODO
def evaluateMa(fname):
    pass


def evaluateOptics(fname):
    f = sdds.SDDS(0)
    f.load(fname)

    res = {}

    res['xi_x'] = getParameter(f,'dnux/dp')[0]
    res['xi_y'] = getParameter(f,'dnuy/dp')[0]
    res['nu_x'] = getParameter(f,'nux')[0]
    res['nu_y'] = getParameter(f,'nuy')[0]

    return res

def get_column(name):
    pass

def get_parameter(name):
    pass

def _plotDa(ids=[0], fname='da4.aper', mark='d--', ax = None, plot_average=False, scale_x=1.0, scale_y=1.0):
     x = sdds.SDDS(0)
     x.load(fname)
     #x.load(input)
     #print(x.__dict__)
     #print(x.columnName)
     print(x.columnName[0], x.columnName[1])

     if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
     ax.set_title("Dynamic Aperture (line search)",fontsize=14)
     ax.set_xlabel("X [mm]",fontsize=12)
     ax.set_ylabel("Y [mm]",fontsize=12)


     #print(x.columnData[0][0], x.columnData[1][0])
     for ipage in ids:
         ax.plot(scale_x*1.e3*np.array(x.columnData[0][ipage]), scale_y*1.e3*np.array(x.columnData[1][ipage]), mark)

     if plot_average:
         xs = 0.0*np.array(x.columnData[0][0])
         ys = 0.0*np.array(x.columnData[1][0])
         for ipage in ids:
             xs += scale_x*1.e3*np.array(x.columnData[0][ipage]) / len(ids)
             ys += scale_y*1.e3*np.array(x.columnData[1][ipage]) / len(ids)
         ax.plot(xs, ys, '-',lw=5)
     return ax


def _plotFma(ipage=0, fname='fma.fma', code='nux',de=-0.02):
     x = sdds.SDDS(0)
     x.load(fname)
     #print(x.__dict__)
     print(x.columnName)
     #idx_x = x.columnName.index('x')
     #idx_y = x.columnName.index('y')
     idx_dif = x.columnName.index(code)
     idx_nux = x.columnName.index('nux')
     idx_nuy = x.columnName.index('nuy')
     idx_x = x.columnName.index('x')
     idx_y = x.columnName.index('y')
     idx_delta = x.columnName.index('delta')

     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_title("Frequency vs. Amplitude",fontsize=14)
     #ax.set_xlabel("X [mm]",fontsize=12)
     #ax.set_ylabel("Y [mm]",fontsize=12)
     ax.set_xlabel(r"$\nu_x$",fontsize=12)
     ax.set_ylabel(r"$\nu_y$",fontsize=12)

     #ax.grid(True,linestyle='-',color='0.75')

     nuxa = np.array(x.columnData[idx_nux][ipage])
     nuya = np.array(x.columnData[idx_nuy][ipage])
     difa = np.array(x.columnData[idx_dif][ipage])
     xa = np.array(x.columnData[idx_x][ipage]) * 1.e3
     ya = np.array(x.columnData[idx_y][ipage]) * 1.e3
     delta = np.array(x.columnData[idx_delta][ipage])

     xa_cut = []
     ya_cut = []
     nux_cut = []
     nuy_cut = []
     dif_cut = []

     for i in range(len(xa)):
         if abs(delta[i] - de) < 1.e-9:
            xa_cut.append(xa[i])
            ya_cut.append(ya[i])
            nux_cut.append(nuxa[i])
            nuy_cut.append(nuya[i])
            dif_cut.append(difa[i])
     xa = np.array(xa_cut)
     ya = np.array(ya_cut)
     nuxa = np.array(nux_cut)
     nuya = np.array(nuy_cut)

     for i in range(len(nuxa)):
           if nuxa[i] > 0.5: nuxa[i] = 1.0 - nuxa[i]
     for i in range(len(nuya)):
           if nuya[i] > 0.5: nuya[i] = 1.0 - nuya[i]

     difa = np.array(dif_cut)

     za = sqrt((xa**2/3. + ya**2/2.))

     p = ax.scatter(nuxa,nuya,s=4,c=za, marker = 's', cmap = cm.jet )
     cbar= fig.colorbar(p, ax=ax)
     ax.set_title(r'$\delta$E/E='+str(de) )
     cbar.set_label(r'$\sqrt{A_X^2 + A_Y^2}$ [mm mrad]')

     #fig = plt.figure()
     #ax = fig.add_subplot(111)
     #ax.plot(delta,'d')

     fig = plt.figure()
     ax = fig.add_subplot(111)

     p = ax.scatter(xa,ya,s=4,c=difa, marker = 's', cmap = cm.jet )
     ax.set_title(r'$\delta$E/E='+str(de) )
     cbar= fig.colorbar(p, ax=ax)
     cbar.set_label('Diffusion')


def _plotAvsE(ipage=0, fname='fma.fma'):
     x = sdds.SDDS(0)
     x.load('c:\workspace\p4-work\hmba\lattice\\v15\elegant\\'+fname)
     #print(x.__dict__)
     print(x.columnName)
     #idx_x = x.columnName.index('x')
     #idx_y = x.columnName.index('y')
     idx_dif = x.columnName.index('diffusion')
     idx_nux = x.columnName.index('nux')
     idx_nuy = x.columnName.index('nuy')
     idx_x = x.columnName.index('x')
     idx_y = x.columnName.index('y')
     idx_delta = x.columnName.index('delta')

     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_xlabel(r'$\delta$E/E (%)',fontsize=12)
     #ax.set_ylabel(r"$A_y$",fontsize=12)

     #ax.grid(True,linestyle='-',color='0.75')

     nuxa = np.array(x.columnData[idx_nux][ipage])
     nuya = np.array(x.columnData[idx_nuy][ipage])
     difa = np.array(x.columnData[idx_dif][ipage])
     xa = np.array(x.columnData[idx_x][ipage]) * 1.e3
     ya = np.array(x.columnData[idx_y][ipage]) * 1.e3
     delta = np.array(x.columnData[idx_delta][ipage])*1.e2

     za = sqrt((xa**2/29. + ya**2/5.))

     p = ax.scatter(delta,za,s=4,c=za, marker = 's', cmap = cm.jet )
     cbar= fig.colorbar(p, ax=ax)
     ax.set_title('Maximum stable amplitude as function of energy offset')
     cbar.set_label(r'$\sqrt{A_X^2 + A_Y^2}$ [mm mrad]')

def _plotDetuning(ipage=0, fname='fma.fma', code='nux'):
     x = sdds.SDDS(0)
     x.load(fname)
     #print(x.__dict__)
     print(x.columnName)
     #idx_x = x.columnName.index('x')
     #idx_y = x.columnName.index('y')
     idx_dif = x.columnName.index(code)
     idx_nux = x.columnName.index('nux')
     idx_nuy = x.columnName.index('nuy')
     idx_x = x.columnName.index('x')
     idx_y = x.columnName.index('y')
     idx_delta = x.columnName.index('delta')

     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_title("Frequency vs. Amplitude",fontsize=14)
     #ax.set_xlabel("X [mm]",fontsize=12)
     #ax.set_ylabel("Y [mm]",fontsize=12)
     ax.set_xlabel(r"$\nu_x$",fontsize=12)
     ax.set_ylabel(r"$\nu_y$",fontsize=12)

     nuxa = np.array(x.columnData[idx_nux][ipage])
     nuya = np.array(x.columnData[idx_nuy][ipage])
     difa = np.array(x.columnData[idx_dif][ipage])
     xa = np.array(x.columnData[idx_x][ipage]) * 1.e3
     ya = np.array(x.columnData[idx_y][ipage]) * 1.e3
     delta = np.array(x.columnData[idx_delta][ipage])

     nux_cut = []
     nuy_cut = []
     # dq/dy
     for i in range(len(xa)):
         if abs(delta[i] - 0.0) < 1.e-9 and abs(xa[i] - 0.0) < 1.e-3 and abs(ya[i] - 0.0) < 3.5:
            nux_cut.append(nuxa[i])
            nuy_cut.append(nuya[i])
     nuxa = np.array(nux_cut)
     nuya = np.array(nuy_cut)
     p1, = ax.plot(nuxa,nuya,'gd')

     nuxa = np.array(x.columnData[idx_nux][ipage])
     nuya = np.array(x.columnData[idx_nuy][ipage])
     difa = np.array(x.columnData[idx_dif][ipage])
     xa = np.array(x.columnData[idx_x][ipage]) * 1.e3
     ya = np.array(x.columnData[idx_y][ipage]) * 1.e3
     delta = np.array(x.columnData[idx_delta][ipage])

     nux_cut = []
     nuy_cut = []
     # dq/dx
     for i in range(len(xa)):
         if abs(delta[i] - 0.0) < 1.e-9 and abs(ya[i] - 0.0) < 1.e-1  and abs(xa[i] - 0.0) < 10.0:
            nux_cut.append(nuxa[i])
            nuy_cut.append(nuya[i])
     nuxa = np.array(nux_cut)
     nuya = np.array(nuy_cut)

     p2,=ax.plot(nuxa,nuya,'rd')

     nuxa = np.array(x.columnData[idx_nux][ipage])
     nuya = np.array(x.columnData[idx_nuy][ipage])
     difa = np.array(x.columnData[idx_dif][ipage])
     xa = np.array(x.columnData[idx_x][ipage]) * 1.e3
     ya = np.array(x.columnData[idx_y][ipage]) * 1.e3
     delta = np.array(x.columnData[idx_delta][ipage])

     nux_cut = []
     nuy_cut = []

     # dq/dx
     for i in range(len(xa)):
         if abs(xa[i] - 0.0) < 1.e-9 and abs(ya[i] - 0.0) < 1.e-2 and delta[i]>=0:
            nux_cut.append(nuxa[i])
            nuy_cut.append(nuya[i])

     nuxa = np.array(nux_cut)
     nuya = np.array(nuy_cut)

     p3,=ax.plot(nuxa,nuya,'bd')

     nuxa = np.array(x.columnData[idx_nux][ipage])
     nuya = np.array(x.columnData[idx_nuy][ipage])
     difa = np.array(x.columnData[idx_dif][ipage])
     xa = np.array(x.columnData[idx_x][ipage]) * 1.e3
     ya = np.array(x.columnData[idx_y][ipage]) * 1.e3
     delta = np.array(x.columnData[idx_delta][ipage])

     nux_cut = []
     nuy_cut = []

     # dq/dx
     for i in range(len(xa)):
         if abs(xa[i] - 0.0) < 1.e-9 and abs(ya[i] - 0.0) < 1.e-2 and delta[i]<=0:
            nux_cut.append(nuxa[i])
            nuy_cut.append(nuya[i])

     nuxa = np.array(nux_cut)
     nuya = np.array(nuy_cut)

     p4,=ax.plot(nuxa,nuya,'md')

     ax.legend([p1,p2,p3,p4],['Y (<3.5mm)','X','dE>0','dE<0'])

     plt.show()

def _plotMmap(ipage=0, fname='ma.mmap', ax=None, mark='rd--'):
     x = sdds.SDDS(0)
     x.load(fname)

     #print(x.__dict__)
     #print(x.columnName)
     idx0 = x.columnName.index('s')
     idx1 = x.columnName.index('deltaNegative')
     idx2 = x.columnName.index('deltaPositive')

     if ax is None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         ax.set_title("Local momentum acceptance",fontsize=14)
         ax.set_xlabel("S (m)",fontsize=12)
         ax.set_ylabel("Dp/p",fontsize=12)

     s = np.array(x.columnData[idx0][ipage])
     d1 = np.array(x.columnData[idx1][ipage])
     d2 = np.array(x.columnData[idx2][ipage])

     print(d1[0])

     idc = np.argsort(s)

     ax.plot(s[idc],d1[idc],mark)
     ax.plot(s[idc],d2[idc],mark)
     return ax

def _plotTunes(ipage=0, fname='fma.fma', code='diffusion'):
     x = sdds.SDDS(0)
     x.load('c:\workspace\p4-work\hmba\lattice\\v15\elegant\\'+fname)
     print(x.__dict__)
     print(x.columnName)
     idx_x = x.columnName.index('nux')
     idx_y = x.columnName.index('nuy')
     idx_dif = x.columnName.index(code)

     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_title("Tunes",fontsize=14)
     ax.set_xlabel(r"$\nu_x$",fontsize=12)
     ax.set_ylabel(r"$\nu_y$",fontsize=12)
     ax.grid(True,linestyle='-',color='0.75')

     xa = np.array(x.columnData[idx_x][ipage])
     ya = np.array(x.columnData[idx_y][ipage])
     za = np.array(x.columnData[idx_dif][ipage])

     # scatter with colormap mapping to z value
     p = ax.scatter(xa,ya,s=2,c=za, marker = 'd', cmap = cm.jet )
     cbar= fig.colorbar(p, ax=ax)
     cbar.set_label(code)
     plt.show()

def _plotEmit(idx=0):
     x = sdds.SDDS(0)
     x.load('c:\workspace\p4-work\hmba\lattice\\v15\elegant\\align4r.twi')
     idx_ex0 = x.parameterName.index('ex0')
     ex0 = 1.e+12 * np.array(x.parameterData[idx_ex0])
     #plt.plot(ex0,'d')
     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_title("Emittance distribution",fontsize=14)
     ax.set_xlabel(r"$\epsilon_{x0}$",fontsize=12)
     ax.set_ylabel(r"Probability density",fontsize=12)

     n, bins, patches = ax.hist(ex0, 10, normed=True,facecolor='g', alpha=0.75)
     #plt.plot(ex0,'d')
     plt.show()
     #print(x.columnName)
     #print(x.columnName[0], x.columnName[1], x.columnName[6], x.columnName[7])

def _plotOptics(idx=0):
     x = sdds.SDDS(0)
     x.load('c:\workspace\p4-work\hmba\lattice\\v15\elegant\\align4r.twi')
     print(x.columnName)
     print(x.columnName[0], x.columnName[1], x.columnName[6], x.columnName[7])
     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_xlabel("s [m]",fontsize=12)
     ax.set_ylabel("[m]",fontsize=12)

     p1,=ax.plot(x.columnData[0][idx], np.array(x.columnData[1][idx]), 'b-', lw=3)
     #plt.plot(x.columnData[0][idx], -10*np.sqrt(30.e-9* np.array(x.columnData[1][idx])), 'b-')

     p2,=ax.plot(x.columnData[0][idx], np.array(x.columnData[7][idx]), 'r-', lw=3)
     ax2 = ax.twinx()
     p3,=ax2.plot(x.columnData[0][idx], np.array(x.columnData[4][idx]), 'g-', lw=3)
     plt.legend([p1,p2,p3],['$\\beta_x$','$\\beta_y$','$\\eta_y$'])
     #plt.plot(x.columnData[0][idx], -10*np.sqrt(30.e-9*np.array(x.columnData[7][idx])), 'r-')

     #plt.plot(x.columnData[0][idx], 1.0*np.array(x.columnData[6][idx]), 'gd')
     #plt.plot(x.columnData[0][idx], -1.0*np.array(x.columnData[6][idx]), 'gd')
     plt.show()


     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_xlabel("s [m]",fontsize=12)
     ax.set_ylabel("[m]",fontsize=12)

     p1,=ax.plot(x.columnData[0][idx], np.array(x.columnData[2][idx]), 'b-', lw=3)
     #plt.plot(x.columnData[0][idx], -10*np.sqrt(30.e-9* np.array(x.columnData[1][idx])), 'b-')

     p2,=ax.plot(x.columnData[0][idx], np.array(x.columnData[8][idx]), 'r-', lw=3)
     ax2 = ax.twinx()
     p3,=ax2.plot(x.columnData[0][idx], np.array(x.columnData[4][idx]), 'g-', lw=3)
     plt.legend([p1,p2,p3],['$\\alpha_x$','$\\alpha_y$','$\\eta_y$'])
     #plt.plot(x.columnData[0][idx], -10*np.sqrt(30.e-9*np.array(x.columnData[7][idx])), 'r-')

     #plt.plot(x.columnData[0][idx], 1.0*np.array(x.columnData[6][idx]), 'gd')
     #plt.plot(x.columnData[0][idx], -1.0*np.array(x.columnData[6][idx]), 'gd')
     plt.show()



     name = "CXY"
     s = []
     mux = []
     muy = []
     for i in range(len(x.columnData[14][idx])):
         if x.columnData[14][idx][i] == name:
             mux.append(x.columnData[3][idx][i])
             muy.append(x.columnData[9][idx][i])
             s.append(x.columnData[0][idx][i])


     p1,=plt.plot(s, np.array(mux)/(2*pi), 'bd--', lw=3)
     p2,=plt.plot(s, np.array(muy)/(2*pi), 'rd--', lw=3)
     plt.legend([p1,p2],['$\\psi_x$','$\\psi_y$'])
     plt.show()


def plotTraj(idx=0, fileName=None, aper=False, aperFile=None, pre_cor=False):

     x = sdds.SDDS(0)
     x.load(fileName)
     if aper:
         x2 = sdds.SDDS(0)
         x2.load(aperFile)

     ax = plt.figure().add_subplot(111)
     if pre_cor:
         ax.plot(x.columnData[0][2*idx], 1.e6*np.array(x.columnData[1][2*idx]), 'b--')
         ax.plot(x.columnData[0][2*idx], 1.e6*np.array(x.columnData[2][2*idx]), 'r--')

     p1, = ax.plot(x.columnData[0][2*idx+1], 1.e6*np.array(x.columnData[1][2*idx+1]), 'bo-', lw=3)
     p2, = ax.plot(x.columnData[0][2*idx+1], 1.e6*np.array(x.columnData[2][2*idx+1]), 'ro-', lw=3)

     ax.set_xlabel('$s \, (m)$')
     ax.set_ylabel('$Orbit$ ($\mu m$)')
     ax.legend([p1,p2],['X','Y'])


     print('rms x [mu m]: ', 1.e6*np.std(np.array(x.columnData[1][2*idx+1])))
     print('rms y [mu m]: ', 1.e6*np.std(np.array(x.columnData[2][2*idx+1])))
     print('max x [mu m]: ', 1.e6*np.max(np.abs(np.array(x.columnData[1][2*idx+1]))))
     print('max y [mu m]: ', 1.e6*np.max(np.abs(np.array(x.columnData[2][2*idx+1]))))

     if aper:
         plt.plot(x2.columnData[0][0][5:], 1.0*np.array(x2.columnData[6][0])[5:], 'gd')
         plt.plot(x2.columnData[0][0][5:], -1.0*np.array(x2.columnData[6][0])[5:], 'gd')
         #plt.ylim(-0.04, 0.04)
     plt.show()

def _plotOrb(idx=0):
     x = sdds.SDDS(0)
     x.load('c:\workspace\p4-work\hmba\lattice\\v15\elegant\\align4r.clo')
     #print(x.columnName)
     x2 = sdds.SDDS(0)
     x2.load('c:\workspace\p4-work\hmba\lattice\\v15\elegant\\align1.twi')

     ax = plt.figure().add_subplot(111)
     p1, = ax.plot(x.columnData[0][idx], 1.e6*np.array(x.columnData[1][idx]), 'b-',lw=3)
     p2, = ax.plot(x.columnData[0][idx], 1.e6*np.array(x.columnData[3][idx]), 'g-',lw=3)
     ax.set_xlabel('$s \, (m)$')
     ax.set_ylabel('Orbit ($\mu m$)')
     ax.legend([p1,p2],['X','Y'])

     print('rms x [mu m]: ', 1.e6*np.std(np.array(x.columnData[1][idx])))
     print('rms y [mu m]: ', 1.e6*np.std(np.array(x.columnData[3][idx])))
     print('max x [mu m]: ', 1.e6*np.max(np.abs(np.array(x.columnData[1][idx]))))
     print('max y [mu m]: ', 1.e6*np.max(np.abs(np.array(x.columnData[3][idx]))))


     #plt.plot(x2.columnData[0][0][5:], 1.e6*np.array(x2.columnData[6][0])[5:], 'gd')
     #plt.plot(x2.columnData[0][0][5:], -1.e6*np.array(x2.columnData[6][0])[5:], 'gd')
     #plt.ylim(-200, 200)
     plt.show()



#todo add chromaticity as constraint
# to be moved to optimization scriptss
def evaluateLattice(values, verbose=False, nlines=5):
    command = elegant+' da.ele  -macro=lattice=p4'
    for k in values.keys():
        command = command + ','+ k + '=' + str(values[k])
    if verbose: print(command)

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if verbose:
        for line in process.stdout:
            print(line)
    process.wait()
    i1,i2,i3, ia =evaluateDa('da.aper',nlines=nlines)
    res = {}
    res['da_x_left'] = i1
    res['da_x_right'] = i3
    res['da_y'] = i2
    res['area'] = ia
    res2 = evaluateOptics('da.twi')


    return {**res, **res2} # merge of two dictionaries

# obsolete, belongs into the optimization code
def fillDict(ar):
    values = {'of1b':ar[0],'of1d':ar[1],'sd1a':ar[2],'sd1d':ar[3],'sf2ah':ar[4],'sf2eh':ar[5]}
    return values

import re

def remove_element(buf, el):
    return re.sub(r',\s*'+el +  '\s*,', ',', buf)

def install_element_after(buf, el_after, el):
    return re.sub(r',\s*'+el_after +  '\s*,', ','+ el_after + ',' + el + ',', buf)

def install_element_before(buf, el_before, el):
    return re.sub(r',\s*'+el_before +  '\s*,', ','+ el + ',' + el_before + ',' , buf)

def make_unique(buf, el):
    ms = [m.start() for m in re.finditer(el, buf)]
    print('replacing occurences: ', len(ms))

    n = len(ms)
    for i in range(n):
        #print('replacing', i)
        buf = re.sub(r',\s*'+el +  '\s*,', ','+ el + str(i)+',' , buf, 1)
    return buf, n
