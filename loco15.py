#! /usr/bin/env python
import sys, time, os
try:
    import sdds
    from elegant import *
except:
    print("Failed to import elegant utils; elegant engine will not work properly")
from pylab import *

try:
    import at
    from at import atpass
    from at import elements
except:
    print("Failed to import AT;  engine will not work")



class Loco:
    def __init__(self,engine):
        if engine == 'elegant':
            self.exec = 'Users/ia/Products/elegant/darwin-x86/elegant'
            self.engine = 'elegant'
            self.verbose = True
        if engine == 'at':
            self.exec = None
            self.engine = 'at'
            self.verbose = True


    def runCommand(self,command):
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if self.verbose:
            for line in process.stdout:
                print(line.decode('UTF-8'))
        process.wait()

    def prepareTwiss(self):
        if self.engine == 'elegant':
            command = self.exec + ' twiss.ele -macro=lattice=fodo'
            self.runCommand(command)
            self.Cx = getOrm(fname='data/twiss.hrm')
            self.Cy = getOrm(fname='data/twiss.vrm')

        if self.engine == 'at':
            print("preparing twiss ..")
            self.opt = at.linopt(self.lattice, refpts=np.arange(1, len(self.lattice)), get_chrom=True)
            tune = self.opt[1]
            chrom = self.opt[2]
            tws = self.opt[3]  # optics functions
            print(f"Tunes={tune}")
            print(f"Chrom={chrom}")
            print(f"calculating ORM ...")
            self.Cx = getOrm_AT_x(self)
            self.Cy = getOrm_AT_y(self)
            #print("finish ORM")
    
    def simulateMachine(self): # simulated errors, compute twiss etc.
        if self.engine == 'elegant':
            command = self.exec+' perturb.ele -macro=lattice=fodo'
            self.runCommand(command)
            
    def measureOrm(self,kickX, kickY, mode):
        if self.engine == 'elegant':
            for corName in self.corNames:
                if self.verbose: print('measuring orm:',corName)
                command = self.exec + ' measure_orm_x.ele -macro=lattice=fodo,cor_name=' + corName + ',cor_var=' + str(kickX)
                self.runCommand(command)
                command = self.exec + ' measure_orm_y.ele -macro=lattice=fodo,cor_name=' + corName + ',cor_var=' + str(kickY)
                self.runCommand(command)




# orbit respoonse matrix
class ORM:
    def __init__(self, nCor, nBpm):
        hcors = []
        vcors = []
        self.Rxx = np.zeros([nCor,nBpm])
        self.Ryy = np.zeros([nCor,nBpm])
        self.Rxy = np.zeros([nCor,nBpm])
        self.Ryx = np.zeros([nCor,nBpm])

# orbit response data
class OrbitResponse:
    def __init__(self):
        self.hcors = []
        self.orbits = {}
        self.scanValues = {}

def getBpmResponse(orbits, refOrbit, iBpm, dkick):

    x = [ orbits[0][iBpm] , refOrbit[0][iBpm]  ]# i is ipage
    y = [ orbits[1][iBpm] ,refOrbit[1][iBpm] ]

    cx = (x[0] - x[1]) / dkick # coefficient, first order
    cy = (y[0] - y[1]) / dkick # coefficient, first order

    return cx,cy


def getScanValues(cor, paramFile, nPages):
    '''
    get array of kicker angles for a single corrector during a scan
    '''
    #cor = fname.replace('orm_x','')

    print('getting scan values for', cor)

    ipage = 0
    n_list = len(paramFile.columnData[0][ipage])
    vx_a = []
    vy_a = []

    for idx in range(n_list):
        for ipage in range(nPages):
            ename = paramFile.columnData[0][ipage][idx] # ElementName
            par = paramFile.columnData[1][ipage][idx] #ElementParameter
            eval = paramFile.columnData[2][ipage][idx] #ParameterValue
            etype = paramFile.columnData[3][ipage][idx] #ElementType
            eocc = paramFile.columnData[4][ipage][idx] #ElementOcuurance
            if etype in ['KICKER'] and ename.lower() == cor.lower():
                if par == "HKICK":
                    vx_a.append(eval)
                if par == "VKICK":
                    vy_a.append(eval)

    return np.array(vx_a), np.array(vy_a)

# build ORM from elegant sdds files, simulation of measurement
# units are [m/rad] or [mm/mrad]
def buildOrmOld(fnames, nPages, plane='x'):
    orb_file = sdds.SDDS(0)
    paramFile  = sdds.SDDS(0)

    ore = OrbitResponse()
    orb_file.load(fnames[0] + '.clo')
    x_, y_ =  getBpms(orb_file, ipage=0)
    n_bpms = len(x_)

    print('building ORM', len(fnames), n_bpms)

    orm = ORM(len(fnames), n_bpms)

    i_cor = 0
    for fname in fnames:
        orb_file.load(fname + '.clo')
        paramFile.load(fname + '.param')
        print(orb_file.columnName)
        ore.hcors.append(fname)
        ore.orbits[fname] = {}
        cname = fname.replace('orm_'+plane+'_','')
        vx, vy = getScanValues(cname, paramFile, nPages)
        if plane == 'x' : ore.scanValues[fname] = vx
        if plane == 'y' : ore.scanValues[fname] = vy

        for ipage in range(nPages):
            ore.orbits[fname][ipage] = getBpms(orb_file, ipage=ipage)

        n_bpms = len(ore.orbits[fname][0][0])
        for i_bpm  in range(n_bpms):
            cx, cy = get_bpm_response(ore.orbits[fname], i_bpm,  ore.scanValues[fname])
            if plane == 'x':
                orm.Rxx[i_cor, i_bpm] = cx
                orm.Rxy[i_cor, i_bpm] = cy
            if plane == 'y':
                orm.Ryx[i_cor, i_bpm] = cx
                orm.Ryy[i_cor, i_bpm] = cy
        i_cor += 1

    return orm, ore

def buildOrm(fnames, refFileName, plane='x', dkick=1.e-4):
    '''
    fnames:        measurement files (closed orbits)
    refFileName: reference orbit file, info to be appended to the orbit response object
    returns:
         orm : contains Rxx, Rxy or Ryx and Ryy
         ore : contains orbits
    '''
    orb_file = sdds.SDDS(0)
    ref_file  = sdds.SDDS(0)

    ore = OrbitResponse()
    #print('opening first measurement file...')
    orb_file.load(fnames[0] + '.clo')
    #print('opening reference file...')
    ref_file.load(refFileName + '.clo')
    x_, y_ =  getBpms(orb_file, ipage=0)
    n_bpms = len(x_)

    print('building ORM', len(fnames), n_bpms)

    orm = ORM(len(fnames), n_bpms)

    ore.orbits['refOrbit'] = getBpms(ref_file, ipage=0)

    i_cor = 0
    for fname in fnames:
        print('reading', fname + '.clo')
        orb_file.load(fname + '.clo')
        #print(orb_file.columnName)
        ore.hcors.append(fname)
        cname = fname.replace('orm_'+plane+'_','')
        #vx, vy = getScanValues(cname, paramFile, nPages)
        #if plane == 'x' : ore.scanValues[fname] = vx
        #if plane == 'y' : ore.scanValues[fname] = vy

        ipage = 0
        ore.orbits[fname] = getBpms(orb_file, ipage=ipage)

        n_bpms = len(ore.orbits[fname][0])
        for i_bpm  in range(n_bpms):
            cx, cy = getBpmResponse(ore.orbits[fname], ore.orbits['refOrbit'], i_bpm, dkick = dkick)
            if plane == 'x':
                orm.Rxx[i_cor, i_bpm] = cx
                orm.Rxy[i_cor, i_bpm] = cy
            if plane == 'y':
                orm.Ryx[i_cor, i_bpm] = cx
                orm.Ryy[i_cor, i_bpm] = cy
        i_cor += 1

    return orm, ore

def getOrm(fname):
    if  not os.path.isfile(fname):
        print("{}: file {} does not exist, quitting".format(getOrm.__name__, fname)) # i'm sorry for that line
        return


    fi = sdds.SDDS(0)
    fi.load(fname)

    n_cor = len(fi.columnName) - 2
    ipage=0
    n_bpm = len(fi.columnData[0][ipage])
    R = np.zeros([n_cor,n_bpm])
    for i in range(n_cor):
        for j in range(n_bpm):
            R[i,j] = fi.columnData[i+2][ipage][j]
    return R


def getOrm_AT_x(self):
    #from pylab import *
    cx = []
    elements = []
    for j in range(len(self.correctors_indexes)):

        # print('+')
        self.lattice[self.correctors_indexes[j]].KickAngle = [self.dkick, 0.00]

        elements.append(self.lattice[self.correctors_indexes[j]])

        lindata0, tune, chrom, lindata = self.lattice.linopt(get_chrom=True, refpts=self.BPM_indexes)
        s_pos = lindata['s_pos']
        closed_orbit = lindata['closed_orbit']
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 1]

        self.lattice[self.correctors_indexes[j]].KickAngle = [0, 0.00]
        cx.append(closed_orbity)
    Cx = np.squeeze(cx) / self.dkick

    for j in Cx:

        idxs = np.argwhere(Cx == j)
        for i in range(len(idxs)):
            file = open("orm_x_CXY_" + str(idxs[i]) + ".txt", "w")
            for k in elements:
                str1 = str(elements[k])
                file.write(str1)
                file.close()

    return Cx




def getOrm_AT_y(self):
    #from pylab import *
    cy = []
    for j in range(len(self.correctors_indexes)):

        # print('+')
        self.lattice[self.correctors_indexes[j]].KickAngle = [self.dkick, 0.00]

        lindata0, tune, chrom, lindata = self.lattice.linopt(get_chrom=True, refpts=self.BPM_indexes)
        s_pos = lindata['s_pos']
        closed_orbit = lindata['closed_orbit']
        closed_orbitx = lindata['closed_orbit'][:, 0]
        closed_orbity = lindata['closed_orbit'][:, 1]




        self.lattice[self.correctors_indexes[j]].KickAngle = [0, 0.00]
        cy.append(closed_orbity)
    #Cy = cy
    Cy = np.squeeze(cy) / self.dkick

    return Cy







def getTheorOrm(idx=0,fname='twiss'):
     f = sdds.SDDS(0)
     f.load(fname + '.twi')
     #print(x.columnName)
     #print(x.columnName[0], x.columnName[1], x.columnName[6], x.columnName[7])
     #print(x.columnData[0])
     s = f.columnData[0][idx]
     betaX = np.array(f.columnData[1][idx])
     betaY = np.array(f.columnData[7][idx])
     muX = np.array(f.columnData[3][idx])
     etaX = np.array(f.columnData[4][idx])
     muY = np.array(f.columnData[9][idx])
     names = np.array(f.columnData[14][idx])
     etypes = np.array(f.columnData[16][idx])

     #print('tunes:', muX[-1], muY[-1])
     Qx = muX[-1] / (2.*pi)
     Qy = muY[-1] / (2.*pi)

     L = s[-1]
     #print('L:', L)

     idx = f.parameterName.index('alphac')
     alphac = float(f.parameterData[idx][0])

     #print('alphac:', alphac)


     idx_cor = []
     idx_bpm = []

     for i in range(len(names)):
         if etypes[i] == 'KICKER': idx_cor.append(i)
         if etypes[i] == 'MONI': idx_bpm.append(i)

     Cx = np.zeros([len(idx_cor), len(idx_bpm)])
     Cy = np.zeros([len(idx_cor), len(idx_bpm)])

     Cxy = np.zeros([len(idx_cor), len(idx_bpm)])
     Cyx = np.zeros([len(idx_cor), len(idx_bpm)])

     for i in range(len(idx_cor)):
         for j in range(len(idx_bpm)):
             ii = idx_cor[i]
             ij = idx_bpm[j]
             Cx[i][j] = np.sqrt( betaX[ii] * betaX[ij]) * cos(abs(muX[ii] - muX[ij]) - pi*Qx) / (2*sin(pi*Qx))
             Cy[i][j] = np.sqrt( betaY[ii] * betaY[ij]) * cos(abs(muY[ii] - muY[ij]) - pi*Qy) / (2*sin(pi*Qy))
             #Cx[i][j] += etaX[ii] * etaX[ij] / (alphac * L )
             #Cy[i][j] += etaX[ii] * etax[ij] / (alphac * L )

     return Cx, Cy
