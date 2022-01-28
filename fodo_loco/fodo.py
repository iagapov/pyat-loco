# pyAT, test CELL

import at
import at.plot
import numpy as np
import matplotlib.pyplot as plt
from math import pi

QF=at.Quadrupole('QF',0.5,1.2)
print(QF)

Dr = at.Drift('Dr', 0.5)
HalfDr = at.Drift('Dr2', 0.25)
QD = at.Quadrupole('QD', 0.5, -1.2)
Bend = at.Dipole('Bend', 1, 2*pi/40)

FODOcell = at.Lattice([HalfDr, Bend, Dr, QF, Dr, Bend, Dr, QD, HalfDr],
                      name='Simple FODO cell', energy=1E9)
print(FODOcell)

FODO = FODOcell*20
print(FODO)

refqf = at.get_cells(FODOcell, 'FamName', 'QF')   # FamName attribute == QF
print(list(FODOcell[refqf]))
refqd = at.get_cells(FODOcell, 'FamName', 'QD')   # FamName attribute == QD
print(list(FODOcell[refqd]))
refbends = at.get_cells(FODOcell, 'BendingAngle') # Existing BendingAngle attribute
print(list(FODOcell[refbends]))

refq1 = at.get_cells(FODOcell, at.checktype(at.Quadrupole))   # class == Quadrupole
print(list(FODOcell[refq1]))
refq2 = at.get_cells(FODOcell, at.checkname('Q[FD]'))         # name matches a pattern
print(list(FODOcell[refq2]))

for elem in FODOcell.select(refqf | refqd):
    print(elem)

if False:
    nturns=200
    Z01 = np.array([.001, 0, 0, 0, 0, 0])
    Z02 = np.array([.002, 0, 0, 0, 0, 0])
    Z03 = np.array([.003, 0, 0, 0, 0, 0])
    Z1=at.lattice_pass(FODO,Z01,nturns)
    Z2=at.lattice_pass(FODO,Z02,nturns)
    Z3=at.lattice_pass(FODO,Z03,nturns)
    plt.plot(Z1[0, 0, 0, :], Z1[1, 0, 0, :],'.')
    plt.plot(Z2[0, 0, 0, :], Z2[1, 0, 0, :],'.')
    plt.plot(Z3[0, 0, 0, :], Z3[1, 0, 0, :],'.')
    plt.show()

[_, beamdata, _] = at.get_optics(FODO, get_chrom=True)

print(beamdata.tune)
print(beamdata.chromaticity)

m44, _ = at.find_m44(FODO,0)
print(m44)

SF = at.Sextupole('SF', 0.1, 1)
SD = at.Sextupole('SD', 0.1, -1)
drs = at.Drift('DRS', 0.2)

FODOcellSext = FODOcell.copy()
FODOcellSext[6:7] = [drs,SD,drs]
FODOcellSext[2:3] = [drs,SF,drs]
FODOSext = FODOcellSext*20
print(FODOSext)

[_, beamdata, _] = at.get_optics(FODOSext, get_chrom=True)
print(beamdata.tune)
print(beamdata.chromaticity)

[elemdata0, beamdata, elemdata] = at.get_optics(FODOcellSext, range(len(FODOcellSext)+1))

if False:
    plt.plot(elemdata.s_pos, elemdata.beta)
    plt.xlabel('s [m]')
    plt.ylabel(r'$\beta$ [m]')


#FODOcellSext.plot_beta()
FODO.plot_beta()
plt.show()