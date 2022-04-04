# pyAT, test CELL

import at
import at.plot
import numpy as np
import matplotlib.pyplot as plt
from math import pi

QF=at.Quadrupole('QF',0.5,1.2)
#print(QF)
#SF = at.Sextupole('SF', 0.1, 1)
#SD = at.Sextupole('SD', 0.1, -1)
#drs = at.Drift('DRS', 0.2)

Dr = at.Drift('Dr', 0.5)
HalfDr = at.Drift('Dr2', 0.25)
QD = at.Quadrupole('QD', 0.5, -1.2)
Bend = at.Dipole('Bend', 1, 2*pi/40)

FODOcell = at.Lattice([HalfDr, Bend, Dr, QF, Dr, Bend, Dr, QD, HalfDr],
                      name='Simple FODO cell', energy=1E9)

FODO = FODOcell*4
#print(FODO)

refqf = at.get_cells(FODO, 'FamName', 'QF')   # FamName attribute == QF
print(f"List of QF quads {list(FODO[refqf])}")
refqd = at.get_cells(FODO, 'FamName', 'QD')   # FamName attribute == QD
print(list(FODO[refqd]))
refbends = at.get_cells(FODO, 'BendingAngle') # Existing BendingAngle attribute
print(list(FODO[refbends]))


#at.shift_elem(FODO[refqf], deltax=0.05, deltaz=0.0, relative=False)
at.shift_elem(FODO[3], deltax=0.01, deltaz=0.0, relative=True)


if False:
    refq1 = at.get_cells(FODOcell, at.checktype(at.Quadrupole))   # class == Quadrupole
    print(list(FODOcell[refq1]))
    refq2 = at.get_cells(FODOcell, at.checkname('Q[FD]'))         # name matches a pattern
    print(list(FODOcell[refq2]))

    for elem in FODOcell.select(refqf | refqd):
        print(elem)

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



[elemdata0, beamdata, elemdata] = at.get_optics(FODO, range(len(FODO)+1))
print(elemdata)
if False:
    plt.plot(elemdata.s_pos, elemdata.beta)
    plt.xlabel('s [m]')
    plt.ylabel(r'$\beta$ [m]')


#FODOcellSext.plot_beta()
FODO.plot_beta()
plt.show()

lindata0, tune, chrom, lindata = FODO.linopt(get_chrom=True, refpts=range(len(FODO)+1))
closed_orbitx = lindata['closed_orbit'][:, 0]

print("closed orbit:")
plt.plot(elemdata.s_pos, elemdata.closed_orbit[:,0])

#plt.plot()

#print(FODO[3])
plt.show()