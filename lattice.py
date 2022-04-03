# init_loco

#from at import atpass
from at import elements
from at import get_refpts
from at.load import load_mat
from matplotlib import pyplot as plt
import at.plot
import numpy as np
from pylab import *
from copy import copy
from copy import deepcopy


B1H = elements.Bend('B1H', 1.3, 0.314159265358979323, 0, EntranceAngle=0, ExitAngle=0)
QF = elements.Quadrupole('QF', 0.1, 1.7)
QD = elements.Quadrupole('QD', 0.1, -2.1)
#QS = elements.Quadrupole('QS', 0.1, 0.0)
#at.tilt_elem(QS,0.7853981633974483)
QS = elements.Quadrupole('QS', 0.1, 0.0, R1=array([[ 0.70710678,  0.        ,  0.70710678,  0.        ,  0.        , 0.        ], [ 0.        ,  0.70710678,  0.        ,  0.70710678,  0.        , 0.        ], [-0.70710678,  0.        ,  0.70710678,  0.        ,  0.        , 0.        ], [ 0.        , -0.70710678,  0.        ,  0.70710678,  0.        , 0.        ], [ 0.        ,  0.        ,  0.        ,  0.        ,  1.        , 0.        ], [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        , 1.        ]]), R2=array([[ 0.70710678,  0.        , -0.70710678,  0.        ,  0.        , 0.        ], [ 0.        ,  0.70710678,  0.        , -0.70710678,  0.        , 0.        ], [ 0.70710678,  0.        ,  0.70710678,  0.        ,  0.        , 0.        ], [ 0.        ,  0.70710678,  0.        ,  0.70710678,  0.        , 0.        ], [ 0.        ,  0.        ,  0.        ,  0.        ,  1.        , 0.        ], [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        , 1.        ]]))
SF = elements.Sextupole('SF', 0.01, 50.0)
BPM = elements.Monitor('BPM')
DD = elements.Drift('DD', 0.1)
m0 = elements.Marker('M0')

CXY = elements.Corrector('CXY', 0, [0, 0])
CXY01 = elements.Corrector('CXY', 0, [0, 0])
CXY02 = elements.Corrector('CXY', 0, [0, 0])
CXY03 = elements.Corrector('CXY', 0, [0, 0])
CXY04 = elements.Corrector('CXY', 0, [0, 0])
CXY05 = elements.Corrector('CXY', 0, [0, 0])
CXY06 = elements.Corrector('CXY', 0, [0, 0])
CXY07 = elements.Corrector('CXY', 0, [0, 0])
CXY08 = elements.Corrector('CXY', 0, [0, 0])
CXY09 = elements.Corrector('CXY', 0, [0, 0])
CXY10 = elements.Corrector('CXY', 0, [0, 0])
CXY11 = elements.Corrector('CXY', 0, [0, 0])
CXY12 = elements.Corrector('CXY', 0, [0, 0])
CXY13 = elements.Corrector('CXY', 0, [0, 0])
CXY14 = elements.Corrector('CXY', 0, [0, 0])
CXY15 = elements.Corrector('CXY', 0, [0, 0])
CXY16 = elements.Corrector('CXY', 0, [0, 0])
CXY17 = elements.Corrector('CXY', 0, [0, 0])
CXY18 = elements.Corrector('CXY', 0, [0, 0])
CXY19 = elements.Corrector('CXY', 0, [0, 0])
CXY20 = elements.Corrector('CXY', 0, [0, 0])
CXY21 = elements.Corrector('CXY', 0, [0, 0])
CXY22 = elements.Corrector('CXY', 0, [0, 0])
CXY23 = elements.Corrector('CXY', 0, [0, 0])
CXY24 = elements.Corrector('CXY', 0, [0, 0])
CXY25 = elements.Corrector('CXY', 0, [0, 0])
CXY26 = elements.Corrector('CXY', 0, [0, 0])
CXY27 = elements.Corrector('CXY', 0, [0, 0])
CXY28 = elements.Corrector('CXY', 0, [0, 0])
CXY29 = elements.Corrector('CXY', 0, [0, 0])
CXY30 = elements.Corrector('CXY', 0, [0, 0])
CXY31 = elements.Corrector('CXY', 0, [0, 0])
CXY32 = elements.Corrector('CXY', 0, [0, 0])
CXY33 = elements.Corrector('CXY', 0, [0, 0])
CXY34 = elements.Corrector('CXY', 0, [0, 0])
CXY35 = elements.Corrector('CXY', 0, [0, 0])
CXY36 = elements.Corrector('CXY', 0, [0, 0])
CXY37 = elements.Corrector('CXY', 0, [0, 0])
CXY38 = elements.Corrector('CXY', 0, [0, 0])
CXY39 = elements.Corrector('CXY', 0, [0, 0])
CXY40 = elements.Corrector('CXY', 0, [0, 0])

ring_SF_off_QS_off = [B1H, CXY01, DD, BPM, QF, DD, B1H, B1H, CXY02, DD, BPM, QD, DD, B1H,
B1H, CXY03, DD, BPM, QF, DD, B1H, B1H, CXY04, DD, BPM, QD, DD, B1H,
B1H, CXY05, DD, BPM, QF, DD, B1H, B1H, CXY06, DD, BPM, QD, DD, B1H,
B1H, CXY07, DD, BPM, QF, DD, B1H, B1H, CXY08, DD, BPM, QD, DD, B1H,
B1H, CXY09, DD, BPM, QF, DD, B1H, B1H, CXY10, DD, BPM, QD, DD, B1H,
B1H, CXY11, DD, BPM, QF, DD, B1H, B1H, CXY12, DD, BPM, QD, DD, B1H,
B1H, CXY13, DD, BPM, QF, DD, B1H, B1H, CXY14, DD, BPM, QD, DD, B1H,
B1H, CXY15, DD, BPM, QF, DD, B1H, B1H, CXY16, DD, BPM, QD, DD, B1H,
B1H, CXY17, DD, BPM, QF, DD, B1H, B1H, CXY18, DD, BPM, QD, DD, B1H,
B1H, CXY19, DD, BPM, QF, DD, B1H, B1H, CXY20, DD, BPM, QD, DD, B1H,
B1H, CXY21, DD, BPM, QF, DD, B1H, B1H, CXY22, DD, BPM, QD, DD, B1H,
B1H, CXY23, DD, BPM, QF, DD, B1H, B1H, CXY24, DD, BPM, QD, DD, B1H,
B1H, CXY25, DD, BPM, QF, DD, B1H, B1H, CXY26, DD, BPM, QD, DD, B1H,
B1H, CXY27, DD, BPM, QF, DD, B1H, B1H, CXY28, DD, BPM, QD, DD, B1H,
B1H, CXY29, DD, BPM, QF, DD, B1H, B1H, CXY30, DD, BPM, QD, DD, B1H,
B1H, CXY31, DD, BPM, QF, DD, B1H, B1H, CXY32, DD, BPM, QD, DD, B1H,
B1H, CXY33, DD, BPM, QF, DD, B1H, B1H, CXY34, DD, BPM, QD, DD, B1H,
B1H, CXY35, DD, BPM, QF, DD, B1H, B1H, CXY36, DD, BPM, QD, DD, B1H,
B1H, CXY37, DD, BPM, QF, DD, B1H, B1H, CXY38, DD, BPM, QD, DD, B1H,
B1H, CXY39, DD, BPM, QF, DD, B1H, B1H, CXY40, DD, BPM, QD, DD, B1H, m0]
ring1 = [deepcopy(x) for x in ring_SF_off_QS_off]




ring_SF_on= [B1H, CXY01, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY02, DD, BPM, QD, DD, B1H,
B1H, CXY03, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY04, DD, BPM, QD, DD, B1H,
B1H, CXY05, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY06, DD, BPM, QD, DD, B1H,
B1H, CXY07, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY08, DD, BPM, QD, DD, B1H,
B1H, CXY09, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY10, DD, BPM, QD, DD, B1H,
B1H, CXY11, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY12, DD, BPM, QD, DD, B1H,
B1H, CXY13, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY14, DD, BPM, QD, DD, B1H,
B1H, CXY15, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY16, DD, BPM, QD, DD, B1H,
B1H, CXY17, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY18, DD, BPM, QD, DD, B1H,
B1H, CXY19, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY20, DD, BPM, QD, DD, B1H,
B1H, CXY21, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY22, DD, BPM, QD, DD, B1H,
B1H, CXY23, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY24, DD, BPM, QD, DD, B1H,
B1H, CXY25, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY26, DD, BPM, QD, DD, B1H,
B1H, CXY27, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY28, DD, BPM, QD, DD, B1H,
B1H, CXY29, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY30, DD, BPM, QD, DD, B1H,
B1H, CXY31, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY32, DD, BPM, QD, DD, B1H,
B1H, CXY33, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY34, DD, BPM, QD, DD, B1H,
B1H, CXY35, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY36, DD, BPM, QD, DD, B1H,
B1H, CXY37, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY38, DD, BPM, QD, DD, B1H,
B1H, CXY39, DD, BPM, QF,SF,QS, DD, B1H, B1H, CXY40, DD, BPM, QD, DD, B1H, m0]

ring2 = [deepcopy(x) for x in ring_SF_on]


ring_SF_off= [B1H, CXY01, DD, BPM, QF,QS, DD, B1H, B1H, CXY02, DD, BPM, QD, DD, B1H,
B1H, CXY03, DD, BPM, QF,QS, DD, B1H, B1H, CXY04, DD, BPM, QD, DD, B1H,
B1H, CXY05, DD, BPM, QF,QS, DD, B1H, B1H, CXY06, DD, BPM, QD, DD, B1H,
B1H, CXY07, DD, BPM, QF,QS, DD, B1H, B1H, CXY08, DD, BPM, QD, DD, B1H,
B1H, CXY09, DD, BPM, QF,QS, DD, B1H, B1H, CXY10, DD, BPM, QD, DD, B1H,
B1H, CXY11, DD, BPM, QF,QS, DD, B1H, B1H, CXY12, DD, BPM, QD, DD, B1H,
B1H, CXY13, DD, BPM, QF,QS, DD, B1H, B1H, CXY14, DD, BPM, QD, DD, B1H,
B1H, CXY15, DD, BPM, QF,QS, DD, B1H, B1H, CXY16, DD, BPM, QD, DD, B1H,
B1H, CXY17, DD, BPM, QF,QS, DD, B1H, B1H, CXY18, DD, BPM, QD, DD, B1H,
B1H, CXY19, DD, BPM, QF,QS, DD, B1H, B1H, CXY20, DD, BPM, QD, DD, B1H,
B1H, CXY21, DD, BPM, QF,QS, DD, B1H, B1H, CXY22, DD, BPM, QD, DD, B1H,
B1H, CXY23, DD, BPM, QF,QS, DD, B1H, B1H, CXY24, DD, BPM, QD, DD, B1H,
B1H, CXY25, DD, BPM, QF,QS, DD, B1H, B1H, CXY26, DD, BPM, QD, DD, B1H,
B1H, CXY27, DD, BPM, QF,QS, DD, B1H, B1H, CXY28, DD, BPM, QD, DD, B1H,
B1H, CXY29, DD, BPM, QF,QS, DD, B1H, B1H, CXY30, DD, BPM, QD, DD, B1H,
B1H, CXY31, DD, BPM, QF,QS, DD, B1H, B1H, CXY32, DD, BPM, QD, DD, B1H,
B1H, CXY33, DD, BPM, QF,QS, DD, B1H, B1H, CXY34, DD, BPM, QD, DD, B1H,
B1H, CXY35, DD, BPM, QF,QS, DD, B1H, B1H, CXY36, DD, BPM, QD, DD, B1H,
B1H, CXY37, DD, BPM, QF,QS, DD, B1H, B1H, CXY38, DD, BPM, QD, DD, B1H,
B1H, CXY39, DD, BPM, QF,QS, DD, B1H, B1H, CXY40, DD, BPM, QD, DD, B1H, m0]

ring3 = [deepcopy(x) for x in ring_SF_off]



ring4_ = [B1H, CXY01, DD, BPM, QF, QS, DD, B1H, B1H, CXY02, DD, BPM, QD, DD, B1H, m0]
ring4 = [deepcopy(x) for x in ring4_]