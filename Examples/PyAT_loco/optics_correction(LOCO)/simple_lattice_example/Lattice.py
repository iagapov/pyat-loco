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
QF1 = elements.Quadrupole('QF1', 0.1, 1.1 )
QF2 = elements.Quadrupole('QF2', 0.1, 1.1 )
QF3 = elements.Quadrupole('QF3', 0.1, 1.1 )
QF4 = elements.Quadrupole('QF4', 0.1, 1.1 )
QF5 = elements.Quadrupole('QF5', 0.1, 1.1 )
QF6 = elements.Quadrupole('QF6', 0.1, 1.1 )
QF7 = elements.Quadrupole('QF7', 0.1, 1.1 )
QF8 = elements.Quadrupole('QF8', 0.1, 1.1 )
QF9 = elements.Quadrupole('QF9', 0.1, 1.1 )
QF10 = elements.Quadrupole('QF10', 0.1, 1.1 )
QF11 = elements.Quadrupole('QF11', 0.1, 1.1 )
QF12 = elements.Quadrupole('QF12', 0.1, 1.1 )
QF13 = elements.Quadrupole('QF13', 0.1, 1.1 )
QF14 = elements.Quadrupole('QF14', 0.1, 1.1 )
QF15 = elements.Quadrupole('QF15', 0.1, 1.1 )
QF16 = elements.Quadrupole('QF16', 0.1, 1.1 )
QF17 = elements.Quadrupole('QF17', 0.1, 1.1 )
QF18 = elements.Quadrupole('QF18', 0.1, 1.1 )
QF19 = elements.Quadrupole('QF19', 0.1, 1.1 )
QF20 = elements.Quadrupole('QF20', 0.1, 1.1 )



QD1 = elements.Quadrupole('QD1', 0.1, -1.8)
QD2 = elements.Quadrupole('QD2', 0.1, -1.8)
QD3 = elements.Quadrupole('QD3', 0.1, -1.8)
QD4 = elements.Quadrupole('QD4', 0.1, -1.8)
QD5 = elements.Quadrupole('QD5', 0.1, -1.8)
QD6 = elements.Quadrupole('QD6', 0.1, -1.8)
QD7 = elements.Quadrupole('QD7', 0.1, -1.8)
QD8 = elements.Quadrupole('QD8', 0.1, -1.8)
QD9 = elements.Quadrupole('QD9', 0.1, -1.8)
QD10 = elements.Quadrupole('QD10', 0.1, -1.8)
QD11 = elements.Quadrupole('QD11', 0.1, -1.8)
QD12 = elements.Quadrupole('QD12', 0.1, -1.8)
QD13 = elements.Quadrupole('QD13', 0.1, -1.8)
QD14 = elements.Quadrupole('QD14', 0.1, -1.8)
QD15 = elements.Quadrupole('QD15', 0.1, -1.8)
QD16 = elements.Quadrupole('QD16', 0.1, -1.8)
QD17 = elements.Quadrupole('QD17', 0.1, -1.8)
QD18 = elements.Quadrupole('QD18', 0.1, -1.8)
QD19 = elements.Quadrupole('QD19', 0.1, -1.8)
QD20 = elements.Quadrupole('QD20', 0.1, -1.8)

QS1 = elements.Quadrupole('QS1', 0.1, 0.0)
at.tilt_elem(QS1,0.7853981633974483)

QS2 = elements.Quadrupole('QS2', 0.1, 0.0)
at.tilt_elem(QS2,0.7853981633974483)

QS3 = elements.Quadrupole('QS3', 0.1, 0.0)
at.tilt_elem(QS3,0.7853981633974483)

QS4 = elements.Quadrupole('QS4', 0.1, 0.0)
at.tilt_elem(QS4,0.7853981633974483)


QS5 = elements.Quadrupole('QS5', 0.1, 0.0)
at.tilt_elem(QS5,0.7853981633974483)

QS6 = elements.Quadrupole('QS6', 0.1, 0.0)
at.tilt_elem(QS6,0.7853981633974483)


QS7 = elements.Quadrupole('QS7', 0.1, 0.0)
at.tilt_elem(QS7,0.7853981633974483)

QS8 = elements.Quadrupole('QS8', 0.1, 0.0)
at.tilt_elem(QS8,0.7853981633974483)

QS9 = elements.Quadrupole('QS9', 0.1, 0.0)
at.tilt_elem(QS9,0.7853981633974483)

QS10 = elements.Quadrupole('QS10', 0.1, 0.0)
at.tilt_elem(QS10,0.7853981633974483)

QS11 = elements.Quadrupole('QS11', 0.1, 0.0)
at.tilt_elem(QS11,0.7853981633974483)

QS12 = elements.Quadrupole('QS12', 0.1, 0.0)
at.tilt_elem(QS12,0.7853981633974483)

QS13 = elements.Quadrupole('QS13', 0.1, 0.0)
at.tilt_elem(QS13,0.7853981633974483)

QS14 = elements.Quadrupole('QS14', 0.1, 0.0)
at.tilt_elem(QS14,0.7853981633974483)

QS15 = elements.Quadrupole('QS15', 0.1, 0.0)
at.tilt_elem(QS15,0.7853981633974483)

QS16 = elements.Quadrupole('QS16', 0.1, 0.0)
at.tilt_elem(QS16,0.7853981633974483)

QS17 = elements.Quadrupole('QS17', 0.1, 0.0)
at.tilt_elem(QS17,0.7853981633974483)

QS18 = elements.Quadrupole('QS18', 0.1, 0.0)
at.tilt_elem(QS18,0.7853981633974483)

QS19 = elements.Quadrupole('QS19', 0.1, 0.0)
at.tilt_elem(QS19,0.7853981633974483)

QS20 = elements.Quadrupole('QS20', 0.1, 0.0)
at.tilt_elem(QS20,0.7853981633974483)



SF = elements.Sextupole('SF', 0.01, 50.0)
BPM = elements.Monitor('BPM', K = 0.0, Systematic_error = [0.0, 0.0])
DD = elements.Drift('DD', 0.1)
m0 = elements.Marker('M0')

CXY = elements.Corrector('CXY', 0, [0, 0])
CXY01 = elements.Corrector('CXY01', 0, [0, 0])
CXY02 = elements.Corrector('CXY02', 0, [0, 0])
CXY03 = elements.Corrector('CXY03', 0, [0, 0])
CXY04 = elements.Corrector('CXY04', 0, [0, 0])
CXY05 = elements.Corrector('CXY05', 0, [0, 0])
CXY06 = elements.Corrector('CXY06', 0, [0, 0])
CXY07 = elements.Corrector('CXY07', 0, [0, 0])
CXY08 = elements.Corrector('CXY08', 0, [0, 0])
CXY09 = elements.Corrector('CXY09', 0, [0, 0])
CXY10 = elements.Corrector('CXY10', 0, [0, 0])
CXY11 = elements.Corrector('CXY11', 0, [0, 0])
CXY12 = elements.Corrector('CXY12', 0, [0, 0])
CXY13 = elements.Corrector('CXY13', 0, [0, 0])
CXY14 = elements.Corrector('CXY14', 0, [0, 0])
CXY15 = elements.Corrector('CXY15', 0, [0, 0])
CXY16 = elements.Corrector('CXY16', 0, [0, 0])
CXY17 = elements.Corrector('CXY17', 0, [0, 0])
CXY18 = elements.Corrector('CXY18', 0, [0, 0])
CXY19 = elements.Corrector('CXY19', 0, [0, 0])
CXY20 = elements.Corrector('CXY20', 0, [0, 0])
CXY21 = elements.Corrector('CXY21', 0, [0, 0])
CXY22 = elements.Corrector('CXY22', 0, [0, 0])
CXY23 = elements.Corrector('CXY23', 0, [0, 0])
CXY24 = elements.Corrector('CXY24', 0, [0, 0])
CXY25 = elements.Corrector('CXY25', 0, [0, 0])
CXY26 = elements.Corrector('CXY26', 0, [0, 0])
CXY27 = elements.Corrector('CXY27', 0, [0, 0])
CXY28 = elements.Corrector('CXY28', 0, [0, 0])
CXY29 = elements.Corrector('CXY29', 0, [0, 0])
CXY30 = elements.Corrector('CXY30', 0, [0, 0])
CXY31 = elements.Corrector('CXY31', 0, [0, 0])
CXY32 = elements.Corrector('CXY32', 0, [0, 0])
CXY33 = elements.Corrector('CXY33', 0, [0, 0])
CXY34 = elements.Corrector('CXY34', 0, [0, 0])
CXY35 = elements.Corrector('CXY35', 0, [0, 0])
CXY36 = elements.Corrector('CXY36', 0, [0, 0])
CXY37 = elements.Corrector('CXY37', 0, [0, 0])
CXY38 = elements.Corrector('CXY38', 0, [0, 0])
CXY39 = elements.Corrector('CXY39', 0, [0, 0])
CXY40 = elements.Corrector('CXY40', 0, [0, 0])

BPM = elements.Monitor('BPM', K = 0.0, Systematic_error = [0.0, 0.0])
BPM01 = elements.Monitor('BPM01', K = 0.0, Systematic_error = [0.0, 0.0])
BPM02 = elements.Monitor('BPM02', K = 0.0, Systematic_error = [0.0, 0.0])
BPM03 = elements.Monitor('BPM03', K = 0.0, Systematic_error = [0.0, 0.0])
BPM04 = elements.Monitor('BPM04', K = 0.0, Systematic_error = [0.0, 0.0])
BPM05 = elements.Monitor('BPM05', K = 0.0, Systematic_error = [0.0, 0.0])
BPM06 = elements.Monitor('BPM06', K = 0.0, Systematic_error = [0.0, 0.0])
BPM07 = elements.Monitor('BPM07', K = 0.0, Systematic_error = [0.0, 0.0])
BPM08 = elements.Monitor('BPM08', K = 0.0, Systematic_error = [0.0, 0.0])
BPM09 = elements.Monitor('BPM09', K = 0.0, Systematic_error = [0.0, 0.0])
BPM10 = elements.Monitor('BPM10', K = 0.0, Systematic_error = [0.0, 0.0])
BPM11 = elements.Monitor('BPM11', K = 0.0, Systematic_error = [0.0, 0.0])
BPM12 = elements.Monitor('BPM12', K = 0.0, Systematic_error = [0.0, 0.0])
BPM13 = elements.Monitor('BPM13', K = 0.0, Systematic_error = [0.0, 0.0])
BPM14 = elements.Monitor('BPM14', K = 0.0, Systematic_error = [0.0, 0.0])
BPM15 = elements.Monitor('BPM15', K = 0.0, Systematic_error = [0.0, 0.0])
BPM16 = elements.Monitor('BPM16', K = 0.0, Systematic_error = [0.0, 0.0])
BPM17 = elements.Monitor('BPM17', K = 0.0, Systematic_error = [0.0, 0.0])
BPM18 = elements.Monitor('BPM18', K = 0.0, Systematic_error = [0.0, 0.0])
BPM19 = elements.Monitor('BPM19', K = 0.0, Systematic_error = [0.0, 0.0])
BPM20 = elements.Monitor('BPM20', K = 0.0, Systematic_error = [0.0, 0.0])
BPM21 = elements.Monitor('BPM21', K = 0.0, Systematic_error = [0.0, 0.0])
BPM22 = elements.Monitor('BPM22', K = 0.0, Systematic_error = [0.0, 0.0])
BPM23 = elements.Monitor('BPM23', K = 0.0, Systematic_error = [0.0, 0.0])
BPM24 = elements.Monitor('BPM24', K = 0.0, Systematic_error = [0.0, 0.0])
BPM25 = elements.Monitor('BPM25', K = 0.0, Systematic_error = [0.0, 0.0])
BPM26 = elements.Monitor('BPM26', K = 0.0, Systematic_error = [0.0, 0.0])
BPM27 = elements.Monitor('BPM27', K = 0.0, Systematic_error = [0.0, 0.0])
BPM28 = elements.Monitor('BPM28', K = 0.0, Systematic_error = [0.0, 0.0])
BPM29 = elements.Monitor('BPM29', K = 0.0, Systematic_error = [0.0, 0.0])
BPM30 = elements.Monitor('BPM30', K = 0.0, Systematic_error = [0.0, 0.0])
BPM31 = elements.Monitor('BPM31', K = 0.0, Systematic_error = [0.0, 0.0])
BPM32 = elements.Monitor('BPM32', K = 0.0, Systematic_error = [0.0, 0.0])
BPM33 = elements.Monitor('BPM33', K = 0.0, Systematic_error = [0.0, 0.0])
BPM34 = elements.Monitor('BPM34', K = 0.0, Systematic_error = [0.0, 0.0])
BPM35 = elements.Monitor('BPM35', K = 0.0, Systematic_error = [0.0, 0.0])
BPM36 = elements.Monitor('BPM36', K = 0.0, Systematic_error = [0.0, 0.0])
BPM37 = elements.Monitor('BPM37', K = 0.0, Systematic_error = [0.0, 0.0])
BPM38 = elements.Monitor('BPM38', K = 0.0, Systematic_error = [0.0, 0.0])
BPM39 = elements.Monitor('BPM39', K = 0.0, Systematic_error = [0.0, 0.0])
BPM40 = elements.Monitor('BPM40', K = 0.0, Systematic_error = [0.0, 0.0])

ring_SF_off_QS_off = [B1H, CXY01, DD, BPM, QF1, DD, B1H, B1H, CXY02, DD, BPM01, QD1, DD, B1H,
B1H, CXY03, DD, BPM02, QF2, DD, B1H, B1H, CXY04, DD, BPM03, QD2, DD, B1H,
B1H, CXY05, DD, BPM04, QF3, DD, B1H, B1H, CXY06, DD, BPM05, QD3, DD, B1H,
B1H, CXY07, DD, BPM06, QF4, DD, B1H, B1H, CXY08, DD, BPM07, QD4, DD, B1H,
B1H, CXY09, DD, BPM08, QF5, DD, B1H, B1H, CXY10, DD, BPM09, QD5, DD, B1H,
B1H, CXY11, DD, BPM10, QF6, DD, B1H, B1H, CXY12, DD, BPM11, QD6, DD, B1H,
B1H, CXY13, DD, BPM12, QF7, DD, B1H, B1H, CXY14, DD, BPM13, QD7, DD, B1H,
B1H, CXY15, DD, BPM14, QF8, DD, B1H, B1H, CXY16, DD, BPM15, QD8, DD, B1H,
B1H, CXY17, DD, BPM16, QF9, DD, B1H, B1H, CXY18, DD, BPM17, QD9, DD, B1H,
B1H, CXY19, DD, BPM18, QF10, DD, B1H, B1H, CXY20, DD, BPM19, QD10, DD, B1H,
B1H, CXY21, DD, BPM20, QF11, DD, B1H, B1H, CXY22, DD, BPM21, QD11, DD, B1H,
B1H, CXY23, DD, BPM22, QF12, DD, B1H, B1H, CXY24, DD, BPM23, QD12, DD, B1H,
B1H, CXY25, DD, BPM24, QF13, DD, B1H, B1H, CXY26, DD, BPM25, QD13, DD, B1H,
B1H, CXY27, DD, BPM26, QF14, DD, B1H, B1H, CXY28, DD, BPM27, QD14, DD, B1H,
B1H, CXY29, DD, BPM28, QF15, DD, B1H, B1H, CXY30, DD, BPM29, QD15, DD, B1H,
B1H, CXY31, DD, BPM30, QF16, DD, B1H, B1H, CXY32, DD, BPM31, QD16, DD, B1H,
B1H, CXY33, DD, BPM32, QF17, DD, B1H, B1H, CXY34, DD, BPM33, QD17, DD, B1H,
B1H, CXY35, DD, BPM34, QF18, DD, B1H, B1H, CXY36, DD, BPM35, QD18, DD, B1H,
B1H, CXY37, DD, BPM36, QF19, DD, B1H, B1H, CXY38, DD, BPM37, QD19, DD, B1H,
B1H, CXY39, DD, BPM38, QF20, DD, B1H, B1H, CXY40, DD, BPM39, QD20, DD, B1H, m0]


ring1 = [deepcopy(x) for x in ring_SF_off_QS_off]


ring_SF_off_QS_on = [B1H, CXY01, DD, BPM, QF1, QS1,DD, B1H, B1H, CXY02, DD, BPM01, QD1, DD, B1H,
B1H, CXY03, DD, BPM02, QF2,DD, DD, B1H, B1H, CXY04, DD, BPM03, QD2, DD, B1H,
B1H, CXY05, DD, BPM04, QF3,QS2, DD, B1H, B1H, CXY06, DD, BPM05, QD3, DD, B1H,
B1H, CXY07, DD, BPM06, QF4,DD, DD, B1H, B1H, CXY08, DD, BPM07, QD4, DD, B1H,
B1H, CXY09, DD, BPM08, QF5, QS3,DD, B1H, B1H, CXY10, DD, BPM09, QD5, DD, B1H,
B1H, CXY11, DD, BPM10, QF6, DD,DD, B1H, B1H, CXY12, DD, BPM11, QD6, DD, B1H,
B1H, CXY13, DD, BPM12, QF7,QS4, DD, B1H, B1H, CXY14, DD, BPM13, QD7, DD, B1H,
B1H, CXY15, DD, BPM14, QF8, DD,DD, B1H, B1H, CXY16, DD, BPM15, QD8, DD, B1H,
B1H, CXY17, DD, BPM16, QF9, QS5,DD, B1H, B1H, CXY18, DD, BPM17, QD9, DD, B1H,
B1H, CXY19, DD, BPM18, QF10, DD,DD, B1H, B1H, CXY20, DD, BPM19, QD10, DD, B1H,
B1H, CXY21, DD, BPM20, QF11, QS6,DD, B1H, B1H, CXY22, DD, BPM21, QD11, DD, B1H,
B1H, CXY23, DD, BPM22, QF12, DD,DD, B1H, B1H, CXY24, DD, BPM23, QD12, DD, B1H,
B1H, CXY25, DD, BPM24, QF13, QS7,DD, B1H, B1H, CXY26, DD, BPM25, QD13, DD, B1H,
B1H, CXY27, DD, BPM26, QF14, DD,DD, B1H, B1H, CXY28, DD, BPM27, QD14, DD, B1H,
B1H, CXY29, DD, BPM28, QF15, QS8,DD, B1H, B1H, CXY30, DD, BPM29, QD15, DD, B1H,
B1H, CXY31, DD, BPM30, QF16, DD,DD, B1H, B1H, CXY32, DD, BPM31, QD16, DD, B1H,
B1H, CXY33, DD, BPM32, QF17, QS9,DD, B1H, B1H, CXY34, DD, BPM33, QD17, DD, B1H,
B1H, CXY35, DD, BPM34, QF18, DD,DD, B1H, B1H, CXY36, DD, BPM35, QD18, DD, B1H,
B1H, CXY37, DD, BPM36, QF19, QS10,DD, B1H, B1H, CXY38, DD, BPM37, QD19, DD, B1H,
B1H, CXY39, DD, BPM38, QF20, DD,DD, B1H, B1H, CXY40, DD, BPM39, QD20, DD, B1H, m0]


ring2 = [deepcopy(x) for x in ring_SF_off_QS_on]



ring_SF_off_QS_on_ = [B1H, CXY01, DD, BPM, QF1, QS1,DD, B1H, B1H, CXY02, DD, BPM01, QD1, DD, B1H,
B1H, CXY03, DD, BPM02, QF2,QS2, DD, B1H, B1H, CXY04, DD, BPM03, QD2, DD, B1H,
B1H, CXY05, DD, BPM04, QF3,QS3, DD, B1H, B1H, CXY06, DD, BPM05, QD3, DD, B1H,
B1H, CXY07, DD, BPM06, QF4,QS4, DD, B1H, B1H, CXY08, DD, BPM07, QD4, DD, B1H,
B1H, CXY09, DD, BPM08, QF5, QS5,DD, B1H, B1H, CXY10, DD, BPM09, QD5, DD, B1H,
B1H, CXY11, DD, BPM10, QF6, QS6,DD, B1H, B1H, CXY12, DD, BPM11, QD6, DD, B1H,
B1H, CXY13, DD, BPM12, QF7,QS7, DD, B1H, B1H, CXY14, DD, BPM13, QD7, DD, B1H,
B1H, CXY15, DD, BPM14, QF8, QS8,DD, B1H, B1H, CXY16, DD, BPM15, QD8, DD, B1H,
B1H, CXY17, DD, BPM16, QF9, QS9,DD, B1H, B1H, CXY18, DD, BPM17, QD9, DD, B1H,
B1H, CXY19, DD, BPM18, QF10, QS10,DD, B1H, B1H, CXY20, DD, BPM19, QD10, DD, B1H,
B1H, CXY21, DD, BPM20, QF11, QS11,DD, B1H, B1H, CXY22, DD, BPM21, QD11, DD, B1H,
B1H, CXY23, DD, BPM22, QF12, QS12,DD, B1H, B1H, CXY24, DD, BPM23, QD12, DD, B1H,
B1H, CXY25, DD, BPM24, QF13, QS13,DD, B1H, B1H, CXY26, DD, BPM25, QD13, DD, B1H,
B1H, CXY27, DD, BPM26, QF14, QS14,DD, B1H, B1H, CXY28, DD, BPM27, QD14, DD, B1H,
B1H, CXY29, DD, BPM28, QF15, QS15,DD, B1H, B1H, CXY30, DD, BPM29, QD15, DD, B1H,
B1H, CXY31, DD, BPM30, QF16, QS16,DD, B1H, B1H, CXY32, DD, BPM31, QD16, DD, B1H,
B1H, CXY33, DD, BPM32, QF17, QS17,DD, B1H, B1H, CXY34, DD, BPM33, QD17, DD, B1H,
B1H, CXY35, DD, BPM34, QF18, QS18,DD, B1H, B1H, CXY36, DD, BPM35, QD18, DD, B1H,
B1H, CXY37, DD, BPM36, QF19, QS19,DD, B1H, B1H, CXY38, DD, BPM37, QD19, DD, B1H,
B1H, CXY39, DD, BPM38, QF20, QS20,DD, B1H, B1H, CXY40, DD, BPM39, QD20, DD, B1H, m0]


ring3 = [deepcopy(x) for x in ring_SF_off_QS_on_]


ring_SF_off_QS_on__ = [B1H, CXY01, DD, BPM, QF1, QS1,DD, B1H, B1H, CXY02, DD, BPM01, QD1, DD, B1H,
B1H, CXY03, DD, BPM02, QF2, DD, B1H, B1H, CXY04, DD, BPM03, QD2, DD, B1H,
B1H, CXY05, DD, BPM04, QF3,QS3, DD, B1H, B1H, CXY06, DD, BPM05, QD3, DD, B1H,
B1H, CXY07, DD, BPM06, QF4, DD, B1H, B1H, CXY08, DD, BPM07, QD4, DD, B1H,
B1H, CXY09, DD, BPM08, QF5, QS5,DD, B1H, B1H, CXY10, DD, BPM09, QD5, DD, B1H,
B1H, CXY11, DD, BPM10, QF6,DD, B1H, B1H, CXY12, DD, BPM11, QD6, DD, B1H,
B1H, CXY13, DD, BPM12, QF7,QS7, DD, B1H, B1H, CXY14, DD, BPM13, QD7, DD, B1H,
B1H, CXY15, DD, BPM14, QF8,DD, B1H, B1H, CXY16, DD, BPM15, QD8, DD, B1H,
B1H, CXY17, DD, BPM16, QF9, QS9,DD, B1H, B1H, CXY18, DD, BPM17, QD9, DD, B1H,
B1H, CXY19, DD, BPM18, QF10,DD, B1H, B1H, CXY20, DD, BPM19, QD10, DD, B1H,
B1H, CXY21, DD, BPM20, QF11, QS11,DD, B1H, B1H, CXY22, DD, BPM21, QD11, DD, B1H,
B1H, CXY23, DD, BPM22, QF12,DD, B1H, B1H, CXY24, DD, BPM23, QD12, DD, B1H,
B1H, CXY25, DD, BPM24, QF13, QS13,DD, B1H, B1H, CXY26, DD, BPM25, QD13, DD, B1H,
B1H, CXY27, DD, BPM26, QF14,DD, B1H, B1H, CXY28, DD, BPM27, QD14, DD, B1H,
B1H, CXY29, DD, BPM28, QF15, QS15,DD, B1H, B1H, CXY30, DD, BPM29, QD15, DD, B1H,
B1H, CXY31, DD, BPM30, QF16,DD, B1H, B1H, CXY32, DD, BPM31, QD16, DD, B1H,
B1H, CXY33, DD, BPM32, QF17, QS17,DD, B1H, B1H, CXY34, DD, BPM33, QD17, DD, B1H,
B1H, CXY35, DD, BPM34, QF18,DD, B1H, B1H, CXY36, DD, BPM35, QD18, DD, B1H,
B1H, CXY37, DD, BPM36, QF19, QS19,DD, B1H, B1H, CXY38, DD, BPM37, QD19, DD, B1H,
B1H, CXY39, DD, BPM38, QF20,DD, B1H, B1H, CXY40, DD, BPM39, QD20, DD, B1H, m0]


ring4 = [deepcopy(x) for x in ring_SF_off_QS_on__]