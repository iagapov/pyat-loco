from matplotlib import pyplot as plt
import numpy as np
from pylab import *

def showQuadsDiff(s_pos, quad_strength , perturped_quad_strength):
    showQuadsDiff = np.array(perturped_quad_strength) / np.array(quad_strength)
    plt.scatter(s_pos, showQuadsDiff)
    plt.show()

def getBetaBeat(beta_x, beta_y, perturped_beta_x, perturped_beta_y):
    beta_x = beta_x
    beta_y = beta_y
    perturped_beta_x = perturped_beta_x
    perturped_beta_y = perturped_beta_y

    bx = np.std((perturped_beta_x - beta_x )/beta_x)
    by = np.std((perturped_beta_y - beta_y)/beta_y)

    return bx, by

def showCorrectors(s_pos, corrector_indexes , lattice):
    correctors_kick_angle = []
    i = 0
    while (i < len(corrector_indexes)):
        corrector_kick_angle = lattice[corrector_indexes[i]].KickAngle
        correctors_kick_angle.append(corrector_kick_angle)
        i += 1
        return correctors_kick_angle



