from matplotlib import pyplot as plt
import numpy as np
from pylab import *

def showQuadsDiff(s_pos, quad_strength , perturped_quad_strength):
    showQuadsDiff = np.array(perturped_quad_strength) / np.array(quad_strength)
    plt.scatter(s_pos, showQuadsDiff)
    plt.show()

def getBetaBeat(beta_x, beta_y, perturped_beta_x, perturped_beta_y):
    beta_x = np.array(beta_x)
    beta_y = np.array(beta_y)
    perturped_beta_x = np.array(perturped_beta_x)
    perturped_beta_y = np.array(perturped_beta_y)

    bx = (beta_x -perturped_beta_x  )/perturped_beta_x
    by = ( beta_y-perturped_beta_y  )/perturped_beta_y

    return bx, by


