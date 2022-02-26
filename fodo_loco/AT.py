from matplotlib import pyplot as plt
import numpy as np
from pylab import *

def showQuadsDiff(s_pos, quad_strength , perturped_quad_strength):
    showQuadsDiff = np.array(perturped_quad_strength) / np.array(quad_strength)
    plt.scatter(s_pos, showQuadsDiff)
    plt.show()
