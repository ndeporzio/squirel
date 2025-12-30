from utils import *
import numpy as np
import matplotlib.pyplot as plt
from classy import Class

# test script for whatever I'm doing
z_array = np.logspace(4,-8,100)
print(z_array,1/(1+z_array)**(-3)*asympt_scaling(1/(1+z_array),0.3,1e3))
