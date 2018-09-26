import numpy as np
import matplotlib.pyplot as plt

def red_noise(mu,sigma,n_elements):
    white_noise = np.random.normal(mu,sigma,n_elements)


    red = np.cumsum(white_noise)
    return red




