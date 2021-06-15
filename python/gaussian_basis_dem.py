import numpy as np
import matplotlib.pyplot as plt

def dem_from_triangles(n,x):
    temp = np.arange(n, dtype='float64')
    dem = np.zeros_like(temp)
    x = np.array(x)
    for i in range(len(x)):
        step = n / (x.shape[0] - 1)
        shift = i * step
        triangle = (1-np.abs(temp-shift)/step)*x[i]
        triangle[triangle<0] = 0
        dem += triangle


    return dem

def dem_from_gaussians(n,x):
    temp = np.arange(n, dtype='float64')
    dem = np.zeros_like(temp)
    x = np.array(x)
    for i in range(len(x)):
        step = n / (x.shape[0] - 1)
        shift = i * step
        print(shift)
        gaussian = x[i] * np.exp(-np.square(temp - shift) / (step*1))
        dem += gaussian
    return dem


if __name__ == '__main__':
    x = [27,22,25,22,22,20,20,20,20,20,20]
    n = 26
    temp = np.arange(n)

    dem = dem_from_triangles(n,x)

    print(dem.shape)
    fig,ax = plt.subplots()
    ax.plot(temp,dem)
    plt.show()

