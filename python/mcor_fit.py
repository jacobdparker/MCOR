import numpy as np
import scipy.io
import pathlib
import matplotlib.pyplot as plt
from kgpy.plot import CubeSlicer
from kgpy.img.coalignment import image_coalignment as img_align
import scipy.optimize
import scipy.signal
import time
import scipy.fft as fft
import scipy.interpolate
import pandas
import ChiantiPy.core as ch


def dem_from_triangles(n, x):
    temp = np.arange(n, dtype='float64')
    dem = np.zeros_like(temp)
    x = np.array(x)
    for i in range(len(x)):
        step = n / (x.shape[0] - 1)
        shift = i * step
        triangle = (1 - np.abs(temp - shift) / step) * x[i]
        triangle[triangle < 0] = 0
        dem += triangle
    return dem


def dem_from_gaussians(n, x):
    temp = np.arange(n, dtype='float64')
    dem = np.zeros_like(temp)
    x = np.array(x)
    for i in range(len(x)):
        step = n / (x.shape[0] - 1)
        shift = i * step
        # print(shift)
        gaussian = x[i] * np.exp(-np.square(temp - shift) / (step * .45))
        dem += gaussian
    return dem


def triangle_dem_merit(x, cube, mcor, n):
    x = dem_from_triangles(n, x)
    cube = cube * 10 ** x[..., None, None, None]
    cube = np.sum(cube, axis=0)
    tcor = mean_row_cc(cube)

    fit = np.sqrt(np.sum((tcor - mcor) ** 2))
    print(fit)

    return fit


def gaussian_dem_merit(x, cube, mcor):
    x = dem_from_gaussians(26, x)
    cube = cube * 10 ** x[..., None, None, None]
    cube = np.sum(cube, axis=0)
    tcor = mean_row_cc(cube)

    fit = np.sqrt(np.sum((tcor - mcor) ** 2))
    print(fit)

    return fit


def merit(x, cube, mcor):
    cube = cube * 10 ** x[..., None, None, None]
    cube = np.sum(cube, axis=0)
    tcor = mean_row_cc(cube)

    fit = np.sqrt(np.sum((tcor - mcor) ** 2))
    print(fit)

    return fit


def mean_row_cc(cube):
    cube /= np.median(cube, axis=(-1, -2), keepdims=True)
    pz = cube[1] - cube[0]
    mz = cube[2] - cube[0]
    imgshp = pz.shape
    pz_std = pz.std(axis=1, keepdims=True)
    pz_normalized = (pz - pz.mean(axis=1, keepdims=True)) / pz_std
    mz_std = mz.std(axis=1, keepdims=True)
    mz_normalized = (mz - mz.mean(axis=1, keepdims=True)) / mz_std
    f_pz = fft.fft(pz_normalized, imgshp[1] * 2 - 1, axis=-1)
    f_mz = fft.fft(mz_normalized, imgshp[1] * 2 - 1, axis=-1)
    cc = fft.ifft(f_pz * f_mz.conj(), imgshp[1] * 2 - 1, axis=-1) / imgshp[1]
    cc = cc.real
    cc = np.roll(cc, imgshp[1] - 1, axis=1)
    mean_cc = np.nanmean(cc, axis=0)
    return mean_cc


if __name__ == '__main__':
    # x = np.array([25.48723033, 24.03133495, 16.45724938, 30.53108401, 29.81265493, 17.77434783,
    #               26.17120737, 31.202132, 24.98253694, 19.3927809, 16.32724126, 21.8978995,
    #               15.23840869, 17.38753519, 24.22666197, 23.9136671, 21.08465408, 24.55670289,
    #               15.038733, 20.94466288, 27.9509734, 12.12873431, 27.53776742, 14.32068209,
    #               22.75790612, 19.13854881])

    x = np.array([25.48723033, 24.03133495, 16.45724938, 30.53108401, 29.81265493, 17.77434783,
                  26.17120737, 31.202132, 24.98253694, 19.3927809, 16.32724126, 21.8978995,
                  15.23840869, 17.38753519, 24.22666197, 23.9136671, 21.08465408, 24.55670289,
                  15.038733, 25, 27.9509734, 12.12873431, 27.55, 15,
                  25, 25])

    # x = np.ones_like(x) * 25
    # x[7]= 30 #he ii basically
    # x[19] = 25.5
    # x[25] = 26

    file_path = pathlib.Path(__file__).parents[1] / 'meit_synth_cube_full.sav'
    sav1 = scipy.io.readsav(file_path)
    synth_cube = sav1.moses_synth_cube
    eit_dem = sav1.eit_dem

    # because of convolving eit with MOSES point spread funtion some zeros need to be removed
    zero_trim = slice(5, -5)
    synth_cube = synth_cube[..., zero_trim, :]

    plot_cube = synth_cube[:, 1] * 10 ** eit_dem[..., None, None]
    # slicer = CubeSlicer(plot_cube, origin='lower', vmax=np.percentile(plot_cube, 99))
    # plt.show()

    file_path3 = pathlib.Path(__file__).parents[1] / 'moses_super.sav'
    sav3 = scipy.io.readsav(file_path3)
    p = sav3.super_plus
    m = sav3.super_minus
    z = sav3.super_zero

    moses_cube = np.array([z, p, m])
    mcor = mean_row_cc(moses_cube[:, zero_trim, ...])

    start = time.time()

    # test_fit = merit(x,synth_cube,mcor)
    # eit_dem = x
    # bounds = [(eit_dem[i] - 5, eit_dem[i] + 5) for i in range(eit_dem.size)]
    #
    # # fit = scipy.optimize.minimize(merit,eit_dem,args=(synth_cube,mcor),bounds=bounds,
    # #                               # options=dict(maxcor=26,eps = 1)
    # #                               )
    #
    # fit = scipy.optimize.differential_evolution(merit, bounds, args=(synth_cube, mcor), workers=4, polish=True)
    # print(fit.x)
    # x = fit.x

    x = [31.08201397, 33.80327896, 21.79443764, 16.7202813, 15.30675487, 24.26310046,
         24.25079284, 14.45231917]
    x = [41.08036155, 36.42910867, 31.77785579, 27.12660291, 22.881341,   19.65105648,
     16.42077197, 13.19048746, 12.0304087,  12.42298426, 12.81555982, 13.20813538,
     11.46301245,  9.36160644,  7.26020044,  5.84213516,  8.52411422, 11.20609329,
     13.88807235, 15.42180299, 15.42453572, 15.42726846, 15.43000119, 13.35222533,
     10.44224604,  7.53226674]

    spread = 1
    bounds = [(x[i] - spread, x[i] + spread) for i in range(len(x))]
    n = 26

    # x = scipy.optimize.brute(triangle_dem_merit, bounds,args=(synth_cube, mcor,n), workers=-1
    #                            , Ns=5)
    # fit = scipy.optimize.minimize(triangle_dem_merit,x,args=(synth_cube,mcor,n),
    #                               bounds=bounds,
    #                               options=dict(maxfun=30000, maxiter=30000)
    #                               )

    # fit = scipy.optimize.differential_evolution(triangle_dem_merit, bounds, args=(synth_cube, mcor, n), workers=6,
    #                                             tol=.001,
    #                                             maxiter=1000,
    #                                             popsize=100,
    #                                             disp=True)
    # print(fit.x)
    # print(fit.success)
    # print(fit.message)
    # x = dem_from_triangles(n, fit.x)

    # x = [28.60054131, 28.80963392, 26.70727614, 20.42273053, 19.48168889, 23.61541075,23.43688185, 19.15181321]
    # x = [30, 28, 27, 26, 24.9, 23, 22, 22, 22, 22.7, 22.3, 22.4, 22, 22.2]
    # x = [30,26,25,23,22.6,22.9,22.6]
    x =[42.08036155, 37.42910867, 31.57393247, 26.12660291, 23.09646747, 19.1234502,
     15.52013682, 12.21448187, 11.0304087, 11.54155543, 11.81555982, 14.20587266,
     12.46299089, 9.6960638, 7.24475971, 5.84213522, 8.49316942, 11.33261624,
     13.45302588, 16.41559511, 15.29675936, 15.62000712, 15.72922864, 12.35222533,
     9.81445679, 7.80544828]

    x = dem_from_triangles(n, x)

    print('Time=', time.time() - start)
    print('Best Fit=', merit(x, synth_cube, mcor))


    fit_cube = synth_cube * 10 ** x[..., None, None, None]
    fit_cube = np.sum(fit_cube, axis=0)

    fit_cube /= np.median(fit_cube, axis=(-1, -2), keepdims=True)
    fit_pz = fit_cube[1] - fit_cube[0]
    fit_mz = fit_cube[2] - fit_cube[0]
    moses_pz = p - z
    moses_mz = m - z

    # fit_pz *= moses_pz.max()/fit_pz.max()
    fit_pz *= moses_pz.std() / fit_pz.std()

    vmax = .25
    fig, ax = plt.subplots(2, 1, sharex=True, sharey=True)
    ax[0].imshow(fit_pz, vmin=-vmax, vmax=vmax, origin='lower')
    ax[0].set_title('Fit (Plus-Zero)')
    ax[1].imshow(moses_pz, vmin=-vmax, vmax=vmax, origin='lower')
    ax[1].set_title('MOSES (Plus-Zero)')

    fig, ax = plt.subplots()
    ax.plot(np.arange(4, 6.6, .1), x)
    ax.set_ylabel('Log DEM (K cm$^{-5}$)')
    ax.set_xlabel('Log T (K)')

    tcor = mean_row_cc(fit_cube)
    lag = np.arange(mcor.shape[0]) - mcor.shape[0] // 2
    fig, ax = plt.subplots()
    ax.plot(lag, mcor, label='MOSES CC Function')
    ax.plot(lag, tcor, color='r', label='Fit CC Function')
    ax.set_ylabel('Cross Correlation')
    ax.set_xlabel('Lag (pix)')
    ax.legend()

    plt.show()
