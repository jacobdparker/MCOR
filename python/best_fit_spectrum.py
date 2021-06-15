import kgpy.chianti
import scipy.interpolate
import numpy as np
import pathlib
import matplotlib.pyplot as plt


def dem_from_triangles(n, x):
    temp = np.arange(n, dtype='float64')
    dem = np.zeros_like(temp)
    x = np.array(x)
    fig, ax = plt.subplots()
    for i in range(len(x)):
        step = n / (x.shape[0] - 1)
        shift = i * step
        triangle = (1 - np.abs(temp - shift) / step) * x[i]
        triangle[triangle < 0] = 0
        dem += triangle
        ax.plot(triangle)
    return dem


if __name__ == '__main__':

    n = 26
    x = [34.49302866, 28.21622703, 30.7782344, 28.30659669, 21.35836151, 19.9022056,
         20.75010404, 25.14890224, 26.25373862, 19.18684885, 24.53881626]
    x = [30, 28, 27, 26, 25, 23, 22, 22, 22, 22, 23, 22, 22,22]
    x = [41.08036155, 23.80427943, 11.80607981, 13.2642176,   5.45899529, 15.42063182,
     15.43078197,  4.62228745]
    x =[42.08036155, 37.42910867, 31.57393247, 26.12660291, 23.09646747, 19.1234502,
     15.52013682, 12.21448187, 11.0304087, 11.54155543, 11.81555982, 14.20587266,
     12.46299089, 9.6960638, 7.24475971, 5.84213522, 8.49316942, 11.33261624,
     13.45302588, 16.41559511, 15.29675936, 15.62000712, 15.72922864, 12.35222533,
     9.81445679, 7.80544828]
    x = dem_from_triangles(n, x)
    # print(x-x_0)

    # Load in MOSES throughput curves
    throughput_dir = pathlib.Path(__file__).parents[1] / 'mosesI_throughput'
    filter1 = np.genfromtxt(throughput_dir / 'filter1.csv')
    filter2 = np.genfromtxt(throughput_dir / 'filter2.csv')
    grating = np.genfromtxt(throughput_dir / 'mosesI_grating.csv')
    secondary = np.genfromtxt(throughput_dir / 'mosesI_secondary.csv')

    filter1_interp = scipy.interpolate.interp1d(filter1[:, 0] * 10, filter1[:, 1])
    filter2_interp = scipy.interpolate.interp1d(filter2[:, 0] * 10, filter2[:, 1])
    grating_interp = scipy.interpolate.interp1d(grating[:, 0], grating[:, 1])
    secondary_interp = scipy.interpolate.interp1d(secondary[:, 0], secondary[:, 1])

    wvl = np.arange(280, 337, .1)

    combined_throughput = (filter1_interp(wvl) * filter2_interp(wvl)) ** 2 * grating_interp(wvl) * secondary_interp(wvl)

    dem_logT = np.arange(4, 6.6, .1)
    temperature = 10 ** dem_logT
    dlnt = np.log(10 ** (dem_logT[1] - dem_logT[0]))
    dt = temperature * dlnt

    em = np.ones_like(x)
    # em = 10**x
    pressure = 1e15
    dens = pressure / temperature
    abund_file = 'sun_coronal_1992_feldman'
    wvl_range = [wvl[0], wvl[-1]]
    # ionlist = ['fe_13', 'si_9', 'si_11', 'mg_8', 'he_2']

    path = pathlib.Path(__file__).parent / 'moses_synthetic_spectrum.pickle'

    if not path.is_file():
        bunch = kgpy.chianti.Bunch(temperature, dens, wvl_range, em=em, abundance=abund_file,
                                   allLines=0,
                                   verbose=True,
                                   minAbund=1e-5,
                                   # ionList=ionlist
                                   )
        bunch.to_pickle(path)
    else:
        bunch = kgpy.chianti.Bunch.from_pickle(path)

    wave = bunch.Intensity['wvl']
    mask = (wave < wvl_range[1]) & (wave > wvl_range[0])
    wave = wave[mask]
    sort = wave.argsort()
    wave = wave[sort]

    unsummed_int = bunch.Intensity['intensity'][:, mask] * 10 ** (x[..., None])
    ions = bunch.Intensity['ionS'][mask]

    print(x)
    integrated = np.trapz(unsummed_int[:, sort], temperature, axis=0)
    throughput_interp = scipy.interpolate.interp1d(wvl, combined_throughput)
    integrated = integrated * throughput_interp(wave)

    integrated = np.trapz(unsummed_int[:, sort], temperature, axis=0)
    throughput_interp = scipy.interpolate.interp1d(wvl, combined_throughput)
    integrated = integrated * throughput_interp(wave)

    relative_int = integrated / integrated.max()

    fig, ax = plt.subplots()
    normalized_throughput = combined_throughput / max(combined_throughput)
    ax.plot(wvl, normalized_throughput, label='Normalized Throughput')
    ax.vlines(wave, np.zeros_like(wave), relative_int, color='r')
    ax.legend()
    ax.set_ylabel('Intensity Relative to He II 303.78 $\AA$')
    ax.set_xlabel('Wavelength ($\AA$)')

    top_lines = 10
    intensity_mask = integrated > np.sort(integrated)[::-1][top_lines]
    relative_int = relative_int[intensity_mask]
    wave = wave[intensity_mask]
    ions = ions[sort][intensity_mask]
    print(ions)
    print(relative_int)
    print(wave)

    # ions = kgpy.chianti.ion_tolatex(ions, use_latex=False)
    # for i in range(wave.shape[0]):
    #     ax.text(wave[i], relative_int[i], ions[i], ha='center', va='bottom', rotation=45)

    plt.show()
