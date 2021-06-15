import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
import pandas
import astropy.units as u
import ChiantiPy

# The goal of this test file is to match the results of Chiantipy to that of ch_ss in IDL to make sure we know how
# everything works.

if __name__ == '__main__':

    # dem_file = '/home/jake/chianti/dem/quiet_sun.dem'
    # dem = pandas.read_csv(dem_file, sep=' ', skipinitialspace=True, skipfooter=9, names=['logT', 'EM'])
    #
    #
    # dem_logT = dem['logT'].to_numpy()
    # temperature = 10 ** dem_logT
    # dlnt = np.log(10 ** (dem_logT[1] - dem_logT[0]))
    # dt = temperature * dlnt
    #
    # em = 10 ** dem['EM'].to_numpy()
    # minabund = 1e-5

    temperature = 1e5
    pressure = 1e15
    dens = pressure / temperature
    em = 10 ** 25
    ion_list = ['o_3', 'o_4', 'o_5', 'he_1', 'mg_10']
    abund_file = 'sun_coronal_2012_schmelz'
    wvl_range = [557, 630]


    test_bunch = ch.bunch(temperature, dens, wvl_range, em=em, abundance=abund_file,
                          allLines=0,
                          verbose=True, ionList=ion_list,
                          # minAbund=minabund,  # Note, minAbund will override ionList
                          )

    wvl = test_bunch.Intensity['wvl']
    print(wvl.shape)
    mask = (wvl < wvl_range[1]) & (wvl > wvl_range[0])
    wvl = wvl[mask]
    print(wvl.shape)




    int = test_bunch.Intensity['integrated']
    unsummed = test_bunch.Intensity['intensity']

    int = int[mask]
    unsummed = unsummed[:,mask]
    sort = int.argsort()
    sort = sort[::-1]
    print(int[sort][:10])
    print(np.sum(unsummed[:,sort]*dt[...,None],axis=0)[:10])
    integrated = np.trapz(unsummed[:, sort], temperature, axis=0)[:10]
    print(integrated/integrated.max())
    print(wvl[sort][:10])




    # test_spectrum = ch.spectrum(temperature, dens, wavelength=wvl_range, abundance=abund_file, verbose=1, allLines=0,
    #                             # ionList=ion_list, #ionList did not work, threw an error
    #                             minAbund=minabund, doContinuum=False)

    # plot_test = test_spectrum.spectrumPlot(integrated=True) # Seems to do nothing
    # plt.show()