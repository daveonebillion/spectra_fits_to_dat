import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.io.fits as fits
from astropy.convolution import convolve, Box1DKernel
import glob

"""Extracts UVES spectra from fits files and turns them into .dat files"""

def filewriter(wavelength, flux, error, save_path, filename):
    """
    write the spectrum to file in savepath/file name in space-separated columns of wavelength, flux, flux_error.
    """
    if not os.path.exists(save_path):
        os.makedirs(save_path)    
    fl=open((save_path+filename),'w') 
    for w, f, e in zip(wavelength, flux, error):
        fl.write('%f %g %g\n'%(w,f,e))
        
def plot_spectrum(wavelength, flux, error, rootname):
    """diagnostic plot. Smooths with a 5pt boxcar. Plots the raw spectrum in grey, then overlays another spectrum with flagged points removed"""
    flux = convolve(flux,Box1DKernel(5))
    plt.figure(rootname)
    plt.plot(wavelength, flux, '0.5', label='spectrum')
    plt.plot(wavelength, error, label='error')  
    plt.legend()
    plt.show()
    plt.close()

    
def heliocentric(w, dv):
    """applies a heliocentric correction to the wavelength axis """
    c = 2.998e8
    w *= (1.0+(dv/c))
    return w

def wavelength_builder(header):
    """builds the wavelength axis using the information in the fits header"""
    w = []
    start, step, length = header['CRVAL1'], header['CDELT1'], header['NAXIS1']
    w = [start+step*i for i in range(length)]
    return np.array(w, dtype=float)

def uves_fits_to_dat(plot=True, file_path=os.getcwd()+'/', save_path=os.getcwd()+'/spectra/'):
    """looks for all datsets in filepath and extracts BLUE, REDL and REDU spectra from each"""
    datasets = os.listdir(file_path)
    colours = ['BLUE', 'REDL', 'REDU']
    for colour in colours:
        for dataset in datasets:
            sci_file = glob.glob(file_path+dataset+'/*FLUXCAL_SCI_POINT_'+colour+'.fits')
            if sci_file != []:
                hdul = fits.open(sci_file[0])
                header = hdul[0].header
                f = hdul[0].data
                hdul.close()
                w = wavelength_builder(header)
                w = heliocentric(w, header['HIERARCH ESO QC VRAD HELICOR'])
                hdul_e = fits.open(glob.glob(file_path+dataset+'/*FLUXCAL_ERRORBAR_SCI_POINT_'+colour+'.fits')[0])
                e = hdul_e[0].data
                hdul_e.close()
                
                filename = header['OBJECT']+'_'+header['INSTRUME']+'_'+colour+'_'+header['DATE-OBS']+'.dat'
                filewriter(w, f, e, save_path+colour+'/', filename)
                

                #make a plot
                if plot == True:
                    plot_spectrum(w, f, e, filename)


                
#test
uves_fits_to_dat(file_path='../csc_storage/wd0137-349/observations/reflex_end_products/2015-10-13T16:42:39/', save_path='../csc_storage/wd0137-349/spectra')
    
    
    

