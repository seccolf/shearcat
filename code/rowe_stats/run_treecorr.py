import treecorr
import numpy as np
from astropy.io import fits
import sys
from os import listdir,environ
###############
#this script assumes you'll star 4 slurm tasks(processes), each of which may have many CPUs
#and each of which will run treecorr on a single band
###############

PROCESS = int(environ('SLURM_PROCID')) 
band_dict = {1:'g', 2:'r', 3:'i', 4:'z'}
band = band_dict[PROCESS]
print('PROCESS %d is running the %s-band'%(PROCESS,band))

#find files
location_of_exposures = '/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/'+band+'/'
all_exposures = listdir(location_of_exposures)
#where to output?
location_of_output = '/home/secco/SHEAR/shearcat/code/rowe_stats/ouput_rhos/'

#combine all outputs into a single array
ra = np.array([])
dec = np.array([])
g1_star = np.array([])
g1_model = np.array([])
g2_star = np.array([])
g2_model = np.array([])
for expame in all_exposures:
	f = fits.open(location_of_exposures+expame)
	ra = np.append(ra, f[1].data['ra'])
	dec = np.append(dec, f[1].data['dec'])
	g1_star = np.append(g1_star, f[1].data['g1_star'])
	g1_model = np.append(g1_model, f[1].data['g1_model'])
	g2_star = np.append(g2_star, f[1].data['g2_star'])
	g2_model = np.append(g2_model, f[1].data['g2_model'])

print('Total number of objects for %d exposures in band %s: %d'%(len(all_exposures),band,len(ra)))

#now create the residual vector and pass everything to treecorr
q1 = g1_star-g1_model
q2 = g2_star-g2_model

cat_q = treecorr.Catalog(g1=q1, g2=q2, ra=ra, dec=dec, ra_units='deg',dec_units='deg')
cat_e = treecorr.Catalog(g1=g1_model, g2=g2_model, ra=ra, dec=dec, ra_units='deg',dec_units='deg')

GG = treecorr.GGCorrelation(nbins=20,min_sep=0.1,max_sep=250.0,sep_units='arcmin',verbose=3,bin_slop=1.0)
GG.process(cat_e)
GG.write(location_of_output+'rho0_%sband_bslop%1._%dexposures.txt'%(band, bin_slop,len(all_exposures)))
print('Done rho0 in band %s'%band)
GG.process(cat_q)
GG.write(location_of_output+'rho1_%sband_bslop%1.2f_%dexposures.txt'%(band, bin_slop,len(all_exposures)))
print('Done rho1 in band %s'%band)
GG.process(cat_e,cat_q)
GG.write(location_of_output+'rho2_%sband_bslop%1.2f_%dexposures.txt'%(band, bin_slop,len(all_exposures)))
print('Done rho2 in band %s'%band)
