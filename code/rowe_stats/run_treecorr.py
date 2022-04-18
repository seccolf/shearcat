import treecorr
import numpy as np
from astropy.io import fits
import sys
from os import listdir,environ
import pandas as pd
###############
#this script assumes you'll star 4 slurm tasks(processes), each of which may have many CPUs
#and each of which will run treecorr on a single band
#but if you run with less than 4 processes, you'll only get some (but not all) of the bands (see band_dict below)
###############

#cut exposures by some threshold in T_eff to see if rho stats go down appreciably:
#query output file:
def remove_exposures_by_teff(list_of_exposures,threshold):
	df = pd.read_csv('/home/secco/SHEAR/shearcat/code/rowe_stats/finalcut_query_withradec_withTeff.csv')
	teff_mask = df['T_EFF']>threshold 
	good_exposures = df['EXPNUM'][teff_mask]
	#now, get the actual exposure number from the input exposures and check if they are in good
	#the input list_of_exposures above gives me items like 'rband_exp273484.fits.fz'
	pass_mask = np.array([], dtype=bool)
	for exp_name in list_of_exposures:
		expnumber = int(re.findall(r'\d+',exp_name)[0])
		print(exp_name, expnumber,'<-- testing')
		true_or_false = expnumber in good_exposures.array 
		print(true_or_false)
		pass_mask=np.append(pass_mask,true_or_false)
	print(pass_mask)
	list_of_exposures = np.array(list_of_exposures)
	return list_of_exposures[pass_mask]

############
#treecorr setup:
############
bslop= 0.01
min_angle = 0.1 #arcmin
max_angle = 250.0 #arcmin
nbins = 25
teff_threshold=0.3
############

#organize processes in nodes
PROCESS = int(environ['SLURM_PROCID']) 
band_dict = {0:'g', 1:'r', 2:'i', 3:'z'}
band = band_dict[PROCESS]
print('PROCESS %d is running the %s-band'%(PROCESS,band),flush=True)

#find files
location_of_exposures = '/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/'+band+'/'
all_exposures_premasking = listdir(location_of_exposures)
all_exposures = remove_exposures_by_teff(all_exposures_premasking,teff_threshold) #removes exposures with Teff<0.3
print('Removed %1.2f%  of %s-band exposures based on Teff>%1.1f'%( 100*(1.0-len(all_exposures)/len(all_exposures_premasking)),band,teff_threshold),flush=True)
#where to output?
location_of_output = '/home/secco/SHEAR/shearcat/code/rowe_stats/output_rhos/'

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

print('Total number of objects for %d exposures in band %s: %d'%(len(all_exposures),band,len(ra)),flush=True)

#now create the residual vector and pass everything to treecorr
q1 = g1_star-g1_model
q2 = g2_star-g2_model

cat_q = treecorr.Catalog(g1=q1, g2=q2, ra=ra, dec=dec, ra_units='deg',dec_units='deg')
cat_e = treecorr.Catalog(g1=g1_model, g2=g2_model, ra=ra, dec=dec, ra_units='deg',dec_units='deg')

GG = treecorr.GGCorrelation(nbins=nbins,min_sep=min_angle,max_sep=max_angle,sep_units='arcmin',verbose=3,bin_slop=bslop)
GG.process(cat_e)
GG.write(location_of_output+'rho0_%sband_bslop%1.3f_%dexposures_teff%1.1f.txt'%(band, bslop,len(all_exposures),teff_threshold))
print('Done rho0 in band %s'%band,flush=True)
GG.process(cat_q)
GG.write(location_of_output+'rho1_%sband_bslop%1.3f_%dexposures_teff%1.1f.txt'%(band, bslop,len(all_exposures),teff_threshold))
print('Done rho1 in band %s'%band,flush=True)
GG.process(cat_e,cat_q)
GG.write(location_of_output+'rho2_%sband_bslop%1.3f_%dexposures_teff%1.1f.txt'%(band, bslop,len(all_exposures),teff_threshold))
print('Done rho2 in band %s'%band,flush=True)
