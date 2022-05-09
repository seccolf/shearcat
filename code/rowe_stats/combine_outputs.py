from astropy.io import fits
from os import environ,listdir
import numpy as np

#organize processes in nodes
PROCESS = int(environ['SLURM_PROCID']) 
band_dict = {0:'g', 1:'r', 2:'i', 3:'z'}
band = band_dict[PROCESS]
print('PROCESS %d is running the %s-band'%(PROCESS,band),flush=True)

#find files
location_of_exposures = '/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/'+band+'/'
all_exposures = listdir(location_of_exposures)
#where to output the full catalog?
location_of_output = '/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/'

#combine all outputs into a single array
focal_x = np.array([])
focal_y = np.array([])
pix_x = np.array([])
pix_y = np.array([])
ra = np.array([])
dec = np.array([])

g1_star = np.array([])
g1_model = np.array([])
T_star = np.array([])
g2_star = np.array([])
g2_model = np.array([])
T_model = np.array([]) 

g1_star_hsm = np.array([])
g1_model_hsm = np.array([])
T_star_hsm = np.array([])
g2_star_hsm = np.array([])
g2_model_hsm = np.array([])
T_model_hsm = np.array([]) 

MAG_AUTO = np.array([])
IMAFLAGS_ISO = np.array([])

N_failed_stars = 0
N_failed_CCDS = 0
N_failed_match = 0 

for expame in all_exposures:
	f = fits.open(location_of_exposures+expame)
	focal_x = np.append(focal_x, f[1].data['focal_x'])
	focal_y = np.append(focal_y, f[1].data['focal_y'])
	pix_x = np.append(pix_x, f[1].data['pix_x'])
	pix_y = np.append(pix_y, f[1].data['pix_y'])
	ra = np.append(ra, f[1].data['ra'])
	dec = np.append(dec, f[1].data['dec'])
	g1_star = np.append(g1_star, f[1].data['g1_star'])
	g1_model = np.append(g1_model, f[1].data['g1_model'])
	T_star = np.append(T_star, f[1].data['T_star'])
	T_model = np.append(T_model, f[1].data['T_model'])
	g2_star = np.append(g2_star, f[1].data['g2_star'])
	g2_model = np.append(g2_model, f[1].data['g2_model'])

	g1_star_hsm = np.append(g1_star_hsm, f[1].data['g1_star_hsm'])
	g1_model_hsm = np.append(g1_model_hsm, f[1].data['g1_model_hsm'])
	T_star_hsm = np.append(T_star_hsm, f[1].data['T_star_hsm'])
	T_model_hsm = np.append(T_model_hsm, f[1].data['T_model_hsm'])
	g2_star_hsm = np.append(g2_star_hsm, f[1].data['g2_star_hsm'])
	g2_model_hsm = np.append(g2_model_hsm, f[1].data['g2_model_hsm'])

	MAG_AUTO = np.append(MAG_AUTO,f[1].data['MAG_AUTO'])
	IMAFLAGS_ISO = np.append(IMAFLAGS_ISO,f[1].data['IMAFLAGS_ISO'])

	N_failed_stars = N_failed_stars + f[1].header['FAILED_stars']
	N_failed_CCDS = N_failed_CCDS + f[1].header['FAILED_ccds']
	N_failed_match = N_failed_match + f[1].header['FAILED_nomatch']

print('Total number of objects for %d exposures in band %s: %d'%(len(all_exposures),band,len(ra)),flush=True)

c1 = fits.Column(name='focal_x', array=focal_x, format='e')
c2 = fits.Column(name='focal_y', array=focal_y, format='e')
c3 = fits.Column(name='pix_x', array=pix_x, format='e')
c4 = fits.Column(name='pix_y', array=pix_y, format='e')
c5 = fits.Column(name='ra', array=ra, format='e')
c6 = fits.Column(name='dec', array=dec, format='e')
c7 = fits.Column(name='g1_star', array=g1_star, format='e')
c8 = fits.Column(name='g2_star', array=g2_star, format='e')
c9 = fits.Column(name='T_star', array=T_star, format='e')
c10 = fits.Column(name='g1_model', array=g1_model, format='e')
c11 = fits.Column(name='g2_model', array=g2_model, format='e')
c12 = fits.Column(name='T_model', array=T_model, format='e')

c13 = fits.Column(name='g1_star_hsm', array=g1_star_hsm, format='e')
c14 = fits.Column(name='g2_star_hsm', array=g2_star_hsm, format='e')
c15 = fits.Column(name='T_star_hsm', array=T_star_hsm, format='e')
c16 = fits.Column(name='g1_model_hsm', array=g1_model_hsm, format='e')
c17 = fits.Column(name='g2_model_hsm', array=g2_model_hsm, format='e')
c18 = fits.Column(name='T_model_hsm', array=T_model_hsm, format='e')

c19 = fits.Column(name='MAG_AUTO', array=MAG_AUTO, format='e')
c20 = fits.Column(name='IMAFLAGS_ISO', array=IMAFLAGS_ISO, format='e')

t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20])
#time2=time()
t.header['FAILED_stars']=(N_failed_stars, 'number of failed ngmix measurements')
t.header['FAILED_ccds']=(N_failed_CCDS, 'number of failed ccds (<100 stars)')
t.header['FAILED_nomatch']=(N_failed_match, 'number of stars not matched to sextractor cat')
#time_it_took = (time2-time1)/60.0
#t.header['RUNTIME'] = ( time_it_took, 'minutes to run this exposure' )
#print('should be done with one exposure')
#pdb.set_trace() 
outputfile_name = location_of_output+band+'band_%dexps.fits.fz'
t.writeto(outputfile_name,overwrite=False)
#print('PROCESS %d DONE: wrote %s to  %s (took %1.2f minutes)'%(PROCESS,expname,outputfile_name,time_it_took))
