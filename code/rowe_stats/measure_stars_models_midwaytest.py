#use this to measure the stars' and psf models' ellipticities

#stuff needed to load PSF files and measure shapes 
import galsim
from galsim.des import DES_PSFEx
import ngmix
from ngmix.fitting import Fitter as LMSimple
from ngmix.admom import AdmomFitter as Admom
from ngmix.flags import get_flags_str as flagstr
#usual stuff
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import logging
import pandas as pd
from os import listdir
from time import time
from re import findall
#
logging.basicConfig(filename='test.log', encoding='utf-8', level=logging.INFO)

#RNG
rng = np.random.RandomState(seed=42)

stampsize=24

#assume given image and psf paths:
#finalcut = pd.read_csv('./finalcut_query.csv')
#prefix='/decade/decarchive/'

def load_psf_and_image(path_to_image, name_of_image, path_to_psf, starlist, psfmodel):
        image = galsim.fits.read(path_to_image+name_of_image)
        weight = galsim.fits.read(path_to_image+name_of_image,hdu=3)
        starlist = fits.open(path_to_psf+starlist)
        des_psfex = galsim.des.DES_PSFEx(path_to_psf+psfmodel,path_to_image+name_of_image)
        image_fits = fits.open(path_to_image+name_of_image)
        wcs_pixel_world = WCS(image_fits[1].header) 
        ccdcoords = findall(r'\d+',image_fits['sci'].header['detsec']) 
        min_x_pix, min_y_pix = int(ccdcoords[0]),int(ccdcoords[2])       
        return image, weight, starlist, des_psfex, wcs_pixel_world, min_x_pix, min_y_pix

def get_psf_stars_index(starlist):
	goodstars= np.where(starlist[2].data['flags_psf']==0)[0] 
	return goodstars


def measure_shear_of_ngmix_obs(obs,prefix,i):
    am = Admom(rng=rng)
    res = am.go(obs, 0.3)
    if res['flags'] != 0:
        logging.info('admom flagged object %d in: '%i,prefix,' ',flagstr(res['flags']))
    lm = LMSimple('gauss')
    lm_res = lm.go(obs, res['pars'])
    #print(lm_res['flags'])
    if lm_res['flags'] == 0:
        g1 = lm_res['pars'][2]
        g2 = lm_res['pars'][3]
        T = lm_res['pars'][4]
        return g1, g2, T
    else: 
        logging.info('lm flagged object %d in: '%i,prefix,' ',flagstr(res['flags']))
        return np.nan, np.nan, np.nan




#########SOME LOOP

rootdir = '/home/secco/project2-kicp-secco/delve/rowe_stats_files/exp145973/'
path_to_image = rootdir+'r/'
path_to_psf = rootdir+'psf_r/'

for name_of_image in listdir(path_to_image):
    time1=time()  
    prefix = name_of_image[0:25] #the prefix containing expnum, band and ccdnum
    print('doing ',prefix)
    outputfile_name =rootdir+'measurements/'+prefix[0:-1]+'.txt'
    outputfile = open(outputfile_name,'w')
    outputfile.write('#focal_x focal_y pix_x pix_y ra dec g1_star g2_star T_star g1_model g2_model T_model\n')
    starlist = prefix+'psfex-starlist.fits'
    psfmodel = prefix+'psfexcat.psf'
    image, weight, starlist, des_psfex, wcs_pixel_world, min_x, min_y = load_psf_and_image(path_to_image,
                                                            name_of_image,
                                                            path_to_psf,
                                                            starlist,
                                                            psfmodel)
    goodstar=get_psf_stars_index(starlist)
    print('found %d stars that pass flags'%len(goodstar))
    for goodstar_index in goodstar:

        X = starlist[2].data['x_image'].astype(int)[goodstar_index] 
        Y = starlist[2].data['y_image'].astype(int)[goodstar_index]
        X_float = starlist[2].data['x_image'][goodstar_index] #getting them as floats also (as in the file itself)
        Y_float = starlist[2].data['y_image'][goodstar_index]

        newbounds = galsim.BoundsI(X-stampsize/2,X+stampsize/2,Y-stampsize/2,Y+stampsize/2)

        image_cutout = image[newbounds].array
        weight_cutout = weight[newbounds].array

        #position where we want the PSF
        psf_pos = galsim.PositionD(X, Y)
        psf_model = des_psfex.getPSF(psf_pos)

        copy_stamp = image[newbounds].copy() #copies the galsim object with wcs and everything
        psf_image=psf_model.drawImage(image=copy_stamp,method='no_pixel')
        psf_cutout = psf_image.array

        #get the wcs at the location 
        psf_wcs = des_psfex.getLocalWCS(psf_pos)

        ra,dec = wcs_pixel_world.pixel_to_world_values(X_float, Y_float)
        ra = ra.item()
        dec=dec.item()
        
        #now create the ngmix observations
        star_obs = ngmix.Observation(
        	image=image_cutout,
        	weight=weight_cutout,
        	jacobian=ngmix.Jacobian(row=stampsize/2 , col=stampsize/2 , wcs=psf_wcs))

        psf_model_obs = ngmix.Observation(
        	image=psf_cutout,
        	weight=np.ones(psf_cutout.shape),
        	jacobian=ngmix.Jacobian(row=stampsize/2 , col=stampsize/2 , wcs=psf_wcs))

        g1_star, g2_star, T_star = measure_shear_of_ngmix_obs(star_obs,prefix,goodstar_index)
        g1_model, g2_model, T_model = measure_shear_of_ngmix_obs(psf_model_obs,prefix,goodstar_index)
        outputfile.write('%f %f %f %f %f %f %f %f %f %f %f %f\n'%(X_float+min_x, Y_float+min_y,
                                                                  X_float, Y_float,
                                                                  ra, dec,
                                                                  g1_star, g2_star, T_star,
                                                                  g1_model, g2_model, T_model))
    outputfile.close()
    time2=time()
    print('wrote ',prefix,'to ',outputfile_name,'(took %1.2f seconds)\n'%(time2-time1))
    








	
