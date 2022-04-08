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
from os import listdir,environ,system
from time import time
from re import findall
#import pdb
#
PROCESS = int(environ['SLURM_PROCID']) #process ID for a single cpu
NTASKS = int(environ['SLURM_NTASKS']) #total number of processes running

#logging.basicConfig(filename='measurement.log', encoding='utf-8', level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p:')
#RNG
rng = np.random.RandomState(seed=666)
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
    #if res['flags'] != 0:
    #        print('admom flagged object %d in: %s with flags %s'%(i,prefix,flagstr(res['flags'])))
    lm = LMSimple('gauss')
    try:
            lm_res = lm.go(obs, res['pars'])
            if lm_res['flags'] == 0:
                    g1 = lm_res['pars'][2]
                    g2 = lm_res['pars'][3]
                    T = lm_res['pars'][4]
                    return g1, g2, T
            else:
                    #print('lm flagged object %d in: %s with flags %s'%(i,prefix,flagstr(lm_res['flags'])))
                    return np.nan, np.nan, np.nan
    except:
            #print("ngmix error in object %d in: %s"%(i,prefix))
            return np.nan, np.nan, np.nan


def get_band_name(subdirectories):
    if len(subdirectories)!=2:
        raise ValueError('Did not find only psf_X/ and X/ subdirectories here')
    if len(subdirectories[0])==1:
        bandname = subdirectories[0] #assumes the images are in a directory called simply 'g' or 'r' for instance!
    else:
        bandname = subdirectories[1]
    return bandname

def delete_all_rsynced_files(location_del,expname_del):
    location_of_stuff = location_del+expname_del+'/'    
    if environ['DELETE_DIR']=='True':
        command_del = 'rm -r '+location_of_stuff+'*'
        command_del_dir = 'rmdir '+location_of_stuff
        print('DELETE: ',command_del)
        print('DELETE: ',command_del_dir)
        system(command_del)
        system(command_del_dir)
    else:
        print('Will not delete inputs in ',location_of_stuff)
    
#########SOME LOOP
location = '/home/secco/project2-kicp-secco/delve/rowe_stats_files/'
output_location = '/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/'
all_exposures = listdir(location) #will get the names of all exposures as strings 'expXXXXXX'
number_of_exps = len(all_exposures)

expnumber_shared = round(number_of_exps/NTASKS +0.5)
exps_for_this_process= all_exposures[ int(PROCESS*expnumber_shared) : int((PROCESS+1)*expnumber_shared) ]
print('PROCESS %d will take care of exposures '%PROCESS,exps_for_this_process)
for expname in exps_for_this_process: #loops over exposures!
    expnum = int(expname[3:])
    if expnum in np.loadtxt(output_location+'DONE_EXPS.txt'):
        print('PROCESS %d will not do exposure %s cause it was already done!'%(PROCESS,expname))
        continue
    #LOOP OF THE TYPE "for expname in exps_for_this_process"
    rootdir = location+expname+'/' #'/home/secco/project2-kicp-secco/delve/rowe_stats_files/exp145973/'
    band = get_band_name(listdir(rootdir)) #finds what band is in this exposure
    print('PROCESS %d doing %s (%s-band)'%(PROCESS,expname,band))  

    path_to_image = rootdir+band+'/' #exp145973/r/ for instance
    path_to_psf = rootdir+'psf_'+band+'/'#exp145973/psf_r/ for instance
    time1=time()
    #now initiate some arrays where the outputs will be appended, then eventually written to fits
    focal_x_out, focal_y_out = np.array([]), np.array([])
    pix_x_out, pix_y_out = np.array([]), np.array([])
    ra_out, dec_out = np.array([]), np.array([])
    g1_star_out, g2_star_out, T_star_out, g1_model_out, g2_model_out, T_model_out = np.array([]), np.array([]),np.array([]), np.array([]),np.array([]), np.array([])
    N_failed_stars = 0
    for name_of_image in listdir(path_to_image): #loops over the CCDs of an exposure!
        prefix = name_of_image[0:25] #the prefix containing expnum, band and ccdnum
        #print('doing ',prefix)
        #outputfile_name =output_location+band+'/'+band+'band_'+prefix[0:-1]+'.txt'
        outputfile_name =output_location+band+'/'+band+'band_'+expname+'.fits.fz'

        #outputfile = open(outputfile_name,'w')
        #outputfile.write('#focal_x focal_y pix_x pix_y ra dec g1_star g2_star T_star g1_model g2_model T_model\n')
        starlist = prefix+'psfex-starlist.fits'
        psfmodel = prefix+'psfexcat.psf'
        image, weight, starlist, des_psfex, wcs_pixel_world, min_x, min_y = load_psf_and_image(path_to_image,
                                                                name_of_image,
                                                                path_to_psf,
                                                                starlist,
                                                                psfmodel)
        goodstar=get_psf_stars_index(starlist)
        Ngoodstar = len(goodstar)
        tmp_focal_x, tmp_focal_y = np.ones(Ngoodstar), np.ones(Ngoodstar)
        tmp_pix_x, tmp_pix_y = np.ones(Ngoodstar), np.ones(Ngoodstar)
        tmp_ra, tmp_dec = np.ones(Ngoodstar), np.ones(Ngoodstar)
        tmp_g1_star, tmp_g2_star, tmp_T_star = np.ones(Ngoodstar), np.ones(Ngoodstar), np.ones(Ngoodstar) 
        tmp_g1_model, tmp_g2_model, tmp_T_model = np.ones(Ngoodstar), np.ones(Ngoodstar), np.ones(Ngoodstar)

        #print('found %d stars that pass flags'%len(goodstar))
        ig = 0
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
            #when g1_star and etc cannot be measured, do not output focal, ra etc
            #but this cannot be the entire problem cause the number of nonzero ras and decs are pretty small
            #print(g1_star, g2_star, T_star, g1_model,g2_model,T_model,ra,dec,'<-- g1,g2,T,g1,g2,T,ra,dec')
            ngmix_outputs = np.array([g1_star, g2_star, T_star, g1_model, g2_model, T_model])
            if np.any(np.isnan(ngmix_outputs)):
                N_failed_stars=N_failed_stars+1

            tmp_focal_x[ig] = X_float+min_x
            tmp_focal_y[ig] = Y_float+min_y
            tmp_pix_x[ig] = X_float
            tmp_pix_y[ig] = Y_float
            tmp_ra[ig] = ra
            tmp_dec[ig] = dec
            tmp_g1_star[ig] = g1_star
            tmp_g2_star[ig] = g2_star
            tmp_T_star[ig] = T_star
            tmp_g1_model[ig] = g1_model
            tmp_g2_model[ig] = g2_model
            tmp_T_model[ig] = T_model
            ig=ig+1
            
            #if ig>300:
            #        pdb.set_trace()

        focal_x_out = np.append(focal_x_out,tmp_focal_x)
        focal_y_out = np.append(focal_y_out,tmp_focal_y)
        pix_x_out = np.append(pix_x_out, tmp_pix_x)
        pix_y_out = np.append(pix_y_out, tmp_pix_y)
        ra_out = np.append(ra_out, tmp_ra)
        dec_out = np.append(dec_out, tmp_dec)
        g1_star_out = np.append(g1_star_out, tmp_g1_star)
        g2_star_out = np.append(g2_star_out, tmp_g2_star)
        T_star_out = np.append(T_star_out, tmp_T_star)
        g1_model_out = np.append(g1_model_out, tmp_g1_model)
        g2_model_out = np.append(g2_model_out, tmp_g2_model)
        T_model_out = np.append(T_model_out, tmp_T_model) 
        #print('should be done with one ccd')
        #pdb.set_trace()
        #outputfile.close()

    c1 = fits.Column(name='focal_x', array=focal_x_out, format='e')
    c2 = fits.Column(name='focal_y', array=focal_y_out, format='e')
    c3 = fits.Column(name='pix_x', array=pix_x_out, format='e')
    c4 = fits.Column(name='pix_y', array=pix_y_out, format='e')
    c5 = fits.Column(name='ra', array=ra_out, format='e')
    c6 = fits.Column(name='dec', array=dec_out, format='e')
    c7 = fits.Column(name='g1_star', array=g1_star_out, format='e')
    c8 = fits.Column(name='g2_star', array=g2_star_out, format='e')
    c9 = fits.Column(name='T_star', array=T_star_out, format='e')
    c10 = fits.Column(name='g1_model', array=g1_model_out, format='e')
    c11 = fits.Column(name='g2_model', array=g2_model_out, format='e')
    c12 = fits.Column(name='T_model', array=T_model_out, format='e')
    t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
    time2=time()
    t.header['FAILED']=(N_failed_stars, 'number of failed ngmix measurements')
    time_it_took = (time2-time1)/60.0
    t.header['RUNTIME'] = ( time_it_took, 'minutes to run this exposure' )
    #print('should be done with one exposure')
    #pdb.set_trace() 
    t.writeto(outputfile_name,overwrite=True)
    print('PROCESS %d DONE: wrote %s to eg. %s (took %1.2f minutes)'%(PROCESS,prefix,outputfile_name,time_it_took))

    track_whats_done = open(output_location+'DONE_EXPS.txt','a')
    track_whats_done.write(str(expnum)+'\n')
    track_whats_done.close()
    
    #NOW DELETE ALL FILES!
    delete_all_rsynced_files(location,expname)






	
