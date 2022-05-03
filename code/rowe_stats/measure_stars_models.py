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
import pdb
#
PROCESS = int(environ['SLURM_PROCID']) #process ID for a single cpu
NTASKS = int(environ['SLURM_NTASKS']) #total number of processes running

#logging.basicConfig(filename='measurement.log', encoding='utf-8', level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p:')
#RNG
rng = np.random.RandomState(seed=666)
stampsize=48

#assume given image and psf paths:
#finalcut = pd.read_csv('./finalcut_query.csv')
#prefix='/decade/decarchive/'

def load_psf_and_image(path_to_image, name_of_image, path_to_psf, starlist, psfmodel,path_to_cat,name_of_cat):
    image = galsim.fits.read(path_to_image+name_of_image)
    weight = galsim.fits.read(path_to_image+name_of_image,hdu=3)
    starlist = fits.open(path_to_psf+starlist)
    des_psfex = galsim.des.DES_PSFEx(path_to_psf+psfmodel,path_to_image+name_of_image)
    image_fits = fits.open(path_to_image+name_of_image)
    wcs_pixel_world = WCS(image_fits[1].header) 
    ccdcoords = findall(r'\d+',image_fits['sci'].header['detsec']) 
    min_x_pix, min_y_pix = int(ccdcoords[0]),int(ccdcoords[2])
    cat = fits.open(path_to_cat+name_of_cat)       
    return image, weight, starlist, des_psfex, wcs_pixel_world, min_x_pix, min_y_pix, cat

def match_star_to_cat(starlist,cat,pixeldistance=1.0):
    good_index = get_psf_stars_index(starlist)
    goodstars = starlist[2].data[good_index]
    #now get the X and Y CCD coordinates of each of the "good" stars
    X,Y = goodstars['x_image'],goodstars['y_image']
    Xcat,Ycat = cat[2].data['x_image'],cat[2].data['y_image'] 
    matched_star_indices=np.array([],dtype=int)
    match_fail=0
    for x,y in zip(X,Y):
        #will match stars in starlist by their X,Y position
        wherex,wherey = np.isclose(x,Xcat,atol=pixeldistance),np.isclose(y,Ycat,atol=pixeldistance)
        product = wherex*wherey
        if np.sum(product)!=1: 
            match_fail=match_fail+1
            continue
        else:
            matched_star_indices = np.append(matched_star_indices,np.where(product)[0])
    print('Failed to find a match in the sextractor fullcat for %1.2f percent of the starlist stars'%(100*len(match_fail)/len(good_index)))
    return matched_star_indices

def get_index_of_star_in_full_catalog(Xstar,Ystar,cat,pixeldistance=1.0):
    #get the X and Y of all objects in catalog and compare with the input
    Xcat,Ycat = cat[2].data['x_image'],cat[2].data['y_image'] 

    wherex,wherey = np.isclose(Xstar,Xcat,atol=pixeldistance),np.isclose(Ystar,Ycat,atol=pixeldistance) 
    product = wherex*wherey 
    if np.sum(product)!=1: 
            return np.nan
    else:
        indout = np.where(product)[0] #this is the index IN THE FULL CATALOG of the star with coords Xstar, Ystar
        if np.isscalar(indout):
            return indout
        else:
            return indout[0]
    
def check_if_star_should_have_been_masked(starlist, cat):
    #get the starlist entry, match it to sextractor catalog imaflags_iso and see if it has a mask
    matched_star_indices = match_star_to_cat(starlist,cat)
    matched_imaflags_iso = cat[2].data['imaflags_iso'][matched_star_indices]
    print('Flags found for matched stars:',np.unique(matched_imaflags_iso))
    return 0

def get_flux_auto_from_cat(starlist,cat):
    matched_star_indices = match_star_to_cat(starlist,cat)
    matched_flux_auto = cat[2].data['flux_auto'][matched_star_indices]
    return matched_flux_auto

def get_psf_stars_index(starlist):
	goodstars= np.where(starlist[2].data['flags_psf']==0)[0] 
	return goodstars


def measure_shear_of_ngmix_obs(obs,prefix,i):
    am = Admom(rng=rng)
    res = am.go(obs, 0.3)
    #pdb.set_trace()#understand what's inside res
    if res['flags'] != 0:
        return np.nan, np.nan, np.nan
    else:
        return res['e1'],res['e2'],res['T']
    #        print('admom flagged object %d in: %s with flags %s'%(i,prefix,flagstr(res['flags'])))
    '''
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
    '''

def measure_hsm_shear():
    



def get_band_name(subdirectories):
    if len(subdirectories)!=3:
        raise ValueError('Did not find only psf_X/ and X/ and cat_X/ subdirectories here')
    if len(subdirectories[0])==1:
        bandname = subdirectories[0] #assumes the images are in a directory called simply 'g' or 'r' for instance!
    if len(subdirectories[1])==1:
        bandname = subdirectories[0]
    else:
        bandname = subdirectories[2]
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
for expname in exps_for_this_process[0:2]: #loops over exposures!
    expnum = int(expname[3:])
    if expnum in np.loadtxt(output_location+'DONE_EXPS.txt'):
        print('PROCESS %d will not do exposure %s cause it was already done!'%(PROCESS,expname))
        continue
    #LOOP OF THE TYPE "for expname in exps_for_this_process"
    rootdir = location+expname+'/' #'/home/secco/project2-kicp-secco/delve/rowe_stats_files/exp145973/'
    band = get_band_name(listdir(rootdir)) #finds what band is in this exposure
    print('PROCESS %d doing %s (%s-band)'%(PROCESS,expname,band))  
    outputfile_name =output_location+band+'/'+band+'band_'+expname+'.fits.fz'

    path_to_image = rootdir+band+'/' #exp145973/r/ for instance
    path_to_psf = rootdir+'psf_'+band+'/'#exp145973/psf_r/ for instance
    path_to_cat = rootdir+'cat_'+band+'/'#exp145973/cat_r/ for instance
    time1=time()
    #now initiate some arrays where the outputs will be appended, then eventually written to fits
    focal_x_out, focal_y_out = np.array([]), np.array([])
    pix_x_out, pix_y_out = np.array([]), np.array([])
    ra_out, dec_out = np.array([]), np.array([])
    g1_star_out, g2_star_out, T_star_out, g1_model_out, g2_model_out, T_model_out = np.array([]), np.array([]),np.array([]), np.array([]),np.array([]), np.array([])
    mag_auto_out, imaflags_iso_out = np.array([]),np.array([])
    N_failed_stars = 0
    for name_of_image in listdir(path_to_image): #loops over the CCDs of an exposure!
        prefix = name_of_image[0:25] #the prefix containing expnum, band and ccdnum
        #print('doing ',prefix)
        #outputfile_name =output_location+band+'/'+band+'band_'+prefix[0:-1]+'.txt'
        

        #outputfile = open(outputfile_name,'w')
        #outputfile.write('#focal_x focal_y pix_x pix_y ra dec g1_star g2_star T_star g1_model g2_model T_model\n')
        starlist = prefix+'psfex-starlist.fits'
        psfmodel = prefix+'psfexcat.psf'
        name_of_cat = prefix+'red-fullcat.fits'

        image, weight, starlist, des_psfex, wcs_pixel_world, min_x, min_y,cat = load_psf_and_image(path_to_image,
                                                                name_of_image,
                                                                path_to_psf,
                                                                starlist,
                                                                psfmodel,
                                                                path_to_cat,
                                                                name_of_cat)
        goodstar=get_psf_stars_index(starlist)
        Ngoodstar = len(goodstar)
        #pdb.set_trace()
        if Ngoodstar<100:
            print('This ccd has less than 100 PSF stars: flag it.')
            flag_bad_ccds = open(output_location+'FLAGGED_CCDS.txt','a')
            flag_bad_ccds.write(name_of_image+'\n')
            flag_bad_ccds.close()

         #returns stars that have psf_flags==0 and for which a match has been found
        tmp_focal_x, tmp_focal_y =np.array([]), np.array([])
        tmp_pix_x, tmp_pix_y = np.array([]), np.array([])
        tmp_ra, tmp_dec = np.array([]), np.array([])
        tmp_g1_star, tmp_g2_star, tmp_T_star = np.array([]), np.array([]), np.array([])
        tmp_g1_model, tmp_g2_model, tmp_T_model = np.array([]), np.array([]), np.array([])
        tmp_mag_auto = np.array([])
        tmp_imaflags_iso = np.array([])
        #print('found %d stars that pass flags'%len(goodstar))
        #ig = 0
        for goodstar_index in goodstar:

            X = starlist[2].data['x_image'].astype(int)[goodstar_index] 
            Y = starlist[2].data['y_image'].astype(int)[goodstar_index]
            X_float = starlist[2].data['x_image'][goodstar_index] #getting them as floats also (as in the file itself)
            Y_float = starlist[2].data['y_image'][goodstar_index]

            #first, check if this star with PSF_FLAGS==0 has a match in the sextractor catalog:
            location_in_catalog = get_index_of_star_in_full_catalog(X_float,Y_float,cat,pixeldistance=1.0)
            #print('location_in_catalog=',location_in_catalog)
            #pdb.set_trace()
            if np.isnan(location_in_catalog):
                print('Did not find a match for this star')
                continue
            #if a match was found, get flags and mag_auto
            IMAFLAG_ISO=cat[2].data['imaflags_iso'][location_in_catalog]
            MAG_AUTO = cat[2].data['mag_auto'][location_in_catalog]
            
            #now continue: re-bound the image around theright location
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

            tmp_focal_x = np.append(tmp_focal_x, X_float+min_x)
            tmp_focal_y = np.append(tmp_focal_y, Y_float+min_y)
            tmp_pix_x = np.append(tmp_pix_x, X_float)
            tmp_pix_y = np.append(tmp_pix_y, Y_float)
            tmp_ra = np.append(tmp_ra, ra)
            tmp_dec = np.append(tmp_dec, dec)
            tmp_g1_star = np.append(tmp_g1_star, g1_star)
            tmp_g2_star = np.append(tmp_g2_star, g2_star)
            tmp_T_star = np.append(tmp_T_star, T_star)
            tmp_g1_model = np.append(tmp_g1_model, g1_model)
            tmp_g2_model = np.append(tmp_g2_model, g2_model)
            tmp_T_model = np.append(tmp_T_model, T_model)
            tmp_mag_auto = np.append(tmp_mag_auto, MAG_AUTO)
            tmp_imaflags_iso = np.append(tmp_imaflags_iso, IMAFLAG_ISO)
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
        imaflags_iso_out = np.append(imaflags_iso_out, tmp_imaflags_iso) 
        mag_auto_out = np.append(mag_auto_out, tmp_mag_auto) 
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
    c13 = fits.Column(name='IMAFLAGS_ISO', array=imaflags_iso_out, format='e')
    c14 = fits.Column(name='MAG_AUTO', array=mag_auto_out, format='e')
    t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14])
    time2=time()
    t.header['FAILED']=(N_failed_stars, 'number of failed ngmix measurements')
    time_it_took = (time2-time1)/60.0
    t.header['RUNTIME'] = ( time_it_took, 'minutes to run this exposure' )
    #print('should be done with one exposure')
    #pdb.set_trace() 
    t.writeto(outputfile_name,overwrite=True)
    print('PROCESS %d DONE: wrote %s to  %s (took %1.2f minutes)'%(PROCESS,expname,outputfile_name,time_it_took))

    track_whats_done = open(output_location+'DONE_EXPS.txt','a')
    track_whats_done.write(str(expnum)+'\n')
    track_whats_done.close()
    
    #NOW DELETE ALL FILES!
    if environ['DELETE_DIR']=='True':
        delete_all_rsynced_files(location,expname)






	
