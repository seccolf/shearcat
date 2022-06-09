#use this to measure the stars' and psf models' ellipticities

#stuff needed to load PSF files and measure shapes 
import galsim
from galsim.des import DES_PSFEx
import ngmix
from ngmix.fitting import Fitter as LMSimple
from ngmix.admom import AdmomFitter as Admom
from ngmix.flags import get_flags_str as flagstr
from ngmix import priors, joint_prior
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
save_individual_CCDs = True #set to True only for debugging purposes
do_flags = False #set to False ONLY FOR DEBUGGING PURPOSES

location = '/home/secco/project2-kicp-secco/delve/rowe_stats_files/' #where to look for exposures
output_location = '/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/problematic_exposures/' #where to write results
#all_exposures = listdir(location) #will get the names of all exposures as strings 'expXXXXXX'
all_exposures = ['exp640670', 'exp830026','exp830031','exp736962']
number_of_exps = len(all_exposures)


#logging.basicConfig(filename='measurement.log', encoding='utf-8', level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p:')
#RNG
rng = np.random.RandomState(seed=500)
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
    fwhm=image_fits['SCI'].header['FWHM']*0.26 
    wcs_pixel_world = WCS(image_fits[1].header) 
    ccdcoords = findall(r'\d+',image_fits['sci'].header['detsec']) 
    min_x_pix, min_y_pix = int(ccdcoords[0]),int(ccdcoords[2])
    cat = fits.open(path_to_cat+name_of_cat)       
    return image, weight, starlist, des_psfex, wcs_pixel_world, min_x_pix, min_y_pix, cat, fwhm

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
    print('Failed to find a match in the sextractor fullcat for %1.2f percent of the starlist stars'%(100*len(match_fail)/len(good_index)),flush=True)
    return matched_star_indices

def get_index_of_star_in_full_catalog(Xstar,Ystar,cat,pixeldistance=4.0):
    #pixel distance is 4 because we don't want any overlap between detected objects (and PSF has size approx 4 pixels on average)
    #get the X and Y of all objects in catalog and compare with the input
    Xcat,Ycat = cat[2].data['x_image'],cat[2].data['y_image'] 

    wherex,wherey = np.isclose(Xstar,Xcat,atol=pixeldistance),np.isclose(Ystar,Ycat,atol=pixeldistance) 
    product = wherex*wherey 
    #pdb.set_trace()
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
    print('Flags found for matched stars:',np.unique(matched_imaflags_iso),flush=True)
    return 0

def get_flux_auto_from_cat(starlist,cat):
    matched_star_indices = match_star_to_cat(starlist,cat)
    matched_flux_auto = cat[2].data['flux_auto'][matched_star_indices]
    return matched_flux_auto

def get_psf_stars_index(starlist):
	goodstars= np.where(starlist[2].data['flags_psf']==0)[0] 
	return goodstars


def measure_shear_of_ngmix_obs(obs,prefix,i,fwhm):
    am = Admom(rng=rng)
    T_guess = (fwhm / 2.35482)**2 * 2.
    res = am.go(obs, T_guess)
    #pdb.set_trace()#understand what's inside res
    pdb.set_trace()
    if res['flags'] != 0:
        return np.nan, np.nan, np.nan
    else:
        g1,g2 = e1e2_to_g1g2(res['e1'],res['e2'])
        if abs(g1) > 0.5 or abs(g2) > 0.5: #bad measurement
            return np.nan, np.nan, np.nan
        else:
            return g1,g2,res['T']
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
def e1e2_to_g1g2(e1, e2): #from ngmix base code
    """
    convert e1,e2 to reduced shear style ellipticity
    parameters
    ----------
    e1,e2: tuple of scalars
        shapes in (ixx-iyy)/(ixx+iyy) style space
    outputs
    -------
    g1,g2: scalars
        Reduced shear space shapes
    """
    ONE_MINUS_EPS=0.99999999999999999
    e = np.sqrt(e1 * e1 + e2 * e2)
    if isinstance(e1, np.ndarray):
        (w,) = np.where(e >= 1.0)
        if w.size != 0:
            raise Error("some e were out of bounds")
        eta = np.arctanh(e)
        g = np.tanh(0.5 * eta)
        np.clip(g, 0.0, ONE_MINUS_EPS, g)
        g1 = np.zeros(g.size)
        g2 = np.zeros(g.size)
        (w,) = np.where(e != 0.0)
        if w.size > 0:
            fac = g[w] / e[w]
            g1[w] = fac * e1[w]
            g2[w] = fac * e2[w]
    else:
        if e >= 1.0:
            raise Error("e out of bounds: %s" % e)
        if e == 0.0:
            g1, g2 = 0.0, 0.0
        else:
            eta = np.arctanh(e)
            g = np.tanh(0.5 * eta)
            if g >= 1.0:
                # round off?
                g = ONE_MINUS_EPS
            fac = g / e
            g1 = fac * e1
            g2 = fac * e2
    return g1, g2

def measure_hsm_shear(im, wt):
    MAX_CENTROID_SHIFT=1.0 #1 pixel recentering at most
    shape_data = im.FindAdaptiveMom(weight=wt, strict=False) #might need to turn on strict
    if shape_data.moments_status == 0:
        dx = shape_data.moments_centroid.x - im.true_center.x
        dy = shape_data.moments_centroid.y - im.true_center.y
        if dx**2 + dy**2 > MAX_CENTROID_SHIFT**2:
            #print('cetroid changed by too much in HSM, will return nan',flush=True)
            return np.nan, np.nan, np.nan
        else:
            e1 = shape_data.observed_shape.e1
            e2 = shape_data.observed_shape.e2
            s = shape_data.moments_sigma
            jac = im.wcs.jacobian(im.true_center)
            M = np.matrix( [[ 1 + e1, e2 ], [ e2, 1 - e1 ]] ) * s*s
            J = jac.getMatrix()
            M = J * M * J.T

            e1 = (M[0,0] - M[1,1]) / (M[0,0] + M[1,1])
            e2 = (2.*M[0,1]) / (M[0,0] + M[1,1])
            T = M[0,0] + M[1,1]

            shear = galsim.Shear(e1=e1, e2=e2)
            g1 = shear.g1
            g2 = shear.g2
            return g1, g2, T
    else:
        return np.nan, np.nan, np.nan 

def get_band_name(subdirectories):
    if len(subdirectories)!=3:
        raise ValueError('Did not find only psf_X/ and X/ and cat_X/ subdirectories here')
    if len(subdirectories[0])==1:
        bandname = subdirectories[0] #assumes the images are in a directory called simply 'g' or 'r' for instance!
    if len(subdirectories[1])==1:
        bandname = subdirectories[1]
    if len(subdirectories[2])==1:
        bandname = subdirectories[2]
    return bandname

def delete_all_rsynced_files(location_del,expname_del):
    location_of_stuff = location_del+expname_del+'/'    
    if environ['DELETE_DIR']=='True':
        command_del = 'rm -r '+location_of_stuff+'*'
        command_del_dir = 'rmdir '+location_of_stuff
        print('DELETE: ',command_del,flush=True)
        print('DELETE: ',command_del_dir,flush=True)
        system(command_del)
        system(command_del_dir)
    else:
        print('Will not delete inputs in ',location_of_stuff,flush=True)
    
def write_measurements_to_fits(OUTFITS_name,focal_x_out,focal_y_out,
                    pix_x_out,pix_y_out,ra_out,dec_out,
                    g1_star_out,g2_star_out,T_star_out,
                    g1_model_out,g2_model_out,T_model_out,
                    g1_star_hsm_out,g2_star_hsm_out,T_star_hsm_out,
                    g1_model_hsm_out,g2_model_hsm_out,T_model_hsm_out,
                    imaflags_iso_out,mag_auto_out,N_failed_stars,N_failed_CCDS,N_bad_match):
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

    c13 = fits.Column(name='g1_star_hsm', array=g1_star_hsm_out, format='e')
    c14 = fits.Column(name='g2_star_hsm', array=g2_star_hsm_out, format='e')
    c15 = fits.Column(name='T_star_hsm', array=T_star_hsm_out, format='e')
    c16 = fits.Column(name='g1_model_hsm', array=g1_model_hsm_out, format='e')
    c17 = fits.Column(name='g2_model_hsm', array=g2_model_hsm_out, format='e')
    c18 = fits.Column(name='T_model_hsm', array=T_model_hsm_out, format='e')

    c19 = fits.Column(name='IMAFLAGS_ISO', array=imaflags_iso_out, format='e')
    c20 = fits.Column(name='MAG_AUTO', array=mag_auto_out, format='e')
    t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20])
    t.header['FAILED_stars']=(N_failed_stars, 'number of failed ngmix measurements')
    t.header['FAILED_ccds']=(N_failed_CCDS, 'number of failed ccds')
    t.header['FAILED_badmatch']=(N_bad_match, 'number of stars not matched to sextractor cat')
    t.writeto(OUTFITS_name,overwrite=True)
    return None

#########START OF ACTUAL LOOPS
expnumber_shared = round(number_of_exps/NTASKS +0.5)
exps_for_this_process= all_exposures[ int(PROCESS*expnumber_shared) : int((PROCESS+1)*expnumber_shared) ]
#print('PROCESS %d will take care of exposures '%PROCESS,exps_for_this_process,flush=True)
for expname in exps_for_this_process: #loops over exposures!
    expnum = int(expname[3:])
    if expnum in np.loadtxt(output_location+'DONE_EXPS.txt'):
        print('PROCESS %d will not do exposure %s cause it was already done!'%(PROCESS,expname),flush=True)
        continue

    rootdir = location+expname+'/' #results in eg. '/home/secco/project2-kicp-secco/delve/rowe_stats_files/exp145973/'
    band = get_band_name(listdir(rootdir)) #finds what band is in this exposure
    print('PROCESS %d doing %s (%s-band)'%(PROCESS,expname,band),flush=True)  
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
    g1_star_hsm_out, g2_star_hsm_out, T_star_hsm_out, g1_model_hsm_out, g2_model_hsm_out, T_model_hsm_out = np.array([]), np.array([]),np.array([]), np.array([]),np.array([]), np.array([])
    mag_auto_out, imaflags_iso_out = np.array([]),np.array([])
    N_failed_stars = 0
    N_failed_CCDS = 0 #number of CCDS for this expnum that did not pass flags
    N_bad_match = 0 #whether the star was not found in the sextractor catalog or if it has a close neighbor
    for name_of_image in listdir(path_to_image): #loops over the CCDs of an exposure!
        prefix = name_of_image[0:25] #the prefix containing expnum, band and ccdnum
        #print('doing ',prefix)
        #outputfile_name =output_location+band+'/'+band+'band_'+prefix[0:-1]+'.txt'
        

        #outputfile = open(outputfile_name,'w')
        #outputfile.write('#focal_x focal_y pix_x pix_y ra dec g1_star g2_star T_star g1_model g2_model T_model\n')
        starlist = prefix+'psfex-starlist.fits'
        psfmodel = prefix+'psfexcat.psf'
        name_of_cat = prefix+'red-fullcat.fits'

        image, weight, starlist, des_psfex, wcs_pixel_world, min_x, min_y,cat,fwhm = load_psf_and_image(path_to_image,
                                                                name_of_image,
                                                                path_to_psf,
                                                                starlist,
                                                                psfmodel,
                                                                path_to_cat,
                                                                name_of_cat)
        

        goodstar=get_psf_stars_index(starlist)
        Ngoodstar = len(goodstar)
        #print('Found %d good stars in CCD %s'%(Ngoodstar,name_of_image),flush=True)
        #pdb.set_trace()
        
        #INITIAL FLAGGING: numbers of stars
        #perform basic flagging of the CCD before running measurements over stars
        if do_flags:
            if Ngoodstar<100: #FLAGGING 1: number of PSF stars is too small
                print('!!!!FLAG CCD %s: has less than 100 PSF stars'%name_of_image,flush=True)
                flag_bad_ccds = open(output_location+'FLAGGED_CCDS.txt','a')
                flag_bad_ccds.write(name_of_image+'\n')
                flag_bad_ccds.close()
                N_failed_CCDS=N_failed_CCDS+1
                continue
            if Ngoodstar>int(0.25*np.sum(cat[2].data['IMAFLAGS_ISO']==0)):#FLAGGING 2: too many PSF stars
                print('!!!!FLAG CCD %s: more than 25%% of the objects in image are stars'%name_of_image,flush=True)
                flag_bad_ccds = open(output_location+'FLAGGED_CCDS.txt','a')
                flag_bad_ccds.write(name_of_image+'\n')
                flag_bad_ccds.close()
                N_failed_CCDS=N_failed_CCDS+1
                continue

        #starting output arrays for this CCD
        tmp_focal_x, tmp_focal_y =np.array([]), np.array([])
        tmp_pix_x, tmp_pix_y = np.array([]), np.array([])
        tmp_ra, tmp_dec = np.array([]), np.array([])
        tmp_g1_star, tmp_g2_star, tmp_T_star = np.array([]), np.array([]), np.array([])
        tmp_g1_model, tmp_g2_model, tmp_T_model = np.array([]), np.array([]), np.array([])
        tmp_g1_star_hsm, tmp_g2_star_hsm, tmp_T_star_hsm = np.array([]), np.array([]), np.array([])
        tmp_g1_model_hsm, tmp_g2_model_hsm, tmp_T_model_hsm = np.array([]), np.array([]), np.array([])
        tmp_mag_auto = np.array([])
        tmp_imaflags_iso = np.array([])

        for goodstar_index in goodstar: #start loop over good stars

            X = starlist[2].data['x_image'].astype(int)[goodstar_index] 
            Y = starlist[2].data['y_image'].astype(int)[goodstar_index]
            X_float = starlist[2].data['x_image'][goodstar_index] #getting them as floats also (as in the file itself)
            Y_float = starlist[2].data['y_image'][goodstar_index]

            #first, check if this star with PSF_FLAGS==0 has a match in the sextractor catalog:
            location_in_catalog = get_index_of_star_in_full_catalog(X_float,Y_float,cat,pixeldistance=4.0)
            #this returns nan if the star was not found or if it has a neighbor detection within 4 pixels (apprx 1 arcsec)

            #NEXT FLAGGING: star is safely un-blended and was found in the sextractor catalog 
            if np.isnan(location_in_catalog):
                N_bad_match=N_bad_match+1
                continue

            #if a good match was found, get imaflags_iso and mag_auto
            IMAFLAG_ISO=cat[2].data['imaflags_iso'][location_in_catalog]
            MAG_AUTO = cat[2].data['mag_auto'][location_in_catalog]
            
            #now continue: re-bound the image around the right location
            newbounds = galsim.BoundsI(X-stampsize/2 +1,X+stampsize/2 +1,Y-stampsize/2 +1,Y+stampsize/2 +1) #addition of 1 here suggested by Mike J.
            image_cutout = image[newbounds].array
            weight_cutout = weight[newbounds].array
            weight_cutout[weight_cutout<0.0]=0.0

            hsm_input_im = image[newbounds]
            hsm_input_wt = weight[newbounds]
            hsm_input_wt.array[hsm_input_wt.array<0.0]=0.0
            #if np.any(weight_cutout<0.0):
            #    pdb.set_trace()
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
            	#weight=np.ones(psf_cutout.shape),
                weight=weight_cutout, #give the same weights found in the image to the model! (suggested by Mike J.) 
            	jacobian=ngmix.Jacobian(row=stampsize/2 , col=stampsize/2 , wcs=psf_wcs))

            g1_star, g2_star, T_star = measure_shear_of_ngmix_obs(star_obs,prefix,goodstar_index,fwhm)
            g1_model, g2_model, T_model = measure_shear_of_ngmix_obs(psf_model_obs,prefix,goodstar_index,fwhm)
            g1_star_hsm, g2_star_hsm, T_star_hsm = measure_hsm_shear(hsm_input_im, hsm_input_wt)
            #g1_model_hsm, g2_model_hsm, T_model_hsm = measure_hsm_shear(psf_image,None)
            g1_model_hsm, g2_model_hsm, T_model_hsm = measure_hsm_shear(psf_image, hsm_input_wt)#give the same weights found in the image to the model! (suggested by Mike J.) 
 
            #print('g1_star, g2_star, T_star = ',g1_star, g2_star, T_star)
            #print('g1_model, g2_model, T_model =',g1_model, g2_model, T_model)
            #print('g1_star_hsm, g2_star_hsm, T_star_hsm =',g1_star_hsm, g2_star_hsm, T_star_hsm)
            #print('g1_model_hsm, g2_model_hsm, T_model_hsm =',g1_model_hsm, g2_model_hsm, T_model_hsm)

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
            
            tmp_g1_star_hsm = np.append(tmp_g1_star_hsm, g1_star_hsm)
            tmp_g2_star_hsm = np.append(tmp_g2_star_hsm, g2_star_hsm)
            tmp_T_star_hsm = np.append(tmp_T_star_hsm, T_star_hsm)
            tmp_g1_model_hsm = np.append(tmp_g1_model_hsm, g1_model_hsm)
            tmp_g2_model_hsm = np.append(tmp_g2_model_hsm, g2_model_hsm)
            tmp_T_model_hsm = np.append(tmp_T_model_hsm, T_model_hsm)
            
            tmp_mag_auto = np.append(tmp_mag_auto, MAG_AUTO)
            tmp_imaflags_iso = np.append(tmp_imaflags_iso, IMAFLAG_ISO)
            #now go back and get the next goodstar to measure the PSF over
        if save_individual_CCDs:
            print('Saving individual CCDS in: ',output_location+expname+'/')
            CCDOUT_name = output_location+expname+'/'+band+name_of_image+'_.fits.fz'
            write_measurements_to_fits(CCDOUT_name,tmp_focal_x,tmp_focal_y,
                    tmp_pix_x,tmp_pix_y,tmp_ra,tmp_dec,
                    tmp_g1_star,tmp_g2_star,tmp_T_star,
                    tmp_g1_model,tmp_g2_model,tmp_T_model,
                    tmp_g1_star_hsm,tmp_g2_star_hsm,tmp_T_star_hsm,
                    tmp_g1_model_hsm,tmp_g2_model_hsm,tmp_T_model_hsm,
                    tmp_imaflags_iso,tmp_mag_auto,N_failed_stars,N_failed_CCDS,N_bad_match)
        #end loop over good stars (meaning the loop over all stars in a CCD)
        
        #NEXT FLAGGING:
        #more than 3% of the goodstars had some nonzero imaflags_iso:
        #more than 3% of the goodstars were not found in the sextractor catalog or had dangerously close neighbors
        #the standard deviation of the sizes of the final PSF stars is >20% of the mean PSF size 
        do_not_write_ccd = 0
        if np.sum(tmp_imaflags_iso==0)/len(tmp_imaflags_iso) < 0.97:
            print('!!!!FLAG CCD %s: more than 3%% of stars have some nonzero IMAFLAGS_ISO'%name_of_image,flush=True)
            flag_bad_ccds = open(output_location+'FLAGGED_CCDS.txt','a')
            flag_bad_ccds.write(name_of_image+'\n')
            flag_bad_ccds.close()
            do_not_write_ccd=do_not_write_ccd+1
        if N_bad_match/Ngoodstar>0.03:
            print('!!!!FLAG CCD %s: more than 3%% of stars were either blended or not found in sextractor cat'%name_of_image,flush=True)
            flag_bad_ccds = open(output_location+'FLAGGED_CCDS.txt','a')
            flag_bad_ccds.write(name_of_image+'\n')
            flag_bad_ccds.close()
            do_not_write_ccd=do_not_write_ccd+1
        if np.nanstd(tmp_T_model)>0.2*np.nanmean(tmp_T_model):
            print('!!!!FLAG CCD %s: stddev of NGMIX size greater than 0.2*NGMIX mean size'%name_of_image,flush=True)
            flag_bad_ccds = open(output_location+'FLAGGED_CCDS.txt','a')
            flag_bad_ccds.write(name_of_image+'\n')
            flag_bad_ccds.close()
            do_not_write_ccd=do_not_write_ccd+1
        if np.nanstd(tmp_T_model_hsm)>0.2*np.nanmean(tmp_T_model_hsm):
            print('!!!!FLAG CCD %s: stddev of HSM size greater than 0.2*HSM mean size'%name_of_image,flush=True)
            flag_bad_ccds = open(output_location+'FLAGGED_CCDS.txt','a')
            flag_bad_ccds.write(name_of_image+'\n')
            flag_bad_ccds.close()
            do_not_write_ccd=do_not_write_ccd+1
        #pdb.set_trace()
        
        #now write the measurements of this CCD to an output array if the flags immediately above passed
        if do_not_write_ccd>0:
            N_failed_CCDS=N_failed_CCDS+1
        else: 
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

            g1_star_hsm_out = np.append(g1_star_hsm_out, tmp_g1_star_hsm)
            g2_star_hsm_out = np.append(g2_star_hsm_out, tmp_g2_star_hsm)
            T_star_hsm_out = np.append(T_star_hsm_out, tmp_T_star_hsm)
            g1_model_hsm_out = np.append(g1_model_hsm_out, tmp_g1_model_hsm)
            g2_model_hsm_out = np.append(g2_model_hsm_out, tmp_g2_model_hsm)
            T_model_hsm_out = np.append(T_model_hsm_out, tmp_T_model_hsm) 

            imaflags_iso_out = np.append(imaflags_iso_out, tmp_imaflags_iso) 
            mag_auto_out = np.append(mag_auto_out, tmp_mag_auto) 
        #go back and start next CCD

    #ended loop over all CCDs (completed the entire expnum), now apply more flags and write to file

    #MISSING: before writing entire exposure to file, remove any CCDs whose mean size is an outlier compared to the other CCDS in the exposure

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

    c13 = fits.Column(name='g1_star_hsm', array=g1_star_hsm_out, format='e')
    c14 = fits.Column(name='g2_star_hsm', array=g2_star_hsm_out, format='e')
    c15 = fits.Column(name='T_star_hsm', array=T_star_hsm_out, format='e')
    c16 = fits.Column(name='g1_model_hsm', array=g1_model_hsm_out, format='e')
    c17 = fits.Column(name='g2_model_hsm', array=g2_model_hsm_out, format='e')
    c18 = fits.Column(name='T_model_hsm', array=T_model_hsm_out, format='e')

    c19 = fits.Column(name='IMAFLAGS_ISO', array=imaflags_iso_out, format='e')
    c20 = fits.Column(name='MAG_AUTO', array=mag_auto_out, format='e')
    t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20])
    time2=time()
    t.header['FAILED_stars']=(N_failed_stars, 'number of failed ngmix measurements')
    t.header['FAILED_ccds']=(N_failed_CCDS, 'number of failed ccds')
    t.header['FAILED_badmatch']=(N_bad_match, 'number of stars not matched to sextractor cat')
    time_it_took = (time2-time1)/60.0
    t.header['RUNTIME'] = ( time_it_took, 'minutes to run this exposure' )
    #print('should be done with one exposure')
    #pdb.set_trace() 
    t.writeto(outputfile_name,overwrite=True)
    print('PROCESS %d DONE: wrote %s to  %s (took %1.2f minutes)'%(PROCESS,expname,outputfile_name,time_it_took),flush=True)

    track_whats_done = open(output_location+'DONE_EXPS.txt','a')
    track_whats_done.write(str(expnum)+'\n')
    track_whats_done.close()
    
    #NOW DELETE ALL FILES!
    if environ['DELETE_DIR']=='True':
        delete_all_rsynced_files(location,expname)






	
