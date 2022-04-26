import pandas as pd 
import os
from time import time
import sys
import logging
###############
#read band, start and finish from command line
band=str(sys.argv[1]) #band
start=int(sys.argv[2]) #start of the loop
end=int(sys.argv[3]) #end of the loop
###############

logging.basicConfig(filename=band+'band.log', encoding='utf-8', level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p:')
logging.info('Running: %s %s %s %s'%(sys.argv[0],sys.argv[1],sys.argv[2],sys.argv[3]))

df = pd.read_csv('./finalcut_query_withradec.csv')
location = '/home/secco/project2-kicp-secco/delve/rowe_stats_files/'

band_select = (df['BAND']==band)
box_select = (df['RADEG']>120.0)*(df['RADEG']<250.0)*(df['DECDEG']>-40.0)*(df['DECDEG']<40.0)
final_select = band_select*box_select
exposures=df[final_select]# apply total masking
expnums = exposures['EXPNUM']

prefix = '/decade/decarchive/'
for expnum in expnums[start:end]:
        time1=time()
        selected = exposures[exposures['EXPNUM']==expnum]
        path_to_im = prefix+selected['PATH'].item()+'/'
        path_to_psf = prefix+selected['PATH.1'].item()+'/'
        path_to_cat = prefix+selected['PATH.1'].item()[0:-3]+'cat/'
        basedir = location+'exp'+str(expnum)
        basedir_image = basedir+'/'+band+'/'
        basedir_psf = basedir+'/psf_'+band+'/'
        basedir_cat = basedir+'/cat_'+band+'/'
        createdirs = 'mkdir '+basedir+' '+basedir_image+' '+basedir_psf+' '+basedir_cat
        #print(createdirs)
        os.system(createdirs)
        logging.info('Downloading exposure %d images (%s-band)...'%(expnum,band))
        #print('Downloading exposure %d images (%s-band)...'%(expnum,band))
        #scp_command_im = 'scp lsecco@deslogin.cosmology.illinois.edu:'+path_to_im+'* '+basedir_image
        #print(scp_command_im)
        #os.system(scp_command_im)        
        rsync_command_im = 'rsync -a lsecco@deslogin.cosmology.illinois.edu:'+path_to_im+' '+basedir_image
        logging.info(rsync_command_im)
        #print(rsync_command_im)        
        os.system(rsync_command_im) 
        time2=time()
        #logging.info('... took %1.2f minutes'%((time2-time1)/60.0))
        #print('... took %1.2f minutes'%((time2-time1)/60.0))
        rsync_command_psf = 'rsync -a lsecco@deslogin.cosmology.illinois.edu:'+path_to_psf+' '+basedir_psf
        os.system(rsync_command_psf)	
        time3=time()
        #logging.info('Total download time including psf files for exposure %d is %1.2f minutes\n'%(expnum, (time3-time1)/60.0))
        #print('Total download time including psf files for exposure %d is %1.2f minutes\n'%(expnum, (time3-time1)/60.0))
        rsync_command_cat = 'rsync -a lsecco@deslogin.cosmology.illinois.edu:'+path_to_cat+' '+basedir_cat
        os.system(rsync_command_cat)    
        time4=time()
        logging.info('Total download time including cat,im,psf files for exposure %d is %1.2f minutes\n'%(expnum, (time4-time1)/60.0))
        #print('Total download time including psf files for exposure %d is %1.2f minutes\n'%(expnum, (time3-time1)/60.0))

