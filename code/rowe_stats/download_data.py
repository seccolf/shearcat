import pandas as pd 
import os
from time import time

df = pd.read_csv('./finalcut_query.csv')
location = '/home/secco/project2-kicp-secco/delve/rowe_stats_files/'

expnums = df['EXPNUM']
rbands = df[df['BAND']=='r']
ibands = df[df['BAND']=='i']
zbands = df[df['BAND']=='z']

#first, let's download a large number of r-band exposures and their PSF files:
prefix = '/decade/decarchive/'
for expnum in rbands['EXPNUM'][0:2]:
	time1=time()
	selected = rbands[rbands['EXPNUM']==expnum]
	path_to_im = prefix+selected['PATH']+'/'
	path_to_psf = prefix+selected['PATH.1']+'/'
	basedir = location+'exp'+str(expnum)
	basedir_image = basedir+'/r/'
	basedir_psf = basedir+'/psf_r/'
	createdirs = 'mkdir '+basedir+' '+basedir_image+' '+basedir_psf
	print(createdirs)
	os.system(createdirs)
	print('Downloading exposure %d images...'%expnum)
	scp_command_im = 'scp lsecco@deslogin.cosmology.illinois.edu:'+path_to_im+'* '+basedir_image
	time2=time()
	print('... took %1.2f minutes'%(time2-time1)/60.0)
	scp_command_psf = 'scp lsecco@deslogin.cosmology.illinois.edu:'+path_to_psf+'* '+basedir_psf
	time3=time()
	print('Total download time including psf files for exposure %d is %1.2f minutes\n'%(expnum, (time3-time1)/60.0))

