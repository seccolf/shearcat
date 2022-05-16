import matplotlib.pyplot as pl 
import numpy as np 
from astropy.io import fits
from os import listdir


do_size_histograms_per_exposure = True
do_both_together= False
do_focalplane = False
do_y1_histograms = False
do_stack_chips = False 

use_hsm=True #look for HSM columns if true

naming = '_May9th_hsm'


if do_size_histograms_per_exposure:
	g_measurements_location ='/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/g'
	r_measurements_location ='/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/r'
	i_measurements_location ='/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/i'
	z_measurements_location ='/home/secco/project2-kicp-secco/delve/rowe_stats_measurements/z'

	for location in [g_measurements_location,r_measurements_location,i_measurements_location,z_measurements_location]:
		filenames = listdir(location)
		bandname = location[-1]
		Ts = np.ones(len(filenames))*-999
		Tm = np.ones(len(filenames))*-999
		i=0
		for filename in filenames:
			g=fits.open(location+'/'+filename)
			if use_hsm:
				Ts_mean = np.nanmean(g[1].data['T_star_hsm'])
				Tm_mean = np.nanmean(g[1].data['T_model_hsm'])
			else:	
				Ts_mean = np.nanmean(g[1].data['T_star'])
				Tm_mean = np.nanmean(g[1].data['T_model'])
			
			Ts[i] = Ts_mean
			Tm[i] = Tm_mean
                        i=i+1
		pl.figure()
		pl.suptitle(bandname+'-band exposures '+naming)
		pl.subplot(121)
		pl.xlabel('T_star')
		pl.title('star sizes')
		pl.hist(Ts_mean,range=[0,3],bins=100)
		pl.subplot(122)
		pl.xlabel('T_model')
		pl.title('model sizes')
		pl.hist(Tm_mean,range=[0,3],bins=100)
		pl.savefig('/home/secco/SHEAR/shearcat/code/rowe_stats/figures/'+bandname+'_size_histograms_exposures_hsm.png',dpi=300)
		pl.tight_layout()
		pl.close()



if do_stack_chips:
	for band in ['g', 'r', 'i', 'z']:		
		if band=='g':
			f = fits.open('gband_500exps.fits.fz')
		if band=='r':
			f = fits.open('rband_500exps.fits.fz')
		if band=='i':
			f = fits.open('iband_500exps.fits.fz')
		if band=='z':
			f = fits.open('zband_500exps.fits.fz')

		g1s,g2s,Ts = f[1].data['g1_star'], f[1].data['g2_star'], f[1].data['T_star']
		g1m,g2m,Tm = f[1].data['g1_model'], f[1].data['g2_model'], f[1].data['T_model']

		modelfail=np.isnan(g1m)*np.isnan(g2m)*np.isnan(Tm)
		starfail=np.isnan(g1s)*np.isnan(g2s)*np.isnan(Ts)
		mask = modelfail+starfail
		x,y = f[1].data['pix_x'][~mask], f[1].data['pix_y'][~mask]
		x_max = 2048
		y_max = 4096

		pl.figure(figsize=(10,5))
		pl.title('Stacking chips in %s-band'%band,fontsize=16)
		pl.hist2d(y,x,bins=[int(y_max/2),int(x_max/2)],range=[[0, y_max], [0, x_max]])
		pl.xlabel('Y',fontsize=16)
		pl.ylabel('X',fontsize=16)
		pl.colorbar(label='Number of flagged (good) stars')
		pl.tight_layout()
		pl.savefig('figures/'+band+'band_stacked_chips+naming.png',dpi=100) 
		pl.show()



if do_both_together:
	alpha=0.6
	emin, emax = -0.25,0.25
	tmin, tmax = 0.0, 2.0
	demin,demax =-0.1, 0.1
	dtmin,dtmax =-0.2, 0.2
	bins=100

	f=fits.open('/Users/secco/Documents/projects/SHEAR/shearcat/code/rowe_stats/DES_PSF_catalogs/psf_y1a1-v13.fits')
	for band in ['g', 'r', 'i', 'z']:		
		mask_des = np.where(f[1].data['filter']==band)
		g1s_des,g2s_des,Ts_des = f[1].data['e1'][mask_des], f[1].data['e2'][mask_des], f[1].data['size'][mask_des]
		g1m_des,g2m_des,Tm_des = f[1].data['psf_e1'][mask_des], f[1].data['psf_e2'][mask_des], f[1].data['psf_size'][mask_des]

		dT_des=Ts_des-Tm_des
		dg1_des = g1s_des-g1m_des
		dg2_des = g2s_des-g2m_des

		if band=='g':
			g = fits.open('gband_500exps.fits.fz')
		if band=='r':
			g = fits.open('rband_500exps.fits.fz')
		if band=='i':
			g = fits.open('iband_500exps.fits.fz')
		if band=='z':
			g = fits.open('zband_500exps.fits.fz')

		if use_hsm:
			g1s,g2s,Ts = g[1].data['g1_star_hsm'], g[1].data['g2_star_hsm'], g[1].data['T_star_hsm']
			g1m,g2m,Tm = g[1].data['g1_model_hsm'], g[1].data['g2_model_hsm'], g[1].data['T_model_hsm']
		else:	
			g1s,g2s,Ts = g[1].data['g1_star'], g[1].data['g2_star'], g[1].data['T_star']
			g1m,g2m,Tm = g[1].data['g1_model'], g[1].data['g2_model'], g[1].data['T_model']

		modelfail=np.isnan(g1m)+np.isnan(g2m)+np.isnan(Tm)
		starfail=np.isnan(g1s)+np.isnan(g2s)+np.isnan(Ts)
		mask = modelfail+starfail
		#~mask is "not mask"
		g1s = g1s[~mask]
		g2s = g2s[~mask]
		Ts = Ts[~mask]
		g1m = g1m[~mask]
		g2m = g2m[~mask]
		Tm = Tm[~mask]

		dT=Ts-Tm
		dg1 = g1s-g1m
		dg2 = g2s-g2m


		pl.figure(figsize=(12,9))
		pl.suptitle('DES Y1 and DECADE (%s-band)'%band)

		pl.subplot(331)
		pl.title(r'$e^{\mathrm{star}}_{1}$')
		pl.hist(g1s,bins=bins,range=[emin,emax],density=True,alpha=alpha)
		pl.hist(g1s_des,bins=bins,range=[emin,emax],density=True,alpha=alpha)
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(332)
		pl.title(r'$e^{\mathrm{star}}_{2}$')
		pl.hist(g2s,bins=bins,range=[emin,emax],density=True,alpha=alpha)
		pl.hist(g2s_des,bins=bins,range=[emin,emax],density=True,alpha=alpha)
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(333)
		pl.title(r'$T^{\mathrm{star}}$')
		pl.hist(Ts,bins=bins,range=[tmin,tmax],density=True,alpha=alpha,label='DECADE')
		pl.hist(Ts_des,bins=bins,range=[tmin,tmax],density=True,alpha=alpha,label='DES Y1')
		pl.legend(loc=0,fontsize=16)
		pl.yscale('log')
		#pl.axvline(x=0.0,ls='--',color='k')


		pl.subplot(334)
		pl.title(r'$e^{\mathrm{model}}_{1}$')
		pl.hist(g1m,bins=bins,range=[emin,emax],density=True,alpha=alpha)
		pl.hist(g1m_des,bins=bins,range=[emin,emax],density=True,alpha=alpha)
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(335)
		pl.title(r'$e^{\mathrm{model}}_{2}$')
		pl.hist(g2m,bins=bins,range=[emin,emax],density=True,alpha=alpha)
		pl.hist(g2m_des,bins=bins,range=[emin,emax],density=True,alpha=alpha)
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(336)
		pl.title(r'$T^{\mathrm{model}}$')
		pl.hist(Tm,bins=bins,range=[tmin,tmax],density=True,alpha=alpha)
		pl.hist(Tm_des,bins=bins,range=[tmin,tmax],density=True,alpha=alpha)
		pl.yscale('log')
		#pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(337)
		pl.title(r'$\delta e_1$')
		pl.hist(dg1,bins=bins,range=[demin,demax],density=True,alpha=alpha)
		pl.hist(dg1_des,bins=bins,range=[demin,demax],density=True,alpha=alpha)
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(338)
		pl.title(r'$\delta e_1$')
		pl.hist(dg2,bins=bins,range=[demin,demax],density=True,alpha=alpha)
		pl.hist(dg2_des,bins=bins,range=[demin,demax],density=True,alpha=alpha)
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(339)
		pl.title(r'$\delta T$')
		pl.hist(dT,bins=bins,range=[dtmin,dtmax],density=True,alpha=alpha)
		pl.hist(dT_des,bins=bins,range=[dtmin,dtmax],density=True,alpha=alpha)
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')
		
		
		pl.tight_layout()
		pl.savefig('figures/'+band+'band_comparison_with_DESY1'+naming+'.png',dpi=300)
		#pl.show()


if do_y1_histograms:
	emin, emax = -0.25,0.25
	tmin, tmax = 0.0, 2.0
	demin,demax =-0.1, 0.1
	dtmin,dtmax =-0.2, 0.2
	bins=100

	f=fits.open('/Users/secco/Documents/projects/SHEAR/shearcat/code/rowe_stats/DES_PSF_catalogs/psf_y1a1-v13.fits')
	for band in ['g', 'r', 'i', 'z']:		
		mask = np.where(f[1].data['filter']==band)
		g1s,g2s,Ts = f[1].data['e1'][mask], f[1].data['e2'][mask], f[1].data['size'][mask]
		g1m,g2m,Tm = f[1].data['psf_e1'][mask], f[1].data['psf_e2'][mask], f[1].data['psf_size'][mask]

		modelfail=np.isnan(g1m)+np.isnan(g2m)+np.isnan(Tm)
		starfail=np.isnan(g1s)+np.isnan(g2s)+np.isnan(Ts)
		mask = modelfail+starfail
		#~mask is "not mask"
		g1s = g1s[~mask]
		g2s = g2s[~mask]
		Ts = Ts[~mask]
		g1m = g1m[~mask]
		g2m = g2m[~mask]
		Tm = Tm[~mask]

		dT=Ts-Tm
		dg1 = g1s-g1m
		dg2 = g2s-g2m

		pl.figure(figsize=(12,9))
		pl.suptitle('DES Y1 %d stars (%s band)'%(len(dT),band))

		pl.subplot(331)
		pl.title(r'$e^{\mathrm{star}}_{1}$')
		pl.hist(g1s,bins=bins,range=[emin,emax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(332)
		pl.title(r'$e^{\mathrm{star}}_{2}$')
		pl.hist(g2s,bins=bins,range=[emin,emax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(333)
		pl.title(r'$T^{\mathrm{star}}$')
		pl.hist(Ts,bins=bins,range=[tmin,tmax])
		pl.yscale('log')
		#pl.axvline(x=0.0,ls='--',color='k')


		pl.subplot(334)
		pl.title(r'$e^{\mathrm{model}}_{1}$')
		pl.hist(g1m,bins=bins,range=[emin,emax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(335)
		pl.title(r'$e^{\mathrm{model}}_{2}$')
		pl.hist(g2m,bins=bins,range=[emin,emax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(336)
		pl.title(r'$T^{\mathrm{model}}$')
		pl.hist(Tm,bins=bins,range=[tmin,tmax])
		pl.yscale('log')
		#pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(337)
		pl.title(r'$\delta e_1$')
		pl.hist(dg1,bins=bins,range=[demin,demax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(338)
		pl.title(r'$\delta e_1$')
		pl.hist(dg2,bins=bins,range=[demin,demax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(339)
		pl.title(r'$\delta T$')
		pl.hist(dT,bins=bins,range=[dtmin,dtmax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')
		
		
		pl.tight_layout()
		pl.savefig('figures/'+band+'band_1Dhists_DESY1'+naming+'.png',dpi=300)
		pl.show()






if do_focalplane:
	emin, emax = -0.25,0.25
	tmin, tmax = 0.0, 2.0
	demin,demax =-0.1, 0.1
	dtmin,dtmax =-0.2, 0.2
	bins=100

	for band in ['g', 'r', 'i', 'z']:		
		if band=='g':
			f = fits.open('gband_500exps.fits.fz')
		if band=='r':
			f = fits.open('rband_500exps.fits.fz')
		if band=='i':
			f = fits.open('iband_500exps.fits.fz')
		if band=='z':
			f = fits.open('zband_500exps.fits.fz')

		g1s,g2s,Ts = f[1].data['g1_star'], f[1].data['g2_star'], f[1].data['T_star']
		g1m,g2m,Tm = f[1].data['g1_model'], f[1].data['g2_model'], f[1].data['T_model']

		modelfail=np.isnan(g1m)*np.isnan(g2m)*np.isnan(Tm)
		starfail=np.isnan(g1s)*np.isnan(g2s)*np.isnan(Ts)
		mask = modelfail+starfail
		#~mask is "not mask"
		g1s = g1s[~mask]
		g2s = g2s[~mask]
		Ts = Ts[~mask]
		g1m = g1m[~mask]
		g2m = g2m[~mask]
		Tm = Tm[~mask]

		dT=Ts-Tm
		dg1 = g1s-g1m
		dg2 = g2s-g2m

		pl.figure(figsize=(12,9))
		pl.suptitle('%d stars (%s band)'%(len(dT),band))

		pl.subplot(331)
		pl.title(r'$e^{\mathrm{star}}_{1}$')
		pl.hist(g1s,bins=bins,range=[emin,emax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(332)
		pl.title(r'$e^{\mathrm{star}}_{2}$')
		pl.hist(g2s,bins=bins,range=[emin,emax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(333)
		pl.title(r'$T^{\mathrm{star}}$')
		pl.hist(Ts,bins=bins,range=[tmin,tmax])
		pl.yscale('log')
		#pl.axvline(x=0.0,ls='--',color='k')


		pl.subplot(334)
		pl.title(r'$e^{\mathrm{model}}_{1}$')
		pl.hist(g1m,bins=bins,range=[emin,emax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(335)
		pl.title(r'$e^{\mathrm{model}}_{2}$')
		pl.hist(g2m,bins=bins,range=[emin,emax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(336)
		pl.title(r'$T^{\mathrm{model}}$')
		pl.hist(Tm,bins=bins,range=[tmin,tmax])
		pl.yscale('log')
		#pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(337)
		pl.title(r'$\delta e_1$')
		pl.hist(dg1,bins=bins,range=[demin,demax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(338)
		pl.title(r'$\delta e_1$')
		pl.hist(dg2,bins=bins,range=[demin,demax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')

		pl.subplot(339)
		pl.title(r'$\delta T$')
		pl.hist(dT,bins=bins,range=[dtmin,dtmax])
		pl.yscale('log')
		pl.axvline(x=0.0,ls='--',color='k')
		
		
		pl.tight_layout()
		pl.savefig('figures/'+band+'band_1Dhists'+naming+'.png',dpi=600)
		pl.show()


	
