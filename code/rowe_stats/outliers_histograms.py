import matplotlib.pyplot as pl 
import numpy as np 
from astropy.io import fits

emin, emax = -0.25,0.25
tmin, tmax = 0.0, 2.0
demin,demax =-0.1, 0.1
dtmin,dtmax =-0.2, 0.2
bins=100

for band in ['g', 'r', 'i', 'z']:		
	if band=='g':
		f = fits.open('gband_1995exps.fits.fz')
	if band=='r':
		f = fits.open('rband_2058exps.fits.fz')
	if band=='i':
		f = fits.open('iband_1992exps.fits.fz')
	if band=='z':
		f = fits.open('zband_1990exps.fits.fz')

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
	pl.hist(g1m,bins=bins,range=[demin,demax])
	pl.yscale('log')
	pl.axvline(x=0.0,ls='--',color='k')

	pl.subplot(338)
	pl.title(r'$\delta e_1$')
	pl.hist(g2m,bins=bins,range=[demin,demax])
	pl.yscale('log')
	pl.axvline(x=0.0,ls='--',color='k')

	pl.subplot(339)
	pl.title(r'$\delta T$')
	pl.hist(Tm,bins=bins,range=[dtmin,dtmax])
	pl.yscale('log')
	pl.axvline(x=0.0,ls='--',color='k')
	
	
	pl.tight_layout()
	pl.savefig(band+'band_1Dhists.png',dpi=600)
	pl.show()


	