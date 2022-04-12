import matplotlib.pyplot as pl 
import numpy as np 
from astropy.io import fits

BINS=150

band='g'

if band=='g':
	f = fits.open('gband_1995exps.fits.fz')
if band=='r':
	f = fits.open('rband_2058exps.fits.fz')
if band=='i':
	f = fits.open('iband_1992exps.fits.fz',memmap=False)
if band=='z':
	f = fits.open('zband_1990exps.fits.fz')

e1min,e1max = -0.05,0.05
e2min,e2max = -0.05,0.05
Tmin,Tmax = 0.2, 0.8
dg1min, dg1max = -0.005, 0.005
dg2min, dg2max = -0.005, 0.005
dTmin, dTmax = -0.02, 0.02


focal_x,focal_y = f[1].data['focal_x'],f[1].data['focal_y'] 
g1s,g2s,Ts = f[1].data['g1_star'], f[1].data['g2_star'], f[1].data['T_star']
g1m,g2m,Tm = f[1].data['g1_model'], f[1].data['g2_model'], f[1].data['T_model']

modelfail=np.isnan(g1m)*np.isnan(g2m)*np.isnan(Tm)
starfail=np.isnan(g1s)*np.isnan(g2s)*np.isnan(Ts)
mask = modelfail+starfail
#~mask is "not mask"
focal_x = focal_x[~mask]
focal_y = focal_y[~mask]
g1s = g1s[~mask]
g2s = g2s[~mask]
Ts = Ts[~mask]
g1m = g1m[~mask]
g2m = g2m[~mask]
Tm = Tm[~mask]


dT=Ts-Tm
dg1 = g1s-g1m
dg2 = g2s-g2m
counts = pl.hist2d(focal_x,focal_y,bins=BINS)[0]

hist_Ts = pl.hist2d(focal_x,focal_y,bins=BINS,weights=Ts)[0]/counts
hist_g1s = pl.hist2d(focal_x,focal_y,bins=BINS,weights=g1s)[0]/counts
hist_g2s = pl.hist2d(focal_x,focal_y,bins=BINS,weights=g2s)[0]/counts

hist_Tm = pl.hist2d(focal_x,focal_y,bins=BINS,weights=Tm)[0]/counts
hist_g1m = pl.hist2d(focal_x,focal_y,bins=BINS,weights=g1m)[0]/counts
hist_g2m = pl.hist2d(focal_x,focal_y,bins=BINS,weights=g2m)[0]/counts

hist_dT = pl.hist2d(focal_x,focal_y,bins=BINS,weights=dT)[0]/counts
hist_dg1 = pl.hist2d(focal_x,focal_y,bins=BINS,weights=dg1)[0]/counts
hist_dg2 = pl.hist2d(focal_x,focal_y,bins=BINS,weights=dg2)[0]/counts

pl.figure(figsize=(12,9))
pl.suptitle('%d stars (%s band)'%(len(dT),band))

pl.subplot(331)
pl.title(r'$e^{\mathrm{star}}_{1}$')
pl.imshow(hist_g1s,cmap='bwr',vmin=e1min,vmax=e1max)
pl.colorbar()

pl.subplot(332)
pl.title(r'$e^{\mathrm{star}}_{2}$')
pl.imshow(hist_g2s,cmap='bwr',vmin=e2min,vmax=e2max)
pl.colorbar()

pl.subplot(333)
pl.title(r'$T^{\mathrm{star}}$')
pl.imshow(hist_Ts,cmap='bwr',vmin=Tmin,vmax=Tmax)
pl.colorbar()

pl.subplot(334)
pl.title(r'$e^{\mathrm{model}}_{1}$')
pl.imshow(hist_g1m,cmap='bwr',vmin=e1min,vmax=e1max)
pl.colorbar()

pl.subplot(335)
pl.title(r'$e^{\mathrm{model}}_{2}$')
pl.imshow(hist_g2m,cmap='bwr',vmin=e2min,vmax=e2max)
pl.colorbar()

pl.subplot(336)
pl.title(r'$T^{\mathrm{model}}$')
pl.imshow(hist_Tm,cmap='bwr',vmin=Tmin,vmax=Tmax)
pl.colorbar()

pl.subplot(337)
pl.title(r'$\delta e_1$')
pl.imshow(hist_dg1,cmap='bwr',vmin=dg1min,vmax=dg1max)
pl.colorbar()

pl.subplot(338)
pl.title(r'$\delta e_2$')
pl.imshow(hist_dg2,cmap='bwr',vmin=dg2min,vmax=dg2max)
pl.colorbar()

pl.subplot(339)
pl.title(r'$\delta T$')
pl.imshow(hist_dT,cmap='bwr',vmin=dTmin,vmax=dTmax)
pl.colorbar()


pl.tight_layout()
pl.savefig(band+'_focalplane.png',dpi=600)
pl.show()