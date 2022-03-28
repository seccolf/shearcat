import matplotlib.pyplot as pl 
import numpy as np 

BINS=40

focal_x,focal_y, g1s,g2s,Ts,g1m,g2m,Tm = np.loadtxt('measurements/D00145973_r.txt',unpack=True,usecols=(0,1,6,7,8,9,10,11))

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
pl.suptitle('%d stars in expnum 145973 (r band)'%len(dT))
pl.subplot(331)
pl.title(r'$e^{\mathrm{star}}_{1}$')
pl.imshow(hist_g1s,cmap='bwr',vmin=-0.1,vmax=0.1)
pl.colorbar()
pl.subplot(332)
pl.title(r'$e^{\mathrm{star}}_{2}$')
pl.imshow(hist_g2s,cmap='bwr',vmin=-0.1,vmax=0.1)
pl.colorbar()
pl.subplot(333)
pl.title(r'$T^{\mathrm{star}}$')
pl.imshow(hist_Ts,cmap='bwr',vmin=0.25,vmax=0.45)
pl.colorbar()

pl.subplot(334)
pl.title(r'$e^{\mathrm{model}}_{1}$')
pl.imshow(hist_g1m,cmap='bwr',vmin=-0.1,vmax=0.1)
pl.colorbar()

pl.subplot(335)
pl.title(r'$e^{\mathrm{model}}_{2}$')
pl.imshow(hist_g2m,cmap='bwr',vmin=-0.1,vmax=0.1)
pl.colorbar()

pl.subplot(336)
pl.title(r'$T^{\mathrm{model}}$')
pl.imshow(hist_Tm,cmap='bwr',vmin=0.25,vmax=0.45)
pl.colorbar()

pl.subplot(337)
pl.title(r'$\delta e_1$')
pl.imshow(hist_dg1,cmap='bwr',vmin=-0.08,vmax=0.08)
pl.colorbar()
pl.subplot(338)
pl.title(r'$\delta e_2$')
pl.imshow(hist_dg2,cmap='bwr',vmin=-0.06,vmax=0.06)
pl.colorbar()
pl.subplot(339)
pl.title(r'$\delta T$')
pl.imshow(hist_dT,cmap='bwr',vmin=-0.06,vmax=0.06)
pl.colorbar()


pl.tight_layout()
pl.savefig('exp145973_rband_focalplane.png',dpi=600)
pl.show()