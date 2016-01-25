import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import astropy
from astropy.io import fits
from photutils import CircularAperture
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5 
from astropy.table import Table, join
from photutils import SkyCircularAperture
from photutils import aperture_photometry
from astropy.table import Column


coords=["210.35953 +54.588268","210.35006 +54.588958","210.28666 +54.613781","210.48034 +54.557275", "210.53029 +54.607541", "210.31118 +54.48097",  "210.17332 +54.513975",  "210.68205 +54.63497"]  
aa = Column(['Nucleus','Hodge602','Searle5','NGC5461','NGC5462','NGC5455','NGC5447','NGC5471'], name='regions')
positions= SkyCoord(coords, FK4, unit=(u.deg))

apertures = SkyCircularAperture(positions, r=10.*u.arcsec)
newtab=Table.read('/Users/Andromeda/Desktop/project/data_mining/m101/fits_table/fuv.fits')

ims=glob.glob('/Users/Andromeda/Desktop/project/data_mining/m101/images/*.fits')

for im in ims:
	hdu = fits.open(im)
 	phot_table = aperture_photometry(hdu, apertures)
 	phot_table.add_column(aa, index=0)
 	phot_table.rename_column('aperture_sum',im)
 	del phot_table['xcenter','ycenter','center_input']
 	newtab=join(newtab,phot_table,keys='regions',join_type='outer')


newtab.rename_column('M101_FUV_Wm-2.fits_1', 'M101_FUV_Wm-2.fits')
del newtab['M101_FUV_Wm-2.fits_2']

newtab.write('/Users/Andromeda/Desktop/project/data_mining/m101/fits_table/m101_phot_data.fits')   
newtab.write('/Users/Andromeda/Desktop/project/data_mining/m101/fits_table/m101_phot_data.csv')   