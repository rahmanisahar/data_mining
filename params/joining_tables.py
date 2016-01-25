import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
from astropy.table import Table, join
import matplotlib
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits



tab_names=glob.glob('~/Desktop/project/data_mining/m31/fits_tables/phot_data/*.fits')

newtab=Table.read('~/Desktop/project/data_mining/m31/fits_tables/phot_data/irac2.fits')

for tab in tab_names:
    jtab=Table.read(tab)
    newtab=join(newtab,jtab,keys=('PUB_ID','RADEG','DECDEG'),join_type='outer')
    
    

newtab.write('~/Desktop/project/data_mining/m31/fits_tables/m31_phot_data.fits')
