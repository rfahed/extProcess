#!/usr/bin/env python

from Phot import catalog,image
import numpy as np
from astropy.io import fits

filter ="g"

imname="../Images/DES{}_sim.fits".format(filter)
incat="input_{}_matchtag.txt".format(filter)
outcat="output_{}_matchtag.txt".format(filter)
mergedcat="merged_{}.txt".format(filter)
realcat="real_{}.txt".format(filter)

inc=catalog.read(incat)
outc=catalog.read(outcat)
merc=catalog.read(mergedcat)
realc=catalog.read(realcat)

ext=5

im=fits.open(imname)
footprint = image.get_footprint(im[ext].header)

# filter bad detection (close to saturation or border)
#merc=merc[merc['FLAGS'] == 0]
#outc=outc[outc['FLAGS'] == 0]
#realc=realc[realc['FLAGS'] == 0]

minb=10
maxb=30
binw=0.2

bins=np.linspace(minb,maxb,int((maxb-minb)/binw))

catalog.P.subplot(2, 1, 1)
n_in,bins,patches_in=catalog.histogramcol(inc,"MAG",label="input catalog", bins=bins)
n_out,bins,patches_out=catalog.histogramcol(outc,"MAG_AUTO", label="output catalog", bins=bins, xlab='Magnitude')
catalog.P.legend(loc=0)

catalog.P.subplot(2, 1, 2)
catalog.P.step(bins[0:-1],n_out/n_in,label="detection fraction ({} band)".format(filter))

catalog.P.xlabel('Magnitude')
catalog.P.ylabel('detection fraction')
catalog.P.tight_layout()
catalog.P.savefig("inout_{}.png".format(filter))
#catalog.P.show()




