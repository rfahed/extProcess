#!/usr/bin/env python

from Phot import catalog
import numpy as np

filter ="g"

incat="input_{}_matchtag.txt".format(filter)
outcat="output_{}_matchtag.txt".format(filter)
mergedcat="merged_{}.txt".format(filter)
realcat="real_{}.txt".format(filter)

inc=catalog.read(incat)
outc=catalog.read(outcat)
merc=catalog.read(mergedcat)
realc=catalog.read(realcat)

# filter bad detection (close to saturation or border)
#merc=merc[merc['FLAGS'] == 0]
#outc=outc[outc['FLAGS'] == 0]
#realc=realc[realc['FLAGS'] == 0]

minb=10
maxb=28
binw=0.2

bins=np.linspace(minb,maxb,int((maxb-minb)/binw))

catalog.P.subplot(2, 1, 1)
n_in,bins,patches_in=catalog.histogramcol(realc,"MAG_AUTO",label="real image", bins=bins)
n_out,bins,patches_out=catalog.histogramcol(outc,"MAG_AUTO", label="sim image", bins=bins, xlab='Magnitude')

catalog.P.subplot(2, 1, 2)
catalog.P.step(bins[0:-1],n_out-n_in,label="sim-real difference ({} band)".format(filter))
catalog.P.xlabel('Magnitude')
catalog.P.ylabel('sim-real difference')
catalog.P.tight_layout()
catalog.P.savefig("realout_{}.png".format(filter))
catalog.P.show()




