#!/usr/bin/env python

from Phot import catalog
import numpy as np

filter="g"

incat="input_{}.txt".format(filter)
outcat="output_{}.txt".format(filter)
mergedcat="merged_{}.txt".format(filter)
realcat="real_{}.txt".format(filter)

inc=catalog.read(incat,format="ext")
outc=catalog.read(outcat)
merc=catalog.read(mergedcat)
realc=catalog.read(realcat)

# filter bad detection (close to saturation or border)
merc=merc[merc['FLAGS'] < 4]
outc=outc[outc['FLAGS'] < 4]
realc=realc[realc['FLAGS'] < 4]

catalog.add_frac_column(outc,'FLUX_MAX','THRESHOLD')
catalog.scattercols(outc,'MAG_AUTO','frac_FLUX_MAX_THRESHOLD',xlab=filter +' mag',ylab='signal to noise ratio',color='g',label="sim data",alpha=0.4)

catalog.add_frac_column(realc,'FLUX_MAX','THRESHOLD')
catalog.scattercols(realc,'MAG_AUTO','frac_FLUX_MAX_THRESHOLD',xlab=filter+' mag',ylab='signal to noise ratio',color='r',label="real data",alpha=0.4)
catalog.p.legend() 
catalog.p.axhline(y=3, linewidth=2, color = 'k')
catalog.p.show()

catalog.scattercols(merc,'MAG','MAG_AUTO',xlab='input mag',ylab='output mag',label=filter+'  band',alpha=0.4)

#func1 = lambda x,y : np.pi*(x**2)*y/4
#func2 = lambda x,y : np.pi*x*y

#mergals = merc[merc["OBJTYPE"] == 2]
#catalog.add_calc_column(mergals,'DLENGTH','DAXISRATIO',func1)
#catalog.add_calc_column(mergals,'A_IMAGE','B_IMAGE',func2)


