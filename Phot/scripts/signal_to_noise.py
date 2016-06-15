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
merc=merc[merc['FLAGS'] == 0]
outc=outc[outc['FLAGS'] == 0]
realc=realc[realc['FLAGS'] == 0]
xlim=[18,25]
ylim=[0,100]


catalog.add_frac_column(outc,'FLUX_MAX','THRESHOLD')
catalog.add_frac_column(realc,'FLUX_MAX','THRESHOLD')

catalog.P.subplot(2, 1, 1)
catalog.scattercols(outc,'MAG_AUTO','frac_FLUX_MAX_THRESHOLD',xlab=filter+' mag',ylab='signal to noise ratio',color='g',label="sim data",alpha=0.4)
catalog.P.axhline(y=3, linewidth=2, color = 'k')
catalog.P.xlim(xlim)
catalog.P.ylim(ylim)
catalog.P.legend() 

catalog.P.subplot(2, 1, 2)
catalog.scattercols(realc,'MAG_AUTO','frac_FLUX_MAX_THRESHOLD',xlab=filter+' mag',ylab='signal to noise ratio',color='r',label="real data",alpha=0.4)
catalog.P.axhline(y=3, linewidth=2, color = 'k')
catalog.P.xlim(xlim)
catalog.P.ylim(ylim)
catalog.P.legend() 

catalog.P.tight_layout()
catalog.P.savefig("s2n_{}.png".format(filter))
catalog.P.show()


#func1 = lambda x,y : np.pi*(x**2)*y/4
#func2 = lambda x,y : np.pi*x*y

#mergals = merc[merc["OBJTYPE"] == 2]
#catalog.add_calc_column(mergals,'DLENGTH','DAXISRATIO',func1)
#catalog.add_calc_column(mergals,'A_IMAGE','B_IMAGE',func2)


