#!/usr/bin/env python

from Phot import catalog
import numpy as np

filters=['g'] #,'r','i','z','Y']
colors=['g'] #,'r'] #,'y','b','k']

for filter,c in zip(filters,colors) :

    incat="input_{}.txt".format(filter)
    outcat="output_{}.txt".format(filter)
    mergedcat="merged_{}.txt".format(filter)
    realcat="real_{}.txt".format(filter)

    #inc=catalog.read(incat,format="ext")
    #outc=catalog.read(outcat)
    merc=catalog.read(mergedcat)
    #realc=catalog.read(realcat)

    # filter bad detection (close to saturation or border)
    merc=merc[merc['FLAGS'] == 0]
    #outc=outc[outc['FLAGS'] == 0]
    #realc=realc[realc['FLAGS'] == 0]

    catalog.scattercols(merc,'MAG','MAG_AUTO',xlab='input mag',ylab='output mag',label=filter+'  band',color=c,alpha=0.4)
    catalog.P.plot([0,40],[0,40],'--')

catalog.P.legend()
catalog.P.xlim([10,30])
catalog.P.ylim([10,30])
catalog.P.tight_layout()
catalog.P.savefig("magmag_{}.png".format('-'.join(filters)))
catalog.P.show()
#func1 = lambda x,y : np.pi*(x**2)*y/4
#func2 = lambda x,y : np.pi*x*y

#mergals = merc[merc["OBJTYPE"] == 2]
#catalog.add_calc_column(mergals,'DLENGTH','DAXISRATIO',func1)
#catalog.add_calc_column(mergals,'A_IMAGE','B_IMAGE',func2)


