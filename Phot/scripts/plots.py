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

catalog.add_frac_column(outc,'FLUX_MAX','THRESHOLD')
catalog.scattercols(outc,'MAG_AUTO','frac_FLUX_MAX_THRESHOLD',color='g',label="sim data",alpha=0.4)

catalog.add_frac_column(realc,'FLUX_MAX','THRESHOLD')
catalog.scattercols(realc,'MAG_AUTO','frac_FLUX_MAX_THRESHOLD',color='r',label="real data",alpha=0.4)
catalog.p.show() 
catalog.p.axhline(y=3, linewidth=2, color = 'k')

func1 = lambda x,y : np.pi*(x**2)*y/4
func2 = lambda x,y : np.pi*x*y

mergals = merc[merc["OBJTYPE"] == 2]
catalog.add_calc_column(mergals,'DLENGTH','DAXISRATIO',func1)
catalog.add_calc_column(mergals,'A_IMAGE','B_IMAGE',func2)


