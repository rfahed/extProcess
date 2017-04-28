#!/usr/bin/env python

from Phot import catalog
import numpy as np
import scipy as sp
import sys

inmag='MAG'
outmag='MAG_BEST'

def magmag(cat,filter="X"):
    xlim=np.array([12,25])
    ylim=np.array([-2,2])
    nbins=xlim[1]-xlim[0]
    cat=catalog.read(cat,format="fits")

    # filter bad detection (close to saturation or border)
    cat=cat[cat['FLAGS'] == 0]
    catalog.add_diff_column(cat,outmag,inmag)
    catalog.P.subplot(211)
    catalog.scattercols(cat,inmag,'diff_{}_{}'.format(outmag,inmag),xlab='input mag',ylab='output mag - input mag',label=filter+'  band',alpha=0.4)
    catalog.P.plot([0,40],[0,0],'--',color='r')

    catalog.P.legend()
    catalog.P.xlim(xlim)
    catalog.P.ylim(ylim)


    catalog.P.subplot(212)
    f = lambda x : np.sqrt(np.mean(x**2))
    mean_error, bins = catalog.histogramcol(cat,inmag,'diff_{}_{}'.format(outmag,inmag), statistic=f, bins=nbins,range=[xlim])
    catalog.P.xlim(xlim)
    catalog.P.ylim(ylim+1)

    catalog.P.tight_layout()
    catalog.P.savefig("magmag_{}.png".format(filter))
    catalog.P.show()
    #func1 = lambda x,y : np.pi*(x**2)*y/4
    #func2 = lambda x,y : np.pi*x*y

    #mergals = merc[merc["OBJTYPE"] == 2]
    #catalog.add_calc_column(mergals,'DLENGTH','DAXISRATIO',func1)
    #catalog.add_calc_column(mergals,'A_IMAGE','B_IMAGE',func2)

if __name__=="__main__":
    magmag(sys.argv[1])
