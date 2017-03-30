#!/usr/bin/env python

from Phot import catalog
import numpy as np
import sys

filters=['g'] #,'r','i','z','Y']
colors=['g'] #,'r'] #,'y','b','k']

def magmag(cat,filter="X"):
    cat=catalog.read(cat,format="fits")
        
    # filter bad detection (close to saturation or border)
    cat=cat[cat['FLAGS'] == 0]
    
    catalog.scattercols(cat,'MAG','MAG_AUTO',xlab='input mag',ylab='output mag',label=filter+'  band',alpha=0.4)
    catalog.P.plot([0,40],[0,40],'--')
    
    catalog.P.legend()
    catalog.P.xlim([10,30])
    catalog.P.ylim([10,30])
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
