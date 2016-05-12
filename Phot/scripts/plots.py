#!/usr/bin/env python

from Phot import catalog

incat="input_g.txt"
outcat="output_g.txt"
mergedcat="merged_g.txt"

inc=catalog.read(incat,format="ext")
outc=catalog.read(outcat)
merc=catalog.read(mergedcat)


