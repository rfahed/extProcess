#!/usr/bin/env bash
run sextractor DESg_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_g.cat
#run sextractor DESr_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_r.cat
#run sextractor DESi_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_i.cat
#run sextractor DESz_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_z.cat
#run sextractor DESY_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_Y.cat

run sextractor DES_g.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_g.cat
#run sextractor DES_r.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_r.cat
#run sextractor DES_i.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_i.cat
#run sextractor DES_z.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_z.cat
#run sextractor DES_Y.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_Y.cat




