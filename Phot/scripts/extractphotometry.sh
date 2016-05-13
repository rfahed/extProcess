#!/usr/bin/env bash
run sextractor DESg_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_g.txt
run sextractor DESr_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_r.txt
run sextractor DESi_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_i.txt
run sextractor DESz_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_z.txt
run sextractor DESY_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_Y.txt

run sextractor DES_g.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_g.txt
run sextractor DES_r.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_r.txt
run sextractor DES_i.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_i.txt
run sextractor DES_z.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_z.txt
run sextractor DES_Y.fits --zeropoint INTEGRATED --zerokey MAGZERO --outputcat real_Y.txt




