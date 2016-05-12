#!/usr/bin/env bash
run sextractor DESg_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_g.txt
run sextractor DESr_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_r.txt
run sextractor DESi_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_i.txt
run sextractor DESz_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_z.txt
run sextractor DESY_sim.fits --zeropoint PER_S --zerokey SIMMAGZP --outputcat output_Y.txt



