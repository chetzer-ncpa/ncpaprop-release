#!/bin/bash 

cd ${NCPAPROP_DIR}/samples

# example 1

    ../bin/WMod --singleprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat \
         --azimuth 90 --freq 0.1 --write_2d_tloss

mv wtloss_1d.nm wtloss_1d.lossless.nm wtloss_2d.nm ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/wmod_gp_files/code_outputs

# example 2

    ../bin/WMod --multiprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat \
         --freq 0.1 --azimuth_start 0 --azimuth_end 360 --azimuth_step 1

mv Nby2D_wtloss_1d.nm Nby2D_wtloss_1d.lossless.nm ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/wmod_gp_files/code_outputs

cd -

