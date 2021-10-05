#!/bin/bash 

cd ${NCPAPROP_DIR}/samples

# example 2

    ../bin/Modess --singleprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat \
         --azimuth 90 --freq 0.1 --write_2d_tloss
         
mv tloss_1d.nm tloss_1d.lossless.nm tloss_2d.nm ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/modess_gp_files/code_outputs/

# example 3

    ../bin/Modess --multiprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat \
         --freq 0.1 --azimuth_start 0 --azimuth_end 360 --azimuth_step 1
mv Nby2D_tloss_1d.nm Nby2D_tloss_1d.lossless.nm ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/modess_gp_files/code_outputs/


# example 4

    ../bin/Modess --singleprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat \
         --azimuth 90 --freq 0.1 --write_2d_tloss --sourceheight_km 20
                 
mv tloss_2d.nm ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/modess_gp_files/code_outputs/tloss_2d_20km.nm

cd -

