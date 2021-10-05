#!/bin/bash 

cd ${NCPAPROP_DIR}/samples

# example 1
../bin/ePape --singleprop --atmosfile NCPA_canonical_profile_trimmed.dat --freq 0.1 --azimuth 90 --maxrange_km 1000

mv tloss_1d.pe ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex1_tloss_1d.pe
mv tloss_2d.pe ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex1_tloss_2d.pe

# example 2
../bin/ePape --singleprop --atmosfile2d toy_profile_2d_summary.dat --freq 0.5 --azimuth 90 --write_2d_tloss --maxrange_km 1800

mv tloss_1d.pe ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex2_tloss_1d.pe
mv tloss_2d.pe ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex2_tloss_2d.pe

# example 3
../bin/ePape --multiprop --atmosfile NCPA_canonical_profile_trimmed.dat \
         --freq 0.5 --azimuth_start 0 --azimuth_end 360 --azimuth_step 2 \
         --maxrange_km 1000 

mv tloss_multiprop.pe ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex3_tloss_multiprop.pe

../bin/ePape --multiprop --atmosfile NCPA_canonical_profile_trimmed.dat \
         --freq 0.5 --azimuth_start 0 --azimuth_end 360 --azimuth_step 2 \
         --maxrange_km 1000 --lossless

mv tloss_multiprop.pe ${NCPAPROP_DIR}/manual/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex3_tloss_multiprop_lossless.pe

cd -
