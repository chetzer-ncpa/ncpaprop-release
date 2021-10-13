#!/bin/bash 

cd ${NCPAPROP_DIR}/samples

# example 1
../bin/ePape --singleprop --atmosfile NCPA_canonical_profile_trimmed.dat --freq 0.1 --azimuth 90 --maxrange_km 1000 --starter self --write_2d_tloss

mv tloss_1d.pe ${NCPAPROP_DIR}/docs/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex1_tloss_1d.pe
mv tloss_2d.pe ${NCPAPROP_DIR}/docs/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex1_tloss_2d.pe

# example 2
../bin/ePape --singleprop --atmosfile2d toy_profile_2d_summary.dat --freq 0.5 --azimuth 90 --write_2d_tloss --maxrange_km 1800 --starter self

mv tloss_1d.pe ${NCPAPROP_DIR}/docs/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex2_tloss_1d.pe
mv tloss_2d.pe ${NCPAPROP_DIR}/docs/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex2_tloss_2d.pe

# example 3
../bin/ePape --multiprop --atmosfile NCPA_canonical_profile_trimmed.dat \
         --freq 0.5 --azimuth_start 0 --azimuth_end 360 --azimuth_step 2 \
         --maxrange_km 1000 --starter self

mv tloss_multiprop.pe ${NCPAPROP_DIR}/docs/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex3_tloss_multiprop.pe

../bin/ePape --multiprop --atmosfile NCPA_canonical_profile_trimmed.dat \
         --freq 0.5 --azimuth_start 0 --azimuth_end 360 --azimuth_step 2 \
         --maxrange_km 1000 --lossless --starter self

mv tloss_multiprop.pe ${NCPAPROP_DIR}/docs/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex3_tloss_multiprop_lossless.pe


# example 4
../bin/ePape --singleprop --topo --starter self --atmosfile2d \
         toy_profile_2d_summary_topo.dat --freq 0.5 --azimuth 90 --maxrange_km \
         1200 --write_2d_tloss --write_topography

mv tloss_2d.pe ${NCPAPROP_DIR}/docs/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex4_tloss_2d_topo.pe
mv topography.pe ${NCPAPROP_DIR}/docs/ncpaprop_run_manual_TeX/figs/ePape_gp_files/code_outputs/ex4_topography.pe

cd -
