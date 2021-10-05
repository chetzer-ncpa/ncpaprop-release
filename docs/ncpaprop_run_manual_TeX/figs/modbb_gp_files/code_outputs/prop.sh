#!/bin/bash

ModBB --pulse_prop_src2rcv myDispersionFile.dat --range_R_km 270 --waveform_out_file mywavf.dat --use_builtin_pulse --max_celerity 320

ModBB --pulse_prop_src2rcv_grid myDispersionFile.dat --R_start_km 220 --DR_km 20 --R_end_km 300 --waveform_out_file mywavf_grid.dat --use_builtin_pulse --max_celerity 320

