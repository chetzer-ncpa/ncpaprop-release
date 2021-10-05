#!/bin/bash

if [ ! -d "files" ]
then
    mkdir files
fi

if [ "$1" = "bib" ]
then
    pdflatex -output-directory files NCPA_prop_manual.tex
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    bibtex files/NCPA_prop_manual
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    pdflatex -output-directory files NCPA_prop_manual.tex
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    pdflatex -output-directory files NCPA_prop_manual.tex
else
    pdflatex -output-directory files NCPA_prop_manual.tex
fi

if [ -e "files/NCPA_prop_manual.pdf" ] 
then
    mv files/NCPA_prop_manual.pdf ..
fi

