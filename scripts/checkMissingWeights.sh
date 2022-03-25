#!/bin/bash

dir=$1
config=$2
zen=$3

if test ! -f ${dir}/${config}_${zen}degzenith_0.5degoffset_2.0to5.0TeV_BDT.weights.xml; then
  echo "Missing file ${dir}/${config}_${zen}degzenith_0.5degoffset_2.0to5.0TeV_BDT.weights.xml, copying from upper energy band"
  cp ${dir}/${config}_${zen}degzenith_0.5degoffset_5.0to100.0TeV_BDT.weights.xml ${dir}/${config}_${zen}degzenith_0.5degoffset_2.0to5.0TeV_BDT.weights.xml
fi

if test ! -f ${dir}/${config}_${zen}degzenith_0.5degoffset_1.0to2.0TeV_BDT.weights.xml; then
  echo "Missing file ${dir}/${config}_${zen}degzenith_0.5degoffset_1.0to2.0TeV_BDT.weights.xml, copying from upper energy band"
  cp ${dir}/${config}_${zen}degzenith_0.5degoffset_2.0to5.0TeV_BDT.weights.xml ${dir}/${config}_${zen}degzenith_0.5degoffset_1.0to2.0TeV_BDT.weights.xml
fi

if test ! -f ${dir}/${config}_${zen}degzenith_0.5degoffset_0.5to1.0TeV_BDT.weights.xml; then
  echo "Missing file ${dir}/${config}_${zen}degzenith_0.5degoffset_0.5to1.0TeV_BDT.weights.xml, copying from upper energy band"
  cp ${dir}/${config}_${zen}degzenith_0.5degoffset_1.0to2.0TeV_BDT.weights.xml ${dir}/${config}_${zen}degzenith_0.5degoffset_0.5to1.0TeV_BDT.weights.xml
fi
    
if test ! -f ${dir}/${config}_${zen}degzenith_0.5degoffset_0.3to0.5TeV_BDT.weights.xml; then
  echo "Missing file ${dir}/${config}_${zen}degzenith_0.5degoffset_0.3to0.5TeV_BDT.weights.xml, copying from upper energy band"
  cp ${dir}/${config}_${zen}degzenith_0.5degoffset_0.5to1.0TeV_BDT.weights.xml ${dir}/${config}_${zen}degzenith_0.5degoffset_0.3to0.5TeV_BDT.weights.xml
fi
       
if test ! -f ${dir}/${config}_${zen}degzenith_0.5degoffset_0.1to0.3TeV_BDT.weights.xml; then
  echo "Missing file ${dir}/${config}_${zen}degzenith_0.5degoffset_0.1to0.3TeV_BDT.weights.xml, copying from upper energy band"
  cp ${dir}/${config}_${zen}degzenith_0.5degoffset_0.3to0.5TeV_BDT.weights.xml ${dir}/${config}_${zen}degzenith_0.5degoffset_0.1to0.3TeV_BDT.weights.xml
fi

if test ! -f ${dir}/${config}_${zen}degzenith_0.5degoffset_0.05to0.1TeV_BDT.weights.xml; then
  echo "Missing file ${dir}/${config}_${zen}degzenith_0.5degoffset_0.05to0.1TeV_BDT.weights.xml, copying from upper energy band"
  cp ${dir}/${config}_${zen}degzenith_0.5degoffset_0.1to0.3TeV_BDT.weights.xml ${dir}/${config}_${zen}degzenith_0.5degoffset_0.05to0.1TeV_BDT.weights.xml
fi
