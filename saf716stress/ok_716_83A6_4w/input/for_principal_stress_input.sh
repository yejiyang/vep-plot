#!/bin/bash

#	$1 layer number
#	$2 x or longi,	$3 y or lat
#	$7 u11, $8 u22, $12 u12, unit Pa, afterProcess unit MPa.
#	use printf not print for format output	
#awk '{ printf "%2s %6s %6s %13s %13s %13s\n",   $1, $2, $3, $7/1e6, $8/1e6, $12/1e6}'  ./sagrid.1_5w.txt > input_stress_for_prin_stress
awk '{ printf "%2s %6s %6s %13s %13s %13s\n",   $1, $2/1000, $3/1000, $7/1e6, $8/1e6, $12/1e6}'  ./sagrid.meter.83A6_4w.txt > input_stress_for_prin_stress
#awk '{ print   $1, $2/1000, $3/1000, $7/1e6, $8/1e6, $12/1e6}' OFS="," ./sagrid.31.txt > input_stress_for_prin_stress.csv
