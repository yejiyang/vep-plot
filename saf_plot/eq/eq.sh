#!/bin/sh

awk -F ' ' '{print $7, $6,$4}' sc_eq_catalog_M5-6_1932_2012.txt > eq5.dat
awk -F ' ' '{print $7, $6,$4}' sc_eq_catalog_M6-7_1932_2012.txt > eq6.dat
awk -F ' ' '{print $7, $6,$4}' sc_eq_catalog_M7-8_1932_2012.txt > eq7.dat
cat eq5.dat eq6.dat eq7.dat > eqall.dat


awk -F ' ' '{print $1, $2, 10,0,4,0,'M5'}' eq5.dat > eq5.gmt
awk -F ' ' '{print $1, $2, 10,0,4,0,'M6'}' eq6.dat > eq6.gmt
awk -F ' ' '{print $1, $2, 10,0,4,0,'M7'}' eq7.dat > eq7.gmt
