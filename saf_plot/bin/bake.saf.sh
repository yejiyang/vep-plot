#!/bin/sh
gmtset FRAME_WIDTH 0.05c TICK_LENGTH -0.075i ANNOT_FONT_SIZE 12p PLOT_DEGREE_FORMAT ddd:mm:ss
title=saftest

input=grid.txt
#prepare the data

psfile=${title}.ps #${psfile}

range="-R-125/-115/32/40"
proj="-JM5.39c"

psbasemap ${range} ${proj} -B4f2eSWn -X1.5c -Y19.0c -P -V -K > ${psfile}
#psbasemap ${range} ${proj} -B4f2 -X1.5c -Y19.0c -P -V -K > ${psfile}

awk -F ' ' '{print $2, $3,sqrt($4*$4+$5*$5)/5}' ${input} > mean.xyz
xyz2grd mean.xyz -Gdata.grd -I0.1 ${range} -V
grdsample data.grd -Gdata1.grd -I0.01 -Qb
#grd2cpt -V data1.grd -CGMT_strainrate.cpt -S0/50/0.5 -L0/50 -D  > velo.cpt
grdimage data1.grd -Cvect_hui.cpt ${range} ${proj} -K -V -P -O >>${psfile}
echo -124.70000     40.00000  > tmp
echo -120.00000     35.40000  >> tmp
echo -118.50000     36.20000  >> tmp
echo -114.90000     32.10000  >>tmp
psxy tmp ${proj} ${range} -W7/0/0/0 -M  -P -O -K >>${psfile}

awk -F ' ' '{if((int($2*10/5)==$2*10/5)&&(int($3*10/5)==$3*10/5)) print $2,$3,$4/15*0.6,$5/15*0.6,0,0,0}' ${input} |psvelo ${proj} ${range} -G100/100/100 -A0.025/0.150/0.05 -Se0.015/0.95/0 -W2/0/0/0 -K -O -V >>${psfile}

pstext <<END ${range} ${proj} -S0.5p -O -K -V >>${psfile}
-121.8  34.6 11 0 0 0 Pacific  
-120.8  33.9 11 0 0 0 Plate    
-119.0  37.3 11 0 0 0 North    
-118.5  36.6 11 0 0 0 American 
-118.0  35.9 11 0 0 0 Plate    
-124.5  32.5 11 0 0 0 (a) 
END

echo  3.185  3.995 > tmp
echo  5.586  3.995 >> tmp
echo  5.586  5.502 >> tmp
echo  3.185  5.502 >> tmp
psxy tmp -R0/30/0/20 -Jx0.38i -G255 -K -L -O -W3/0/0/0 -V >> ${psfile}

pstext <<END -R -Jx -P -S0.5p -O -K -V>>${psfile}
3.35 4.94 8 0 0 0 Velocity (mm/yr) 
END

gmtset  ANNOT_FONT_SIZE 6p
psscale -D4.235c/4.550c/1.96c/0.14ch -Cvect_hui.cpt  -I0.3 -B10  -O -K -V>> ${psfile}
gmtset  ANNOT_FONT_SIZE 12p


#rem ------------------------------------------------------------------------------
gmtset  ANNOT_FONT_SIZE 9p

psbasemap $range  $proj -B4f2eSWn -X6.5c -K -O >>$psfile

awk -F ' ' '{print $2,$3,$9*1e7}' $input |blockmean -R -I0.1 > mean.xyz
xyz2grd mean.xyz -Gdata.grd -I0.1 -R
grdsample data.grd -Gdata1.grd -I0.01 -Qb
grdimage data1.grd -Ccfs_hui.cpt $proj -R -K -V -P -O >>$psfile

echo -124.70000     40.00000  > tmp
echo -120.00000     35.40000  >> tmp
echo -118.50000     36.20000  >> tmp
echo -114.90000     32.10000  >>tmp
psxy tmp $proj -R -W7/0/0/0 -M  -P -O -K >>$psfile

pstext <<END -R $proj -S0.5p -O -K >>$psfile
-121.8  34.6 11 0 0 0 Pacific  
-120.8  33.9 11 0 0 0 Plate    
-119.0  37.3 11 0 0 0 North    
-118.5  36.6 11 0 0 0 American 
-118.0  35.9 11 0 0 0 Plate    
-124.5  32.5 11 0 0 0 (b)      
END

echo  3.185  4.095 > tmp
echo  5.586  4.095 >> tmp
echo  5.586  5.502 >> tmp
echo  3.185  5.502 >> tmp
psxy tmp -R0/30/0/20 -Jx0.38i -G255 -K -L -O -W3/0/0/0  >> $psfile

pstext <<END -R -Jx -P -S0.5p -O -K >>$psfile
3.30 4.94  7 0 4 0 Strain Rate (10@+-7@+/yr)
END

gmtset  ANNOT_FONT_SIZE 6p
psscale -D4.235c/4.550c/1.96c/0.14ch -Ccfs_hui.cpt -B10 -S -O -K>> $psfile

#rem ------------------------------------------------------------------------------
gmtset  ANNOT_FONT_SIZE 9p

psbasemap $range  $proj -B4f2eSWn -X6.5c  -K -O >>$psfile

awk -F ' ' '{print $2,$3,$7/1e6}' $input |blockmean -R -I0.1 > mean.xyz
xyz2grd mean.xyz -Gdata.grd -I0.1 -R
grdsample data.grd -Gdata1.grd -I0.01 -Qb
grdimage data1.grd -Cshear_hui.cpt $proj -R -K -V -P -O >>$psfile

echo -124.70000     40.00000  > tmp
echo -120.00000     35.40000  >> tmp
echo -118.50000     36.20000  >> tmp
echo -114.90000     32.10000  >>tmp
psxy tmp $proj -R -W7/0/0/0 -M  -P -O -K >>$psfile

pstext <<END -R $proj -S0.5p -O -K >>$psfile
-121.8  34.6 11 0 0 0 Pacific  
-120.8  33.9 11 0 0 0 Plate    
-119.0  37.3 11 0 0 0 North    
-118.5  36.6 11 0 0 0 American 
-118.0  35.9 11 0 0 0 Plate    
-124.5  32.5 11 0 0 0 (c)      
END

echo  3.185  4.095 > tmp
echo  5.586  4.095 >> tmp
echo  5.586  5.502 >> tmp
echo  3.185  5.502 >> tmp
psxy tmp -R0/30/0/20 -Jx0.38i -G255 -K -L -O -W3/0/0/0  >> $psfile

pstext <<END -R -Jx -P -S0.5p -O -K >>$psfile
3.465 5.175 7 0 0 0 Maximum Shear 
3.675 4.895 7 0 0 0 Stress (MPa)  
END

gmtset  ANNOT_FONT_SIZE 6p
psscale -D4.235c/4.550c/1.96c/0.14ch -Cshear_hui.cpt -B15 -S -O -K >> $psfile

rm *.grd mean.xyz

ps2pdf ${psfile} ${title}.pdf
ps2raster $psfile -A -P -Tg

gs ${psfile}
#rcp ${psfile} geo01:/home/jiyang/Dropbox/temp

