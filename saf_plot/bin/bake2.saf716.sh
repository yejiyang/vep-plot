#!/bin/sh
gmtset FRAME_WIDTH 0.05c TICK_LENGTH -0.075i ANNOT_FONT_SIZE 12p PLOT_DEGREE_FORMAT ddd:mm:ss
#title=safpap5
title=saf716_5

input=cgrid716_5.txt
#prepare the data

psfile=${title}.ps #${psfile}

#range="-R-125/-115/32/40"
range="-R-121.6/-113.2/31.4/38.1"
proj="-JM8.39c"
ticks="-B2f1eSWn"

psbasemap ${range} ${proj} $ticks -X1.5c -Y19.0c -P -V -K > ${psfile}
#psbasemap ${range} ${proj} -B4f2 -X1.5c -Y19.0c -P -V -K > ${psfile}

awk -F ' ' '{print $2, $3,sqrt($4*$4+$5*$5)/5}' ${input} > mean.xyz
xyz2grd mean.xyz -Gdata.grd -I0.1 ${range} -V
grdsample data.grd -Gdata1.grd -I0.01 -Qb
#grd2cpt -V data1.grd -CGMT_strainrate.cpt -S0/50/0.5 -L0/50 -D  > velo.cpt
grdimage data1.grd -Cvect_hui.cpt ${range} ${proj} -K -V -P -O >>${psfile}
#plot the faults
psclip saf627bddegree.txt $range $proj -m -O -K >>$psfile
#psxy fault_data.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psxy saf1.txt.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psclip -C -O -K >>$psfile
psxy saf627bddegree.txt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile

awk -F ' ' '{if((int($2*10/5)==$2*10/5)&&(int($3*10/5)==$3*10/5)) print $2,$3,$4/15*0.6,$5/15*0.6,0,0,0}' ${input} |psvelo ${proj} ${range} -G100/100/100 -A0.025/0.150/0.05 -Se0.015/0.95/0 -W2/0/0/0 -K -O -V >>${psfile}

pstext <<END ${range} ${proj} -S0.5p -O -K -V >>${psfile}
#-121.8  34.6 11 0 0 0 Pacific  
#-120.8  33.9 11 0 0 0 Plate    
#-119.0  37.3 11 0 0 0 North    
#-118.5  36.6 11 0 0 0 American 
#-118.0  35.9 11 0 0 0 Plate    
-121.0  32.0 11 0 0 0 (a) 
END

echo  3.885  4.50 > tmp
echo  5.686  4.50 >> tmp
echo  5.686  5.502 >> tmp
echo  3.885  5.502 >> tmp
psxy tmp -R0/30/0/20 -Jx0.58i -G255 -K -L -O -W3/0/0/0 -V >> ${psfile}

pstext <<END -R -Jx -P -S0.5p -O -K -V>>${psfile}
4.09 5.2 8 0 0 0 Velocity (mm/yr) 
END

gmtset  ANNOT_FONT_SIZE 6p
psscale -D7.035c/7.3c/2.06c/0.18ch -Cvect_hui.cpt  -I0.3 -B10  -O -K -V>> ${psfile}
gmtset  ANNOT_FONT_SIZE 12p


#rem ------------------------------------------------------------------------------
gmtset  ANNOT_FONT_SIZE 9p

psbasemap $range  $proj $ticks -X9.5c -K -O >>$psfile

awk -F ' ' '{print $2,$3,$9*1e7}' $input |blockmean -R -I0.1 > mean.xyz
xyz2grd mean.xyz -Gdata.grd -I0.1 -R
grdsample data.grd -Gdata1.grd -I0.01 -Qb
grdimage data1.grd -Ccfs_hui.cpt $proj -R -K -V -P -O >>$psfile

#plot the faults
#psxy fault_data.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psclip saf627bddegree.txt $range $proj -m -O -K >>$psfile
psxy saf1.txt.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psclip -C -O -K >>$psfile
psxy saf627bddegree.txt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile

pstext <<END -R $proj -S0.5p -O -K >>$psfile
#-121.8  34.6 11 0 0 0 Pacific  
#-120.8  33.9 11 0 0 0 Plate    
#-119.0  37.3 11 0 0 0 North    
#-118.5  36.6 11 0 0 0 American 
#-118.0  35.9 11 0 0 0 Plate    
-121.0  32.0 11 0 0 0 (b) 
END

echo  3.885  4.50 > tmp
echo  5.686  4.50 >> tmp
echo  5.686  5.502 >> tmp
echo  3.885  5.502 >> tmp
psxy tmp -R0/30/0/20 -Jx0.58i -G255 -K -L -O -W3/0/0/0 -V >> ${psfile}

pstext <<END -R -Jx -P -S0.5p -O -K -V>>${psfile}
3.95 5.2 8 0 0 0 Strain Rate (10@+-7@+/yr) 
END

gmtset  ANNOT_FONT_SIZE 6p
psscale -D7.035c/7.3c/2.06c/0.18ch -Ccfs_hui.cpt  -S -B10  -O -K -V>> ${psfile}
#psscale -D4.235c/4.550c/1.96c/0.14ch -Ccfs_hui.cpt -B10 -S -O -K>> $psfile

#rem ------------------------------------------------------------------------------
gmtset  ANNOT_FONT_SIZE 9p

psbasemap $range  $proj $ticks -X-9.5c -Y-9.5c  -K -O >>$psfile

awk -F ' ' '{print $2,$3,$7/1e6}' $input |blockmean -R -I0.1 > mean.xyz
xyz2grd mean.xyz -Gdata.grd -I0.1 -R
grdsample data.grd -Gdata1.grd -I0.01 -Qb
grdimage data1.grd -Cshear_hui.cpt $proj -R -K -V -P -O >>$psfile

#plot the faults
#psxy fault_data.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psclip saf627bddegree.txt $range $proj -m -O -K >>$psfile
psxy saf1.txt.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psclip -C -O -K >>$psfile
psxy saf627bddegree.txt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile

pstext <<END -R $proj -S0.5p -O -K >>$psfile
#-121.8  34.6 11 0 0 0 Pacific  
#-120.8  33.9 11 0 0 0 Plate    
#-119.0  37.3 11 0 0 0 North    
#-118.5  36.6 11 0 0 0 American 
#-118.0  35.9 11 0 0 0 Plate    
-121.0  32.0 11 0 0 0 (c) 
END

echo  3.885  4.50 > tmp
echo  5.686  4.50 >> tmp
echo  5.686  5.502 >> tmp
echo  3.885  5.502 >> tmp
psxy tmp -R0/30/0/20 -Jx0.58i -G255 -K -L -O -W3/0/0/0 -V >> ${psfile}

pstext <<END -R -Jx -P -S0.5p -O -K -V>>${psfile}
4.09 5.3 8 0 0 0 Maximum Shear 
4.09 5.1 8 0 0 0 Stress (Mpa)
END

gmtset  ANNOT_FONT_SIZE 6p
psscale -D7.035c/7.3c/2.06c/0.18ch -Cshear_hui.cpt   -B15 -S  -O -K -V>> ${psfile}
#psscale -D4.235c/4.550c/1.96c/0.14ch -Cshear_hui.cpt -B15 -S -O -K >> $psfile

#rem ------------------------------------------------------------------------------
gmtset  ANNOT_FONT_SIZE 9p

psbasemap $range  $proj $ticks -X9.5c -K -O >>$psfile

#plot topo
#gmtset COLOF_MODEL hsv
#makecpt -Cno_green -Qo -T1/1000/1 -D -Z > strain.cpt
grd2cpt topo.grd -Crelief -Z > ca.cpt
grdgradient topo.grd -A0 -M -Nt -Gtopo_gradient.grd
#makecpt -Cgray -T-1.4/.7/.1 -Z > grad.cpt
grdimage topo.grd -Itopo_gradient.grd  $proj $range -Cca.cpt  -P -O -K -V >>$psfile
grdcontour topo.grd -C500 -R -J -O -V -K >> $psfile
psscale -Cca.cpt -D4/-1.2/6/0.5h -B2000/:meters: -O -K >>$psfile
psxy fault_data.gmt -m $proj $range -W2 -O -V -K >> $psfile

#plot the faults
#psxy fault_data.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psclip saf627bddegree.txt $range $proj -m -O -K >>$psfile
psxy saf1.txt.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psclip -C -O -K >>$psfile
psxy saf627bddegree.txt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile

pstext <<END -R $proj -S0.5p -O -K >>$psfile
-121.0  32.0 11 0 0 0 (d) 
END


###########
#rm *.grd mean.xyz

ps2pdf ${psfile} ${title}.pdf
ps2raster $psfile -A -P -Tg

#gs ${psfile}
#rcp ${psfile} geo01:/home/jiyang/Dropbox/temp
tar -czvf ${title}.fig.tgz ${psfile} ${title}.pdf ${title}.png
rcp ${title}.fig.tgz geo01:/home/jiyang/Dropbox/SAF/res

#rm ${psfile} ${title}.pdf ${title}.png
