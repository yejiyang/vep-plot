#!/bin/sh

# about beachball.txt  
# this file includes the data to plot beachball, calculated from Yang's method.
# the original file is from '/home/jiyang/opensource/vep-plot/app3stress/ori_results/app32', calculated by 'p12345.f'

m=4  
input=./input/sagrid.meter.83A6_4w.txt
inputybb=./input/beachball.meter.83A6_4w.txt
inputprinstr=./input/prin_stress.txt
#input=./ores.$1/sagrid.$1.txt
#prepare the data
#title=$1.beachball.vep
title=83A6_4w.$m.beachball.vep
psfile=${title}.ps

#range="-R-121.6/-113.2/31.4/38.1"
#proj="-JM18.39c"
#ticks="-B2f1eSWn"
range="-R-220/280/-350/300"
proj="-JX6i/8i"
#ticks="-B:."Only SAF, 15km Depth":50f25eSWn"
let depth=5\*$m-5
title1=Only_SAF_${depth}_km
ticks="-B50f25eSWn:."$title1":"
V="-V"
#V=""

para_path=./para/

psbasemap ${range} ${proj} $ticks -X3.0c -Y6c -P $V -K > ${psfile}
#plot the faults
psxy ./para/saf.km  $proj $range -W3,black -m -O -V -K>>$psfile
#psxy ${para_path}saf627bddegree.txt $proj $range -W4,grey,- -m -O $V -K>>$psfile


# beachball from hui
#awk '{if(($1=='$m')&& (int($2*10/5)==($2*10/5)) && int($3*10/5)==($3*10/5))  print $2,$3,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20}' $input > bb.txt
#psmeca bb.txt  -R -J -Sx0.09i -Ggrey -Ewhite  -O  -W2  -C -V -K>>$psfile

# beachball from yang
awk '{if(($3=='$m')&& int($1*10/400000)==($1*10/400000) && int($2*10/400000)==($2*10/400000))  print $1/1000,$2/1000,5*($3-1),$4,$5,$6,$7,$8,$9,$10,$11/1000,$12/1000}' $inputybb > yang.txt 
psmeca yang.txt -R  -J  -C -O -V  -Sd0.19i/1 -Ggrey50 -Ewhite -K>> $psfile

# about the direction of 2D principal stress
# our matlab output, extension--positive, compression--negative, but the sita equals to I1 counter-clockwise to x-axis.
# in gmt, psvelo, -Sx option, extension and compression definition the same. But the angle defition is different. 
# the angle in gmt, is azimuth of II2 in degrees CW (clockwise) from North.
# So we need transform the matlab output sita to GMT angle. I1=sita, I2=sita+90, angle=90-sita=90-sita-90=-sita, this's why -$9
awk '{if(($1=='$m')&& int(($2)*10/400)==(($2)*10/400) && int($3*10/400)==($3*10/400))  print $2,$3,$7,$8,-$9}' $inputprinstr > prin.txt
# plot extension, red color
awk '{if(($3>0))  print $1,$2,$3,0,$5}' prin.txt > prin1.txt
awk '{if(($4>0))  print $1,$2,0,$4,$5}' prin.txt > prin2.txt
psvelo prin1.txt -A0.01/0.6/0.1 -W0.25p,red  -R  -J   -O -V  -Sx0.012  -K>> $psfile
psvelo prin2.txt -A0.01/0.6/0.1 -W0.25p,red  -R  -J   -O -V  -Sx0.012  -K>> $psfile
# plot compression, blue color
awk '{if(($3<0))  print $1,$2,$3,0,$5}' prin.txt > prin3.txt
awk '{if(($4<0))  print $1,$2,0,$4,$5}' prin.txt > prin4.txt
psvelo prin3.txt -A0.01/0.6/0.1 -W0.25p,blue  -R  -J   -O -V  -Sx0.012  -K>> $psfile
psvelo prin4.txt -A0.01/0.6/0.1 -W0.25p,blue  -R  -J   -O -V  -Sx0.012  -K>> $psfile
# -A default -A0.03/0.12/0.09


#---------end of plot ---------
psxy -J -R -T -O >> $psfile
### clean the trash
#rm bb.txt prin.txt yang.txt prin0.txt

#gnome-open ${title}.ps
open ${title}.ps
#mv *.ps ./output
