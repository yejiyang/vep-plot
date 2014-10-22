#!/bin/sh

# about beachball.txt  
# this file includes the data to plot beachball, calculated from Yang's method.
# the original file is from '/home/jiyang/opensource/vep-plot/app3stress/ori_results/app32', calculated by 'p12345.f'

m=1  
input=./input/sagrid.83A6_2w.txt
inputybb=./input/beachball.83A6_2w.txt
inputprinstr=./input/prin_stress.txt
#input=./ores.$1/sagrid.$1.txt
#prepare the data
#title=$1.beachball.vep
title=83A6_2w.$m.beachball.vep
psfile=${title}.ps

range="-R-121.6/-113.2/31.4/38.1"
proj="-JM18.39c"
ticks="-B2f1eSWn"
V="-V"
#V=""

para_path=./para/

psbasemap ${range} ${proj} $ticks -X1.5c -Y2c -P $V -K > ${psfile}
#plot the faults
psxy ./para/716_83A6_3w.fault.gmt  $proj $range -W3,black -m -O -V -K>>$psfile
psxy ${para_path}saf627bddegree.txt $proj $range -W4,grey,- -m -O $V -K>>$psfile

awk '{if(($1==1)&& (int($2*10/5)==($2*10/5)) && int($3*10/5)==($3*10/5))  print $2,$3,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20}' $input > bb.txt

# beachball from hui
#psmeca bb.txt  -R -J -Sx0.09i -Ggrey -Ewhite  -O  -W2  -C -V -K>>$psfile

awk '{if(($3==1)&& int($1*10/2)==($1*10/2) && int($2*10/2)==($2*10/2))  print $1,$2,5*$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $inputybb > yang.txt

# beachball from yang
psmeca yang.txt -R  -J  -C -O -V  -Sd0.07i/1 -Gred -Ewhite -K>> $psfile

# about the direction of 2D principal stress
# our matlab output, extension--positive, compression--negative, but the sita equals to I1 counter-clockwise to x-axis.
# in gmt, psvelo, -Sx option, extension and compression definition the same. But the angle defition is different. 
# the angle in gmt, is azimuth of II2 in degrees CW (clockwise) from North.
# So we need transform the matlab output sita to GMT angle. I1=sita, I2=sita+90, angle=90-sita=90-sita-90=-sita, this's why -$9
awk '{if(($1==1)&& int($2*10/2)==($2*10/2) && int($3*10/2)==($3*10/2))  print $2,$3,$7,$8,-$9}' $inputprinstr > prin.txt
#awk '{if(($1==1) && int($2*10/2)==($2*10/2) )  print $2,$3,$7,$8,-$9}' $inputprinstr > prin0.txt
#awk '{if(int($2*10/2)==($2*10/2))   print $1,$2,$3,$4,$5}' prin0.txt > prin.txt
# -A default -A0.03/0.12/0.09
psvelo prin.txt -A0.01/0.6/0.1 -W0.25p,blue  -R  -J   -O -V  -Sx0.004  -K>> $psfile


#---------end of plot ---------
psxy -J -R -T -O >> $psfile
### clean the trash
rm bb.txt prin.txt yang.txt prin0.txt

#gnome-open ${title}.ps
open ${title}.ps
#mv *.ps ./output
