#!/bin/sh
##gmtset FRAME_WIDTH 0.05c TICK_LENGTH -0.075i ANNOT_FONT_SIZE 12p PLOT_DEGREE_FORMAT ddd:mm:ss

input=./input/sagrid.32.txt
#input=./ores.$1/sagrid.$1.txt
#prepare the data
#title=$1.beachball.vep
title=32.beachball.vep
psfile=${title}.ps

#range="-R-125/-115/32/40"
#range="-R-121.6/-113.2/31.4/38.1"
#range="-R-116.5/-115/32.5/34"
#range="-R-116.25/-115.25/32.75/33.5"
#proj="-JM15.39c"
#ticks="-B0.25f0.125eSWn"

range="-R75/225/100/300"
proj="-JX6i/8i"
ticks="-B25f25eSWn"


psbasemap ${range} ${proj} $ticks -X3.5c -Y3.0c -P -V -K > ${psfile}
#plot the faults
psxy ./para/ss3.fault.gmt  $proj $range -W9/0/0/0 -m -O -V -K>>$psfile
psxy ./para/profile.gmt  $proj $range -W0.1p,- -m -O -V -K  >>$psfile
#pscoast ${proj} ${range} -Dh -A10  -I8 -N1 -N2 -W0.5p,black -O  -K -V>>${psfile}
#psxy ../para/fault_data.gmt -m $proj $range -W6,blue -K -V -P -O >> $psfile
#psxy ../para/fault_BSZ.gmt -m $proj $range -W6,blue,- -K -V -P -O >> $psfile
awk '{if(($1==1)&& int(($2+5000)*10/100000)==(($2+5000)*10/100000) && int($3*10/100000)==($3*10/100000))  print $2/1000,$3/1000,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20}' $input > bb.txt
#awk '{if(($1==1)&&((int(($2*100)/5)*10000==($2*100)/5*10000)&&(int(($3*100)/5)*10000==($3*100)/5*10000))  print $2,$3,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20}' $input > bb.txt
# looks like a bug here, try $3*100/5
#awk '{if(($1==1)&&(int($2*100/5)==$2*100/5)&&(int($3*100/50)==($3*100/50)))  
###awk '{if(($1==1)&&(int($2*100/5)==$2*100/5)&&(int($3*1000/50)==($3*1000/50)))  
###	print $2,$3,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20,$2,$3 ;
###	else
###	print $2,$3,int($2*100/5), $2*100/5,  int(($3*1000)/50), $3*1000/50;
###	}' $input > bb.txt
psmeca bb.txt  -R -J -Sx0.09i -Gblack -Ewhite  -O  -W2  -C -V>>$psfile
#psmeca bb.txt  -R -J -Sx0.60 -G0/0/0 -K -O -o-1 -W2 -L -C >>$psfile

### clean the trash
rm bb.txt

ps2pdf ${psfile} ${title}.pdf
ps2raster $psfile -A -P -Tg
mv *.pdf ./output
mv *.ps ./output
mv *.png ./output
#gnome-open ${title}.pdf
#cp $psfile ~/Dropbox/temp/
