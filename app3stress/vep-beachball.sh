#!/bin/sh

# about beachball.txt  
# this file includes the data to plot beachball, calculated from Yang's method.
# the original file is from '/home/jiyang/opensource/vep-plot/app3stress/ori_results/app32', calculated by 'p12345.f'

input=./input/sagrid.32.txt
#input=./ores.$1/sagrid.$1.txt
#prepare the data
#title=$1.beachball.vep
title=32.beachball.vep
psfile=${title}.ps
psfile1=test.ps

range="-R75/225/100/300"
proj="-JX6i/8i"
ticks="-B25f25eSWn"


psbasemap ${range} ${proj} $ticks -X3.5c -Y3.0c -P -V -K > ${psfile}
#plot the faults
psxy ./para/ss3.fault.gmt  $proj $range -W9/0/0/0 -m -O -V -K>>$psfile
psxy ./para/profile.gmt  $proj $range -W0.1p,- -m -O -V -K  >>$psfile
awk '{if(($1==1)&& int(($2+5000)*10/100000)==(($2+5000)*10/100000) && int($3*10/100000)==($3*10/100000))  print $2/1000,$3/1000,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20}' $input > bb.txt
#awk '{if(($1==1)&&((int(($2*100)/5)*10000==($2*100)/5*10000)&&(int(($3*100)/5)*10000==($3*100)/5*10000))  print $2,$3,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20}' $input > bb.txt
# looks like a bug here, try $3*100/5
#awk '{if(($1==1)&&(int($2*100/5)==$2*100/5)&&(int($3*100/50)==($3*100/50)))  
###awk '{if(($1==1)&&(int($2*100/5)==$2*100/5)&&(int($3*1000/50)==($3*1000/50)))  
###	print $2,$3,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20,$2,$3 ;
###	else
###	print $2,$3,int($2*100/5), $2*100/5,  int(($3*1000)/50), $3*1000/50;
###	}' $input > bb.txt

psmeca bb.txt  -R -J -Sx0.09i -Ggrey -Ewhite  -O  -W2  -C -V -K>>$psfile
#psmeca bb.txt  -R -J -Sx0.60 -G0/0/0 -K -O -o-1 -W2 -L -C >>$psfile

awk '{if(int(($1+5000)*10/100000)==(($1+5000)*10/100000) && int($2*10/100000)==($2*10/100000))  print $1/1000,$2/1000,$3,$4,$5,$6,$7,$8,$9,$10,$11/1000,$12/1000+5}' beachball.txt > yang.txt

psmeca yang.txt -R  -J  -C -O -V  -Sd0.15i/1 -G255/0/0 >> $psfile

#
#psbasemap ${range} ${proj} $ticks -X3.5c -Y3.0c -P -V -K > $psfile1
##plot the faults
#psxy ./para/ss3.fault.gmt  $proj $range -W9/0/0/0 -m -O -V -K>>$psfile1
#psxy ./para/profile.gmt  $proj $range -W0.1p,- -m -O -V -K  >>$psfile1
#awk '{if(int(($1+5000)*10/100000)==(($1+5000)*10/100000) && int($2*10/100000)==($2*10/100000))  print $1/1000,$2/1000,$3,$4,$5,$6,$7,$8,$9,$10,$11/1000,$12/1000}' beachball.txt > yang.txt
#
#psmeca yang.txt -R  -J  -C -O -V  -Sd0.15i/1 -G255/0/0 >> $psfile1

### clean the trash
rm bb.txt

gnome-open ${title}.ps
#gnome-open $psfile1
#mv *.ps ./output
