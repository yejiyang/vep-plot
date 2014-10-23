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

range="-R75/225/100/300"
proj="-JX6i/8i"
title2=1saf_vs_3sjf
ticks="-B25f25eSWn:.$title2:"


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

# beachball from hui
#psmeca bb.txt  -R -J -Sx0.09i -Ggrey -Ewhite  -O  -W2  -C -V -K>>$psfile

awk '{if(int(($1+5000)*10/100000)==(($1+5000)*10/100000) && int($2*10/100000)==($2*10/100000))  print $1/1000,$2/1000,$3,$4,$5,$6,$7,$8,$9,$10,$11/1000,$12/1000}' beachball.txt > yang.txt

# beachball from yang
psmeca yang.txt -R  -J  -C -O -V  -Sd0.25i/1 -Gred -Ewhite -K>> $psfile

awk '{if(int(($2+5)*10/100)==(($2+5)*10/100) && int($3*10/100)==($3*10/100))  print $2,$3,$7,$8,-$9}' prin_stress1.txt > prin.txt
# -A default -A0.03/0.12/0.09
psvelo prin.txt -A0.01/0.6/0.1 -W0.25p,blue  -R  -J   -O -V  -Sx0.06  -K>> $psfile


#---------end of plot ---------
psxy -J -R -T -O >> $psfile
### clean the trash
rm bb.txt prin.txt yang.txt

#gnome-open ${title}.ps
open ${title}.ps
#gnome-open $psfile1
#mv *.ps ./output
