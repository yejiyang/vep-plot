#!/bin/sh
#gmtset FRAME_WIDTH 0.05c TICK_LENGTH -0.075i ANNOT_FONT_SIZE 12p PLOT_DEGREE_FORMAT ddd:mm:ss

input=cgrid.txt
#prepare the data
title=beachball.vep
psfile=${title}.ps

#range="-R-125/-115/32/40"
#range="-R-121.6/-113.2/31.4/38.1"
#range="-R-116.5/-115/32.5/34"
range="-R-116.25/-115.25/32.75/33.5"
proj="-JM15.39c"
ticks="-B0.25f0.125eSWn"

psbasemap ${range} ${proj} $ticks -X3.5c -Y3.0c -P -V -K > ${psfile}
#awk '{if(($1==1))  print $2,$3,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20}' $input > bb.txt
#awk '{if(($1==1)&&((int(($2*100)/5)*10000==($2*100)/5*10000)&&(int(($3*100)/5)*10000==($3*100)/5*10000))  print $2,$3,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20}' $input > bb.txt
# looks like a bug here, try $3*100/5
#awk '{if(($1==1)&&(int($2*100/5)==$2*100/5)&&(int($3*100/50)==($3*100/50)))  
awk '{if(($1==1)&&(int($2*100/5)==$2*100/5)&&(int($3*1000/50)==($3*1000/50)))  
	print $2,$3,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20,$2,$3 ;
	else
	print $2,$3,int($2*100/5), $2*100/5,  int(($3*1000)/50), $3*1000/50;
	}' $input > bb.txt
psmeca bb.txt  -R -J -Sx0.09i -Gblack -Ewhite  -O  -W2  -C >>$psfile
#psmeca bb.txt  -R -J -Sx0.60 -G0/0/0 -K -O -o-1 -W2 -L -C >>$psfile


psbasemap ${range} ${proj} $ticks -X3.5c -Y3.0c -P -V -K > yang.ps
awk '{if((int($1*100/5)==$1*100/5)&&(int($2*1000/50)==($2*1000/50)))  
	print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10,$11,$12 ;
	else
	print $2,$3,int($2*100/5), $2*100/5,  int(($3*1000)/50), $3*1000/50;
	}' beachball.txt >  yang.txt

psmeca yang.txt -R  -J  -C -O  -Sd0.15i/1 -G255/0/0 >> yang.ps

#ps2pdf ${psfile} ${title}.pdf
#ps2raster $psfile -A -P -Tg
gnome-open ${title}.pdf
gnome-open yang.ps
#cp $psfile ~/Dropbox/temp/
