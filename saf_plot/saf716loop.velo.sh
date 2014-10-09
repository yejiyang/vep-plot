#!/bin/sh
gmtset FRAME_WIDTH 0.05c TICK_LENGTH -0.075i ANNOT_FONT_SIZE 12p PLOT_DEGREE_FORMAT ddd:mm:ss
#title=safpap5
title=$1
efffault=$2
eb=$3
nv=$4
m=$5

#cpt_path=/home/jiyang/vepmpich2-June302013/results/cpt/
#bin_path=/home/jiyang/vepmpich2-June302013/results/bin/
#para_path=/home/jiyang/vepmpich2-June302013/results/para/
#eq_path=/home/jiyang/vepmpich2-June302013/results/eq/
#ores_path=/home/jiyang/vepmpich2-June302013/results/ores/
cpt_path=../cpt/
bin_path=../bin/
para_path=../para/
eq_path=../eq/
ores_path=../ores/

input=cgrid.$m.txt
#prepare the data

psfile=${title}.velo.ps #${psfile}

#range="-R-125/-115/32/40"
range="-R-121.6/-113.2/31.4/38.1"
#range="-R-121.6/-113.2/31.3/38.1"
proj="-JM8.39c"
ticks="-B2f1eSWn"
#V="-V"
V=""

psbasemap ${range} ${proj} $ticks -X1.5c -Y19.0c -P $V -K > ${psfile}
#psbasemap ${range} ${proj} -B4f2 -X1.5c -Y19.0c -P $V -K > ${psfile}

awk -F ' ' '{print $2, $3,sqrt($4*$4+$5*$5)/5}' ${input} > mean.xyz
#xyz2grd mean.xyz -Gdata.grd -I0.1 ${range} $V
#xyz2grd mean.xyz -Gdata.grd -I0.1 `grdinfo ${para_path}scatopo_gradient1.grd -I0.1/0.1` $V
xyz2grd mean.xyz -Gdata.grd -I0.05 $range $V  -NNaN
grdsample data.grd -Gdata1.grd -I0.01 -Ql.1 $range
#surface mean.xyz -I0.01 $range -T.5 -Li0. -Gsurf.grd $V
#grdmath surf.grd data1.grd OR = data2.grd
#grdsample data2.grd $range -I0.00833333/0.00833333 -F -Gdata3.grd
#mv data3.grd data1.grd 
grdsample data1.grd $range -I0.00833333/0.00833333 -F -Gdata2.grd
mv data2.grd data1.grd 

#apath=/home/jiyang/vepmpich2-June302013/results/sc_strain/

grdmath ${para_path}scatopo.grd -2000 GT 0. NAN = mask.grd $V
mv data1.grd tmp.grd
grdmath tmp.grd mask.grd OR = data1.grd

grdimage ${para_path}scatopo_gradient1.grd ${proj} ${range} -C${para_path}grad.cpt -K $V -O>>${psfile}
grdimage data1.grd -I${para_path}scatopo_gradient1.grd -C${cpt_path}vect_hui.cpt ${range} ${proj} -K $V  -O >>${psfile}
#grdimage data1.grd  -C${cpt_path}vect_hui.cpt ${range} ${proj} -K $V  -O >>${psfile}
#plot the faults
psclip ${para_path}saf627bddegree.txt $range $proj -m -O -K >>$psfile
#psxy fault_data.gmt $proj $range -W7/0/0/0 -m -O $V -K>>$psfile
psxy ${para_path}fault_data.gmt -m $proj $range -W1 -O $V -K >> $psfile
psxy ${para_path}$efffault $proj $range -W5/0/0/0 -m -O $V -K>>$psfile
##for eb, if within plot dash line
if [ "$eb" = "1" ]; then
        psxy ${para_path}eb.gmt $proj $range -W1.3p,- -m -O $V -K  >>$psfile
fi
if [ "$nv" = "1" ]; then
        psxy ${para_path}nv.gmt $proj $range -W1.3p,- -m -O $V -K  >>$psfile
#        echo "nv = $nv **********************************************"
fi
##for eb
psclip -C -O -K >>$psfile
#psxy ${para_path}saf627bddegree.txt $proj $range -W7/0/0/0 -m -O $V -K>>$psfile
#psxy ${para_path}saf627bddegree.txt $proj $range -W4,snow,- -m -O $V -K>>$psfile

awk -F ' ' '{if((int($2*10/5)==$2*10/5)&&(int($3*10/5)==$3*10/5)) print $2,$3,$4/15*0.6,$5/15*0.6,0,0,0}' ${input} |psvelo ${proj} ${range} -G100/100/100 -A0.025/0.150/0.05 -Se0.045/0.95/0 -W2/0/0/0 -K -O $V >>${psfile}

pstext <<END ${range} ${proj} -S0.5p -O -K $V >>${psfile}
#-121.8  34.6 11 0 0 0 Pacific  
#-120.8  33.9 11 0 0 0 Plate    
#-119.0  37.3 11 0 0 0 North    
#-118.5  36.6 11 0 0 0 American 
#-118.0  35.9 11 0 0 0 Plate    
#-121.0  32.0 11 0 0 0 (a) 
END

echo  3.885  4.50 > tmp
echo  5.686  4.50 >> tmp
echo  5.686  5.502 >> tmp
echo  3.885  5.502 >> tmp
psxy tmp -R0/30/0/20 -Jx0.58i -G255 -K -L -O -W3/0/0/0 $V >> ${psfile}

pstext <<END -R -Jx -P  -O -K $V>>${psfile}
4.09 5.2 8 0 0 0 Velocity (mm/yr) 
END

#draw the whole title
pstext <<END -R/0/4/0/4 -Jx1i -P -S0.5P -O -K $V >>$psfile
0.00  3.05 10 0 0 0 $title 
END

gmtset  ANNOT_FONT_SIZE 6p TICK_LENGTH -0.1905c TICK_PEN 0.05p
gmtset FRAME_WIDTH 0.5c  FRAME_PEN 0.35p
psscale -D7.035c/7.3c/2.06c/0.18ch -C${cpt_path}vect_hui.cpt -E -I0.3  -B8  -O -K $V>> ${psfile}
gmtset  ANNOT_FONT_SIZE 12p

psxy -J -R -T -O >> $psfile
###########
rm *.grd mean.xyz tmp

#ps2pdf ${psfile} ${title}.velo.pdf
#ps2raster $psfile -A -P -Tg

#gnome-open  ${psfile}
