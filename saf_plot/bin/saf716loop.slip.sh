#!/bin/sh
gmtset FRAME_WIDTH 0.05c TICK_LENGTH -0.075i ANNOT_FONT_SIZE 12p PLOT_DEGREE_FORMAT ddd:mm:ss
#title=safpap5
title=$1
efffault=$2
eb=$3
nv=$4

cpt_path=/home/jiyang/vepmpich2-June302013/results/cpt/
bin_path=/home/jiyang/vepmpich2-June302013/results/bin/
para_path=/home/jiyang/vepmpich2-June302013/results/para/
eq_path=/home/jiyang/vepmpich2-June302013/results/eq/
ores_path=/home/jiyang/vepmpich2-June302013/results/ores/


input=cgrid.txt
#prepare the data

psfile=${title}.slip.ps #${psfile}

#range="-R-125/-115/32/40"
range="-R-121.6/-113.2/31.4/38.1"
proj="-JM8.39c"
ticks="-B2f1eSWn"

psbasemap ${range} ${proj} $ticks -X1.5c -Y19.0c -P -V -K > ${psfile}
#psbasemap ${range} ${proj} -B4f2 -X1.5c -Y19.0c -P -V -K > ${psfile}

grdimage ${para_path}scatopo_gradient1.grd  $proj $range -C${para_path}grad.cpt   -O -K -V >>$psfile

psclip ${para_path}saf627bddegree.txt $range $proj -m -O -K >>$psfile
grdimage ${para_path}scatopo.grd -I${para_path}scatopo_gradient.grd  $proj $range -C${cpt_path}ca1.cpt  -P -O -K -V >>$psfile
#grdimage ./para/scatopo_gradient.grd  $proj $range -C./cpt/ca1.cpt  -P -O -K -V >>$psfile
psclip -C -O -K >>$psfile
psxy ${para_path}fault_data.gmt -m $proj $range -W1 -O -V -K >> $psfile

#plot the faults
#psxy fault_data.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psclip ${para_path}saf627bddegree.txt $range $proj -m -O -K >>$psfile
psxy ${para_path}$efffault $proj $range -W7/0/0/0   -m -O -V -K>>$psfile
##for eb, if within plot dash line
if [ "$eb" = "1" ]; then
        psxy ${para_path}eb.gmt $proj $range -W1.3p,- -m -O -V -K  >>$psfile
fi
if [ "$nv" = "1" ]; then
        psxy ${para_path}nv.gmt $proj $range -W1.3p,- -m -O -V -K  >>$psfile
        echo "nv = $nv **********************************************"
fi
##for eb
psclip -C -O -K >>$psfile
psxy ${para_path}saf627bddegree.txt $proj $range -W6/0/0/0 -m -O -V -K>>$psfile

#label the slip rate
awk -F ' ' '{print $2, $3, 8, 0, 0, 0, $4}' slip > slip.gmt
#pstext slip.gmt  -R $proj -S0.5p -O -K -Gblue -Sthinnest >>$psfile
pstext slip.gmt  -R $proj -Gred -O -K  >>$psfile

pstext <<END -R $proj -S0.5p -O -K >>$psfile
#-121.0  32.0 11 0 0 0 (d) 
END

#draw the whole title
pstext <<END -R/0/4/0/4 -Jx1i -P -S0.5P -O -K -V >>$psfile
0.00  3.05 10 0 0 0 $title 
END

psxy -J -R -T -O >> $psfile
###########
rm *.grd mean.xyz tmp #slip.gmt

#gs $psfile
ps2pdf ${psfile} ${title}.slip.pdf
ps2raster $psfile -A -P -Tg
#rcp *.png geo01:~/Dropbox/temp
