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

psfile=${title}.j2.ps #${psfile}

#range="-R-125/-115/32/40"
range="-R-121.6/-113.2/31.4/38.1"
proj="-JM8.39c"
ticks="-B2f1eSWn"

psbasemap ${range} ${proj} $ticks -X1.5c -Y19.0c -P -V -K > ${psfile}
#psbasemap ${range} ${proj} -B4f2 -X1.5c -Y19.0c -P -V -K > ${psfile}
#rem3 -j2-----------------------------------------------------------------------------

awk -F ' ' '{print $2,$3,$7/1e6}' $input  > mean.xyz
xyz2grd mean.xyz -Gdata.grd -I0.05 -R
grdsample data.grd -Gdata1.grd -I0.01 -Ql.1

grdsample data1.grd $range -I0.00833333/0.00833333 -F -Gdata2.grd
mv data2.grd data1.grd
grdmath ${para_path}scatopo.grd -2000 GT 0. NAN = mask.grd -V
mv data1.grd tmp.grd
grdmath tmp.grd mask.grd OR = data1.grd

grdimage ${para_path}scatopo_gradient1.grd ${proj} ${range} -C${para_path}grad.cpt -K -V -O>>${psfile}
grdimage data1.grd -I${para_path}scatopo_gradient1.grd -C${cpt_path}shear_hui.cpt $proj -R -K -V -P -O >>$psfile

#plot the faults
#psxy fault_data.gmt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
psclip ${para_path}saf627bddegree.txt $range $proj -m -O -K >>$psfile
psxy   ${para_path}fault_data.gmt -m $proj $range -W2 -O -V -K >> $psfile
psxy   ${para_path}$efffault $proj $range -W7/0/0/0 -m -O -V -K>>$psfile
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
psxy ${para_path}saf627bddegree.txt $proj $range -W7/0/0/0 -m -O -V -K>>$psfile

pstext <<END -R $proj -S0.5p -O -K >>$psfile
#-121.8  34.6 11 0 0 0 Pacific  
#-120.8  33.9 11 0 0 0 Plate    
#-119.0  37.3 11 0 0 0 North    
#-118.5  36.6 11 0 0 0 American 
#-118.0  35.9 11 0 0 0 Plate    
#-121.0  32.0 11 0 0 0 (c) 
END

echo  3.885  4.50 > tmp
echo  5.686  4.50 >> tmp
echo  5.686  5.502 >> tmp
echo  3.885  5.502 >> tmp
psxy tmp -R0/30/0/20 -Jx0.58i -G255 -K -L -O -W3/0/0/0 -V >> ${psfile}

pstext <<END -R -Jx -P  -O -K -V>>${psfile}
4.09 5.3 8 0 0 0 Maximum Shear 
4.09 5.1 8 0 0 0 Stress (Mpa)
END

#draw the whole title
pstext <<END -R/0/4/0/4 -Jx1i -P -S0.5P -O -K -V >>$psfile
0.00  3.05 10 0 0 0 $title 
END
#psscale -C./cpt/ca.cpt -D4/-1.2/6/0.5h -B2000/:meters: -O -K >>$psfile

gmtset  ANNOT_FONT_SIZE 6p
psscale -D7.035c/7.3c/2.06c/0.18ch -C${cpt_path}shear_hui.cpt   -B15 -S  -O -K -V>> ${psfile}
#psscale -D4.235c/4.550c/1.96c/0.14ch -Cshear_hui.cpt -B15 -S -O -K >> $psfile

#draw earthquake
#eq_data=./eq/caeq.dat
psxy ${eq_path}eq5.gmt $range $proj -Sc0.1  -W0.8/black -O -P -V -K>>$psfile
psxy ${eq_path}eq6.gmt $range $proj -Sc0.2  -W0.8/black -O -P -V -K>>$psfile
psxy ${eq_path}eq7.gmt $range $proj -Sc0.3 -W0.8/black -O -P -V -K>>$psfile

psxy -J -R -T -O >> $psfile
###########
rm *.grd mean.xyz tmp

ps2pdf ${psfile} ${title}.j2.pdf
ps2raster $psfile -A -P -Tg

#gs ${psfile}
#rcp ${psfile} geo01:/home/jiyang/Dropbox/temp
#tar -czvf ${title}.fig.tgz ${psfile} ${title}.pdf ${title}.png
#rcp ${title}.fig.tgz geo01:/home/jiyang/Dropbox/SAF/res

#rm ${psfile} ${title}.pdf ${title}.png
