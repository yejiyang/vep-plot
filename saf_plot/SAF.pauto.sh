#!/bin/sh
# pauto means I1 plot with auto color scale
gmtset FRAME_WIDTH 0.05c TICK_LENGTH -0.075i ANNOT_FONT_SIZE 12p PLOT_DEGREE_FORMAT ddd:mm:ss
gmtset COLOR_NAN white
#title=safpap5
title=$1
efffault=$2
eb=$3
nv=$4
m=$5

cpt_path=../cpt/
bin_path=../bin/
para_path=../para/
eq_path=../eq/
ores_path=../ores/

input=cgrid.$m.txt
#prepare the data

psfile=${title}.pauto.ps #${psfile}

range="-R-121.6/-113.2/31.4/38.1"
proj="-JM8.39c"
ticks="-B2f1eSWn"
#V="-V"
V=""

# p- auto----------------------------------------------------------------------------
psbasemap ${range} ${proj} $ticks -X1.5c -Y19.0c -P $V -K > ${psfile}

awk -F ' ' '{print $2,$3,$8/1e6}' $input  > mean.xyz
xyz2grd mean.xyz -Gdata.grd -I0.05 -R
grdsample data.grd -Gdata1.grd -I0.01 -Ql.1

grd2cpt data1.grd -Crainbow $V -D>$title.p.cpt
mv $title.p.cpt ${cpt_path}

grdsample data1.grd $range -I0.00833333/0.00833333 -F -Gdata2.grd
mv data2.grd data1.grd
grdmath ${para_path}scatopo.grd -2000 GT 0. NAN = mask.grd $V
mv data1.grd tmp.grd
grdmath tmp.grd mask.grd OR = data1.grd

#grdimage ${para_path}scatopo_gradient1.grd ${proj} ${range} -C${para_path}grad.cpt -K $V -O>>${psfile}
grdimage data1.grd -I${para_path}scatopo_gradient1.grd  -C${cpt_path}$title.p.cpt $proj -R -K $V -P -O >>$psfile

#plot the faults
#psxy fault_data.gmt $proj $range -W7/0/0/0 -m -O $V -K>>$psfile
psclip ${para_path}saf627bddegree.txt $range $proj -m -O -K >>$psfile
psxy   ${para_path}fault_data.gmt -m $proj $range -W0.5 -O $V -K >> $psfile
psxy   ${para_path}$efffault $proj $range -W3/0/0/0 -m -O $V -K>>$psfile
##for eb, if within plot dash line
if [ "$eb" = "1" ]; then
        psxy ${para_path}eb.gmt $proj $range -W3,- -m -O $V -K  >>$psfile
fi
if [ "$nv" = "1" ]; then
        psxy ${para_path}nv.gmt $proj $range -W3,- -m -O $V -K  >>$psfile
fi
##for eb
psclip -C -O -K >>$psfile
psxy ${para_path}saf627bddegree.txt $proj $range -W5,grey,- -m -O $V -K>>$psfile

echo  3.885  4.50 > tmp
echo  5.686  4.50 >> tmp
echo  5.686  5.502 >> tmp
echo  3.885  5.502 >> tmp
psxy tmp -R0/30/0/20 -Jx0.58i -G255 -K -L -O -W3/0/0/0 $V >> ${psfile}

pstext <<END -R -Jx -P  -O -K $V>>${psfile}
4.00 5.3 8 0 0 0 Mean normal stress 
4.09 5.1 8 0 0 0 Stress (Mpa)
END

#draw the whole title
pstext <<END -R/0/4/0/4 -Jx1i -P -S0.5P -O -K $V >>$psfile
0.00  3.05 10 0 0 0 $title 
END

gmtset  ANNOT_FONT_SIZE 6p TICK_LENGTH -0.1905c TICK_PEN 0.05p
gmtset FRAME_WIDTH 0.5c  FRAME_PEN 0.35p
psscale -D7.035c/7.3c/2.06c/0.18ch -C${cpt_path}$title.p.cpt -E  -B15 -S  -O -K $V>> ${psfile}
gmtset  ANNOT_FONT_SIZE 12p

#draw earthquake
#eq_data=./eq/caeq.dat
psxy ${eq_path}eq5.gmt $range $proj -Sc0.1  -W0.8/grey -O -P $V -K>>$psfile
psxy ${eq_path}eq6.gmt $range $proj -Sc0.2  -W0.8/grey -O -P $V -K>>$psfile
psxy ${eq_path}eq7.gmt $range $proj -Sc0.3  -W0.8/grey -O -P $V -K>>$psfile

psxy -J -R -T -O >> $psfile
###########
rm *.grd mean.xyz tmp

