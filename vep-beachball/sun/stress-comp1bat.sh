#!/bin/sh
psfile=sun_beachball.ps
#makecpt -Csealand -D-7500/2500/100 -Z > topography.cpt
psbasemap -R119.0/123.0/20.5/25.5 -JM3i -B1/1WSeN -P -K > $psfile
#grdimage ETOPO2v2g_f4.nc -IETOPO2v2g_f4.grd -R119.0/123.0/20.5/25.5 -JM3i -Ctopography.cpt -B -O -K -Ei >> $psfile
#  psxy -R -J -O -M faultR.gen -Gwhite -Sf0.4i/0.05ilt -W0.01i/red -B -K >> $psfile
#  psxy -R -J -O -M faultL.gen -Gwhite -Sf0.4i/0.05irt -W0.01i/red -B -K >> $psfile
psmeca -R -J -O earthmech4-m.txt -Sy0.13i -W -Gblack -Ewhite -C -K >> $psfile
echo 122.0 20.5 > box.d
echo 123.0 20.5 >> box.d
echo 123.0 21.1 >> box.d
echo 122.0 21.1 >> box.d
psxy -R -J -O -K -G240 -L -Wthicker box.d>> $psfile
psmeca -R -J -O epicenter.txt -Sy0.13i -W -C -K >> $psfile
echo 122.7	20.6	10	0	4	2 10MPa | pstext -R -J -O -K >> $psfile
echo 122.5	20.9	12	0	4	2 10km | pstext -R -J -O -K >> $psfile
pscoast -R -JM -B -Dc  -W0/0/255 -Di -P -O -K >> $psfile
echo 119.0 20.5 > box.d
echo 119.6 20.5 >> box.d
echo 119.6 21.0 >> box.d
echo 119.0 21.0 >> box.d
psxy -R -J -O -K -G240 -L -Wwhite box.d >> $psfile
echo 119.1 20.8 16 0.0 1 5 \(a\) | pstext -R -J -O -K >> $psfile
#####
# pause
psbasemap -R119.0/123.0/20.5/25.5 -JM3i -B1/1wSEN -P -X3.4i -O -K >> $psfile
#grdimage ETOPO2v2g_f4.nc -IETOPO2v2g_f4.grd -R119.0/123.0/20.5/25.5 -JM3i -Ctopography.cpt -B -O -K -Ei >> $psfile
# psxy -R -J -O -M faultR.gen -Gwhite -Sf0.4i/0.05ilt -W0.01i/red -B -K >> $psfile
# psxy -R -J -O -M faultL.gen -Gwhite -Sf0.4i/0.05irt -W0.01i/red -B -K >> $psfile
#psmeca -R -J -O batsmech15.txt -Sa0.16i -W -C -K >> $psfile
echo 122.0 20.5 > box.d
echo 123.0 20.5 >> box.d
echo 123.0 21.1 >> box.d
echo 122.0 21.1 >> box.d
psxy -R -J -O -K -G240 -L -Wthicker box.d >> $psfile
#psmeca -R -J -O batssign.txt -Sa0.16i -W -C -K >> $psfile
echo 122.7	20.6 10 0.0 4 2 Mw=6 | pstext -R -J -O -K >> $psfile
echo 122.5	20.9	12	0	4	2 0~15km | pstext -R -J -O -K >> $psfile
pscoast -R -JM -B -Dc  -W0/0/255 -Di -P -O -K >> $psfile 
echo 119.0 20.5 > box.d
echo 119.6 20.5 >> box.d
echo 119.6 21.0 >> box.d
echo 119.0 21.0 >> box.d
psxy -R -J -O -K -G240 -L -Wwhite box.d >> $psfile
echo 119.1 20.8 16 0.0 1 5 \(b\) | pstext -R -J -O -K >> $psfile
#pause
cp $psfile ~/Dropbox/temp/
