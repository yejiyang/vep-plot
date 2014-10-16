#!/bin/sh

psfile=beachball.ps

psbasemap -R119.0/123.0/20.5/25.5 -JM3i -B1/1WSeN -P -K > $psfile

#psmeca -R -J -O earthmech4-m.txt -Sy0.13i -W -Gblack -Ewhite -C -K >> $psfile
psmeca -R -J -O earthmech4-m.txt -Sy0.13i -W -Gblack -Ewhite -C  >> $psfile

#psxy -R -J -O -K -G240 -L -Wthicker box >> $psfile
#psmeca -R -J -O epicenter.txt -Sy0.13i -W -C -K >> $psfile
#echo 122.7	20.6	10	0	4	2	10MPa | pstext -R -J -O -K >> $psfile
#echo 122.5	20.9	12	0	4	2	10km  | pstext -R -J -O -K >> $psfile

gnome-open $psfile
