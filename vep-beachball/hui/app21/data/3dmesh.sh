#!/bash/sh

# /*-----------------------------------------------------------------------------------------*/
# /*-----------------------------------------------------------------------------------------*/
# /*-----------------------------------------------------------------------------------------*/
output=mesh.ps
perspective=-E-135/30
region=-R-1/400/-1/600/0/60
project=-JX4i/6i

psxyz 3dgrid.txt $region $project -JZ1.0i -W3/150/150/150 -C3dgrid.cpt	-Be100/100/10 -H -L -G255 -M $perspective -K > $output

psxyz <<END $region $project -JZ -W8/0/0/0 $perspective -O -K >> $output
150  300  60
150  600  60 
END
psxyz <<END $region $project -JZ -W8/0/0/0 $perspective -O >> $output
250   0   0
250   0  60
250 300  60 
END

#ps2pdf $output
ps2raster $output -A -P -Tg

exit
