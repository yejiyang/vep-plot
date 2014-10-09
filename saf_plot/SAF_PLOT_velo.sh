#!/bin/sh

plotAll=1
plotVelo=1
plotE=1
plotESlip=1
plotJ2=1
plotJ2auto=1
plotPauto=1
plotSlip=1
plotUplift=1


#for k in `cat templist`
while read line
do
	input=$line
	i=$(echo  $input | awk '{ print $1 }')
	eb=$(echo $input | awk '{ print $2 }')
	nv=$(echo $input | awk '{ print $3 }')
done<toP

##plot muitilayer results
#for m in 1 2 3 4 5 6 7 8
for m in 1 
do
	if [ "$plotVelo" = "1" ]; then  # plot velo
	cd res.saf$i
	echo "running SAF.velo.sh ..."
	../SAF.velo.sh saf$i $i.fault.gmt $eb $nv $m
	mv saf$i.velo.ps  saf$i.velo.$m.ps
	cd ..
	fi

	if [ "$plotUplift" = "1" ]; then  # plot uplift
	cd res.saf$i
	echo "running SAF.uplift.sh ..."
	../SAF.uplift.sh saf$i $i.fault.gmt $eb $nv $m
	mv saf$i.uplift.ps  saf$i.uplift.$m.ps
	cd ..
	fi

	if [ "$plotE" = "1" ]; then  # plot plastic strain
	cd res.saf$i
	echo "running SAF.e.sh ..."
	../SAF.e.sh saf$i $i.fault.gmt $eb $nv $m
	mv saf$i.e.ps  saf$i.e.$m.ps
	cd ..
	fi

	if [ "$plotJ2" = "1" ]; then  # plot maximum stress
	cd res.saf$i
	echo "running SAF.j2.sh ..."
	../SAF.j2.sh saf$i $i.fault.gmt $eb $nv $m
	mv saf$i.j2.ps  saf$i.j2.$m.ps
	cd ..
	fi

	if [ "$plotJ2auto" = "1" ]; then  # plot maximum stress, auto scale
	cd res.saf$i
	echo "running SAF.j2auto.sh ..."
	../SAF.j2auto.sh saf$i $i.fault.gmt $eb $nv $m
	mv saf$i.j2auto.ps  saf$i.j2auto.$m.ps
	cd ..
	fi

	if [ "$plotPauto" = "1" ]; then  # plot mean normal stress, auto scale
	cd res.saf$i
	echo "running SAF.pauto.sh ..."
	../SAF.pauto.sh saf$i $i.fault.gmt $eb $nv $m
	mv saf$i.pauto.ps  saf$i.pauto.$m.ps
	cd ..
	fi
done
#
