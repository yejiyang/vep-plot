#!/bin/sh

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
for m in 1 6
do
		cd res.saf$i
		../saf716loop.velo.sh saf$i $i.fault.gmt $eb $nv $m
		echo "******* end of  saf716loop.velo.sh *************"
		mv saf$i.velo.ps  saf$i.velo.$m.ps
		cd ..
done
#
