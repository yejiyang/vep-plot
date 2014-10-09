#!/bin/sh
#
bin_path=./bin/
para_path=./para/
cpt_path=./cpt/
eq_path=./eq/
ores_path=./ores/

#input="716_1_50000     saf"
#	echo " input=:  $input"
# read res.list and remove # item
cat paper.res.list | sed '/#/d' >templist

#for k in `cat templist`
while read line
#for2
do
input=$line

#
#----------step 1:  check the input file------------------
#----------1.1 check the case original data---------------
#
	i=$(echo $input | awk '{ print $1 }')
#	# i, case name, should equal to 716_1_50000
#	# determine whether the ores.saf$i.tgz exists or not, 
#	# if doesn't exist, quit current loop, output error info
	if [ -e ores.saf$i.tgz ]; then
		echo "ores.saf$i.tgz exists, continue ..."
	else 
		## echo "ores.saf$i.tgz doesn't exists, exit current process, jump to next one"
		## continue
		echo "ores.saf$i.tgz doesn't exists, exit the  process!"
		exit 1;
		#continue;
	fi

#----------1.2 check the effective faults---------------
#--- get the effective fault name, especially for NV and ECSZ, will plot with dash line;
#--- create the fault.gmt file for later plot.
#
#	# eb, for determine whether plot eb on figure
#	# nv, for determine whether plot navada block on figure
	eb=0;	nv=0
#
#	# j store the effective faults in the model
	j=$(echo $input | awk '{ print $2 }')
#
#
#	#if eb and nv, using dash line
	if echo $j | grep -Fq "eb"; then
		eb=1;
	fi
	if echo $j | grep -Fq "nv"; then nv=1; fi;

#	# the other faults plot solid line
#	# including, saf,sjf,gf,pf,tf
	rm $i.fault.gmt
	if echo $j | grep -Fq "saf"; then
		cat ${para_path}saf.gmt >>$i.fault.gmt;
	fi 
	if echo $j | grep -Fq "sjf"; then cat ${para_path}sjf.gmt >>$i.fault.gmt; fi; 
	if echo $j | grep -Fq "gf"; then cat ${para_path}gf.gmt >>$i.fault.gmt; fi; 
	if echo $j | grep -Fq "pf"; then cat ${para_path}pf.gmt >>$i.fault.gmt; fi; 
	if echo $j | grep -Fq "tf"; then cat ${para_path}tf.gmt >>$i.fault.gmt; fi; 		
	cp $i.fault.gmt ${para_path}
#	echo "eb= $eb"
#	echo "nv= $nv"
#
#----------1.3 rearrange data-----------------------------
#	# datapost
#	cp ores.saf$i.tgz ./
	tar -xzvf ores.saf$i.tgz
#	rm ores.saf$i.tgz
	gfortran  -o pp123 ${bin_path}pp123.f
	pp123 >& pp123saf$i.log
	s_mkdir res.saf$i
	mv  res.fl* pp123saf$i.log  res_xy.txt ./res.saf$i
	cp cgrid.m.txt ./res.saf$i/.
	rm -f amend.flavia.msh coor0 data1 elem0 munod munod1
#
#
#----------2 prepare the plot data--------------
#
#
##plot muitilayer results
#for m in 1 2 3 4 5 6 7 8
for m in 1 6
#start for1
do
	echo $m
	awk '{if($1=='$m') print; }' cgrid.m.txt > cgrid.txt
	awk '{if($1=='$m') print; }' slip.m > slip
#
#------@@@-process slip, plot fault slip rate only on the effected faults
	cp slip slip.bak
# 	#j the faults in the model
#	#j=$(echo $input | awk '{ print $2 }')
#	echo "j= $j"
	rm slip
#	#if f1 match with some fault, than output this fault's slip rate
	if echo $j | grep -Fq "saf"; then
		awk  'NR==1, NR==7 {print;}' slip.bak >> slip;
	fi 
	if echo $j | grep -Fq "sjf"; then awk  'NR==13,  NR==14 {print;}' slip.bak >> slip; fi;
	if echo $j | grep -Fq "gf";  then awk  'NR==8,   NR==9  {print;}' slip.bak >> slip; fi;
	if echo $j | grep -Fq "pf";  then awk  'NR==15,  NR==18 {print;}' slip.bak >> slip; fi;
	if echo $j | grep -Fq "tf";  then awk  'NR==10,  NR==12 {print;}' slip.bak >> slip; fi;
	if echo $j | grep -Fq "eb";  then awk  'NR==19,  NR==20 {print;}' slip.bak >> slip; fi;
#	#always plot the slip rate on eb
#	awk  'NR==19, NR==20 {print;}' slip.bak >> slip
#
#-------@@@process slip
		mv cgrid.txt ./res.saf$i/cgrid.$m.txt
		mv slip ./res.saf$i/slip.$m
		mv slip.bak ./res.saf$i/slip.$m.bak
#end for1
done
#
#
	mv slip.m ./res.saf$i/
	rm  slip*
	mv $i.fault.gmt ./res.saf$i
	rm  cgrid.m.txt
#endfor2
done<templist

echo $i $eb $nv > toP
#mv toP ./res.saf$i/

rm templist
