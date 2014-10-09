#!/bin/sh
#
bin_path=/home/jiyang/vepmpich2-June302013/results/bin/
para_path=/home/jiyang/vepmpich2-June302013/results/para/
cpt_path=/home/jiyang/vepmpich2-June302013/results/cpt/
eq_path=/home/jiyang/vepmpich2-June302013/results/eq/
ores_path=/home/jiyang/vepmpich2-June302013/results/ores/
#
# clear figure folder, then creat one
#	s_mkdir  figure
#	s_mkdir  out_data
#
#
#input="716_1_50000     saf"
	echo " input=:  $input"
# read res.list and remove # item
#cat all.list | sed '/#/d' >templist
cat paper.res.list | sed '/#/d' >templist

#for k in `cat templist`
while read line
#for2
do
input=$line

#
#
#----------step 1:  check the input file------------------
#----------1.1 check the case original data---------------
#
	i=$(echo $input | awk '{ print $1 }')
#	# i, case name, should equal to 716_1_50000
#	# determine whether the ores.saf$i.tgz exists or not, 
#	# if doesn't exist, quit current loop, output error info
	if [ -e ${ores_path}ores.saf$i.tgz ]; then
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
	cp ${ores_path}ores.saf$i.tgz ./
	tar -xzvf ores.saf$i.tgz
#	rm ores.saf$i.tgz
	${bin_path}pp123 >& pp123saf$i.log
	s_mkdir res.saf$i
	mv  res.fl* pp123saf$i.log  res_xy.txt ./res.saf$i
	cp cgrid.m.txt ./res.saf$i/.
	rm -f amend.flavia.msh coor0 data1 elem0 munod munod1
#
#
#
#----------step 2:  plot the results------------------
#
##plot muitilayer results
for m in 1 2 3 4 5 6 7 8
#for m in 1
#start for1
do
#----------2.1 prepare the plot data--------------
	echo $m
	awk '{if($1=='$m') print; }' cgrid.m.txt > cgrid.txt
	awk '{if($1=='$m') print; }' slip.m > slip
#
#-------process slip, plot fault slip rate only on the effected faults
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
#-------process slip
#
#----------2.2  plot ------------------------
	plotAll=1
	plotVelocity=1
	plotE=1
	plotESlip=1
	plotJ2=1
	plotJ2auto=1
	plotPauto=1
	plotSlip=1
	plotUplift=1
#
#	plot all
#	echo "plotAll = $plotAll"
	if [ "$plotAll" = "1" ]; then	
		echo "************running saf716loop.sh ************************************"
		${bin_path}saf716loop.sh saf$i $i.fault.gmt $eb $nv
		echo "************end of  saf716loop.sh ************************************"
	#	mv cgrid.txt ./res.saf$i
		mv saf$i.pdf saf$i.$m.pdf
		mv saf$i.png saf$i.$m.png
		mv saf$i.ps  saf$i.$m.ps
		cp saf$i.$m.p* ./figure
		mv saf$i.$m.p* ./res.saf$i
	fi
#	plot velocity
	if [ "$plotVelocity" = "1" ]; then
		echo "************running saf716loop.velo.sh ************************************"
		${bin_path}saf716loop.velo.sh saf$i $i.fault.gmt $eb $nv
		echo "************end of  saf716loop.velo.sh ************************************"
	#	mv cgrid.txt ./res.saf$i
		mv saf$i.velo.pdf saf$i.velo.$m.pdf
		mv saf$i.velo.png saf$i.velo.$m.png
		mv saf$i.velo.ps  saf$i.velo.$m.ps
		cp saf$i.velo.$m.p* ./figure
		mv saf$i.velo.$m.p* ./res.saf$i
	fi
	#  plot e, plastic strain
	if [ "$plotE" = "1" ]; then
		echo "************running saf716loop.e.sh ************************************"
		${bin_path}saf716loop.e.sh saf$i $i.fault.gmt $eb $nv
		echo "************end of  saf716loop.e.sh ************************************"
	#	mv cgrid.txt ./res.saf$i
		mv saf$i.e.pdf saf$i.e.$m.pdf
		mv saf$i.e.png saf$i.e.$m.png
		mv saf$i.e.ps  saf$i.e.$m.ps
		cp saf$i.e.$m.p* ./figure
		mv saf$i.e.$m.p* ./res.saf$i
	fi
	#  plot e-slip, plastic strain
	if [ "$plotESlip" = "1" ]; then
		echo "************running saf716loop.e-slip.sh ************************************"
		${bin_path}saf716loop.e-slip.sh saf$i $i.fault.gmt $eb $nv
		echo "************end of  saf716loop.e-slip.sh ************************************"
	#	mv cgrid.txt ./res.saf$i
		mv saf$i.e-slip.pdf saf$i.e-slip.$m.pdf
		mv saf$i.e-slip.png saf$i.e-slip.$m.png
		mv saf$i.e-slip.ps  saf$i.e-slip.$m.ps
		cp saf$i.e-slip.$m.p* ./figure
		mv saf$i.e-slip.$m.p* ./res.saf$i
	fi
	#  plot j2
	if [ "$plotJ2" = "1" ]; then
		echo "************running saf716loop.j2.sh ************************************"
		${bin_path}saf716loop.j2.sh saf$i $i.fault.gmt $eb $nv
		echo "************end of  saf716loop.j2.sh ************************************"
	#	mv cgrid.txt ./res.saf$i/cgrid.$m.txt
		mv saf$i.j2.pdf saf$i.j2.$m.pdf
		mv saf$i.j2.png saf$i.j2.$m.png
		mv saf$i.j2.ps  saf$i.j2.$m.ps
		cp saf$i.j2.$m.p* ./figure
		mv saf$i.j2.$m.p* ./res.saf$i
	fi
	#  plot j2-auto
	if [ "$plotJ2auto" = "1" ]; then
		echo "************running saf716loop.j2auto.sh ************************************"
		${bin_path}saf716loop.j2auto.sh saf$i $i.fault.gmt $eb $nv
		echo "************end of  saf716loop.j2auto.sh ************************************"
	#	mv cgrid.txt ./res.saf$i/cgrid.$m.txt
		mv saf$i.j2auto.pdf saf$i.j2auto.$m.pdf
		mv saf$i.j2auto.png saf$i.j2auto.$m.png
		mv saf$i.j2auto.ps  saf$i.j2auto.$m.ps
		cp saf$i.j2auto.$m.p* ./figure
		mv saf$i.j2auto.$m.p* ./res.saf$i
	fi
	#  plot p-auto
	if [ "$plotPauto" = "1" ]; then
		echo "************running saf716loop.pauto.sh ************************************"
		${bin_path}saf716loop.pauto.sh saf$i $i.fault.gmt $eb $nv
		echo "************end of  saf716loop.pauto.sh ************************************"
	#	mv cgrid.txt ./res.saf$i/cgrid.$m.txt
		mv saf$i.pauto.pdf saf$i.pauto.$m.pdf
		mv saf$i.pauto.png saf$i.pauto.$m.png
		mv saf$i.pauto.ps  saf$i.pauto.$m.ps
		cp saf$i.pauto.$m.p* ./figure
		mv saf$i.pauto.$m.p* ./res.saf$i
	fi
	#  plot uplift
	if [ "$plotUplift" = "1" ]; then
		echo "************running saf716loop.uplift.sh ************************************"
		${bin_path}saf716loop.uplift.sh saf$i $i.fault.gmt $eb $nv
		echo "************end of  saf716loop.uplift.sh ************************************"
		mv cgrid.txt ./res.saf$i/cgrid.$m.txt
		mv saf$i.uplift.pdf saf$i.uplift.$m.pdf
		mv saf$i.uplift.png saf$i.uplift.$m.png
		mv saf$i.uplift.ps  saf$i.uplift.$m.ps
		cp saf$i.uplift.$m.p* ./figure
		mv saf$i.uplift.$m.p* ./res.saf$i
	fi
	#  plot slip
	if [ "$plotSlip" = "1" ]; then
		echo "************running saf716loop.slip.sh ************************************"
		${bin_path}saf716loop.slip.sh saf$i $i.fault.gmt $eb $nv
		echo "************end of  saf716loop.slip.sh ************************************"
		mv cgrid.txt ./res.saf$i/cgrid.$m.txt
		mv slip  ./res.saf$i/slip.$m
		mv slip.bak  ./res.saf$i/slip.$m.bak
		mv slip.gmt  ./res.saf$i/slip.$m.gmt
#		mv slip.m  ./res.saf$i/slip.$m.m
		mv saf$i.slip.pdf saf$i.slip.$m.pdf
		mv saf$i.slip.png saf$i.slip.$m.png
		mv saf$i.slip.ps  saf$i.slip.$m.ps
		cp saf$i.slip.$m.p* ./figure
		mv saf$i.slip.$m.p* ./res.saf$i
	fi
#end for1
done
#
#
	mv slip.m ./res.saf$i/
	rm  slip.f*
	mv $i.fault.gmt ./res.saf$i
	mv ores.saf$i.tgz ./res.saf$i
	rm  res.saf$i.tgz #clean first
#	tar -czvf res.saf$i.tgz res.saf$i
	#rm -fr res.saf$i cgrid.m.txt
	#rm -fr res.saf$i cgrid.m.txt
	rm  cgrid.m.txt
#endfor2
done<templist

#rm -fr figureagu1uplift #clean first
#rm -fr figureagu1uplift.tgz #clean first
#mv figure figureagu1uplift
#tar -czvf figureagu1uplift.tgz figureagu1uplift
#rm pp123.f pp123 saf716loop.sh 
#rm saf716loop.velo.sh 
#rm saf716loop.e.sh 
#rm saf716loop.e-slip.sh 
#rm saf716loop.j2.sh 
#rm saf716loop.slip.sh 
#rm saf716loop.uplift.sh 
#rm saf716loop.pauto.sh
#rm saf716loop.j2auto.sh 


#rm -fr figure
