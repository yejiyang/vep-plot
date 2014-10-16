C
C Include file for memory allocation in each master sub-program and 
C slave sub-programs, all the parameters add togeter to define the 
C whole memory usage of one computation. These parameters are list as
C following:
C
C maxaa--maximum dimension of double-precision data common block  
C maxia--maximum dimension of integer data common block
C maxrpoolm--maximum dimension of double-precision temp data and files
C            assigned to one-dimension memory common block in master node
C maxipoolm--maximum dimension of integer temp data and files assigned 
C            to one-dimension memory common block in master node
C maxrpools--maximum dimension of double-precision temp data and files
C            assigned to one-dimension memory common block in slave mode
C maxipools--maximum dimension of integer temp data and files assigned
C            to one-dimension memory common block in slave node
C
        parameter (maxaa =     35000000)
	parameter (maxia =     85000000)
	parameter (maxrpoolm = 25000000)
	parameter (maxipoolm = 25000000)
	parameter (maxrpools = 25000000)
	parameter (maxipools = 25000000)
C End of memory allocation data file	
