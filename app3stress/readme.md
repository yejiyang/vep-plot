# README for app3stress, including beachball and 2-D principal stress direction plot
a sucessful example, see ./ok_app31/

## 1. Beach ball plot

There are three ways to plot beach ball using model stress results. Hui and Yang use fortran code, Sun use Matlab code.


### a. Hui beach ball

Check app3stress/ori_results/app32/p12345.f, and the output file grid.txt

The plot bash file, ./vep-beachball.sh.
```bash
	awk '{if(($1==1)&& int(($2+5000)*10/100000)==(($2+5000)*10/100000) && int($3*10/100000)==($3*10/100000))  print $2/1000,$3/1000,10,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20}' $input > bb.txt
	psmeca bb.txt  -R -J -Sx0.09i -Ggrey -Ewhite  -O  -W2  -C -V -K>>$psfile
```

### b. Yang beach ball

Check app3stress/ori_results/app32/p12345.f, and the output file beachball.txt
The plot bash file, ./vep-beachball.sh.
```bash
	awk '{if(int(($1+5000)*10/100000)==(($1+5000)*10/100000) && int($2*10/100000)==($2*10/100000))  print $1/1000,$2/1000,$3,$4,$5,$6,$7,$8,$9,$10,$11/1000,$12/1000+5}' beachball.txt > yang.txt
	psmeca yang.txt -R  -J  -C -O -V  -Sd0.15i/1 -G255/0/0 >> $psfile
```

## 2. 2-D principal stress

Check the function principal_2d.m and test_principal_2d.m

Plot file: ./vep-beachball.sh
```bash
	# about the direction of 2D principal stress
	# our matlab output, extension--positive, compression--negative, but the sita equals to I1 counter-clockwise to x-axis.
 	# in gmt, psvelo, -Sx option, extension and compression definition the same. But the angle defition is different. 
	# the angle in gmt, is azimuth of II2 in degrees CW (clockwise) from North.
	# So we need transform the matlab output sita to GMT angle. I1=sita, I2=sita+90, angle=90-sita=90-sita-90=-sita, this's why -$9
	awk '{if(int(($2+5)*10/100)==(($2+5)*10/100) && int($3*10/100)==($3*10/100))  print $2,$3,$7,$8,-$9}' $inputprinstr > prin.txt
	# -A default -A0.03/0.12/0.09
	psvelo prin.txt -A0.01/0.6/0.1 -W0.25p,blue  -R  -J   -O -V  -Sx0.06  -K>> $psfile
```

