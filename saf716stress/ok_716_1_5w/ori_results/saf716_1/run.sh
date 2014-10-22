#!/bin/sh

# p12345.f for salton sea app32, output both hui and yang's beach ball
# p12345.f.ok1.clean salton sea app32, only output hui's beach ball

# pp123.f for saf716 model, output p, e, no beach ball
# pp1234.f , for saf716 model, output both hui and yang's beach ball, in long and lat
# pp12345.f, for saf716 model, output both hui and yang's beach ball, in x, y

#ifort -o p12345 p12345.f
#./p12345
#mv grid.txt sagrid.31.txt
#mv beachball.txt beachball.31.txt
#cp sagrid.31.txt ../../input/
#cp beachball.31.txt ../../input/

ifort -o pp1234 pp1234.f
./pp1234
mv cgrid.m.txt sagrid.1_5w.txt
mv beachball.txt beachball.1_5w.txt
cp sagrid.1_5w.txt ../../input/
cp beachball.1_5w.txt ../../input/
