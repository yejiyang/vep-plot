#!/bin/sh

ifort -o p12345 p12345.f
./p12345
mv grid.txt sagrid.31.txt
mv beachball.txt beachball.31.txt
cp sagrid.31.txt ../../input/
cp beachball.31.txt ../../input/
