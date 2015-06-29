#!/bin/bash

for ((num=1285; num<=1350 ; num = num + 1))
do
#copy file from a directory to this
cd Input/
cp Ion2.in ../Ion2.in 
cp Ion1.restart ../Ion1.restart
cd ..

#change certain parameter in data and in file 
sed -i -re "s/111/$num/"   Ion2.in

#run lammps
mpirun -np 1 ../lammps-26Jan2014/src/lmp_openmpi < Ion2.in

#grab data
line_num=125
string1=$(cat log.lammps | sed ${line_num}\!d)
line_num=165
string2=$(cat log.lammps | sed ${line_num}\!d)

string=$( echo "scale = 5;($string2 - $string1)/332.06371" | bc )
echo "$num      $string" >> test.txt
rm Ion* log*
done