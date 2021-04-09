#!/bin/bash

num=$(awk 'BEGIN{for(i=2;i<=12;i+=0.5) printf "%.1f ", i}')

for ph in $num
do
/usr/bin/python3 ./tksamc2-python3.py -f $1'.pdb' -s $2 -ph $ph >> output.txt
G=`grep 13-Total_dG= Output_''$2''_''$1''_pH_''$ph''_T_300.0.dat | awk '{print $NF}'`
echo $ph $G >> Gqq-ph_$1''.dat
done

/usr/bin/python3 ./plot_Gqq-ph.py -f Gqq-ph_$1''.dat >> output_plot.txt 2> output_2plot.txt


zip ${3}.zip *.jpg *.dat *.pdb *. README.txt

find . -type l -delete
rm tksamc.exe

touch finished
