#!/bin/bash


source /home/pulodogato/venvs/ntksamc/bin/activate

num=$(awk 'BEGIN{for(i=2;i<=12;i+=0.5) printf "%.1f ", i}')

#for ph in 2.0 3.0 4.0
for ph in $num
#0.58 1.08 1.68 2.22 2.71 3.31 3.85 4.28 4.95 5.46 6.04 6.72 7.12 8.25 9.25 10.45 ## $num
#for ph in 0.65 1.63 2.18 2.76 3.31 3.86 4.39 4.95 5.5 6.04 6.61 7.31 8.25 9.34 10.45
#for ph in 2.5 4.5 7.5 10.5
do
/home/pulodogato/venvs/ntksamc/bin/python3.6 ./tksamc2-python3.py -f $1'.pdb' -s $2 -ph $ph >> output.txt
G=`grep 13-Total_dG= Output_''$2''_''$1''_pH_''$ph''_T_300.0.dat | awk '{print $NF}'`
echo $ph $G >> Gqq-ph_$1''.dat
done


python2 ./plot_Gqq-ph.py -f Gqq-ph_$1''.dat >> output_plot.txt 2> output_2plot.txt


zip ${3}.zip *.jpg *.dat *.pdb *. README.txt

find . -type l -delete
rm tksamc.exe

touch finished
