#!/bin/bash

/usr/bin/python3 ./tksamc2-python3.py -T $1 -ph $2 -s MC -f $3 > output.txt

sed '/#/Id' Output*.dat | awk '{print $3" "$12}' > dG_Energy_`expr "$3" : "processed_\(.*\).pdb"`.dat

zip ${4}.zip *.jpg *.dat *.pdb *. README.txt

find . -type l -delete
rm tksamc.exe

touch finished
