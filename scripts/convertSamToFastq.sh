#!/bin/sh

sam=$1
fqF=$2
fqR=$3

cat $sam | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > $2 
cat $sam | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > $3


