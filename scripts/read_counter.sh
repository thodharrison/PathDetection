for f1 in $1*_R1_*
do 
	i=$(($i+1))
	echo ${f1}
	f2="${f1/R1/R2}"
	size=$(wc -l < $f1)
	#expr $size / 4
	echo ${size}
done
