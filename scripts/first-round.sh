i = 0
for f1 in $1*_R1_*
do 
	i=$(($i+1))
	f2="${f1/R1/R2}"
	echo "Doing sample ${f1}"
	echo ${f1}---------${f2}
	python scripts/pathogen_detection.py -p/s/fir/a/nobackup/data/VIRAL_DB/pathogens.fa -rresults/tmp/refs/felCat5.fa -f${f1} -b${f2} -oresults
	echo "----------------------- file name ----------------------- ${i}.cov"
	python scripts/reporter.py -p/s/fir/a/nobackup/data/VIRAL_DB/pathogens.fa -rresults/tmp/refs/felCat5.fa -f${f1} -b${f2} -oresults > results/report/${i}.cov
done
