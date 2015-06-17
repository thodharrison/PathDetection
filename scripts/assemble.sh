for f1 in $1*_R1_*.fastq
do
    f2="${f1/_R1_/_R2_}"
	echo "Merging the reads ... ${f1}    ${f2}"
	fa="${f1/.fastq/.fasta}"
	echo "fasta file ${fa}"
	bin/fq2fa --merge --filter $f1 $f2 $fa 
	fname=`basename ${f1}` 
	direct=results/assembly/viral/"${fname/.fastq/ }"
    echo "Assembling sample ${direct}"
	mkdir -p $direct
	bin/idba_ud -o $direct -r $fa
done

