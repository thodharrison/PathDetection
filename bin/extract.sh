samtools view $1 | awk '{OFS="\t"; print ">"$1"\n"$10}' - > $2

