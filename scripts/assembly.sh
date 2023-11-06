

for f in 02.QC/*R1_trim.fq.gz; do
    n=`echo ${f##*/} | cut -d _ -f 1-2` &&
    megahit -1 $f -2 $n'_R2_trim.fq.gz --min-contig-len $MIN_CONTIG_SIZE -m 0.85 -o 03.ASSEMBLY/ -t $NUM_THREADS
done
