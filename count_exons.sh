for fn in HTT_alignments/*merged*HTT.bam
do
    BASE=`basename $fn .bam`
    echo $fn
    #htseq-count -i exon_id -r pos -f bam -s no -t exon -m intersection-nonempty $fn HTT.gff > counts/${BASE}_exon_counts.txt
    python count_exons.py $fn > counts/${BASE}_exon_counts.txt
done

#for fn in HTT_alignments/ml*HTT.bam
#do
#    BASE=`basename $fn .bam`
#    echo $fn
#    #htseq-count -i exon_id -r pos -f bam -s no -t exon -m intersection-nonempty $fn HTT.gff > counts/${BASE}_exon_counts.txt
#    python count_exons.py $fn > counts/${BASE}_exon_counts.txt
#done
