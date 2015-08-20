for bamfn in HTT_alignments/*.bam
do
  BASE=`basename $bamfn .bam`
  python3 SeqSpliceView.py -b $BASE.bed -r ${BASE}_reads.bam -s -m 1000 -c 10 $bamfn HTT.gff > $BASE.svg
  samtools index $bamfn
done
