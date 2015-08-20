#!/usr/bin/env python

import sys
from optparse import OptionParser
from collections import Counter
from csv import reader, writer
import re

import pandas as pd

import pysam

parser = OptionParser()

if __name__ == '__main__' :

  opts, args = parser.parse_args()

  num_reads = 0
  bpos = Counter()

  bam_f = pysam.Samfile(args[0],'rb')

  unique_cigars = set()

  wtf = None
  for r in bam_f :
    num_reads += 1
    for rp, gp in r.aligned_pairs :
      if rp is not None and gp is not None :
        bpos[gp] += 1

  cnts = pd.DataFrame(list(bpos.items()),columns=['pos','count'])

  htt = reader(open('HTT_full.gff'),delimiter='\t')

  out_f = writer(sys.stdout,delimiter="\t")
  for r in htt :
    chrm, src, typ, st, en, name, strand, phase, attr = r
    st, en = map(int,(st,en))

    bin_cnts = cnts[(cnts.pos>=st) & (cnts.pos<en)].sum()

    out_f.writerow([chrm,st,en,name, bin_cnts['count']])
