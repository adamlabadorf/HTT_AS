# given a gff formatted file of exons in a gene and a BAM file, create an SVG
# image showing all splicing events in proportion to the number of reads
# supporting each intronic splice pattern

# parse gff file, get bounds for reads
# extract reads from BAM according to gff annotation
# generate SVG

import os
import sys

from collections import Counter, defaultdict
from csv import reader, writer
from optparse import OptionParser

from math import log10

cigar_d = {
'M': 0, # BAM_CMATCH
'I': 1, # BAM_CINS
'D': 2, # BAM_CDEL
'N': 3, # BAM_CREF_SKIP
'S': 4, # BAM_CSOFT_CLIP
'H': 5, # BAM_CHARD_CLIP
'P': 6, # BAM_CPAD
'=': 7, # BAM_CEQUAL
'X': 8, # BAM_CDIFF
}

import pysam

svg_props = {'width':3200,
             'model_margin': {
                'top': 20
             }
            }
svg_width = 3200
tmpl = """\
<?xml version="1.0"?>
<!DOCTYPE svg PUBLIC '-//W3C//DTD SVG 1.0//EN'
          'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'>
<svg width="%(width)d" height="%(height)d" xmlns:xlink="http://www.w3.org/1999/xlink"
     style="fill-opacity:1; stroke:black; stroke-linecap:square; stroke-miterlimit:10; stroke-opacity:1; shape-rendering:auto; fill:black; stroke-dasharray:none; font-weight:normal; stroke-width:0; font-family:&apos;Dialog&apos;; font-style:normal; stroke-linejoin:miter; font-size:12;" xmlns="http://www.w3.org/2000/svg">
<!--
<g id="grid" style="padding:10px;fill:rgb(200,200,255);opacity:0.8;" transform="translate(0,0)">
%(grid)s
</g>
-->
<g id="model" style="padding:10px;fill:rgb(200,200,255)" transform="translate(0,0)">
%(model)s
</g>
<g id="spliced_reads" style="padding:10px;fill:rgb(100,100,178)" transform="translate(0,20)">
%(rects)s
</g>
</svg>
"""

rect_tmpl = """\
<g style="text-align:center;font-size:%(fontsize)fpt;text-anchor:middle;">
  <rect x="%(x)f" y="%(y)f" width="%(w)f" height="%(h)f"/>
  <text x="%(text_x)d" y="%(text_y)d" fill="rgb(50,50,50)">%(c)d</text>
</g>
"""
def rect(**kwargs) :
  kwargs['text_x'] = x+w/2
  kwargs['text_y'] = y+h/2
  kwargs['fontsize'] = min(12,min(w,h)*.5)
  return rect_tmpl%kwargs

grid_tmpl = """\
<g style="font-size:9pt;">
  <rect x="%(x)d" y="0" width="2" height="%(h)d" fill="rgb(90,90,90)"/>
  <text x="%(x)d" y="0" fill="rgb(50,50,50)" style="writing-mode:tb">%(c)d</text>
</g>
"""
def grid(**kwargs):
  return grid_tmpl%kwargs

parser = OptionParser()
parser.add_option('-m','--margin',type='int',default=300,
    help="margin around gene boundaries to extract reads for")
parser.add_option('-f','--flag',type='int',default=255,
    help="require these bits to be set on each read's flag, e.g. 4 == secondary only")
parser.add_option('-F','--flag-not',type='int',default=0,
    help=("require these bits to be NOT set on each read's flag, "
          "e.g. 4 == primary only"))
parser.add_option('-s','--strict',action='store_true',default=False,
    help=("include reads only if their alignments are strictly within the gene "
          "bounds, otherwise any partially overlapping read will be included"))
parser.add_option('-c','--min-count',type='int',default=10,
    help=("exclude splice events supported by fewer than this number of reads"))
parser.add_option('-w','--min-width',type='int',default=9,
    help=("exclude splice events shorter than this number of bases"))
parser.add_option('-b','--bed-fn',default=None,help="write out introns to bed file")
parser.add_option('-r','--read-fn',default="splice_reads.bam",
    help="write out spliced reads for each intron to bam file [default: %default]")

if __name__ == '__main__' :

  opts, args = parser.parse_args()

  bam_f = pysam.Samfile(args[0],'rb')
  gff_fn = args[1]
  gff_r = reader(open(gff_fn),delimiter="\t")

  gff_coords = []
  gff_chrms = set()
  for r in gff_r :
    gff_chrms.add(r[0])
    gff_coords.append(tuple(map(int,r[3:5])))

  if len(gff_chrms) != 1 :
    parser.error("The GFF file contains records from more than one chromosome."
        " Ensure the annotation is correct.")

  gene_chrm = gff_chrms.pop()
  gene_coords = (min(gff_coords)[0],max(_[1] for _ in gff_coords))
  read_window = gene_coords[0]-opts.margin, gene_coords[1]+opts.margin

  unique_cigars = set()

  juncs = Counter()
  junc_reads = defaultdict(list)

  read_iter = bam_f.fetch(gene_chrm,*gene_coords)
  for r in read_iter :

    # filter based on flag opts
    if (opts.flag & r.flag)==0 or (opts.flag_not & r.flag)!=0 :
        #print('filtering read based on flag: %d'%r.flag,file=sys.stderr)
        continue

    # filter based on strict coordinates
    if opts.strict and \
        (r.reference_start < gene_coords[0] or \
         r.reference_end > gene_coords[1] ) :
        #print('filtering read based on coords: %s, %s'%(gene_coords,
        #    (r.reference_start,r.reference_end)),file=sys.stderr)
        continue
        
    ct = r.cigartuples
    splice_coords = [None,None]
    for i in range(len(ct)-1) :
      c1,c2 = ct[i],ct[i+1]
      run_sum = sum(_[1] for _ in ct[:i] if _[0] not in (cigar_d['I'],cigar_d['S']))
      if c1[0] == cigar_d['M'] and c2[0] == cigar_d['N'] :
        splice_coords[0] = r.reference_start+c1[1]+run_sum
      elif c1[0] == cigar_d['N'] and c2[0] == cigar_d['M'] :
        splice_coords[1] = r.reference_start+c1[1]+run_sum
      if None not in splice_coords :
        juncs[tuple(splice_coords)] += 1
        junc_reads[tuple(splice_coords)].append(r)
        splice_coords = [None,None]

  if opts.bed_fn :
    fout = writer(open(opts.bed_fn,'w'),delimiter="\t")
    fout.writerows([[gene_chrm,_[0],_[1],'intron_%d'%i,juncs[_]] for i,_ in enumerate(sorted(juncs))])

  coords = list(juncs.keys())
  #reads_out = open(opts.read_fn,'w')
  reads_out = pysam.AlignmentFile(opts.read_fn+"_tmp",'wb',template=bam_f)
  intron_i = 0
  for k in sorted(coords) :
    count = juncs[k]
    w = abs(k[1]-k[0])
    if count < opts.min_count or w <= opts.min_width :
      del juncs[k]
      del junc_reads[k]
    else :
      #reads_out.write('#intron_%d %s:%d-%d %d\n'%(intron_i,gene_chrm,k[0],k[1],juncs[k]))
      #reads_out.write("".join([str(_)+"\n" for _ in junc_reads[k]]))
      for _ in junc_reads[k] :
          reads_out.write(_)
      intron_i += 1
  reads_out.close()

  # sort the stupid bam
  pysam.sort("-f",opts.read_fn+"_tmp",opts.read_fn)
  os.remove(opts.read_fn+"_tmp")

  min_start, max_start = gene_coords
  for k in juncs.keys() :
    min_start = min(min_start,*k)
    max_start = max(max_start,*k)

  rects = []
  curr_height = 0
  max_height = 0
  coord_stack = [(min(juncs),0)]
  for coords in sorted(juncs) :
    count = juncs[coords]
    w = (abs(coords[1]-coords[0]))/(max_start-min_start)*svg_width
    x = (coords[0]-min_start)/(max_start-min_start)*svg_width

    while len(coord_stack) > 0 and coords[0]>coord_stack[-1][0][1] :
      coord_stack.pop()

    if len(coord_stack) == 0 :
      curr_height = 0
    else :
      curr_height = coord_stack[-1][1]

    y = curr_height+2
    h = 3+count/10
    curr_height += h+2

    coord_stack.append((coords,curr_height))

    max_height = max(curr_height,max_height)

    rects.append(rect(x=x,y=y,w=w,h=h,c=count))

  
  # do the model rects
  gff_rects = []
  for i,(st,en) in enumerate(gff_coords) :
    w = max(1,(en-st)/(max_start-min_start)*svg_width)
    x = (st-min_start)/(max_start-min_start)*svg_width
    h = max_height
    gff_rects.append(rect(x=x,y=0,w=w,h=h,c=i+1))

  # grid
  grid_elems = []
  num_lines = 30
  for i in range(num_lines) :
    grid_elems.append(grid(x=(svg_width*i/num_lines),h=max_height,c=(min_start+(max_start-min_start)*i/num_lines)))

  print(tmpl%{
    "rects":"".join(rects),
    "model":"".join(gff_rects),
    "grid":"".join(grid_elems),
    "width":svg_props['width']+w,
    "height":max_height+svg_props['model_margin']['top']})

