#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *31.01.2012

"""

import os,sys
import io
import gzip
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-o", "--newpath", dest="newpath", help="PATH for the split EMF and EMF index files (def ./split/)",default="./split/")
(options, args) = parser.parse_args()

if options.newpath == "": options.newpath = "./"
elif not os.path.isdir(options.newpath):
  options.newpath = options.newpath.rstrip('/')+"/"
  os.makedirs(options.newpath)

translations = {'Homo_sapiens': 'Hsap',
                'Pan_troglodytes': 'Ptro',
                'Gorilla_gorilla': 'Ggor',
                'Pongo_pygmaeus' : 'Ppyg',
                'Macaca_mulatta' : 'Mmul',
                'Callithrix_jacchus' : 'Cjac',
                'Ancestral_sequences' : 'Aseq'}

repeated = set()
outfile_hsa = open(options.newpath+"hsa_emf.index",'w')
outfile_hsa.write("#Chr\tStart\tEnd\tFile\tByte\tLine\n")
outfile_ptr = open(options.newpath+"ptr_emf.index",'w')
outfile_ptr.write("#Chr\tStart\tEnd\tFile\tByte\tLine\n")

print("Reading EMF files...")
for filename in args:
  fcount = 1
  if os.path.isfile(filename) and filename.endswith('.gz'):
    print (filename)
    f = gzip.open(filename,'rb')
    file_content = f.read().splitlines()
    f.close()
  #elif os.path.exists(filename):
    #print filename
    #file_content = open(filename)
  else:
    continue

  start_data = False
  start_seq = False
  start_offset = None
  start_loffset = None
  skip = False
  seq_order = []
  humans = []
  chimps = []
  newfilename = filename.replace(".emf",".%d.emf"%(fcount))
  outfile = gzip.open(options.newpath+newfilename, 'wb')
  offset = 0
  lcount = 0
  content = ''

  for line in file_content:
    if line.startswith('#'): continue # FILE HEADER & COMMENTS
    else:
      #outfile.write(line)
      content += line+"\n"

    if line.startswith('SCORE'): pass
    elif line.startswith('//'): # END OF ENTRY
      if start_data and not skip:
        for elem in humans: outfile_hsa.write('\t'.join(map(str,list(seq_order[elem][1:4])+[newfilename,start_offset,start_loffset]))+"\n")
        for elem in chimps: outfile_ptr.write('\t'.join(map(str,list(seq_order[elem][1:4])+[newfilename,start_offset,start_loffset]))+"\n")
        if lcount > 2500000:
          outfile.write(content)
          content = ''
          outfile.close()
          fcount += 1
          newfilename = filename.replace(".emf",".%d.emf"%(fcount))
          outfile = gzip.open(options.newpath+newfilename, 'wb')
          offset = -len(line)
          lcount = -1

      start_data = False
      start_seq = False
      start_offset = None
      start_loffset = None
      skip = False
      seq_order = []
      humans = []
      chimps = []
    elif line.startswith('SEQ'):
      if not start_seq:
        start_offset = offset
        start_loffset = lcount
        start_seq = True

      fields = line.split()
      if fields[2] in translations: fields[2] = translations[fields[2]]
      # species, chromosome, start, end, strand, length
      cid = (fields[1][:1].upper()+fields[1].split('_')[1][:3],fields[2],int(fields[3])-1,int(fields[4]),('-' if (fields[5][0] == '-') else '+'),int(fields[6].split('=')[1][:-1]))

      if cid[0] == 'Hsap': humans.append(len(seq_order))
      if cid[0] == 'Ptro': chimps.append(len(seq_order))

      seq_order.append(cid)
      if '_'.join(map(str,cid[:-2])) in repeated: skip = True
      else: repeated.add('_'.join(map(str,seq_order[0][:-2])))

    elif line.startswith('TREE'): pass
    elif line.startswith('DATA'): start_data = True
    offset += len(line)
    lcount += 1

  outfile.write(content)
  outfile.close()

outfile_hsa.close()
outfile_ptr.close()
print("Wrote: %shsa_emf.index"%options.newpath)
print("Wrote: %sptr_emf.index"%options.newpath)