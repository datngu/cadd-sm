#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *08.05.2012

"""

import os,sys
import string
import gzip
import shutil
from optparse import OptionParser
from bx.intervals.intersection import Intersecter, Interval
from collections import defaultdict
import traceback

# REVERSE COMPLEMENT DNA
table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv')
def revcompl(seq):
  global table
  return seq.translate(table)[::-1]

# COMPLEMENT DNA
def compl(seq):
  global table
  return seq.translate(table)

# SHORTEN SEQUENCE (USEFUL FOR ABBREVIATING LONG INDELS)
def shorten_seq(seq):
  if len(seq) > 12: 
    sub_seq = seq.split(",")
    if len(sub_seq) == 1:
      return "%s..%d..%s"%(seq[:4],len(seq)-8,seq[-4:])
    else:
      res = []
      for seq in sub_seq:
        if len(seq) > 12: res.append("%s..%d..%s"%(seq[:4],len(seq)-8,seq[-4:]))
        else: res.append(seq)
      return ",".join(res)
  else: return seq

# READ FASTA INDEX TO MEMORY
def read_fasta_index(filename):
  infile = open(filename)
  res = {}
  for line in infile:
    fields = line.split()
    if len(fields) == 5:
      cname,length,start,line,cline = fields
      res[cname]=int(length),int(start),int(line),int(cline)
    else:
      sys.stderr.write('Error: Unexpected line in fasta index file: %s\n'%(line.strip()))
      sys.exit()
  infile.close()
  return res

# TRANSLATION OF EPO SPECIES IDENTIFIERS
translations = {'Homo_sapiens': 'Hsap',
                'Pan_troglodytes': 'Ptro',
                'Gorilla_gorilla': 'Ggor',
                'Pongo_pygmaeus' : 'Ppyg',
                'Macaca_mulatta' : 'Mmul',
                'Callithrix_jacchus' : 'Cjac',
                'Ancestral_sequences' : 'Aseq'}

# RETURN A DICTIONARY (CHROM) OF DICTIONARIES (POS) WITH INFORMATION ABOUT FILE, BYTES AND LINES
def init_index(filename):
  if os.path.isfile(filename):
    infile = open(filename)
    res = {}
    for line in infile:
      if line.startswith("#"): continue
      fields = line.split() #Chr\tStart\tEnd\tFile\tByte\tLine
      if len(fields) == 6:
        chrom,start,end,cfile,bytes,lines = fields
        start,end = int(start)+1,int(end)+1 #MAKE COORDINATES ONE BASED HALF-OPEN
        if chrom not in res: 
          res[chrom] = {}
          res[chrom]["coords"] = Intersecter()
        res[chrom]["coords"].add_interval(Interval(start, end))
        res[chrom][(start,end)] = cfile,int(bytes),int(lines)
      else:
        sys.stderr.write("Unexpected line [%s]: %s"%(filename,line))
    infile.close()
    return res
  return None

# RETURN EMF BLOCK TYPE
def emf_get_type():
  global current
  if current["btype"] == None: return "NA"
  else: return current["btype"]

def emf_get_base_at_position(species,chrom,pos,outspecies="Hsap-Ptro",silent=False):
  global current
  global hsa_index
  global ptr_index
  global translations
  if species == "Hsap": cindex = hsa_index
  elif species == "Ptro": cindex = ptr_index
  else: return None

  if chrom.startswith('chr'): chrom = chrom[3:]
  if chrom.startswith('gl'): chrom = chrom.upper()
  if chrom == "M": chrom = "MT"
  if chrom == "23": chrom = "X"
  result = []

  if current['species'] != species or current['chrom'] != chrom or (pos < current['bstart']) or (pos >= current['bend']):
    if chrom in cindex:
      res = cindex[chrom]["coords"].find(pos,pos+1)
      #print chrom,pos,res
      if len(res) == 1:
        rstart = res[0].start
        rend = res[0].end
        filename,bytes,loffset = cindex[chrom][(rstart,rend)]
        #bytes,loffset=int(bytes),int(loffset)
        if current['filename'] != filename:
          sys.stderr.write("Reading file: %s\n"%(options.ancestor_path+"/"+filename))
          if os.path.isfile(options.ancestor_path+"/"+filename) and filename.endswith('.gz'):
            f = gzip.open(options.ancestor_path+"/"+filename,'rb')
            del current['data']
            current['data'] = f.read().splitlines()
            f.close()
            current['filename'] = filename
          else:
            sys.stderr.write("Error: EMF input file (%s) is not available\n"%(filename))
        current['species'] = species
        current['chrom'] = chrom
        current['btype'] = None
        start_data = False
        is_compl = False
        found_pos = False
        current['block'] = []
        current['bline'] = None
        current['bstart'] = rstart
        current['bend'] = rend
        current['bseqs'] = []
        current['bpos'] = []
        current['bstrands'] = []
        current['bgaps'] = []
        current['bancestors'] = {}
        current['bref'] = None
        sys.stderr.write("Initiating block (%s %d-%d)...\n"%(chrom,rstart,rend))
        for line in current['data'][loffset:]:
          if line.startswith('#'): continue # HEADER
          elif line.startswith('SCORE'): continue
          elif line.startswith('//'): # END OF ENTRY
            break
            #sys.stderr.write("Unexpected error: hit block end...\n")
          elif line.startswith('SEQ'):
            fields = line.split()
            if fields[2] in translations: fields[2] = translations[fields[2]]
            # species, chromosome, start, end, strand, length
            cid = (fields[1][:1].upper()+fields[1].split('_')[1][:3],fields[2],int(fields[3])-1,int(fields[4]),('-' if (fields[5][0] == '-') else '+'),int(fields[6].split('=')[1][:-1]))
            if species == cid[0] and chrom == cid[1] and rstart == cid[2]+1 and rend == cid[3]+1:
              current['bref'] = len(current['bseqs'])
              current['bline'] = 0
              is_compl = (fields[5][0]=='-')
            current['bstrands'].append(fields[5][0]!='-')
            current['bgaps'].append(0)
            if current['bstrands'][-1]: current['bpos'].append(int(fields[3])-1)
            else: current['bpos'].append(int(fields[4])+1)
            current['bseqs'].append(cid)

          elif line.startswith('TREE'):
            refname = "%s_%s_%d_%d[%s]"%(current['bseqs'][current['bref']][0],current['bseqs'][current['bref']][1],current['bseqs'][current['bref']][2]+1,current['bseqs'][current['bref']][3],current['bseqs'][current['bref']][4])
            tree = line
            keller = []
            hname = ''
            hpos = 0
            while hpos < len(tree):
              elem = tree[hpos]
              if elem == '(':
                if hname != '': 
                  keller.append(hname)
                  hname = ''
                keller.append(elem)
              elif elem == ',':
                if hname != '':
                  keller.append(hname)
                  hname = ''
              elif elem == ')':
                if hname != '':
                  keller.append(hname)
                  hname = ''
                ancestor = ''
                hpos += 1
                while (hpos < len(tree)) and (tree[hpos] != ',')  and (tree[hpos] != ')'):
                  elem = tree[hpos]
                  ancestor += elem
                  hpos += 1
                ancestor = ancestor.split(':')[0]
                if (hpos < len(tree)): hpos -= 1
                #print '--------------'
                #print "K-I",keller
                #print 'Ancestor',ancestor
                value = keller.pop()
                reorder = []
                while value != '(':
                  if (value.split(':')[0] == refname): value = '!'+value.replace(species,"SELF")
                  elif value.startswith('SELF'): value = '#'+value
                  elif value.startswith(species): value = '$'+value
                  else: value = '-'+value
                  reorder.append(value.split(':')[0])
                  value = keller.pop()
                reorder.sort()
                val1 = reorder[0].split('_')[0][1:].split('-')[0]
                val2 = reorder[-1].split('_')[0][1:].split('-')[0]
                if val1 == val2: save = val1+'_'+'_'.join(ancestor.split('_')[1:][:-2])
                else: save = val1+'-'+val2+'_'+'_'.join(ancestor.split('_')[1:][:-2])
                keller.append(save)
                current['bancestors']['_'.join(save.split('_')[1:])] = save.split('_')[0]
                #print "K-O",keller
              else:
                hname += elem
              hpos += 1
          elif line.startswith('DATA'):
            start_data = True
            current['btype'] = ''.join(map(lambda x: x[0][0] if not x[0].startswith('Pabe') else 'O', filter(lambda x:not x[0].startswith('Aseq'),current['bseqs'])))
            #print current['bancestors']
            #print current['bseqs']
            for hpos,elem in enumerate(current['bseqs']):
              if elem[1] in current['bancestors']:
                helper = list(elem)
                helper[0] = current['bancestors'][elem[1]]
                current['bseqs'][hpos] = tuple(helper)
            #print current['bseqs']
          elif start_data and len(line) >= len(current['bseqs']):
            line = "".join(line.split(' ')) # FORMAT SPECIFICATION ALLOWS FOR SPACES
            if is_compl: line = compl(line)
            current['block'].append(line)
            # COUNT BASES:
            if not found_pos:
              for hpos,cid in enumerate(current['bseqs']):
                if line[hpos] == "-": current['bgaps'][hpos]+=1
                elif line[hpos] == "~": continue
                else:
                  current['bgaps'][hpos] = 0
                  if current['bstrands'][hpos]: current['bpos'][hpos]+=1
                  else: current['bpos'][hpos]-=1
              if current['bpos'][current['bref']] == pos:
                found_pos = True
                current['bline'] = len(current['block'])-1
                #hhcount = 0
            #if found_pos and hhcount < 20:
              #print line.strip()
              #hhcount += 1
      elif len(res) == 0:
        if not silent: sys.stderr.write("Error: Coordinate (%s:%d) not available from EMF alignments\n"%(chrom,pos))
      else:
        sys.stderr.write("Error: Coordinate (%s:%d) seems present in several EMF alignments\n"%(chrom,pos))
    else:
      if not silent: sys.stderr.write("Error: Chromosome (%s) not available from EMF alignments\n"%(chrom))

  if (len(result)==0) and (current['species'] == species) and (current['chrom'] == chrom) and (pos >= current['bstart']) and  (pos < current['bend']):
    # GO TO POSITION
    # FORWARD
    if (current['bstrands'][current['bref']] and current['bpos'][current['bref']] <= pos) or (not current['bstrands'][current['bref']] and current['bpos'][current['bref']] >= pos):
      while (current['bline'] < (len(current['block'])-1)) and (current['bpos'][current['bref']] != pos):
        current['bline']+=1
        line = current['block'][current['bline']]
        for hpos,cid in enumerate(current['bseqs']):
          if line[hpos] == "-": current['bgaps'][hpos]+=1
          elif line[hpos] == "~": continue
          else:
            current['bgaps'][hpos] = 0
            if current['bstrands'][hpos]: current['bpos'][hpos]+=1
            else: current['bpos'][hpos]-=1
    # REVERSE
    else:
      while (current['bline'] > 0) and (current['bpos'][current['bref']] != pos):
        current['bline']-=1
        line = current['block'][current['bline']]
        for hpos,cid in enumerate(current['bseqs']):
          if line[hpos] == "-": current['bgaps'][hpos]-=1
          elif line[hpos] == "~": continue
          else:
            current['bgaps'][hpos] = 0
            if current['bstrands'][hpos]: current['bpos'][hpos]-=1
            else: current['bpos'][hpos]+=1
    # ARE AT THE POSITION
    if current['bpos'][current['bref']] == pos:
      # SAVE CURRENT POSITION
      helper = (list(current['bgaps']),list(current['bpos']),current['bline'])
      if current['bstrands'][current['bref']]: # PLUS STRAND
        while (current['bline'] < len(current['block'])) and (current['bpos'][current['bref']] == pos): # WHILE WE ARE AT THE SAME REFERENCE POSITION WE COLLECT ALL BASES FROM OTHERS
          line = current['block'][current['bline']]
          if outspecies == "SELF":
            if len(result) == 0: result.append(line[current['bref']])
            else: result[0]+=line[current['bref']]
          elif outspecies == "ALL":
            if len(result) == 0:
              for hpos,cid in enumerate(current['bseqs']):
                if not cid[0].startswith('Aseq') and not cid[0].startswith('SELF-'): result.append(line[hpos])
            else:
              ind = 0
              for hpos,cid in enumerate(current['bseqs']):
                if not cid[0].startswith('Aseq') and not cid[0].startswith('SELF-'):
                  if line[hpos]!="-": result[ind]+=line[hpos]
                  ind += 1
          else:
            if len(result) == 0:
              for hpos,cid in enumerate(current['bseqs']):
                if cid[0] == outspecies: result.append(line[hpos])
            else:
              ind = 0
              for hpos,cid in enumerate(current['bseqs']):
                if cid[0] == outspecies:
                  if line[hpos]!="-": result[ind]+=line[hpos]
                  ind += 1
          current['bline']+=1
          if current['bline'] < len(current['block']):
            line = current['block'][current['bline']]
            for hpos,cid in enumerate(current['bseqs']):
              if line[hpos] == "-": current['bgaps'][hpos]+=1
              elif line[hpos] == "~": continue
              else:
                current['bgaps'][hpos] = 0
                if current['bstrands'][hpos]: current['bpos'][hpos]+=1
                else: current['bpos'][hpos]-=1
      else: # MINUS STRAND
        while (current['bline'] >= 0) and (current['bpos'][current['bref']] == pos): # WHILE WE ARE AT THE SAME REFERENCE POSITION WE COLLECT ALL BASES FROM OTHERS
          line = current['block'][current['bline']]
          if outspecies == "SELF":
            if len(result) == 0: result.append(line[current['bref']])
            else: result[0]+=line[current['bref']]
          elif outspecies == "ALL":
            if len(result) == 0:
              for hpos,cid in enumerate(current['bseqs']):
                if not cid[0].startswith('Aseq') and not cid[0].startswith('SELF-'): result.append(line[hpos])
            else:
              ind = 0
              for hpos,cid in enumerate(current['bseqs']):
                if not cid[0].startswith('Aseq') and not cid[0].startswith('SELF-'): 
                  if line[hpos]!="-": result[ind]+=line[hpos]
                  ind += 1
          else:
            if len(result) == 0:
              for hpos,cid in enumerate(current['bseqs']):
                if cid[0] == outspecies: result.append(line[hpos])
            else:
              ind = 0
              for hpos,cid in enumerate(current['bseqs']):
                if cid[0] == outspecies:
                  if line[hpos]!="-": result[ind]+=line[hpos]
                  ind += 1
          current['bline']-=1
          if (current['bline'] >= 0):
            line = current['block'][current['bline']]
            for hpos,cid in enumerate(current['bseqs']):
              if line[hpos] == "-": current['bgaps'][hpos]-=1
              elif line[hpos] == "~": continue
              else:
                current['bgaps'][hpos] = 0
                if current['bstrands'][hpos]: current['bpos'][hpos]-=1
                else: current['bpos'][hpos]+=1
      # SET BACK POSITION TO FIRST HIT
      current['bgaps'],current['bpos'],current['bline'] = helper
  return result

def format_ancestor(hlist,last=None,unique=True):
  if last != None:
    if len(hlist) == 0 and last == "N/A": return "N/A"
    elif len(hlist) == 0: return last
    elif len(hlist) == 1 and last == "N/A": return hlist[0].upper()
    elif len(hlist) == 1: return last+hlist[0].upper()
    elif len(hlist) > 1 and last == "N/A":
      if unique: return ",".join(list(set(map(lambda x:x.upper(),hlist))))
      else: return ",".join(map(lambda x:x.upper(),hlist))
    else:
      seqs = last.split(',')
      if len(seqs) == hlist: 
        hlist = map(lambda x:x.upper(),hlist)
      else:
        hlist = list(set(map(lambda x:x.upper(),hlist)))
        if len(seqs) != hlist: return "N/A"
      for elem,ind in seqs:
        hlist[ind] = elem+hlist[ind]
      return ",".join(hlist)
  else:
    if len(hlist) == 0: return "N/A"
    elif len(hlist) == 1: return hlist[0].upper()
    else:
      if unique: return ",".join(list(set(map(lambda x:x.upper(),hlist))))
      else: return ",".join(map(lambda x:x.upper(),hlist))


parser = OptionParser()
parser.add_option("-c","--chroms", dest="chromosomes", help="Chromosomes considered (e.g. '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22')",default='')
parser.add_option("-i", "--infile", dest="infile", help="Infile name of genome (prefix for fasta index, def /net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g/whole_genome.fa)",default="/net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g/whole_genome.fa")
parser.add_option("-o", "--outfile", dest="outfile", help="Name of output file (def STDOUT)",default=None)
parser.add_option("-a", "--ancestor_path", dest="ancestor_path", help="PATH of the EMF and EMF index files (def /net/shendure/vol1/home/mkircher/sequencedb/wga/epo_6_primate_v66/split_mod)",default="/net/shendure/vol1/home/mkircher/sequencedb/wga/epo_6_primate_v66/split_mod")
parser.add_option("-r", "--reference", dest="reference", help="Reference species (def Hsap)",default="Hsap")
(options, args) = parser.parse_args()

# OPEN OUTPUT FILE IF GIVEN, OTHERWISE PRINT ON SCREEN
if options.outfile != None:
  output = open(options.outfile,'w')
else:
  output = sys.stdout

# OUTPUT ERROR MESSAGES AND EXIST IF THERE ARE MISSING FILES
if not os.path.isdir(options.ancestor_path) or not os.path.exists(options.ancestor_path+"/hsa_emf.index") or not os.path.exists(options.ancestor_path+"/ptr_emf.index"):
  sys.stderr.write("Require valid path to EMF and EMF index files.\n")
  sys.exit()

if not os.path.exists(options.infile) or not os.path.exists(options.infile+".fai"):
  sys.stderr.write("Invalid path for genome and genome fasta index file.\n")
  sys.exit()

fastaindex = read_fasta_index(options.infile+".fai")

# MAKE AN INDEX OF HUMAN EMF FILES
hsa_index = init_index(options.ancestor_path+"/hsa_emf.index")
# MAKE AN INDEX OF CHIMP EMF FILES
ptr_index = init_index(options.ancestor_path+"/ptr_emf.index")
current = { "data":None, "filename":None, "species":None, 
            "chrom":None, "pos":None,
            "block": None, "bline":None, "bstart":None, "bend":None,
            "bseqs": None, "bpos": None, "bstrands": None, "bgaps": None,
            "bancestors": None, "bref":None, "btype":None }

if options.reference == "Hsap":
  CAnc = "SELF-Ptro"
else:
  CAnc = "SELF-Hsap"

if options.chromosomes != '':
  chroms = options.chromosomes.split(',')
else:
  chroms = filter(lambda x: x.isdigit() or x in ['X','Y'], fastaindex.keys())

try:
  for chrom in chroms:
    output.write("###CHROM %s\n"%chrom)
    length,sblock,bline,cline = fastaindex[chrom]

    totalrefA = 0
    totalrefC = 0
    totalrefG = 0
    totalrefT = 0
    total = 0
    mut = 0
    totalCpG = 0 # Total bases in CpG sites
    mutCpG = 0 # Total number of mutations CpG sites

    ACn = 0
    AGn = 0
    ATn = 0
    CAn = 0
    CGn = 0
    CTn = 0
    GAn = 0
    GCn = 0
    GTn = 0
    TAn = 0
    TCn = 0
    TGn = 0

    CA = 0
    CG = 0
    CT = 0
    GA = 0
    GC = 0
    GT = 0

    insertionsizes,deletionsizes = defaultdict(int),defaultdict(int)

    lref,lhp,lhg,lho = "N","N","N","N"
    ref,hp,hg,ho = "N","N","N","N"
    insertion = 0
    deletion = 0
    window_start = None
    window_lmut,window_ltotal,window_lmutCpG,window_ltotalCpG = 0,0,0,0
    window_lA,window_lC,window_lG,window_lT = 0,0,0,0
    count = 0
    for pos in xrange(1,length+1):
      lref,lhp,lhg,lho = ref,hp,hg,ho
      ref = format_ancestor(emf_get_base_at_position(options.reference,chrom,pos,"SELF",True))
      if ref != "N/A":

        hp = format_ancestor(emf_get_base_at_position(options.reference,chrom,pos,CAnc,True))
        hg = format_ancestor(emf_get_base_at_position(options.reference,chrom,pos,"SELF-Ggor",True))
        ho = format_ancestor(emf_get_base_at_position(options.reference,chrom,pos,"SELF-Pabe",True))
        nref = format_ancestor(emf_get_base_at_position(options.reference,chrom,pos+1,"SELF",True))
        nhp = format_ancestor(emf_get_base_at_position(options.reference,chrom,pos+1,CAnc,True))

        if (hp != "N/A") and ((hp == hg) or (hp == ho)) and (ref[0] in ["A","C","G","T"]) and (nref[0] in ["A","C","G","T"]) and (lref[0] in ["A","C","G","T"]):
          if window_start == None: window_start = pos
          if hp == "-" and ((insertion > 0) or (insertion == 0 and lhp[0] in ["A","C","G","T"])):
            insertion += 1
          elif hp[0] in ["A","C","G","T"]:
            if insertion > 0:
              insertionsizes[insertion]+=1
              insertion = 0
            if deletion > 0:
              deletionsizes[deletion]+=1
              deletion = 0
            if len(ref) > 1:
              trhp = hp[1:].replace("-","")
              if ref[1:].replace("-","") == "" and len(trhp) > 0:
                deletion = len(trhp)

            if hp[0] == "A": totalrefA+=1
            elif hp[0] == "C": totalrefC+=1
            elif hp[0] == "G": totalrefG+=1
            else: totalrefT+=1

            if (hp[0] == 'G' and lref == 'C' and lhp == 'C') or (hp[0] == 'C' and nref == 'G' and nhp == 'G'):
              totalCpG += 1
              if ref[0] != hp[0]:
                mutCpG += 1
                if hp[0] == 'C' and ref[0] == 'A': CA += 1
                elif hp[0] == 'C' and ref[0] == 'G': CG += 1
                elif hp[0] == 'C' and ref[0] == 'T': CT += 1
                elif hp[0] == 'G' and ref[0] == 'A': GA += 1
                elif hp[0] == 'G' and ref[0] == 'C': GC += 1
                elif hp[0] == 'G' and ref[0] == 'T': GT += 1
            else:
              total += 1
              if ref[0] != hp[0]:
                mut += 1
                if hp[0] == 'A' and ref[0] == 'C': ACn += 1
                elif hp[0] == 'A' and ref[0] == 'G': AGn += 1
                elif hp[0] == 'A' and ref[0] == 'T': ATn += 1
                elif hp[0] == 'C' and ref[0] == 'A': CAn += 1
                elif hp[0] == 'C' and ref[0] == 'G': CGn += 1
                elif hp[0] == 'C' and ref[0] == 'T': CTn += 1
                elif hp[0] == 'G' and ref[0] == 'A': GAn += 1
                elif hp[0] == 'G' and ref[0] == 'C': GCn += 1
                elif hp[0] == 'G' and ref[0] == 'T': GTn += 1
                elif hp[0] == 'T' and ref[0] == 'A': TAn += 1
                elif hp[0] == 'T' and ref[0] == 'C': TCn += 1
                elif hp[0] == 'T' and ref[0] == 'G': TGn += 1

            if total % 100000 == 0:
              cmut = mut-window_lmut
              ctotal = total-window_ltotal
              cmutCpG = mutCpG-window_lmutCpG
              ctotalCpG = totalCpG-window_ltotalCpG

              cA,cC,cG,cT = totalrefA-window_lA,totalrefC-window_lC,totalrefG-window_lG,totalrefT-window_lT
              output.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n"%(chrom,window_start,pos,cmut,ctotal,cmutCpG,ctotalCpG,cA,cC,cG,cT,0 if ctotal == 0 else cmut/float(ctotal),0 if ctotalCpG == 0 else cmutCpG/float(ctotalCpG),0 if ctotal == 0 else cA/float(ctotal),0 if ctotal == 0 else cC/float(ctotal),0 if ctotal == 0 else cG/float(ctotal),0 if ctotal == 0 else cT/float(ctotal)))
              window_start=pos+1
              window_lmut=mut
              window_ltotal=total
              window_lmutCpG=mutCpG
              window_ltotalCpG=totalCpG
              window_lA,window_lC,window_lG,window_lT = totalrefA,totalrefC,totalrefG,totalrefT

          #else:
            #print "Uups",chrom,pos,ref,hp,hg,ho
      else:
        insertion = 0
        deletion = 0
        ref,ho,hg,ho = "N","N","N","N"
    if total > 0:
      if (total-window_ltotal)/100000. > 0.1:
        cmut = mut-window_lmut
        ctotal = total-window_ltotal
        cmutCpG = mutCpG-window_lmutCpG
        ctotalCpG = totalCpG-window_ltotalCpG
        cA,cC,cG,cT = totalrefA-window_lA,totalrefC-window_lC,totalrefG-window_lG,totalrefT-window_lT
        output.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n"%(chrom,window_start,pos,cmut,ctotal,cmutCpG,ctotalCpG,cA,cC,cG,cT,0 if ctotal == 0 else cmut/float(ctotal),0 if ctotalCpG == 0 else cmutCpG/float(ctotalCpG),0 if ctotal == 0 else cA/float(ctotal),0 if ctotal == 0 else cC/float(ctotal),0 if ctotal == 0 else cG/float(ctotal),0 if ctotal == 0 else cT/float(ctotal)))

      output.write("\n##STATS\n")
      output.write("#A\tC\tG\tT\tCpGs\n")
      output.write("\t".join(map(str,[totalrefA,totalrefC,totalrefG,totalrefT,totalCpG]))+"\n")
      output.write("#y\tN\tAC\tAG\tAT\tCA\tCG\tCT\tGA\tGC\tGT\tTA\tTC\tTG\n")
      output.write("\t".join(map(str,[mut,total,ACn,AGn,ATn,CAn,CGn,CTn,GAn,GCn,GTn,TAn,TCn,TGn]))+"\n")
      output.write("#yCpG\tNCpG\tCA\tCG\tCT\tGA\tGC\tGT\n")
      output.write("\t".join(map(str,[mutCpG,totalCpG,CA,CG,CT,GA,GC,GT]))+"\n")
      output.write("##INSERTIONS\n")
      output.write("#len\tcount\n")
      for key,value in insertionsizes.iteritems():
        output.write("%d\t%d\n"%(key,value))
      output.write("##DELETIONS\n")
      output.write("#len\tcount\n")
      for key,value in deletionsizes.iteritems():
        output.write("%d\t%d\n"%(key,value))
except:
  exc_type, exc_value, exc_traceback = sys.exc_info()
  sys.stderr.write("%s\n"%str(exc_value))
  traceback.print_tb(exc_traceback)
  sys.stderr.write('Script terminated early. Printing current values.\n')

  if total > 0:
    if (total-window_ltotal)/100000. > 0.1:
      cmut = mut-window_lmut
      ctotal = total-window_ltotal
      cmutCpG = mutCpG-window_lmutCpG
      ctotalCpG = totalCpG-window_ltotalCpG
      cA,cC,cG,cT = totalrefA-window_lA,totalrefC-window_lC,totalrefG-window_lG,totalrefT-window_lT
      output.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n"%(chrom,window_start,pos,cmut,ctotal,cmutCpG,ctotalCpG,cA,cC,cG,cT,cmut/float(ctotal),cmutCpG/float(ctotalCpG),cA/float(ctotal),cC/float(ctotal),cG/float(ctotal),cT/float(ctotal)))

    output.write("\n##STATS\n")
    output.write("#A\tC\tG\tT\tCpGs\n")
    output.write("\t".join(map(str,[totalrefA,totalrefC,totalrefG,totalrefT,totalCpG]))+"\n")
    output.write("#y\tN\tAC\tAG\tAT\tCA\tCG\tCT\tGA\tGC\tGT\tTA\tTC\tTG\n")
    output.write("\t".join(map(str,[mut,total,ACn,AGn,ATn,CAn,CGn,CTn,GAn,GCn,GTn,TAn,TCn,TGn]))+"\n")
    output.write("#yCpG\tNCpG\tCA\tCG\tCT\tGA\tGC\tGT\n")
    output.write("\t".join(map(str,[mutCpG,totalCpG,CA,CG,CT,GA,GC,GT]))+"\n")
    output.write("##INSERTIONS\n")
    output.write("#len\tcount\n")
    for key,value in insertionsizes.iteritems():
      output.write("%d\t%d\n"%(key,value))
    output.write("##DELETIONS\n")
    output.write("#len\tcount\n")
    for key,value in deletionsizes.iteritems():
      output.write("%d\t%d\n"%(key,value))