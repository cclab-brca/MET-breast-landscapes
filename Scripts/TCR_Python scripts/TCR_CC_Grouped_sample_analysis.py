## Script for TCR repertoire paper:
## The integrated genomic and immune landscapes of lethal metastatic breast cancer
## Written by Dr R. Bashford-Rogers, Department of Medicine, University of Cambridge, UK

## To run
# python TCR_CC_Grouped_sample_analysis.py <SAMPLE_GROUPING_FILE> <SAMPLE_GROUP_INDEX> 
# where: 
#   The sequence file is in /DIRECTORY_FOR_DATA/NETWORKS/Fully_reduced_<SAMPLE_ID>.fasta
#   The IMGT annotation files are in  /DIRECTORY_FOR_DATA/ANNOTATIONS/IMGT_ISOTYPER/IMGT_<SAMPLE_ID>_3_Nt-sequences.txt
#   The cluster identities file (output from code from Bashford-Rogers et al, 2013) is in /DIRECTORY_FOR_DATA/NETWORKS/Cluster_identities_<SAMPLE_ID>.txt
#   And an output directory created /DIRECTORY_FOR_DATA/CDR3_TEMPORAL/
#   Example <SAMPLE_GROUPING_FILE> given in: Samples_grouping_example_file.txt

#!/usr/bin/python
import math
import sys
from collections import defaultdict
import commands
from operator import itemgetter
from itertools import izip
from operator import add
import numpy as np
import random
from numpy import median, mean
import networkx as nx
from itertools import izip


def fasta_iterator(fh):
  while True:
    line = fh.readline()
    if line.startswith('>'): break	
  while True:
    header = line[1:-1].rstrip()
    sequence = fh.readline().rstrip()
    while True:
      line = fh.readline()
      if not line: break
      if line.startswith('>'): break
      sequence += line.rstrip()
    yield(header, sequence)
    if not line: return

def Get_sequences(file):
  fh = open (file, "r")
  seqs = {}
  for header,sequence in fasta_iterator(fh):
    header = header.replace(":","").split()[0].replace("-","")
    header = header.split("#")[0]
    seqs[header]=sequence
  fh.close()
  return(seqs)

class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

def CDR3_defined_sequences_IMGT(grouped, pat, dir, combined_sequences_file,header_file,CDR3_info_file,merge_cluster_file):
  for f in [CDR3_info_file,combined_sequences_file,merge_cluster_file]:
    fh=open(f,"w")
    fh.close()
  seqs = Tree()
  ids,ind,chains_use = [],0,{}
  fh=open(CDR3_info_file,"w")
  fh.close() 
  out, ind = '#sample\tid\tv\tj\tCDR3 region\tV:J:CDR3\n',0
  cluster_codes,codes_clusters=Tree(),Tree()
  all_clusters,total_seq = Tree(),0
  all_seqs={} 
  for info in grouped[pat]:
    sample = info[1]
    print sample
    seq_file = info[4]+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+sample+".fasta"
    annot_file = info[4]+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+sample+"_3_Nt-sequences.txt"
    cluster_file = info[4]+"ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"+sample+".txt"
    fh=open(seq_file,"r")
    for header,sequence in fasta_iterator(fh):
      all_seqs[header.split("__")[0]] = [header+"|"+sample, sequence]
    fh.close()
    fh=open(annot_file,"r")
    out, ind = '#sample\tid\tV\tJ\tCDR3\tcode\n',0
    codes, inv_codes={},Tree()
    for l in fh:
      if(l[0]!="#"):
        l=l.strip().split("\t")
        if(l[2].count("unknown")==0):
          if(len(l)>=16):
            cdr3,v,j = l[15],l[3].split()[1],l[4].split()[1]
            v,j = v.split("*")[0], j.split("*")[0]
            if(len(cdr3)<=7):cdr3 = "-"
            out=out+sample+"\t"+l[1]+"\t"+v+"\t"+j+"\t"+cdr3+"\t"+v+":"+j+":"+cdr3+"\n"
            if(cdr3 != "-"):
              codes[l[1].split("__")[0]] = v+":"+j+":"+cdr3
              inv_codes[v+":"+j+":"+cdr3][l[1].split("__")[0]].value = 1
            ind = ind+1
            if(ind>500):
              Write_output(out, CDR3_info_file)
              out, ind = '',0
    Write_output(out, CDR3_info_file)
    out, ind = '',0
    fh.close()
    fh=open(cluster_file,"r")
    for l in fh:
      if(l[0]!="#"):
        l=l.strip().split()
        clust, id = l[1]+":"+sample,l[2].split("__")[0]
        all_clusters[clust][id].value = 1
        total_seq = total_seq+1
        if(id in codes):
          cluster_codes[clust][codes[id]].value = 1
          codes_clusters[codes[id]][clust].value = 1
    fh.close()
  cluster_sizes = []
  for c in all_clusters:
    cluster_sizes.append([c,len(all_clusters[c])])
  cluster_sizes.sort(key=lambda x: x[1], reverse = True)
  clusters_found,codes_found = {},{}
  clusters_merged = Tree()
  per_cluster_codes = {}
  for c in cluster_sizes:
    c = c[0]
    if(c not in clusters_found):
      cluster_merge_id = "M:"+c
      new_codes = []
      clusters_found[c] = cluster_merge_id
      for code in cluster_codes[c]:
        per_cluster_codes[code] = cluster_merge_id
        if(len(codes_clusters[code])>1):
          for c1 in codes_clusters[code]:
            if(c1 != c):
              if(c1 not in clusters_found):
                clusters_found[c1]=cluster_merge_id
                for code1 in cluster_codes[c1]:
                  if(code1 not in per_cluster_codes):
                    if(code1 not in new_codes):
                      new_codes.append(code1)
                      per_cluster_codes[code1] = cluster_merge_id
      while(len(new_codes)>0):
        additional_clusters,additional_codes = [],[]
        for c1 in new_codes:
          for clust1 in codes_clusters[c1]:
            if(clust1 not in clusters_found):
              clusters_found[clust1]=cluster_merge_id
              additional_clusters.append(clust1)
              for codes1 in cluster_codes[clust1]:
                if(codes1 not in per_cluster_codes):
                  per_cluster_codes[codes1] = cluster_merge_id
                  additional_codes.append(codes1)
        new_codes = additional_codes
  print "CODES:",len(codes_clusters),len(per_cluster_codes)
  print "CLUSTERS:",len(clusters_found), len(cluster_codes)
  del new_codes, per_cluster_codes,cluster_codes,codes_clusters
  out,ind="#merge cluster\tid\toriginal cluster\n",0
  ids_found = 0
  clustered_sequence = Tree()
  for c in all_clusters:
    if(c in clusters_found):
      for id in all_clusters[c]:
        out=out+clusters_found[c]+"\t"+id+"\t"+c+"\n"
        clustered_sequence[clusters_found[c]][id].value = 1
        ids_found = ids_found+1
        ind = ind+1
        if(ind>500):
          Write_output(out,merge_cluster_file)
          out, ind = '',0
    else:
      for id in all_clusters[c]:
        clustered_sequence["U:"+c][id].value = 1
        out=out+"U:"+c+"\t"+id+"\t"+c+"\n"
        ind = ind+1
        if(ind>500):
          Write_output(out,merge_cluster_file)
          out, ind = '',0
  Write_output(out,merge_cluster_file)
  print "ORIG seqs:",total_seq, "FOUND IDS",ids_found
  out,ind = '',0
  for c in clustered_sequence:
    for id in clustered_sequence[c]:
      out=out+">"+c+"|"+all_seqs[id][0]+"\n"+all_seqs[id][1]+"\n"
      ind = ind+1
      if(ind>400):
        Write_output(out,combined_sequences_file)
        out,ind = '',0
  Write_output(out,combined_sequences_file)
  return()


def Get_overlapping_istoype_cluster_summaries (grouped,merge_cluster_file,dir,overlapping_subsampled,isotype_overlapping_subsampled):
  subsample_names, subsample_ids = {},{}
  totals,freqs,cluster_inv = {},{},{}
  freq_inv,classifications = {},{}
  threshold = 10000000 ##################### remove
  for id in grouped[pat]:
    subsample_cluster_file = dir+"ISOTYPER/Subsample_"+id[1]+".txt"
    fh = open(subsample_cluster_file, "r")
    subsample_names[id[1]] =[]
    totals[id[1]] =0
    ind = 0
    l1 = ''
    for l in fh:
      if(l[0]!="#"):
        l=l.strip().split()
        if(l[0] not in subsample_names[id[1]]):
          ind = ind+1
          subsample_names[id[1]] = subsample_names[id[1]]+[l[0]]
          totals[l[0]] =0
        if(totals[l[0]]<threshold):
          l1 = l[2]
          id1,subsample = l[2].split("__")[0], l[0]
          f = map(int, l[2].split("__")[1].split("|")[0].split("_"))
          if(id1 in subsample_ids):subsample_ids[id1] = subsample_ids[id1] + [l[0]]
          else:subsample_ids[id1] = [l[0]]
          totals[l[0]] =totals[l[0]] + sum(f)
          freqs[id1] = f
          cluster_inv[subsample+":"+id1] = l[1]
          n = subsample+":"+id1
          freq_inv[n] = f
    fh.close()
    if(len(l1)>0):
      classes = []
      for c in l1.split("|")[1].split("_"):
        classes.append(c)
      for subsample in subsample_names[id[1]]:
        classifications[subsample] = classes
  fh=open(merge_cluster_file, "r")
  clusters = Tree()
  clusters_subsample = Tree()
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[1] in subsample_ids):
        clusters[l[0]][l[2].split(":")[1]][l[1]].value = 1
        for subsample in subsample_ids[l[1]]:
          clusters_subsample[subsample][l[0]].value = 1
  fh.close()
  overlapping1, overlapping2, overlapping_clusters = {},{},{}
  overlap_cluster_isotypes = {}
  fh=open(isotype_overlapping_subsampled,"w")
  fh.close()
  chains = ["IGHA1/2","IGHD","IGHE","IGHG1/2","IGHG3","IGHG4","IGHM"]
  out = "#sample1\tsample2\trepID\tCluster_ID\treads1\t"+"\t".join(chains)+"\treads2\t"+"\t".join(chains)+"\n"
  ind = 0
  for c in clusters:
    if(len(clusters[c])>=1):
      overlap_ids = []
      overlap_frequencies = {}
      sample_array = []
      for id in grouped[pat]:
        if(id[1] in clusters[c]):
          sample_array.append(id[1])
          array_subsamples = {}
          total_sub_cluster = []
          for id1 in clusters[c][id[1]]:
            for subsample in subsample_ids[id1]:
              total_sub_cluster = freq_inv[subsample+":"+id1] 
              fx = [0]*len(chains)
              for i in range(0,len(classifications[subsample])):
                if(classifications[subsample][i] in chains):
                  index = chains.index(classifications[subsample][i])
                  fx[index] = fx[index]+total_sub_cluster[i]
              overlap_cluster_isotypes[subsample] = fx
              if(subsample in array_subsamples):array_subsamples[subsample] = array_subsamples[subsample]+sum(freqs[id1])
              else:array_subsamples[subsample]=sum(freqs[id1])
          subs = []
          for subsample in array_subsamples:
            subs.append(subsample)
            overlap_frequencies[subsample] = array_subsamples[subsample]
          overlap_ids.append(subs)
      array_subsamples= overlap_frequencies
      for subsample1 in range(0,len(sample_array)):
        for subsample2 in range(subsample1,len(sample_array)):
          if(subsample1<subsample2):
            for i1 in range(0,len(overlap_ids[subsample1])):
              for i2 in range(0,len(overlap_ids[subsample2])):
                name = sample_array[subsample1]+"\t"+sample_array[subsample2]+"\t"+overlap_ids[subsample1][i1]+":"+overlap_ids[subsample2][i2]
                if(name in overlapping1):overlapping1[name] = overlapping1[name]+array_subsamples[overlap_ids[subsample1][i1]]
                else:overlapping1[name]=array_subsamples[overlap_ids[subsample1][i1]]
                if(name in overlapping2):overlapping2[name] = overlapping2[name]+array_subsamples[overlap_ids[subsample2][i2]]
                else:overlapping2[name]=array_subsamples[overlap_ids[subsample2][i2]]
                if(name in overlapping_clusters):overlapping_clusters[name] = overlapping_clusters[name] +1
                else:overlapping_clusters[name] =1
                out=out+name+"\t"+"\t".join(map(str, [c, sum(overlap_cluster_isotypes[overlap_ids[subsample1][i1]])]+overlap_cluster_isotypes[overlap_ids[subsample1][i1]]+[sum(overlap_cluster_isotypes[overlap_ids[subsample2][i2]])]+ overlap_cluster_isotypes[overlap_ids[subsample2][i2]]  ))+"\n"
                ind = ind+1
                if(ind>300):
                  Write_output(out, isotype_overlapping_subsampled)
                  out, ind = '',0
  Write_output(out, isotype_overlapping_subsampled)
  ind =0
  out = "#sample1\tsample2\trepID\treads1\treads2\tclusters1\tclusters2\toverlap_reads1\toverlap_reads2\toverlap_clusters\n"
  for ind1 in range(0,len(grouped[pat])):
    for ind2 in range(ind1,len(grouped[pat])):
      if(ind1<ind2):
        for subsample1 in subsample_names[grouped[pat][ind1][1]]:
          for subsample2 in subsample_names[grouped[pat][ind2][1]]:
            name =grouped[pat][ind1][1]+"\t"+grouped[pat][ind2][1]+"\t"+subsample1+":"+subsample2
            s1, s2 =subsample1,subsample2
            if(name in overlapping_clusters):
              out=out+name+"\t"+"\t".join(map(str, [totals[s1], totals[s2], len(clusters_subsample[s1]),len(clusters_subsample[s2]),overlapping1[name], overlapping2[name], overlapping_clusters[name] ]))+"\n"
            else:
              out=out+name+"\t"+"\t".join(map(str,[totals[s1], totals[s2],len(clusters_subsample[s1]),len(clusters_subsample[s2]),0,0,0]))+"\n"
  fh=open(overlapping_subsampled, "w")
  fh.write(out)
  fh.close()
  return()

def Get_overlap_between_samples(cluster_file, header_file, overlap_file,pat):
  cols = []
  fh = open(header_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cols = cols+[l[1]]
  fh.close()
  fh=open(cluster_file, "r")
  clusters,per_seq_overlap, per_cluster_overlap,col_length = {},{},{},len(cols)
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      clust, freqs = int(l[1]), map(int, l[3:3+col_length])
      if(clust in clusters):clusters[clust] = map(add, clusters[clust], freqs)
      else:clusters[clust] =freqs
      nz = [i for i in range(len(freqs)) if freqs[i]!=0]
      if(len(nz)>1):
        for i in range(0,len(nz)):
          for j in range(i,len(nz)):
            if(i<j):
              name = cols[nz[i]]+"\t"+cols[nz[j]]
              if(name in per_seq_overlap):per_seq_overlap[name] = map(add, per_seq_overlap[name], [freqs[nz[i]],freqs[nz[j]]])
              else:per_seq_overlap[name] = [freqs[nz[i]],freqs[nz[j]]]
  fh.close()
  total = [0]*col_length
  for c in clusters:
    freqs = clusters[c]
    nz = [i for i in range(len(freqs)) if freqs[i]!=0]
    total = map(add, total, freqs)
    if(len(nz)>1):
      for i in range(0,len(nz)):
        for j in range(i,len(nz)):
          if(i<j):
            name = cols[nz[i]]+"\t"+cols[nz[j]]
            if(name in per_cluster_overlap):per_cluster_overlap[name] = map(add, per_cluster_overlap[name],[freqs[nz[i]],freqs[nz[j]]])
            else:per_cluster_overlap[name]=[freqs[nz[i]],freqs[nz[j]]]
  out="#ID\tsample1\tsample2\tper_cluster_overlap1\tper_cluster_overlap2\tper_sequence_overlap1\tper_sequence_overlap2\ttotal_sample1\ttotal_sample2\n"
  print total, cols
  for name in per_cluster_overlap:
    if(name in per_seq_overlap):overlap_seq = "\t".join(map(str,per_seq_overlap[name]))
    else:overlap_seq = "0\t0"
    totals = [total[cols.index(name.split("\t")[0])], total[cols.index(name.split("\t")[1])]]
    out=out+pat+"\t"+name+"\t"+"\t".join(map(str,per_cluster_overlap[name]))+"\t"+overlap_seq+"\t"+"\t".join(map(str, totals))+"\n"
  fh=open(overlap_file,"w")
  fh.write(out)
  fh.close()
  return()

def Write_output(out, file):
  fh=open(file, "a")
  fh.write(out)
  fh.close()
  return ()

def Get_info(file):
  fh=open(file,"r")
  ids, grouped,previous = {},{},''
  sample_group_RNA=[]
  for l in fh:
     if(l[0]!="#"):
       l=l.replace("/lustre/scratch108/viruses/rbr1/","/lustre/scratch118/infgen/team146/rbr1/")
       l1 = l.strip()
       l=l.split()
       #if(len(l)<8 or l[7]!="No"):
       if(len(l)>=7):
         pat, sample_id,index,dir,outdir = l[0],l[1],l[2],l[5], l[6]
         if(pat!=previous):
           if(len(sample_group_RNA)>=1):grouped[previous] = sample_group_RNA
           else:
             sample_group_RNA.append([pat, sample_id,index,dir,outdir])
             grouped[pat] = sample_group_RNA
           sample_group_RNA=[]
           sample_group_RNA.append([pat, sample_id,index,dir,outdir])
         else:sample_group_RNA.append([pat, sample_id,index,dir,outdir])
         previous = pat
  if(len(sample_group_RNA)>=1):grouped[pat] = sample_group_RNA
  fh.close()
  return(grouped)


###########################
file = sys.argv[1]
index_required = int(sys.argv[2])
grouped=Get_info(file)
batch = file.split("Samples_")[1].split(".")[0]

########################### Files for QC and filtering
index = 0
for pat in grouped:
  pat_sub=pat
  index = index+1
  if(index == index_required):
    for info in grouped[pat]:
      print pat, index
      id,dir,outdir = info[1],info[4]+"ORIENTATED_SEQUENCES/", info[4]+"ORIENTATED_SEQUENCES/TEMPORAL/"
      outdir = info[4]+"ORIENTATED_SEQUENCES/CDR3_TEMPORAL/"
      combined_sequences_file =outdir+"Combined_fully_reduced_sequences_"+pat+".fasta"
      combined_sequences_file_for_total_analysis = outdir+"Fully_reduced_"+pat+".fasta"
      CDR3_info_file = outdir+"CDR3_information_"+pat+".txt"
      merge_cluster_file = outdir+"Merge_clustering_"+pat+".txt"
      header_file = outdir+"Headers_"+pat+".txt"
      overlapping_subsampled = outdir+"Subsampled_overlapping_values_"+pat+".txt"
      isotype_overlapping_subsampled = outdir+"Subsampled_overlapping_isotypevalues_"+pat+".txt"
      read_number_division = "__"
######################### Commands
      CDR3_defined_sequences_IMGT(grouped, pat, dir, combined_sequences_file,header_file,CDR3_info_file,merge_cluster_file)
      Get_overlapping_istoype_cluster_summaries (grouped,merge_cluster_file,dir,overlapping_subsampled,isotype_overlapping_subsampled)
      break

      
     
