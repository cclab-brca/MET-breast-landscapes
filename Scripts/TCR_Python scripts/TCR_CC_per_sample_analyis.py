## Script for TCR repertoire paper:
## The integrated genomic and immune landscapes of lethal metastatic breast cancer
## Written by Dr R. Bashford-Rogers, Department of Medicine, University of Cambridge, UK

## To run
# python TCR_CC_per_sample_analyis.py <SAMPLE_ID> <SAMPLE_ID> <DIRECTORY_FOR_DATA> <SPECIES> 
# where: <SPECIES>  = HOMO_SAPIENS
#   The sequence file is in /DIRECTORY_FOR_DATA/NETWORKS/Fully_reduced_<SAMPLE_ID>.fasta
#   The IMGT annotation files are in  /DIRECTORY_FOR_DATA/ANNOTATIONS/IMGT_ISOTYPER/IMGT_<SAMPLE_ID>_1_Summary.txt
#                                     /DIRECTORY_FOR_DATA/ANNOTATIONS/IMGT_ISOTYPER/IMGT_<SAMPLE_ID>_3_Nt-sequences.txt
#                                     /DIRECTORY_FOR_DATA/ANNOTATIONS/IMGT_ISOTYPER/IMGT_<SAMPLE_ID>_5_AA-sequences.txt
#                                     /DIRECTORY_FOR_DATA/ANNOTATIONS/IMGT_ISOTYPER/IMGT_<SAMPLE_ID>_7_V-REGION-mutation-and-AA-change-table.txt
#                                     /DIRECTORY_FOR_DATA/ANNOTATIONS/IMGT_ISOTYPER/IMGT_<SAMPLE_ID>_8_V-REGION-nt-mutation-statistics.txt
# The cluster identities file (output from code from Bashford-Rogers et al, 2013) is in /DIRECTORY_FOR_DATA/NETWORKS/Cluster_identities_<SAMPLE_ID>.txt
# And an output directory created /DIRECTORY_FOR_DATA/ISOTYPER/

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

class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

def fasta_iterator(fh):
  while True:
    line = fh.readline()
    if line.startswith('>'): break	
    elif(len(line)==0):break
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

def Get_annotations(annot_file):
  CDR3s,mutations = Tree(),{}
  fh=open(annot_file,"r")
  for l in fh:
    l=l.strip().split("\t")
    if(l[0]!='Sequence number'):
      if(len(l)>=21 and l[2].count('productive')!=0):
        id,v_muts, j_muts,cdr3 = l[1].split("__")[0],l[6], l[12], l[20].replace("#","")
        v,j = l[3], l[9]
        if(v_muts.count('nt')!=0 and j_muts.count('nt')!=0 and len(cdr3)>1):
          v = v.split()[1].split("*")[0]
          j = j.split()[1].split("*")[0]
          cdr3 = cdr3.split()[0]+":"+v+":"+j
          v_muts, j_muts = map(int,v_muts.split()[0].split("/")), map(int,j_muts.split()[0].split("/"))
          v_muts, j_muts = v_muts[1]-v_muts[0], j_muts[1]-j_muts[0]
          mutations[id] = [v_muts, j_muts,cdr3, v+"\t"+j]
          CDR3s[cdr3][id].value = 1
  fh.close()
  return(CDR3s,mutations)

def Get_sequences(seq_file):
  fh=open(seq_file,"r")
  seqs,freqs,alias = {},{},{}
  for header,sequence in fasta_iterator(fh):
    seqs[header.split("__")[0]]=sequence
    f = map(int, header.split("__")[1].split("|")[0].split("_"))
    freqs[header.split("__")[0]] = f
    alias[header.split("__")[0]] = header
  fh.close()
  return(seqs, freqs, alias)

def Get_annotation(annot_file):
  fh=open(annot_file,"r")
  annots = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(len(l)>=19):
        v,j,v_mm,j_mm = l[1],l[13],int(l[len(l)-4]), int(l[len(l)-3])
        if(v_mm+j_mm>60):print l
        annots[l[0].split("__")[0]] = [v,j,v_mm,j_mm]
  fh.close()
  return(annots)

def Get_cluster_occupancy(merged_cluster_file):
  clusters, inv_clusters= Tree(),{}
  fh=open(merged_cluster_file, "r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split() 
      cluster, id = l[1],l[2].split("__")[0]
      clusters[cluster][id].value = 1
      inv_clusters[id]=cluster
  fh.close()
  return(clusters, inv_clusters)

def Classify_sequences_into_developmental_stages_per_cluster(sample, annot_file_IMGT, cluster_file,per_cluster_developmental_classification_file,reverse_primer_group):
  clusters, inv_clusters = Get_cluster_occupancy(cluster_file)
  CDR3s,mutations = Get_annotations(annot_file_IMGT)
  seqs,freqs,alias = Get_sequences(seq_file)
  total = 0 
  for id in freqs:
    if(id in inv_clusters):
      l=len(freqs[id])
      chains= alias[id].split("|")[1].split("_")
      total = total+sum(freqs[id])
  print total
  for i in range(0,len(chains)):
    chains[i] = chains[i].split("*")[0]
    if(chains[i].count("TRBC")==1):chains[i] = "TRBC"
  cluster_freqs,chains = {}, np.array(chains)
  fh=open(per_cluster_developmental_classification_file,"w")
  fh.close()
  out,ind = '#ID\tsequence\tclassifiation\tall_classes\tV\tJ\tmutation\tCDR3\tnode size\tnode %\tcluster ID\tcluster size\n',0
  print len(alias)
  for c in clusters: 
    if(len(clusters[c])>0):
      freq,mm,vj = [0]*l,[],{}
      for id in clusters[c]:
        freq = map(add, freq, freqs[id])
        if(id in mutations):
          mm = mm+[mutations[id][0]+mutations[id][1]]
          vj_class = mutations[id][3]
          if(vj_class in vj):vj[vj_class] = vj[vj_class]+sum(freqs[id])
          else:vj[vj_class]=sum(freqs[id])
      max_vj_f,max_vj = 0,'-\t-'
      for i in vj:
        if(max_vj_f<vj[i]):max_vj_f,max_vj = vj[i], i
      if(len(mm)==0):mm=-1
      else:mm=mean(mm)
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      chain = ",".join(np.sort(np.unique(chains[nz])))
      chains_used = np.sort(np.unique(chains[nz]))
      classifiations=[chain]
      MD = 0
      if(len(chains_used)>1):
        for i in range(0,len(chains_used)):
          classifiations.append(chains_used[i])
      for id in clusters[c]:
        out=out+sample+"\t"+id+"\t"+"|".join(classifiations)+"\t"+chain+"\t"+max_vj+"\t"+str(mm)+"\t"+"CDR3"+"\t"+str(sum(freq))+"\t"+str(sum(freqs[id])*100/total)+"\t"+c+"\t"+str(sum(freq)*100.0/total)+"\n"
        ind = ind+1
        if(ind>500):
          Write_out(out,per_cluster_developmental_classification_file)
          out, ind = '',0
  Write_out(out,per_cluster_developmental_classification_file)
  return()

def Classify_sequences_into_developmental_stages(sample, annot_file, cluster_file, developmental_classification_file,reverse_primer_group):
  CDR3s,annots = Get_annotations(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  total = 0
  for id in freqs:
    chains= alias[id].split("|")[1].split("_")
    total = total+sum(freqs[id])
  for i in range(0,len(chains)):
    chains[i] = chains[i].split("*")[0]
  chains = np.array(chains)
  out,ind = '#ID\tsequence\tclassifiation\tall_classes\tV\tJ\tmutation\tCDR3\tnode size\tnode %\n',0
  fh=open(developmental_classification_file, "w")
  fh.close()
  for id in freqs:
    mm,info = -1,"-\t-\t-1"
    if(id in annots):
      #info = "\t".join(map(str,annots[id]))
      info = annots[id][3]+"\t"+str(annots[id][0])+"\t"+annots[id][2].split(":")[0]
      mm = annots[id][0]
    freq = freqs[id]
    nz = [i for i in range(len(freq)) if freq[i]!=0]
    chain = ",".join(np.sort(np.unique(chains[nz])))
    chains_used = np.sort(np.unique(chains[nz]))
    classifiations=[chain]
    if(len(chains_used)>1):
      for i in range(0,len(chains_used)):
        classifiations.append(chains_used[i])
    out=out+sample+"\t"+id+"\t"+"|".join(classifiations)+"\t"+chain+"\t"+info+"\t"+str(sum(freq))+"\t"+str(sum(freq)*100.0/total)+"\n"
    ind = ind+1
    if(ind>500):
      Write_out(out, developmental_classification_file)
      out, ind = '',0
  Write_out(out, developmental_classification_file)
  return()

def Gini_index(cpoints,cvdf, vpoints,vvdf): 
  (vgini)=Get_Gini(vpoints,vvdf)
  (cgini)=Get_Gini(cpoints,cvdf)
  return(vgini, cgini)

def Get_Gini(n,v):
  values=[]
  for i in range(0,len(n)):
    for j in range(0,v[i]):
      values.append(n[i])
  n = len(values)
  assert(n > 0), 'Empty list of values'
  sortedValues = sorted(values) #Sort smallest to largest
  cumm = [0]
  for i in range(n):
    cumm.append(sum(sortedValues[0:(i + 1)]))
  LorenzPoints = [[], []]
  sumYs = 0           #Some of all y values
  robinHoodIdx = -1   #Robin Hood index max(x_i, y_i)
  for i in range(1, n + 2):
    x = 100.0 * (i - 1)/n
    y = 100.0 * (cumm[i - 1]/float(cumm[n]))
    LorenzPoints[0].append(x)
    LorenzPoints[1].append(y)
    sumYs += y
    maxX_Y = x - y
    if maxX_Y > robinHoodIdx: robinHoodIdx = maxX_Y   
  giniIdx = 100 + (100 - 2 * sumYs)/n #Gini index 
  return(giniIdx/100)

def Get_network_parameters_per_classificiation_subsample(samplex,per_cluster_developmental_classification_file,per_cluster_developmental_classification_network_parameters,reverse_primer_group,cluster_file,developmental_classification_file):
  fh = open(developmental_classification_file,"r")
  muts1,muts = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id, mut = l[1], int(l[6])
      muts1[id] = mut
  fh.close()
  fh = open(per_cluster_developmental_classification_file,"r")
  classification_ids,classification_clusters = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id, cluster, mut, types = l[1],l[10], float(l[6]),l[2].replace("IGHA","IGHA1/2").replace("IGHA1/21/2","IGHA1/2").split("|")
      if(id in muts1):
        muts[id] = [muts1[id], cluster, types]
      else:
        print id,"FAIL"
  fh.close()
  del muts1
  print "READ DEVELOP FILE"
  fh=open(cluster_file,"r")
  v_sizes,c_sizes = {},Tree()
  index = 0
  sample_tree,reads_per_sample = Tree(),{}
  cluster_size_cumul = {}
  inv_cluster ={}
  cluster_per_chain ={}
  for l in fh:
    index=index+1
    if (index>1):
      l=l.strip().split()
      id,freq, id_short,clust = l[2],map(int, l[2].split("__")[1].split("|")[0].split("_")),l[2].split("__")[0], l[1]
      clust = l[1]
      classification = l[2].split("|")[1].split("_")
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      classes = {}
      inv_cluster[id_short]=clust
      for i in nz:
        c = classification[i].split("*")[0]
        if(c.count("TRBC")==1):c = "TRBC"
        if(c in classes):classes[c] = classes[c]+freq[i]
        else:classes[c]=freq[i]
      for c1 in classes:
        c = c1
        if(c in v_sizes):
          v_sizes[c] = v_sizes[c]+[classes[c1]]
          reads_per_sample[c] = reads_per_sample[c] +classes[c1]
        else:
          reads_per_sample[c]=classes[c1]
          v_sizes[c]=[classes[c1]]
        if(c in cluster_per_chain):cluster_per_chain[c]=cluster_per_chain[c]+[id_short]*classes[c1]
        else:cluster_per_chain[c]=[id_short]*classes[c1]
      classes1 = muts[id_short][2]
      for c1 in classes1:
        c = c1
  fh.close()
  subsample_depth = {}
  for c in classification:
    c = classification[i].split("*")[0]
    if(c.count("TRBC")==1):c = "TRBC"
    subsample_depth[c] = 1000
  out="#Id\tIsotype\tN reads\tN vertices\tVertex Gini Index\tCluster Unique Gini Index\tCluster Total Gini Index\tLargest Cluster (total %)\t2nd Largest Cluster (total %)\tLargest Cluster (unique %)\t2nd Largest Cluster (unique %)\tnumber of unique clusters\ttotal reads per sample\tsample\tmean_unique_cluster_size\tmean_total_cluster_size\n"
  repeats =20
  for c in v_sizes:
    if(c in subsample_depth):
      if(subsample_depth[c]>0):
        depth = subsample_depth[c]
        v = v_sizes[c]
        print sum(v), depth
        if(sum(v)>depth):
          prob = np.array(v)/(sum(v)*1.0)
          xlen = len(v)
          for r in range(0,repeats):
            sample = np.random.choice(xlen, depth,p=prob)
            unique, counts = np.unique(sample, return_counts=True)
            vpoints,vvdf = np.unique(counts, return_counts=True)
            vgini=Get_Gini(vpoints,vvdf)

            if(c in cluster_per_chain):
              cluster_sizes_sub = cluster_per_chain[c]
              sample_cs = np.random.choice(cluster_sizes_sub, depth)
              unique, counts = np.unique(sample_cs, return_counts=True)
              unique_BCRs_per_cluster, total_BCRs_per_cluster = {},{}
              unique_cluster_sizes, total_cluster_sizes = [],[]
              for i in range(0,len(unique)):
                id = unique[i]
                c1 = inv_cluster[id]
                if(c1 in unique_BCRs_per_cluster):
                  unique_BCRs_per_cluster[c1] = unique_BCRs_per_cluster[c1]+1
                  total_BCRs_per_cluster[c1] = total_BCRs_per_cluster[c1]+counts[i]
                else:
                  unique_BCRs_per_cluster[c1]=1
                  total_BCRs_per_cluster[c1]=counts[i]
              for c1 in unique_BCRs_per_cluster:
                unique_cluster_sizes = unique_cluster_sizes+[unique_BCRs_per_cluster[c1]]
                total_cluster_sizes=total_cluster_sizes+[total_BCRs_per_cluster[c1]]
              mean_unique_cluster_size,mean_total_cluster_size = mean(unique_cluster_sizes), mean(total_cluster_sizes)
              unique_uniq_cs, counts_uniq_cs = np.unique(unique_cluster_sizes, return_counts=True)
              unique_total_cs, counts_total_cs = np.unique(total_cluster_sizes, return_counts=True)
              uniq_clust_gini=Get_Gini(unique_uniq_cs,counts_uniq_cs)
              total_clust_gini=Get_Gini(unique_total_cs, counts_total_cs)
              max_pop, max_1_pop = unique_total_cs[len(unique_total_cs)-1]*100.0/sum(total_cluster_sizes), unique_total_cs[len(unique_total_cs)-2]*100.0/sum(total_cluster_sizes)
              max_pop_uniq, max_1_pop_uniq = unique_uniq_cs[len(unique_uniq_cs)-1]*100.0/sum(unique_cluster_sizes), unique_uniq_cs[len(unique_uniq_cs)-2]*100.0/sum(unique_cluster_sizes)
            else:
              max_pop, max_1_pop,max_pop_uniq, max_1_pop_uniq  = -1,-1,-1,-1
              uniq_clust_gini,total_clust_gini = -1,-1
            out = out+"REP"+str(r)+"\t"+c+"\t"+str(len(sample))+"\t"+str(len(unique))+"\t"+str(vgini)+"\t"+str(uniq_clust_gini)+"\t"+str(total_clust_gini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\t"+str(max_pop_uniq)+"\t"+str(max_1_pop_uniq)+"\t"+str(len(unique))+"\t"+str(sum(v))+"\t"+samplex+"\t"+str(mean_unique_cluster_size)+"\t"+str(mean_total_cluster_size)+"\n"
  fh=open(per_cluster_developmental_classification_network_parameters, "w")
  fh.write(out)
  fh.close()
  return()

def Get_expanded_cluster_summary(cluster_SUBSAMPLE_information_file,cluster_isotype_expansion_SUBSAMPLE_file,sample,developmental_classification_file,per_cluster_developmental_classification_file):
  fh = open(developmental_classification_file,"r")
  muts,muts_ids = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,mut = l[1],float(l[6])
      muts_ids[id] = mut
  fh.close()
  fh = open(per_cluster_developmental_classification_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,clust = l[1],l[10]
      if(clust in muts):muts[clust] =muts[clust] +[muts_ids[id]]
      else:muts[clust] = [muts_ids[id]]
  fh.close()
  del muts_ids
  fh=open(cluster_SUBSAMPLE_information_file,"r")
  total_BCRs_per_isotype,unique_BCRs_per_isotype, singletons_per_isotype = {},{},{}
  expanded_clonal_range, expanded_clonal_SD={},{}
  perc_exp_3, perc_exp_5, perc_exp_10, perc_exp_20 = {},{},{},{}
  unmutated_singletons, unmutated_perc_exp_3 = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      rep, isotype =l[1],l[3]
      reads_per_isotype = int(l[7])
      range_muts, n_unique_BCRs,sd_muts, total_cluster_size = int(l[5]), int(l[8]),float(l[7]),int(l[9])
      typ = rep+"\t"+isotype
      mut = -1
      clust = l[2]
      if(clust in muts):mut = mean(muts[clust])
      if(typ in total_BCRs_per_isotype):
        total_BCRs_per_isotype[typ] = total_BCRs_per_isotype[typ] + reads_per_isotype
        unique_BCRs_per_isotype[typ] = unique_BCRs_per_isotype[typ] + n_unique_BCRs
      else:
        total_BCRs_per_isotype[typ] = reads_per_isotype
        unique_BCRs_per_isotype[typ] = n_unique_BCRs
      if(total_cluster_size==reads_per_isotype and n_unique_BCRs==1):
        if(typ in singletons_per_isotype):singletons_per_isotype[typ] = singletons_per_isotype[typ]+1
        else:singletons_per_isotype[typ] = 1
        if(mut!=-1 and mut<=4):
          if(typ in unmutated_singletons):unmutated_singletons[typ] = unmutated_singletons[typ]+1
          else:unmutated_singletons[typ] = 1
      if(n_unique_BCRs>=3):
        if(mut!=-1 and mut<=4):
          if(typ in unmutated_perc_exp_3):unmutated_perc_exp_3[typ] = unmutated_perc_exp_3[typ]+reads_per_isotype
          else:unmutated_perc_exp_3[typ]=reads_per_isotype
        if(typ in expanded_clonal_range):
          expanded_clonal_range[typ] = expanded_clonal_range[typ]+[range_muts*1.0/n_unique_BCRs]
          expanded_clonal_SD[typ] = expanded_clonal_SD[typ]+[sd_muts]
          perc_exp_3[typ] = perc_exp_3[typ]+reads_per_isotype
        else:
          expanded_clonal_range[typ]=[range_muts*1.0/n_unique_BCRs]
          expanded_clonal_SD[typ]=[sd_muts]
          perc_exp_3[typ] = reads_per_isotype
        if(n_unique_BCRs>=5):
          if(typ in perc_exp_5):perc_exp_5[typ] = perc_exp_5[typ]+reads_per_isotype
          else:perc_exp_5[typ]=reads_per_isotype
        if(n_unique_BCRs>=10):
          if(typ in perc_exp_10):perc_exp_10[typ] = perc_exp_10[typ]+reads_per_isotype
          else:perc_exp_10[typ]=reads_per_isotype
        if(n_unique_BCRs>=20):
          if(typ in perc_exp_20):perc_exp_20[typ] = perc_exp_20[typ]+reads_per_isotype
          else:perc_exp_20[typ]=reads_per_isotype
  fh.close()
  out = "#sample\tsubsample\tisotype\t% singletons (total BCRs)\t% singletons (unique)\tnumber of expanded clusters\tmean_cluster mutational range (bp/unique sequences)\tSD cluster mutational range\ttotal_BCRs_per_isotype\tunique_BCRs_per_isotype\tn_in_expanded_3_clusters\tn_in_expanded_5_clusters\tn_in_expanded_10_clusters\tn_in_expanded_20_clusters\tunmutated_singletons\tunmutated_perc_expanded_3_clusters\n"
  for typ in total_BCRs_per_isotype:
    expanded_clonal_ranges, expanded_clonal_SDs,length = -1,-1,0
    singletons1,singletons2 = 0,0
    if(typ in expanded_clonal_range):
      length = len(expanded_clonal_range[typ])
      expanded_clonal_ranges = mean(expanded_clonal_range[typ])
      expanded_clonal_SDs = mean(expanded_clonal_SD[typ])
    if(typ in singletons_per_isotype):
      singletons1 = singletons_per_isotype[typ]*100.0/total_BCRs_per_isotype[typ]
      singletons2 = singletons_per_isotype[typ]*100.0/unique_BCRs_per_isotype[typ]
    array = []
    for x in [perc_exp_3, perc_exp_5, perc_exp_10, perc_exp_20,unmutated_singletons, unmutated_perc_exp_3 ]:
      if(typ in x):array.append(x[typ])
      else:array.append(0)
    out=out+sample+"\t"+typ+"\t"+"\t".join(map(str, [singletons1, singletons2, length, expanded_clonal_ranges, expanded_clonal_SDs,total_BCRs_per_isotype[typ],unique_BCRs_per_isotype[typ]]+array))+"\n"
  fh=open(cluster_isotype_expansion_SUBSAMPLE_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_network_parameters_total_subsample(sample,developmental_classification_file, subsample_file, reverse_primer_group,clonality_file_SUBSAMPLED,isotype_usages_SUBSAMPLED,cluster_SUBSAMPLE_information_file,cluster_isotype_expansion_SUBSAMPLE_file,per_cluster_developmental_classification_file):
  Get_raw_files(sample,developmental_classification_file, subsample_file, reverse_primer_group,clonality_file_SUBSAMPLED,isotype_usages_SUBSAMPLED,cluster_SUBSAMPLE_information_file,cluster_isotype_expansion_SUBSAMPLE_file)
  Get_expanded_cluster_summary(cluster_SUBSAMPLE_information_file,cluster_isotype_expansion_SUBSAMPLE_file,sample,developmental_classification_file,per_cluster_developmental_classification_file)
  return()

def Get_raw_files(sample,developmental_classification_file, subsample_file, reverse_primer_group,clonality_file_SUBSAMPLED,isotype_usages_SUBSAMPLED,cluster_SUBSAMPLE_information_file,cluster_isotype_expansion_SUBSAMPLE_file):
  fh=open(cluster_SUBSAMPLE_information_file,"w")
  fh.close()
  fh = open(developmental_classification_file,"r")
  classification_ids,classification_clusters = {},{}
  muts = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id, mut = l[1].split("__")[0], int(l[6])
      muts[id] = mut
  fh.close()
  print "READ DEVELOP FILE"
  fh=open(subsample_file,"r")
  v_sizes,c_sizes,reads_per_sample = {},Tree(),{}
  percent_isotypes,clusters, cluster_types = {},Tree(),Tree()
  index = 0
  muts_per_cluster = {}
  cluster_sizes = {}
  total_cluster_sizes = {}
  for l in fh:
    index=index+1
    if (index>1):
      l=l.strip().split()
      id,freq, id_short,clust = l[2].split("__")[0],map(int, l[2].split("__")[1].split("|")[0].split("_")),l[2].split("__")[0], l[1]
      clust,f = l[1],sum(freq)
      classification = l[2].split("|")[1].split("_")
      subsample = l[0]
      if(subsample in reads_per_sample):reads_per_sample[subsample] = reads_per_sample[subsample]+sum(freq)
      else:reads_per_sample[subsample]=sum(freq)
      clust1=clust+":"+subsample
      if(clust1 in total_cluster_sizes):total_cluster_sizes[clust1] = total_cluster_sizes[clust1] +f
      else:total_cluster_sizes[clust1] = f
       ## overal network params
      c = subsample
      if(c in v_sizes):v_sizes[c] = v_sizes[c]+[f]
      else:v_sizes[c]=[f] 
      c_sizes[c][clust][id].value = 1
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      classes = {}
      for i in nz:
        c = classification[i].split("*")[0]
        if(c.count("TRBC")==1):c = "TRBC"
        if(c in classes):classes[c] = classes[c]+freq[i]
        else:classes[c]=freq[i]
      for c1 in classes:
        c = subsample+"\t"+c1
        if(c in percent_isotypes):percent_isotypes[c] = percent_isotypes[c]+classes[c1]
        else:percent_isotypes[c]=classes[c1]
        cluster_types[subsample][clust][c1].value = 1
        c = subsample+":"+clust+":"+c1
        if(c in cluster_sizes):cluster_sizes[c] = cluster_sizes[c]+f
        else:cluster_sizes[c]=f
        if(c in muts_per_cluster):muts_per_cluster[c] = muts_per_cluster[c]+[muts[id]]*classes[c1]
        else:muts_per_cluster[c]=[muts[id]]*classes[c1]
      clusters[subsample][clust][f][id].value = 1
  fh.close()
  cluster_size_cumul = {}
  for c in cluster_sizes:
    subsample = c.split(":")[0]
    if(cluster_sizes[c]>1):
      if(subsample in cluster_size_cumul):cluster_size_cumul[subsample] = cluster_size_cumul[subsample]+[cluster_sizes[c]]
      else:cluster_size_cumul[subsample]=[cluster_sizes[c]]
  out = "#Id\tIsotype\tN reads isotype\tN reads sampled\t% reads per isotype\tN clusters per isotype\tN clusters\t% clusters per isotype\tN clusters per isotype expanded\tN clusters expanded\t% clusters per isotype expanded\tmean mutations per cluster\tmean mutations per expanded cluster\tsample\n"
  out1 = "#Id\tsubsample\tIsotype\tcluster\tmean_mutations\trange_mutations\tsd mutations\tN reads per isotype\tunique reads per isotype\ttotal cluster size\n"
  ind1 = 0
  for subsample in reads_per_sample:
    expanded_clusters_per_isotype = {}
    for clust in clusters[subsample]:
      if(len(clusters[subsample][clust])>=3):
        expanded_clusters_per_isotype[clust] = 1
      else:
        f = 0
        for f1 in clusters[subsample][clust]:
          f = f+f1*len(clusters[subsample][clust][f1])
          if(f>=3):
            expanded_clusters_per_isotype[clust] = 1
    classes_per_expanded_cluster, classes_per_cluster = {},{}
    muts_per_clusters, muts_per_clusters_expanded = {},{}
    for clust in cluster_types[subsample]:
      for classes in cluster_types[subsample][clust]:
        c = subsample+":"+clust+":"+classes
        clust1=clust+":"+subsample
        if(classes.count("IGH")==1):
          out1 = out1+sample+"\t"+subsample+"\t"+clust+"\t"+classes+"\t"+str(mean(muts_per_cluster[c]))+"\t"+str(max(muts_per_cluster[c])-min(muts_per_cluster[c]))+"\t"+str(np.std(muts_per_cluster[c]))+"\t"+str(cluster_sizes[c])+"\t"+str(len(muts_per_cluster[c]))+"\t"+str(total_cluster_sizes[clust1])+"\n"
          ind1=ind1+1
        if(ind1>500):
          Write_out(out1, cluster_SUBSAMPLE_information_file)
          out1, ind1 ='',0
        if(classes in classes_per_cluster):
          classes_per_cluster[classes] = classes_per_cluster[classes]+1
          muts_per_clusters[classes] = muts_per_clusters[classes]+muts_per_cluster[c]
        else:
          classes_per_cluster[classes]=1
          muts_per_clusters[classes]=muts_per_cluster[c]
        if(clust in expanded_clusters_per_isotype):
          if(classes in classes_per_expanded_cluster):
            classes_per_expanded_cluster[classes]=classes_per_expanded_cluster[classes]+1
            muts_per_clusters_expanded[classes] =muts_per_clusters_expanded[classes]+muts_per_cluster[c]
          else:
            classes_per_expanded_cluster[classes] = 1
            muts_per_clusters_expanded[classes]=muts_per_cluster[c]
    for classes in classes_per_cluster: 
      c1 = subsample+"\t"+classes
      n_expanded = 0
      mutations = -1
      if(classes in classes_per_expanded_cluster):
        n_expanded = classes_per_expanded_cluster[classes]
        mutations = mean(muts_per_clusters_expanded[classes])
      if(c1 in percent_isotypes):
        if(len(expanded_clusters_per_isotype)>0):perc_expanded = n_expanded*100.0/len(expanded_clusters_per_isotype)
        else:perc_expanded = -1
        out=out+c1+"\t"+str(percent_isotypes[c1])+"\t"+str(reads_per_sample[subsample])+"\t"+str(percent_isotypes[c1]*100.0/reads_per_sample[subsample])+"\t"+str(classes_per_cluster[classes])+"\t"+str(len(cluster_types[subsample]))+"\t"+str(classes_per_cluster[classes]*100.0/len(cluster_types[subsample]))+"\t"+str(len(expanded_clusters_per_isotype))+"\t"+str(n_expanded)+"\t"+str(perc_expanded)+"\t"+str(mean(muts_per_clusters[classes]))+"\t"+str(mutations)+"\t"+sample+"\n"
  Write_out(out1, cluster_SUBSAMPLE_information_file)
  out1, ind1 ='',0
  fh=open(isotype_usages_SUBSAMPLED, "w")
  fh.write(out)
  fh.close()
  out="#Id\tIsotype\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\tnumber of unique clusters\ttotal reads per sample\tsample\tmean_cluster_size\tmedian_cluster_size\n"
  for subsample in reads_per_sample:
    c_sizes_classification = []
    for clust in c_sizes[subsample]:
      c_sizes_classification.append(len(c_sizes[subsample][clust]))
    if(len(v_sizes[subsample])>1 and len(c_sizes_classification)>1):
      (vpoints,vvdf)=VDF(v_sizes[subsample])
      (cpoints,cvdf)=VDF(c_sizes_classification)
      vgini, cgini=Gini_index(cpoints,cvdf, vpoints,vvdf)
      cluster_sizes = cluster_size_cumul[subsample]
      cluster_sizes.sort(reverse = True)
      max_pop, max_1_pop =cluster_sizes[0]*100.0/reads_per_sample[subsample], cluster_sizes[1]*100.0/reads_per_sample[subsample]
      mean_cluster_size, median_cluster_size = mean(c_sizes_classification), median(c_sizes_classification)
      #max_pop, max_1_pop = cpoints[len(cpoints)-1]*100.0/sum(cpoints), cpoints[len(cpoints)-2]*100.0/sum(cpoints)
      out = out+subsample+"\tALL\t"+str(sum(v_sizes[subsample]))+"\t"+str(len(v_sizes[subsample]))+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\t"+str(len(c_sizes_classification))+"\t"+str(reads_per_sample[subsample])+"\t"+sample+"\t"+str(mean_cluster_size)+"\t"+str(median_cluster_size)+"\n"
  fh=open(clonality_file_SUBSAMPLED, "w")
  fh.write(out)
  fh.close()
  return()

def Uniq(v):
  C=set(v)
  return list(C)

def VDF (n):
  points=sorted(Uniq(n))
  vdf=[]
  for i in range(0,len(points)):
    vdf.append(n.count(points[i]))
  return (points,vdf)

def Subsample_repertoire(subsample_file, sample_depth, seq_file,samplex,cluster_file):
  n_repeats = 40
  clusters = {}
  fh=open(cluster_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      clusters[l[2].split("__")[0]] = l[1]
  fh.close()
  fh=open(seq_file,"r")
  sequences,id_seqs = [],{}
  for header,sequence in fasta_iterator(fh):
    id = header.split("__")[0]
    if(id in clusters):
      freq,classes = map(int, header.split("__")[1].split("|")[0].split("_")), header.split("|")[1].split("_")
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      for i in nz: 
        for j in range(0,freq[i]):
          sequences.append([sequence,classes[i],id])
  fh.close()
  fh=open(subsample_file,"w")
  fh.write("#repeat\tcluster\tid\ttotal frequency\n")
  fh.close()
  if(sample_depth<len(sequences)):
    classes_dict = {}
    for i in range(0,len(classes)):
      classes_dict[classes[i]] = i
    classes_lab = "_".join(classes)
    for r in range(0,n_repeats):
      rep_id = "REP"+str(r)+"_"+samplex
      print "\t",sample_depth
      sample = random.sample(sequences, sample_depth)
      seqs = Tree()
      out,ind='',0
      for s in sample:
        seqs[s[0]][s[1]][s[2]][ind].value = 1
        ind = ind+1
      for s in seqs:
        f=[0]*len(classes)
        for clas in seqs[s]:
          for header in seqs[s][clas]:
            f[classes_dict[clas]] = f[classes_dict[clas]]+len(seqs[s][clas][header])
            break
        cluster = clusters[header]
        header = header+"__"+"_".join(map(str, f))+"|"+classes_lab
        out=out+rep_id+"\t"+cluster+"\t"+header+"\t"+str(sum(f))+"\n"
      Write_out(out, subsample_file)
      out=''
  return()

def Write_out(out, file):
  fh = open (file,"a")
  fh.write(out)
  fh.close()
  return()

###########################
id = sys.argv[1]
group = sys.argv[2]
input_dir = sys.argv[3]+"ORIENTATED_SEQUENCES/"
output_dir = input_dir
species = sys.argv[4]
reverse_primer_group = "STANDARD"
###########################
id1 = id
###########################
seq_file = input_dir+"NETWORKS/Fully_reduced_"+id+".fasta"
annot_file = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_1_Summary.txt"
annot_file2 = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_8_V-REGION-nt-mutation-statistics.txt"
annot_file3 = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_3_Nt-sequences.txt"
annot_file4 = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_5_AA-sequences.txt"
annot_file7 = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_7_V-REGION-mutation-and-AA-change-table.txt"
cluster_file = input_dir+"NETWORKS/Cluster_identities_"+id+".txt"
developmental_classification_file = output_dir+"ISOTYPER/Developmental_classification_"+id+".txt"
per_cluster_developmental_classification_file = output_dir+"ISOTYPER/Developmental_per_cluster_classification_"+id+".txt"
per_cluster_developmental_classification_network_parameters = output_dir+"ISOTYPER/per_cluster_network_parameters_"+id+".txt"
########################## subsampled
subsample_file = output_dir+"ISOTYPER/Subsample_"+id+".txt"
clonality_file_SUBSAMPLED = output_dir+"ISOTYPER/Clonality_per_sample_SUBSAMPLED_"+id+".txt"
isotype_usages_SUBSAMPLED = output_dir+"ISOTYPER/Isotype_usages_SUBSAMPLED_"+id+".txt"
per_cluster_developmental_classification_network_parameters_subsampled = output_dir+"ISOTYPER/per_cluster_network_parameters_SUBSAMPLED"+id+".txt"
cluster_SUBSAMPLE_information_file = output_dir+"ISOTYPER/Cluster_SUBSAMPLE_information_"+id+".txt"
cluster_isotype_expansion_SUBSAMPLE_file = output_dir+"ISOTYPER/Expansion_cluster_isotype_SUBSAMPLE_"+id+".txt"
###########################
sample = id
########################### 
subsample_depth = 1000

command="STANDARD"

if(command == "STANDARD"):
  Classify_sequences_into_developmental_stages(sample, annot_file, cluster_file, developmental_classification_file,reverse_primer_group)
  Classify_sequences_into_developmental_stages_per_cluster(sample, annot_file, cluster_file,per_cluster_developmental_classification_file,reverse_primer_group)
  ### Subsampled for clonality
  Subsample_repertoire(subsample_file, subsample_depth, seq_file,sample,cluster_file)
  Get_network_parameters_per_classificiation_subsample(sample,per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters_subsampled,reverse_primer_group,cluster_file,developmental_classification_file)
  Get_network_parameters_total_subsample(sample,developmental_classification_file, subsample_file, reverse_primer_group,clonality_file_SUBSAMPLED,isotype_usages_SUBSAMPLED,cluster_SUBSAMPLE_information_file,cluster_isotype_expansion_SUBSAMPLE_file,per_cluster_developmental_classification_file)
  
