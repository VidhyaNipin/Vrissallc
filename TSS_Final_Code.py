import pyranges as pr
import pandas as pd
import pybedtools
from fuc import pybed
 
 # *********To Generate TSS Sites from the known gene ref file********** #

# gr = pr.read_gtf("/Users/neha/Desktop/Gene_Mapping/Reference_genomes/hg38.knownGene.gtf") #read the known gene ref file using pyranges

# print(gr.columns)

# gr = gr[["Chromosome",'Feature', "Start", 'End', 'Strand','gene_id']] # Filter the columns

# print(gr)

# gr = gr.features.tss() #to find the tss site using pyranges


# s = (gr.Chromosome == "chr21") #Filter Chr21

# gr[s].df.to_csv('tss_sites_knowngene.csv', header=True, index=False) # Saving tss sites to a csv file


#  # *********To convert csv file to bed file  ********** #

# df = pd.read_csv('tss_sites_knowngene.csv', usecols = ['Chromosome','Start', 'End','Strand',"gene_id"])
# df.columns = ['Chromosome', 'Start', 'End', 'Strand', "gene_id"]

# bf = pybed.BedFrame.from_frame(meta=[], data=df)
# bf.to_file('tss_sites_knowngene.bed')


#  # *********To extend TSS sites 5kb upstream & 5kb downstream ********** #


tss_sites = pybedtools.BedTool('tss_sites_knowngene.bed') #get the tss sites bed file

# tss_extend_5kb= tss_sites.slop(b=5000, genome='hg38', output='tss_extend-5kb.bed') #extend 5kb upstream and downstream

#  # *********To find intersection of sample with 5kb upstream & 5kb downstream region of tss  ********** #

# bt_sample = pybedtools.BedTool("BT_sample.bed")   #read sample file


# #Fag - wo :- Write the original A and B entries plus the number of base pairs of overlap between the 
# #two features. Only A features with overlap are reported
# sample_intersect_5kb = tss_extend_5kb.intersect(bt_sample, wo=True, output="sample_intersect_5kb_wo.bed")

## **********Calculating the bin from TSS site Upstream and Downstream ## *********
# tss_Minus_1000 = tss_sites.slop(l=1000, r=0, s=True, genome='hg38', output="TSSMinus1000.bed")
# tss_Minus_2000 = tss_Minus_1000.slop(l=1000, r=-1000, s=True, genome='hg38', output="TSSMinus2000.bed")
# tss_Minus_3000 = tss_Minus_2000.slop(l=1000, r=-1000, s=True, genome='hg38', output="TSSMinus3000.bed")
# tss_Minus_4000 = tss_Minus_3000.slop(l=1000, r=-1000, s=True, genome='hg38', output="TSSMinus4000.bed")
# tss_Minus_5000 = tss_Minus_4000.slop(l=1000, r=-1000, s=True, genome='hg38', output="TSSMinus5000.bed")

# tss_plus_1000 = tss_sites.slop(l=0, r=1000, s=True, genome='hg38', output="TSSPlus1000.bed")
# tss_plus_2000 = tss_plus_1000.slop(l=-1000, r=1000, s=True, genome='hg38', output="TSSPlus2000.bed")
# tss_plus_3000 = tss_plus_2000.slop(l=-1000, r=1000, s=True, genome='hg38', output="TSSPlus3000.bed")
# tss_plus_4000 = tss_plus_3000.slop(l=-1000, r=1000, s=True, genome='hg38', output="TSSPlus4000.bed")
# tss_plus_5000 = tss_plus_4000.slop(l=-1000, r=1000, s=True, genome='hg38', output="TSSPlus5000.bed")

## **********Calculating the coverage count of reads from TSS site and Upstream/Downstream region in each bin ## *********

# bt_sample = pybedtools.BedTool("BT_sample.bed")
# bt_sample1 = pr.read_bed("BT_sample.bed")

# TSSMinus1000_count = tss_Minus_1000.coverage(bt_sample, counts=True, output="TSSMinus1000_coverage.bed")
# TSSMinus2000_count = tss_Minus_2000.coverage(bt_sample, counts=True, output="TSSMinus2000_coverage.bed")
# TSSMinus3000_count = tss_Minus_3000.coverage(bt_sample, counts=True, output="TSSMinus3000_coverage.bed")
# TSSMinus4000_count = tss_Minus_4000.coverage(bt_sample, counts=True, output="TSSMinus4000_coverage.bed")
# TSSMinus5000_count = tss_Minus_5000.coverage(bt_sample, counts=True, output="TSSMinus5000_coverage.bed")

# TSSPlus1000_count = tss_plus_1000.coverage(bt_sample, counts=True, output="TSSPlus1000_coverage.bed")
# TSSplus2000_count = tss_plus_2000.coverage(bt_sample, counts=True, output="TSSPlus2000_coverage.bed")
# TSSPlus3000_count = tss_plus_3000.coverage(bt_sample, counts=True, output="TSSPlus3000_coverage.bed")
# TSSPlus4000_count = tss_plus_4000.coverage(bt_sample, counts=True, output="TSSPlus4000_coverage.bed")
# TSSPlus5000_count = tss_plus_5000.coverage(bt_sample, counts=True, output="TSSPlus5000_coverage.bed")

## **********Calculating the coverage fraction of reads from TSS site and Upstream/Downstream region in each bin ## *********

# bt_sample = pybedtools.BedTool("BT_sample.bed")
# bt_sample1 = pr.read_bed("BT_sample.bed")

# TSSMinus1000_count = tss_Minus_1000.coverage(bt_sample, output="TSSMinus1000_coverage_f.bed")
# TSSMinus2000_count = tss_Minus_2000.coverage(bt_sample, output="TSSMinus2000_coverage_f.bed")
# TSSMinus3000_count = tss_Minus_3000.coverage(bt_sample, output="TSSMinus3000_coverage_f.bed")
# TSSMinus4000_count = tss_Minus_4000.coverage(bt_sample, output="TSSMinus4000_coverage_f.bed")
# TSSMinus5000_count = tss_Minus_5000.coverage(bt_sample, output="TSSMinus5000_coverage_f.bed")

# TSSPlus1000_count = tss_plus_1000.coverage(bt_sample, output="TSSPlus1000_coverage_f.bed")
# TSSplus2000_count = tss_plus_2000.coverage(bt_sample, output="TSSPlus2000_coverage_f.bed")
# TSSPlus3000_count = tss_plus_3000.coverage(bt_sample, output="TSSPlus3000_coverage_f.bed")
# TSSPlus4000_count = tss_plus_4000.coverage(bt_sample, output="TSSPlus4000_coverage_f.bed")
# TSSPlus5000_count = tss_plus_5000.coverage(bt_sample, output="TSSPlus5000_coverage_f.bed")

### ******* Coverage Total Sum Upstream ***********

# gr_1 = pr.read_bed("TSSMinus1000_coverage.bed")
# gr_1.df.to_csv("TSSMinus1000_coverage.csv")
# # print(gr_1)

# coverage_tssminus1000 = pd.read_csv("TSSMinus1000_coverage.csv", usecols = [6])
# # # print(coverage_tssminus1000)
# list_tssminus1000 = coverage_tssminus1000["Strand"].to_list()
# # print(list_tssminus1000)
# list_tssminus1000 = (sum(list_tssminus1000))
# print(list_tssminus1000)

# gr_2 = pr.read_bed("TSSMinus2000_coverage.bed")
# gr_2.df.to_csv("TSSMinus2000_coverage.csv")
# # print(gr_2)

# coverage_tssminus2000 = pd.read_csv("TSSMinus2000_coverage.csv", usecols = [6])
# # print(coverage_tssminus2000)
# list_tssminus2000 = coverage_tssminus2000["Strand"].to_list()
# # print(list_tssminus2000)
# list_tssminus2000 = (sum(list_tssminus2000))
# print(list_tssminus2000)

# gr_3 = pr.read_bed("TSSMinus3000_coverage.bed")
# gr_3.df.to_csv("TSSMinus3000_coverage.csv")
# # # print(gr_3)

# coverage_tssminus3000 = pd.read_csv("TSSMinus3000_coverage.csv", usecols = [6])
# # print(coverage_tssminus3000)
# list_tssminus3000 = coverage_tssminus3000["Strand"].to_list()
# # print(list_tssminus3000)
# list_tssminus3000 = (sum(list_tssminus3000))
# print(list_tssminus3000)

# gr_4 = pr.read_bed("TSSMinus4000_coverage.bed")
# gr_4.df.to_csv("TSSMinus4000_coverage.csv")
# # print(gr_4)

# coverage_tssminus4000 = pd.read_csv("TSSMinus4000_coverage.csv", usecols = [6])
# # print(coverage_tssminus4000)
# list_tssminus4000 = coverage_tssminus4000["Strand"].to_list()
# # print(list_tssminus4000)
# list_tssminus4000 = (sum(list_tssminus4000))
# print(list_tssminus4000)

# gr_5 = pr.read_bed("TSSMinus5000_coverage.bed")
# gr_5.df.to_csv("TSSMinus5000_coverage.csv")
# # print(gr_5)

# coverage_tssminus5000 = pd.read_csv("TSSMinus5000_coverage.csv", usecols = [6])
# # print(coverage_tssminus5000)
# list_tssminus5000 = coverage_tssminus5000["Strand"].to_list()
# # print(list_tssminus5000)
# list_tssminus5000 = (sum(list_tssminus5000))
# print(list_tssminus5000)

Upstream_List = [86302, 55212, 33950, 45580, 48944]
print(Upstream_List)


### ******* Coverage Total Sum Downstream ***********

# gr_1_d = pr.read_bed("TSSPlus1000_coverage.bed")
# gr_1_d.df.to_csv("TSSPlus1000_coverage.csv")
# print(gr_1_d)

# coverage_tssplus1000 = pd.read_csv("TSSPlus1000_coverage.csv", usecols = [6])
# print(coverage_tssplus1000)
# list_tssplus1000 = coverage_tssplus1000["Strand"].to_list()
# print(list_tssplus1000)
# list_tssplus1000 = (sum(list_tssplus1000))
# print(list_tssplus1000)

# gr_2_d = pr.read_bed("TSSPlus2000_coverage.bed")
# gr_2_d.df.to_csv("TSSPlus2000_coverage.csv")
# print(gr_2_d)

# coverage_tssplus2000 = pd.read_csv("TSSPlus2000_coverage.csv", usecols = [6])
# print(coverage_tssplus2000)
# list_tssplus2000 = coverage_tssplus2000["Strand"].to_list()
# print(list_tssplus2000)
# list_tssplus2000 = (sum(list_tssplus2000))
# print(list_tssplus2000)

# gr_3_d = pr.read_bed("TSSPlus3000_coverage.bed")
# gr_3_d.df.to_csv("TSSPlus3000_coverage.csv")
# print(gr_3_d)

# coverage_tssplus3000 = pd.read_csv("TSSPlus3000_coverage.csv", usecols = [6])
# # print(coverage_tssplus3000)
# list_tssplus3000 = coverage_tssplus3000["Strand"].to_list()
# print(list_tssplus3000)
# list_tssplus3000 = (sum(list_tssplus3000))
# print(list_tssplus3000)

# gr_4_d = pr.read_bed("TSSPlus4000_coverage.bed")
# gr_4_d.df.to_csv("TSSPlus4000_coverage.csv")
# print(gr_4_d)

# coverage_tssplus4000 = pd.read_csv("TSSPlus4000_coverage.csv", usecols = [6])
# # print(coverage_tssplus4000)
# list_tssplus4000 = coverage_tssplus4000["Strand"].to_list()
# print(list_tssplus4000)
# list_tssplus4000 = (sum(list_tssplus4000))
# print(list_tssplus4000)

# gr_5_d = pr.read_bed("TSSPlus5000_coverage.bed")
# gr_5_d.df.to_csv("TSSPlus5000_coverage.csv")
# print(gr_5_d)

# coverage_tssplus5000 = pd.read_csv("TSSPlus5000_coverage.csv", usecols = [6])
# # print(coverage_tssplus5000)
# list_tssplus5000 = coverage_tssplus5000["Strand"].to_list()
# print(list_tssplus5000)
# list_tssplus5000 = (sum(list_tssplus5000))
# print(list_tssplus5000)

Downstream_List = [35050, 38668, 51121, 46646, 50926]
print(Downstream_List)

Final_List_Count = Upstream_List + Downstream_List
print ("Concatenated list: " + str(Final_List_Count)) 









