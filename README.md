# SNP-cluster
This is an algorithm program about SNP clustering
## The SNP_cluster.py contains a series of functions for SNP clustering
### density():
For simple density analysis of mutant SNPs on the genome.                                                                       
### SNP_anno():
Store the annotations for all SNPs and annotateing the target SNP list. 
### Operon_manage(): 
handle the operon file including:profile the gene operons,define the range of each gene operon,re-define the gene operons for target SNPs due to the fact that some operons may include several genes,profile the snp list for each operon.
### clust_initaliz():
initialize the targetClust_dic
### SNP_dist():
Calculate the distance between SNPs
### cluster_merge(): 
Merge clusters with neighboring clusters,sorte the clusters.
### cluster_out():
Output the result of SNP clustering
### cluster_Pval():
Statistical test was conducted on the clustering of SNPs to determine whether it had statistical significance.
## The MakeTest.py is the main program
### Useage:python MakeTest.py [arguments]
### Parameters: [required]
#### -AnalTpye
The type of SNP analysis to be performed:density or cluster or Pval
#### -input
the vcf file or clustering file.
#### -anno
tab-delimited file containing the annotations for the SNPs, with four columns, SNP locations, SNP coding types, the gene name where the SNP is located, and the corresponding gene ID. SNP coding types include: “CDS_synon” for synonymous SNPs; “CDS_nonSynon” for non-synonymous SNPs; “pseudo-gene” for SNPs in pseudogene; “inter-gene” for inter-genic SNPs.  Among the four types, non-synonymous SNPs should be defined as “CDS_nonSynon” because the keyword “CDS_nonSynon” will be used in SNP clustering. Other types could be defined as what you like. gene name and gene ID could be blank if the SNP is inter-genic; or alternatively, could be defined specifically, such as -gene1-gene2 indicating the presence of the SNP in the inter-genic region between gene1 and gene2.An example looks like below:260680  CDS_nonSynon oppD  AP53_241.The first row should be Loc  Type  Gene_Name Gene_ID
#### -operon
tab-delimited file defining the gene operons, with one line for each gene. Each gene contains five columns for gene ID, the start location of the gene, the end location of the gene, the orientation of the gene, and the operon name the gene belongs to. If the gene operons of the genome have not yet defined well, the column of operon name could be replaced by the gene ID. However, the well-defined gene operons will be helpful for SNP clustering. An example looks like below:AP53_107  121551  122138  + Operon_30.
#### -recomb_lg
 the estimated average recombination tract length for the genome. This value should be estimated before the running of SNP clustering, using ClonalFrame or other algorithms. 
#### -output
the output file
#### -rate
It is a pre-estimated average mutation rate under the null hypothesis that the SNPs are independently and randomly distributed.
### Parameters: [optional]:
#### -total: 
the total length of the genome. This value is useful for calculation of the SNP density at the boundaries. If not defined, the location of the last SNP in the input file will be used.
#### -step:
 the step size of the sliding window method. If not defined, default value of 300 will be used.
 #### -ws:
 the window size of the sliding window method. If not defined, default value of 2000 will be used.
 #### -scale:
  the value for normalization of the SNP density. The SNP density is calculated as the total number of SNPs in each window normalized to a unit length. If not defined, default value of 1000 will be used.
  #### -scan_loop:
  the number of scanning times the program will perform for the SNPs. The higher number of scanning will include as many as qualified but ambiguous boundary SNPs into the clusters.  Default value is set to 100, which has been able to perform well.
  #### -max_dist：
  maximum inter-SNP distance allowed within a cluster. This parameter is used to optimize the identified clusters. The SNPs falling into the same gene/gene operon are initially defined to be in the same cluster, which may be longer than the estimated recombination tract length “recomb_lg”. However, multiple gene sweeping events may occur in the same gene/gene operon.  This parameter will split the cluster if any inter-SNP distance is greater than the max_dist_cluster. This value may depend on the specific species. Default value is set to 4000 bp, which should be a good value for bacterial genomes of length ~ 1-3 Mbp.
  #### -min_num：
   minimum number of SNPs per cluster. The default value is set to 2, meaning that each cluster should contain at least 2 SNPs.
##The Rscript is used to graph the results of SNP clustering,is referenced in recomb_lg_test.sh
##The recomb_lg_test.sh is used to test the effect of different recombination length SNP clustering results
We find that there are many factors that influence the clustering results.Recombination length has the greatest effect, which directly affects the number of cluster
### Useage:./recomb_lg_test.sh [arguments]
### Parameters: [required]
#### -vcf:
the vcf file
#### -anno:
tab-delimited file containing the annotations for the SNPs, with four columns, SNP locations, SNP coding types, the gene name where the SNP is located, and the corresponding gene ID. SNP coding types include: “CDS_synon” for synonymous SNPs; “CDS_nonSynon” for non-synonymous SNPs; “pseudo-gene” for SNPs in pseudogene; “inter-gene” for inter-genic SNPs.  Among the four types, non-synonymous SNPs should be defined as “CDS_nonSynon” because the keyword “CDS_nonSynon” will be used in SNP clustering. Other types could be defined as what you like. gene name and gene ID could be blank if the SNP is inter-genic; or alternatively, could be defined specifically, such as -gene1-gene2 indicating the presence of the SNP in the inter-genic region between gene1 and gene2.An example looks like below:260680  CDS_nonSynon oppD  AP53_241.The first row should be Loc  Type  Gene_Name Gene_ID
#### -operon:
tab-delimited file defining the gene operons, with one line for each gene. Each gene contains five columns for gene ID, the start location of the gene, the end location of the gene, the orientation of the gene, and the operon name the gene belongs to. If the gene operons of the genome have not yet defined well, the column of operon name could be replaced by the gene ID. However, the well-defined gene operons will be helpful for SNP clustering. An example looks like below:AP53_107  121551  122138  + Operon_30.
#### -rate:
It is a pre-estimated average mutation rate under the null hypothesis that the SNPs are independently and randomly distributed.
#### -start:
Recombination length used at the beginning of clustering.
#### -end :
Recombination length used at the ending of clustering.
#### -size :
The step size of recombination length.
#### -interv:
intervalue is used to obtain the first and second derivatives of the clustering model by using the finite difference method.
#### -scan_loop:
the number of scanning times the program will perform for the SNPs. The higher number of scanning will include as many as qualified but ambiguous boundary SNPs into the clusters.  Default value is set to 100, which has been able to perform well.
#### -max_dist:
maximum inter-SNP distance allowed within a cluster. This parameter is used to optimize the identified clusters. The SNPs falling into the same gene/gene operon are initially defined to be in the same cluster, which may be longer than the estimated recombination tract length “recomb_lg”. However, multiple gene sweeping events may occur in the same gene/gene operon.  This parameter will split the cluster if any inter-SNP distance is greater than the max_dist_cluster. This value may depend on the specific species. Default value is set to 4000 bp, which should be a good value for bacterial genomes of length ~ 1-3 Mbp.
#### -min_num:
 minimum number of SNPs per cluster. The default value is set to 2, meaning that each cluster should contain at least 2 SNPs.
