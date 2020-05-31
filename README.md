# SweepCluster
SweepCluster is a python library and toolkit for implementation of SNP clustering and significance estimation based on the anchor-extension method.

## Running SweepCluster
### SweepCluster.py [-h] {Density,Cluster,Pval,Dbscan} 

#### Command line usage                        
    -h, --help          show the help message
    Density             Calculate the SNP density on the genome
    Cluster             Perform SNP clustering
    Pval                Calculate the significance p-value of SNP clustering
    Dbscan              The machine learning python module DBSACN was used to
                        cluster SNPs

### Calculate SNP density using the sliding window method
#### SweepCluster.py Density [-h] -vcf VCF -out OUT  [-scale SCALE]  [-step STEP]  [-win WINDOW]  [-length LENGTH]  

#### {arguments}
    -h, --help      show this help message and exit
    -vcf VCF        the vcf file
    -out OUT        the output file
    -scale SCALE    the value for normalization of SNP density. The SNP density is calculated as the total number of SNPs in each window normalized to a unit length. If not defined, default value of 1000 will be used.
    -length LENGTH  the total length of the genome. This value is useful for calculation of the SNP density at the boundaries. If not defined, the location of the last SNP in the input file will be used.
    -step STEP      the step size of the sliding window method. If not defined, default value of 300 will be used.
    -win WINDOW     the window size of the sliding window method. If not defined, default value of 2000 will be used.
 


### Perform SNP clustering using anchor-extension method
#### SweepCluster.py Cluster [-h] -vcf VCF -out OUT -anno ANNO -operon OPERON [-sweep_lg SWEEP_LG]  [-scan_loop SCAN_LOOP]  [-min_num MIN_NUM]  [-max_dist MAX_DIST]

#### {arguments}
    -h, --help            show this help message and exit
    -vcf VCF              the vcf file
    -cluster CLUSTER      the output file of clustering result
    -snp   SNP            the output file of snps after clustering
    -anno ANNO            tab-delimited file containing the annotations for the SNPs
    -operon OPERON        tab-delimited file defining the gene operons, 
    -recomb_lg RECOMB_LG  the estimated average recombination tract length for the genome. This value should be estimated before the running of SNP clustering, using ClonalFrame or other algorithms. 
    -scan_loop SCAN_LOOP   the number of scanning times the program will perform for the SNPs. The higher number of scanning will include as many as qualified but ambiguous boundary SNPs into the clusters.  Default value is set to 100, which has been able to perform well.
    -min_num MIN_NUM      minimum number of SNPs per cluster. The default value is set to 2, meaning that each cluster should contain at least 2 SNPs.
     -max_dist MAX_DIST   maximum inter-SNP distance allowed within a cluster. This parameter is used to optimize the identified clusters. The SNPs falling into the same gene/gene operon are initially defined to be in the same cluster, which may be longer than the estimated recombination tract length “recomb_lg”. However, multiple gene sweeping events may occur in the same gene/gene operon.  This parameter will split the cluster if any inter-SNP distance is greater than the max_dist_cluster. This value may depend on the specific species. Default value is set to 4000 bp, which should be a good value for bacterial genomes of length ~ 1-3 Mbp.
    
####    The SNP annotation file should contains four columns, i.e., SNP locations, SNP coding types, the gene name where the SNP is located, and the corresponding gene ID. SNP coding types include: “CDS_synon” for synonymous SNPs; “CDS_nonSynon” for non-synonymous SNPs; “pseudo-gene” for SNPs in pseudogene; “inter-gene” for inter-genic SNPs.  Among the four types, non-synonymous SNPs should be defined as “CDS_nonSynon” because the keyword “CDS_nonSynon” will be used in SNP clustering. Other types could be defined as what you like. gene name and gene ID could be blank if the SNP is inter-genic; or alternatively, could be defined specifically, such as -gene1-gene2 indicating the presence of the SNP in the inter-genic region between gene1 and gene2.An example looks like below:260680  CDS_nonSynon oppD  AP53_241.The first row should be Loc  Type  Gene_Name Gene_ID
####   The operon file defines the genes in each operon with one line for each gene. Each gene contains five columns for gene ID, the start location of the gene, the end location of the gene, the orientation of the gene, and the operon name the gene belongs to. If the gene operons of the genome have not yet defined well, the column of operon name could be replaced by the gene ID. However, the well-defined gene operons will be helpful for SNP clustering. An example looks like below:AP53_107  121551  122138  + Operon_30.
####    The average recombination tract length can be estimated using the independent tools, such as ClonalFrame, or can be optimized using the script "recomb_lg_simulation.sh" in the package. 
####    The maximum inter-SNP distance is used to optimize the identified clusters. The SNPs falling into the same gene/gene operon are initially defined to be in the same cluster, which may be longer than the estimated recombination tract length “recomb_lg”. However, multiple gene sweeping events may occur in the same gene/gene operon.  This parameter will split the cluster if any inter-SNP distance is greater than the max_dist. This value may depend on the specific species. Default value is set to 4000 bp, which should be a good value for bacterial genomes of length ~ 1-3 Mbp.

### subcommands[Pval]:
#### usage: Maketest.py Pval [-h] -cluster CLUSTER -out OUT -rate RATE

#### {optional arguments}
    -h, --help        show this help message and exit
    -cluster CLUSTER  the snp clustering file
    -out OUT          the output file .
    -rate RATE        it is a pre-estimated average mutation rate under the null
                    hypothesis that the SNPs are independently and randomly
                    distributed.



### subcommands[Dbscan]:
#### usage: Maketest.py Dbscan [-h] -vcf VCF -eps EPS -min_sample MIN_SAMPLE -out
                          OUT

#### optional arguments:
    -h, --help            show this help message and exit
    -vcf VCF              the vcf file
    -eps EPS              ϵ - neighborhood distance threshold, and the sample is
                        more than ϵ sample points is not in the neighborhood.
    -min_sample MIN_SAMPLE
                        Sample points to be a core object need ϵ -
                        neighborhood sample threshold.
    -out OUT              the output file


## The recomb_lg_simulation.sh is used to test the effect of different recombination length SNP clustering results
We find that there are many factors that influence the clustering results.Recombination length has the greatest effect, which directly affects the number of cluster
### Useage:./recomb_lg_simulation.sh [arguments]
#### {optional arguments}
    -vcf:     the vcf file
    -anno:    tab-delimited file containing the annotations for the SNPs, with four columns, SNP locations, SNP coding types, the gene name where the SNP is located, and the corresponding gene ID. SNP coding types include: “CDS_synon” for synonymous SNPs; “CDS_nonSynon” for non-synonymous SNPs; “pseudo-gene” for SNPs in pseudogene; “inter-gene” for inter-genic SNPs.  Among the four types, non-synonymous SNPs should be defined as “CDS_nonSynon” because the keyword “CDS_nonSynon” will be used in SNP clustering. Other types could be defined as what you like. gene name and gene ID could be blank if the SNP is inter-genic; or alternatively, could be defined specifically, such as -gene1-gene2 indicating the presence of the SNP in the inter-genic region between gene1 and gene2.An example looks like below:260680  CDS_nonSynon oppD  AP53_241.The first row should be Loc  Type  Gene_Name Gene_ID
    -operon:    tab-delimited file defining the gene operons, with one line for each gene. Each gene contains five columns for gene ID, the start location of the gene, the end location of the gene, the orientation of the gene, and the operon name the gene belongs to. If the gene operons of the genome have not yet defined well, the column of operon name could be replaced by the gene ID. However, the well-defined gene operons will be helpful for SNP clustering. An example looks like below:AP53_107  121551  122138  + Operon_30.
    -start:     Recombination length used at the beginning of clustering.
    -end :      Recombination length used at the ending of clustering.
    -size :     The step size of recombination length.
    -interv:    intervalue is used to obtain the first and second derivatives of the clustering model by using the finite difference method.
    -scan_loop:   the number of scanning times the program will perform for the SNPs. The higher number of scanning will include as many as qualified but ambiguous boundary SNPs into the clusters.  Default value is set to 100, which has been able to perform well.
    -max_dist:    maximum inter-SNP distance allowed within a cluster. This parameter is used to optimize the identified clusters. The SNPs falling into the same gene/gene operon are initially defined to be in the same cluster, which may be longer than the estimated recombination tract length “recomb_lg”. However, multiple gene sweeping events may occur in the same gene/gene operon.  This parameter will split the cluster if any inter-SNP distance is greater than the max_dist_cluster. This value may depend on the specific species. Default value is set to 4000 bp, which should be a good value for bacterial genomes of length ~ 1-3 Mbp.
    -min_num:   minimum number of SNPs per cluster. The default value is set to 2, meaning that each cluster should contain at least 2 SNPs.
 
## The DBSCAN_simulation.py is used to test the effect of different eps and minimum samples SNP clustering results
### usage: python DBSCAN_simulation.py [arguments]
#### {optional arguments}
    -h, --help            show this help message and exit
    -vcf VCF              the vcf file
    -eps_start EPS_START  the start value of eps,The default value is the
                        minimum of the average distance between SNPS
    -eps_end EPS_END      the end value of eps,The default value is the maximum
                        of the average distance between SNPS
    -sample_start SAMPLE_START
                        the start value of eps,The default value is 5
    -sample_end SAMPLE_END
                        the end value of eps,The default value is 50

    
