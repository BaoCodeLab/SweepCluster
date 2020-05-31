# SweepCluster
SweepCluster is a python library and toolkit for implementation of SNP clustering and significance estimation based on the anchor-extension method.

## Install
### Requirements
  Python 2.7+, but Python 3.7+ is recommended。
  scipy
  numpy 
  pandas  
  vcf  
  scikit-learn
  multiprocessing

####  
    git clone https://github.com/BaoCodeLab/SweepCluster
    cd SweepCluster
    pip install .
    

## Running SweepCluster
### Usage:  SweepCluster.py [-h] {Density,Cluster,Pval,Dbscan}

#### Command line usage                        
    -h, --help          Whow the help message.
    Density             Calculate the SNP density on the genome.
    Cluster             Perform SNP clustering.
    Pval                Calculate the significance p-value of SNP clustering.
    Dbscan              The machine learning python module DBSACN was used to
                        cluster SNPs.

### Calculate SNP density using the sliding window method
#### Usage:  SweepCluster.py Density [-h] -vcf VCF -out OUT  [-scale SCALE]  [-step STEP]  [-win WINDOW]  [-length LENGTH]  

#### {arguments}
    -h, --help      Show the help message and exit.
    -vcf VCF        The vcf file.
    -out OUT        The output file.
    -scale SCALE    The value for normalization of SNP density. The SNP density is calculated as the total number of SNPs in each window normalized to a unit length. If not defined, default value of 1000 will be used.
    -length LENGTH  The total length of the genome. This value is useful for calculation of the SNP density at the boundaries. If not defined, the location of the last SNP in the input file will be used.
    -step STEP      The step size of the sliding window method. If not defined, default value of 300 will be used.
    -win WINDOW     The window size of the sliding window method. If not defined, default value of 2000 will be used.
 


###  Perform SNP clustering based on anchor-extension method
#### Usage:  SweepCluster.py Cluster [-h] -vcf VCF -out OUT -anno ANNO -operon OPERON [-sweep_lg SWEEP_LG]  [-scan_loop SCAN_LOOP]  [-min_num MIN_NUM]  [-max_dist MAX_DIST]

#### {arguments}
    -h, --help            Show the help message and exit.
    -vcf VCF              The vcf file.
    -cluster CLUSTER      The output file of clustering result.
    -snp   SNP            The output file of snps after clustering.
    -anno ANNO            Tab-delimited file containing the annotation of the SNPs.
    -operon OPERON        Tab-delimited file defining the gene operons.
    -sweep_lg SWEEP_LG    The estimated average sweep length in the genome. 
    -scan_loop SCAN_LOOP  The number of iterations the program will perform for merging clusters. Default value is set to 100.
    -min_num MIN_NUM      Minimum number of SNPs per cluster. Default value is set to 2.
    -max_dist MAX_DIST    Maximum inter-SNP distances allowed within a cluster.
    
####   The SNP annotation file should contain four columns, i.e., SNP locations, SNP coding types, Gene name, and Gene ID. SNP coding types include: “CDS_synon” for synonymous SNPs; “CDS_nonSynon” for non-synonymous SNPs; “pseudo-gene” for SNPs in pseudogene; “inter-gene” for inter-genic SNPs.  Among the four types, non-synonymous SNPs should be defined as “CDS_nonSynon” because the keyword “CDS_nonSynon” will be used in SNP clustering. Other types could be defined as what you like. Gene name and Gene ID could be blank if the SNP is inter-genic; or alternatively, could be defined specifically, such as gene1-gene2, indicating the location of the SNP in the inter-genic region between gene1 and gene2. An example looks like below:
       260680    CDS_nonSynon    oppC   gene241
       261541    CDS_synon       oppD   gene242
       261592    CDS_synon       oppD   gene242
       261883    inter-gene      -      gene242-gene243
       261896    inter-gene      -      gene242-gene243
       262321    CDS_synon       oppF   gene243
       
####   The operon file defines the genes in each operon with one line for each gene. Each gene contains five columns, i.e., gene ID, the start location of the gene, the end location of the gene, the orientation of the gene, and the operon ID the gene belongs to. If the gene operons of the genome have not yet defined well, the column of operon ID could be replaced by the gene ID. However, the well-defined gene operons will be helpful for SNP clustering. 
       gene107   121551   122138   +   Operon_30
       gene108   122299   123504   ‐   Operon_30
       gene109   123889   126090   +   Operon_31
       gene110   126356   126913   ‐   Operon_31
       gene111   127281   128687   +   Operon_31
       gene112   128757   129668   ‐   Operon_32

####   The average sweeping length can be simulated and optimized using the script "sweep_lg_simulation.sh" in the package. 
####   The maximum inter-SNP distance is used to optimize the identification of clusters. Multiple gene sweeping events with large distances may occur in the same gene/gene operon. This parameter will split the cluster if any inter-SNP distance is greater than the specified threshold. This value may depend on the specific species. Default value is set to 1000 bp, which should be a good value for bacterial genomes.
    

### Estimate the significance p-value of each cluster based on a gamma distribution model of SNPs
#### Usage:  SweepCluster.py Pval [-h] -cluster CLUSTER -out OUT -rate RATE

#### {arguments}
    -h, --help        Show the help message and exit.
    -cluster CLUSTER  The snp clustering result file.
    -out OUT          The output file.
    -rate RATE        The pre-estimated mutation rate under the null hypothesis that the SNPs are independently and randomly distributed.


### Use python module Dbscan to perform blind SNP clustering without considering SNP annotation
#### Usage:  SweepCluster.py Dbscan [-h] -vcf VCF -eps EPS -min_sample MIN_SAMPLE -out OUT

#### {arguments}
    -h, --help          Show the help message and exit.
    -vcf VCF            The vcf file.
    -eps EPS            ϵ - The neighborhood distance threshold for any two SNPs.
    -min_num MIN_NUM    Minimum number of SNPs in a neighborhood for a SNP to be considered as a core SNP.
    -out OUT            The output file.

    
### The shell script sweep_lg_simulation.sh is a driver program for simulating the effect of different sweep length on SNP clustering results.
We found that many factors influence the clustering results. Sweep length is the most influential one. We evaluate the effect by simulating the dynamics of the number of clusters against the sweep length based on generalized additive model. 
#### Usage:  sweep_lg_simulation.sh -vcf VCF -anno ANNO -operon OPERON -start START -end END -step STEP -delta DELTA [-scan_loop SCAN_LOOP] [-max_dist MAX_DIST]  [-min_num MIN_NUM] 

#### {arguments}
    -vcf VCF              The vcf file
    -anno ANNO            The tab-delimited file containing the annotation of the SNPs.
    -operon OPERON        The tab-delimited file defining the gene operons.
    -start START          The lower-bound of the range of sweep length used for simulation.
    -end END              The upper-bound of the range of sweep length used for simulation.
    -step STEP            The step size of sweep length ranging from START to END in simulation.
    -delta DELTA          The small difference used for approximation of the first and second derivatives in the finite difference method.
    -scan_loop SCAN_LOOP  The number of iterations the clustering program will perform for merging clusters. Default value is set to 100.
    -max_dist MAX_DIST    Maximum inter-SNP distances allowed within a cluster. Default value is set to 1000.
    -min_num MIN_NUM      Minimum number of SNPs per cluster. The default value is set to 2.
 
### The DBSCAN_simulation.py is a driver script for simulating the effect of eps and min_num on SNP clustering using Dbscan method.
#### Usage:  DBSCAN_simulation.py -vcf VCF [-eps_start EPS_START] [-eps_end EPS_END] [-min_start MIN_START] [-min_end MIN_END]
#### {arguments}
    -h, --help             Show the help message and exit.
    -vcf VCF               The vcf file.
    -eps_start EPS_START   The lower-bound of the range of eps. Default value is the minimum of the average distance between SNPS.
    -eps_end EPS_END       The upper-bound of the range of eps. Default value is the maximum of the average distance between SNPS.
    -min_start MIN_START   The lower-bound of the range of min_num. Default value is set to 5.
    -min_end MIN_END       The upper-bound of the range of min_num. Default value is set to 50.

    
