import argparse
from library.Density import density
from library.cluster import cluster
from library.Pval import Pval
from library.Dbscan import DBSCA


if __name__=="__main__":
    ### Defines the parameters of the subprocess for calculating SNP density
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="subcommands",description="Density or Cluster or Pval or Dbscan",help="config subscommand")
    
    density_parser=subparsers.add_parser("Density",help="Calculate SNP density using the sliding window method") 
    density_parser.add_argument("-vcf", required=True, help="the vcf file")
    density_parser.add_argument("-out", required=True, help="the output file")
    density_parser.add_argument("-scale", help="The value for normalization of SNP density. The SNP density is calculated as the total number of SNPs in each window normalized to a unit length. If not defined, default value of 1000 will be used.")
    density_parser.add_argument("-length", help="The total length of the genome. This value is useful for calculation of the SNP density at the boundaries. If not defined, the location of the last SNP in the input file will be used.")
    density_parser.add_argument("-step", help="The step size of the sliding window method. If not defined, default value of 300 will be used.")
    density_parser.add_argument("-win", help="The window size of the sliding window method. If not defined, default value of 2000 will be used.")
    density_parser.set_defaults(func=density)

    
    ### Defines the parameters of the subprocess for SNP clustering
    cluster_parser=subparsers.add_parser("Cluster",help="Perform SNPs clustering based on anchor-extension method")
    cluster_parser.add_argument("-vcf", required=True, help="The vcf file")
    cluster_parser.add_argument("-cluster", required=True, help="The output file of clustering result")
    cluster_parser.add_argument("-snp", required=True, help="The output file of snps after clustering")
    cluster_parser.add_argument("-anno", required=True,help="Tab-delimited file containing the annotations for the SNPs")
    cluster_parser.add_argument("-operon", required=True,help="Tab-delimited file defining the gene operons")
    cluster_parser.add_argument("-sweep_lg", required=True,help="The estimated average sweep length in the genome")
    cluster_parser.add_argument("-scan_loop",help="The number of iterations the program will perform for merging clusters. Default value is set to 100.")
    cluster_parser.add_argument("-max_dist",help="Maximum inter-SNP distances allowed within a cluster. Default value is set to 4000 bp")
    cluster_parser.add_argument("-min_num", help="Minimum number of SNPs per cluster. Default value is set to 2.")
    cluster_parser.set_defaults(func=cluster)

    
    ### Defines the parameters of the subprocess for calculating SNPs cluster p-value
    Pval_parser=subparsers.add_parser("Pval",help="Estimate the significance p-value of each cluster based on a gamma distribution model of SNPs")   
    Pval_parser.add_argument("-cluster", required=True, help="The snp clustering result file.")
    Pval_parser.add_argument("-out", required=True,help="the output file")
    Pval_parser.add_argument("-rate",required=True,The pre-estimated mutation rate under the null hypothesis that the SNPs are independently and randomly distributed.")
    Pval_parser.set_defaults(func=Pval)

    
    ###### Defines the parameters of the subprocess for clusterting the SNPs by DBSCAN
    DBSCAN_parser=subparsers.add_parser("Dbscan",help="Use python module Dbscan to perform blind SNP clustering without considering SNP annotation")
    DBSCAN_parser.add_argument("-vcf", required=True, help="The vcf file")
    DBSCAN_parser.add_argument("-eps",required=True,help="ϵ - The neighborhood distance threshold for any two SNPs.")
    DBSCAN_parser.add_argument("-min_num",required=True,help="Minimum number of SNPs in a neighborhood for a SNP to be considered as a core SNP.")
    DBSCAN_parser.add_argument("-out",required=True,help="the output file")
    DBSCAN_parser.set_defaults(func=DBSCA)
    args = parser.parse_args()
    args.func(args)
