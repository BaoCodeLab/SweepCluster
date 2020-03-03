import argparse
from library.Density import density
from library.cluster import cluster
from library.Pval import Pval
from library.Dbscan import DBSCA


if __name__=="__main__":
    ###Defines the parameters of the subprocess for calculating the SNPs density
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="subcommands",description="density or cluster or Pval or DBSCAN",help="config subscommand help")
    density_parser=subparsers.add_parser("density",help="Calculate the density of SNPs")
    density_parser.add_argument("-vcf", required=True, help="the vcf file")
    density_parser.add_argument("-out", required=True, help="the output file")
    density_parser.add_argument("-scale", help="the value for normalization of the SNP density")
    density_parser.add_argument("-total", help="the total length of the genome")
    density_parser.add_argument("-step", help="the step size of the sliding window method")
    density_parser.add_argument("-ws", help="the window size of the sliding window method")
    density_parser.set_defaults(func=density)

    ###Defines the parameters of the subprocess for  the SNPs clustering
    cluster_parser=subparsers.add_parser("cluster",help="Perform SNPs clustering")
    cluster_parser.add_argument("-vcf", required=True, help="the vcf file")
    cluster_parser.add_argument("-cluster", required=True, help="the output file of clustering result")
    cluster_parser.add_argument("-snp", required=True, help="the output file of snp after clustering")
    cluster_parser.add_argument("-anno", required=True,help="tab-delimited file containing the annotations for the SNPs with four columns")
    cluster_parser.add_argument("-operon", required=True,help="tab-delimited file defining the gene operons, with one line for each gene.")
    cluster_parser.add_argument("-recomb_lg", required=True,help="the estimated average recombination tract length for the genome.")
    cluster_parser.add_argument("-scan_loop",help="the number of scanning times the program will perform for the SNPs.Default value is set to 100.")
    cluster_parser.add_argument("-max_dist",help="maximum inter-SNP distance allowed within a cluster.Default value is set to 4000 bp.")
    cluster_parser.add_argument("-min_num", help=" minimum number of SNPs per cluster.The default value is set to 2.")
    cluster_parser.set_defaults(func=cluster)

    ###Defines the parameters of the subprocess for calculating the SNPs cluster P value
    Pval_parser=subparsers.add_parser("Pval",help="Calculate the P value of SNP clustering")
    Pval_parser.add_argument("-cluster", required=True, help="the snp clustering file")
    Pval_parser.add_argument("-out", required=True,help="the output file .")
    Pval_parser.add_argument("-rate",required=True,help="it is a pre-estimated average mutation rate under the null hypothesis that the SNPs are independently and randomly distributed.")
    Pval_parser.set_defaults(func=Pval)

    ######Defines the parameters of the subprocess for clusterting the SNPs by DBSCAN
    DBSCAN_parser=subparsers.add_parser("DBSCAN",help="The DBSACN method of machine learning was used to cluster SNPS")
    DBSCAN_parser.add_argument("-vcf", required=True, help="the vcf file")
    DBSCAN_parser.add_argument("-eps",required=True,help="ϵ - neighborhood distance threshold, and the sample is more than ϵ sample points is not in the neighborhood.")
    DBSCAN_parser.add_argument("-min_sample",required=True,help="Sample points to be a core object need ϵ - neighborhood sample threshold.")
    DBSCAN_parser.add_argument("-out",required=True,help="the output file")
    DBSCAN_parser.set_defaults(func=DBSCA)
    args = parser.parse_args()
    args.func(args)
