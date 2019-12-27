import argparse
import sys
from snp_cluster import density
from snp_cluster import SNP_anno
from snp_cluster import Operon_manage
from snp_cluster import clust_initalize
from snp_cluster import SNP_dist
from snp_cluster import cluster_merge
from snp_cluster import cluster_out
from snp_cluster import cluster_Pval
if __name__=="__main__":
    pars = argparse.ArgumentParser()
    pars.add_argument("-AnalType", required=True, help="The type of SNP analysis to be performed:density or cluster or Pval")
    pars.add_argument("-input", required=True, help="the vcf file or cluster file")
    pars.add_argument("-anno", required=True,help="tab-delimited file containing the annotations for the SNPs with four columns")
    pars.add_argument("-operon", required=True,help="tab-delimited file defining the gene operons, with one line for each gene.")
    pars.add_argument("-recomb_lg", required=True,help="the estimated average recombination tract length for the genome.")
    pars.add_argument("-output", required=True,help="the output file .")
    pars.add_argument("-rate",required=True,help="it is a pre-estimated average mutation rate under the null hypothesis that the SNPs are independently and randomly distributed.")
    ### Parameters optional:
    pars.add_argument("-total", help="the total length of the genome")
    pars.add_argument("-step", help="the step size of the sliding window method")
    pars.add_argument("-ws", help="the window size of the sliding window method")
    pars.add_argument("-scale", help="the value for normalization of the SNP density")
    pars.add_argument("-scan_loop",help="the number of scanning times the program will perform for the SNPs.Default value is set to 100.")
    pars.add_argument("-max_dist",help="maximum inter-SNP distance allowed within a cluster.Default value is set to 4000 bp.")
    pars.add_argument("-min_num", help=" minimum number of SNPs per cluster.The default value is set to 2.")
    args = pars.parse_args()
    if(args.AnalType != "density" and args.AnalType != "cluster" and args.AnalType != "Pval"):
        raise Exception("Please provide correct analysis type:density or cluster or Pval")
    if(args.AnalType == "density"):
        if args.input:
            input_dx = sys.argv.index("-input")
            input_file = sys.argv[input_dx + 1]
        else:
            raise Exception("Please input the vcf file!")
        if args.output:
            output_dx = sys.argv.index("-output")
            output_file =sys.argv[output_dx + 1]
        else:
            raise Exception("Please provide the output file!")
        if args.total:
            totLg_dx = sys.argv.index("-total")
            tot_lg = int(sys.argv[totLg_dx + 1])
        else:
            tot_lg = 0
            print("No total genome size is provided; the maximum location of SNPs will be used!")

        if args.step:
            stp_dx = sys.argv.index("-step")
            stp = int(sys.argv[stp_dx + 1])
        else:
            stp = 300
            print("No step size is provided; default value of 300 will be used!")

        if args.ws:
            wind_dx = sys.argv.index("-ws")
            wind = int(sys.argv[wind_dx + 1])
        else:
            wind = 2000
            print("No window size is provided; default value of 2000 will be used!")

        if args.scale:
            scale_dx = sys.argv.index("-scale")
            scale = int(sys.argv[scale_dx + 1])
        else:
            scale = 1000
            print("No scale factor is provided; default value of 1000 will be used!")
        density(vcf_file=input_file,out_file=output_file,genome_length=tot_lg,step_size=stp,window_size=wind,scale=scale)
    if (args.AnalType == "cluster"):
        if args.input:
            input_dx = sys.argv.index("-input")
            input_file = sys.argv[input_dx + 1]
        #		snp_list=snp_file.read().rstrip()
        else:
            raise Exception("Please input the vcf file!")

        if args.anno:
            anno_dx = sys.argv.index("-anno")
            anno_file =sys.argv[anno_dx + 1]
        #		annotation=anno_file.read().rstrip()
        else:
            raise Exception("Please input the gene annotation file!")

        if args.operon:
            operon_dx = sys.argv.index("-operon")
            operon_file = sys.argv[operon_dx + 1]
        #		operon=operon_file.read().rstrip()
        else:
            raise Exception("Please input the operon file!")

        if args.recomb_lg:
            recombLg_dx = sys.argv.index("-recomb_lg")
            recomb_lg = int(sys.argv[recombLg_dx + 1])
        else:
            raise Exception("Please input the estimated recombination tract length.")

        if args.output:
            outputClust_dx = sys.argv.index("-output")
            output_file = sys.argv[outputClust_dx + 1]
        else:
            raise Exception("Please designate output file for cluster list")

        if args.scan_loop:
            scanLoop_dx = sys.argv.index("-scan_loop")
            scan_loop = int(sys.argv[scanLoop_dx + 1])
        else:
            scan_loop = 100
            print("No scan_loop number is provided, default value of 100 will be used.")

        if args.max_dist:
            maxDist_dx = sys.argv.index("-max_dist")
            max_dist_cluster = int(sys.argv[maxDist_dx + 1])
        else:
            max_dist_cluster = 5000
            print("No max_dist_cluster is provided, default value of 5000 bp will be used.")

        if args.min_num:
            minNumSnp_dx = sys.argv.index("-min_num")
            min_num_snp = int(sys.argv[minNumSnp_dx + 1])
        else:
            min_num_snp = 2
            print("No min_num_snp is provided, default value of 2 will be used.")
        targetAnno_dic,targetClustRev_dic,targetClust_dic,dist_dic,distRev_dic,targetClustRev_dic={},{},{},{},{},{}
        targetAnno_dic,targetOperon_dic=SNP_anno(vcf_file=input_file,anno_file=anno_file)
        ### store the annotations for all SNPs and annotate the target SNP list, caution: time-consuming
        snp_AllPos,targetClustRev_dic=Operon_manage(operon_file=operon_file,targetOperon_dic=targetOperon_dic,targetAnno_dic=targetAnno_dic)
        ### handle the operon file including:profile the gene operons,define the range of each gene operon
        #re-define the gene operons for target SNPs due to the fact that some operons may include several genes
        # profile the snp list for each operon
        # initialize the snp clusters if containing non-synonymous snp
        targetClust_dic=clust_initalize(targetClustRev_dic=targetClustRev_dic)### initialize the targetClust_dic
        dist_dic,distRev_dic=SNP_dist(snp_AllPos=snp_AllPos)
        ### loop and update targetClustRev and targetClust
        Clust_lst_single = list(targetClustRev_dic.keys())
        Clust_lst = Clust_lst_single * scan_loop
        ClustLg = len(Clust_lst)
        clust_order,NewClust3Lg,targetClustRev_dic,NewClust2Lg=cluster_merge(Clust_Length=ClustLg,Clust_list=Clust_lst,targetClustRev_dic=targetClustRev_dic,dist_dic=dist_dic,distRev_dic=distRev_dic,snp_AllPos=snp_AllPos,targetClust_dic=targetClust_dic,recomb_lg=recomb_lg,scan_loop=scan_loop,targetOperon_dic=targetOperon_dic,max_dist_cluster=max_dist_cluster,min_num_snp=min_num_snp,targetAnno_dic=targetAnno_dic)
        cluster_out(clust_order, NewClust3Lg, targetClustRev_dic, NewClust2Lg, recomb_lg, output_file)
    if (args.AnalType == "Pval"):
        if args.input:
            input_dx = sys.argv.index("-input")
            input_file = sys.argv[input_dx + 1]
        else:
            raise Exception("Please input the snp cluster file!")
        if args.output:
            output_dx = sys.argv.index("-output")
            output_file =sys.argv[output_dx + 1]
        else:
            raise Exception("Please provide the output file!")
        if args.rate:
            mutRat_dx = sys.argv.index("-rate")
            mutRate = float(sys.argv[mutRat_dx + 1])
        else:
            raise Exception("Please provide the mutation rate!")
        cluster_Pval(input_file, output_file, mutRate )



