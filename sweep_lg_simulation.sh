#!/bin/bash

dir=$(pwd)
scan_loop=100
max_dist=4000
min_num=2

echo
while [ -n "$1" ]
do
	case "$1" in
		-vcf)	param="$2"
			vcf_file=$2
			echo "The vcf file is $2"
			shift;;
		-anno)  param="$2"
			anno_file=$2
			echo "The annotation file is $2"
			shift ;;
		-operon) param="$2"
			operon_file=$2
			echo "The operon file is $2"
			shift ;;
		-start) param="$2"
			start_num=$2
			echo "The start recombination length is $2"
			shift ;;
		-end) param="$2"
			end_num=$2
			echo "The end recombination length is $2"
			shift ;;
		-step) param="$2"
			step_size=$2
			echo "The step size of recombination length is $2"
			shift ;;
		-delta) param="$2"
			interv=$2
			echo "The gam model interv value is $2"
			shift ;;
		-scan_loop) param="$2"
			scan_loop=$2
			shift ;;
		-max_dist) param="$2"
			max_dist=$2
			shift ;;
		-min_num) param="$2"
			min_num=$2
			shift ;;
	esac
	shift
done

echo "the number of scanning times the program will perform for the SNPs is $scan_loop"
echo "maximum inter-SNP distance allowed within a cluster is $max_dist"
echo "minimum number of SNPs per cluster is $min_num"

mkdir $dir/snp_cluster_result
#
count=1
for param in "$@"
do
	echo "Parameter #$count: $param"
	count=$[ $count + 1 ] 
done

length=$(($end_num-$start_num))
echo $length

number=$(($length/$step_size))
echo "the toatl number of recombination length is $number "

starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo "Start running SNP clustering in parallel"
echo $(seq $start_num $step_size $end_num)

func()
{
	python SweepCluster.py Cluster -vcf $2 -anno $3 -operon $4 -sweep_lg $1 -cluster $5/snp_cluster_result/clust_out_$1 -snp $5/snp_cluster_result/snp_out_$1  -scan_loop $6 -max_dist $7 -min_num $8
}
export -f func

parallel -j 8 func ::: $(seq $start_num $step_size $end_num) ::: $vcf_file ::: $anno_file ::: $operon_file ::: $dir ::: $scan_loop ::: $max_dist ::: $min_num
echo "The results of the clustering are saved in the snp_cluster_result"

grep "recomb_lg" log.txt | awk -F "=" '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > recomb.txt
rm log.txt
echo "The result of the relationship between recombination length and the cluster number are saved in recomb_result.txt"

echo "Start to run R program and use generalized additive model and linear model to fit the data"
Rscript ./library/sweep_lg_Rscript $interv
rm recomb.txt
echo "The description of the R program is stored in Rcode_log.txt"
echo "The results of the R program run are saved in Rplots.pdf"

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s)
end_seconds=$(date --date="$endtime" +%s)
echo "The total time used is:： "$((end_seconds-start_seconds))"s"
