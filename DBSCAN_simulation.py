import multiprocessing as mp
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
import vcf
import argparse
import sys
import os


def DBSCAN_test(args):
    pars = argparse.ArgumentParser()
    pars.add_argument("-vcf", required=True, help="The vcf file.")
    pars.add_argument("-eps_start",help="The lower-bound of the range of eps. Default value is the minimum of the average distance between SNPS.")
    pars.add_argument("-eps_end", ,help="The upper-bound of the range of eps. Default value is the maximum of the average distance between SNPS.")
    pars.add_argument("-min_start", help="The lower-bound of the range of min_num. Default value is set to 5.")
    pars.add_argument("-min_end", help=" The upper-bound of the range of min_num. Default value is set to 50.")
    args = pars.parse_args()
    
    if args.vcf:
        input_dx = sys.argv.index("-vcf")
        vcf_file = sys.argv[input_dx + 1]
    else:
        raise Exception("Please input the vcf file!")
        
    if args.eps_start:
        eps_start_dx = sys.argv.index("-eps_start")
        eps_start = int(sys.argv[eps_start_dx + 1])
    else:
        print("No eps_start is provided; the minimum of the average distance between SNPS will be used!")
        eps_start=0
        
    if args.eps_end:
        eps_end_dx = sys.argv.index("-eps_end")
        eps_end = int(sys.argv[eps_end_dx + 1])
    else:
        print("No eps_end is provided; the maximum of the average distance between SNPS will be used!")
        eps_end=0
        
    if args.min_start:
        sample_start_dx = sys.argv.index("-min_start")
        sample_start = int(sys.argv[sample_start_dx + 1])
    else:
        print("No min_sample start is provided; default value of 5 will be used!")
        sample_start=5
        
    if args.min_end:
        sample_end_dx = sys.argv.index("-min_end")
        sample_end = int(sys.argv[sample_end_dx + 1])
    else:
        print("No min_sample end is provided; default value of 50 will be used!")
        sample_end=50
        
    loc=[]
    try:
        vcf_reader = vcf.Reader(filename=vcf_file)
        for record in vcf_reader:
            loc.append(int(record.POS))
        loc.sort()
        
    except:
        for pos in open(vcf_file):
            POS = int(pos)
            loc.append(POS)
        loc.sort()
    
    
    #Calculate the distance between the two SNPS
    one_list=[]
    two_list=[]
    
    for i in range(1,len(loc),2):
        one_list.append(loc[i])
        
    for i in range(0,len(loc),2):
        two_list.append(loc[i])
        
    if(len(one_list)>len(two_list)):
        one_list=one_list[0:-1]
    else:
        two_list=two_list[0:-1]
        
    snp_dist=[]
    for i in range(0,len(one_list)):
        dist=one_list[i]-two_list[i]
        snp_dist.append(dist)
        
        
    #Calculate the mean distance between SNPS, the maximum distance and the minimum distance
    size=2000
    num_list=[]
    
    for i in range(0,len(snp_dist),size):
        num_list.append(i)
        
    aver_list=[]
    for i in range(0,len(num_list)-1):
        dist_sum=0
        for index in range(num_list[i],num_list[i+1]):
            dist_sum+=snp_dist[index]
        aver=dist_sum*1.0/size
        aver_list.append(aver)
        
    if(eps_start==0):
        eps_start=int(min(aver_list))
        print("The minimum of the average distance of SNPS is:",min(aver_list))
        
    if(eps_end==0):
        eps_end=int(max(aver_list))
        print("The maximum of the average distance of SNPs is :",max(aver_list))
        
    total_sum=0
    for i in aver_list:
        total_sum+=i
    print("The average distance of SNPs",total_sum/len(aver_list))
    size=int((eps_end-eps_start)/50)
    loc=np.array(loc)
    x=loc.reshape(-1,1)
    
    
    def cal_std(eps,min_sample,loc):
        y_pred = DBSCAN(eps= eps, min_samples = min_sample).fit_predict(x)
        cluster_num=(max(y_pred))
        loc=list(loc)
        y_pred=y_pred.tolist()
        
        std_sum=0
        data_dic={}
        data_dic["loc"]=loc
        data_dic["index"]=y_pred
        data_df=pd.DataFrame(data_dic)
        data_df=data_df.groupby("index")
        
        for i in range(0,cluster_num+1):
            std_sum+=(data_df.get_group(i)["loc"].std())
            
        if(cluster_num==0):
            aver_std=0
        else:
            aver_std=(std_sum/cluster_num)
            
        with open(r"log.txt","a+") as f:
            print("%d\t%d\t%d\t%.3f" %(eps,min_sample,cluster_num,aver_std),file=f)
 

    pool = mp.Pool(mp.cpu_count())
    print(pool)
    for eps in range(eps_start,eps_end,size):
        pool.starmap(cal_std,[(eps,min_sample,loc) for min_sample in range(sample_start,sample_end,5)])
    pool.close()
    os.system('Rscript ./library/DBSCAN_simulation_Rscript log.txt')
    
    
