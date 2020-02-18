import os
from scipy import stats


###Defines a function that computes the p-value of the cluster and passes the parameters of the main program
def Pval(args):
    args_dic = vars(args)
    if (args_dic["cluster"] == None):
        raise Exception("Please input the snp cluster file!")
    else:
        cluster_file = args_dic["cluster"]
    if (args_dic["out"] == None):
        raise Exception("Please provide the output file!")
    else:
        out_file = args_dic["out"]
    if (args_dic["rate"] == None):
        raise Exception("Please provide the mutation rate!")
    else:
        mutRate = float(args_dic["rate"])
    print("The snp cluster file is:",os.path.abspath(cluster_file))
    print("The output file providing the p-value of clustering for each SNP cluster is :",os.path.abspath(out_file))
    print("the pre-estimated average mutation rate is :",mutRate)
    mutRate = float(mutRate)
    clust_file=open(cluster_file,'r')
    next(clust_file)
    #mutRate=0.0001
    big_sz=100000000
    buff0,buff1,buff='','',''
    clustID_tab,clustSiz_tab,clustLg_tab=0,3,4
    outStr=''
    while True:
        buff1=clust_file.read(big_sz).rstrip()
        if buff1 == '':
            buff=buff0+'\n'
        else:
            buff=buff0+buff1

        st=0
        while True:
            ed=buff.find('\n',st+1)
            if ed==-1:
                buff0=buff[st:]
                break
            the_buff=buff[st:ed].lstrip().rstrip('\n')
            the_buff_tmp=the_buff.split('\t')
            clust_size=int(the_buff_tmp[clustSiz_tab])
            clust_lg=int(the_buff_tmp[clustLg_tab])
            clust_ID=the_buff_tmp[clustID_tab]

            pVal=stats.gamma.cdf(clust_lg,clust_size,loc=0,scale=1.0/mutRate)
            outStr+='%s\t%d\t%d\t%.4f\n' % (clust_ID,clust_size,clust_lg,pVal)
            st=ed
        if buff1=='':
            break
    clust_file.close()
    output_file=open(out_file,'w')
    output_file.write("ID\tSNP_num\tlength\tp-value\n")
    output_file.write(outStr)
    output_file.close()