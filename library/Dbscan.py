import vcf
import numpy
from sklearn.cluster import DBSCAN


###DBSCAN algorithm in machine learning was used for SNP clustering
def DBSCA(args):
    args_dic = vars(args)
    
    if (args_dic["vcf"] == None):
        raise Exception("Please provide the vcf file!")
    else:
        vcf_file = args_dic["vcf"]
        
    if (args_dic["eps"] == None):
        raise Exception("Please provide the value of eps!")
    else:
        eps=int(args_dic["eps"])
        
    if(args_dic["min_sample"]==None):
        raise Exception("Please provide the value of min_sample!")
    else:
        min_sample=int(args_dic["min_sample"])
        
    if (args_dic["out"] == None):
        raise Exception("Please provide the output file!")
    else:
        out_file = args_dic["out"]
        
        
    vcf_reader = vcf.Reader(filename=vcf_file)
    loc = []
    
    for record in vcf_reader:
        loc.append(record.POS)
        
    loc.sort()
    loc = numpy.array(loc)
    x = loc.reshape(-1, 1)
    y_pred = DBSCAN(eps=eps, min_samples = min_sample).fit_predict(x)
    cluster_num = max(y_pred)
    
    y_pred = y_pred.tolist()
    loc = list(loc)
    
    with open(out_file, 'w') as f:
        print("ClusterID\tstart\tend\tsize\tlength", file=f)
        
        for i in range(0, cluster_num + 1):
            pos_list = []
            for index in range(0, len(loc)):
                if y_pred[index] == i:
                    pos_list.append(loc[index])
                    
            cluster_ID = i + 1
            cluster_start = min(pos_list)
            cluster_end = max(pos_list)
            size = len(pos_list)
            length = cluster_end - cluster_start
            
            print("%d\t%d\t%d\t%d\t%d" % (cluster_ID, cluster_start, cluster_end, size, length), file=f)
