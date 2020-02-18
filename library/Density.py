import vcf
import os
import math

###Define a function for calculating the density of SNPS
def calu_density(vcf_file,out_file,genome_length=0,step_size=300,window_size=2000,scale=1000):
    print("the vcf file is :",os.path.abspath(vcf_file))
    print("the snp density file is :",os.path.abspath(out_file))
    loc_dic = {}
    vcf_reader=vcf.Reader(filename =vcf_file)
    loc=[]
    for record in vcf_reader:
        loc.append(record.POS)
    loc.sort()
    for i in loc:
        loc_dic[i]=1
    allKys = loc_dic.keys()
    allKys=list(allKys)
    allKys.sort()
    if(genome_length<=allKys[-1]):
        tot_lg = allKys[-1]
    else:
        tot_lg=genome_length
    print("the total genome size is : ",tot_lg)
    stp=step_size
    print("the step size is : ",stp)
    wind =window_size
    print("the window size is : ",wind)
    print("the value for normalization of the SNP density is : " ,scale)
    eff_lg=math.trunc((tot_lg-wind)/stp)*stp
    print("the effective length is %d: " % eff_lg)
    print("the total length is %d: " % tot_lg)
    wind_dic=[]

    for ii in range(0,eff_lg,stp):
        if ii  not in loc_dic:
            loc_dic[ii]=0
        if (ii+wind)  not in loc_dic:
            loc_dic[ii+wind]=0
        wind_dic += [[ii,wind+ii]]

    if eff_lg  not in loc_dic:
        loc_dic[eff_lg]=0

    if (eff_lg+wind) not in loc_dic:
        loc_dic[eff_lg+wind]=0
    wind_dic += [[eff_lg,eff_lg+wind]]

    allKys=loc_dic.keys()
    allKys=list(allKys)
    allKys.sort()

    out_str=''
    for ech in wind_dic:
        left=allKys.index(ech[0])
        right=allKys.index(ech[1])
        rat=0
        for dd in range(left,right):
            rat+=loc_dic[allKys[dd]]
        out_str+='%d\t%d\t%.4f\n' % (ech[1],rat,rat*1.0*scale/wind)
    output_file=open(out_file,'w')
    output_file.write("Win_start\tSNP_num\tDensity\n")
    output_file.write(out_str)
    output_file.close()


###Defines a function that calls the function:calu_density, and passes the arguments to the main program
def density(args):
    args_dic=vars(args)
    if(args_dic["vcf"]==None):
        raise Exception("Please provide the vcf file!")
    else:
        vcf_file = args_dic["vcf"]
    if(args_dic["out"]==None):
        raise Exception("Please provide the output file!")
    else:
        out_file = args_dic["out"]
    if(args_dic["total"]==None):
        tot_lg = 0
        print("No total genome size is provided; the maximum location of SNPs will be used!")
    else:
        tot_lg=int(args_dic["total"])
    if(args_dic["step"]==None):
        stp = 300
        print("No step size is provided; default value of 300 will be used!")
    else:
        stp=int(args_dic["step"])
    if(args_dic["ws"]==None):
        wind = 2000
        print("No window size is provided; default value of 2000 will be used!")
    else:
        wind=int(args_dic["ws"])
    if(args_dic["scale"]==None):
        scale = 1000
        print("No scale factor is provided; default value of 1000 will be used!")
    else:
        scale=int(args_dic["scale"])
    calu_density(vcf_file=vcf_file, out_file=out_file, genome_length=tot_lg, step_size=stp, window_size=wind,scale=scale)
