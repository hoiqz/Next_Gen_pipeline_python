import re
import gzip
import glob
import sys

folder=sys.argv[1]
gene=sys.argv[2]

print (folder)
print (gene)
#read the genelist 
genedict=dict()
readgene=open(gene,"r")
for gene in readgene:
    if re.search('^#',gene):
        next
    else:
        gene=gene.replace("\n","")
        gene_arr=gene.split("\t")
        gene_name=gene_arr[0]+" "+gene_arr[1]
        coor=gene_arr[1]
        if gene_name in genedict:
            next
        else:
            genedict[gene_name] =0
# PRINT THE HEADER
header=["file_name"]
for key,value in genedict.iteritems():
#    print (key,value)
    header.append(key)
header="\t".join(header)
print header

def checkposition(check_chr,check_pos):
    for ori_keys,value in genedict.iteritems():
#        print (keys)
        keys=str(ori_keys)
#        print ("check key")
#        print (keys)
        key_arr=keys.split()
#        print (key_arr)
        chr_pos=key_arr[1]
#        print(chr_pos)
        gene_chr,range=chr_pos.split(":")
        gene_chr="chr"+gene_chr
        start,end=range.split("-")
        start=int(start)
#        print(start,end)
        end=int(end)
        
        #print("checking this %s %d against chromosome %s %d-%d"%(check_chr,check_pos,gene_chr,start,end))
        if gene_chr == check_chr:
            #print("%s match %s" %(gene_chr,check_chr))
            if check_pos >= start and check_pos <= end:
                genedict[ori_keys] += 1
            else:
                next




# now i take all the .gz files in the folder and process them
# assumption is that the ion torrent files are in gz form

outfile=open(folder+"/result.tsv","w")
outfile.write(header+"\n")
get_all=glob.glob(folder+"/*.gz")
for file in get_all:
    print file
    with gzip.open(file,'r') as infile:
        #set the dict values to zero before starting
        genedict=dict.fromkeys(genedict,0)
        for line in infile:
            line=line.replace("\n","")
            if re.search('^#',line):
                next
            else:
                split_line=line.split("\t")
                chr_now=split_line[0]
                position_now=int(split_line[1])
                #print(chr_now,position_now)
                checkposition(chr_now,position_now)
        #once finish reading file print the value into the  file
        toprint=[file]
        for key,value in genedict.iteritems():
            print (key,value)
            toprint.append(str(value))
        outfile.write("\t".join(toprint)+"\n")




          #  print (line)


print ("end")
outfile.close()

    
