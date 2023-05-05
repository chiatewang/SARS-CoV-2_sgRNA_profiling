# generate sample list contain fastq path and variant name
import os
import glob
import pandas
cwd = os.getcwd()
fastq_path = glob.glob(cwd+"/fastq/*")
with open("sample.txt", "w") as f:
    f.write("id\tfq1\tfq2\tvariant_name\n")

for i in fastq_path:
    prefix = i.split('/')
    variant_name = prefix[-1].split('_')
    path = glob.glob(i+"/*_1*")
    for k in path:        
        a = k.split('/')
        sample_name = a[-1].split('_')
        sample_name = sample_name[0]
        with open("sample.txt", "a",encoding='utf-8') as f:
            f.write(str(sample_name))
            f.write('\t')
            fq1 = i+"/"+sample_name+"_1.fastq"
            f.write(str(fq1))
            f.write('\t')
            fq2 = i+"/"+sample_name+"_2.fastq"
            f.write(str(fq2))
            f.write('\t')
            variant=os.path.basename(os.path.abspath(os.path.join(k, os.path.pardir)))
            f.write(str(variant)+'\n')
            f.close