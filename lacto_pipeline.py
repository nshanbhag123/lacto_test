#%%
import glob
import os
import random
import subprocess
import pandas as pd
import numpy as np
import argparse
#%%
def count_reads(fastq):
    f = open(fastq, "r")
    lines = f.readlines()
    f.close()
    return (len(lines))

def init_folders(output_folder, sample_no):
    folders = [f"{output_folder}/{sample_no}", f"{output_folder}/{sample_no}/ss", f"{output_folder}/{sample_no}/trim_ss", f"{output_folder}/{sample_no}/assemble", f"{output_folder}/{sample_no}/ani"]
    # Create paths to each directory
    for fold in folders:
        if not os.path.exists(fold):
            os.makedirs(fold)

def subsample_trim_assemble(data_folder, output_folder, subsamples, no_reads_in_subsample):
    ordered_fastqs = sorted([fastq for fastq in glob.glob(f"{data_folder}/*.fastq")]) #sorted fastqs
    for i in range(0, len(ordered_fastqs), 2): #for each paired fastq file (you need to subsample 50 times)
        fastq1_path = ordered_fastqs[i] #fastq1
        fastq2_path = ordered_fastqs[i+1] #fastq2
        file_name = os.path.basename(fastq1_path) #file name
        sample_no = file_name[0:file_name.find("_")] #sample no
                
        random.seed(1) #setting seed
        no_reads_in_sample = count_reads(fastq1_path) #counting number of reads
        randomList=random.sample(range(0,no_reads_in_sample),subsamples) #50 random numbers from 0 to no_reads range - your seeds for seqtk

        init_folders(output_folder, sample_no) #initialize folders

        for index,seed in enumerate(randomList): #subsampling reads 50 times
            fastq1 = f"{sample_no}_1"
            fastq2 = f"{sample_no}_2"
            os.chdir(f"{output_folder}/{sample_no}/ss")

            seqtk_fastq1 = f"seqtk sample -s{seed} {fastq1_path} {no_reads_in_subsample} > {fastq1}_ss{index}.fastq"
            seqtk_fastq2 = f"seqtk sample -s{seed} {fastq2_path} {no_reads_in_subsample} > {fastq2}_ss{index}.fastq"
            subprocess.run(f'conda run -n lacto {seqtk_fastq1}', shell = True)
            subprocess.run(f'conda run -n lacto {seqtk_fastq2}', shell = True)

            bbduk = f"bbduk.sh -Xmx1G overwrite=t in1={fastq1}_ss{index}.fastq in2={fastq2}_ss{index}.fastq out1={output_folder}/{sample_no}/trim_ss/{fastq1}_ss{index}_trim.fastq out2={output_folder}/{sample_no}/trim_ss/{fastq2}_ss{index}_trim.fastq qtrim=rl ftl=15 ftr=135 maq=20 maxns=0 stats={output_folder}/{sample_no}/trim_ss/{sample_no}_ss{index}_trim.stats statscolumns=5 trimq=20"
            subprocess.run(f'conda run -n lacto {bbduk}', shell = True)

            os.chdir(f"{output_folder}/{sample_no}/trim_ss")
            spades_path = "/home/nshanbhag/software/SPAdes-3.15.5-Linux/bin"
            spades = f"python3 {spades_path}/spades.py --only-assembler -t 12 -1 {fastq1}_ss{index}_trim.fastq -2 {fastq2}_ss{index}_trim.fastq -o {output_folder}/{sample_no}/assemble/{sample_no}_ss{index}_assembly"
            os.system(spades)

            os.chdir(f"{output_folder}/{sample_no}/assemble")
            copy = f"cp {sample_no}_ss{index}_assembly/contigs.fasta {sample_no}_ss{index}_contigs.fasta"
            os.system(copy)

        remove = f"rm -r {output_folder}/{sample_no}/assemble/*assembly" #removing the extra junk from spades
        os.system(remove)

        calc_anis(output_folder, sample_no)

def calc_anis(output_folder, sample_no):
    assemblies = [assembly for assembly in glob.glob(f"{output_folder}/{sample_no}/assemble/*_contigs.fasta")]
    pairwise_sim_anis = open(f"{output_folder}/{sample_no}/ani/pairwise_sim_anis.tsv", "w")
    i = 1
    while len(assemblies) != 0:
        random.seed(1)
        pair = random.sample(assemblies,2)

        fastani = f"fastANI -q {pair[0]} -r {pair[1]} -o {output_folder}/{sample_no}/ani/ani{i}.out"
        subprocess.run(f'conda run -n ani {fastani}', shell = True)

        f = open(f"{output_folder}/{sample_no}/ani/ani{i}.out")
        ani = f.readline()
        pairwise_sim_anis.write(f"{ani}")
        f.close()

        assemblies.remove(pair[0])
        assemblies.remove(pair[1])
        i+=1 

    pairwise_sim_anis.close()

def collate_anis(output_folder):
    dicty = {}
    for sample_fold in glob.glob(f"{output_folder}/*"):
        print(sample_fold)
        sample = os.path.basename(sample_fold)
        ani = pd.read_csv(f"{sample_fold}/ani/pairwise_sim_anis.tsv", delimiter= "\t", header = None).iloc[:,2].tolist()
        dicty[sample] = ani
    df = pd.DataFrame(dicty)
    df.to_csv(f"{output_folder}/coallate_ani.csv")
    return(df)

#ARG PARSE
parser = argparse.ArgumentParser(description='Simulate Distribution of ANI values for a Given Species')
parser.add_argument('--subsamples','-sub', type = int, default = 50, help='The number of subsamples per paired-end fastq')
parser.add_argument('--reads','-r', type=int, default = 150000, help='Number of reads in each subsample - can adjust based on coverage')
parser.add_argument('--output', '-o', required = True, type = str, help = 'Destination path of output file')
parser.add_argument('--input', '-i', required = True, help = 'Path to paired end fastqs folder')
args = parser.parse_args()

# START RUNNING FUNCTIONS
if glob.glob(f"{args.input}"):
    print(args.input, args.output, args.subsamples, args.reads)
    # subsample_trim_assemble(args.input, args.output, args.subsamples, args.reads)
    # df = collate_anis(args.output)
    # os.system(f'Rscript /home/nshanbhag/lacto_compare/stats.r')
else:
    print("Data was not found in the given input folder. Please try again")
# %%
