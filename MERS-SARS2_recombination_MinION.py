#!/bin/python3
import pandas as pd
import os
import subprocess

wd = "/home/denison-thelio/Current_projects/MinION/MERS-SARS2_recombination_2020/"
od = wd
sample_file = wd + "samples.txt"
samples = [line.rstrip('\n') for line in open(wd + "samples.txt")]

# Set new folder for chimeric junction files
if not os.path.exists(od + '/Chimeric_Junctions'):
    os.makedirs(od + '/Chimeric_Junctions')

# Filter chimeric junctions by inner join
for name in samples:
    mers_bed = pd.read_csv(wd + name + "/" + name + "_MERS.bed",
                           sep="\t",
                           names=['Genome',
                           'Read_Start',
                           'Read_End',
                           'Name',
                           'Qual',
                           'strand',
                           'Start',
                           'Stop',
                           'RGB',
                           'Blocks',
                           'Block_Length',
                           'Block_Start'])
    sars2_bed = pd.read_csv(wd + name + "/" + name + "_SARS2.bed",
                           sep="\t",
                           names=['Genome',
                                  'Read_Start',
                                  'Read_End',
                                  'Name',
                                  'Qual',
                                  'strand',
                                  'Start',
                                  'Stop',
                                  'RGB',
                                  'Blocks',
                                  'Block_Length',
                                  'Block_Start'])
    mers_bed = mers_bed.drop(['Qual', 'Start', 'Stop', 'RGB'], axis=1)
    sars2_bed = sars2_bed.drop(['Qual', 'Start', 'Stop', 'RGB'], axis=1)
    recombined_bed = mers_bed.merge(sars2_bed, how="inner", on="Name", suffixes=("_m", "_s"))
    # command = 'samtools fasta -n -@ 32 ' + name + "/" + name + '_SARS2.sort.bam > ' + name + "/" + name + '.fasta'
    # process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    # output, error = process.communicate()
    recombined_reads = recombined_bed['Name']
    recombined_bed.to_csv(od + "Chimeric_Junctions/" + name + "_recombined_reads.txt", sep="\t", index=False)
    recombined_reads.to_csv(od + name + "_read_names.txt", header=False, index=False, sep="\t")
    command2 = "sed 's/..$//' " + wd + name + "_read_names.txt > " + od + name + "_read_names.txt"
    command3 = "seqtk subseq " + name + "/" + name + "_SARS2.fasta " + name + "_read_names.txt > " + od + "Chimeric_Junctions/" + name + "_recombined.fasta"
    process = subprocess.Popen(command2.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    process = subprocess.Popen(command3.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()