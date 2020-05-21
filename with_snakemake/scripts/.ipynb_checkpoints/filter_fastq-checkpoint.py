from Bio import SeqIO
import gzip
import subprocess



if snakemake.input[0].endswith(".fq"):
    filteredPath = snakemake.input[0][:-3] + "_filtered.fq"
if snakemake.input[0].endswith(".fq.gz"):
    filteredPath = snakemake.input[0][:-6] + "_filtered.fq"
if snakemake.input[0].endswith(".fastq"):
    filteredPath = snakemake.input[0][:-6] + "_filtered.fastq"
if snakemake.input[0].endswith(".fastq.gz"):
    filteredPath = snakemake.input[0][:-9] + "_filtered.fastq"
if snakemake.input[0].endswith(".fq.gz") or snakemake.input[0].endswith(".fastq.gz"):
    with gzip.open(r"{}".format(snakemake.input[0]), "rt") as fo, open(snakemake.output[0], "w") as filter_fo:
        record_all = []
        seq = []
        for record in SeqIO.parse(fo, "fastq"):
            if record.seq.count("N") < 10: # reads with less than 10 "N" will be passed on
                count = 0
                for qual in record.letter_annotations["phred_quality"]:
                    if qual < 20:
                        count += 1
                poor_perc = count / len(record.letter_annotations["phred_quality"]) # this will calculate the percentage of bases with a score less than 20 (which means a 1% chance of being wrong)
                if poor_perc < 0.6:       
                    SeqIO.write(record, filter_fo , "fastq")
if snakemake.input[0].endswith(".fq") or snakemake.input[0].endswith(".fastq"):
    with open(r"{}".format(snakemake.input[0]), "rt") as fo, open(snakemake.output[0], "w") as filter_fo:
        record_all = []
        seq = []
        for record in SeqIO.parse(fo, "fastq"):
            if record.seq.count("N") < 10: # reads with less than 10 "N" will be passed on
                count = 0
                for qual in record.letter_annotations["phred_quality"]:
                    if qual < 20:
                        count += 1
                poor_perc = count / len(record.letter_annotations["phred_quality"]) # this will calculate the percentage of bases with a score less than 20 (which means a 1% chance of being wrong)
                if poor_perc < 0.6:       
                    SeqIO.write(record, filter_fo, "fastq")
                    

