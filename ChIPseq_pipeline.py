#!/usr/bin/env python3

from Bio import SeqIO
import os, sys, re, glob, shutil
import subprocess
import gzip
import argparse
import pandas as pd
from datetime import datetime

# other things to add:
## read trimming
## ngs plot
## deeptools QC 
###correlation
###coverage
###pca
## peak calling - homer and macs2?
## make the log separate for each sample??

# dependencies - picard, fastqc, bwa, bowtie2, samtools, deeptools

def filterFastq(fastqPath):
    if fastqPath.endswith(".fq"):
        filteredPath = fastqPath[:-3] + "_filtered.fq"
    if fastqPath.endswith(".fastq"):
        filteredPath = fastqPath[:-6] + "_filtered.fastq"
    with open(r"{}".format(fastqPath), "r") as fo, open(r"{}".format(filteredPath), "w") as filter_fo:
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
    return(filteredPath)
     
        
def downloadFasta(species):
    # download 'soft-masked' primary assembly genomes per this thread - https://bioinformatics.stackexchange.com/questions/540/what-ensembl-genome-version-should-i-use-for-alignments-e-g-toplevel-fa-vs-p
    possible_genomes = ["hg38", "hg19", "mm10"]
    if species not in possible_genomes:
        print("FAILED:'genome' must be 'hg38', 'hg19', or 'mm10'")
        exit(3)
    if species == "hg38":
        fileGrab = "wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
    elif species == "hg19":
        fileGrab = "wget ftp://ftp.ensembl.org/pub/grch37/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz"
    elif species == "mm10":
        fileGrab = "wget ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
    fastaGrab_run = subprocess.run(fileGrab, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    if fastaGrab_run.returncode != 0:
        print("FAILED! Problem retrieving fasta from Ensembl FTP.")
        exit(2)
    else:
        contents=os.listdir()
        for files in contents:
            if files.endswith('dna_sm.primary_assembly.fa.gz'):
                fastaGz = files
                with gzip.open(fastaGz,'rb') as f_in:
                    with open(fastaGz[:-3],'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
        newName = fastaGz[:-3]
        os.remove(fastaGz)
        print('\n\n{} was successfully downloaded.'.format(newName) )
    return(newName)
       

def main():
    parser = argparse.ArgumentParser(description="""
        Pipeline for ChIPseq 
        """)
    parser.add_argument('-f','--fasta', help='The aligner to use, bwa-mem will be used by default')
    parser.add_argument('-a','--aligner', choices = ['bwa-mem', 'bowtie2'], help='The name of the fasta file for the reference genome (eg hg19.fa). This can be a full path if the file is not in the current working directory.')
    parser.add_argument('-q','--pathToFastqs', help='The path to the folder that contains the fastq files.')
    parser.add_argument('-g', '--genome' , choices = ['hg38', 'hg19', 'mm10'], help="This is required if no fasta reference file is provided in the '-f' argument. This will determine which reference genome fasta is downloaded. This is also required if peak calling is desired")
    parser.add_argument('-i', '--index', choices = ['True', 'False'], help='If a reference fasta is provided, please indicate whether you need it to be indexed with BWA')
    parser.add_argument('-d', '--remove_duplicates', choices = ['True', 'False'], help='This will determine if duplicates should be removed, by defualt they will be flagged in the bam file, but not removed.')
    parser.add_argument('-s', '--sampleSheet', help='If peaks should be called, the path to the peak sample sheet should be provided here. If no path is included, peaks will not be called.')
    parser.add_argument('-p', '--peakCaller', choices = ['macs2', 'homer'], help='The peak caller to be used.')
    parser.add_argument('-x', '--extraPeakArgs', help="A string to be added onto either macs2 or homer command in addition to the standard arguments. Make sure to put in quotes. E.g. -x '--bw 300 -q 0.01'")
    parser.add_argument('-c', '--pcaCorr', choices = ['True', 'False'], help="Perform PCA oand correlation analysis with deepTools")
    

    args = parser.parse_args()
    
    pathtoFastqs = args.pathToFastqs
    
    if args.genome:
        genome = args.genome
        
    Log_path = pathtoFastqs + "/StdError_Stdout_log.txt"
    #outLog_path = pathtoFastqs + "/Standard_out_log"
    
    # set the aligner to be bwa-mem if nothing was input
    if not args.aligner:
        aligner = 'bwa-mem'
        indexer = "bwa"
    elif args.aligner == 'bwa-mem':
        aligner = 'bwa-mem'
        indexer = "bwa"
    elif args.aligner == 'bowtie2':
        aligner = 'bowtie2'
        indexer = "bowtie2-build"
    else:
        print("Alinger must be 'bwa-mem' or 'bowtie2'")
    
    with open(Log_path, "w") as Log:
    
        if args.fasta:
            fasta = args.fasta
            if not os.path.exists(fasta):
                print(fasta, "does not exist!")
                exit(1)
            if not args.index:
                print("Since you included your own fasta reference file with the '-f' argument, you must indicate if you need it indexed with bwa using the '-i' argument")
                exit(5)

            if args.index == 'True':
                print("Indexing user provided reference fasta file with " + indexer + "...")
                if indexer == "bwa":
                    index_command = "bwa index " + fasta
                elif indexer == "bowtie2-build":
                    index_command = "bowtie2-build " + fasta + " " + fasta
                index_run = subprocess.run(index_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
                Log.write(indexer + " log:\n" + index_run.stderr.decode() + "\n" + index_run.stdout.decode() + "\n" ) # dont think I actually need the stdout section as both out and error are same thing it seems
                #outLog.write(indexer + " standard out:\n" + index_run.stdout.decode() + "\n")
                if index_run.returncode != 0:
                    print("Problem indexing with " + indexer)
                    exit(2)
        else:
            if not args.genome:
                print("No reference genome or preferred genome was provided. Either use the '-f' argument if a reference fasta is present, or the '-g' argument to download a specific fasta from Ensembl.")
                exit(4)
            print("No fasta reference was provided, so the", genome, "will be downloaded from Ensembl.")
            print("Indexing reference fasta file from Ensembl with bwa...")
            fasta = downloadFasta(genome) # this will download and return the path to the fasta file
            if indexer == "bwa":
                index_command = "bwa index " + fasta
            elif indexer == "bowtie2-build":
                index_command = "bowtie2-build " + fasta + " " + fasta
            index_run = subprocess.run(index_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write(indexer + " log:\n" + index_run.stderr.decode() + "\n" + index_run.stdout.decode() + "\n")
            #outLog.write(indexer + " standard out:\n" + index_run.stdout.decode() + "\n")
            if index_run.returncode != 0:
                print("Problem indexing with " + indexer)
                exit(2)

        
        
        print("Generating FastQC reports for fastq file and filtering low quality reads...")
        os.mkdir(args.pathToFastqs + '/fastQC_results')
        filtered_paths = []
        fastqs = os.listdir(args.pathToFastqs)
        for file in fastqs:
            if file.endswith(".fq") or file.endswith("fastq"):
                # run fastqc
                fastqc_command = "~/Tools2/FastQC/fastqc " + args.pathToFastqs + '/' + file + ' -q -o ' +  args.pathToFastqs + '/fastQC_results/'
                fastqc_run = subprocess.run(fastqc_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
                Log.write("FastQC log:\n" + fastqc_run.stderr.decode() + "\n" + fastqc_run.stdout.decode() + "\n")
                #outLog.write("FastQC standard out:\n" + fastqc_run.stdout.decode() + "\n")
                if fastqc_run.returncode != 0:
                    print("Problem with FastQC")
                    exit(2)
                # use function above to filter fastq
                filtered_fastq = filterFastq(args.pathToFastqs + '/' + file) # this will filter and return the path to the filtered fastq
                filtered_paths.append(filtered_fastq)
            
        
        print("Aligning with " + aligner + "...")
        os.mkdir(args.pathToFastqs + '/bam_files')
        if args.remove_duplicates  == 'TRUE' or args.remove_duplicates == 'True':
                print('Duplicates will be removed...')
        else:
                print("Duplicates will NOT be removed, but will be flagged in the bam file...")
        bam_files = []
        for fastq in filtered_paths:
            fastq_basename = os.path.basename(fastq)
            fastq_basename_noend = ".".join(fastq_basename.split(".")[0:-1]) # this will split on period, remove the last thing, then rejoin (in case there is more than one period)
            
            if aligner ==  'bwa-mem':
                align_command = "bwa mem " + fasta + " " + fastq + " | samtools view -b - > " + args.pathToFastqs + '/bam_files/' + fastq_basename_noend + ".bam"
            else:
                align_command = "bowtie2 -x " + fasta + " -U " + fastq + " | samtools view -b - > " + args.pathToFastqs + '/bam_files/' + fastq_basename_noend + ".bam"
            align_run = subprocess.run(align_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write(aligner + " log:\n" + align_run.stderr.decode() + "\n" + align_run.stdout.decode() + "\n")
            #outLog.write(aligner + " standard out:\n" + align_run.stdout.decode() + "\n")
            if align_run.returncode != 0:
                print("Problem aligning with " + aligner + "!")
                exit(2)

            samtools_sort_command = "samtools sort " + args.pathToFastqs + '/bam_files/' + fastq_basename_noend + ".bam -o " + args.pathToFastqs + '/bam_files/' + fastq_basename_noend + "_sorted.bam"
            samtools_sort_run = subprocess.run(samtools_sort_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write("samtools sort log:\n" + samtools_sort_run.stderr.decode() + "\n" + samtools_sort_run.stdout.decode() + "\n")
            #outLog.write("samtools sort bamCoverage standard out:\n" + samtools_sort_run.stdout.decode() + "\n")
            if samtools_sort_run.returncode != 0:
                print("Problem sorting BAM")
                exit(2)
            os.remove(args.pathToFastqs + '/bam_files/' + fastq_basename_noend + ".bam")
            
            
            # finding the number of duplicates and removing if specified               
               
            if args.remove_duplicates  == 'TRUE' or args.remove_duplicates == 'True':
                picard_dup_command = "java -jar ~/Tools2/picard.jar MarkDuplicates I=" + args.pathToFastqs + '/bam_files/' + fastq_basename_noend + "_sorted.bam O=" + args.pathToFastqs + '/bam_files/'+ fastq_basename_noend + "_sorted_rmdup.bam M=" + args.pathToFastqs + '/bam_files/' + fastq_basename_noend + "_sorted.metrics REMOVE_DUPLICATES=TRUE"
                picard_dup_run = subprocess.run(picard_dup_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
                Log.write("MarkDuplicates log:\n" + picard_dup_run.stderr.decode() + "\n" + picard_dup_run.stdout.decode() + "\n")
                #outLog.write("MarkDuplicates standard out:\n" + picard_dup_run.stdout.decode() + "\n")
                if picard_dup_run.returncode != 0:
                    print("Problem flagging duplicates")
                    exit(2)
                bam_files.append(args.pathToFastqs + '/bam_files/'+ fastq_basename_noend + "_sorted_rmdup.bam")
            else:
                picard_dup_command = "java -jar ~/Tools2/picard.jar MarkDuplicates I=" + args.pathToFastqs + '/bam_files/' + fastq_basename_noend + "_sorted.bam O=" + args.pathToFastqs + '/bam_files/'+ fastq_basename_noend + "_sorted_mdup.bam M=" + args.pathToFastqs + '/bam_files/' + fastq_basename_noend + "_sorted.metrics"
                picard_dup_run = subprocess.run(picard_dup_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
                Log.write("MarkDuplicates log:\n" + picard_dup_run.stderr.decode() + "\n" + picard_dup_run.stdout.decode() + "\n")
                #outLog.write("MarkDuplicates standard out:\n" + picard_dup_run.stdout.decode() + "\n")
                if picard_dup_run.returncode != 0:
                    print("Problem flagging duplicates")
                    exit(2)
                bam_files.append(args.pathToFastqs + '/bam_files/' + fastq_basename_noend + "_sorted_mdup.bam")
                

            samtools_index_command = "samtools index " + bam_files[len(bam_files) - 1]
            samtools_index_run = subprocess.run(samtools_index_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write("samtools index log:\n" + samtools_index_run.stderr.decode() + "\n" + samtools_index_run.stdout.decode() + "\n")
            #outLog.write("samtools index standard out:\n" + samtools_index_run.stdout.decode() + "\n")
            if samtools_index_run.returncode != 0:
                print("Problem indexing BAM")
                exit(2)

        '''
        print("Making bigWigs normalized to counts per million...")
        os.mkdir(args.pathToFastqs + '/bigwigs')
        for bam_file in bam_files:
            bam_basename = os.path.basename(bam_file)
            bam_basename_noend = ".".join(bam_basename.split(".")[0:-1]) # this will split on period, remove the last thing, then rejoin (in case there is more than one period)
            
            bigwig_command = "bamCoverage -b " + bam_file + " -o " + args.pathToFastqs + '/bigwigs/' + bam_basename_noend + ".bw --normalizeUsing CPM"
            bigwig_run = subprocess.run(bigwig_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write("deepTools bamCoverage log:\n" + bigwig_run.stderr.decode() + "\n" + bigwig_run.stdout.decode() + "\n")
            #outLog.write("deepTools bamCoverage standard out:\n" + bigwig_run.stdout.decode() + "\n")
            if bigwig_run.returncode != 0:
                print("Problem making bigwig")
                exit(2)
        '''
        
        if args.pcaCorr == 'TRUE' or args.pcaCorr == "True":
            # PCA and correlations with deeptools
            print("Performing multiBamSummary, PCA, correlation with BAM files...")
            os.mkdir(args.pathToFastqs + '/bam_deepTools_pcaCorr')

            bam_files_string = " ".join(bam_files)
            now = datetime.today().isoformat()
            now = now.replace("-", "_")
            now = now.replace(":", ".")
            multiBamSummary_command = "multiBamSummary bins -b " + bam_files_string + " -o " + args.pathToFastqs + "/bam_deepTools_pcaCorr/multiBamSummary" + now + ".npz" 
            multiBamSummary_run = subprocess.run(multiBamSummary_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write("deepTools multiBamSummary log:\n" + multiBamSummary_run.stderr.decode() + "\n" + multiBamSummary_run.stdout.decode() + "\n")
            if multiBamSummary_run.returncode != 0:
                print("Problem with deeptools multiBamSummary")
                exit(2)

            plotPCA_command = "plotPCA --corData " + args.pathToFastqs + "/bam_deepTools_pcaCorr/multiBamSummary" + now + ".npz" + " --plotFile " + args.pathToFastqs + "/bam_deepTools_pcaCorr/multiBamSummary" + now + "_PCA.pdf" 
            plotPCA_run = subprocess.run(plotPCA_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write("deepTools plotPCA log:\n" + plotPCA_run.stderr.decode() + "\n" + plotPCA_run.stdout.decode() + "\n")
            if plotPCA_run.returncode != 0:
                print("Problem with deeptools plotPCA")
                exit(5)

            plotCorrelation_command = "plotCorrelation --corData " + args.pathToFastqs + "/bam_deepTools_pcaCorr/multiBamSummary" + now + ".npz" + " --corMethod spearman --whatToPlot heatmap --plotFile " + args.pathToFastqs + "/bam_deepTools_pcaCorr/multiBamSummary" + now + "_corr.pdf" 
            plotCorrelation_run = subprocess.run(plotCorrelation_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write("deepTools plotCorrelation log:\n" + plotCorrelation_run.stderr.decode() + "\n" + plotCorrelation_run.stdout.decode() + "\n")
            if plotCorrelation_run.returncode != 0:
                print("Problem with deeptools plotCorrelation")
                exit(6)

        
        
        # peak calling
        if args.sampleSheet:
            print("Calling peaks...")
            if not args.genome:
                print("No reference genome or preferred genome was provided. Either use the '-f' argument if a reference fasta is present, or the '-g' argument to download a specific fasta from Ensembl.")
                exit(4)
            os.mkdir(args.pathToFastqs + '/peaks')
            ss = pd.read_csv(args.sampleSheet)
            if args.peakCaller == "macs2":
                if genome == "hg19" or genome == "hg38":
                    macs2_genome = "hs"
                if genome == "mm9" or genome == "mm10":
                    macs2_genome = "mm"
                        
                for index, row in ss.iterrows():
                    input_noend = ".".join(row[2].split(".")[0:-1])
                    for bam_file in bam_files:
                        if re.search(input_noend, bam_file): # this will be a probelm when the fastq names are just longer versions of other fastqs. Just need to make sure the names are all unique without any substrings of other fastq names
                            input_bam = bam_file
                    peakfq_noend = ".".join(row[1].split(".")[0:-1])
                    for bam_file in bam_files:
                        if re.search(peakfq_noend, bam_file): # this will be a probelm when the fastq names are just longer versions of other fastqs. Just need to make sure the names are all unique without any substrings of other fastq names
                            peak_bam = bam_file
                    
                    if not row[3]:
                        os.mkdir(args.pathToFastqs + '/peaks/' + row[0] + "_macs2_narrow/")
                        output_path = args.pathToFastqs + '/peaks/' + row[0] + "_macs2_narrow/"
                        macs_command = "macs2 callpeak -t " + peak_bam + " -c " + input_bam + " -g " + macs2_genome + " -n " + row[0] + " --outdir " + output_path + " " + args.extraPeakArgs

                    else:
                        os.mkdir(args.pathToFastqs + '/peaks/' + row[0] + "_macs2_broad/")
                        output_path = args.pathToFastqs + '/peaks/' + row[0] + "_macs2_broad/"
                        macs_command = "macs2 callpeak -t " + peak_bam + " -c " + input_bam + " -g " + macs2_genome + " -n " + row[0] + " --outdir " + output_path + " --broad " + args.extraPeakArgs
                    macs_run = subprocess.run(macs_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
                    Log.write("samtools index log:\n" + macs_run.stderr.decode() + "\n" + samtools_index_run.stdout.decode() + "\n")
                    if macs_run.returncode != 0:
                        print("Problem calling peaks with macs2 for peak file", peak_bam, "and input file", input_bam)
            #if args.peakCaller == "homer":
                # want to build this in, but is kind of annoying becuase it seems like it works best if you convert to bed, then make the tag directory, then call peaks. Not hard, just a bit annoying and this isnt necessary right away. Will do more QC stuff first. 
                
                
        else:
            print("No sample sheet included, peaks will not be called.")

        
            
if __name__ == '__main__':
	main()
    
    
    
