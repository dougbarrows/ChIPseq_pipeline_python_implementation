#!/usr/bin/env python3

from Bio import SeqIO
import os, sys, re, glob, shutil
import subprocess
import gzip
import argparse
import pandas as pd
from datetime import datetime
import multiprocessing as mp
from itertools import chain

# other things to add:
## read trimming
## input gene list for ngsplot?
## peak calling - homer and macs2?
## make the log separate for each sample??

# dependencies - picard, fastqc, bwa, bowtie2, samtools, deeptools, macs2


# this function is largely from Menas function from CSHL course, only thing I changed was to make the negative strand genes TSS be the 'end'
def gff_to_TSS(gff3filename):
    # Determine the output file name
    outputfilename = "./"+os.path.splitext(gff3filename)[0]+".TSS.bed"

    # fields in a gff3 files
    field_str = "seqid source genetype start end score strand phase attributes"
    fields = field_str.split(' ') # makes a list of the columns in the gff file


    # open gff3 file, read in lines with gene information
    with open(gff3filename, 'r') as fo, open(outputfilename, 'w') as outputfile:

        for line in fo:
            line = line.strip()

            # skip header lines
            if line.startswith('#'):
                continue
            # get info from non-header lines
            else:
                data = dict(zip(fields,line.split('\t')))	# dict(zip()) is a clever way to make a dictionary from two lists with the keys as list 1 and values as list 2	
                # only keep lines that are of type "gene"
                # print relevant columns to file

                if data['genetype'] == 'gene':

                    attributes = data['attributes'].split(";")
                    gene_name = attributes[0]+";"+attributes[1]
                    strand = data['strand']

                    bed_line = "chr{}\t{}\t{}\t{}\t{}\t{}\n"
                    if data['strand'] == "-":
                        outputfile.write(bed_line.format(data['seqid'], int(data['end']), int(data['end']) + 1, gene_name, '', strand))
                    else:
                        outputfile.write(bed_line.format(data['seqid'], int(data['start']), int(data['start']) + 1, gene_name, '', strand))
        return(outputfilename)
    
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
       
def filterFastq(fastqPath):
    if fastqPath.endswith(".fq"):
        filteredPath = fastqPath[:-3] + "_filtered.fq"
    if fastqPath.endswith(".fq.gz"):
        filteredPath = fastqPath[:-6] + "_filtered.fq"
    if fastqPath.endswith(".fastq"):
        filteredPath = fastqPath[:-6] + "_filtered.fastq"
    if fastqPath.endswith(".fastq.gz"):
        filteredPath = fastqPath[:-9] + "_filtered.fastq"
    if fastqPath.endswith(".fq.gz") or fastqPath.endswith(".fastq.gz"):
        with gzip.open(r"{}".format(fastqPath), "rt") as fo, open(r"{}".format(filteredPath), "w") as filter_fo:
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
    if fastqPath.endswith(".fq") or fastqPath.endswith(".fastq"):
        with open(r"{}".format(fastqPath), "rt") as fo, open(r"{}".format(filteredPath), "w") as filter_fo:
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
     
def fastqc_and_filter(path_to_each_fastq):
    fastqc_command = "~/Tools2/FastQC/fastqc " + path_to_each_fastq + ' -q -o ' +  os.path.dirname(path_to_each_fastq) + '/fastQC_results/'
    fastqc_run = subprocess.run(fastqc_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    for_log = "FastQC log for " + os.path.basename(path_to_each_fastq) + ":\n" + fastqc_run.stderr.decode() + "\n" + fastqc_run.stdout.decode() + "\n"
    #outLog.write("FastQC standard out:\n" + fastqc_run.stdout.decode() + "\n")
    if fastqc_run.returncode != 0:
        print("Problem with FastQC")
        #return("FastQC log:\n" + fastqc_run.stderr.decode() + "\n" + fastqc_run.stdout.decode() + "\n")
    # use function above to filter fastq
    filtered_fastq = filterFastq(path_to_each_fastq) # this will filter and return the path to the filtered fastq
    return(filtered_fastq, for_log) # just for clarity sake I make another line to show I am returning the path to the filtered fastq
    
    
def align_sort_markDups_index_bam(fastq, aligner, fasta, path_to_fastq_dir, remove_dups, cores):

    #for fastq in filtered_paths:
    fastq_basename = os.path.basename(fastq)
    fastq_basename_noend = ".".join(fastq_basename.split(".")[0:-1]) # this will split on period, remove the last thing, then rejoin (in case there is more than one period)
            
    if aligner ==  'bwa-mem':
        align_command = "bwa mem " + fasta + " " + fastq + " -t " + str(cores) + " | samtools view -b - > " + path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + ".bam"
    else:
        align_command = "bowtie2 -x " + fasta + " -U " + fastq + " -p " + str(cores) + " | samtools view -b - > " + path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + ".bam"
    align_run = subprocess.run(align_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    for_log_align = aligner + " log for " + fastq_basename + ":\n" + align_run.stderr.decode() + "\n" + align_run.stdout.decode() + "\n"
    if align_run.returncode != 0:
        print("Problem aligning with " + aligner + ": \n" + align_run.stderr.decode() + "\n" + align_run.stdout.decode() + "\n")
        exit(2)

    samtools_sort_command = "samtools sort " + path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + ".bam -o " + path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + "_sorted.bam -@ " + str(cores)
    samtools_sort_run = subprocess.run(samtools_sort_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    for_log_samtools_sort = "samtools sort log for " + fastq_basename + ":\n" + samtools_sort_run.stderr.decode() + "\n" + samtools_sort_run.stdout.decode() + "\n"
    if samtools_sort_run.returncode != 0:
        print("Problem sorting BAM : \n" + samtools_sort_run.stderr.decode() + "\n" + samtools_sort_run.stdout.decode() + "\n")
        exit(2)
    os.remove(path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + ".bam")
            
            
    # finding the number of duplicates and removing if specified               
               
    if remove_dups:
        picard_dup_command = "java -jar ~/Tools2/picard.jar MarkDuplicates I=" + path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + "_sorted.bam O=" + path_to_fastq_dir + '/bam_files/'+ fastq_basename_noend + "_sorted_rmdup.bam M=" + path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + "_sorted.metrics REMOVE_DUPLICATES=TRUE"
        picard_dup_run = subprocess.run(picard_dup_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
        for_log_mark_dups = "MarkDuplicates log for " + fastq_basename + ":\n" + picard_dup_run.stderr.decode() + "\n" + picard_dup_run.stdout.decode() + "\n"
        if picard_dup_run.returncode != 0:
            print("Problem flagging duplicates \n" + picard_dup_run.stderr.decode() + "\n" + picard_dup_run.stdout.decode() + "\n")
            exit(2)
        bam_files = path_to_fastq_dir + '/bam_files/'+ fastq_basename_noend + "_sorted_rmdup.bam"
    else:
        picard_dup_command = "java -jar ~/Tools2/picard.jar MarkDuplicates I=" + path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + "_sorted.bam O=" + path_to_fastq_dir + '/bam_files/'+ fastq_basename_noend + "_sorted_mdup.bam M=" + path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + "_sorted.metrics"
        picard_dup_run = subprocess.run(picard_dup_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
        for_log_mark_dups = "MarkDuplicates log for " + fastq_basename + ":\n" + picard_dup_run.stderr.decode() + "\n" + picard_dup_run.stdout.decode() + "\n"
        if picard_dup_run.returncode != 0:
            print("Problem flagging duplicates \n" + picard_dup_run.stderr.decode() + "\n" + picard_dup_run.stdout.decode() + "\n")
            exit(2)
        bam_files = path_to_fastq_dir + '/bam_files/' + fastq_basename_noend + "_sorted_mdup.bam"
                

    samtools_index_command = "samtools index " + bam_files
    samtools_index_run = subprocess.run(samtools_index_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    for_log_samtools_index = "samtools index log for " + fastq_basename + ":\n" + samtools_index_run.stderr.decode() + "\n" + samtools_index_run.stdout.decode() + "\n"
    if samtools_index_run.returncode != 0:
        print("Problem indexing BAM")
        exit(2)
                
    return(bam_files, for_log_align, for_log_samtools_sort, for_log_mark_dups, for_log_samtools_index)

def make_bigwig(bam_file, path_to_fastq_dir):
    bam_basename = os.path.basename(bam_file)
    bigwigs = []
    bam_basename_noend = ".".join(bam_basename.split(".")[0:-1]) # this will split on period, remove the last thing, then rejoin (in case there is more than one period)
    bigwigs.append(path_to_fastq_dir + '/bigwigs/' + bam_basename_noend + ".bw") # make a list of bigwig path to use later on in deeptools computeMatrix
    bigwig_command = "bamCoverage -b " + bam_file + " -o " + path_to_fastq_dir + '/bigwigs/' + bam_basename_noend + ".bw --normalizeUsing CPM"
    bigwig_run = subprocess.run(bigwig_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    for_log_bamCoverage = "deepTools bamCoverage log for " + bam_basename + ":\n" + bigwig_run.stderr.decode() + "\n" + bigwig_run.stdout.decode() + "\n"
    if bigwig_run.returncode != 0:
        print("Problem making bigwig")
        exit(2)
    return(bigwigs, for_log_bamCoverage)
    
def main():
    parser = argparse.ArgumentParser(description="""
        Pipeline for ChIPseq 
        """)
    parser.add_argument('-f','--fasta', help='The name of the fasta file for the reference genome (eg hg19.fa). This can be a full path if the file is not in the current working directory.')
    parser.add_argument('-a','--aligner', choices = ['bwa-mem', 'bowtie2'], help='The aligner to use, bwa-mem will be used by default')
    parser.add_argument('-q','--pathToFastqs', help='The path to the folder that contains the fastq files.')
    parser.add_argument('-g', '--genome' , choices = ['hg38', 'hg19', 'mm10'], help="This is required if no fasta reference file is provided in the '-f' argument, if peaks are to be called, or if '-t' is True. This will determine which reference genome fasta is downloaded. This is also required if peak calling is desired")
    parser.add_argument('-i', '--index', choices = ['True', 'False'], help='If a reference fasta is provided, please indicate whether you need it to be indexed with BWA')
    parser.add_argument('-d', '--remove_duplicates', choices = ['True', 'False'], help='This will determine if duplicates should be removed, by defualt they will be flagged in the bam file, but not removed.')
    parser.add_argument('-s', '--sampleSheet', help='If peaks should be called, the path to the peak sample sheet should be provided here. If no path is included, peaks will not be called.')
    parser.add_argument('-p', '--peakCaller', choices = ['macs2', 'homer'], help='The peak caller to be used.')
    parser.add_argument('-P', '--extraPeakArgs', help="A string to be added onto either macs2 or homer command in addition to the standard arguments. Make sure to put in quotes. E.g. -x '--bw 300 -q 0.01'")
    parser.add_argument('-c', '--pcaCorr', help="Perform PCA and correlation analysis with deepTools")
    parser.add_argument('-t', '--tssPlot', choices = ['True', 'False'], help="Output a signal plot over TSS.")
    parser.add_argument('-n', '--ngsPlot', choices = ['tss', 'tes', 'genebody', 'exon'], nargs = '*', help="Output a signal plot over the specified region type.")
    parser.add_argument('-N', '--ngsPlotExtra', help="A string to be added onto 'ngs.plot.r' command within the ngsplot package. The arguments (ngs.plot.r arguments, not from this pipeline) already included are -G, -R, -C, -O.  Make sure to put in quotes. E.g. -N ' -D ensembl -FL 300'")
    parser.add_argument('-w', '--numCores', type=int, help="Number of cores/workers to use in parallelizing code")

    args = parser.parse_args()
    
    path_to_fastq_dir = args.pathToFastqs
    
    if args.remove_duplicates:
        remove_dups = True
        
    if args.genome:
        genome = args.genome
    
    if args.numCores:
        cores = args.numCores
    else:
        cores = 1
        
    Log_path = path_to_fastq_dir + "/StdError_Stdout_log.txt"
    
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
            if index_run.returncode != 0:
                print("Problem indexing with " + indexer)
                exit(2)

        
        fastqs = os.listdir(args.pathToFastqs)
        fastqs_to_analyze = []
        for file in fastqs:
            if file.endswith("fq") or file.endswith("fastq") or file.endswith("fq.gz") or file.endswith("fastq.gz"):
                fastqs_to_analyze.append(file)
        if len(fastqs_to_analyze) == 0:
            print("no fastq files detected, make sure the fastq files end with either 'fq', 'fastq', 'fq.gz', 'fastq.gz'")
            exit(3)
        else:
            print("fastqs put through pipeline: ")
            for fastq in fastqs_to_analyze:
                print(fastq)
                
        print("\nGenerating FastQC reports for fastq files and filtering low quality reads...")
        os.mkdir(args.pathToFastqs + '/fastQC_results')
        path_to_each_fastq = []
        for fastq in fastqs_to_analyze:
            path_to_each_fastq.append( args.pathToFastqs + "/" + fastq)
            
        filtered_paths = []
        
        pool = mp.Pool(cores)
        filtered_paths, for_log = zip(*pool.map(fastqc_and_filter, path_to_each_fastq))# this outputs two tuples, one with the path to filtered fastqs and one with the log for each
        pool.close()
        
        Log.write("\n".join(for_log))                
        
        print("Aligning with " + aligner + "...")
        os.mkdir(args.pathToFastqs + '/bam_files')
        if args.remove_duplicates:
                print('Duplicates will be removed...')
        else:
                print("Duplicates will NOT be removed, but will be flagged in the bam file...")
        
        '''
        filtered_paths = list(filtered_paths) # this is a tuple from the multiprocess before, so convert to a list
        aligner_list = ([aligner] * len(filtered_paths))
        fasta_list = ([fasta] * len(filtered_paths))
        path_to_fastq_dir_list = ([path_to_fastq_dir] * len(filtered_paths))
        remove_dups_list = ([remove_dups] * len(filtered_paths))
        cores_list = ([cores] * len(filtered_paths))
        
        # the arguments for starmap must be a list of lists, with each list contiang the arguments for each iteration through the function used in starmap
        input_list = []
        for i in range(len(filtered_paths)):
            input_list.append([filtered_paths[i], 
                             aligner_list[i], 
                             fasta_list[i], 
                             path_to_fastq_dir_list[i], 
                             remove_dups_list[i],
                             cores_list[i]])
        
        pool = mp.Pool(1)
        bam_files, for_log_align, for_log_samtools_sort, for_log_mark_dups, for_log_samtools_index = zip(*pool.starmap(align_sort_markDups_index_bam, input_list))
        pool.close()
        '''
        
        bam_files_list = []
        '''
        for_log_align_list = []
        for_log_samtools_sort_list = []
        for_log_mark_dups_list = []
        for_log_samtools_index_list = []
        '''
            
        for fastq in filtered_paths:
            bam_files, for_log_align, for_log_samtools_sort, for_log_mark_dups, for_log_samtools_index = align_sort_markDups_index_bam(fastq, aligner, fasta, path_to_fastq_dir, remove_dups, cores)

            bam_files_list.append(bam_files)
            Log.write(for_log_align + "\n") 
            Log.write(for_log_samtools_sort + "\n")
            Log.write(for_log_mark_dups + "\n")
            Log.write(for_log_samtools_index + "\n")
            
            
            
            '''
            for_log_align_list = for_log_align_list.append()
            for_log_samtools_sort_list = for_log_samtools_sort_list.append()
            for_log_mark_dups_list = for_log_mark_dups_list.append()
            for_log_samtools_index_list = for_log_samtools_index_list.append()
            
       
        Log.write("\n".join(for_log_align_list)) 
        Log.write("\n".join(for_log_samtools_sort_list))
        Log.write("\n".join(for_log_mark_dups_list))
        Log.write("\n".join(for_log_samtools_index_list)) 

        
        bam_files_list = list(bam_files)
        bam_files_list = list(chain.from_iterable(bam_files_list))
        # this comes out as as tuple from the above multiprocess command, so convert to a list to iterate below
        '''

        print("Making bigWigs normalized to counts per million...")
        os.mkdir(args.pathToFastqs + '/bigwigs')

        
        path_to_fastq_dir_list = ([path_to_fastq_dir] * len(filtered_paths))
        
        input_list_bw = []
        for i in range(len(bam_files_list)):
            input_list_bw.append([bam_files_list[i],
                                  path_to_fastq_dir_list[i]])
                    
        pool = mp.Pool(cores)
        bigwigs, for_log_bamCoverage = zip(*pool.starmap(make_bigwig, input_list_bw))
        pool.close()
        
        Log.write("\n".join(for_log_bamCoverage)) 
        
        if args.pcaCorr == 'TRUE' or args.pcaCorr == "True":
            # PCA and correlations with deeptools
            print("Performing multiBamSummary, PCA, correlation with BAM files...")
            os.mkdir(args.pathToFastqs + '/bam_deepTools_pcaCorr')

            bam_files_string = " ".join(bam_files_list)
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
        
        
        if args.tssPlot == 'True' or args.tssPlot == 'TRUE':
            print("Generating plots of signal over TSS with deeptools...")
            os.mkdir(args.pathToFastqs + '/deeptools_TSSplot')
            gff3filename = "Mus_musculus.GRCm38.98.gff3" # this needs to be downloaded eventually (or check if its already there), amybe make an agument with the path?
            tss_bed = gff_to_TSS(gff3filename) # this will write the bed file with the TSS's and also return the path to that file
            
            bigwig_string = " ".join(bigwigs)
            now = datetime.today().isoformat()
            now = now.replace("-", "_")
            now = now.replace(":", ".")
            computeMatrix_command = "computeMatrix reference-point -R " + tss_bed + " -S " + bigwig_string + " -o " + args.pathToFastqs + '/deeptools_TSSplot/computematrix' + now + ".mat" + " --beforeRegionStartLength 1000 --afterRegionStartLength 1000"
            computeMatrix_run = subprocess.run(computeMatrix_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write("deepTools computeMatrix log:\n" + computeMatrix_run.stderr.decode() + "\n" + computeMatrix_run.stdout.decode() + "\n")
            if computeMatrix_run.returncode != 0:
                print("Problem with deeptools computeMatrix")
                exit(6)
            
            plotProfile_indiv_command = "plotProfile -m " + args.pathToFastqs + '/deeptools_TSSplot/computematrix' + now + ".mat" + " -o " + args.pathToFastqs + '/deeptools_TSSplot/individual_plotProfile' + now + ".pdf" 
            plotProfile_indiv_run = subprocess.run(plotProfile_indiv_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write("deepTools plotProfile-individual log:\n" + plotProfile_indiv_run.stderr.decode() + "\n" + plotProfile_indiv_run.stdout.decode() + "\n")
            if plotProfile_indiv_run.returncode != 0:
                print("Problem with deeptools plotProfile (individual)")
                exit(10)   
            
            plotProfile_group_command = "plotProfile -m " + args.pathToFastqs + '/deeptools_TSSplot/computematrix' + now + ".mat" + " -o " + args.pathToFastqs + '/deeptools_TSSplot/group_plotProfile' + now + ".pdf --perGroup" 
            plotProfile_group_run = subprocess.run(plotProfile_group_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            Log.write("deepTools plotProfile-group log:\n" + plotProfile_group_run.stderr.decode() + "\n" + plotProfile_group_run.stdout.decode() + "\n")
            if plotProfile_group_run.returncode != 0:
                print("Problem with deeptools plotProfile (group)")
                exit(10)

                
        if args.ngsPlot:
            print("Generating plot of signal over specified regions with ngsplot...")
            os.mkdir(args.pathToFastqs + '/ngsplots')
            regions = args.ngsPlot
            config_paths = []
            for region in regions:
                config_path = args.pathToFastqs + '/ngsplots/' + region + '_config.txt'
                config_paths.append(config_path)
                with open(config_path, "w") as config_file:
                    for bam_file in bam_files_list: # iterate through the 'bam_files' variable made above
                        bam_base = os.path.basename(bam_file)
                        bam_base_noend = ".".join(bam_base.split(".")[0:-1])                   
                        config_file.write(bam_file + "\t-1\t" + bam_base_noend + "\n")
                if not args.ngsPlotExtra:
                    args.ngsPlotExtra = " "
                ngs_group_command = "ngs.plot.r -G " + genome + " -R " + region + " -C " +  config_path + " -O " + args.pathToFastqs + '/ngsplots/' + region + "_ngsPlot_group " + args.ngsPlotExtra
                ngs_group_run = subprocess.run(ngs_group_command, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
                Log.write("ngsplot-group (" + region + ") log:\n" + ngs_group_run.stderr.decode() + "\n" + ngs_group_run.stdout.decode() + "\n")
                if ngs_group_run.returncode != 0:
                    print("Problem with ngsplot (group)")
                    exit(10)
            
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
    
    
    
