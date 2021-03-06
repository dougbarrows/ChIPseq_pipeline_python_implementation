import pandas as pd

ss = pd.read_csv(snakemake.config['peak_samplesheet_path'])
if snakemake.config['genome'] == "hg19" or snakemake.config['genome'] == "hg38":
    macs2_genome = "hs"
if snakemake.config['genome'] == "mm9" or snakemake.config['genome'] == "mm10":
    macs2_genome = "mm"

for index, row in ss.iterrows():
    input_noend = ".".join(row[2].split(".")[0:-1])
    for bam_file in snakemake.input[0]:
        if re.search(input_noend, bam_file): # this will be a probelm when the fastq names are just longer versions of other fastqs. Just need to make sure the names are all unique without any substrings of other fastq names
            input_bam = bam_file
    peakfq_noend = ".".join(row[1].split(".")[0:-1])
    for bam_file in bam_files_list:
        if re.search(peakfq_noend, bam_file): # this will be a probelm when the fastq names are just longer versions of other fastqs. Just need to make sure the names are all unique without any substrings of other fastq names
            peak_bam = bam_file

    if not row[3]:
        os.mkdir(args.pathToFastqs + '/peaks/' + row[0] + "_macs2_narrow/")
        output_path = args.pathToFastqs + '/peaks/' + row[0] + "_macs2_narrow/"
        macs_command = "macs2 callpeak -t " + peak_bam + " -c " + input_bam + " -g " + macs2_genome + " -n " + row[0] + " --outdir " + output_path + " " + args.extraPeakArgs

    else:
        os.mkdir(args.pathToFastqs + '/peaks/' + row[0] + "_macs2_broad/")
