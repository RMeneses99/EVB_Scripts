# Snakefile

import os

# Define reference genome paths
Ref_gens = { 
    "resident_fasta" : "/home/Reference_genomes/CP054662.1_IS_annotated.fasta",
    "resident_gbk" : "/home/Reference_genomes/CP054662.1_IS_annotated.gbk",
    "invader_fasta" : "/home/Reference_genomes/NC_000913.2.fasta",
    "invader_gbk" : "/home/Reference_genomes/NC_000913.2.gbk",
    "rp1_fasta" : "/home/Reference_genomes/CP054663.fasta",
    "rp1_gbk" : "/home/Reference_genomes/CP054663.gbk",
    "rp2_fasta" : "/home/Reference_genomes/CP054664.fasta",
    "rp2_gbk" : "/home/Reference_genomes/CP054664.gbk"
}

# Define input file paths
raw_reads = [
    "IL10KO_SPF_M8_D371_RI_I_L8_R1.fastq.gz", 
    "IL10KO_SPF_M8_D371_RI_I_L8_R2.fastq.gz"
]

# Rule for read trimming and filtering with fastp
rule fastp:
    input:
        raw_reads
    output:
        r1 = "fastp/{basename}".format(basename=os.path.basename(raw_reads[0])),
        r2 = "fastp/{basename}".format(basename=os.path.basename(raw_reads[1]))
    shell:
        """
        fastp -q 20 -u 50 --length_required 100 --dedup 1 --detect_adapter_for_pe -p 3 -5 -M 20 -W 4 -c \
            -i {input[0]} -I {input[1]} -o {output.r1} -O {output.r2} --html /dev/null --json /dev/null
        """

# Rule for contaminant removal with bbsplit
rule bbsplit:
    input:
        r1 = "fastp/{basename}".format(basename=os.path.basename(raw_reads[0])),
        r2 = "fastp/{basename}".format(basename=os.path.basename(raw_reads[1]))
    output:
        directory("decontam")
    params:
        ref_genomes = ",".join(Ref_gens[key] for key in ["resident_fasta", "invader_fasta", "rp1_fasta", "rp2_fasta"]),
        sample_name = lambda wildcards: "_".join(os.path.basename(raw_reads[0]).split("_")[:5])
        
    shell: 
        """
        bbsplit.sh in1={input.r1} in2={input.r2} ambig2=toss \
            ref={params.ref_genomes} basename={output}/{params.sample_name}_%.fastq.gz
        """