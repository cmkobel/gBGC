from gwf import *



goodbye = "echo JOBID $SLURM_JOBID; jobinfo $SLURM_JOBID"



def smalt_index(target_dir, group, reference_stem):

    inputs = 'reference_genomes/' + reference_stem + '.fasta'
    outputs = ['reference_genomes/' + reference_stem + '.smi',
               'reference_genomes/' + reference_stem + '.sma',
               'reference_genomes/' + reference_stem + '.cds.gff3']
    options = {'nodes': 1, 'cores': 1, 'memory': '8g', 'walltime': '1:00:00',  'account': 'clinicalmicrobio'}
    
    spec = f"""

# mkdir {target_dir}/group
# cd {target_dir}/group

cd reference_genomes/

smalt index -k 14 -s 8 {reference_stem} {reference_stem}.fasta

# Later, when we are going to extract the CDS sequences from each isolate, we are only interested in the CDS.
cat {reference_stem}.gff3 | awk '$3 == "CDS" {{print}}' > {reference_stem}.cds.gff3




{goodbye}
    """
    return inputs, outputs, options, spec



def smalt_map(group, sample, forward, reverse, reference_stem):
    inputs = ['reference_genomes/' + reference_stem + '.smi',
               'reference_genomes/' + reference_stem + '.sma']
    outputs = ['output/' + group + '/' + sample + '.sorted.bam']
    options = {'nodes': 1, 'cores': 8, 'memory': '32g', 'walltime': '2:00:00',  'account': 'clinicalmicrobio'}
    spec = f"""

mkdir -p output/{group}
cd output/{group}


# map reads to reference and sort
smalt map -f bam -n 4 ../../reference_genomes/{reference_stem} ../../{forward} ../../{reverse} | \
samtools sort > {sample}.sorted.bam



{goodbye}
    """
    return inputs, outputs, options, spec



def coverage(group, sample, forward, reverse, reference_stem):
    inputs = ['output/' + group + '/' + sample + '.sorted.bam']
    outputs = ['output/' + group + '/' + sample + '.sorted.filtered.bam']
    options = {'nodes': 1, 'cores': 8, 'memory': '32g', 'walltime': '2:00:00',  'account': 'clinicalmicrobio'}
    spec = f"""

mkdir -p output/{group}
cd output/{group}


#cat {sample}.bam | /home/cmkobel/software/sambamba71/sambamba-0.7.1-linux-static view -f bam -F "proper_pair" -S -t 8 /dev/stdin > {sample}.proppair.bam
#todo: overvej at bruge cram (sambamba view -f bam -T reference.fasta > out.bam)


# Mark duplicates in the alignment
/home/cmkobel/software/sambamba71/sambamba-0.7.1-linux-static markdup -t 8 {sample}.sorted.bam --tmpdir=/scratch/$GWF_JOBID/ {sample}.markdup.bam

# Remove low-quality content.
sambamba view -F "not (duplicate or secondary_alignment or unmapped)" -f bam {sample}.markdup.bam > {sample}.sorted.filtered.bam


samtools depth {sample}.sorted.filtered.bam > {sample}_cov.tab
# TODO use samtools bedcov and gff2bed instead of this whole R data wrangling nightmare


# clean up
#rm {sample}.markdup.bam
#rm {sample}.markdup.bam.bai



{goodbye}
    """
    return inputs, outputs, options, spec





def consensus(group, sample, reference_stem):
    inputs = ['output/' + group + '/' + sample + '.sorted.bam']
    outputs = ['output/' + group + '/' + sample + '.fasta']
    options = {'nodes': 1, 'cores': 4, 'memory': '16g', 'walltime': '2:00:00',  'account': 'clinicalmicrobio'}
    spec = f"""

mkdir -p output/{group}
cd output/{group}


samtools mpileup -f ../../reference_genomes/{reference_stem}.fasta {sample}.sorted.bam | \
GenomicConsensus --fna_file ../../reference_genomes/{reference_stem}.fasta > {sample}.fasta 


{goodbye}
    """
    return inputs, outputs, options, spec



def mcorr_bam_fit(group, sample, forward, reverse, reference_stem):
    # todo: split op mellem mcorr_bam og mcorr_fit. Mange ting kan laves uden at fitte. 
    inputs = ['reference_genomes/' + reference_stem + '.cds.gff3',
              'output/' + group + '/' + sample + '.sorted.bam']
    outputs = ['output/' + group + '/' + sample + '.csv',  # mcorr_bam
               'output/' + group + '/' + sample + '.json', # mcorr_bam
               'output/' + group + '/' + sample + '_best_fit.svg',  # mcorr_fit
               'output/' + group + '/' + sample + '_fit_reports.txt',  # mcorr_fit
               'output/' + group + '/' + sample + '_fit_results.csv']  # mcorr_fit
               #'output/' + group + '/' + sample + '_parameter_histograms.svg']  # mcorr_fit # fejl muligvis relateret til -t -T i samtools sort
    options = {'nodes': 1, 'cores': 4, 'memory': '64g', 'walltime': '2:00:00',  'account': 'clinicalmicrobio'}
    spec = f"""
mkdir -p output/{group}
cd output/{group}

mcorr-bam ../../reference_genomes/{reference_stem}.cds.gff3 {sample}.sorted.bam {sample}

mcorr-fit {sample}.csv {sample} || echo "probably not enough memory"



{goodbye}
"""
    return inputs, outputs, options, spec




def extract_fasta_from_bed(group, sample, reference_stem):
    # todo: overvej at skrive denne target sammen med 'consensus'
    inputs = ['output/' + group + '/' + sample + '.fasta',
              'reference_genomes/' + reference_stem + '.gff3']
    outputs = ['output/' + group + '/' + sample + '.extracted.fasta',
               'output/' + group + '/' + sample + '.gc.tab']
    options = {'nodes': 1, 'cores': 4, 'memory': '8g', 'walltime': '2:00:00',  'account': 'clinicalmicrobio'}
    spec = f"""
mkdir -p output/{group}
cd output/{group}

# The index is automatically generated
bedtools getfasta -fi {sample}.fasta -bed ../../reference_genomes/{reference_stem}.gene.gff3 > {sample}.extracted.fasta

infoseq -only -name -pgc {sample}.extracted.fasta | awk '{{print $1 "\t" $2}}' > {sample}.gc.tab


{goodbye}
"""
    return inputs, outputs, options, spec



