configfile: "config.yaml"

import os
from os.path import basename, join

WD = config['wd']
TMP_DIR = config['tmp_dir']

FASTQ_DIR = join(WD, config['fastq_dir'])
REFGEN_DIR = join(WD, config['refgen_dir'])
GENCODE_DIR = join(WD, config['gencode_dir'])
VCF_DIR = join(WD, config['vcf_dir'])

SNP_DIR = join(WD, config['snp_dir'])
STAR_GENOME_DIR = join(WD, config['star_genome_dir'])
SAM_DIR = join(WD, config['sam_dir'])
READ_GROUPS_DIR = join(WD, config['add_rg_dir'])
MARK_DUPS_DIR = join(WD, config['mark_dups_dir'])
FIND_SNPS_DIR = join(WD, config['find_snps_dir'])
REMAP_DIR = join(WD, config['remap_dir'])
REMAP_RG_DIR = join(WD, config['remap_add_rg_dir'])
REMAP_MD_DIR = join(WD, config['remap_mark_dups_dir'])
FILT_REMAPPED_READS_DIR = join(WD, config['filter_remapped_reads'])
MERGED_DIR = join(WD, config['merged_dir'])
SORTED_DIR = join(WD, config['sorted_dir'])
READ_COUNTS_DIR = join(WD, config['read_counts_dir'])


WC_fastqs = glob_wildcards(join(FASTQ_DIR, '{sample}.{pair}.fq.gz'))
WC_refgens = glob_wildcards(join(REFGEN_DIR, '{genome}.fa'))
WC_gencodes = glob_wildcards(join(GENCODE_DIR, '{gencode}.gtf'))

SAMPLES = set(WC_fastqs.sample)
PAIRS = ['R1', 'R2']

GENOMES = set(WC_refgens.genome)
GENCODES = set(WC_gencodes.gencode)


rule all:
    input:
        expand('{read_counts_dir}/{{sample}}.{{genome}}.tsv'.format(read_counts_dir=READ_COUNTS_DIR),
        sample=SAMPLES, genome=GENOMES)
    run:
        print("ASECOUNT FINISHED WITH NO EXCEPTIONS!")


rule star_genome:
    log: "{star_genome_dir}/{{genome}}.log".format(star_genome_dir=STAR_GENOME_DIR)
    input:
        refgen="{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR),
        gencode="{gencode_dir}/{{genome}}.gtf".format(gencode_dir=GENCODE_DIR)
    output:
        dir="{star_genome_dir}/{{genome}}".format(star_genome_dir=STAR_GENOME_DIR),
        sa="{star_genome_dir}/{{genome}}/SA".format(star_genome_dir=STAR_GENOME_DIR)
    threads: config['star_genome_threads']
    params:
        overhang=config['read_length']-1
    shell:
        "mkdir -p {output.dir} &>> {log};"
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.dir} "
        "--genomeFastaFiles {input.refgen} --sjdbGTFfile {input.gencode} --sjdbOverhang {params.overhang} &>> {log}"


rule star_align:
    log: '{sam_dir}/{{sample}}.{{genome}}.log'.format(sam_dir=SAM_DIR)
    input:
        genome_dir="{star_genome_dir}/{{genome}}".format(star_genome_dir=STAR_GENOME_DIR),
        f1='{fastq_dir}/{{sample}}.R1.fq.gz'.format(fastq_dir=FASTQ_DIR),
        f2='{fastq_dir}/{{sample}}.R2.fq.gz'.format(fastq_dir=FASTQ_DIR),
        gencode="{gencode_dir}/{{genome}}.gtf".format(gencode_dir=GENCODE_DIR)
    output:
        dir='{sam_dir}/{{sample}}.{{genome}}'.format(sam_dir=SAM_DIR),
        sam='{sam_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.Aligned.out.sam'.format(sam_dir=SAM_DIR)
    threads: config['star_align_threads']
    params:
        overhang=config['read_length']-1,
        out_prefix='{sam_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.'.format(sam_dir=SAM_DIR)
    shell:
        "mkdir -p {output.dir} &>> {log}; "
        "STAR --readFilesIn {input.f1} {input.f2} --outFileNamePrefix {params.out_prefix} "
        "--genomeDir {input.genome_dir} --readFilesCommand zcat --runThreadN {threads} "
        "--genomeLoad NoSharedMemory --outFilterType BySJout "
        "--outSAMunmapped Within "
        "--outSAMattributes NH HI AS NM MD NM "
        "--twopassMode Basic "
        "--sjdbOverhang {params.overhang} "
        "--sjdbGTFfile {input.gencode} &>> {log}"


rule add_read_groups:
    log: '{add_rg_dir}/{{sample}}.{{genome}}.log'.format(add_rg_dir=READ_GROUPS_DIR)
    input:
        '{sam_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.Aligned.out.sam'.format(sam_dir=SAM_DIR)
    output:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=READ_GROUPS_DIR)
    params:
        tmp_dir=TMP_DIR
    run:
        command = "{picard} AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=id " \
                  "RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample TMP_DIR={params} &>> {log}; ".format(picard=config['picard_path'],
                  input=input, output=output, params=params.tmp_dir, log=log)

        print(command)
        shell(command)


rule mark_duplicates:
    log: '{mark_dups_dir}/{{sample}}.{{genome}}.log'.format(mark_dups_dir=MARK_DUPS_DIR)
    input:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=READ_GROUPS_DIR)
    output:
        bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=MARK_DUPS_DIR),
        metrics='{mark_dups_dir}/{{sample}}.{{genome}}.metrics'.format(mark_dups_dir=MARK_DUPS_DIR)
    params:
        tmp_dir=TMP_DIR
    run:
        command = "{picard} MarkDuplicates I={input} O={output_bam} CREATE_INDEX=true " \
                  "VALIDATION_STRINGENCY=SILENT M={output_metrics} TMP_DIR={params} &>> {log};".format(picard=config['picard_path'],
                  input=input, output_bam=output.bam, output_metrics=output.metrics, params=params.tmp_dir, log=log)

        print(command)
        shell(command)


rule make_snp_dir:
    input:
        '{vcf_dir}/{{sample}}.{{genome}}.vcf'.format(vcf_dir=VCF_DIR)
    output:
        '{snp_dir}/{{sample}}.{{genome}}'.format(snp_dir=SNP_DIR)
    run:
        import vcf, gzip
        os.makedirs(str(output))
        reader = vcf.Reader(filename=str(input))


        snp_files = {}
        for rec in reader:
            if rec.CHROM not in snp_files.keys():
                snp_files[rec.CHROM] = gzip.open(join(str(output), rec.CHROM+".snps.txt.gz"), 'wt')

            snp_files[rec.CHROM].write(' '.join(map(str, [rec.POS, rec.REF, rec.ALT[0]]))+'\n')

        for chrom, out in snp_files.items():
            out.close()


rule find_intersecting_snps:
    log: '{find_snps_dir}/{{sample}}.{{genome}}.log'.format(find_snps_dir=FIND_SNPS_DIR)
    input:
        snps='{snp_dir}/{{sample}}.{{genome}}'.format(snp_dir=SNP_DIR),
        bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=MARK_DUPS_DIR)
    output:
        dir='{find_snps_dir}/{{sample}}.{{genome}}'.format(find_snps_dir=FIND_SNPS_DIR),
        keep_bam='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.keep.bam'.format(find_snps_dir=FIND_SNPS_DIR),
        remap_bam='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.to.remap.bam'.format(find_snps_dir=FIND_SNPS_DIR),
        fq1='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.fq1.gz'.format(find_snps_dir=FIND_SNPS_DIR),
        fq2='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.fq2.gz'.format(find_snps_dir=FIND_SNPS_DIR),
        single='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.single.fq.gz'.format(find_snps_dir=FIND_SNPS_DIR)
    conda:
        "envs/wasp.yaml"
    shell:
        "python scripts/find_intersecting_snps.py --is_paired_end --is_sorted "
        "--output_dir {output.dir} --snp_dir {input.snps} {input.bam} &>> {log}"


rule star_remap:
    log: '{remap_dir}/{{sample}}.{{genome}}.log'.format(remap_dir=REMAP_DIR)
    input:
        genome_dir="{star_genome_dir}/{{genome}}".format(star_genome_dir=STAR_GENOME_DIR),
        gencode="{gencode_dir}/{{genome}}.gtf".format(gencode_dir=GENCODE_DIR),
        f1='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.fq1.gz'.format(find_snps_dir=FIND_SNPS_DIR),
        f2='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.fq2.gz'.format(find_snps_dir=FIND_SNPS_DIR)
    output:
        dir='{remap_dir}/{{sample}}.{{genome}}'.format(remap_dir=REMAP_DIR),
        sam='{remap_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.Aligned.out.sam'.format(remap_dir=REMAP_DIR)
    threads: config['star_align_threads']
    params:
        overhang=config['read_length']-1,
        out_prefix = '{remap_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.'.format(remap_dir=REMAP_DIR)
    shell:
        "mkdir -p {output.dir}; "
        "STAR --readFilesIn {input.f1} {input.f2} --outFileNamePrefix {params.out_prefix} "
        "--genomeDir {input.genome_dir} --readFilesCommand zcat --runThreadN {threads} "
        "--genomeLoad NoSharedMemory --outFilterType BySJout "
        "--outSAMunmapped Within "
        "--outSAMattributes NH HI AS NM MD NM "
        "--twopassMode Basic "
        "--sjdbOverhang {params.overhang} "
        "--sjdbGTFfile {input.gencode} &>> {log}"


rule remap_add_read_groups:
    log: '{add_rg_dir}/{{sample}}.{{genome}}.log'.format(add_rg_dir=REMAP_RG_DIR)
    input:
        '{remap_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.Aligned.out.sam'.format(remap_dir=REMAP_DIR)
    output:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=REMAP_RG_DIR)
    params:
        tmp_dir=TMP_DIR
    run:
        command = "{picard} AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=id " \
                  "RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample TMP_DIR={params} &>> {log};".format(picard=config['picard_path'],
                  input=input, output=output, params=params.tmp_dir, log=log)

        print(command)
        shell(command)


rule remap_mark_duplicates:
    log: '{mark_dups_dir}/{{sample}}.{{genome}}.log'.format(mark_dups_dir=REMAP_MD_DIR),
    input:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=REMAP_RG_DIR)
    output:
        bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=REMAP_MD_DIR),
        metrics='{mark_dups_dir}/{{sample}}.{{genome}}.metrics'.format(mark_dups_dir=REMAP_MD_DIR)
    params:
        tmp_dir=TMP_DIR
    run:
        command = "{picard} MarkDuplicates I={input} O={output_bam} CREATE_INDEX=true " \
                  "VALIDATION_STRINGENCY=SILENT M={output_metrics} TMP_DIR={params} &>> {log}".format(picard=config['picard_path'],
                  input=input, output_bam=output.bam, output_metrics=output.metrics, params=params.tmp_dir, log=log)

        print(command)
        shell(command)


rule filter_remapped_reads:
    log: '{filt_remapped_reads_dir}/{{sample}}.{{genome}}.log'.format(filt_remapped_reads_dir=FILT_REMAPPED_READS_DIR)
    input:
        to_remap_bam='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.to.remap.bam'.format(find_snps_dir=FIND_SNPS_DIR),
        remap_bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=REMAP_MD_DIR),
    output:
        '{filt_remapped_reads_dir}/{{sample}}.{{genome}}.bam'.format(filt_remapped_reads_dir=FILT_REMAPPED_READS_DIR)
    conda:
        "envs/wasp.yaml"
    shell:
        'python scripts/filter_remapped_reads.py {input.to_remap_bam} {input.remap_bam} {output} &>> {log}'


rule samtools_merge:
    log: '{merged_dir}/{{sample}}.{{genome}}.log'.format(merged_dir=MERGED_DIR)
    input:
        filtered_bam='{filt_remapped_reads_dir}/{{sample}}.{{genome}}.bam'.format(filt_remapped_reads_dir=FILT_REMAPPED_READS_DIR),
        keep_bam='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.keep.bam'.format(find_snps_dir=FIND_SNPS_DIR)
    output:
        '{merged_dir}/{{sample}}.{{genome}}.bam'.format(merged_dir=MERGED_DIR)
    shell:
        'samtools merge {output} {input.filtered_bam} {input.keep_bam} &>> {log}'


rule samtools_sort:
    log: '{sorted_dir}/{{sample}}.{{genome}}.log'.format(sorted_dir=SORTED_DIR)
    input:
        '{merged_dir}/{{sample}}.{{genome}}.bam'.format(merged_dir=MERGED_DIR)
    output:
        '{sorted_dir}/{{sample}}.{{genome}}.bam'.format(sorted_dir=SORTED_DIR)
    params:
        tmp_dir=TMP_DIR,
        sample='{sample}'
    shell:
        'samtools sort -T {params.tmp_dir}/{params.sample} {input} 1> {output} 2> {log}'


rule samtools_index:
    log: '{sorted_dir}/{{sample}}.{{genome}}.bam.bai.log'.format(sorted_dir=SORTED_DIR)
    input:
        '{sorted_dir}/{{sample}}.{{genome}}.bam'.format(sorted_dir=SORTED_DIR)
    output:
        '{sorted_dir}/{{sample}}.{{genome}}.bam.bai'.format(sorted_dir=SORTED_DIR)
    shell:
        'samtools index {input} &>> {log}'


rule index_genome:
    input:
        "{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR)
    output:
        "{refgen_dir}/{{genome}}.fa.fai".format(refgen_dir=REFGEN_DIR)
    shell:
        "samtools faidx {input}"


rule create_genome_dictionary:
    log: "{refgen_dir}/{{genome}}.dict.log".format(refgen_dir=REFGEN_DIR)
    input:
        "{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR)
    output:
        "{refgen_dir}/{{genome}}.dict".format(refgen_dir=REFGEN_DIR)
    run:
        command = "{picard} CreateSequenceDictionary R={input} O={output} &>> {log}".format(picard=config['picard_path'],
                  input=input, output=output, log=log)

        print(command)
        shell(command)


rule ase_read_counter:
    log: '{read_counts_dir}/{{sample}}.{{genome}}.log'.format(read_counts_dir=READ_COUNTS_DIR)
    input:
        bam='{sorted_dir}/{{sample}}.{{genome}}.bam'.format(sorted_dir=SORTED_DIR),
        idx='{sorted_dir}/{{sample}}.{{genome}}.bam.bai'.format(sorted_dir=SORTED_DIR),
        ref='{refgen_dir}/{{genome}}.fa'.format(refgen_dir=REFGEN_DIR),
        refidx="{refgen_dir}/{{genome}}.fa.fai".format(refgen_dir=REFGEN_DIR),
        refdict="{refgen_dir}/{{genome}}.dict".format(refgen_dir=REFGEN_DIR),
        vcf='{vcf_dir}/{{sample}}.{{genome}}.vcf'.format(vcf_dir=VCF_DIR)
    output:
        '{read_counts_dir}/{{sample}}.{{genome}}.tsv'.format(read_counts_dir=READ_COUNTS_DIR)
    run:
        command = "{gatk} -T ASEReadCounter -I {input} -o {output} -sites {vcf} -R {refgen} -U ALLOW_N_CIGAR_READS -minDepth -1 " \
                  "--minBaseQuality -1 --minMappingQuality -1 -drf DuplicateRead &>> {log}".format(
                  gatk=config['gatk_path'], input=str(input.bam), vcf=input.vcf, output=str(output),
                  refgen=str(input.ref), log=log)

        print(command)
        shell(command)
