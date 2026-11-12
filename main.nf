nextflow.enable.dsl = 2

process FASTQC {
    container 'biocontainers/fastqc:v0.11.9_cv8'

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path("*.{html,zip}")

    script:
    """
    fastqc ${read1} ${read2}
    """
}

process BWA_ALIGN {
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0'

    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path reference
    path reference_index, stageAs: '*'
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
	"""
	bwa mem -t 2 -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \
    	${reference} ${read1} ${read2} | \
	samtools sort -@ 2 -m 512M -o ${sample_id}.sorted.bam -

	"""
    
}


process SAMTOOLS_FLAGSTAT {
    container 'quay.io/biocontainers/samtools:1.19--h50ea8bc_0'

    publishDir "${params.outdir}/alignment_qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}.flagstat.txt")

    script:
    """
    samtools flagstat ${bam} > ${sample_id}.flagstat.txt
    """
}

process MOSDEPTH {
    container 'quay.io/biocontainers/mosdepth:0.3.14--h87be163_2'

    publishDir "${params.outdir}/depth_qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}.mosdepth.summary.txt")
    path("${sample_id}.mosdepth.global.dist.txt")

    script:
    """
    mosdepth -t 2 -c ${params.interval} ${sample_id} ${bam}
    """
}

process INDEX_BAM {

    tag "${sample_id}"

    container 'quay.io/biocontainers/samtools:1.19--h50ea8bc_0'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai")

    script:
    """
    samtools index ${bam}
    """
}

process MARK_DUPLICATES {
    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/dedup", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.metrics.txt")

    script:
    """
    gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample_id}.dedup.bam \
        -M ${sample_id}.dedup.metrics.txt
    """
}

process HAPLOTYPE_CALLER {
    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(metrics)
    path reference
    path reference_index

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    """
    samtools index ${bam}
    gatk HaplotypeCaller \
        -R ${reference} \
        -I ${bam} \
        -L ${params.interval} \
        -O ${sample_id}.vcf.gz
    """
}

process FILTER_VARIANTS {
    tag "${sample_id}"

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/filtered_variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz")

    script:
    """
    gatk IndexFeatureFile \
        -I ${vcf}

    gatk VariantFiltration \
        -V ${vcf} \
        -O ${sample_id}.filtered.vcf.gz \
        --filter-name "LowQUAL" \
        --filter-expression "QUAL < 30.0"
    """

}

process MULTIQC {
    container 'multiqc/multiqc:v1.21'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc .
    """
}


process BENCHMARK_HAPPY {
    tag "${sample_id}"

    container 'genonet/hap-py:0.3.15'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}/benchmark", mode: 'copy'

    cpus 2

    input:
    tuple val(sample_id), path(vcf)
    path truth_vcf
    path truth_vcf_index
    path truth_bed
    path reference
    path reference_index


    output:
    path("happy.summary.csv")
    path("happy.runinfo.json")

    script:
    """
    /opt/hap.py/bin/hap.py \
        ${truth_vcf} \
        ${vcf} \
        -f ${truth_bed} \
        -r ${reference} \
        -o happy \
        --engine xcmp \
        --pass-only \
        --threads=${task.cpus}
    """
}

workflow {

    reads_ch = channel.of(
        [params.sample_id, file(params.reads_r1), file(params.reads_r2)]
    )

    reference_ch = file(params.reference)

    reference_index_ch = channel.fromPath([
        "${params.reference}.*",
        params.reference_dict
    ]).collect()

    // 1. Initial sequencing quality control
    FASTQC(reads_ch)

    // 2. Read alignment
    BWA_ALIGN(reads_ch, reference_ch, reference_index_ch)

    // 3. Duplicate marking
    MARK_DUPLICATES(BWA_ALIGN.out)

    // 4. Prepare duplicate-marked BAM for quality assessment
    dedup_bam_ch = MARK_DUPLICATES.out.map { sample_id, bam, metrics ->
        tuple(sample_id, bam)
        }

    // 5. BAM indexing and alignment quality assessment
    INDEX_BAM(dedup_bam_ch)

    SAMTOOLS_FLAGSTAT(INDEX_BAM.out)

    // 6. Sequencing depth and coverage assessment
    MOSDEPTH(INDEX_BAM.out)
    
    // 7. Germline variant calling
    HAPLOTYPE_CALLER(
        MARK_DUPLICATES.out,
        reference_ch,
        reference_index_ch
    )

    // 8. Variant quality filtering
    FILTER_VARIANTS(HAPLOTYPE_CALLER.out)

    // 9. Optional GIAB benchmarking
    if (params.benchmark) {

        if (!params.truth_vcf || !params.truth_bed) {
            error "Benchmarking requires --truth_vcf and --truth_bed"
        }

        def truthVcf = file(params.truth_vcf)
        def truthBed = file(params.truth_bed)

        if (!truthVcf.exists() || !truthBed.exists()) {
            error "Benchmarking truth VCF or BED file does not exist"
        }

        BENCHMARK_HAPPY(
            FILTER_VARIANTS.out,
            truthVcf,
            file("${params.truth_vcf}.tbi"),
            truthBed,
            reference_ch,
            file("${params.reference}.fai")
        )



    }

    // 10. Consolidated quality control reporting
    MULTIQC(
        FASTQC.out
            .mix(MARK_DUPLICATES.out.map { _sample_id, _bam, metrics -> metrics })
            .mix(SAMTOOLS_FLAGSTAT.out)
            .mix(MOSDEPTH.out[0])
            .mix(MOSDEPTH.out[1])
            .collect()
    )
}


