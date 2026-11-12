# Reproducible Germline Variant Calling and Benchmarking with Nextflow

**Nextflow DSL2 · Docker · BWA-MEM · GATK · SAMtools · Mosdepth · MultiQC · hap.py**

## 1. Problem

A bioinformatics pipeline can complete successfully without producing reliable biological results. This is particularly important in low-coverage whole-genome sequencing, where insufficient sequencing depth can lead to missed genetic variants despite successful read alignment.

This project investigated two questions:

1. **How does sequencing coverage affect the accuracy of germline SNP and INDEL detection?**
2. **Do additional quality-control and filtering steps change benchmark performance, and why?**

The objective was to develop reproducible germline variant-calling workflows, evaluate their performance against a trusted reference and understand the distinction between workflow quality assessment and variant-calling accuracy.

## 2. Approach

### Pipeline development

I developed and evaluated two containerised Nextflow DSL2 workflows using paired-end whole-genome sequencing data from NA12878 (HG001).

**Pipeline 1 - Original Nextflow workflow**

The original pipeline contained five processes:

1. **FastQC:** Sequencing read quality assessment
2. **BWA-MEM:** Alignment to the GRCh38 reference genome
3. **GATK MarkDuplicates:** Duplicate-read marking
4. **GATK HaplotypeCaller:** Germline SNP and INDEL calling
5. **MultiQC:** Consolidated quality-control reporting

**Pipeline 2 - Sarek-informed extended workflow**

After reviewing [nf-core/sarek](https://nf-co.re/sarek), I extended the original workflow with five additional processes:

- **SAMtools index:** BAM indexing
- **SAMtools flagstat:** Alignment and mapping statistics
- **Mosdepth:** Sequencing depth and breadth assessment
- **GATK variant filtering:** Quality-based filtering using a QUAL threshold of 30
- **hap.py:** Automated benchmarking against Genome in a Bottle (GIAB) reference calls

The extended pipeline contained ten processes and integrated additional quality assessment, filtering and benchmarking capabilities.

*nf-core/sarek was used as a workflow design reference. Both pipelines were independently developed and are not executions of the official nf-core/sarek pipeline.*

### Workflow overview

```
Paired-end FASTQ
      |
      +----> FastQC
      |
      v
 BWA-MEM Alignment
      |
      v
 GATK MarkDuplicates
      |
      +----> SAMtools Index / Flagstat
      |
      +----> Mosdepth
      |
      v
 GATK HaplotypeCaller
      |
      v
    Raw VCF
      |
      v
 QUAL-based Filtering
      |
      v
    hap.py
      |
      v
 Precision / Recall / F1
```

### Dataset provenance and preparation

Paired-end whole-genome sequencing data from **NA12878 (HG001), SRA accession SRR8454587**, were used to evaluate both workflows.

Three low-coverage datasets were prepared using **seqtk v1.4-r122**, applying the same random seed (`42`) to both FASTQ files.

| Target coverage | Sampling fraction | Random seed |
|---|---|---|
| 0.65× | 0.021326 | 42 |
| 1.3× | 0.042652 | 42 |
| 2.6× | 0.085304 | 42 |

The 0.65× subset was initially prepared for the original Nextflow workflow and subsequently reused in the Sarek-informed pipeline. The 1.3× and 2.6× subsets were prepared from the same sequencing accession.

Using the same source dataset and sampling strategy allowed variant-calling performance to be compared across increasing sequencing coverage levels.

### Experimental design and benchmarking

| Component | Description |
|---|---|
| Biological sample | NA12878 / HG001 |
| Sequencing accession | SRR8454587 |
| Reference genome | GRCh38 |
| Evaluation region | Chromosome 21 |
| Coverage conditions | 0.65×, 1.3×, 2.6× |
| Variant caller | GATK HaplotypeCaller |
| Benchmark truth set | GIAB HG001 GRCh38 v4.2.1 |
| Benchmarking tool | hap.py |
| Comparison engine | xcmp |
| Evaluation metrics | Precision, recall, F1 score |

Variant calls were compared against GIAB chromosome 21 truth variants using the corresponding confident-region BED file.

**Precision** measures the proportion of assessable positive calls matching the reference, **recall** measures the proportion of known truth variants detected, and **F1** represents the harmonic mean of precision and recall.

SNPs and INDELs were evaluated separately.

### Software and container versions

| Component | Version / container |
|---|---|
| Nextflow | 26.04.6 |
| seqtk | 1.4-r122 |
| FastQC | `biocontainers/fastqc:v0.11.9_cv8` |
| BWA-MEM | 0.7.17-r1188 |
| SAMtools (alignment environment) | 1.12 |
| SAMtools (QC/index container) | `quay.io/biocontainers/samtools:1.19--h50ea8bc_0` |
| Mosdepth | `quay.io/biocontainers/mosdepth:0.3.14--h87be163_2` |
| GATK | `broadinstitute/gatk:4.5.0.0` |
| MultiQC | `multiqc/multiqc:v1.21` |
| hap.py | `genonet/hap-py:0.3.15` |

Nextflow was used for workflow orchestration, Docker for containerised execution and Linux/WSL2 as the computational environment.

## 3. Results

### 3.1 Sequencing coverage assessment

Mosdepth analysis confirmed increasing chromosome 21 sequencing depth across the three datasets.

| Nominal coverage | Observed mean depth | Approximate breadth ≥1× |
|---|---:|---:|
| 0.65× | 0.60× | 37% |
| 1.3× | 1.21× | 57% |
| 2.6× | 2.40× | 74% |

At 2.6× coverage, approximately 26% of chromosome 21 reference positions remained uncovered.

However, part of this uncovered fraction reflects **assembly gaps (N regions)**, particularly on chromosome 21's short arm, which cannot be covered by aligned sequencing reads. Therefore, chromosome-wide breadth should not be interpreted as a direct measurement of insufficient sequencing depth at every uncovered position.

Coverage measurements and GIAB confident-region benchmarking also evaluate different sets of genomic positions.

### 3.2 SNP benchmarking

| Coverage | Precision | Recall | F1 score |
|---|---:|---:|---:|
| 0.65× | 67.61% | 5.90% | 10.86% |
| 1.3× | 74.20% | 16.20% | 26.60% |
| 2.6× | 83.34% | 38.41% | 52.58% |

**Interpretation**

Increasing coverage from 0.65× to 2.6× substantially improved SNP detection.

SNP F1 increased from **10.86% to 52.58%**, while recall increased from **5.90% to 38.41%**.

However, recall remained considerably lower than precision. At 2.6×, the pipeline detected **17,413 true-positive SNPs** but missed **27,927 known SNPs**.

This demonstrates that relatively high precision does not necessarily indicate comprehensive variant detection.

### 3.3 INDEL benchmarking

| Coverage | Precision | Recall | F1 score |
|---|---:|---:|---:|
| 0.65× | 72.07% | 3.95% | 7.49% |
| 1.3× | 76.77% | 11.21% | 19.57% |
| 2.6× | 82.12% | 27.50% | 41.21% |

**Interpretation**

INDEL detection also improved as sequencing coverage increased.

INDEL F1 increased from **7.49% to 41.21%**, while recall increased from **3.95% to 27.50%**.

At 2.6×, the pipeline recovered **1,886 true-positive INDELs** but missed **4,971 known INDELs**.

INDEL recall remained lower than SNP recall across all three coverage levels, highlighting the challenges of detecting insertions and deletions under the tested low-coverage conditions.

### 3.4 Original vs Sarek-informed pipeline

Both workflows produced matching benchmark metrics at all three sequencing coverage levels.

| Evaluation component | Original pipeline | Sarek-informed pipeline |
|---|---|---|
| Read quality assessment | FastQC | FastQC |
| Read alignment | BWA-MEM | BWA-MEM |
| Duplicate marking | GATK | GATK |
| Germline variant calling | GATK HaplotypeCaller | GATK HaplotypeCaller |
| Alignment statistics | Not integrated | SAMtools flagstat |
| Coverage assessment | Not integrated | Mosdepth |
| Variant filtering | Not integrated | QUAL-based filtering |
| Automated benchmarking | Not integrated | GIAB / hap.py |
| Benchmark accuracy | Baseline | Matching baseline metrics |

The additional QC processes measure workflow outputs rather than modifying variant calls. Therefore, improvements in variant-calling accuracy were not expected from these processes alone.

The only additional variant-processing step capable of changing the evaluated calls was the QUAL-based filter.

#### Investigating the filtering result

The extended pipeline applied a variant quality threshold of **QUAL ≥30**.

To investigate why the benchmark metrics remained unchanged, I examined the raw HaplotypeCaller VCF using:

```bash
bcftools view -H -i 'QUAL<30' \
  results/coverage_2.6x/variants/NA12878.vcf.gz | wc -l
```

**Observed result:**

```text
0
```

The raw 2.6× VCF contained no variants with QUAL below 30.

This is consistent with GATK HaplotypeCaller's standard minimum calling-confidence threshold of 30, although the effective threshold also depends on the caller configuration.

At 2.6× coverage, the raw HaplotypeCaller VCF contained no variants with QUAL below 30. Therefore, the additional QUAL-based filtering criterion did not exclude any variants at this coverage level. The identical benchmarking results between the original and extended pipelines indicate that the added filtering step did not improve measured variant-calling performance under the tested conditions.

**Key finding:** A filtering criterion that is already satisfied by the caller's output provides no additional discrimination under the tested conditions. More informative filtering experiments would require appropriately selected criteria, such as quality-by-depth (QD), strand-bias (FS) or mapping-quality (MQ) annotations.

The Sarek-informed workflow therefore improved **quality assessment, coverage reporting, traceability and automated benchmarking**, rather than measured variant-calling accuracy.

### 3.5 Overall findings

**Finding 1 - Sequencing coverage strongly influenced variant detection.**

Increasing sequencing depth improved precision, recall and F1 scores for both SNPs and INDELs.

**Finding 2 - Recall remained the primary performance limitation.**

Even at 2.6×, many known variants were not recovered, demonstrating the limitations of low-coverage germline variant calling.

**Finding 3 - Filtering must provide information beyond existing calling criteria.**

The raw 2.6× VCF contained no variants below QUAL 30, explaining why applying that threshold did not improve benchmark performance.

**Finding 4 - Successful workflow execution does not guarantee reliable biological results.**

High read-mapping rates and completed pipeline processes cannot substitute for independent variant-call benchmarking.

## 4. Limitations

**Low sequencing coverage:** All three datasets were below conventional high-coverage germline sequencing depths, limiting sensitivity.

**Single reference sample:** The evaluation used only NA12878 (HG001), so performance across additional biological samples was not assessed.

**Chromosome-specific evaluation:** Benchmarking focused on chromosome 21 rather than genome-wide variant-calling performance.

**Single variant caller:** GATK HaplotypeCaller was used throughout. Alternative callers were not evaluated.

**Limited filtering experiment:** The QUAL ≥30 filter did not provide additional discrimination in the inspected 2.6× dataset. More informative filtering strategies were not tested.

**Benchmark confidence regions:** Calls outside confidently evaluated regions cannot be fully classified, limiting interpretation of precision across all emitted variants.

**No high-coverage comparison:** The project did not evaluate performance at substantially higher sequencing depths.

Future work could investigate higher-coverage sequencing, alternative variant callers or additional filtering criteria while controlling experimental variables.

## 5. Conclusion

This project demonstrated that **sequencing coverage substantially influences germline variant-calling performance**, with increasing depth improving SNP and INDEL precision, recall and F1 scores.

However, recall remained limited across all three low-coverage datasets, highlighting the importance of independent benchmarking rather than relying solely on successful workflow execution.

Extending the original Nextflow pipeline with alignment QC, coverage assessment, variant filtering and automated benchmarking provided a more comprehensive framework for evaluating analytical results.

The investigation also demonstrated that **additional filtering does not necessarily improve accuracy when the applied threshold is already satisfied by the variant caller's output**.

Overall, this project demonstrates practical skills in **Nextflow DSL2, Docker, Linux, reproducible workflow development, sequence alignment, GATK variant calling, quality assessment, GIAB benchmarking and critical interpretation of computational results**.

---

## Reproducibility

### Requirements

- Linux or WSL2
- Nextflow and compatible Java runtime
- Docker
- Paired-end FASTQ sequencing files
- GRCh38 reference genome and required indexes
- GIAB truth VCF and confident-region BED for benchmarking

### Dataset preparation

The three sequencing subsets were prepared from SRR8454587 using `seqtk sample` with seed 42.

Example for the 2.6× dataset:

```bash
seqtk sample -s42 SRR8454587_1.fastq.gz 0.085304 \
  | gzip > NA12878_R1.fastq.gz

seqtk sample -s42 SRR8454587_2.fastq.gz 0.085304 \
  | gzip > NA12878_R2.fastq.gz
```

The corresponding sampling fractions for 0.65× and 1.3× were `0.021326` and `0.042652`, respectively.

### Running the workflow

Configure the appropriate input FASTQ files, reference genome, interval and output directory in `nextflow.config`.

Run:

```bash
nextflow run main.nf
```

To enable GIAB benchmarking, set:

```groovy
params.benchmark = true
```

and configure the matching truth VCF and confident-region BED file.

To resume an interrupted analysis:

```bash
nextflow run main.nf -resume
```

Nextflow reuses eligible cached tasks, reducing unnecessary recomputation.

Large sequencing datasets, reference files and intermediate analysis outputs are excluded from version control.

### Related project

[Original Nextflow Germline Variant-Calling Pipeline](https://github.com/Adeola29-star/nextflow-based-variant-calling-pipeline)

### References

- [nf-core/sarek](https://nf-co.re/sarek)
- [Genome in a Bottle — NIST](https://www.nist.gov/programs-projects/genome-bottle)
- [hap.py Benchmarking Toolkit](https://github.com/Illumina/hap.py)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
