# UCLA UDN RNA-seq pipeline instructions for DNAnexus

This is our UDN RNA-seq pipeline (loosely adapted from those used by the GTEx Consortium and TOPMed).

There are three sections:
1. Software used and where to retrieve them (hyperlinks).
2. Steps to generate references/annotation.
3. Pipeline steps and execution.


## Section 1. Software
* [Picard tools v2.9.0](https://github.com/broadinstitute/picard/releases/download/2.9.0/picard.jar) - Requires Java 1.8
* [STAR v2.6.0a](https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz) - Used pre-compiled binaries in bin/Linux_x86_64
* [RNASeQC v1.1.9](https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz) - Requires Java 1.7
* [RSEM v1.3.0](https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz) 

## Section 2. References and annotation

### Genome 
WGS alignments for the UCLA UDN are based on the hs37d5.fa (GRCh37/hg19 + decoy sequences).
For RNA-seq, we will use GRCh37 without the decoy sequence or alternate/haplotype contigs (primary assembly).

*[GRCh37 fasta](http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz)
*[GRCh37 index](http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai)

### Gene annotation
We use GENCODE v19, which is the last annotation build native to GRCh37. We need to change the chromosome naming convention to match our genomic reference.
```bash
wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
zcat gencode.v19.annotation.gtf.gz | sed 's/^chrM/chrMT/1;s/^chr//1' > gencode.v19.annotation.hs37d5.gtf
```
### STAR index
The majority of our sequencing was performed using 75bp reads; we generate a STAR index using a 74bp overhang.
```bash
STAR \
--runMode genomeGenerate \
--genomeDir star_oh74/ \
--genomeFastaFiles hs37d5_nodecoy.fasta \
--sjdbGTFfile gencode.v19.annotation.hs37d5.gtf \
--sjdbOverhang 74 \
--runThreadN 8
```
### RSEM reference
```bash
rsem-prepare-reference \
    --num-threads 8 \
    --gtf gencode.v19.annotation.hs37d5.gtf \
    hs37d5_nodecoy.fasta \
    rsem/rsem_reference
```
### RNASeQC GTF reference

Any GTF file will work; we use here a flattened gene model version and exclude certain intron rentention transcripts on the gTEx blacklist.

## Section 3. Pipeline execution

### 1. STAR alignment and generation of UDN BAM file for upload

This will run the STAR aligner, generate a coordinate sorted, duplicate marked, compressed BAM file for upload to the UDN gateway, and produce other files needed for quantification.

Runs STAR and Picard tools (I guess this could be strung together from the existing apps on DNAnexus; I did not test with the versions that are available but it should work).

This step requires the following input variables to be set (with example values):

```bash
sample_id="UDN369194-Fibroblast"		#sample name
star_index="/references/star_oh74"		#location of the STAR index
fastq1="/reads/UDN369194-Fibroblast_R1.fastq.gz"		#R1 fastq file
fastq2="/reads/UDN369194-Fibroblast_R2.fastq.gz"		#R2 fastq file
outdir="/output/star"		#prefix output DIRECTORY
```
BASH script:

```bash

## define output file prefix
prefix="${outdir}/${sample_id}."

## run STAR alignment
/opt/star2.6.0a/bin/Linux_x86_64/STAR \
    --runMode alignReads \
    --runThreadN 6 \
    --genomeDir ${star_index} \
    --twopassMode Basic \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFilterType BySJout \
    --outFilterScoreMinOverLread 0.33 \
    --outFilterMatchNminOverLread 0.33 \
    --limitSjdbInsertNsj 1200000 \
    --readFilesIn ${fastq1} ${fastq2} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${prefix} \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs None \
    --alignSoftClipAtReferenceEnds Yes \
    --quantMode TranscriptomeSAM GeneCounts \
    --quantTranscriptomeBAMcompression 9 \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMcompression 1 \
    --outSAMunmapped Within \
    --genomeLoad NoSharedMemory \
    --chimSegmentMin 15 \
    --chimJunctionOverhangMin 15 \
    --chimOutType Junctions WithinBAM SoftClip \
    --chimMainSegmentMultNmax 1 \
    --outSAMattributes NH HI AS nM NM ch \
    --outSAMattrRGline ID:rg1 PL:Illumina LB:${sample_id} SM:${sample_id}

## Picard mark duplicates

java -jar /udn_rnaseq/bin/picard.jar MarkDuplicates \
	I=${prefix}Aligned.sortedByCoord.out.bam \
    O=${prefix}bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M=${prefix}rmdup.metrics

## Compress STAR output files
gzip ${outdir}/*.junction ${outdir}/*.tab

```
Output files to keep (based on supplied inputs):

```bash
/output/star/UDN369194-Fibroblast.bam		#coordinate-sorted rmdup BAM for UDN upload
/output/star/UDN369194-Fibroblast.bai		#BAM index
/output/star/UDN369194-Fibroblast.ReadsPerGene.out.tab.gz		# Genecount table
/output/star/UDN369194-Fibroblast.Chimeric.out.junction.gz		# CHimera table
/output/star/UDN369194-Fibroblast.SJ.out.tab.gz		#Splice table
/output/star/UDN369194-Fibroblast.Aligned.toTranscriptome.out.bam		#Transcriptome BAM

```
### 2. RSEM quantification

Perform RSEM quantification (Requires Java 1.7).

Required inputs:
```bash
transcriptbam="/output/star/UDN369194-Fibroblast.Aligned.toTranscriptome.out.bam		#input transcriptome bam (from STAR)"
prob=0		#strand-specific setting (default to 0 for TruSeq stranded RNA)
rsem_reference="/references/rsem/rsem_reference"		#prefix of rsem reference
prefix="/output/rsem/UDN369194-Fibroblast"		#prefix of output files
```
BASH code:
```bash

## Run quantification
/udn_rnaseq/opt/RSEM-1.3.0/rsem-calculate-expression \
    --num-threads 4 \
    --fragment-length-max 1000 \
    --no-bam-output \
    --paired-end \
    --estimate-rspd \
    --forward-prob 0 \
    --calc-ci
    --bam output/star/${SAMPLE}.Aligned.toTranscriptome.out.bam \
    ${rsem_reference} \
    ${prefix}
    
## Compress outputs
gzip ${prefix}.*.results

```
Output files to keep:
```bash
/output/rsem/UDN369194-Fibroblast.isoforms.results.gz		#isoform quant
/output/rsem/UDN369194-Fibroblast.genes.results.gz		#gene quant
```


### 3. RNASeQC

Inputs:
```bash
sample_id="UDN369194-Fibroblast"
input_bam="/output/star/UDN369194-Fibroblast.bam"
reference_fasta="/reference/hs37d5_nodecoy.fasta"
gtf_file="/reference/hs37d5_genes.gtf"
outdir="/output/rnaseqc"
```

BASH:
```bash
## run program
java -Xmx6g -jar /opt/RNA-SeQC.jar \
	-n 1000 \
    -s ${sample_id},${input_bam},${sample_id} \
    -t ${gtf_file} \
    -r ${reference_fasta} \
    -strictMode \
    -noDoC \
    -gatkFlags --allow_potentially_misencoded_quality_scores \
    -o $outdir
    
## clean output
mv  ${outdir}/genes.rpkm.gct ${outdir}/${sample_id}/${sample_id}.genes.rpkm.gct
cd ${outdir} && tar -czvf ../${sample_id}.rnaseqc.tar.gz *
```

Output file to keep:
```bash
/output/rnaseqc/UDN369194-Fibroblast.rnaseqc.tar.gz
```
Sample plots from RNASeQC output (note that these are test cases; many are known problematic samples):

![alt tag](https://github.com/furbelows/udn_rnaseq/blob/master/qc1.png "qc1")
![alt tag](https://github.com/furbelows/udn_rnaseq/blob/master/qc2.png "qc2")

## Additional thoughts
* We should add proportion of globin RNA as a QC metric for blood samples.
* Include raw counts for genes and exons (featurecounts or htseq-count).
* WASP filtering?
* Docker image?
