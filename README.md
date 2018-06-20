# UCLA UDN RNA-seq Pipeline
The UDN RNA-seq pipeline is loosely adapted from those used by the GTEx Consortium and TOPMed. Once finalized, this pipeline will be converted into a docker image for deployment on DNAnexus.
## References and annotation
### Genome 
WGS alignments for the UCLA UDN are based on the hs37d5.fa (GRCh37/hg19 + decoy sequences).
For RNA-seq, we remove the decoy sequences from this reference, which makes it essentially identical to GRCh37 without alternate/haplotype contigs (aka primary assembly).
```bash
head -n 51696830 human_g1k_hs37d5.fasta > hs37d5_nodecoy.fasta
```
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
## Pipeline execution
### 1. Running STAR alignment
This is the excerpt from the  BASH script that runs the STAR aligner. The only input that needs to be set is the ```$SAMPLE``` variable.

I have decided on the following changes:
* Alignment settings have been modified slightly to match GTEx/TOPMed. 
* Chimeric alignments are output to a separate splicejunction file (will be useful for contrasting with SVs). 
* ```ID:rg1``` and ```SM:${SAMPLE}``` are now input at this stage (removes downstream processing step). I will just use ```rg1``` for all readgroup IDs, I think putting a full sample identifier here makes the BAM file unecessarily big.
* TO DOs:
	1. Test if duplicate marking using STAR is suitable for RNASeQC (removes an additional downstream step).
	2. Generate stranded BIGWIG files (useful for future genome browser / database application).
	3. Parse out RG information from fastq file??
```bash
SAMPLE="UDN369194-Fibroblast"

star_index=$HOME/udn_rnaseq/references/star_oh74
fastq1=$HOME/udn_rnaseq/reads/${SAMPLE}_R1.fastq.gz
fastq2=$HOME/udn_rnaseq/reads/${SAMPLE}_R2.fastq.gz
sample_id=${SAMPLE}
prefix="$HOME/udn_rnaseq/output/star/${SAMPLE}."

$HOME/bin/STAR \
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
    --quantTranscriptomeBAMcompression -1 \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMcompression -1 \
    --outSAMunmapped Within \
    --genomeLoad NoSharedMemory \
    --chimSegmentMin 15 \
    --chimJunctionOverhangMin 15 \
    --chimOutType Junctions WithinBAM SoftClip \
    --chimMainSegmentMultNmax 1 \
    --outSAMattributes NH HI AS nM NM ch \
    --outSAMattrRGline ID:rg1 PL:Illumina LB:${sample_id} SM:${sample_id} CN:UNGC 
```
### 2. RSEM quantification
Below is the excerpt that performs RSEM quantification.
We are counting reads in a stranded fashion (```--forward-prob 0```). We may wish to set this parameter to ```0.5``` to enable a more equitable comparison with unstranded RNA-seq (e.g. GTEx).
```bash
 ~/RSEM-1.3.0/rsem-calculate-expression \
    --num-threads 8 \
    --fragment-length-max 1000 \
    --no-bam-output \
    --paired-end \
    --estimate-rspd \
    --forward-prob 0 \
    --bam output/star/${SAMPLE}.Aligned.toTranscriptome.out.bam \
    ~/udn_rnaseq/references/rsem/rsem_reference \
    ~/udn_rnaseq/output/rsem/${SAMPLE}
```
### 3. Output generation
Mark duplicates and creates final BAM file and index for upload to UDN server. Compression level is set to maximum and this actually takes FOREVER; will save on storage cost but result in slower processing downstream (but we probably won't touch the file much after this).
```bash
java -Xmx24g -jar ~/picard.jar MarkDuplicates \
	I=${SAMPLE}.Aligned.sortedByCoord.out.bam \
    O=${SAMPLE}.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M=${SAMPLE}.markdup.metrics \
    COMPRESSION_LEVEL=9
```
### 4. RNASeQC

Still checking parameters for this part...

```bash
java -Xmx24g -jar RNA-SeQC.jar \
	-n 1000 \
    -s ${SAMPLE},${SAMPLE}.bam,${SAMPLE} \
    -t gencode.v19.annotation.hs37d5.collapsed.gtf \
    -r hs37d5_nodecoy.fasta \
    -strictMode \
    -gatkFlags ---allow_potentially_misencoded_quality_scores
```
