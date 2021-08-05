# Viral Variant Detection
Detection, annotation, and quantification of viral variants in short-read Illumina RNA-seq datasets.

## Requirements
Datasets must be pre-aligned to the indexed viral genome. Please choose an appropriate aligner. Bowtie2 and STAR are recommended for accuracy and sensitivity. BBMap is a flexible tool for alignment and data parsing. See https://www.ecseq.com/support/ngs/best-RNA-seq-aligner-comparison-of-mapping-tools for more info.

Example alignment workflow:

```
cd target directory/
input="./samples.txt"
while IFS= read -r line
do
	java -jar /home/denison-thelio/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 ${line}_R1.fastq.gz ${line}_R2.fastq.gz ${line}_R1_paired.fastq ${line}_R1_unpaired.fastq ${line}_R2_paired.fastq ${line}_R2_unpaired.fastq ILLUMINACLIP:/home/denison-thelio/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	bowtie2 -p 32 -q -x SARSCoV2 -1 ${line}_R1_paired.fastq -2 ${line}_R2_paired.fastq -U ${line}_R1_unpaired.fastq,${line}_R2_unpaired.fastq -S ${line}_bowtie2.sam
	samtools view -b -@ 32 ${line}_bowtie2.sam > ${line}_bowtie2.bam
	samtools sort -@ 32 -o ${line}_bowtie2.sort.bam ${line}_bowtie2.bam
	samtools index -@ 32 -b ${line}_bowtie2.sort.bam ${line}_bowtie2.sort.bam.bai
	/home/denison-thelio/bbmap/pileup.sh in=${line}_bowtie2.sam basecov=${line}_bowtie2_coverage.txt delcoverage=f 32bit=t -Xmx64g
	lofreq call-parallel --pp-threads 32 -f SARSCoV2_virema.fasta -d 100000 -o ${line}.vcf ${line}_bowtie2.sort.bam
done < "$input"
```

Python scripts require packages `pandas` `os` and `fnmatch`

Script was written under Python version 3.9.

**Required input files:**
1) VCF (v4.0) named {sample}.vcf
```
##fileformat=VCFv4.0
##fileDate=20201030
##source=lofreq call -f SARSCoV2_virema.fasta -d 100000 -o SARS2-05-1.vcf SARS2-05-1_bowtie2.sort.bam 
##reference=SARSCoV2_virema.fasta
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias at this position">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description="Indicates that the variant is a consensus variant (as opposed to a low frequency variant).">
##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer length to the right of report indel position">
##FILTER=<ID=min_dp_10,Description="Minimum Coverage 10">
##FILTER=<ID=sb_fdr,Description="Strand-Bias Multiple Testing Correction: fdr corr. pvalue > 0.001000">
##FILTER=<ID=min_snvqual_68,Description="Minimum SNV Quality (Phred) 68">
##FILTER=<ID=min_indelqual_20,Description="Minimum Indel Quality (Phred) 20">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
MT020881.1	821	.	G	A	108	PASS	DP=955;AF=0.014660;SB=0;DP4=454,486,7,7
MT020881.1	884	.	C	T	74	PASS	DP=743;AF=0.012113;SB=13;DP4=392,342,8,1
MT020881.1	3927	.	C	T	68	PASS	DP=468;AF=0.017094;SB=5;DP4=278,181,3,5
MT020881.1	4162	.	G	A	92	PASS	DP=599;AF=0.020033;SB=1;DP4=224,363,5,7
MT020881.1	5150	.	G	A	69	PASS	DP=603;AF=0.016584;SB=6;DP4=279,310,7,3
MT020881.1	5457	.	C	T	459	PASS	DP=412;AF=0.070388;SB=0;DP4=203,179,16,13
MT020881.1	11083	.	G	T	171	PASS	DP=417;AF=0.035971;SB=25;DP4=186,212,13,2
MT020881.1	13422	.	T	C	86	PASS	DP=585;AF=0.023932;SB=32;DP4=239,328,0,14
MT020881.1	14637	.	T	C	68	PASS	DP=658;AF=0.012158;SB=5;DP4=280,370,5,3
MT020881.1	14679	.	T	C	101	PASS	DP=544;AF=0.027574;SB=26;DP4=241,285,1,14
MT020881.1	17142	.	T	C	84	PASS	DP=496;AF=0.022177;SB=23;DP4=341,143,3,8
```
2) BBMap coverage file named {sample}_coverage.txt
Currently, only coverage input files with column headers accepted. A future version will include coverage files that do not have headers for user flexibility.
```
#RefName	Pos	Coverage
MT020881.1	0	2
MT020881.1	1	2
MT020881.1	2	5
MT020881.1	3	7
MT020881.1	4	7
MT020881.1	5	7
MT020881.1	6	8
MT020881.1	7	9
MT020881.1	8	13
MT020881.1	9	13
MT020881.1	10	15
MT020881.1	11	28
MT020881.1	12	63
MT020881.1	13	85
MT020881.1	14	113
MT020881.1	15	273
MT020881.1	16	573
MT020881.1	17	800
MT020881.1	18	1611
MT020881.1	19	1886
MT020881.1	20	2574
MT020881.1	21	3100
```

3) Text file of samples to run in an experiment. No underscores, or spaces. Must match {sample} name of VCF and coverage files.
```
SARS2-DMSO-1
SARS2-DMSO-2
SARS2-DMSO-3
SARS2-0125-1
SARS2-0125-2
SARS2-0125-3
SARS2-05-1
SARS2-05-2
SARS2-05-3
```

## Usage

LoFreq_VCF_analysis_CLI.py is run via command line. There is a "debug" version that can be run in PyCharm and is used for updating functionalities.

```
python3 ./LoFreq_VCF_analysis_CLI.py Sample_List Virus Working_Directory Experiment [options]

Sample_List = path to .txt file listing samples in an experiment (1 sample name per line)
Virus = "MHV", "MERS", "SARS2", or "other"
Working_Directory = path to directory where VCF and coverage files are found
Experiment = name of experiment for outout files

Options:
--freq Frequency cutoff (0 to 1) for filtering variants
--tag File naming tag to match filter applied by --freq flag
```
Future versions will have the ability to input feature files (GTF) to provide annotations for genomes other than MHV (AY910861.1), MERS-CoV (JX869059.2), and SARS-CoV-2 (MT020881.1). If you select "other" for virus field, gene in output file will be listed as "unknown"

## Output files
For each sample in an experiment, a tab-delineated {sample}_variants.txt file and a {experiment}_summary.txt file are created.
1) Variants file
```genome	position	reference	variant	frequency	raw_depth	ref_f_count	ref_r_count	variant_f_count	variant_r_count	variant_total	variant_type	gene
MT020881.1	4162	G	A	0.020033000000000002	599	224	363	5	7	12	transition	nsp3
MT020881.1	5150	G	A	0.016584	603	279	310	7	3	10	transition	nsp3
MT020881.1	5457	C	T	0.070388	412	203	179	16	13	29	transition	nsp3
MT020881.1	11083	G	T	0.035970999999999996	417	186	212	13	2	15	transversion	nsp6
MT020881.1	13422	T	C	0.023932	585	239	328	0	14	14	transition	nsp10
MT020881.1	14637	T	C	0.012158	658	280	370	5	3	8	transition	nsp12
MT020881.1	14679	T	C	0.027574	544	241	285	1	14	15	transition	nsp12
MT020881.1	17142	T	C	0.022177000000000002	496	341	143	3	8	11	transition	nsp13
MT020881.1	19275	T	C	0.01566	894	447	430	4	10	14	transition	nsp14
MT020881.1	19633	G	A	0.007402	1486	757	715	8	3	11	transition	nsp15
MT020881.1	19896	T	C	0.011918000000000002	1762	978	756	6	15	21	transition	nsp15
MT020881.1	21747	T	G	0.017797999999999998	2416	1322	1049	14	29	43	transversion	S protein
MT020881.1	21784	T	A	0.021172	2645	1290	1277	29	27	56	transversion	S protein
MT020881.1	22009	C	T	0.006667	3300	1489	1785	8	14	22	transition	S protein
MT020881.1	22020	T	C	0.007547	3180	1402	1749	12	12	24	transition	S protein
MT020881.1	22114	T	C	0.022702	1806	1023	694	10	31	41	transition	S protein
MT020881.1	22350	C	T	0.009745	2976	1383	1560	15	14	29	transition	S protein
MT020881.1	23242	G	A	0.0063880000000000004	4070	2211	1825	12	14	26	transition	S protein
MT020881.1	23403	A	G	0.010151	3448	1523	1881	17	18	35	transition	S protein
```
Future releases will contain the nonsynonymous amino acid changes based on annotated genes.

2) Summary file
```
sample	unique_variants	variant_nts	total_nts	transition_nts	transversion_nts	AtoG_nts	GtoA_nts	CtoT_nts	TtoC_nts	AtoT_nts	TtoA_nts	AtoC_nts	CtoA_nts	CtoG_nts	GtoC_nts	GtoT_nts	TtoG_nts	mutation_freq	transition_freq	transversion_freq	AtoG_freq	GtoA_freq	CtoT_freq	TtoC_freq	AtoT_freq	TtoA_freq	AtoC_freq	CtoA_freq	CtoG_freq	GtoC_freq	GtoT_freq	TtoG_freq
SARS2-DMSO-1	65	34588	1572038080	15317	19271	2192	805	8836	3484	7314	6304	1475	46	0	0	108	4024	2.20E-05	9.74E-06	1.23E-05	1.39E-06	5.12E-07	5.62E-06	2.22E-06	4.65E-06	4.01E-06	9.38E-07	2.93E-08	0	0	6.87E-08	2.56E-06
SARS2-DMSO-2	62	35948	1950550128	14779	21169	1817	280	9859	2823	7999	6491	2330	67	0	0	131	4151	1.84E-05	7.58E-06	1.09E-05	9.32E-07	1.44E-07	5.05E-06	1.45E-06	4.10E-06	3.33E-06	1.19E-06	3.43E-08	0	0	6.72E-08	2.13E-06
SARS2-DMSO-3	56	27333	2044217597	12982	14351	1788	270	8649	2275	250	6548	1594	49	0	0	530	5380	1.34E-05	6.35E-06	7.02E-06	8.75E-07	1.32E-07	4.23E-06	1.11E-06	1.22E-07	3.20E-06	7.80E-07	2.40E-08	0	0	2.59E-07	2.63E-06
```
