# Bulk_RNAseq_Analysis
This repository contains scripts, workflows, and documentation for performing bulk RNA sequencing (RNA-seq) analysis. All steps were implemented hands-on, including cloud-based setup, raw data processing, alignment, QC, and downstream analysis, demonstrating practical proficiency in RNA-Seq workflows.

# Table of Contents
1. [Dataset Overview](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#introduction) 
2. [Google Cloud Setup](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#setting_up_google_cloud_console) 
3. [Installation of Tools](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#setting_up_the_terminal_and_installation_of_tools ) 
4. [Data Downloading](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#1-dataset_downloading) 
5. [Organizing Data](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#2-organizing_dataset_concatenate__renaming) 
6. [Quality Check](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#3-quality_check) 
7. [Trimming](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#4-trimming_if_needed) 
8. [Alignment & Indexing](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#5-indexing_and_alignment) 
9. [BAM File QC](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#6-quality_check_of_bam_files) 
10. [Gene Quantification](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#7-gene_expression_quantification_using_featurecounts) 
11. [Count Matrix Generation](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#8-gene-level_count_matrix_generation_python__r) 
12. [Downstream Analysis](https://github.com/Bidya122/Bulk_RNAseq_Analysis/blob/main/README.md#8-gene-level_count_matrix_generation_python__r)

[Find Scripts](https://github.com/Bidya122/Bulk_RNAseq_Analysis/tree/main/Data)
[Find Output files](https://github.com/Bidya122/Bulk_RNAseq_Analysis/tree/main/Output)   


# Introduction  
This repository contains scripts, workflows, and documentation for performing bulk RNA sequencing (RNA-seq) analysis. All steps were implemented hands-on, including cloud-based setup, raw data processing, alignment, QC, and downstream analysis, demonstrating practical proficiency in RNA-Seq workflows.

**ONECUT2 Signaling Role in Normoxia and Hypoxia Condition for Neuroendocrine Prostate Cancer**  
Experimental Design:  
This repository contains datasets and analysis scripts for studying the role of ONECUT2 (OC2) signaling in prostate cancer in Normoxia and Hypoxia conditions. The experiments were performed using LNCaP and PC3 prostate cancer cell lines under normoxia and hypoxia conditions.

**RNA-Seq**  
LNCaP cells: OC2 overexpression vs. control
PC3 cells: OC2 knockdown vs. control
Goal: Identify OC2-regulated gene expression changes under different oxygen conditions.

**ChIP-Seq**  
PC3 cells with Flag-OC2 fusion protein: anti-Flag ChIP to map OC2 binding sites
HIF1A ChIP-Seq: AR-negative PC3 cells under hypoxia ± ONECUT2 or SMAD3 knockdown
SMAD3 and HIF2A ChIP-Seq: PC3 cells under hypoxia
Goal: Map transcription factor binding sites and identify OC2 targets in hypoxia response.

**ChIP-reChIP-Seq**  
Sequential ChIP experiments to study co-binding of transcription factors:
Primary ChIP using anti-SMAD3 antibody
ReChIP using anti-HIF1A or anti-HIF2A antibodies
Goal: Identify genomic regions co-bound by multiple transcription factors.

**Dataset Metadata**  
Field	Value
SRA Study	SRP122890
BioProject	PRJNA416257
Consent	public
Center Name	GEO
Datastore filetype	fastq, run.zq, sra
Datastore provider	gs, ncbi, s3
Datastore region	gs.us-east1, ncbi.public, s3.us-east-1
Library Layout	SINGLE
Organism	Homo sapiens
Platform	ILLUMINA
Release Date	2019-02-27
Version	1
Tumor Type	prostate cancer

Purpose:
Understand the transcriptional network driven by ONECUT2 in neuroendocrine prostate cancer
Integrate RNA-Seq and ChIP-Seq data to elucidate hypoxia-mediated gene regulation
Explore transcription factor interactions under hypoxia conditions

<img width="455" height="535" alt="image" src="https://github.com/user-attachments/assets/9d8f5783-3134-40fe-acbc-3e6cf1e449fe" />

<img width="545" height="356" alt="image" src="https://github.com/user-attachments/assets/348328ac-c1cf-4ad1-83bf-d9f37a6483ef" />  

Link for Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305

# Setting_up_Google_Cloud_Console

**Lets Start with the Workflow, Since my laptop has very low RAM and cannot support the analysis, I chose to perform the work on Google Cloud Platform (GCP). This not only solves the hardware limitation but also provides me with valuable experience in cloud computing.**
*Set up Google Cloud Console*  
- Login with your gmail-id.
- Click on **console** on the top right near your account info.  
  <img width="91" height="57" alt="image" src="https://github.com/user-attachments/assets/b7cf3086-55e8-409a-820a-76354fcd8531" />

- On left corner click the menu (three lines)  
  <img width="194" height="50" alt="image" src="https://github.com/user-attachments/assets/521a04e2-b23f-40de-a8a9-ccf3bda93064" />

- click on Compute Engine and then on the sublist click VM Instances.(This is going to be your Virtual Machine)
- Click **Enable** Computer Engine API
- Then, set the billing account.(paid or free account-as you wish!)
- Click on **+ Create Instance**. You will see this kind of window,  
  <img width="1013" height="587" alt="image" src="https://github.com/user-attachments/assets/611d7316-98b0-4775-bbab-d654f9f49b75" />

- On left, as shown in the screenshot, There is an option called **Machine Configuration** make sure you are there to see these options
- You can name your project as you like, and fill out the machine details as required by your specifications.
  Name: bulk-rna-seq-vm  
  Region and Zone: Which ever is nearest to you.
  Select machine type as E2  
  <img width="734" height="154" alt="image" src="https://github.com/user-attachments/assets/4bf5443d-7f8e-4bbd-a3b2-bc7237165b2b" />
  Scroll down, you will see
  <img width="745" height="273" alt="image" src="https://github.com/user-attachments/assets/99af91f1-6b7e-4527-91ba-67874e28910f" />
  click Standard > your choice of machine you want
  <img width="683" height="257" alt="image" src="https://github.com/user-attachments/assets/e6cec3de-2a32-49a1-af50-681b07ef07a7" />

- On left, click **OS and Storage**
  click **change**  
  <img width="111" height="61" alt="image" src="https://github.com/user-attachments/assets/593c315b-e18b-4a1e-9675-a5bf6ce18af0" />
  For Operating System
  Select Ubuntu and other requirements, it should look something like this,
  <img width="456" height="189" alt="image" src="https://github.com/user-attachments/assets/b2f82550-19df-4e1e-885f-638be142dfca" />

- On left, click **Networking**  
<img width="655" height="480" alt="image" src="https://github.com/user-attachments/assets/9400867c-2462-40c6-b46c-d946c7a25302" />  

-Just click **Create** and its done. Voila!! you have your VM. It will take a couple of minutes and the window will open. You can see your VM name and a green tick near it which means its running. Click on **SSH** button and your terminal will open.  


## Setting_up_the_terminal_and_Installation_of_Tools  
I have written the installation steps first, rather than writing it during the workflow. You can do it one by one during each step of the analysis or however you wish. My intention was to have a clear workflow with clear steps.

**Update and upgrade system packages**
```bash 
sudo apt update && sudo apt upgrade -y
```
Install essential tools
```bash 
sudo apt install -y build-essential wget curl git unzip htop
```
**R Installation, it installs the base R interpreter and standard R libraries**  
```bash 
sudo apt install -y -r base
```
R-server Installation
```bash
sudo apt install gdebi-core -y
```
```bash
wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2025.09.0-387-amd64.deb
```
```bash
sudo gdebi rstudio-server-2025.09.0-387-amd64.deb
```
(alternative without gdebi)
```bash
sudo apt install -y ./rstudio-server-2025.09.0-387-amd64.deb
```
```bash
sudo rstudio-server start
```
After this, you have to set up the firewall rule for your VM. 
Again go to Google Cloud Console and click the menu:  
 <img width="194" height="50" alt="image" src="https://github.com/user-attachments/assets/521a04e2-b23f-40de-a8a9-ccf3bda93064" />  
**On MENU, click VPC Network > Firewall > + create Firewall Rule**
Name: rstudio-server  
Network: default  
Targets: All instances in the network  
Source\Prangs: 0.0.0.0/0  
Protocols and Ports: click specified protocols and ports  
                     enter 8787   
Now on your Browser.. type (external IP you can find on VM Instances)
```bash
http://your_external_IP:8787
```
Then, type in your username and password. R-Studio could be used now. 

**SRA Toolkit Installation - necessary to access, download, and convert sequencing data from the NCBI Sequence Read Archive (SRA) into usable formats for analysis.**  
```bash
sudo apt update && sudo apt upgrade -y
```
```bash
sudo apt install sra-toolkit -y
```
```bash
which prefetch
```
```bash
which fastq-dump
```
**Fastq Installation.**  
```bash
sudo apt install fastqc -y
```
```bash
fastqc --version
```
**Multiqc Installation.**  
```bash
python3 -m pip install --upgrade pip
```
```bash
pip install multiqc
```
```bash
multiqc --version
```
**Trimmomatic Installation**
(on home directory)
```bash
mkdir -p ~/tools
```
```bash
cd ~/tools
```
```bash
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
```
```bash
unzip Trimmomatic-0.39.zip
```
```bash
ls
```
```bash
cd Trimmomatic-0.39
```
you should see this  
<img width="337" height="22" alt="image" src="https://github.com/user-attachments/assets/7d19209d-fbe0-4dca-a66d-86031af7d4cd" />  

**Hisat2 Installation**
```bash
sudo apt install hisat2 -y
```
```bash
hisat2 --version
```
**Samtools Installation**
```bash
sudo apt install samtools
```
**Qualimap Installation**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  
bash Miniconda3-latest-Linux-x86_64.sh
```
```bash
source ~/.bashrc
```
```bash
conda --version
```
```bash
conda install -c bioconda qualimap
```
```bash
qualimap --version 
```
to come out of conda base
```bash
conda config --set auto_activate_base false
```
**Install Subread**
```bash
sudo apt install subread
```
## Workflow Starts:
Finally we have come to the actual workflow of this Bulk RNA Seq Analysis. Few points to remember according to my understanding:   
- To be careful about the directories you will be making. Be careful about where you are saving your results.  
- In any script, the most important thing to keep in mind is the path of the file that the script is fetching and the output directory of the output files, so yes, the path matters.  
-  Being mindful of the files you will be downloading, before downloading double check the working dir.  
- The path and the file names could be different, also the tools if you choose to work with anything different.  
- Let's START!!
   
# **1. DATASET_DOWNLOADING**

I have chosen to download the dataset using a script because of automation, reproducibility, and reliability instead of manually doing it from the database. It would also cause errors. It also helps in proper organization of the data into folders and helps in error handling and resuming the process after already present data. This dataset was publicly available!! 

For that I have used fastq_download.py FASTQ is the raw sequencing reads from a sequencing machine (like Illumina). Each read is basically the DNA or RNA fragment that was sequenced. FASTQ files contain both the sequence and the quality scores for each base. In the downstream I have attached an example of the how the file looks.  
*Command/Script Explanation: Download raw SRA files from NCBI > Convert them to FASTQ files ready for RNA-Seq analysis > Compress and organize the reads automatically > Reports the time taken for each step.*

```bash
#for a new folder named bulk_RNA_analysis
mkdir -p ~/bulk_RNA_analysis 
```
```bash
#if you are working on a VM (GCP), upload the .py file first using the upload option on top right and move it into this dir  
mv ~/fastq_download.py ~/bulk_RNA_analysis/
```
**Please find the fastq_download.py file on this repository in the Data folder**  
[fastq_download.py](/Data)

```bash
#to run the .py file
python3 fastq_download.py
```

<img width="940" height="105" alt="image" src="https://github.com/user-attachments/assets/f7e17d06-d30c-4911-848e-fbaae3289af9" />  
<img width="940" height="70" alt="image" src="https://github.com/user-attachments/assets/23cc73e5-6d96-4855-ab27-ee278f16cae0" />  


# **2. ORGANIZING_DATASET_(CONCATENATE_&_RENAMING)**  

As I am working with publicly available dataset it was done as a preprocessing step before actually running the analysis on it. As we saw it from the data base SRR7179504, SRR7179505, SRR7179506, SRR7179507 are all pieces of the same biological sample (LNCAP_Normoxia_S1) so I just joined it together so the results are not messy later on. 
*Command/Script Explanation: cat joins all files and puts into whatever is written after > . and mv moves the matter of one file to the other eg. matter from SRR7179536_pass.fastq.gz to PC3_Normoxia_S1.fastq.gz*
```bash
#concatenating/joining to combine multiple sequencing runs into one file per sample
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz  > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz  > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz  > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz  > LNCAP_Hypoxia_S2.fastq.gz
```
```bash
#renaming to meaningful biological names
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
```
after carefully viewing and making sure that we have all the required concatenated and renamed files we can remove the SRR files. 
```bash
#removing SRR files for a clean folder
rm SRR*
```
<img width="940" height="49" alt="image" src="https://github.com/user-attachments/assets/ca66da82-5628-460b-b034-76089e479718" />  

This is how the fastq file looks. It has four lines.. 1.Read Id starts with @ 2.Read sequence 3.somewhat same as first line starts with +  4.ASCII characters encoding phred quality scores  

<img width="940" height="85" alt="image" src="https://github.com/user-attachments/assets/f557eb5c-41e1-4267-bfec-7573d8442700" />


# **3. QUALITY_CHECK**
We perform quality checks (QC) on FASTQ files because raw sequencing data is noisy and error-prone. If we don’t check quality early, errors can propagate into alignment, quantification, and differential expression results. Fastqc tells us about read quality across the entire length, Detect adapter contamination, Identify overrepresented sequences, Helps detect technical artifacts, rRNA contamination, or PCR duplicates, Check GC content distribution, Ensure enough sequencing depth & uniformity and helps to Catch issues before wasting compute. 
*Command/Script Explanation: fastq is the tool which I used to do this*
```bash
mkdir -p fastqc_results
```
```bash
#runs fastqc on the fastq.gz files and puts ouput in fastqc_results  
fastqc *.fastq.gz -o fastqc_results/ --threads 8
```
<img width="940" height="145" alt="image" src="https://github.com/user-attachments/assets/eda30212-0af2-4ef2-aeb2-61e758368264" />  

```bash
mkdir -p multiqc_report
```
```bash
#runs multiqc on the files in fastqc_results and puts the output in multiqc_report
multiqc fastqc_results/ -o multiqc_report/
```
<img width="940" height="204" alt="image" src="https://github.com/user-attachments/assets/2b4fb960-a3f0-48fc-a65d-effd3f6b3c80" />  

Quality assessment was performed using FastQC on one representative sample (LNCAP_Hypoxia_S1) to demonstrate the QC workflow. My objective was to get familiarize with the tool and interpret the major metrics (per-base sequence quality, GC content, and adapter contamination). The FastQC report indicated high-quality reads across all positions with no adapter contamination. Therefore, trimming was not required for this or other samples. This decision ensured data integrity was preserved and unnecessary processing was avoided.  
The summary reports for all samples confirmed consistent sequencing quality, allowing the pipeline to proceed directly to alignment.  

**Please find the fastqc and multiqc html files on this repository in the Output > Qualitycheck folder** 
[Results of QC](/Output/Qualitycheck)  


# **4. TRIMMING_IF_NEEDED**
In RNA-Seq pipelines, trimming is often applied to remove adapter contamination and low-quality bases at the read ends, which can otherwise reduce alignment accuracy. To gain familiarity with preprocessing tools, I performed trimming on one FASTQ file (LNCAP_Hypoxia_S1.fastq.gz) using Trimmomatic. The trimmed file was then re-evaluated with FastQC to confirm the improvement in quality metrics. 
```bash
#runs trimmomatic on the specified sample
java -jar ~/tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 \
  ~/bulk_RNA_analysis/fastq/LNCAP_Hypoxia_S1.fastq.gz \
  ~/bulk_RNA_analysis/fastq/LNCAP_Hypoxia_S1_trimmed.fastq.gz \
  TRAILING:10 -phred33
```
<img width="940" height="107" alt="image" src="https://github.com/user-attachments/assets/2583d58b-edd1-4ee6-b601-8f748e1635e0" />  


After trimming, the trimmed output was then re-evaluated with FastQC, which confirmed the expected improvement in quality metrics.

```bash
#runs fastqc on the trimmed fastq file
fastqc LNCAP_Hypoxia_S1_trimmed.fastq.gz
```
<img width="940" height="173" alt="image" src="https://github.com/user-attachments/assets/78a2d0f1-b35e-4ab0-8bc1-d01a890a0f64" />  

**Please find the trimmedqc.sh file to perform QC with fastqc on the trimmed fastq files on this repository in the Data folder**  
[trimmedqc.sh](/Data)  

Since overall FastQC/MultiQC reports showed consistently high Phred scores and negligible adapter contamination, trimming was not applied to the remaining samples, in order to preserve read length and maximize mapping efficiency.  

[Results of QC](/Output/Qualitycheck)  


# **5. INDEXING_AND_ALIGNMENT**  

As the raw reads are processed, its time to align it to a reference genome. Because the genome is very huge its converted to genome indexes using hisat2-build. Before aligning RNA-Seq reads, HISAT2 requires a genome index. Instead of building it from scratch (which is time-consuming), I used a prebuilt GRCh38 human genome index provided by the HISAT2 team.It will align the single-end RNA-Seq FASTQ files against the HISAT2 genome index, then sort and index the BAM outputs with samtools. HISAT2 was used because of Lower memory requirement, Works well on modest VMs and has pre built indices.

Prebuilt indexes available → Saves hours of index building.

Easy integration
In the same directory only,
```bash
#downloads the genome index
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
```
```bash
#unzips the index
tar -xvzf grch38_genome.tar.gz

<img width="936" height="26" alt="image" src="https://github.com/user-attachments/assets/b24050d3-00be-406a-a68c-8b57c8a7d035" />



```
On home dir,
```bash
#unzips the index
cd ~/tools
```
```bash
#Downloads Gene annotation file
wget ftp://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
```
```bash
#unzips the index
gunzip Homo_sapiens.GRCh38.115.gtf.gz
ls -lh Homo_sapiens.GRCh38.115.gtf
```

<img width="940" height="45" alt="image" src="https://github.com/user-attachments/assets/9ab0e67c-c4f0-48d4-ad34-9474116bef06" />  


hisat2alignment.sh file was uploaded on GCP  
**Please find the hisat2alignment.sh file to perform alignment on the processed fastq files on this repository in the Data folder**  
[hisat2alignment.sh](/Data)  

```bash
#change the paths and folder names or any specifics  
nano hisat2alignment.sh
```
```bash
#make the.sh file executable and run
chmod +x hisat2alignment.sh
./hisat2alignment.sh
```
After running HISAT2, I evaluated the mapping efficiency for each sample using the alignment summary reports. All samples aligned efficiently to the reference genome, confirming the integrity and quality of the RNA-seq libraries. In a standard workflow, any sample with a mapping rate below 85–90% would be rechecked for adapter contamination or reference mismatch.  
**FASTQ → (HISAT2) → SAM → (Samtools sort) → BAM → (Samtools index) → BAM + BAI**  

<img width="940" height="148" alt="image" src="https://github.com/user-attachments/assets/a2530284-2fc2-4857-93a4-acccd54e50b5" />  

Aligned RNA-Seq reads from the PC3 Normoxia sample (PC3_Normoxia_S1.bam) to the human reference genome (GRCh38) using HISAT2, and visualized the alignments in IGV. The IGV snapshot shows the coverage and individual read alignments for a selected region of chromosome 1, allowing verification of mapping quality and coverage depth. Junction tracks were also visualized to inspect potential splicing events.  
Interpretation, In this region:  
- Only a few reads align (very low coverage).  
- There are no junction reads, suggesting it’s probably a continuous (non-spliced) region.  
- You can visually confirm that your RNA-seq reads are aligning correctly to the reference genome.

Tracks:
There are three tracks loaded from the BAM file (PC3_Normoxia_S1.bam):
Coverage Track (PC3_Normoxia_S1.bam Coverage) - Shows the depth of sequencing coverage at each position and The grey bars indicate how many reads cover each base.  
Junction Track (PC3_Normoxia_S1.bam Junction) - Shows spliced alignments, i.e., reads that span introns. Often empty for non-spliced regions.  
BAM Track (PC3_Normoxia_S1.bam) - Shows individual aligned reads. Each horizontal grey bar is a read aligned to that region.  

<img width="1366" height="707" alt="image" src="https://github.com/user-attachments/assets/a9ad1ec6-63a2-4e0f-b12d-935212f10b77" />  

**Read Info Pop-up**  
When you click a read, IGV shows detailed info:  
Read name: SRR7179536.8906826.1 – Unique identifier from sequencing.  
Read length: 58 bp – The read spans 58 bases.  
Flags: 256 – Indicates the read is a secondary alignment.  
Mapping: Secondary @ MAPQ 1 – MAPQ (mapping quality) is low (1), secondary alignment.  
CIGAR: 58M – The read aligns perfectly for 58 bases.  
Location: chr1:124,439,623 – Position of the read on the reference.  
Base and Quality: C @ QV 38 – The base at this position is C with a high quality score (38).  
Other fields (NH, NM, etc.) provide technical details about alignment, mismatches, or splicing.  


# **6. QUALITY_CHECK_OF_BAM_FILES**  

After aligning RNA-Seq reads to the GRCh38 reference genome and generating sorted BAM files, I performed quality assessment using Qualimap. This tool evaluates mapping quality, coverage uniformity, strand specificity, duplication rates, and GC bias across the dataset. Running Qualimap ensures that the BAM files are of high quality before downstream analysis, such as gene quantification and differential expression. For RNA-Seq, the BAM files were analyzed against the Ensembl GTF annotation to generate comprehensive reports, including alignment statistics and coverage plots, facilitating early detection of potential issues in the sequencing or alignment process.
**If qualimap is run immediately after installing it using conda it be can run using the following command immediately or first create conda base and then only it could be run**
**Please find the script in the data folder**
[qualimapQC.sh](/Data)  

```bash
#to run qualimap using .sh script
nano qualimapQC.sh
```
```bash
#to run qualimap using .sh script
chmod +x qualimapQC.sh
./qualimapQC.sh
```
If this doesn't work, conda environment doesn't work,
```bash
source ~/miniconda3/etc/profile.d/conda.sh
```
```bash
conda activate base
```
```bash
bash qualimapQC.sh
```

<img width="940" height="259" alt="image" src="https://github.com/user-attachments/assets/f0b6d1d4-547e-4bab-9d2c-382bd05599c9" />  

<img width="940" height="314" alt="image" src="https://github.com/user-attachments/assets/22fea255-48e8-46d3-a9a9-30fa6d976655" />  

For each sample, the pipeline:
- Extracts the sample name from the alignment file.  
- Runs Qualimap RNA-Seq to evaluate mapping quality, gene coverage, and transcript-level metrics using the human GTF     annotation (GRCh38.115).  
- Stores detailed QC reports in a dedicated output folder per sample.  
- Memory optimization ensured efficient processing with 12 GB allocated to Java (required by Qualimap).

  [Results of QC](/Output/Qualimap_Reports)

  Alignment summaries confirmed that all samples mapped efficiently to the genome, and reads genomic origin analysis showed that 89–90% of reads mapped to exonic regions, ~9% to introns, and <2% to intergenic regions. No reads were detected in rRNA regions, indicating clean library preparation. Samples with low exon coverage (<70%) would typically be flagged for reanalysis, such as checking alignment parameters, trimming, or re-examining library preparation. In this dataset, all samples passed the expected coverage thresholds, so no reanalysis was required. 

# **7. GENE_EXPRESSION_QUANTIFICATION_USING_FEATURECOUNTS**  

- Performed gene-level read quantification from RNA-Seq alignments (BAM files) using FeatureCounts.  
- Assigned sequencing reads to genes based on exon overlap using the reference GTF annotation (GRCh38.115).  
- Generated a read count matrix across all samples, forming the input for downstream differential expression analysis. 
- Ensured accurate handling of paired-end and stranded RNA-Seq libraries, producing a reliable dataset for expression profiling.
```bash
mkdir -p quants
```
**Upload the .sh file on GCP and change the path and file name accordingly.**
[featurecounts.sh](/Data)  

```bash
nano featurecounts.sh 
```
```bash
#to run featurecounts using .sh script
chmod +x featurecounts.sh
./featurecounts.sh
```
<img width="940" height="611" alt="image" src="https://github.com/user-attachments/assets/b557857c-d55f-48e6-8e49-0ce962c56c4b" />  
<img width="940" height="94" alt="image" src="https://github.com/user-attachments/assets/6627cfc1-1a85-4231-af47-2257d34fef91" />  

# **8. GENE-LEVEL_COUNT_MATRIX_GENERATION_(Python_&_R)** 

For the downstream analysis, and generating a count matrix I cumulated the featurecounts data for a count matrix and I tried my hands on both Python and R. The script ( Rscript_for_merging_featurecounts.R & countsmatrix_wholedata.py) for both is on Data folder in this repository. 

*Python Implementation:*  
- Used pandas to read each FeatureCounts output file, extracted gene IDs and counts, and iteratively merged all   samples into a single CSV.  
- Automated logging of processing time and number of genes per sample for efficiency tracking.

*R Implementation:*
- Used dplyr to read and process each FeatureCounts file, keeping only the gene ID and counts.
- Merged all samples using full_join to produce a complete counts matrix, ready for DESeq2 or edgeR.

**Note: The quants folder containing the FeatureCounts output files was transferred from the VM to my local machine to reduce cloud computing costs. Since the heavy-duty alignment and counting steps were already completed, all downstream analysis, including merging counts and differential expression, was performed locally in R.**  
[countsmatrix_wholedata.py](/Data)  

```bash
#to install python and run the script 
python3 --version
pip3 install pandas
sudo apt update
sudo apt install python3-pip -y
python3 countsmatrix_wholedata.py
```
<img width="1346" height="467" alt="image" src="https://github.com/user-attachments/assets/76fd9be0-5d67-407a-a037-cd08f4011b62" />  

To run R, I used the Google cloud console only for the generation of read count matrix using the R script Rscript_for_merging_featurecounts.R  

[Rscript_for_merging_featurecounts.R](/Data)  

(ON bash)
```bash
#to activate R on VM 
sudo systemctl status rstudio-server #checks status
```
```bash
#if shows inactive 
sudo systemctl start rstudio-server
```
```bash
#if shows inactive 
sudo systemctl status rstudio-server 
```
<img width="940" height="238" alt="image" src="https://github.com/user-attachments/assets/2a78ee95-d92a-4551-8735-2acce8fcec84" />

**Then , Browser and type URL - http://your_external_IP:8787 > then type in ur ID and password and rn the script**  

**-------- Downstream Processing for Differential Gene Expression Analysis Begins now on R -----------**
# **9. DOWNSTREAM_ANALYSIS**

```bash
#to install and load the libraries and set working directory
install.packages("BiocManager")
BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)
BiocManager::install(version = "devel")
BiocManager::install("DESeq2")
library("DESeq2")
install.packages("tidyverse")
library("tidyverse")
setwd("D:/BIDYA")
```

```bash
#to load the counts matrix file
raw_counts <- read.csv("GSE106305_counts_matrix.csv", header = TRUE, row.names = "Geneid", stringsAsFactors = FALSE)  
head(raw_counts)
```
<img width="972" height="286" alt="image" src="https://github.com/user-attachments/assets/8aa101a4-bc3b-4747-861c-af773974715a" />


```bash
#to give a sum total of all the read counts in each sample
raw_counts <- raw_counts[,sort(colnames(raw_counts))]
colSums(raw_counts)
```
<img width="1050" height="118" alt="image" src="https://github.com/user-attachments/assets/f0f06ac2-aea6-43e6-90e0-3d66be0f8a84" />


```bash
#to define experimental conditions for RNA-seq analysis
condition <- c(rep("LNCAP_Hypoxia", 2), rep("LNCAP_Normoxia", 2), rep("PC3_Hypoxia", 2), rep("PC3_Normoxia", 2))
print(condition)
```
<img width="1011" height="96" alt="image" src="https://github.com/user-attachments/assets/6f742808-62b5-4501-9719-bca1210b93a2" />


```bash
#to define Rows are named after your samples and Column named after condition to later work with DeSeq2
my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(raw_counts)
head(my_colData)
```
<img width="403" height="166" alt="image" src="https://github.com/user-attachments/assets/ec27a164-296d-42fc-ae7d-f592caf6b5c0" />

- To perform differential expression analysis using DESeq2, I first created a DESeqDataSet object using raw count data and sample metadata. This object serves as the foundation for downstream analysis.
```bash
#create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = my_colData,
                              design = ~condition)
dds  ##inspect dataset 
head(counts(dds)) ##preview and check dataset dimensions
dim(counts(dds))
```
<img width="1494" height="573" alt="image" src="https://github.com/user-attachments/assets/3f91073d-5734-4830-9d00-ea6f43df2905" /> 

- Before proceeding with normalization and differential expression testing, I examined how many genes have zero counts in each sample. This step helps identify lowly or non-expressed genes that may be removed to reduce noise. 

```bash
#Quality check
count_matrix <- counts(dds) ##Extracts the raw count matrix from the DESeq2 object dds. Rows = genes, Columns = samples.  
dim(count_matrix) ##Returns the dimensions of the count matrix.
zero_counts_per_gene <- rowSums(count_matrix == 0) ##This line calculates the number of zero counts per gene across all samples.  
count_matrix <- as.data.frame(count_matrix) ##Converts the matrix to a data frame
zero_summary <- table(zero_counts_per_gene) ##Creates a frequency table showing how many genes have a given number of zero counts.
print(zero_summary)
```
<img width="994" height="223" alt="image" src="https://github.com/user-attachments/assets/ac4109cb-c688-41f9-9cdd-4185014861fb" />  

- To enhance interpretability of the raw RNA-seq counts, I merged Ensembl gene IDs with metadata from a GRCh38 annotation file. I removed version numbers from Ensembl IDs to ensure accurate matching, then used a left join to attach gene symbols and biotypes to each row in the count matrix. This allowed me to generate a streamlined, annotated dataset suitable for downstream analysis and visualization.

```bash
#Gene annotation join and Cleanup
annotation_file <- "GRCh38annotation.csv" ##Loads gene annotation file that we made
annotation <- fread(annotation_file, stringsAsFactors = FALSE) 
counts_gse <- read.csv("GSE106305_counts_matrix.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE) ##Loads the raw count matrix 

counts_gse$Geneid <- sub("\\..*$", "", counts_gse$Geneid)
annotation$Geneid <- sub("\\..*$", "", annotation$Geneid) ##removes the version numbers so both files have matching geneids  
annotated_counts <- left_join(counts_gse, annotation, by = "Geneid") %>%  ##performs the joining of both files
  select(Geneid, Genesymbol, Genebiotype, 
         LNCAP_Hypoxia_S1, LNCAP_Hypoxia_S2, LNCAP_Normoxia_S1, LNCAP_Normoxia_S2, 
         PC3_Hypoxia_S1, PC3_Hypoxia_S2, PC3_Normoxia_S1, PC3_Normoxia_S2)  ##keeps only these columns
```
<img width="1648" height="561" alt="image" src="https://github.com/user-attachments/assets/0419c09d-7210-4faf-9586-190c53ca8fe1" />  

- To improve the quality and interpretability of downstream analyses, I applied a two-step filtering approach on the annotated count matrix:  
Biotype filtering:  
I retained genes belonging to biologically relevant categories including protein-coding genes and immune-related gene types (immunoglobulin and T-cell receptor genes). This reduces noise from uninformative or poorly annotated gene biotypes.
Expression filtering:  
Genes with zero counts across all samples were removed to focus the analysis on genes that are at least minimally expressed. This step enhances statistical power and reduces computational burden.  
```bash
#Two step filtering
biotypes_to_keep <- c("protein_coding", "IG_J_gene", "IG_V_gene", "IG_C_gene", "IG_D_gene", "TR_D_gene", "TR_C_gene", "TR_V_gene", "TR_J_gene") ##gene biotypes I want to keep  

filtered_counts <- annotated_counts %>%
  filter(Genebiotype %in% biotypes_to_keep)  ##keeps only rows that matches the biotypes above

filtered_counts$Geneid <- sub("\\..*$", "", filtered_counts$Geneid)
head(filtered_counts, n = 3)  ##Strips Ensembl version suffixes from Geneid again, just in case!!

output_file <- "9biotype_count_matrix.csv"  ##Saves the filtered count matrix (only chosen biotypes) to a CSV file
fwrite(filtered_counts, file = output_file, sep = ",", row.names = FALSE)
zero_counts1 <- rowSums(filtered_counts[, 4:11] == 0)  ##Counts how many samples have zero counts per gene.
zero_summary2 <- table(zero_counts1)
print(zero_summary2) 

############## filtering 2 steps
keep_genes <- zero_counts1 < 7  ##Keeps genes that have counts in at least one sample (i.e., less than 7 zeros across 7 samples)
filtered_counts_nozero <- filtered_counts[keep_genes, ]   
cat("Number of genes after filtering (zeros in <7 samples):", nrow(filtered_counts_nozero), "\n")

new_zero_counts <- rowSums(filtered_counts_nozero[, 4:11] == 0)  ##Recalculates zero count summary after filtering.
cat("New zero counts distribution:\n")
print(table(new_zero_counts))

output_file <- "filtered_biotype_nozero_count_matrix.csv"  ##Exports the dataset after both filters: biotype and zero counts.
fwrite(filtered_counts_nozero, file = output_file, sep = ",", row.names = FALSE)

head(filtered_counts_nozero, n = 3)

dds_filtered <- dds[rownames(dds) %in% filtered_counts_nozero$gene_id, ]  ##Keeps only rows (genes) in the original DESeq2 object dds that passed filters. rownames(dds) are gene IDs.

cat("Dimensions of filtered DESeqDataSet:", dim(dds_filtered), "\n")

removed_genes <- filtered_counts[!keep_genes, ]  ## Get genes filtered out due to zeros in all samples
cat("Biotype distribution of removed genes:\n")
removed_biotype_dist <- table(removed_genes$Genebiotype)
print(removed_biotype_dist)
```
<img width="801" height="171" alt="image" src="https://github.com/user-attachments/assets/bfc70355-d8df-40c5-b05e-85b5df885ffe" />  
Summarizes the frequency of genes having 0, 1, 2, ... zeros across the samples.  

<img width="1124" height="116" alt="image" src="https://github.com/user-attachments/assets/11b6f28a-d1de-4fc0-8b7b-d3aeb7394420" />  
After removing genes that are zero (not expressed) in 7 or more samples, 17,478 genes remain.

<img width="1239" height="189" alt="image" src="https://github.com/user-attachments/assets/f94bd9bf-8556-45a4-9e61-dd588a5a2e79" />  
These are the biotypes of removed genes

- I filtered the gene count data by removing genes with zero counts in 7 or more samples to improve data quality. The filtered dataset was saved, and the DESeq2 object was updated accordingly to ensure consistency for downstream analysis. This step helps focus on reliably expressed genes and reduces noise. The first filter (e.g., by biotype) narrows down genes to relevant types, but it doesn't guarantee all those genes are sufficiently expressed. The second filter targets expression levels, removing genes with too many zeros, improving data quality further.
```bash
print(colnames(filtered_counts)) ##outputs the column names of filtered_counts data frame.
zero_counts <- rowSums(filtered_counts[, 4:11] == 0) ##For each gene(each row), counts how many samples(columns 4 to 11) have a zero count.
zero_summary <- table(zero_counts)
print(zero_summary)
keep_genes <- zero_counts < 7
filtered_counts_nozero <- filtered_counts[keep_genes, ]  ##Keeps genes expressed(non-zero) in at least 2 samples(since fewer than 7 zeros).
cat("Number of genes after filtering (zeros in <7 samples):", nrow(filtered_counts_nozero), "\n")

new_zero_counts <- rowSums(filtered_counts_nozero[, 4:11] == 0  ##Counts zeros again in the filtered dataset to check the new distribution
print(table(new_zero_counts))

output_file <- "filtered_biotype_6.csv"  ##filtered count matrix to a CSV file
fwrite(filtered_counts_nozero, file = output_file, sep = ",", row.names = FALSE)

head(filtered_counts_nozero, n = 3) 

dds_filtered <- dds[rownames(dds) %in% filtered_counts_nozero$Geneid, ]  ##Filters the original DESeq2 dataset (dds) to keep only genes that passed the zero-count filtering.
```

<img width="1527" height="234" alt="image" src="https://github.com/user-attachments/assets/a0e78145-979d-4b78-bc7e-393461615d1e" />  
This shows the filtered read counts where all genes are expressed and contains no zeroes. 

- I analyzed the distribution of gene biotypes in the filtered dataset by calculating the proportion of genes belonging to each biotype category. Using ggplot2, I created a clear and informative bar plot visualizing these proportions, which highlights the dominant gene types retained after filtering. The plot was saved as a high-resolution image for documentation and presentation purposes.
```bash
biotype_counts <- filtered_counts_nozero %>%  ##Using the zero-count-cleaned dataset (filtered_counts_nozero), it counts how many genes belong to each gene biotype (Genebiotype). Then it calculates the proportion and percentage of each biotype relative to the total genes.
  count(Genebiotype) %>%
  mutate(Proportion = n / sum(n),
         Percentage = Proportion * 100) %>%
  rename(Biotype = Genebiotype)  # Rename for clarity in plot
print(biotype_counts)

p <- ggplot(biotype_counts, aes(x = reorder(Biotype, -Proportion), y = Proportion, fill = Biotype)) +    ##create plot with ggplo2
    geom_bar(stat = "identity") +
    labs(title = "Proportion of Genes by Biotype",
        x = "Gene Biotypes",
          y = "Proportion") +
     scale_y_continuous(labels = scales::percent_format(scale = 100)) +  # Show proportions as percentages
    theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Legend is shown now
    scale_fill_brewer(palette = "Set2")  # Use distinct colors 
p
output_plot <- "genebiotype_proportions1.png"
ggsave(output_plot, plot = p, width = 8, height = 6, dpi = 300)
```

<img width="1351" height="693" alt="image" src="https://github.com/user-attachments/assets/c2d84e33-eacb-4471-90d8-1e233e7192bc" />  
As shown in the plot, the most abundant gene biotype in the experiment is protein_coding, indicating that the majority of retained genes after filtering are protein-coding genes.

-  I conducted Principal Component Analysis (PCA) on variance-stabilized gene expression data to examine sample relationships and variability. The PCA plot visualizes the first two principal components, with samples colored by experimental condition. Clear labels were added using ggrepel to avoid overlap. The plot was saved as a high-quality image for reporting and interpretation.
```bash
vsd <- vst(dds_filtered, blind = TRUE)  # blind=TRUE for exploratory PCA  ##transforms the filtered DESeq2 data for variance stabilization, which is good for PCA  
plot_PCA = function (vsd.obj) {
  pcaData <- plotPCA(vsd.obj,  intgroup = c("condition"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
         title = "PCA Plot colored by condition") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}

png(filename = "pcab.png", 
    width = 2000, height = 2000, res = 300)  # adjust width/height as needed
plot_PCA(vsd)
dev.off()
using ntop=500 top features by variance  ##selecting the 500 genes (features) with the highest variability across samples
```
<img width="612" height="607" alt="image" src="https://github.com/user-attachments/assets/6c25f946-c855-40bc-917c-4a98d1da89f4" /> 
The PCA plot reveals distinct clustering of samples according to their experimental conditions, indicating that the treatment strongly influences gene expression profiles. For example, samples under hypoxia cluster separately from normoxia, reflecting biological differences captured by the data. The close grouping of replicates confirms good experimental consistency. No obvious outliers were observed, suggesting reliable data quality.It explains nearly 99% of the variance, clearly separates samples based on the hypoxia condition. This strong separation indicates that hypoxia has a major impact on gene expression, driving most of the variability in the dataset.

- I performed differential expression analysis on the filtered dataset using the DESeq2 pipeline. After fitting the model with DESeq(), I extracted the normalized count data to correct for library size and sequencing depth differences. The normalized counts were saved as a CSV file for downstream analyses and reporting.
```bash
dds <- DESeq(dds_filtered)  ##fits the negative binomial model to the count data using the DESeq() function, which estimates size factors and dispersions necessary for normalization and statistical testing
dds
normalized_counts <- counts(dds, normalized = T)  ##extracted the normalized counts
normalized_counts_df <- as.data.frame(normalized_counts)
write.csv(normalized_counts_df, file = "normalized_counts.csv", row.names = TRUE)
```
<img width="1382" height="356" alt="image" src="https://github.com/user-attachments/assets/26b8def1-abd9-4797-ad6d-22f5a32fe528" /> 
This only carried the gene ids, if needed can be joined with gene annotation file like earlier

- I generated a sample-to-sample distance heatmap using variance-stabilized data to assess the similarity between samples. The heatmap shows hierarchical clustering of samples based on Euclidean distances, confirming expected grouping by experimental condition (e.g., hypoxia vs. normoxia).
```bash
vsd <- vst(dds, blind = TRUE) ##Applies variance-stabilizing transformation
plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj))) ##Computes Euclidean distances between samples using the transformed data (columns = samples)
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd.obj$condition)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(55)
  pheatmap::pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,  clustering_distance_cols = sampleDists, col = colors, fontsize_row = 4, fontsize_col = 4, fontsize_legend = 4, fontsize = 4)
}
png(filename = "sampleheatmap1.png", width = 1000, height = 900, res = 300)  # adjust width/height as needed
plotDists(vsd)
dev.off()
```
<img width="676" height="609" alt="image" src="https://github.com/user-attachments/assets/51a5333a-8ee7-4cda-8357-9b6932c1f25b" />  
A sample-to-sample distance heatmap was generated using variance-stabilized transformed (VST) expression data. The Euclidean distances between samples reveal strong clustering by cell line (LNCaP vs. PC3) and condition (Hypoxia vs. Normoxia). Replicates within each group (e.g., LNCAP_Hypoxia_S1 and S2) show high similarity (darker blue), confirming consistency. Additionally, the dendrogram clearly separates LNCaP and PC3 groups, validating distinct transcriptional profiles across cell lines.

```bash
variable_gene_heatmap <- function (vsd.obj, num_genes = 40, annotation, title = "") {
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)

  # Color palette
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]

  # Get top variable genes
  stabilized_counts <- assay(vsd.obj)
  row_variances <- rowVars(stabilized_counts)
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing=TRUE)[1:num_genes],]
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm=TRUE)

  # Map gene names
  gene_names <- annotation$Genesymbol[match(rownames(top_variable_genes), annotation$Geneid)]
  gene_names[is.na(gene_names)] <- rownames(top_variable_genes)[is.na(gene_names)]
  rownames(top_variable_genes) <- gene_names

  # Get metadata
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL

  # Plot
  pheatmap::pheatmap(
    top_variable_genes,
    color = mr,
    annotation_col = coldata,
    fontsize_col = 8,
    fontsize_row = max(6, 250 / num_genes),
    border_color = NA,
    main = title
  )
}
png(filename = "variable_gene_heatmap.png", width = 1400, height = 2000, res = 300)
variable_gene_heatmap(vsd, num_genes = 40, annotation = annotation)
dev.off()
```
<img width="425" height="609" alt="image" src="https://github.com/user-attachments/assets/542aa7a1-764b-4d99-8317-c6b17d7de0aa" />  

The heatmap displays the expression patterns of the top variable genes across different samples or conditions. Each row corresponds to a highly variable gene, and each column represents a sample. The color intensity reflects the normalized expression level of each gene in each sample, with the color scale indicating relative expression (e.g., blue for low expression, red for high expression). Clustering of rows (genes) and columns (samples) reveals groups of genes with similar expression profiles and samples with similar gene expression patterns, respectively. This visualization helps identify distinct gene expression signatures associated with specific conditions or sample groups, highlighting potentially important genes for further biological interpretation.

- The MA-plot visualizes differential expression results by plotting the log2 fold change (M) of each gene against its average normalized expression (A). It helps to identify genes that are significantly up- or down-regulated between the hypoxia and normoxia conditions in the LNCAP samples. Genes with significant adjusted p-values are highlighted to show the most relevant changes.
```bash
dds_lncap <- dds[, grepl("LNCAP", colnames(dds))]   ##subset the DESeqDataSet (dds) to keep only samples with "LNCAP" in their column names
dds_lncap
dds_lncap$condition <- droplevels(dds_lncap$condition)  ##Removes unused factor levels in the condition variable, so only relevant ones remain.
dds_lncap$condition <- relevel(dds_lncap$condition, ref = "LNCAP_Normoxia")   ##Sets the reference level for the condition factor to "LNCAP_Normoxia" which means means comparisons will be made relative to this baseline.
dds_lncap <- DESeq(dds_lncap)

res_lncap <- results(dds_lncap, contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))  ##Extracts the results comparing LNCAP_Hypoxia vs LNCAP_Normoxia conditions.
res_lncap
summary(res_lncap)
# res0.01 <- results(dds, alpha = 0.5)
# summary(res0.01)
sum(reslncapOrdered$padj < 0.05, na.rm = TRUE)
reslncapOrdered <- res_lncap[order(res_lncap$padj), ]
head(reslncapOrdered)
summary(reslncapOrdered)
write.csv(as.data.frame(reslncapOrdered), file = "DEGs_lncap.csv")
plotMA(res_lncap)
```
<img width="1017" height="661" alt="image" src="https://github.com/user-attachments/assets/27c212ca-3d50-46cd-b14b-abbe23809017" />
The MA-plot visualizes the relationship between the mean expression (average counts) and the log2 fold changes of genes when comparing LNCAP cells under hypoxia versus normoxia conditions. The x-axis (A) represents the average normalized expression level of each gene across all samples. The y-axis (M) shows the log2 fold change in expression between hypoxia and normoxia conditions. Each point corresponds to a gene. Genes with significant differential expression (adjusted p-value < 0.05) are typically highlighted in blue, indicating upregulated or downregulated genes under hypoxia.

-This step highlights genes that are significantly up- or downregulated between hypoxia and normoxia in LNCAP cells. The volcano plot combines fold change and significance to help quickly spot important genes. This aids in focusing on candidates for further study and understanding the biological response to hypoxia.

```bash
# Assuming res_lncap is DESeq2 results object:
reslncapOrdered <- res_lncap[order(res_lncap$padj), ]

# Convert to dataframe
res_df <- as.data.frame(reslncapOrdered)
res_df <- na.omit(res_df)
res_df$gene <- rownames(res_df)

# Categorize genes by regulation status
res_df$regulation <- "Not Significant"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# Load ggplot2
library(ggplot2)

# Plot volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "#FEA405", 
                                "Downregulated" = "purple", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  annotate("text", x = min(res_df$log2FoldChange), y = -log10(0.05) + 0.5,
           label = "padj = 0.05", hjust = 0, size = 3) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression in LNCAP Cells",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(plot.title = element_text(hjust = 0.5))

# Save plot
ggsave("vp_lncap.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
```
<img width="817" height="613" alt="image" src="https://github.com/user-attachments/assets/05acf681-9c4a-48ef-8dab-0a3801c7f053" />  
The volcano plot visualizes gene expression changes between hypoxia and normoxia conditions. The x-axis shows the log2 fold change, indicating how much a gene’s expression increases or decreases. The y-axis shows the statistical significance (adjusted p-value) of those changes.  
Genes with large positive fold changes and low p-values (top right) are significantly upregulated.  
Genes with large negative fold changes and low p-values (top left) are significantly downregulated.  
Genes near the center or bottom are not significantly changed.  

- To investigate the biological pathways affected by hypoxia in LNCaP cells, I prepared gene sets from the MSigDB Hallmark collection using the msigdbr R package. This step is essential to perform pathway-level enrichment analysis with tools like FGSEA, which helps identify key molecular processes altered under hypoxic conditions.
















 
















