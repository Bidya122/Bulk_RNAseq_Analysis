# Bulk_RNAseq_Analysis
This repository contains scripts, workflows, and documentation for performing bulk RNA sequencing (RNA-seq) analysis.

**ONECUT2 and Hypoxia Signaling in Neuroendocrine Prostate Cancer**  
Experimental Design:  
This repository contains datasets and analysis scripts for studying the role of ONECUT2 (OC2) in hypoxia signaling in prostate cancer. The experiments were performed using LNCaP and PC3 prostate cancer cell lines under normoxia and hypoxia conditions.

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

# Setting up Google Cloud Console

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


## Setting up the terminal and Installation of Tools.  

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












