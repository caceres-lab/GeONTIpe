# GeONType
A Snakemake workflow to genotype inversions and to discover structural variants using Oxford Nanopore Long Reads 

## Pipeline Overview
GeONType is a Snakemake workflow designed to genotype inversions, specifically those mediated by inverted repeats, using Oxford Nanopore long reads. The pipeline encompasses various steps, from sample downloading to inversion genotyping. It employs a range of software tools and custom scripts. The core concept in genotyping involves capturing reads that span the inversion breakpoints and determining their orientation through probe blasting. 
Despite the pipeline's initial purpose, it can also be used for:
  - Genotyping inversions not mediated by inverted repeats (previous creation of artificial inverted repeats of X bp around the breakpoints of the inversion)
  - Using PacBio reads (adjusting parameters to improve performance)
  - Using assemblies (adjusting parameters to improve performance)

For the proper functioning of this workflow, a config.yaml file is included for users to adjust each parameter, along with a yml file for use with conda, which includes all the necessary software. The pipeline can be used in its entirety or only one section (parameter option:complete/genotype), depending on whether you need to download and map your reads or only genotype inversions.

## Installation
To run snakefile locally you must have execute:
```
git clone https://github.com/RMoreiraP/GeONType
cd GeONType
conda env create -f environment/Inversiones.yml
conda activate Inversiones
```
## Inputs
### Inputs for Running the Complete Pipeline
If you need to execute the entire pipeline (from downloading to genotyping), the required inputs are:
  - coords.txt: Specifies the inversions to be tested. Each column must be tab-separated and include: Inv_ID, chromosome, BP1_start, BP1_end, BP2_start, BP2_end (chromosome should be 1, 2, ..., X, Y).
  - ListaTodos.txt: Contains all individuals with their respective URLs for downloading the necessary data in FASTQ or FASTA format. Must be tab-separated with the format: Individual, URL (if there are multiple URLs per individual, each URL should be associated with the same individual).
  - gender: Provides the sex of each sample, separated by a tab: Sample, sex (Women/Men).
  - conversion.txt: Converts chromosome names, e.g., chr1, chr2, ... (modifiable based on the reference genome used and its chromosome names).

### References:
  - Std_ref.fa (reference genome to use, e.g., hg38)
  - t2t_ref.fa (T2T or other secondary reference genome)
  - All_snp.vcf.gz (SNPs used in the analysis)
  - Ref.fa (reference used when cram file is used)

### Inputs for Using Pre-Mapped Files
If mapped files are already available, either locally or on an FTP server, you do not need ListaTodos.txt. Instead, provide:
  - Individuals.txt: List of all the samples to be tested.

### Directory Structure
  - Place all input files in the Infor directory.
  - Place all necessary references in the Infor/Reference directory.

## Output
To track the workflow steps, you can enter in the tracking_pipeline directory, where files for each step will be created. The main outputs will be found in the Genotyping directory and the directories of each analyzed sample.

In the Genotyping directory, for each analyzed inversion, you will find all the generated probes. Note that if fewer than three probes are found, the specific directory for that inversion will not be created. This information can be checked in tracking_pipeline/created_probes.txt.

Within each individual’s directory, there will again be a Genotyping directory, containing each inversion in its directory. Inside each inversion directory, you will find:
  - Inversion.png: An image displaying the reads and the generated probes.
  - Inversion_Genotype.txt: A file with information on each read and its respective genotype.
Both outputs are shown here as examples.

| Imagen | Datos |
|--------|-------|
| <img src="https://github.com/RMoreiraP/GeONType/blob/main/example/HG00096/Genotyping/INV00000035/INV00000035.png" alt="Inversion.png" width="400"/> | ```plaintext
Read	Dist	SignGT	ProbGeno
0ee43437-b1c9-4d73-9f57-0419ff7525ff	Std	Std	
10a3d390-28f0-4291-8553-a5ae61c7a7f9	Std	Std	
118bd44d-a00b-4bd6-85c6-bfea8d297d33	Std	Std	
1c262ba9-6b6c-48c9-9322-04bc5ea7426e	Std	Std	
2333d74a-fa95-49ac-8637-cca73188082c	Std	Std	
2347ce7d-e9cf-417e-b38b-91316657d4ec	Std	Std	
27af50dd-e8fe-4d1c-b599-e7020e63dd5b	Std	Std	
37a38fca-7cc2-40e8-95f7-760ac25742e6	Std	Std	
41a97261-02e1-422e-b8d5-4a397c8a0b59	Std	Std	
4d1fcc3a-00e2-48b5-b733-a72e59d167bc	Std	Std	
5a829f2c-bc5f-4519-817a-d6b97a94fac6	Std	Std	
5e345c18-9218-4c42-89ca-0b3d553dfadc	Std	Std	
5e7b1dbd-7c1d-421a-9ebe-8d6992353472	Std	Std	
6182a2d8-1cf1-4bf6-ab19-11b85a54784a	Std	Std	
7354255c-1411-42dc-8044-2dd94dc63121	Std	Std	
76218d15-365f-46ea-a8f8-0a8250ac4e9e	Std	Std	
86e55b6f-f901-49b1-8745-517def763d99	Std	Std	
8f4d9b4d-6e6a-457c-9978-519eb01e16cc	Std	Std	
9bbff0f0-b523-41e2-b7ca-72517a1db3b1	Std	Std	
a86db99a-6a1d-40a9-b2e3-af349ea8d149	Std	Std	
abd10e7d-2a52-45f6-b9e6-799cb23fceae	Std	Std	
bc882571-11c2-43c4-bc8f-7fe8c3c13a84	Std	Std	
c254af39-5881-4030-a926-ef7914e9c177	Std	Std	
c88a592f-d59e-4867-808c-56c79e689b22	Std	Std	
d5a96459-d0ed-4180-8222-3f9fe1fce7e3	Std	Std	
dc4de944-e2b1-4a3e-9730-531410e2a47f	Std	Std	
e508c28c-0c94-49c0-b995-43f394baf744	Std	Std	
e6437356-54b2-4c8c-9263-4e47632d580e	Std	Std	
eb44b2db-c36d-4ed6-b282-971809e85478	Std	Std	
f44ee506-7173-4d36-9816-782aaf516a79	Std	Std	
f647c213-ba9e-4477-ab6d-4f1cd30dc1b1	Std	Std	
FinalGenotype	Std		1
``` |

## Usage
To execute the workflow, the first step is to access the config file to modify the main directory paths:

```
{
    ## Execution options
    "option": "genotype",

    ## Directories
    "directory": "/main/route",
    "Ref": "/main/route/Infor/Reference/Std_ref.fa",
    
    ## Split of files
    "number_of_splits": 30,
    
    ## Filtering reads on statistics
    "Quality_reads": 7,
    "Read_length": 5000,
    
    ## Mapping of reads
    "Kind_reads": "map-ont",
    "Inversion_detection": "400,100",
    "SV_detection": "100,1000",
    "Mapq_filter": 20,
    
    ## Creating input to genotypes and probes
    "flanking_region": 100000,
    "extra_region": 500000,
    "size_probes": 300,
    "overlap_probes": 290,
    "seq_ref": "/main/route/Infor/Reference/chm13v2.0.fa",
    "extra_route": "/main/route/Infor/Reference",
    "min_uni_seq": 50,

    ## Genotyping inversions
    "local": "Y",
    "ftp_route": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38",
    "task": "megablast",
    "coverage_of_probes": 90,
    "identity_probe": 85,
    "Ref_to_cram": "/main/route/Infor/Reference/reference_used_to_cram_creation",

    ## Genotype classification
    "Minor_allele_proportion": 0.15,
    "Low_confidence_proportion": 0.05,

    ## Alternative SVs classification
    "Difference_respect_expected": 0.05,
    "Minimum_size_SVs": 200,
    "Ratio_to_split_SV": 0.1,
    
    ## SNP analysis
    "SNP_selection_route": "/main/route/Infor/Reference/All_snp.vcf.gz",
    "Quality_by_base": 1,
    "Min_freq_snps": 0.15,
    "Min_coverage": 1,
    "Distance_of_tree": 0.5,
}
```
If you have the mapping files in a local directory or perform them with workflow:
  - Parameter local:"Y"
  - Parameter ftp_route:"X" 

However, if you wish to work from mapping files located in FTPs:
  - Parameter local:"N"
  - Parameter ftp_route:"/route/ftp" (The mapping files must have the following naming convention: individual.bam; If this is not possible, you can modify the script located at Scripts/Genotyping_scripts/runall.sh to adapt the naming)

If you want to keep all parameters after modifying the main paths, proceed with the following steps:
```
snakemake --cores all
```

## Notes
Here is an explanation of the function of some adjustable parameters:
  - Option: Determines the steps of the workflow. Options:
    - "complete": Full workflow from download to genotype.
    - "genotype": Genotype inversions using existing BAM files.
    - "genotype_wt_probes": Genotype inversions using existing BAM files and using designed probes.
  - Number_of_splits: This parameter determines how many fragments each downloaded file will be divided. If you have limited RAM and prefer to launch many jobs simultaneously, it's preferable to increase this number. Conversely, if you have ample RAM, you can set this parameter to 1 to map the entire file without splitting it (Default: 30).
  - Quality_reads: Minimum quality threshold required for reads to pass filtering. (Default: 7)
  - Read_length: Minimum desired length of reads in base pairs. (Default: 5000)
  - Kind_reads: Equivalent to the -ax parameter of minimap2. (Default: map-ont)
  - Inversion_detection: Equivalent to the -z parameter of minimap2. (Default: 400,100)
  - SV_detection: Equivalent to the -r parameter of minimap2. (Default: 100,1000)
  - Mapq_filter: Filtering of reads based on mapping quality. (Default: 20)
  - Flanking_region: Additional sequence added upstream and downstream using the breakpoint coordinates as reference. For high-complexity regions or inversiones mediated by inverted repeats, values of 50,000-100,000 are recommended. For low-complexity regions, values of 10,000-50,000 are sufficient. (Default: 100,000)
  - Extra_region: Additional sequence used to assess the uniqueness of generated probes. Should be adjusted based on read type (e.g., ONT, PacBio Hi-fi). Longer reads require more extra region. (Default: 500000)
  - Size_probes: Size of the generated probes. Adjust based on base quality; higher size is needed in reads with low-quality bases. Range of 250-500 bp for ONT, 100-300 bp for PacBio Hi-fi. (Default: 300)
  - Overlap_probes: Overlap between each tested probe. (Default: 275)
  - Min_uni_seq: Amount of sequence not classified as mobile elements. Lower values increase the likelihood of non-specific probes affecting genotyping. (Default: 70)
  - Task: Parameter -task for blast. (Default: megablast)
  - Coverage_of_probes: Parameter -qcov_hsp_perc for blast. (Default: 90)
  - Identity_probe: Parameter -perc_identity for blast. (Default: 85)
  - Minor_allele_proportion: Proportion of reads of the minor allele required to be considered a real allele. (Default: 0.15)
  - Low_confidence_proportion: Proportion of reads of the minor allele to consider inconclusive genotypes. Used together with Minor_allele_proportion for inconclusive genotype thresholds. If low_confidence_proportion is the same value as minor_allele_proportion, not low-confident genotypes will be generated (Default: 0.05)
  - Difference_respect_expected: Minimum percentage difference to consider an additional structural variant. Higher values are needed for reads with more errors.(Default: 0.05)
  - Minimum_size_SVs: Minimum size of detection of structural variants. (Default: 200)
  - Ratio_to_split_SV: Percentage difference to determine if two similar found structural variants should be classified as one or two. (Default: 0.1)
  - Quality_by_base: Minimum quality at each base when using samtools mpileup. Higher values extract fewer but more likely correct bases. (Default: 1)
  - Min_freq_snps: Minimum frequency of the minor allele found by SNPs. (Default: 0.15)
  - Min_coverage: Minimum depth required at each SNP position found. (Default: 1)
  - Distance_of_tree: Minimum distance in the generated read cluster to consider reads as different. (Default: 0.5)

These parameters can be adjusted in the config file according to your specific data and analysis needs.

## How to cite
Ricardo Moreira-Pinhal, Konstantinos Karakostis, Illya Yakymenko, Oscar Conchillo, Maria Díaz-Ros, Andrés Santos, Miquel Àngel Senar, Jaime Martínez-Urtaza, Marta Puig, Mario Cáceres, Complete characterization of human polymorphic inversions and other complex variants from long read data.

