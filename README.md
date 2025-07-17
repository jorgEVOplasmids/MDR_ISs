# MDR_ISs

This repository summarizes the code developed for the work **"Plasmids promote antimicrobial resistance through Insertion Sequence-mediated gene inactivation"**.

## Statistical analysis

### Phenotypic resistant mutation rate analysis

To analyze the phenotypic resistant mutation rate of the diverse *Klebsiella pneumoniae* strains tested in the project , we used the R package RSalvador, developed by [Qi Zheng](https://academic.oup.com/g3journal/article/7/12/3849/6027424). We wrote the R script [RSalvador_COL_RIF_mutrate.R](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Statistical_analysis/Mut_rate_analysis/RSalvador_COL_RIF_mutrate.R) to calculate the mutation rate, and later compare between each version (pOXA-48- *vs* pOXA-48+) per strain using maximum likelihood ratio tests, adjusting by Bonferroni-Holm. We followed the same reasoning for comparing the strain containing the IS1-less pOXA-48 (pOXA-48ΔΔIS1). The results of these analyses are shown in **Fig. 1B** and **S. Fig. 1**, summarized in **S. Table 2** and can be plotted using the R script [mut_rates_plot.R](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Statistical_analysis/Mut_rate_analysis/mut_rates_plot.R).

### High-throughput AMR analysis

To analyze whether the strains analyzed during the high-throughput antimicrobial susceptibility assays (**Fig. 2A**) survived better when carrying pOXA-48 against a wide range of antibiotics, we built a logistic regression model for each antibiotic and strain. We used the *glm* R function implemented in the [high_throughput_AMR.R](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Statistical_analysis/AMR_HT_analysis/high_throughput_AMR.R) script to predict the survival probability and compare between genotypes, which also includes the code to plot the inset model plots shown in **Fig. 2A**.

## Genomic Analysis

### Assembly of reference genomes

We assembled the reference genomes of the clinical strains ancestral genomes using hybrid assembly of ONT and Illumina data. We used Unicycler v.0.5.0 and checked the genome completeness using Bandage. Then, we annotated the assembled genomes using Bakta v.1.9.3. The commands for assembling and annotating the reference genomes can be found in the script [reference_assembly_and_annotation.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Genome_Analysis/reference_assembly_and_annotation.sh).

### Variant calling

We analyzed the mutations in the resistant colonies to diverse antibiotics isolated from fluctuation assays, high-throughput antibiotic susceptibility assays, and community assays. To perform the analysis, we used breseq v.39.0, using as input the Illumina reads of the resistant clones against the reference genome of the corresponding ancentral strain. The commands used can be found in [breseq_loop.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Genome_Analysis/breseq_loop.sh).

## Databases analysis

### Distribution of KO-AMR genes by plasmid-borne ISs

To analyze the distribution of KOs in AMR related genes by plasmid-encoded ISs, we analyzed the genomes available at the [BV-BRC](https://www.bv-brc.org/). To retrieve these, we downloaded the complete metadata of the genomes available from the specific [database section](https://www.bv-brc.org/view/Bacteria/2#view_tab=genomes), filtering by genome completeness (*complete genomes*), and quality (*good*). We did so through FTP, and then built the plasmid, IS and AMR determinants local ABRicate databases (using the command *abricate --setupdb*). Then, we looped over the genomes to detect IS and partial or complete AMR determinants using the code specified in [BVBRC_parsing.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/AMR_determinant_distribution/BVBRC_parsing.sh).

Then, we merged the information from the ABRicate results and identified those IS disrupting AMR genes using the R script **DBSCRIPT1**. The statistical analyses and representation of the BV-BRC analyses plots are summarized in **DBSCRIPT2**.

### Analysis of KOs by ISs impact on AMR phenotype
To investigate the relationship between the disruption of AMR-related genes and their effects on resistance phenotypes, we analyzed the genomes available from the NCBI Pathogen Detection Database (https://www.ncbi.nlm.nih.gov/pathogens/ast/#). Specifically, we examined this correlation for our experimental antibiotics — Chloramphenicol, Ciprofloxacin, Colistin (and other polymyxins), Kanamycin — as well as for fosfomycin and all betalactams due to their clinical relevance.

1. We downloaded the complete metadata, selecting those from the *Enterobacteriaceae* family and using the NCBI FTP to download the whole genomes. The code can be found in [Download_Genomes.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/Download_Genomes.sh).
2. We filtered out inconsistent data by removing genomes without phenotypic information, genomes labeled as both Resistant and Susceptible without corresponding M.I.C. data and genomes where the Susceptible M.I.C. value was higher than the Resistant M.I.C. one. In cases with several coincidences for the same genome and phenotype, we kept only the highest M.I.C. values for Susceptible data and the lowest for Resistant data. The code can be found in [General_Filtering.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/General_Filtering.sh).
3. We used Abricate to detected ISs and disrupted genes that could confer resistance. To ensure the detection of both considering the nature of our samples, we set the minimum alignment coverage to 5 and the minimum nucleotide identity to 50. The code can be found in [Annotation_Genomes.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/Annotation_Genomes.sh).
4. With the results of the ISs, the genes in the Megares database (plus our experimental targets) and the main database, we determined if there was a coincidence in a broken gene with an adjacent IS. In order to be considered a broken gene the `%Coverage` must be under 100%.  We ran this code for every antibiotic individually and then we compiled the result of the targets (both the gene per se and the `Element` or family in which is included). The code can be found in [Main_analysis.R](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/Main_Analysis.R).
5. Finally, we performed Fisher's test for all antibiotics (with more than 100 samples and more than 5 samples per phenotypic condition) to determine whether or not there is a significant relationship between the phenotype and the presence of KO AMR-related genes. The code can be found in [KO_analysis.R](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/KO_analysis.R).

