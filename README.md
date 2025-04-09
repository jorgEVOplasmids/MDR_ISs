# MDR_ISs

This repository summarizes the code developed for the work **"Plasmids promote multidrug resistance in bacteria through IS-mediated gene inactivation"**.

## Statistical analysis

### Phenotypic resistant mutation rate analysis

To analyze the phenotypic resistant mutation rate of the diverse *Klebsiella pneumoniae* strains tested in the project (**Fig. 1**), we used the R package RSalvador, developed by [Qi Zheng](https://academic.oup.com/g3journal/article/7/12/3849/6027424). We wrote the R script **RSALVADORSCRIPT** to calculate the mutation rate, and later compare between each version (pOXA-48- *vs* pOXA-48+) per strain using maximum likelihood ratio tests, adjusting by Bonferroni-Holm. We followed the same reasoning for comparing the strain containing the IS1-less pOXA-48 (pOXA-48ΔΔIS1).

### Generalized Linear Model of survival during evolutionary rescue

To analyze whether the strains analyzed during the evolutionary rescue (**Fig. 2A**) survived better when carrying pOXA-48 against a wide range of antibiotics, we built a logistic regression model for each antibiotic and strain. We used the *glm* R function implemented in the **GLMSCRIPT** to predict the survival probability and compare between genotypes.

### Log-rank test of survival during experimental evolution

To compare the survival of the KPN08 samples with and without pOXA-48 (**Fig. 2B**) we used the *survival* and *survminer* R packages implemented in the **EESCRIPT**. We used a Log-Rank test per antibiotic comparing the pOXA-48- and pOXA-48+ strains survival against the different antibiotics. We specifically used the modified Log-Rank test with the Fleming-Harrington method, which assigns larger weights to final times of the experiment.

## Databases analysis

### Distribution of KOs in AMR determinants by plasmid-borne ISs

To analyze the distribution of KOs in AMR related genes by plasmid-encoded ISs, we analyzed the genomes available at the [BV-BRC](https://www.bv-brc.org/). To retrieve these, we downloaded the complete metadata of the genomes available from the specific [database section](https://www.bv-brc.org/view/Bacteria/2#view_tab=genomes), filtering by genome completeness (*complete genomes*), quality (*good*) and family (*Enterobacteriaceae*). 

### Analysis of KOs by ISs impact on AMR
To investigate the relationship between the disruption of AMR-related genes and their effects on resistance phenotypes, we analyzed the genomes available from the NCBI Pathogen Detection Database (https://www.ncbi.nlm.nih.gov/pathogens/ast/#). Specifically, we examined this correlation for our experimental antibiotics — Chloramphenicol, Ciprofloxacin, Colistin (and other polymyxins), Kanamycin, and Fosfomycin — as well as for betalactams due to their clinical relevance. 

1. We downloaded the complete metadata, selecting those from the *Enterobacteriaceae* family and using the NCBI FTP to download the whole genomes. The code can be found in [Download_Genomes.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/Download_Genomes.sh).
2. We filtered out inconsistent data by removing genomes without phenotypic information, genomes labeled as both Resistant and Susceptible without corresponding M.I.C. data and genomes where the Susceptible M.I.C. value was higher than the Resistant M.I.C. one. In cases with several coincidences for the same genome and phenotype, we kept only the highest M.I.C. values for Susceptible data and the lowest for Resistant data. The code can be found in [General_Filtering.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/General_Filtering.sh).
3. ANNNOTATION
4. R CODE PRINCIPAL
5. Finally, we performed Fisher's test for all antibiotics (with more than 100 samples and more than 5 samples per phenotypic condition) to determine whether or not there is a significant relationship between the phenotype and the presence of KO AMR-related genes. The code can be found in [KO_analysis.R](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/KO_analysis.R).
