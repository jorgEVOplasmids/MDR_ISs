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

1. We downloaded the complete metadata, selecting those from the *Enterobacteriaceae* family ([Databases_analysis/KO_impact_Phenotypic_AMR/Download_Genomes.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/Download_Genomes.sh)).
2. We filtered out inconsistent data. We deleted genomes without phenotypic information, labeled as both Resistant and Susceptible without M.I.C. information or genomes with higher values of Susceptible M.I.C. than Resistant M.I.C. If there are multiple data for a single genome and a single phenotype, we only keep the highest M.I.C values for Susceptible data and the lowest for Resistant data. ([Databases_analysis/KO_impact_Phenotypic_AMR/General_Filtering.sh](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/General_Filtering.sh)).
3. ANOTATE
4. R CODE PRINCIPAL
5. The final plotting and statistical analysis can be found in [Databases_analysis/KO_impact_Phenotypic_AMR/KO_analysis.R](https://github.com/jorgEVOplasmids/MDR_ISs/blob/main/Databases_analysis/KO_impact_Phenotypic_AMR/KO_analysis.R).
