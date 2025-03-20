# MDR_ISs

This repository summarizes the code developed for the work **"Plasmids promote multidrug resistance in bacteria through IS-mediated gene inactivation"**.

### Phenotypic resistant mutation rate analysis

To analyze the phenotypic resistant mutation rate of the diverse *Klebsiella pneumoniae* strains tested in the project (**Fig. 1**), we used the R package RSalvador, developed by [Qi Zheng](https://academic.oup.com/g3journal/article/7/12/3849/6027424). We wrote the R script **RSALVADORSCRIPT** to calculate the mutation rate, and later compare between each version (pOXA-48- *vs* pOXA-48+) per strain using maximum likelihood ratio tests, adjusting by Bonferroni-Holm. We followed the same reasoning for comparing the strain containing the IS1-less pOXA-48 (pOXA-48ΔΔIS1).

### Generalized Linear Model of survival during evolutionary rescue

To analyze whether the strains analyzed during the evolutionary rescue (**Fig. 2A**) survived better when carrying pOXA-48 against a wide range of antibiotics, we built a logistic regression model for each antibiotic and strain.
