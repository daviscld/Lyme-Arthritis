# Lyme-Arthritis
This repository contains the code corresponding to the publication: "Borrelia-specific antibody profiles and antibody-dependent complement deposition in joint fluid distinguish antibiotic-refractory Lyme arthritis from antibiotic-responsive arthritis" (submitted).

Abstract: Lyme disease, caused by the spirochete Borrelia burgdorferi, is the most common vector-borne disease in the United States Lyme arthritis is the most common feature of late disseminated disease, seen in approximately 60% of untreated individuals. While most Lyme arthritis resolves with oral or IV antibiotics, termed “antibiotic-responsive” arthritis, a small percentage of individuals develop progressive synovitis despite both oral and IV antibiotic therapy, called “antibiotic refractory” Lyme arthritis (LA). Patients with antibiotic-refractory LA are usually treated with immunosuppressive disease-modifying antirheumatic drugs (DMARDs). The primary drivers behind antibiotic refractory disease are likely multifactorial and remain incompletely understood. More specifically, it remains unclear whether antibodies, either specific for Borrelia or autoantibodies, may act as biomarkers or play a mechanistic role in driving pathogenesis in the joint. We performed a matched, cross-compartmental comparison of the antibody profiles from the blood and joint fluid of individuals with antibiotic responsive (n = 11) or antibiotic refractory arthritis (n= 31). We found that, while serum antibody profiles poorly discriminated responsive from refractory LA patients, a discrete profile of B.burgdorferi-specific antibodies in the joint fluid was able to discriminate antibiotic responsive from refractory LA patients. Cross compartmental comparison of antibody glycosylation, IgA1, and antibody-dependent complement deposition (ADCD) revealed differences across compartments, with more poorly coordinated Lyme-specific humoral responses and increased antibody-dependent complement deposition in antibiotic-refractory disease. These data point to B.burgdorferi-specific serological markers that may support the early stratification and clinical management of LA. Moreover, these findings point to immune complex driven complement activation as a key mechanism underlying persistent disease.

Last updated: 2023-11-20

This repository is administered by Christine Wiggins (daviscld@mit.edu)

Required Software:

R version 4.2.3

Required Packages:

pROC_1.18.0
gridExtra_2.3
forcats_1.0.0
caret_6.0-94
lattice_0.21-8
ggraph_2.1.0
igraph_1.4.2
corrr_0.4.4
ropls_1.30.0
inlmisc_0.5.5 dplyr_1.1.2
reshape2_1.4.4
robustbase_0.95-1
gplots_3.1.3
purrr_1.0.1
tidyr_1.3.0
stringr_1.5.0
ggpubr_0.6.0
ggplot2_3.4.2
summarytools_1.0.1
RColorBrewer_1.1-3
viridis_0.6.2
viridisLite_0.4.2
pheatmap_1.0.12
readxl_1.4.2
systemsseRology_1.1

This is the tested environment, but program may run with other specifications.

Analysis was carried out under macOS Ventura 13.2.1 with R run via R Studio Version 2023.03.0+386

To generate any figures of interest, run the appropriate Data_load.R first (either the one for Figs 1, 3, 4, 5 & S1/2/5/6, or the one for Figs 2, 6, S3/4). Then run Load_functions.R, followed by the file for the figure you would like to generate.
