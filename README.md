# OC_CPTAC_iTRAQ_SWATH
This repository provides R scripts on the analysis of the following publication: 


The analysis contains following steps: 
A) Processing of CPTAC iTRAQ DDA data downloaded from CPTAC data portal to a complete and filtered protein matrix.
    SCRIPT: 1_iTRAQ_protein_processing.R
    INPUT:  Downloaded protein level data from CPTAC data portal analyzed with Common Data Analysis Pipeline (CDAP)
    OUTPUT: CPTAC_DDA_protein_level_CDAP.Rdata containing DDA and DDA_missing0
    
B) Processing the output from the described OpenSWATH, pyProphet, TRIC pipeline with SWATH2stats and mapDIA for protein inference
    SCRIPT: 2_SWATH_preprocessing.R 3_SWATH_norm_batch.R 4_mapDIA_parameter
    INPUT:  feature_alignment.tsv can be downloaded from ProteomeXchange
    OUTPUT: fragments_for_protein_quantification.txt
    
C) Processing the mapDIA output of SWATH data to a final filtered and imputed protein matrix
    SCRIPT: 5_SWATH_processing.R
    INPUT:  fragments_for_protein_quantification.txt
    OUTPUT: CPTAC_SWATH_peptide_level.Rdata contains several intermediate steps, data_peptide, data_pepratio, data_pepratio_log
            CPTAC_SWATH_protein_level.Rdata contains the full unfiltered SWATH protein data data_protein
            CPTAC_SWATH_protein_level_filtered.Rdata contains the SWATH data after filtering (with gene names) SWATH_filtered
            CPTAC_SWATH_imputed.Rdata contains SWATH data with imputed values SWATH_imputed
            Figure S1C
            
D) Analyze SWATH data and iTRAQ data for technical details like peptide and protein variability, completeness and correlation
    SCRIPT: 6a_CPTAC_completeness.R 6b_CPTAC_correlation.R 6c_CPTAC_peptide_variability.R
    INPUT:  CPTAC_SWATH_protein_level.Rdata (data_protein)
            CPTAC_DDA_protein_level_CDAP.Rdata (DDA, DDA_missing0)
            CPTAC_SWATH_imputed.Rdata (SWATH_imputed)
            CPTAC_SWATH_peptide_level.Rdata (data_pepratio)
            PSM files from CPTAC data portal
    OUTPUT: Figure S1A and B
            Figure 2A-C
            CPTAC_SWATH_overlap.Rdata contains matrix with common proteins SWATH_overlap
            CPTAC_DDA_overlap.Rdata contains matrix with common proteins DDA_overlap
            SWATH_iTRAQ_correlations.Rdata
            
E) Comparison of molecular subtype classification
    SCRIPT: 7a_CPTAC_Classification.R 7b_CPTAC_WGCNA.R 7c_CPTAC_heatmap.R
    INPUT:  CPTAC_SWATH_overlap.Rdata
            CPTAC_DDA_overlap.Rdata
            CPTAC_DDA_protein_level_CDAP.Rdata (DDA_missing0)
            CPTAC_SWATH_protein_level.Rdata (data_protein)
            classification_result_template.csv
            clinical_data.csv
    OUTPUT: SWATH_classification_workflow.Rdata
            DDA_classification_workflow.Rdata
            classification_result.Rdata
            Figure 3
            Figure S2A
            
F) Bootstrapping algorithm for classification stability and mesenchymal stability
    SCRIPT: 7d_CPTAC_stability.R 
    INPUT:  classification_result.Rdata
            CPTAC_DDA_protein_level_CDAP.Rdata
            CPTAC_SWATH_imputed.Rdata
            SWATH_classification_workflow.Rdata
            DDA_classification_workflow.Rdata
    OUTPUT: Figure 4A-D
            Figure S2B

G) Differential expression analysis of previously identified molecular subgroups
    SCRIPT: 7e_CPTAC_class_comparison.R
    INPUT:  CPTAC_DDA_WGCNA.Rdata
            CPTAC_SWATH_WGCNA.Rdata
            SWATH_iTRAQ_correlations.Rdata
            CPTAC_SWATH_imputed.Rdata
            CPTAC_DDA_protein_level_CDAP.Rdata
    OUTPUT: Figure 5
            Figure 6
    
H) Analysis of HRD vs. non-HRD subgroups
    SCRIPT: 8_
    INPUT:
    OUTPUT:
