# MultiNeuroDisease_TimeCourseAnalysis

Workflow: Preprocessing.R -> Time Course DEG Analysis.R -> Functional Analysis.R

**Preprocessing.R**
1. Initializes and formats annotation matrix with the correct sample IDs and column data

    *NOTE: Code has been optimized for "Ismael-QC-SHU_14219_B01_EXS_RNA_batchIdentifiers.xlsx"*
    
2. Initializes and formats count matrix

    *NOTE: Code has been optimized for "results_Sep2120_cmp_same_age/featureCounts_count_matrix.txt"*

3. Creates DESeq Data Set from count matrix --> estimate size factors (each column is divided by the geometric means of the rows. The median of these ratios is the size factor for this column) -> normalization 

4. Normalized counts are then exported for use in "Time Course DEG Analysis.R"


**Time Course DEG Analysis.R**

**Functional Analysis.R**
