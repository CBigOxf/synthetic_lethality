# Synthetic lethality analysis for AGO2
Repository for the experiments performed.

### Prerequisites
R version 4.0 and above should be working as normal. Required packages are loading within the functions themselves.
Data-sets should be downloaded from their corresponding websites.
For TCGA_BRCA and TCGA_METABRIC please visit https://www.cbioportal.org/, for CCLE: https://depmap.org/portal/

## Running the functions

Some preprocessing might be required to use the data-sets, such as transposition of the expression matrices. 
Please see the functions to gain insight how the mutational and expression matrices have been used. You should
create a subfolder within the parent directory where the functions will reside named : "inputs" (or change the 
way the functions use the directory to your custom selection). The following are files expected to be found for
the functions to work:

* CELL_LINES (folder including a file containing .csv files which correspond to cell-lines mapped to specific cancer types
* brca_CNA.txt (can be downloded from cBioPortal)
* brca_expression.Rdata (expression matrix transposed from cBioPortal, samples as rows, genes as columns)
* brca_mutations.txt (can be downloded from cBioPortal)
* calls.txt (GISTIC calls for AGO2, MYC from cBioPortal on CCLE)
* CCLE20Q2.Rdata (expression,mutation matrix from depmap, transposed as required)
* cell_line_names.csv (a file containing the mapping of BROAD_ID cell line names to their actual names)
* metabric_CNA.txt (can be downloded from cBioPortal)
* METABRIC_expression.Rdata (expression matrix transposed from cBioPortal, samples as rows, genes as columns)
* metabric_mutations.txt (can be downloded from cBioPortal)

MYC neutral cases or not (MYC diploidy) can be run using the same functions and respectivelly commenting out the line of code each time that
incudes the filter, such as "dplyr::filter(MYC == 0)". 

### References

