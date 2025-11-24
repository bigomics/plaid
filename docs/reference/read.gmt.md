# Read GMT File

Read data from a GMT file (Gene Matrix Transposed). The GMT format is
commonly used to store gene sets or gene annotations.

## Usage

``` r
read.gmt(gmt.file, dir = NULL, add.source = FALSE, nrows = -1)
```

## Arguments

- gmt.file:

  Path to GMT file.

- dir:

  (Optional) The directory where the GMT file is located.

- add.source:

  (optional) Include the source information in the gene sets' names.

- nrows:

  (optional) Number of rows to read from the GMT file.

## Value

A list of gene sets: each gene set is represented as a character vector
of gene names.

## Examples

``` r
# \donttest{
# Read GMT file (requires file to exist)
gmt_file <- system.file("extdata", "hallmarks.gmt", package = "plaid")
if (file.exists(gmt_file)) {
  gmt <- read.gmt(gmt_file)
  print(names(gmt))
  print(head(gmt[[1]]))
  
  # Read with source information
  gmt_with_source <- read.gmt(gmt_file, add.source = TRUE)
  print(head(names(gmt_with_source)))
}
#>  [1] "HALLMARK_TNFA_SIGNALING_VIA_NFKB"          
#>  [2] "HALLMARK_HYPOXIA"                          
#>  [3] "HALLMARK_CHOLESTEROL_HOMEOSTASIS"          
#>  [4] "HALLMARK_MITOTIC_SPINDLE"                  
#>  [5] "HALLMARK_WNT_BETA_CATENIN_SIGNALING"       
#>  [6] "HALLMARK_TGF_BETA_SIGNALING"               
#>  [7] "HALLMARK_IL6_JAK_STAT3_SIGNALING"          
#>  [8] "HALLMARK_DNA_REPAIR"                       
#>  [9] "HALLMARK_G2M_CHECKPOINT"                   
#> [10] "HALLMARK_APOPTOSIS"                        
#> [11] "HALLMARK_NOTCH_SIGNALING"                  
#> [12] "HALLMARK_ADIPOGENESIS"                     
#> [13] "HALLMARK_ESTROGEN_RESPONSE_EARLY"          
#> [14] "HALLMARK_ESTROGEN_RESPONSE_LATE"           
#> [15] "HALLMARK_ANDROGEN_RESPONSE"                
#> [16] "HALLMARK_MYOGENESIS"                       
#> [17] "HALLMARK_PROTEIN_SECRETION"                
#> [18] "HALLMARK_INTERFERON_ALPHA_RESPONSE"        
#> [19] "HALLMARK_INTERFERON_GAMMA_RESPONSE"        
#> [20] "HALLMARK_APICAL_JUNCTION"                  
#> [21] "HALLMARK_APICAL_SURFACE"                   
#> [22] "HALLMARK_HEDGEHOG_SIGNALING"               
#> [23] "HALLMARK_COMPLEMENT"                       
#> [24] "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"        
#> [25] "HALLMARK_PI3K_AKT_MTOR_SIGNALING"          
#> [26] "HALLMARK_MTORC1_SIGNALING"                 
#> [27] "HALLMARK_E2F_TARGETS"                      
#> [28] "HALLMARK_MYC_TARGETS_V1"                   
#> [29] "HALLMARK_MYC_TARGETS_V2"                   
#> [30] "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
#> [31] "HALLMARK_INFLAMMATORY_RESPONSE"            
#> [32] "HALLMARK_XENOBIOTIC_METABOLISM"            
#> [33] "HALLMARK_FATTY_ACID_METABOLISM"            
#> [34] "HALLMARK_OXIDATIVE_PHOSPHORYLATION"        
#> [35] "HALLMARK_GLYCOLYSIS"                       
#> [36] "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY"  
#> [37] "HALLMARK_P53_PATHWAY"                      
#> [38] "HALLMARK_UV_RESPONSE_UP"                   
#> [39] "HALLMARK_UV_RESPONSE_DN"                   
#> [40] "HALLMARK_ANGIOGENESIS"                     
#> [41] "HALLMARK_HEME_METABOLISM"                  
#> [42] "HALLMARK_COAGULATION"                      
#> [43] "HALLMARK_IL2_STAT5_SIGNALING"              
#> [44] "HALLMARK_BILE_ACID_METABOLISM"             
#> [45] "HALLMARK_PEROXISOME"                       
#> [46] "HALLMARK_ALLOGRAFT_REJECTION"              
#> [47] "HALLMARK_SPERMATOGENESIS"                  
#> [48] "HALLMARK_KRAS_SIGNALING_UP"                
#> [49] "HALLMARK_KRAS_SIGNALING_DN"                
#> [50] "HALLMARK_PANCREAS_BETA_CELLS"              
#> [1] "JUNB"    "CXCL2"   "ATF3"    "NFKBIA"  "TNFAIP3" "PTGS2"  
#> [1] "HALLMARK_TNFA_SIGNALING_VIA_NFKB (http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_TNFA_SIGNALING_VIA_NFKB)"      
#> [2] "HALLMARK_HYPOXIA (http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_HYPOXIA)"                                      
#> [3] "HALLMARK_CHOLESTEROL_HOMEOSTASIS (http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_CHOLESTEROL_HOMEOSTASIS)"      
#> [4] "HALLMARK_MITOTIC_SPINDLE (http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_MITOTIC_SPINDLE)"                      
#> [5] "HALLMARK_WNT_BETA_CATENIN_SIGNALING (http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_WNT_BETA_CATENIN_SIGNALING)"
#> [6] "HALLMARK_TGF_BETA_SIGNALING (http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_TGF_BETA_SIGNALING)"                
# }
```
