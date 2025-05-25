## Followings steps of "R Packages" book
##

library(usethis)
library(devtools)

setwd("~/Playground/projects/plaid")
create_package("~/Playground/projects/plaid")
use_git()

use_r("plaid")
use_gpl3_license()

usethis::use_vignette("plaid")

use_testthat()
use_test("plaid")

## List package dependencies
use_package("Matrix")
use_package("Rfast")
use_package("matrixStats")
use_package("sparseMatrixStats")
use_package("metap")

## For vignette
#use_package("Seurat", type = "Suggests")
#use_package("SeuratData", type = "Suggests")

use_readme_md()

check()
install()

## Reload
load_all()




