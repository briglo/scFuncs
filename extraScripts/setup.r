#from http://r-pkgs.had.co.nz/intro.html and anne

install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
install.packages("rstudioapi")

#installed RStudio 1.2.1541 
devtools::install_github("r-lib/devtools")
library(devtools)
has_devel()

create()
#tweak DESCRIPTION
document()
 #GoodToGoexit :)