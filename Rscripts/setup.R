## @knitr setup
suppressWarnings(library(plyr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape))
suppressWarnings(library(gridExtra))
suppressWarnings(library(MPTinR))
suppressWarnings(library(stringr))
suppressWarnings(library(snowfall))
options(width=60)
#listing <- function(x, options) {
#	paste("\\begin{lstlisting}[basicstyle=\\ttfamily,breaklines=true]\n", x, "\\end{lstlisting}\n", sep = "")
#}
#knit_hooks$set(source=listing, output=listing)