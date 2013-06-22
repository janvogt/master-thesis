## @knitr setup
suppressMessages(library(plyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(gridExtra))
suppressMessages(library(MPTinR))
suppressMessages(library(stringr))
options(width=60)
#listing <- function(x, options) {
#	paste("\\begin{lstlisting}[basicstyle=\\ttfamily,breaklines=true]\n", x, "\\end{lstlisting}\n", sep = "")
#}
#knit_hooks$set(source=listing, output=listing)