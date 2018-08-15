################################################################################################
#
#           5. Contingency table analysis                      
#
#
################################################################################################

# This is using data from the contingency tables (Tables 2 and 3 in the manuscript)
# This script was lightly modified from a script provided by David W. Crowder

library(vcd)
library(DescTools)

# Set up overall matrix for Tacotalpa
# Include all data regardless of subset (Upreg, Fixed, etc.)
everything_taco<-matrix(c(34,326,304,5697),nrow=2)

# Fisher exact test on overall data for Tacotalpa
fisher.test(everything_taco)

# Enter data with the six groupings for Tacotalpa (don't include everything due to nesting)
Input =(
"                  Group Nonsyn Syn Upreg Downreg Fixed Notfixed
Type     Type2
DE       Under           13     21   19   15      17    17
Not             347    339  153  173     343   343
NotDE    Under           77     227  319  323     150   154
Not             5924   5774 5870 5850    5851  5847
")

# Read this in as a table
Taco = as.table(read.ftable(textConnection(Input)))

# View the table
ftable(Taco) 

# Generate individual Fisher exact tests for each grouping
n = dim(Taco)[3]

for(i in 1:n){
  Name = dimnames(Taco)[3]$Group[i]
  P.value = fisher.test(Taco[,,i])$p.value
  cat(Name, "\n")
  cat("Fisher test p-value: ", P.value, "\n")
  cat("\n")
}

# Get log odds ratios for everything and each grouping
oddsratio(everything_taco)
oddsratio(Taco)

# Mantel-Haenzel test of overall differences, considering stratification
mantelhaen.test(Taco)

# Woolf Test, and Breslow-Day test, for differences among the strata
woolf_test(Taco)
BreslowDayTest(Taco)

# Set up overall matrix for Puyacatengo
everything_puya<-matrix(c(18,256,111,5976),nrow=2)

# Fisher exact test on overall data for Puyacatengo
fisher.test(everything_puya)

# Enter data with the six groupings for Puyacatengo (don't include everything due to nesting)
Input2 =(
"                  Group Nonsyn Syn  Upreg Downreg Fixed Notfixed
Type     Type2
DE       Under           9      9    12    6       2     16
Not                      265    265  117   139     272   258
NotDE    Under           37     74   117   123     9     100
Not                      6050   6013 6115  6093    6078  5987
")

# Read this in as a table
Puya = as.table(read.ftable(textConnection(Input2)))

# View the table
ftable(Puya) 

# Generate individual Fisher exact tests for each grouping
n = dim(Puya)[3]

for(i in 1:n){
  Name = dimnames(Puya)[3]$Group[i]
  P.value = fisher.test(Puya[,,i])$p.value
  cat(Name, "\n")
  cat("Fisher test p-value: ", P.value, "\n")
  cat("\n")
}

# Mantel-Haenzel test of overall differences, considering stratification
mantelhaen.test(Puya)

# Woolf Test, and Breslow-Day test, for differences among the strata
woolf_test(Puya)
BreslowDayTest(Puya)