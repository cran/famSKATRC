\name{sample.ped.geno}
\alias{sample.ped.geno}
\docType{data}
\title{
Sample Pedigree Genotype Data
}
\description{
A sample pedigree file with SNP data, already processed by the process_data()
function.  The raw data is also included in the package.
}
\usage{data("sample.ped.geno")}
\format{
  A data frame with 20 observations on the following 36 variables.
  \describe{
    \item{\code{FID}}{a numeric vector}
    \item{\code{IID}}{a character vector}
    \item{\code{FA}}{a character vector}
    \item{\code{MO}}{a character vector}
    \item{\code{SEX}}{a numeric vector}
    \item{\code{Phenotype}}{a numeric vector}
    \item{\code{rs1}}{a numeric vector}
    \item{\code{rs2}}{a numeric vector}
    \item{\code{rs3}}{a numeric vector}
    \item{\code{rs4}}{a numeric vector}
    \item{\code{rs5}}{a numeric vector}
    \item{\code{rs6}}{a numeric vector}
    \item{\code{rs7}}{a numeric vector}
    \item{\code{rs8}}{a numeric vector}
    \item{\code{rs9}}{a numeric vector}
    \item{\code{rs10}}{a numeric vector}
    \item{\code{rs11}}{a numeric vector}
    \item{\code{rs12}}{a numeric vector}
    \item{\code{rs13}}{a numeric vector}
    \item{\code{rs14}}{a numeric vector}
    \item{\code{rs15}}{a numeric vector}
    \item{\code{rs16}}{a numeric vector}
    \item{\code{rs17}}{a numeric vector}
    \item{\code{rs18}}{a numeric vector}
    \item{\code{rs19}}{a numeric vector}
    \item{\code{rs20}}{a numeric vector}
    \item{\code{rs21}}{a numeric vector}
    \item{\code{rs22}}{a numeric vector}
    \item{\code{rs23}}{a numeric vector}
    \item{\code{rs24}}{a numeric vector}
    \item{\code{rs25}}{a numeric vector}
    \item{\code{rs26}}{a numeric vector}
    \item{\code{rs27}}{a numeric vector}
    \item{\code{rs28}}{a numeric vector}
    \item{\code{rs29}}{a numeric vector}
    \item{\code{rs30}}{a numeric vector}
  }
}
\details{
A pedigree file with 20 individuals from two families.  The first column is a family ID,
the second an individual ID, the 3rd gives their father, the 4th their mother, and the
5th their sex.  The 6th column gives their phenotype and columns 7 through 36 give their
genotype, 30 SNP loci.
}
\source{
This is synthetically generated data.
}
\examples{
data(sample.ped.geno)
library(kinship2)
sample.ped.geno$FA[sample.ped.geno$FA == "0"] <- NA
sample.ped.geno$MO[sample.ped.geno$MO == "0"] <- NA
sample.ped.geno$Phenotype[sample.ped.geno$Phenotype >= 0] <- 1
sample.ped.geno$Phenotype[sample.ped.geno$Phenotype < 0] <- 0
ped.list <- pedigree(id = sample.ped.geno$IID, dadid = sample.ped.geno$FA,
                     momid = sample.ped.geno$MO, sex = sample.ped.geno$SEX,
                     famid = sample.ped.geno$FID,
                     affected = sample.ped.geno$Phenotype)
ped1 <- ped.list[1]
ped2 <- ped.list[2]
plot(ped1)
plot(ped2)
}
\keyword{datasets}
