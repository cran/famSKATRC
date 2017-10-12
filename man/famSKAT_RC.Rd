\name{famSKAT_RC}
\alias{famSKAT_RC}
\title{
famSKAT_RC
}
\description{
FamSKAT-RC is a family-based association kernel test for both rare and common variants.
}
\usage{
famSKAT_RC (PHENO, genotypes, id, fullkins, covariates=NULL, sqrtweights_c,
            sqrtweights_r, binomialimpute=FALSE, acc=NULL, maf, phi)
}
\arguments{
  \item{PHENO}{
  The vector of the phenotype values. The missing values must be represented by NA.
  }
  \item{genotypes}{
  The genotype matrix. Its size should be N * P, where N is the number of individuals and P is the number of SNPs.
  }
  \item{id}{
  The vector of individual IDs to be included in the test. These IDs are present in the kinship matrix.
  }
  \item{fullkins}{
  The full kinship matrix that contains all individuals.
  }
  \item{covariates}{
  A matrix of possible covariates.
  }
  \item{sqrtweights_c}{
  The weight function to be assigned for common variants. An example is: function(maf) ifelse(maf>0, dbeta(maf,0.5,0.5), 0).
  }
  \item{sqrtweights_r}{
  The weight function to be assigned for rare variants. An example is: function(maf) ifelse(maf>0, dbeta(maf,1,25), 0).
  }
  \item{binomialimpute}{
  If TRUE, then impute missing genotypes using a binomial distribution (e.g. rbinom(10, 2, MAF), if 10 genotypes are missing).
  }
  \item{acc}{
   The accuracy used in the Davies approximation. For example, acc= 1e-06.
  }
  \item{maf}{
  A MAF threshold used to define rare and common variants (e.g. maf=0.01).
  }
  \item{phi}{
  This parameter indicates the contribution portion of rare variants. For example, a value of phi=0.5 means that the contribution of rare and common variants to the test is equaland a value of phi=1 means that only rare variants contribute to the test. You can vary the phi values as you wish and you can also provide a grid of different values (e.g. phi = c(0,0.2,0.5,0.9)). In this case, four p-values will be obtained, one for each value of phi.
  }
}
\details{
FamSKAT-RC is a family-based association kernel test for both rare and common variants. This test is general and several special cases are known as other methods: famSKAT, which only focuses on rare variants in family-based data, SKAT, which focuses on rare variants in population-based data (unrelated individuals), and SKAT-RC, which focuses on both rare and common variants in population-based data. When one applies famSKAT-RC and sets the value of phi to 1, famSKAT-RC becomes famSKAT. When one applies famSKAT-RC and set the value of phi to 1 and the kinship matrix to the identity matrix, famSKAT-RC becomes SKAT. When one applies famSKAT-RC and set the kinship matrix (fullkins) to the identity matrix (and phi is not equal to 1), famSKAT-RC becomes SKAT-RC.
}

\value{
  The test statistic p-value.
  \item{temp_pvalues}{Some Description}
}

\examples{
library(kinship2)
sample.ped.geno <- process_data()
KIN = kinship(sample.ped.geno$IID, sample.ped.geno$FA, sample.ped.geno$MO)
IID = sample.ped.geno$IID
wuweights_r <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)
wuweights_c <- function(maf) ifelse(maf>0, dbeta(maf,0.5,0.5), 0)
P_VALUES <- famSKAT_RC(PHENO=sample.ped.geno[,"Phenotype"],genotypes=as.matrix(
            sample.ped.geno[,7:ncol(sample.ped.geno)]), binomialimpute=TRUE,
            id=IID,fullkins=KIN,maf=0.05, sqrtweights_c=wuweights_c,
            sqrtweights_r=wuweights_r, phi = c(0,0.2,0.5,0.9))
print(P_VALUES)
}
