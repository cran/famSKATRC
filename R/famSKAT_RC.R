## Mohamad Saad, July, 2014
## mhsaad@uw.edu
##
## University of Washington, Seattle, USA
## Statistical genetics: Prof. Ellen wijsman's Group.
## Reference: Saad M and Wijsman EM, 2014. Combining family- and population-based imputation data for association analysis of rare and common variants in large pedigrees. Genet Epi (in press).


famSKAT_RC<-function(PHENO, genotypes, id, fullkins, covariates=NULL, sqrtweights_c, sqrtweights_r,  binomialimpute=FALSE, acc=NULL, maf, phi) {
MAF=colMeans(genotypes)/2
SET1 = which( MAF <= maf)   # group of SNPs with maf is smaller than maf
SET2 = which( MAF > maf)    # group of SNPs with maf is greater than maf
if (length(SET2)>1){ NAMES = colnames(genotypes[,c(SET2)])} else{NAMES=names(SET2)}
C = SET2 ; R=SET1

        phenotype = PHENO
        if(class(phenotype)!="numeric" && class(phenotype)!="integer") stop("phenotype should be a numeric vector!")
        n<-length(phenotype)
        if(is.data.frame(genotypes)) genotypes<-as.matrix(genotypes)
        if(!is.matrix(genotypes)) stop("genotypes should be a matrix!")
        if(nrow(genotypes)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
        if(!is.vector(id)) stop("id should be a vector!")
        if(any(duplicated(id))) stop("Duplicated id exists. Check your data...")
        if(length(id)!=n) stop("Number of individuals inconsistent between phenotype and id. Check your data...")
        #Add new is.bdsmatrix check here
        #fullkins <- as.bdsmatrix(fullkins)
        #if(class(fullkins) != "bdsmatrix") stop("fullkins should be a kinship matrix from makekinship function!")
        if(!is.null(covariates)) {
                covariates<-as.matrix(covariates)
                if(nrow(covariates)!=n) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")
        }
        method="Davies"
        acc = 1e-06
        missidx<-is.na(phenotype)
        if(!is.null(covariates)) {
                missidx<-missidx | apply(is.na(covariates), 1, any)
        }
        phenotype<-phenotype[!missidx]
        genotypes<-genotypes[!missidx,]
        id<-id[!missidx]
        if(!is.null(covariates)) covariates<-covariates[!missidx,]
        n<-length(phenotype)
        Z<-genotypes
        nZ<-ncol(Z)
        MAF<-colMeans(Z, na.rm = TRUE)/2
############## Weight for common variants
        if(is.function(sqrtweights_c)) {
                weights_c<-sqrtweights_c(MAF[SET2])
        } else {
                if(!is.vector(sqrtweights_c)) stop("Check the class of your variable sqrtweights_c: should be function or vector!")
                if(length(sqrtweights_c)!=nZ) stop("Number of variants inconsistent between genotypes and sqrtweights_c. Check your data...")
                weights_c<-sqrtweights_c
        }
############## Weight for rare variants
        if(is.function(sqrtweights_r)) {
                weights_r<-sqrtweights_r(MAF[SET1])
        } else {
                if(!is.vector(sqrtweights_r)) stop("Check the class of your variable sqrtweights_r: should be function or vector!")
                if(length(sqrtweights_r)!=nZ) stop("Number of variants inconsistent between genotypes and sqrtweights_r. Check your data...")
                weights_r<-sqrtweights_r
        }
######################################
        tmpidx<-!is.na(match(dimnames(fullkins)[[1]], id))
        tmpkins<-fullkins[tmpidx, tmpidx]
        tmpidx<-match(id, dimnames(tmpkins)[[1]])
        if(any(is.na(tmpidx))) stop("Some id not exist in fullkins. Check your full pedigree data...")
        kins<-as.matrix(tmpkins)[tmpidx, tmpidx]
        for(pp in 1:nZ) {
                IDX<-which(is.na(Z[,pp]))
                if(length(IDX)>0) {
                        if(binomialimpute) { # random imputation based on binomial distribution
                                Z[IDX,pp]<-rbinom(length(IDX), 2, MAF[pp])
                        } else { # default: set missing to 0
                                Z[IDX,pp]<-0
                        }
                }
        }
        Gcommon<-t(t(Z[,C])*weights_c) ;
        Grare<-t(t(Z[,R])*weights_r)
         # let lmekin estimate both variance components for the null model
                if(is.null(covariates)) {
                        X<-matrix(rep(1,n),n,1)
                        tmpdata<-data.frame(phenotype,id)
                        nullmod<-lmekin(phenotype~1 + (1|id),tmpdata,varlist=list(tmpkins))
                        #nullmod<-lme(phenotype~1 + (1|id),tmpdata,varlist=list(tmpkins))
                } else {
                        X<-cbind(rep(1,n),as.matrix(covariates))
                        tmpdata<-data.frame(phenotype,id,covariates)
                        nc<-ncol(covariates)
                        if(nc==1) {
                                exprs<-"phenotype ~ covariates"
                        } else {
                                exprs<-paste("phenotype ~", paste(names(tmpdata)[-c(1,2)],collapse=" + "))
                        }
                        nullmod<-lmekin(as.formula(exprs),tmpdata,random=~1|id,varlist=list(tmpkins))
                }
           res<-nullmod$res
           SIGMA<-as.numeric(nullmod$vcoef)*as.matrix(kins)+as.numeric(nullmod$sigma^2)*diag(n)
	### Scale the Gcommon and Grare:
	#Gcommon_t = Gcommon - X%*%solve(t(X)%*%X)%*%t(X)%*%Gcommon
  Gcommon_t = Gcommon - X%*%(1/(t(X)%*%X))%*%t(X)%*%Gcommon
	#Grare_t = Grare - X%*%solve(t(X)%*%X)%*%t(X)%*%Grare
	Grare_t = Grare - X%*%(1/(t(X)%*%X))%*%t(X)%*%Grare
	Gcommon_t = t(Gcommon_t) %*% Gcommon_t
	Grare_t = t(Grare_t) %*% Grare_t
	mean_c = sum(diag(Gcommon_t))
	mean_r = sum(diag(Grare_t))
	var_c = sum(Gcommon_t*Gcommon_t)
	var_r = sum(Grare_t*Grare_t)
	Gcommon = Gcommon/sqrt(sqrt(var_c))
	Grare = Grare/sqrt(sqrt(var_r))
  SIGMAi<-solve(SIGMA)
	#P=SIGMA - X %*% solve(t(X) %*% SIGMAi %*% X) %*% t(X)
  P=SIGMA - X %*% (1/(t(X) %*% SIGMAi %*% X)) %*% t(X)
	K1 = Grare %*% t(Grare)
        K2 = Gcommon %*% t(Gcommon)
	K3 = Z%*%t(Z)
	temp_pvalues = vector("numeric")
	for (dd in 1:length(phi))
	{
	Q = t(res) %*% SIGMAi %*% (phi[dd]*K1+(1-phi[dd])*K2) %*% SIGMAi %*% res
	#FIX
	pfamskat<-0.5
# 	eig<-eigen(  SIGMAi %*% (phi[dd]*K1+(1-phi[dd])*K2) %*% SIGMAi %*% P, symmetric=F, only.values=T)
#         evals<-eig$values
#         tmpout<-davies(Q, evals, acc=acc)
#         pfamskat<-tmpout$Qq
#         pfamskat<-min(pfamskat,1)
#         pfamskat<-max(pfamskat,0)
#         message<-tmpout$ifault
	temp_pvalues = c(temp_pvalues, pfamskat)
	}
return(temp_pvalues)
}
