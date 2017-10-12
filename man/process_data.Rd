\name{process_data}
\alias{process_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process Data
}
\description{
Processes the raw data included to produce data identical to the (also included)
sample.ped.geno data set.  This can be used for processing your own pedigrees with
SNP data.
}
\usage{
process_data(Data = read.table(system.file("extdata", "data",
                                            package = "famSKATRC"), header = TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
A string, the path to the location of the data file you are processing, formatted as
the included example, which can be loaded with:
\code{read.table(system.file("extdata", "data", package = "famSKATRC"), header = TRUE)}
and can be found in your filesystem at: system.file("extdata", "data",
                                                                 package = "famSKATRC")
}
}
\value{
Returns the data frame with completed preprocessing changes for famSKATRC.
Mainly reworking IDs so there are not duplicates.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{famSKATRC}}
}
\examples{
sample.ped.geno <- process_data()
## The function is currently defined as
function(Data = read.table(system.file("extdata", "data",
                                                package = "famSKATRC"), header = TRUE))
{
  Data[ , "IID"] = paste(Data[ , "FID"]  , Data[ , "IID"]  ,sep=".")
  Data[Data[,"FA"]!=0 , "FA"] = paste(Data[Data[,"FA"]!=0 , "FID"], Data[Data[,"FA"]!=0,
                                                                        "FA"]  ,sep=".")
  Data[Data[,"FA"]!=0 , "MO"] = paste(Data[Data[,"FA"]!=0 , "FID"], Data[Data[,"FA"]!=0,
                                                                        "MO"]  ,sep=".")
  return(Data)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pedigree }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~data }% __ONLY ONE__ keyword per line
