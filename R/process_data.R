process_data<-function(Data = read.table(system.file("extdata", "data",
                                                package = "famSKATRC"), header = TRUE))
{
  Data[ , "IID"] = paste(Data[ , "FID"]  , Data[ , "IID"]  ,sep=".")
  Data[Data[,"FA"]!=0 , "FA"] = paste(Data[Data[,"FA"]!=0 , "FID"], Data[Data[,"FA"]!=0,
                                                                        "FA"]  ,sep=".")
  Data[Data[,"FA"]!=0 , "MO"] = paste(Data[Data[,"FA"]!=0 , "FID"], Data[Data[,"FA"]!=0,
                                                                        "MO"]  ,sep=".")
  return(Data)
}
