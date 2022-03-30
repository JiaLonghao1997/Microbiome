## -- load XML --
library(XML);

## -- load file --
doc <- xmlTreeParse(file = "F:\\Zhaolab\\MEDIC\\dataset\\PRJNA398089.xml");
#doc
### -- first get document root --
sampleSet = xmlRoot(doc);

length(sampleSet)

df <- data.frame( sam.biosample = character(), 
                  sam.Sample_name = character(), 
                  sam.SRA  = character()
                  )
for(i in 1 : length(sampleSet) ){
  
  ## -- get biosample --
  biosamplesNodeSet <- getNodeSet( sampleSet[[i]], "//BioSample/Ids/Id[@db='BioSample']" );
  as.character( xmlValue( biosamplesNodeSet[[1]] ))
  if(length(biosamplesNodeSet) == 0){
    biosample = 'NA';
  }else{
    biosample <- as.character( xmlValue( biosamplesNodeSet[[1]] ));
  }
  
  
  ## -- get SRA --
  sample_namesNodeSet <- getNodeSet( sampleSet[[i]], "//BioSample/Ids/Id[@db_label='Sample name']" );
  if(length(sample_namesNodeSet) == 0){
    sample_name = 'NA';
  }
  else{
    sample_name <- as.character( xmlValue( sample_namesNodeSet[[1]] ) );
  }
  
  
  ## -- get SRA --
  IdsNodeSet <- getNodeSet( sampleSet[[i]], "//BioSample/Ids/Id[@db='SRA']" );
  if(length(IdsNodeSet) == 0){
    sra = 'NA';
  }
  else{
    sra <- as.character( xmlValue( IdsNodeSet[[1]] ) );
  }
  
  
  data_return <-  data.frame( sam.biosample = biosample, 
                              sam.Sample_name = sample_name, 
                              sam.SRA = sra ) 
  df <- rbind(df , data_return)
} 
write.csv(df ,file = "F:\\Zhaolab\\MEDIC\\dataset\\PRJNA398089.csv",row.names = F)
