## -- load XML --
library(XML);

## -- load file --
doc <- xmlTreeParse(file = "F:\\Zhaolab\\MEDIC\\dataset\\PRJNA389280.xml");
doc
### -- first get document root --
sampleSet = xmlRoot(doc);

sampleSet[[1]]

df <- data.frame( sam.biosample = character(), 
                  sam.Sample_name = character(), 
                  sam.SRA  = character()
                  )
for(i in 1 : length(sampleSet) ){
  
  ## -- get biosample --
  biosamplesNodeSet <- getNodeSet( sampleSet[[i]], "//BioSample/Ids/Id[@db='BioSample']" );
  biosample <- as.character( xmlValue( biosamplesNodeSet[[1]] ));
  
  ## -- get SRA --
  sample_namesNodeSet <- getNodeSet( sampleSet[[i]], "//BioSample/Ids/Id[@db_label='Sample name']" );
  sample_name <- as.character( xmlValue( sample_namesNodeSet[[1]] ) );
  
  
  ## -- get SRA --
  IdsNodeSet <- getNodeSet( sampleSet[[i]], "//BioSample/Ids/Id[@db='SRA']" );
  sra <- as.character( xmlValue( IdsNodeSet[[1]] ) );
  
  
  data_return <-  data.frame( sam.biosample = biosample, 
                              sam.Sample_name = sample_name, 
                              sam.SRA = sra ) 
  df <- rbind(df , data_return)
} 
write.csv(df ,file = "F:\\Zhaolab\\MEDIC\\dataset\\PRJNA389280.csv",row.names = F)
