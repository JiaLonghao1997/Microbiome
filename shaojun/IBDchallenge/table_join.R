sample_name <- read.csv(file = "F:\\Zhaolab\\MEDIC\\dataset\\PRJNA398089.csv", header = T)
metadata <- read.csv(file = "F:\\Zhaolab\\MEDIC\\dataset\\hmp2_metadata_deduplicated.csv", header = T)
head(sample_name)
head(metadata)
m1 <- merge(sample_name, metadata, by.x = "Sample_name", by.y = "Project")
write.csv(m1 ,file = "F:\\Zhaolab\\MEDIC\\dataset\\PRJNA398089_metadata.csv",row.names = F)
