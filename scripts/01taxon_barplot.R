library(ggplot2)
library(doBy)
library(RColorBrewer)
library(scales)
library(plyr)
library(reshape)

inputdir="F:\\Zhaolab2020\\virome\\TGScompare\\Results0228"
taxon<-read.table(file.path(inputdir, "merged_abundance_table_species.txt"), sep="\t", header=TRUE, row.names = 1)
print(taxon[1:5, ])
taxon<-t(taxon)
taxon_rel<-taxon/rowSums(taxon)
taxon_rel<-as.data.frame(t(taxon_rel))
print(taxon_rel[1:5, ])                 
print(colSums(taxon_rel))


taxon_rel$sum <- rowSums(taxon_rel)
# taxon_rel1 <- cbind(taxon_rel, sum=sum)
# taxon_rel <- taxon_rel1
print(taxon_rel[1:5, ])
taxon_rel <- taxon_rel[order(taxon_rel$sum, decreasing = TRUE), ]
print(taxon_rel[1:5, ])
taxon_top10 <- taxon_rel[1:10, -ncol(taxon_rel) ]
taxon_top10['others', ] <- 1 - colSums(taxon_top10)

#names(taxon_top10) <- c("Taxonomy", "sample", "abundance") 

print(taxon_top10)


taxon_top10$Taxonomy <- factor(rownames(taxon_top10), levels = rev(rownames(taxon_top10)))
taxon_top10 <- melt(taxon_top10, id = 'Taxonomy')
names(taxon_top10) <- c("Taxonomy", "sample", "abundance") 
print(taxon_top10)

group <- read.table(file.path(inputdir,'sample_group.csv'), sep = ',', header=TRUE, stringsAsFactors = FALSE)
print(group)
taxon_top10 <- merge(taxon_top10, group, by.x = 'sample', by.y='Sample')

print(taxon_top10)

taxon_top10$Taxonomy = factor(taxon_top10$Taxonomy)
taxon_top10$Taxonomy = factor(taxon_top10$Taxonomy, levels = unique(taxon_top10$Taxonomy[order(taxon_top10$abundance, decreasing = T)]))



colors = c(brewer.pal(11, 'Set3'))
print(colors)

png(filename = file.path(inputdir, "taxon_barplot_species.png"),width = 1200,height = 900,res = 300)
p <- ggplot(taxon_top10, aes(sample, abundance, fill = Taxonomy)) +
  geom_col(position = 'stack', width = 0.6) + scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  facet_wrap(~group, scales = 'free_x') + 
  labs(x = 'Samples', y = 'Relative abundance', title = "Relative abundance in species") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一???
print(p)
dev.off()

