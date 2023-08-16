#install.packages(file.path(workdir, "fishtaco-plot-1.0.2.tar.gz"), repos = NULL, type="source")
require(FishTacoPlot); 
require(ggplot2); 
require(scales); 
require(grid)
require(ComplexHeatmap)

rm(list=ls())
workdir="D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\04FishTaco"
setwd("D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\04FishTaco")



#?MultiFunctionTaxaContributionPlots
stages = c("AD", "MCI")
#stages = c("SCS")
module_types = c("GBMs", "GMMs")

##计算分组均值。
# for (module_type in module_types){
#     ##
#     #module_type = "GMMs"
#     profiledir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\00profile"
#     module_df = read.csv(file.path(profiledir, paste0("metadata_", module_type, ".csv")), header=1, row.names=1)
#     module_df = module_df[module_df$Study %in% c("CHN", "CHN2"), ]
#     print(unique(module_df$Group))
#     module_df$Group = factor(module_df$Group, levels=c("AD", "MCI", "SCD", "SCS","NC" ))
#     module_meta = module_df[,1:14]
#     module_feature = module_df[, 15:ncol(module_df)]
#     module_feature_rel = module_feature / rowSums(module_feature)
#     module_feature_rel[is.na(module_feature_rel)] = 0
#     module_feature_scale = as.data.frame(scale(module_feature_rel))
#     #module_meta_rel = join(module)
#     module_feature_scale_mean = aggregate(module_feature_scale, by=list(module_meta$Group), FUN=mean)
#     row.names(module_feature_scale_mean) = module_feature_scale_mean[, 1]
#     diff_modules = list()
#     ###<=============================>###
#     inputdir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\04ReporterScore"
#     infile = file.path(inputdir, paste0(module_type, "out"), "merge_reporter_score.csv")
#     reporter_df = read.csv(infile, header=1)
#     diff_function_df = reporter_df[rowSums(abs(reporter_df[, -c(1,2,3)])>1.68)>0, ]
#     print(paste0(nrow(diff_function_df), " modules in ", module_type))
#     diff_modules = diff_function_df$ID
#     
#         
#     module_feature_scale_mean_filter = module_feature_scale_mean[, colnames(module_feature_scale_mean) %in% diff_modules]
#     module_feature_scale_mean_filter = t(module_feature_scale_mean_filter)
#     #img = Heatmap(module_feature_scale_mean_filter, cluster_rows=FALSE)
#     
#     pdf(paste0(module_type, "_heatmap.pdf"), width=4, height=4+0.2*nrow(module_feature_scale_mean_filter))
#     img = Heatmap(module_feature_scale_mean_filter, cluster_columns=FALSE)
#     draw(img)
#     dev.off()
# }

###<---------------------利用FishTaco------------------------>###
#stages = c("SCS")         
for (module_type in module_types){
    for(stage in stages){
        #module_type = "GBMs"
        #stage = "AD"
        print(paste0("deal with NCvs", stage, ": ", module_type))
        source(file.path(workdir, "MultiFunctionTaxaContributionPlots.R"))
        inputdir = file.path(workdir, paste0("NCvs", stage), module_type)
        inputfile = file.path(inputdir, "fishtaco_out_de_novo_inf_STAT_taxa_contributions_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab")
        if (file.exists(inputfile)){
            bars <- customMultiFunctionTaxaContributionPlots(input_dir=inputdir, 
                                                    input_prefix="fishtaco_out_de_novo_inf",
                                                    input_suffix=".tab", 
                                                    input_permutation = "single_taxa",
                                                    show_only_taxa_with_function = TRUE,
                                                    min_cont_as_separate = 0.025, 
                                                    input_taxa_taxonomy=file.path(workdir, paste0("NCvs", stage), 'species_taxon.txt'), 
                                                    sort_by="list", plot_type="bars",
                                                    #input_function_filter_list=c("ko00020", "ko00540","ko02040"), 
                                                    add_predicted_da_markers=TRUE, add_original_da_markers=TRUE)
            

            y_labels = ggplot_build(bars)$layout$panel_params[[1]]$y$get_labels()
            
            infile = file.path(inputdir, "fishtaco_out_de_novo_inf_STAT_DA_function_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab")
            diff_function_df = read.table(infile, header=1, sep="\t")
            #diff_modules = c(diff_modules, unique(diff_function_df$module))
            diff_function_df$module = paste0(diff_function_df$Function, ": ", diff_function_df$Description)
            diff_function_df = diff_function_df[match(y_labels, diff_function_df$Function), ]
            module_num = length(diff_function_df$module)
            print(length(diff_function_df$module))
            print(bars)
            
            ##<---------------------------------------->###
            bars <- bars + scale_x_continuous(breaks=seq(1, module_num, 1), 
                                        labels=diff_function_df$module) +
                guides(fill=guide_legend(ncol=6)) + ylab("Enrichment score (W)") +
                labs(title=paste0(module_type, " enriched in ", stage)) + 
                theme(plot.title=element_text(size=14,colour="black",face="plain"), 
                      axis.title.x=element_text(size=14,colour="black",face="plain"),
                      axis.text.x=element_text(size=12,colour="black",face="plain"), 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(size=12,colour="black",face="plain"), 
                      axis.ticks.y=element_blank(),
                      axis.ticks.x=element_blank(), 
                      panel.grid.major.x = element_line(colour="light gray"), 
                      panel.grid.major.y = element_line(colour="light gray"),
                      panel.grid.minor.x = element_line(colour="light gray"), 
                      panel.grid.minor.y = element_line(colour="light gray"),
                      panel.background = element_rect(fill="transparent",colour=NA), 
                      panel.border = element_rect(fill="transparent",colour="black"),
                      legend.background=element_rect(colour="black"), 
                      legend.title=element_text(size=12), 
                      legend.text=element_text(size=8,face="plain"),
                      legend.key.size=unit(0.8,"line"), 
                      legend.margin=unit(0.1,"line"), 
                      legend.position="bottom")
            if(stage=="AD"){
                width=8; height=3+0.3*module_num
            }else if(stage=="MCI"){
                width=8.5; height=3.5+0.3*module_num
            }else{
                width=8; height=3+0.3*module_num
            }
            pdf(paste0("NCvs", stage, '_', module_type, ".pdf"), width=width, height=height)
            print(bars)
            dev.off()
        }
    }
}

