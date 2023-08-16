## FishTacoPlot原理和使用说明

>   FishTacoPlot:  https://github.com/borenstein-lab/fishtaco-plot/blob/master/R/MultiFunctionTaxaContributionPlots.R

### 1. MultiFunctionTaxaContributionPlots.R

#### 1.1 读取输入文件(Line219-707)

-   Line219: 检查输入的一致性。不能够指定展示相互矛盾的内容。

-   Line244: 读取物种贡献。输入文件为:

    ```R
    ##参考案例: https://github.com/borenstein-lab/fishtaco-plot/blob/be0f9b86e6a2535b15879a8bd9859b0ba132ea62/examples/HMP_fishtaco_STAT_DA_taxa_SCORE_wilcoxon_ASSESSMENT_multi_taxa.tab
    input_file_name = paste(input_dir,"/",input_prefix,"_STAT_",input_contribution,"_SCORE_",
                              input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep="")
    
    ##行名为各个物种，列名为各个KOs, pathways或者modules。
    ```

-   Line258: 如果指定flip_case_control, 需要对物种贡献矩阵乘以-1。

-   Line270: 保存对每个代谢功能原始的正和负贡献之和。对于正贡献求和时，所有负贡献被替换为0；对于负贡献，所有正的贡献被替换为0。

-   Line285: 读取功能过滤文件input_function_filter_fil，该文件包含一个要绘图的特定功能列表文件，默认情况下时NULL。类似的是input_function_filter_list，直接传递一个列表。

-   Line311: 读取物种分类文件input_taxa_taxonomy,  该文件第一列是物种分类名称，后面的7列，分别为'kingdom','phylum','class','order','family','genus','species'。

```
s__Atopobium_vaginae	k__Bacteria	p__Actinobacteria	c__Actinobacteria	o__Coriobacteriales	f__Coriobacteriaceae	g__Atopobium	s__vaginae
s__Deinococcus_unclassified	k__Bacteria	p__Thermi	c__Deinococci	o__Deinococcales	f__Deinococcaceae	g__Deinococcus	s__unclassified
s__Peptoniphilus_unclassified	k__Bacteria	p__Firmicutes	c__Clostridia	o__Clostridiales	f__Clostridiales_Family_XI_Incertae_Sedis	g__Peptoniphilus	s__unclassified
s__Arthrobacter_unclassified	k__Bacteria	p__Actinobacteria	c__Actinobacteria	o__Actinomycetales	f__Micrococcaceae	g__Arthrobacter	s__unclassified
```

-   Line322: 读取功能计数文件input_function_count[The median count of functions across samples]，用于过滤低计数的功能, 默认是NULL。仅保留中位数Median大于min_function_counts的代谢功能。
-   Line334: 读取功能统计文件input_function_stats[Statistics about each function across the taxa],  默认为NULL。
-   Line344: 读取taxa-to-function文件，如果我们想只展示包含功能的物种，就将不包含功能的物种贡献度设置0。
-   Line365: 读取预测的丰度一致文件input_predicted_abundance_agreement[The agreement between the taxa- and metagenome based functional abundance], 默认是predicted_function_agreement。

```R
##参考文件:HMP_fishtaco_STAT_predicted_function_agreement_SCORE_wilcoxon_ASSESSMENT_multi_taxa.tab
##包括7列，分别是:KO,R^2,PearsonCorr,PearsonPval,SpearmanCorr,SpearmanPval,MeanAbsDiff,StdAbsDiff
predicted_abundance_agreement = read.table(paste(input_dir,"/",input_prefix,"_STAT_",
                                                     input_predicted_abundance_agreement,"_SCORE_",
                                                     input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                                               sep = "\t", header=TRUE, stringsAsFactors=FALSE, quote="")
```

-   Line378: 读取预测的差异丰度(DA)/转变(shift)文件input_predicted_da[The predicted (taxa-based) shift score], 默认值是predicted_DA_value。该文件包括2列，第一列是KOs, 第二列是wilcoxon打分。
-   Line392:  读取原始功能水平的差异丰度/shift数据input_original[The original (metagenome-based) calculated shift (default: "original_value")]，然后据此排序。该文件包括2列，第一列是KOs, 第二列是wilcoxon打分。【<font color="red">Line378和Line392读取的文件内容几乎相同</font>】随后按照wilcoxcon打分降序排列，默认是采用`predicted_da`方式进行排序。

```R
##参考示例文件分别为:
Line378: HMP_fishtaco_STAT_predicted_DA_value_SCORE_wilcoxon_ASSESSMENT_multi_taxa.tab
Line392: HMP_fishtaco_STAT_original_value_SCORE_wilcoxon_ASSESSMENT_multi_taxa.tab
```

-   Line445:  计算每个功能DA/shift被解释的百分比，过滤掉预测转变非常差的功能。如果预测input_predicted_abundance_agreement[<font color="red">Line365读取</font>]不为NULL, 那么要去预测的丰度一致性的KOs属于Functions。
-   Line471: 如果functions的数目大于max_function_to_show，则去掉差异丰度打分较小的functions。
-   Line500: 读取物种差异丰度数据input_taxa_da[The calculated taxonomic composition shift score (default: "DA_taxa")], 其具体内容如下:

```R
# PyCode/FiShTaCo/fishtaco/compute_contribution_to_DA.py {'functional_profile_already_corrected_with_musicc': True, 'map_function_level': 'pathway', 'taxa_abun_file': 'HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/METAPHLAN_taxa_vs_SAMPLE.tab', 'number_of_shapley_orderings_per_taxa': '5', 'max_da_functions_cases_controls': '200', 'function_abun_file': 'HMP_DATA/WGS_KO_vs_SAMPLE_MUSiCC.tab', 'number_of_permutations': '100', 'normalization_mode': 'scale_permuted', 'taxa_to_function_file': 'HMP_DATA/METAPHLAN_taxa_vs_KO.tab', 'write_log': True, 'da_result_file': 'HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/Output_METAPHLAN_REMOVE_RESIDUAL_SCALE_PERMUTED_MAPPING_IMG_MEAN/Parallel/DA_of_PATHWAY.tab', 'taxa_assessment_method': 'multi_taxa', 'max_score_cutoff': '100', 'na_rep': '0', 'case_label': '1', 'single_function_filter': 'ko00020', 'map_function_file': '/net/gs/vol1/home/ohadm/MUSiCC/Matrices/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15', 'control_label': '0', 'perform_inference_on_ko_level': False, 'score_to_compute': 'wilcoxon', 'class_file': 'HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/SAMPLE_vs_CLASS.tab', 'apply_inference': False, 'output_pref': 'HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/Output_METAPHLAN_REMOVE_RESIDUAL_SCALE_PERMUTED_MAPPING_IMG_MEAN/Parallel/ko00020/pathway_with_t2f', 'residual_mode': 'remove_residual', 'multiple_hypothesis_correction': 'Bonf', 'permutation_mode': 'blocks'}
Taxa	meanCases	meanControls	StatValue	pval	singLogP	Bonf	FDR-0.01	FDR-0.05	FDR-0.1	Metadata
s__Abiotrophia_defectiva	0.0002202359291866977	0.0046370172362582	-6.131326395190045	8.714939985429503e-10	-9.0597355992696	1.0	1.0	1.0	1.0	s__Abiotrophia_defectiva
s__Actinomyces_odontolyticus	0.016424583202763625	0.002093538658915435	11.446624943767306	2.4448004779184698e-30	29.611756578206933	1.0	1.0	1.0	1.0	s__Actinomyces_odontolyticus
s__Actinomyces_oris	0.00018681609787312662	0.004893592465039642	-10.223694928407552	1.5530571564900069e-24	-23.80881256082557	1.0	1.0	1.0	1.0	s__Actinomyces_oris
```

-   Line535: 去掉在所有功能中贡献非常低的taxa。此处的df来自于第244行读取的taxa contribution score。
-   ~~Line552: 如果绘图类型是inferred_copy_number_heatmap并且input_inference_copy_number[The inferred copy number for each function (default: "taxa_learned_copy_num")]不为空，那么将读取推断的拷贝数==>示例文件中不存在，不太常用。~~
-   ~~Line565: 读取功能的宏信息input_function_meta[Metadata on the given functions (default: NULL).], 重新命名功能的meta-names。~~
-   ~~Line617: 读取真实的功能丰度数据input_function_abundance[The metagenome-based functional abundance (default: NULL)], 即基于宏基因组的功能丰度。【<font color="red">作用不明确。</font>】~~
-   Line631: 读取预测的功能丰度数据input_predicted_function [The predicted (taxa-based) functional abundance (default: "predicted_function_abundance").], 对应的测试文件是: HMP_fishtaco_STAT_predicted_function_abundance_SCORE_wilcoxon_ASSESSMENT_multi_taxa.tab, 该文件的行为KOs, 列为样本名称。
-   Line654: 读取预测功能丰度的残差input_predicted_residual[The residual functional abundance (default: "residual_function_abundance")],  对应的测试文件为: HMP_fishtaco_STAT_residual_function_abundance_SCORE_wilcoxon_ASSESSMENT_multi_taxa.tab, 该文件的行为KOs, 列为样本名称，几乎所有的值都为0。
-   Line707： <b>FINISHED READING INPUT FILES</b>

也用于后续绘图的数据包括:

1.    Line244行读取，经过各个步骤过滤的物种贡献打分，后续代码中的df。

![image-20230815151555062](https://jialh.oss-cn-shanghai.aliyuncs.com/img2/image-20230815151555062.png)

2.   functions列表，包括如下内容:

```R
> functions
[1] "MF0051" "MF0035" "MF0054" "MF0103" "MF0040" "MF0034"
```

#### 1.2 准备用于绘图的数据(Line720~1243)

-   Line726: 将每个功能划分为正的和负的, 正负的判断是根据物种贡献打分来判断的。

-   Line750: 如果input_taxa_taxonomy[Phylogenetic assignment for each taxon (default: NULL)]不为空，对应的示例文件为HMP_TAXONOMY.tab。

    -   找到贡献最大的分类群，保持它们的原样，然后尝试找到他们的分类名称。
    -   在列举最大贡献的分类群时，如果有unknow的物种，即便其贡献比较小，也应当展示它。
    -   获得贡献的和的排序列表，如果我们希望按照最前的贡献物种排序。
    -   如果strong_contributors的长度大于0，那么对于每一个关键物种，如果"species"水平的名称不包含缺失值，如果"species"水平的名称都不等于"s__", 那么名称改为phylum+genus+species【类似地，根据界门纲目科属种的比较，可以根据各个水平，对名称进行优化】。

-   Line842: 对于所有其他物种(贡献比较少)， 将它们合并到phylum水平。如果需要将富集和删除的物种分开separate_enriched_depleted_taxa[Separate case- and control-associated taxa (default: TRUE)]。

-   Line988: 如果我们想按最强排序，那么请根据最强排序对病例和对照分别排序。经过上述处理，最终将得到两个变量，分别是df_pos和df_neg.

    ![image-20230815155912337](https://jialh.oss-cn-shanghai.aliyuncs.com/img2/image-20230815155912337.png)

    ![image-20230815155848034](https://jialh.oss-cn-shanghai.aliyuncs.com/img2/image-20230815155848034.png)

-   <b>Line1048: set the color pallete[设置调色盘]</b>

-   Line1075: 调色盘的设置是通过expand.grid()函数实现的，该函数的作用是Create a data frame from all combinations of the supplied vectors or factors. See the description of the return value for precise details of the way this is done.   cPalette是一个data.frame, 其内容如下图所示: ![image-20230815160819821](https://jialh.oss-cn-shanghai.aliyuncs.com/img2/image-20230815160819821.png) 

-   Line1135:  如果我们有物种的差异丰度信息，我们可以将这个作为另一个置于底部的bar。
-   <b>Line1243:  FINISHED PREPARING DATA FOR PLOT</b>

```R
if (verbose) {
        print("Finished Preparing data for plots")
        print(head(df_pos.m))
        print(head(df_neg.m))
        print(cPalette)
    }
```

>   <font color="red"><b>问题1 : 设置调色板时，cPalette中的h,s,v与最终的颜色之间的关系是什么？</b></font>
>
>   -   Line1205完成cPalette定义以后，直到Line1405才开始使用。
>   -   hsv(h = 1, s = 1, v = 1, alpha)的作用是通过指定色调、饱和度和数值的向量创建颜色向量。
>       -   Hue（色调、色相): 颜色圆环上所有的颜色都是光谱上的颜色，从红色开始按逆时针方向旋转，Hue=0 表示红色，Hue=120 表示绿色，Hue=240 表示蓝色等等。
>       -   Saturation（饱和度、色彩纯净度）: 饱和度为0表示纯白色。取值范围为0～100%，值越大，颜色越饱和。
>       -   Value（明度）: 竖直方向表示明度，决定颜色空间中颜色的明暗程度，明度越高，表示颜色越明亮，范围是 0-100%。明度为0表示纯黑色（此时颜色最暗）。

![img](https://jialh.oss-cn-shanghai.aliyuncs.com/img2/v2-e9f9c843e7d60e8f7aa7de1cd61d1818_720w.webp)

```R
##除去8个指定的pylum以外，还剩下3个phylum和ZOther taxa。
> print(taxa_from_other_phyla)
[1] "p__Basidiomycota g__Puccinia s__Puccinia_striiformis"                                                    
[2] "p__Candidatus_Thermoplasmatota g__Methanomassiliicoccus s__Candidatus_Methanomassiliicoccus_intestinalis"
[3] "p__Verrucomicrobia g__Akkermansia s__Akkermansia_muciniphila"                                            
[4] "ZOther taxa(127)" 

> apply(cPalette[,-1], 1, function(x) hsv(x[1],x[2],x[3]))
 [1] "#A1E6BC" "#00E65C" "#AEE6A1" "#9CE68A" "#8AE673" "#77E65C" "#65E645" "#53E62E" "#40E617" "#2EE600" "#D8A1E6" "#D697E6"
[13] "#D48DE6" "#D283E6" "#D078E6" "#CE6EE6" "#CC64E6" "#CA5AE6" "#C850E6" "#C646E6" "#C43CE6" "#C232E6" "#C028E6" "#BE1EE6"
[25] "#BC14E6" "#BA0AE6" "#B800E6" "#D8A1E6" "#E0A1E6" "#A1E6E6" "#86E6E6" "#6BE6E6" "#50E6E6" "#36E6E6" "#1BE6E6" "#00E6E6"
[37] "#A1BCE6" "#8AAEE6" "#73A1E6" "#5C93E6" "#4585E6" "#2E77E6" "#176AE6" "#005CE6" "#AEA1E6" "#8E78E6" "#6E50E6" "#4E28E6"
[49] "#2E00E6" "#E6A1E3" "#808080"
        
```

#### 1.3 开始绘图(Line1255~Line2415)

最常用的绘图代码位于Line1261~1347：

```R
BarPlot = ggplot()
if (separate_enriched_depleted_taxa) {
    if (!use_facets_for_separation) {
        no_facet_df_pos.m = df_pos.m
        no_facet_df_pos.m$x_val = 0
        for (i in 1:length(levels(df_pos.m$variable))) {
            no_facet_df_pos.m$x_val[no_facet_df_pos.m$variable == 
                                    levels(df_pos.m$variable)[i]] = i
        }
        no_facet_df_pos.m$x_val[no_facet_df_pos.m$Shift == 
                                "ENRICHED"] = no_facet_df_pos.m$x_val[no_facet_df_pos.m$Shift == 
                                                                      "ENRICHED"] + 0.175
        no_facet_df_pos.m$x_val[no_facet_df_pos.m$Shift == 
                                "DEPLETED"] = no_facet_df_pos.m$x_val[no_facet_df_pos.m$Shift == 
                                                                      "DEPLETED"] - 0.175
        no_facet_df_neg.m = df_neg.m
        no_facet_df_neg.m$x_val = 0
        for (i in 1:length(levels(df_neg.m$variable))) {
            no_facet_df_neg.m$x_val[no_facet_df_neg.m$variable == 
                                    levels(df_neg.m$variable)[i]] = i
        }
        no_facet_df_neg.m$x_val[no_facet_df_neg.m$Shift == 
                                "ENRICHED"] = no_facet_df_neg.m$x_val[no_facet_df_neg.m$Shift == 
                                                                      "ENRICHED"] + 0.175
        no_facet_df_neg.m$x_val[no_facet_df_neg.m$Shift == 
                                "DEPLETED"] = no_facet_df_neg.m$x_val[no_facet_df_neg.m$Shift == 
                                                                      "DEPLETED"] - 0.175
        if (!show_only_neg) {
            ##如果不是只展示负值。
            BarPlot = BarPlot + geom_bar(data = no_facet_df_pos.m, 
                                         aes(x = x_val, fill = Taxa, y = value), stat = "identity", 
                                         width = bar_width) + geom_bar(data = no_facet_df_pos.m, 
                                                                       aes(x = x_val, fill = Taxa, y = value), stat = "identity", 
                                                                       colour = "black", width = bar_width, show_guide = FALSE)
        }
        if (!show_only_pos) {
            ##如果不是只展示正值。
            BarPlot = BarPlot + geom_bar(data = no_facet_df_neg.m, 
                                         aes(x = x_val, fill = Taxa, y = value), stat = "identity", 
                                         width = bar_width) + geom_bar(data = no_facet_df_neg.m, 
                                                                       aes(x = x_val, fill = Taxa, y = value), stat = "identity", 
                                                                       colour = "black", width = bar_width, show_guide = FALSE)
        }
        BarPlot = BarPlot + scale_x_continuous(breaks = 1:length(levels(df_pos.m$variable)), 
                                               labels = levels(df_pos.m$variable))
    }else {
        if (!show_only_neg) {
            BarPlot = BarPlot + geom_bar(data = df_pos.m, 
                                         aes(x = Shift, fill = Taxa, y = value), stat = "identity", 
                                         width = bar_width) + geom_bar(data = df_pos.m, 
                                                                       aes(x = Shift, fill = Taxa, y = value), stat = "identity", 
                                                                       colour = "black", width = bar_width, show_guide = FALSE)
        }
        if (!show_only_pos) {
            BarPlot = BarPlot + geom_bar(data = df_neg.m, 
                                         aes(x = Shift, fill = Taxa, y = value), stat = "identity", 
                                         width = bar_width) + geom_bar(data = df_neg.m, 
                                                                       aes(x = Shift, fill = Taxa, y = value), stat = "identity", 
                                                                       colour = "black", width = bar_width, show_guide = FALSE)
        }
        BarPlot = BarPlot + geom_abline(data = df_pos.m, 
                                        intercept = 0, slope = 0, size = 1) + scale_y_continuous(breaks = round(seq(min(original_sum_of_neg_contributions$NegSum), 
                                                                                                                    max(original_sum_of_pos_contributions$PosSum), 
                                                                                                                    by = 1), 1)) + facet_wrap(~variable, ncol = 1) + 
        theme(panel.margin = unit(0, "cm"), plot.title = element_text(face = "bold", 
                                                                      size = 20), axis.title = element_text(face = "bold", 
                                                                                                            size = 20), axis.text.x = element_text(face = "bold", 
                                                                                                                                                   size = 15), panel.grid = element_blank(), 
              panel.border = element_blank(), axis.line = element_line(colour = "black", 
                                                                       size = 1))
        if (!add_facet_labels) {
            BarPlot = BarPlot + theme(strip.text = element_blank(), 
                                      strip.background = element_rect(fill = "transparent", 
                                                                      colour = NA))
        }else {
            BarPlot = BarPlot + theme(strip.text = element_text(size = 20))
        }
    }
}
```

