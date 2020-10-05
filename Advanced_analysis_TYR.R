#####Advanced analysis and visualization
####Part I: More details about the top drivers####
ms_tab <- ms_tab_combined_edit ## get the master table data frame
ms_tab <- ms_tab[which(ms_tab$Size_Isabel>=30 & ms_tab$Size_Isabel <=1000),] 
###filtered ms_tabl by individual z value > 1.96, 
#that is make sure each drivers are up-regulated and significant in individual assay
ms_tab_TYR<-ms_tab[which(ms_tab$'Z-statistics.TYR.vs.Others_DA_Isabel'>1.96 & ms_tab$'Z-statistics.TYR.Vs.others_DA_Xu'>1.96),]
#ranked ms_tab_TYR
ms_tab_TYR<-ms_tab_TYR[order(ms_tab_TYR$'Z-statistics_combineDataDA_TYR_vs_others',decreasing=T),]
out2excel(ms_tab_TYR,"/Volumes/XuYang/ATRT/analysis/ATRT_combined_Isabel_Xu/ATRT_combined_2020-09-22/DATA/TYR_driver_table.xlsx")
ms_tab_TYR<-ms_tab_TYR[!duplicated(ms_tab_TYR$gene_label),]
sig_driver<-ms_tab_TYR[1:800,]
#for sig_gene
ms_tab_TYR_gene<-ms_tab[which(ms_tab$'Z-statistics_TYR_vs_others_DE_Isabel'>1.96 & ms_tab$'Z-statistics_TYR_vs_others_DE_Xu'>1.96),]
ms_tab_TYR_gene<-ms_tab_TYR_gene[order(sig_gene$'Z-statistics_combinateDataDE_TYR_vs_others',decreasing=T),]
ms_tab_TYR_gene<-ms_tab_TYR_gene[!duplicated(ms_tab_TYR_gene$gene_label),]
sig_gene<-ms_tab_TYR_gene[1:800,]

comp_name <- 'TYR.Vs.others' ## get the comparison name
#1e-18 538
sig_driver <- draw.volcanoPlot(dat=ms_tab_TYR,label_col='gene_label',logFC_col="logFC.TYR.Vs.others_DA_Xu",
                               Pv_col="P.Value_combineDataDA_TYR_vs_others",logFC_thre=0,Pv_thre=1e-18,
                               main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                               pdf_file=sprintf('%s/vocalno_label_DA.pdf',analysis.par$out.dir.PLOT),label_cex = 1)
sig_driver<-ms_tab_TYR[1:800,]

#1e-13
sig_gene <- draw.volcanoPlot(dat=ms_tab_TYR,label_col='geneSymbol',logFC_col="logFC_combinateDataDE_TYR_vs_others",
                             Pv_col="P.Value_combinateDataDE_TYR_vs_others",logFC_thre=1,Pv_thre=1e-3,
                             main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                             pdf_file=sprintf('%s/vocalno_label_DE.pdf',analysis.par$out.dir.PLOT),label_cex = 1)
sig_gene<-ms_tab_TYR_gene[1:800,]

####QI.2: How to interpret the significance of top DA drivers ?####
# Get the DE data frame of target genes
DE <- analysis.par$DE[[comp_name]]
driver_list_full <- rownames(sig_driver) # The rownames is the originalID_label
driver_list<-sig_driver[driver_list_full,"originalID"]


draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab_TYR[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab_TYR[driver_list,"Z-statistics_combineDataDA_TYR_vs_others"],
                 driver_DE_Z=ms_tab_TYR[driver_list,"Z-statistics_combinateDataDE_TYR_vs_others"],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in TYR',right_annotation = 'high in others',
                 main=comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_TYR.pdf',analysis.par$out.dir.PLOT))

####QI.3: What is the expression/activity pattern of these top DA drivers across sample subtypes?####    

exp_mat <- exprs(analysis.par$cal.eset) # expression matrix, the rownames must be the originalID
ac_mat <- exprs(analysis.par$merge.ac.eset) # activity matrix, the rownames must be the originalID_label
phe_info <- pData(analysis.par$cal.eset) # phenotype data frame
#DA of DA
phe_info2<-phe_info[order(phe_info$Subgroup_new),]
ac_mat2<-ac_mat[,rownames(phe_info2)]
driver<-driver_list_full[1:50]
draw.heatmap(mat=ac_mat2,use_genes=driver,use_gene_label=sig_driver[driver,'gene_label'],
             use_samples=colnames(ac_mat2),
             phenotype_info=phe_info2,use_phe=c("Subgroup_new"),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=F,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_DAbyDA.pdf',analysis.par$out.dir.PLOT))

#DE of DA
exp_mat2<-exp_mat[,rownames(phe_info2)]
driver<-driver_list_full[1:50]
draw.heatmap(mat=exp_mat2,use_genes=ms_tab_TYR[driver,'geneSymbol'],use_gene_label=ms_tab_TYR[driver,'gene_label'],
             use_samples=colnames(ac_mat2),
             phenotype_info=phe_info2,use_phe=c("Subgroup_new"),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=F,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_DAbyDE.pdf',analysis.par$out.dir.PLOT))

####QI.4: What are the biological functions of these top DA drivers ?####
# Download gene sets from MSigDB and save as RData, creat a global variable all_gs2gene
gs.preload(use_spe='Homo sapiens',update=FALSE)
print(all_gs2gene_info)
# Gene Set Enrichment Analysis
res_up <- funcEnrich.Fisher(input_list=sig_driver[driver_list_full,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs=c('H','CP:REACTOME','BP',"CP:KEGG","CGP"), Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
res_up_H <- funcEnrich.Fisher(input_list=sig_driver[driver_list_full,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs='H', Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
res_up_kegg <- funcEnrich.Fisher(input_list=sig_driver[driver_list_full,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs='CP:KEGG', Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
# Save gene set enrichment analysis results as EXCEl
out2excel(res_up,out.xlsx=sprintf('%s/fisher_res_TYR.xlsx',analysis.par$out.dir.PLOT))

# Gene set enrichment analysis Barplot
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_up.pdf',analysis.par$out.dir.PLOT))  # with gene name

draw.funcEnrich.bar(funcEnrich_res= res_up_H,top_number=30,main='Function Enrichment for Top drivers',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_up_H.pdf',analysis.par$out.dir.PLOT))  # with gene name
draw.funcEnrich.bar(funcEnrich_res= res_up_kegg,top_number=30,main='Function Enrichment for Top drivers',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_up_kegg.pdf',analysis.par$out.dir.PLOT))  # with gene name


# Gene set enrichment analysis Function Cluster Plot
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_cluster_DA_up.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95,Pv_thre=0.05)
draw.funcEnrich.cluster(funcEnrich_res= res_up_H,top_number=30,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_cluster_DA_up_H.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95,Pv_thre=0.05)
draw.funcEnrich.cluster(funcEnrich_res= res_up_kegg,top_number=30,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_cluster_DA_up_kegg.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95,Pv_thre=0.05)

# Gene Set Enrichment Analysis of DE
DE_list_up <- rownames(sig_gene) # up
DE_res_up <- funcEnrich.Fisher(input_list=ms_tab_TYR[DE_list_up,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs=c('H','CP:REACTOME','BP',"CP:KEGG","CGP"), Pv_thre=1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
DE_res_up_H <- funcEnrich.Fisher(input_list=ms_tab_TYR[DE_list_up,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs='H', Pv_thre=1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
DE_res_up_kegg <- funcEnrich.Fisher(input_list=ms_tab_TYR[DE_list_up,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs='CP:KEGG', Pv_thre=1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)


# Save gene set enrichment analysis results as EXCEl
out2excel(list(up=DE_res_up),out.xlsx=sprintf('%s/fisher_DE.xlsx',analysis.par$out.dir.PLOT))

# Gene set enrichment analysis Barplot
#draw.funcEnrich.bar(funcEnrich_res= DE_res_up,top_number=30,main='Function Enrichment for Top drivers',pdf_file=sprintf('%s/funcEnrich_bar_nogene.pdf',analysis.par$out.dir.PLOT)) #without gene name
draw.funcEnrich.bar(funcEnrich_res= DE_res_up,top_number=30,main='Function Enrichment for Top DE',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_DE_up.pdf',analysis.par$out.dir.PLOT))  # with gene name

draw.funcEnrich.bar(funcEnrich_res= DE_res_up_H,top_number=30,main='Function Enrichment for Top DE',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_DE_up_H.pdf',analysis.par$out.dir.PLOT))  # with gene name

draw.funcEnrich.bar(funcEnrich_res= DE_res_up_kegg,top_number=30,main='Function Enrichment for Top DE',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_DE_up_kegg.pdf',analysis.par$out.dir.PLOT))  # with gene name

# Gene set enrichment analysis Function Cluster Plot
draw.funcEnrich.cluster(funcEnrich_res= DE_res_up,top_number=30,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_cluster_DE_up.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95)

###QI.5:the biological functions of the target genes of these top DA drivers#########
# Get ID conversion table
# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)

# Get all comparison names
transfer_tab<-analysis.par$transfer_tab

# Bubble Plot to show target genes enriched biological functions

draw.bubblePlot(driver_list= driver_list,show_label=ms_tab_TYR[driver_list_full,'gene_label'],
                Z_val=ms_tab_TYR[driver_list_full,"Z-statistics_combineDataDA_TYR_vs_others"],
                driver_type=NULL,
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=5,max_gs_size=500,use_gs=c('H','CP:REACTOME','BP',"CP:KEGG","CGP"),
                top_geneset_number=15,top_driver_number=15,
                pdf_file = sprintf('%s/bubblePlot.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')

##Hallmark
draw.bubblePlot(driver_list= driver_list,show_label=ms_tab_TYR[driver_list,'gene_label'],
                Z_val=ms_tab_TYR[driver_list,"Z-statistics_combineDataDA_TYR_vs_others"],
                driver_type=NULL,
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=5,max_gs_size=500,use_gs='H',
                top_geneset_number=15,top_driver_number=15,
                pdf_file = sprintf('%s/bubblePlot_H.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')
##Hallmark
draw.bubblePlot(driver_list= driver_list,show_label=ms_tab_TYR[driver_list,'gene_label'],
                Z_val=ms_tab_TYR[driver_list,"Z-statistics_combineDataDA_TYR_vs_others"],
                driver_type=NULL,
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=5,max_gs_size=500,use_gs='CP:KEGG',
                top_geneset_number=15,top_driver_number=15,
                pdf_file = sprintf('%s/bubblePlot_kegg.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')


#####GS analysis by DA###### this part should be the same as myself data, bc it is only need ac_matix
# with larger target size for duplicate drivers
#first filtered ac_mat to make one gene only have one probe.
#using ac_mat_iqr_final and
ms_tab_xu<-analysis.par[["final_ms_tab"]]
ac_mat_iqr<-as.data.frame(ac_mat)
iqr<-data.matrix(apply(ac_mat,1,IQR))
colnames(iqr)<-"IQR"
ac_mat_iqr$IQR<-iqr[,1]
ac_mat_iqr<-ac_mat_iqr[order(ac_mat_iqr$IQR,decreasing = T),]  #23078
ac_mat_iqr$external_gene_name<-ms_tab_xu[rownames(ac_mat_iqr),"geneSymbol"]
ac_mat_iqr_final<-ac_mat_iqr[!duplicated(ac_mat_iqr$external_gene_name),] #10331
ac_mat_iqr_final$PGnames<-rownames(ac_mat_iqr_final)
rownames(ac_mat_iqr_final)<-ac_mat_iqr_final$external_gene_name
ac_mat_gene<-ac_mat_iqr_final[,1:21]

use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,use_gs=c('H','CP:REACTOME','BP',"CP:KEGG","CGP"))
use_gs2gene_H <- merge_gs(all_gs2gene=all_gs2gene,use_gs='H')
use_gs2gene_kegg <- merge_gs(all_gs2gene=all_gs2gene,use_gs='CP:KEGG')

driver_ac_gs <- cal.Activity.GS(use_gs2gene = use_gs2gene,
                                cal_mat = ac_mat_gene)    
driver_ac_gs_H<- cal.Activity.GS(use_gs2gene = use_gs2gene_H,
                                 cal_mat = ac_mat_gene)
driver_ac_gs_kegg<- cal.Activity.GS(use_gs2gene = use_gs2gene_kegg,
                                 cal_mat = ac_mat_gene)

# Calculate DA 
phe_info <- pData(analysis.par$cal.eset)

G1  <- rownames(phe_info)[which(phe_info$`Subgroup_new`=='TYR')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup_new`!='TYR')] # Combine other groups as the Control group

DA_gs_bid_TYR <- getDE.BID.2G(eset=generate.eset(driver_ac_gs),G1=G1,G0=G0,G1_name='TYR',G0_name='others')
DA_gs_bid_TYR_H<- getDE.BID.2G(eset=generate.eset(driver_ac_gs_H),G1=G1,G0=G0,G1_name='TYR',G0_name='others')
DA_gs_bid_TYR_kegg<- getDE.BID.2G(eset=generate.eset(driver_ac_gs_kegg),G1=G1,G0=G0,G1_name='TYR',G0_name='others')

DA_gs_bid_TYR_up<-DA_gs_bid_TYR[which(DA_gs_bid_TYR[,4]>0),]
DA_gs_bid_TYR_down<-DA_gs_bid_TYR[which(DA_gs_bid_TYR[,4]<0),]

out2excel(list(up=DA_gs_bid_TYR_up,down=DA_gs_bid_TYR_down),out.xlsx=sprintf('%s/gs_DA_TYR.xlsx',analysis.par$out.dir.PLOT))

# Draw volcano plot for top significant gene sets
DA_gs_TYR <- draw.volcanoPlot(dat= DA_gs_bid_TYR,label_col='ID',logFC_col='logFC',
                              Pv_col='P.Value',logFC_thre=0,Pv_thre=1e-3,
                              main='Volcano Plot for gene sets',show_label=TRUE,label_type = 'distribute',label_cex = 0.5,
                              pdf_file=sprintf('%s/vocalno_GS_DA_TYR.pdf',analysis.par$out.dir.PLOT)) 


DA_gs_TYR_H <- draw.volcanoPlot(dat= DA_gs_bid_TYR_H,label_col='ID',logFC_col='logFC',
                                Pv_col='P.Value',logFC_thre=0,Pv_thre=1e-3,
                                main='Volcano Plot for gene sets',show_label=TRUE,label_type = 'distribute',label_cex = 0.5,
                                pdf_file=sprintf('%s/vocalno_GS_DA_TYR_H.pdf',analysis.par$out.dir.PLOT)) 

DA_gs_TYR_kegg <- draw.volcanoPlot(dat= DA_gs_bid_TYR_kegg,label_col='ID',logFC_col='logFC',
                                Pv_col='P.Value',logFC_thre=0,Pv_thre=1e-3,
                                main='Volcano Plot for gene sets',show_label=TRUE,label_type = 'distribute',label_cex = 0.5,
                                pdf_file=sprintf('%s/vocalno_GS_DA_TYR_kegg.pdf',analysis.par$out.dir.PLOT)) 


DA_gs_TYR_up<-DA_gs_TYR[which(DA_gs_TYR[,2]>0),]
DA_gs_TYR_down<-DA_gs_TYR[which(DA_gs_TYR[,2]<0),]
DA_gs_TYR_H_up<-DA_gs_TYR_H[which(DA_gs_TYR_H[,2]>0),]
DA_gs_TYR_kegg_up<-DA_gs_TYR_kegg[which(DA_gs_TYR_kegg[,2]>0),]

# Draw GSEA plot for top significant gene sets
DA <- analysis.par$DA[[comp_name]]
DA_name <- ac_mat_iqr_final$PGnames
DA<-DA[DA_name,]
rownames(DA)<-ac_mat_iqr_final$external_gene_name
DA$ID<-rownames(DA)

# up-regulated
draw.GSEA.NetBID.GS(DE=DA,name_col='ID',profile_col='logFC',profile_trend='pos2neg',
                    sig_gs_list = DA_gs_TYR_up$ID,
                    gs_DA_Z= DA_gs_bid_TYR[DA_gs_TYR_up$ID,'Z-statistics'],
                    use_gs2gene = use_gs2gene,
                    top_gs_number=20,target_col='RdBu',
                    left_annotation = 'high in mutant',right_annotation = 'high in WT',
                    main= comp_name,Z_sig_thre=1.64,profile_sig_thre = 0,
                    pdf_file=sprintf('%s/NetBID_GSEA_GS_DA_up.pdf',analysis.par$out.dir.PLOT))   

draw.GSEA.NetBID.GS(DE=DA,name_col='ID',profile_col='logFC',profile_trend='pos2neg',
                    sig_gs_list = DA_gs_TYR_H_up$ID,
                    gs_DA_Z= DA_gs_bid_TYR_H[DA_gs_TYR_H_up$ID,'Z-statistics'],
                    use_gs2gene = use_gs2gene,
                    top_gs_number=20,target_col='RdBu',
                    left_annotation = 'high in mutant',right_annotation = 'high in WT',
                    main= comp_name,Z_sig_thre=1.64,profile_sig_thre = 0,
                    pdf_file=sprintf('%s/NetBID_GSEA_GS_DA_H.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID.GS(DE=DA,name_col='ID',profile_col='logFC',profile_trend='pos2neg',
                    sig_gs_list = DA_gs_TYR_kegg_up$ID,
                    gs_DA_Z= DA_gs_bid_TYR_kegg[DA_gs_TYR_kegg_up$ID,'Z-statistics'],
                    use_gs2gene = use_gs2gene,
                    top_gs_number=20,target_col='RdBu',
                    left_annotation = 'high in mutant',right_annotation = 'high in WT',
                    main= comp_name,Z_sig_thre=1.64,profile_sig_thre = 0,
                    pdf_file=sprintf('%s/NetBID_GSEA_GS_DA_kegg.pdf',analysis.par$out.dir.PLOT))





save.image("/Volumes/XuYang/ATRT/analysis/ATRT_combined_Isabel_Xu/combined.RData")









