###immune_targets_GS analysis####
load("/Volumes/XuYang/ATRT/analysis/ATRT_combined_Isabel_Xu/geneset_lists.RData")
use_gs2gene <- merge_gs(all_gs2gene=curated,use_gs=c("kelli", "immunotx", "genesig"))
use_gs2gene2 <- merge_gs(all_gs2gene=all_gs2gene,use_gs='H')
use_gs2gene3<-c(use_gs2gene,use_gs2gene2)
ac_gs <- cal.Activity.GS(use_gs2gene = use_gs2gene3,cal_mat = ac_mat_gene)

#MYC
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup_new`=='MYC')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup_new`!='MYC')] # Combine other groups as the Control group
DA_gs_bid_MYC <- getDE.BID.2G(eset=generate.eset(ac_gs),G1=G1,G0=G0,G1_name='MYC',G0_name='others')
sig_gs <- draw.volcanoPlot(dat= DA_gs_bid_MYC,label_col='ID',logFC_col='logFC',
                           Pv_col='P.Value',logFC_thre=0.0,Pv_thre=1e-3,
                           main='Volcano Plot for gene sets',show_label=TRUE,label_type = 'distribute',label_cex = 0.5,
                           pdf_file=sprintf('%s/volcano_GS_DA_MYC.pdf',analysis.par$out.dir.PLOT))

phe_info_2 <- phe_info[order(phe_info$Subgroup_new),]
draw.heatmap(mat=ac_gs[sig_gs$ID,phe_info_2$Identifier],
             pdf_file=sprintf('%s/heatmap_kelli.pdf',analysis.par$out.dir.PLOT),
             scale='row',
             phenotype_info=phe_info_2,use_phe='Subgroup_new',
             cluster_columns = F, cluster_rows = T)

####immune_target heatmap####
imm_gene<-list()
imm_gene<-c(curated$kelli,curated$immunotx)
imm_gene<-as.matrix(imm_gene)
for(i in 1:length(imm_gene)){
  if (i==1){
  imm_genes<-imm_gene[[i]]
  }
  else{
    imm_genes<-c(imm_genes,imm_gene[[i]]) 
  }
}

imm_genes<-unique(imm_genes)
#for heatmap, first order the samples followed by subgroup_new
phe_info2<-phe_info[order(phe_info$Subgroup_new),]
ac_mat2<-ac_mat[,rownames(phe_info2)]
#find the drug target gene in ms_tab and unique the name
driver_probename<-ms_tab[which(ms_tab$geneSymbol %in% imm_genes),] #446
driver_probename<-driver_probename[order(driver_probename$Size_Xu,decreasing = T),] #446
driver_probename<-driver_probename[!duplicated(driver_probename$geneSymbol),] #181
driver<-driver_probename$originalID_label
draw.heatmap(mat=ac_mat2,use_genes=driver,use_gene_label=driver_probename[driver,'gene_label'],
             use_samples=colnames(ac_mat2),
             phenotype_info=phe_info2,use_phe=c("Subgroup_new","Metastatic.Status_new"),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=F,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_immunothreapy.pdf',analysis.par$out.dir.PLOT))

#shh immuno drug target  #9
#DA  9
immuno_SHH<-ms_tab_SHH[which(ms_tab_SHH$originalID_label %in% driver),]  
immuno_SHH$immune_target<-"SHH"
#tyr immuno drug target  #19
#DA 19
immuno_TYR<-ms_tab_TYR[which(ms_tab_TYR$originalID_label %in% driver),] 
immuno_TYR$immune_target<-"TYR"
#MYC immuno drug target   #77
#DA  77
immuno_MYC<-ms_tab_MYC[which(ms_tab_MYC$originalID_label %in% driver),] 
immuno_MYC$immune_target<-"MYC"

immune_tot<-rbind(immuno_MYC,immuno_TYR)
immune_tot<-rbind(immune_tot,immuno_SHH)

out2excel(immune_tot,"/Volumes/XuYang/ATRT/analysis/ATRT_combined_Isabel_Xu/ATRT_combined_2020-09-22/DATA/immunetarget.xlsx")





