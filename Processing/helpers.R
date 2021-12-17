check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg, update = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

needed_packages <- c("Seurat","devtools","Matrix","scales","ensembldb","biomaRt","gghighlight","patchwork","ggplot2","DropletUtils","ggrepel","ggforce","cluster","ggplotify","tidyverse","SingleCellExperiment","scater","purrr","stringr","scran", "cowplot","RColorBrewer","BiocSingular","irlba","wesanderson","clustree","ggpubr","scDblFinder","writexl","zeallot","ComplexHeatmap")

check.packages(needed_packages)

options(width = 900)
options(future.globals.maxSize= (10050*1024^2))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
date <- as.character(Sys.Date())

kneeplot <- function(bc_rank,upper,lower){
  
  hl_ls <- upper
  ll_ls <- lower
  if(!exists("ll_ls")){ll_ls=800}
  if(!exists("hl_ls")){hl_ls=100000}
  print(ll_ls)
  print(hl_ls)
  
  bc_rank_knee_plot <- qplot(bc_rank$rank, bc_rank$total, geom = "point") +
    geom_hline(yintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
    geom_hline(yintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
    geom_hline(yintercept = c(ll_ls,hl_ls), color = "red", linetype = 2) +
    annotate("text", x = 1000, y = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection, ll_ls,hl_ls),
             label = c(paste0("Knee: ",metadata(bc_rank)$knee), paste0("Inflection: ",metadata(bc_rank)$inflection),"Lower Limit", "Upper Limit"), color = c("blue", "green","red","red")) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank", y = "Total UMI counts") +
    theme_bw()
  bc_rank_knee_plot
  
}


Import_STAR_Velocyto <- function(folder_path){

  print("Loading RNA matrix...")
  RNA_mat <- readMM(paste0(folder_path,"Gene/raw/matrix.mtx"))
  barcodes <- read.table(file = paste0(folder_path,"Gene/raw/barcodes.tsv"), sep = '\t', header = FALSE)[,1]
  features <- read.table(file = paste0(folder_path,"Gene/raw/features.tsv"), sep = '\t', header = FALSE)[,2]
  colnames(RNA_mat) <- barcodes
  rownames(RNA_mat) <- features

  print("Done. Loading Velocyto matrix...")
  t <- readr::read_delim(paste0(folder_path,"Velocyto/raw/matrix.mtx"), delim = ' ', skip = 3, col_names = FALSE)
  barcodes <- read.table(file = paste0(folder_path,"Velocyto/raw/barcodes.tsv"), sep = '\t', header = FALSE)[,1]
  features <- read.table(file = paste0(folder_path,"Velocyto/raw/features.tsv"), sep = '\t', header = FALSE)[,2]
  dim(t)

  spliced_mat <- Matrix::sparseMatrix(
    i = t$X1,
    j = t$X2,
    x = t$X3,
    dims=c(length(features),length(barcodes))
  )


  unspliced_mat <- Matrix::sparseMatrix(
    i = t$X1,
    j = t$X2,
    x = t$X4,
    dims=c(length(features),length(barcodes))
  )

  colnames(spliced_mat) <- barcodes
  rownames(spliced_mat) <- features

  colnames(unspliced_mat) <- barcodes
  rownames(unspliced_mat) <- features


  mats <- list(RNA_mat,spliced_mat,unspliced_mat)
  names(mats) <- c("RNA","spliced","unspliced")
  print("Done.")
  return(mats)
}

Cells_per_Cluster <- function(seurat=NA,sample_name=opt$sample,sample.ident="SAMPLE",cluster.ident="seurat_clusters",level.order=c("control_1","control_2","cKO_1","cKO_2")){

  if(class(seurat@meta.data[[cluster.ident]]) != "factor"){
    seurat@meta.data[[cluster.ident]] <- as.factor(seurat@meta.data[[cluster.ident]])
  }

  table_samples_by_clusters <- as.data.frame(seurat@meta.data) %>%
    dplyr::rename("sample_name":= sample.ident,"clusters_used" :=cluster.ident) %>%
    dplyr::group_by(sample_name, clusters_used) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    spread(clusters_used, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c('sample_name', 'total_cell_count', everything())) %>%
    arrange(factor(sample_name, levels = levels(seurat@meta.data[[sample.ident]])))

  table_samples_by_clusters_percent <- as.data.frame(seurat@meta.data)%>%
    dplyr::rename("sample_name":= sample.ident,"clusters_used" :=cluster.ident)%>%
    dplyr::group_by(sample_name, clusters_used) %>%
    dplyr::summarize(count = dplyr::n())%>%
    spread(clusters_used, count, fill = 0)%>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    pivot_longer(-c(sample_name,total_cell_count),names_to="clusters") %>%
    dplyr::mutate(percent=(value/total_cell_count)*100) %>%
    dplyr::select(-value) %>% pivot_wider(names_from=clusters,values_from=percent)


  table_clusters_by_samples <- seurat@meta.data %>%
    dplyr::rename("sample_name":= sample.ident,"clusters_used" :=cluster.ident) %>%
    group_by(clusters_used, sample_name) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    spread(sample_name, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c('clusters_used', 'total_cell_count', everything())) %>%
    arrange(factor(clusters_used, levels = levels(seurat@meta.data[[cluster.ident]])))

  table_clusters_by_samples_percent <- seurat@meta.data %>%
    dplyr::rename("sample_name":= sample.ident,"clusters_used" :=cluster.ident) %>%
    group_by(clusters_used, sample_name) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    spread(sample_name, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    pivot_longer(-c(clusters_used,total_cell_count),names_to="clusters") %>%
    dplyr::mutate(percent=(value/total_cell_count)*100) %>%
    dplyr::select(-value) %>% pivot_wider(names_from=clusters,values_from=percent)


  p1 <- table_samples_by_clusters %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'sample_name') %>%
    mutate(sample_name = factor(sample_name, levels = c(level.order) ) ) %>%
    ggplot(aes(sample_name, value, fill = variable)) +
    geom_bar(position = 'stack', stat = 'identity') +
    scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
    scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
    theme_bw() +
    theme(
      legend.position = 'left',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

  p1_percent <- table_samples_by_clusters %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'sample_name') %>%
    mutate(sample_name = factor(sample_name, levels = c(level.order) ) ) %>%
    ggplot(aes(sample_name, value, fill = variable)) +
    geom_bar(position = 'fill', stat = 'identity') +
    scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
    scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
    theme_bw() +
    theme(
      legend.position = 'left',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

  p2 <- table_clusters_by_samples %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'clusters_used') %>%
    mutate(clusters_used = factor(clusters_used, levels = levels(seurat@meta.data[[cluster.ident]])),variable = factor(variable, levels = c(level.order) ) ) %>%
    ggplot(aes(clusters_used, value, fill = variable)) +
    geom_bar(position = 'stack', stat = 'identity') +
    scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
    scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
    theme_bw() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

  p2_percent <- table_clusters_by_samples %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'clusters_used') %>%
    mutate(clusters_used = factor(clusters_used, levels = levels(seurat@meta.data[[cluster.ident]])),variable = factor(variable, levels = c(level.order) ) ) %>%
    ggplot(aes(clusters_used, value, fill = variable)) +
    geom_bar(position = 'fill', stat = 'identity') +
    scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
    scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
    theme_bw() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )



  absolute_plot <- p1*p2+plot_layout(ncol = 2, widths = c(
    seurat@meta.data[[sample.ident]] %>% unique() %>% length(),
    seurat@meta.data[[cluster.ident]] %>% unique() %>% length()
  ))

  percent_plot <- p1_percent*p2_percent+plot_layout(ncol = 2, widths = c(
    seurat@meta.data[[sample.ident]] %>% unique() %>% length(),
    seurat@meta.data[[cluster.ident]] %>% unique() %>% length()
  ))

  return(list(absolute_plot,percent_plot))
}



`%notin%` <- Negate(`%in%`)


Enrichr_Plot <- function(markers.list,enrichr.database,top.n=50,export=FALSE,export.path,group.name,width=10,height=20){

  enrichr_results <- enrichr(markers.list,enrichr.database)

  top_20 <- enrichr_results[[1]] %>% dplyr::top_n(n=top.n,wt= desc(Adjusted.P.value)   ) %>% dplyr::arrange(Adjusted.P.value)

  top_20$Term <- factor(top_20$Term,levels = rev( as.vector(top_20$Term)) )

  plot <- ggplot(top_20,aes(x=-log10(Adjusted.P.value),y=Term,fill=Adjusted.P.value)) + geom_col() + geom_vline(xintercept=-log10(0.01),color="red",linetype="dotted") + theme_cowplot()

  if(export == TRUE){

    ggsave(paste0(export.path,"Enrichr_BarPlot_Gene_Sets_",enrichr.database,"_",group.name,'.pdf'),plot, width=width, height=height, limitsize=F)
    writexl::write_xlsx( enrichr_results[[1]], path=paste0(export.path,"Enrichr_Gene_Sets_",enrichr.database,"_",group.name,'.xlsx'),col_names = TRUE )
  }

  return(plot)
}



plot_pct_genes <- function(mat, top_n = 20) {
  pct_tx <- rowSums(mat)
  gs <- rownames(mat)[order(-pct_tx)]
  df <- as.data.frame(t(mat[gs[1:20],]))
  df <- df %>%
    mutate_all(function(x) x/colSums(mat)) %>%
    pivot_longer(everything(), names_to = "gene")
  
  df %>%
    mutate(gene = fct_reorder(gene, value, .fun = median)) %>%
    ggplot(aes(gene, value)) +
    geom_boxplot() +
    labs(x = "", y = "Proportion of total counts") +
    coord_flip()
}


UMAP_Plot <- function(seurat,pt.size=0.1,group.by="final_cell_type",reduction="umap",split.by=NA,facet_labels=NA,border_col=NA,colors=cell_type_main_cols,legend=TRUE,tick_labels=TRUE,n.row=1,n.col=NA,x.axis_label=TRUE,y.axis_label=TRUE,x.axis_line=TRUE,y.axis_line=TRUE){
  
  umap_data <- as.data.frame(seurat@reductions[[reduction]]@cell.embeddings) %>% rownames_to_column(var="cell_name")
  colnames(umap_data) <- c("cell_name","x","y")
  meta_data <- as.data.frame(seurat@meta.data) %>% rownames_to_column(var="cell_name")
  
  if(!is.na(split.by) ){
    
    meta_data <- meta_data[,colnames(meta_data) == group.by | colnames(meta_data) == split.by | colnames(meta_data) == "cell_name" ]
    
    umap_data <- umap_data %>% left_join(meta_data,by="cell_name")
    
    umap_data[[split.by]] <- factor(umap_data[[split.by]], levels=names(facet_labels) )
    
  }else{
    
    meta_data <- meta_data[,colnames(meta_data) == group.by | colnames(meta_data) == "cell_name" ]
    
    umap_data <- umap_data %>% left_join(meta_data,by="cell_name")
    
  }
  
  
  
  umap_plot <- ggplot(umap_data,aes(x=x,y=y,fill= !!sym(group.by),color= !!sym(group.by) ))
  
  if( !is.na(border_col) ){
    umap_plot <- umap_plot + geom_point(size=pt.size,pch=21,color=border_col) + scale_fill_manual("Cell Type", values = colors)
  }else{
    umap_plot <- umap_plot + geom_point(size=pt.size) + scale_color_manual("Cell Type", values = colors)
  }
  
  umap_plot <- umap_plot  + theme_cowplot() +
    labs(x="UMAP 1",y="UMAP 2") +
    theme(axis.title.x = element_text(hjust=0.2,vjust=1),axis.title.y = element_text(hjust=0.2,vjust=1)) + guides(fill="legend")
  
  if(legend == FALSE){umap_plot <- umap_plot + theme(legend.position = "none")}
  
  if(tick_labels == FALSE){umap_plot <- umap_plot + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())}
  
  if(!is.na(split.by) ){umap_plot <- umap_plot +
    facet_wrap(c(split.by),labeller=as_labeller(facet_labels),nrow=n.row,ncol=n.col )+
    theme(strip.background = element_blank())
  
  }
  
  if(x.axis_label == FALSE){umap_plot <- umap_plot + theme(axis.title.x = element_text(colour = "white"))}
  if(y.axis_label == FALSE){umap_plot <- umap_plot + theme(axis.title.y = element_text(colour = "white"))}
  if(x.axis_line == FALSE){umap_plot <- umap_plot + theme(axis.line.x=element_blank())}
  if(y.axis_line == FALSE){umap_plot <- umap_plot + theme(axis.line.y=element_blank())}
  
  
  umap_plot
}

Feature_Plot <- function(seurat,pt.size=0.1,scale.label="Expression",features=c("Acan"),split.by=NA,legend=TRUE,tick_labels=TRUE,use_assay="RNA",use_slot="data",n.row=2,n.col=3,color_lower="grey85",color_upper="deepskyblue3",max.cutoff="q90",min.cutoff=NA,reduction="umap"){
  
  feature_data <- as.data.frame(seurat@reductions[[reduction]]@cell.embeddings) %>% rownames_to_column(var="cell_name")
  colnames(feature_data) <- c("cell_name","x","y")
  meta_data <- as.data.frame(seurat@meta.data) %>% rownames_to_column(var="cell_name")
  read_data <- as.data.frame( t( GetAssayData(seurat,assay=use_assay,slot=use_slot) )  ) %>% rownames_to_column(var="cell_name")
  read_data <- read_data[,c("cell_name",features)] %>% pivot_longer(cols=-cell_name,names_to="Feature",values_to="Expression")
  
  
  feature_data_final <- feature_data %>% left_join(meta_data,by="cell_name") %>% full_join(read_data,by="cell_name")
  feature_data_final[["Feature"]] <- factor(feature_data_final[["Feature"]],levels = features)
  
  
  if(!is.na(max.cutoff) & startsWith( as.character(max.cutoff) ,"q")){
    quantile <- as.numeric( paste0(".",str_extract(max.cutoff,"\\d{2}") ) )
    limit_max <- as.numeric( quantile(feature_data_final$Expression, c(quantile) ) )
    
  }else{
    limit_max <- ifelse( !is.na(max.cutoff),max.cutoff,max(feature_data_final$Expression) )
  }
  
  if(!is.na(min.cutoff) & startsWith( as.character(min.cutoff),"q")){
    quantile <- as.numeric( paste0(".",str_extract(min.cutoff,"\\d{2}") ) )
    limit_min <- as.numeric( quantile(feature_data_final$Expression, c(quantile) ) )
    
  }else{
    limit_min <- ifelse( !is.na(min.cutoff),min.cutoff,min(feature_data_final$Expression) )
  }
  
  
  
  feature_plot <- ggplot(feature_data_final,aes(x=x,y=y,color= Expression,group=split.by )) +
    geom_point(size=pt.size)  + theme_cowplot()
  
  if(any(feature_data_final$Expression <0)){
    feature_plot <- feature_plot +
      scale_colour_gradient2(low=color_lower,high=color_upper,mid="white",midpoint=0,limits=c(limit_min,limit_max),oob=squish)
  } else{
    feature_plot <- feature_plot +
      scale_colour_gradient(low=color_lower,high=color_upper,limits=c(limit_min,limit_max),oob=squish)
  }
  
  feature_plot <- feature_plot +
    labs(x="UMAP 1",y="UMAP 2",color=scale.label) +
    theme(axis.title.x = element_text(hjust=0.2,vjust=1),axis.title.y = element_text(hjust=0.2,vjust=1))
  
  if(legend == FALSE){feature_plot <- feature_plot + theme(legend.position = "none")}
  
  if(tick_labels == FALSE){feature_plot <- feature_plot + theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())}
  
  if(length(features) > 1 & is.na(split.by)){feature_plot <- feature_plot +
    facet_wrap(~Feature,nrow=n.row,ncol=n.col )+
    theme(strip.background = element_blank())
  
  }
  
  
  feature_plot
}
