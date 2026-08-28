
####################################
## Utility functions for scRNA-seq
## data aanlysis
####################################

get_ncells <- function(list_of_sobjs){
  ncells_list <- list()
  for (i in 1:length(list_of_sobjs)){
    sobj <- list_of_sobjs[[i]]
    ncells <- data.frame(sample = sobj@project.name, ncells = ncol(sobj@assays$RNA@layers$counts))
    ncells_list[[i]] <- ncells
  }
  ncells <- dplyr::bind_rows(ncells_list)  
  return(ncells)
}

###################################
###################################
###################################

filter_sobjs <- function(list_of_sobjs,
                         read_count,
                         feature_count,
                         feature_count_max,
                         mt_percent,
                         scrubdir){
  res_list <- list()
  for (i in 1:length(list_of_sobjs)){
    s_obj <- list_of_sobjs[[i]]
    # filter
    filtered_sobj <- subset(s_obj, subset = nCount_RNA >= read_count & nFeature_RNA >=feature_count & nFeature_RNA < feature_count_max & percent.mt < mt_percent)
    
    # scrublet
    id <- filtered_sobj@meta.data$orig.ident[1]
    f = paste0(scrubdir, id, ".tsv.gz")
    
    scrubbed <- read.csv(f, sep="\t", header=TRUE)
    to_remove <- scrubbed[scrubbed$scrub_predicted_doublets == "True", ]$barcode
    
    flog.info(paste0("removing ", length(to_remove), " based on scrublet"))
    
    filtered_sobj <- subset(filtered_sobj, cells = setdiff(rownames(filtered_sobj@meta.data), to_remove))
    
    res_list[[i]] <- filtered_sobj  
  }
  return(res_list)
}

###################################
###################################
###################################

plot_props <- function(sobj, metadata, annotation_colors){
  
  annos <- sobj@meta.data
  annos$annotation <- Idents(sobj)
  
  annos <- annos %>% 
    group_by(orig.ident, annotation) %>%
    summarize(n = n())
  annos$disease_type <- metadata[annos$orig.ident,]$disease_type
  
  ggplot(annos, aes(x=orig.ident, y=n, fill=annotation)) +
    geom_col(position = "fill") +
    theme_classic() +
    scale_fill_manual(values = annotation_colors) +
    facet_wrap(~disease_type, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

###################################
###################################
###################################
feature_plot <- function(obj, feature = "Alpk1", high="red3"){
  
  embeddings <- as.data.frame(obj@reductions$pca@cell.embeddings)[,c(1,2)]
  expression <- unlist(obj@assays$RNA$data[feature,])
  embeddings$expression <- expression
  colnames(embeddings) <- c("PC_1", "PC_2", "expression")
  embeddings_sub <- embeddings[embeddings$expression > 0,]
  
  ggplot(embeddings, aes(x=PC_1, y=PC_2, color=expression)) +
    geom_point(size = 0.5) +
    theme_classic() +
    scale_color_gradient(low = "snow2", high = high) +
    geom_point(size = 0.5, data = embeddings_sub, aes(x=PC_1, y=PC_2, color=expression), inherit.aes=FALSE)
}

###################################
###################################
###################################

pathway_fold_enrichment <- function(res, nsig, pathways, universe){
  
  # iterate over pathways and estimate fold changes
  for (i in 1:nrow(res)){
    pathway <- res[i, "pathway"]
    pathway_set <- unlist(pathways[[pathway]])
    background_overlap <- length(intersect(universe, pathway_set))
    overlap <- res[i, "overlap"]
    fe <- (overlap/nsig)/(background_overlap/length(universe))
    res[i, "foldEnrichment"] <- fe
  }
  return(res)
}

###################################
###################################
###################################

volcano_plot <- function(x, use="pvalue", value=0.05, labels = NULL){
  
  if (is.null(labels)){
    x$label <- ifelse(x[,use] < value &!(is.na(x[,use])), rownames(x), NA)
  }else{
    x$label <- ifelse(rownames(x) %in% labels, rownames(x), NA)
  }
  
  x$color <- ifelse(x[,use] < value &!(is.na(x[,use])), "Sig", NA)
  ggplot(x, aes(x=log2FoldChange, y = -log10(pvalue), color=color)) +
    geom_point(aes(size=log2(baseMean)), alpha=0.5) +
    geom_vline(xintercept=c(-1,1), linetype = "dashed") +
    ggrepel::geom_text_repel(label = x$label, max.overlaps=20) +
    theme_classic() +
    scale_color_manual(values = c("blue4")) 
}


###################################
###################################
###################################

get_ppos <- function(obj, gene, type = "orig.ident"){
  
  counts <- as.data.frame(GetAssayData(obj))
  ppos <- data.frame(individual = obj@meta.data[,type],
                     counts = unlist(counts[gene,]))
  ppos <- ppos %>%
    group_by(individual) %>%
    summarise(npos = sum(counts > 0), n = n()) %>%
    mutate(ppos = (npos/n)*100)
  return(ppos)
  
}


