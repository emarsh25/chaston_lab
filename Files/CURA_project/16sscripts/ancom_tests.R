make_ancom2_plots <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, random.formula=NULL, color = color, breaks = breaks, labs = labs) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)

  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% 
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  
  #  taxonomic_level <- "X.OTU.ID"
    ## Cluster rows by taxonomy hierarchically 
  taxon_vector <- c("kingdom","phylum","class","order","family","genus","species","X.OTU.ID")
  taxon_number <- which(taxon_vector==taxonomic_level)
  taxon_vector2 <- c()
  for (i in 1:taxon_number) {
    taxon_vector2 <- c(taxon_vector2,taxon_vector[i])
  }

  if(taxonomic_level!="X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      tidyr::unite(clustered_taxonomy, (taxon_vector[1:length(taxon_vector2)])) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy))
  } else if (taxonomic_level == "X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      mutate(OTUID2 = X.OTU.ID) %>%
      tidyr::unite(clustered_taxonomy, c(taxon_vector[1:7],"OTUID2")) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      dplyr::select(-X.OTU.ID.1)
        
    } 

  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% 
    dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))

  phyl2 <- otu_table %>% 
    inner_join(otu2, by = "X.OTU.ID") %>%
    group_by(clustered_taxonomy) %>% 
    summarize_at(ref_names,sum, na.rm=T)
  rownames(phyl2) <- as.character(unlist(phyl2$clustered_taxonomy))
  
  phyl3 <- phyl2 %>% dplyr::select(-clustered_taxonomy) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  phyl5[1,]
  map3[1,]
   
  comparison_test = ANCOM(otu_data = phyl5, 
              meta_data = map3, 
              main_var = main.var,  
              zero_cut = .6, 
              p_adjust_method = multcorr, 
              alpha = sig, 
              adj_formula = adj.formula, 
              rand_formula = random.formula)
  
  write.csv(comparison_test, paste("ancomexcel_",main.var,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test) %>% filter(detected_0.9==T) %>% dplyr::select(otu_id)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu_id, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  otu3 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu3$OTU <- as.character(otu3$X.OTU.ID)
  
  if(taxonomic_level!="OTU") {
   otu3[,paste(taxonomic_level)] = gsub(" ","",otu3[,paste(taxonomic_level)])
  }

  otu5 <- otu3 %>% inner_join(otu2, by = "X.OTU.ID") %>% filter(clustered_taxonomy%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)

  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean, fill = detected_taxa[i])) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             scale_fill_manual(values=color[i]) +
             scale_x_discrete(breaks=breaks, labels=labs) +
             theme_cowplot() +
             theme(axis.text=element_text(size=14),
                   axis.text.x = element_text(angle = 45, hjust = 1),
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=20),
                   legend.position = "none") +
             labs(y="Relative Abundance",x=str_to_title(var1),title=tail(strsplit(detected_taxa[i], '__')[[1]], n = 1)))
  }
  
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 100)
  do.call("grid.arrange", c(plist, ncol=length(plist)))
  dev.off()
}

library()


# breaks=c("iPLA2-VIA Control", "iPLA2-VIA[Delta174]", "iPLA2-VIA[Delta192]"), labels=c("iPLA2-VIA (Control)", "iPLA2-VIA[Delta174] (Mutant)", "iPLA2-VIA[Delta192] (Mutant)"



#breaks=c("pink1 RV FM6 (Control)", "pink1 B9 FM6 Mutants"), labels=c("pink1 RV FM6 (Control)", "pink1 B9 FM6 (Mutant)"))

# "park25 (Mutant)", "parkrvA (Control)"