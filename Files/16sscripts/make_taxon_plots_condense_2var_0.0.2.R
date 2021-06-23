# 0.0.2	update Aug 2019 ot make compatible with varying file headers and names in the mapper file and taxonomy folders


make_taxon_plots_condense_2var <- function(file_path, map_path="", mapper_file=mapper_file, var1, var2, taxonomic_level="genus", the_two_var = "twovar", the_id="X.SampleID", rpa_in_chart = 0.05, plot_color = "", plot_order = "", x_val=the_id, read_depth = read_depth, legend_position = legend_position) {

  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t") %>% select_(.dots = list(the_id, var1, var2))

  print(colnames(otu_table))
  ## make melted OTU table
  otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species","OTU")]
  otu2$OTU <- as.character(otu2$OTU)
  melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
  for (i in 2:(dim(otu_table)[2]-7-1)) {
    rm(new_table)
    new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
    melted_table <- rbind(melted_table, new_table)
  }

  ## sanity check
  dim(melted_table)[1]/(59)==dim(new_table)[1]

  ## start merging
  melted_table_no0 <- melted_table %>% filter(Count>-1)
  mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))
  mt0[1,]
  mt0$sample2 <- gsub("\\.","", mt0$Sample)
  mt0$sample2 <- gsub("_","", mt0$sample2)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  mt0[1,]
  mt1 <- mt0 %>% inner_join(map2)
  mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
  mt2$sort_order <- c(1:dim(mt2)[1])
  mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
  mt2$order <- paste(mt2$class, mt2$order,sep="_")
  mt2$family <- paste(mt2$order, mt2$family,sep="_")
  mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
  mt2$species <- paste(mt2$genus, mt2$species,sep="_")
  mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")

  ## make the interactive variable
  mt2$twovar <- paste(mt2[,var1], mt2[,var2],sep="-")
  table(list(mt2$twovar))

  ## find the rare taxa
  rare_taxa <- mt2 %>% group_by_(.dots = taxonomic_level) %>% summarize(phylum.abun=sum(Count)) %>% mutate(rpa = phylum.abun / sum(phylum.abun)) %>% filter(rpa<=rpa_in_chart) %>% select_(.dots = taxonomic_level) %>% unlist

  rt2 <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames("paste('other')", taxonomic_level))

  ## find the abundant taxa
  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) %>% data.frame()

  abundant_taxa <- names(table(list(abun_taxa[,paste(taxonomic_level)])))

  #abundant_taxa <- abundant_taxa[!abundant_taxa %in% c("__")]
  print(abundant_taxa)
  i
  for(i in 1:length(abundant_taxa)) {
    try(space_split_vector <- strsplit(abundant_taxa[i], split = " "),F)
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      try(space_split_vector[[1]] <- head(space_split_vector[[1]], -1),F)#; space_split_vector[[1]]
    }, F)
    try(abundant_taxa[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""),F)
  }
  print(abundant_taxa)

  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames(paste("abundant_taxa"), taxonomic_level)) %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level))

  ## merge the two into a new file
  rta_median_table <- rt2 %>% bind_rows(abun_taxa)	%>% droplevels
  rta_median <- median(table(list(rta_median_table$X.SampleID)))

  rta_time <- data.frame(table(list(rta_median_table$twovar))) %>% mutate(Freq=Freq/rta_median)

  rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=c("twovar" = "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))

  #	rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=structure(names=var1, "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))


  sum(rta$rpa)


  if(length(plot_order) == 1) {
    plot_order <- abundant_taxa
  }

  rta$genus2 <- unlist(rta[,paste0(taxonomic_level)])
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))

  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }
  x_val
  p <- 	 ggplot(rta, aes(x = get(x_val), y = rpa, fill = genus2)) +
    facet_grid(.~get(var1) * get(var2) , drop = TRUE, space = "free", scales = "free") +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(name = "Legend", values=plot_color) + theme_cowplot() + theme(legend.position=legend_position)


  jpeg(h=800, w=1600, paste("taxon_plot_condensed_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  plot(p)
  dev.off()

  return(p)
}

make_taxon_plots_condense_2var_nolegend <- function(file_path, map_path="", var1, var2, mapper_file = mapper_file, taxonomic_level="genus", the_two_var = "twovar", the_id="X.SampleID", rpa_in_chart = 0.05, plot_color = "", plot_order = "", x_val=the_id, read_depth = read_depth, legend_position=legend_position) {

  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))

  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t")# %>% select_(.dots = list(the_id, var1, var2))

  print(colnames(otu_table))
  ## make melted OTU table
  otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species","OTU")]
  otu2$OTU <- as.character(otu2$OTU)
  melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
  for (i in 2:(dim(otu_table)[2]-7-1)) {
    rm(new_table)
    new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
    melted_table <- rbind(melted_table, new_table)
  }

  ## sanity check
  dim(melted_table)[1]/(59)==dim(new_table)[1]

  ## start merging
  melted_table_no0 <- melted_table %>% filter(Count>-1)
  mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))
  mt0[1,]
  mt0$sample2 <- gsub("\\.","", mt0$Sample)
  mt0$sample2 <- gsub("_","", mt0$sample2)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  map2$sample2 <- gsub("\\.","",map2$sample2)

  mt1 <- mt0 %>% inner_join(map2)
  mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
  mt2$sort_order <- c(1:dim(mt2)[1])
  mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
  mt2$order <- paste(mt2$class, mt2$order,sep="_")
  mt2$family <- paste(mt2$order, mt2$family,sep="_")
  mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
  mt2$species <- paste(mt2$genus, mt2$species,sep="_")
  mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")

  ## make the interactive variable
  mt2$twovar <- paste(mt2[,var1], mt2[,var2],sep="-")
  table(list(mt2$twovar))

  ## find the rare taxa
  rare_taxa <- mt2 %>% group_by_(.dots = taxonomic_level) %>% summarize(phylum.abun=sum(Count)) %>% mutate(rpa = phylum.abun / sum(phylum.abun)) %>% filter(rpa<=rpa_in_chart) %>% select_(.dots = taxonomic_level) %>% unlist

  rt2 <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames("paste('other')", taxonomic_level))

  ## find the abundant taxa
  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) %>% data.frame()

  abundant_taxa <- names(table(list(abun_taxa[,paste(taxonomic_level)])))

  #abundant_taxa <- abundant_taxa[!abundant_taxa %in% c("__")]
  print(abundant_taxa)
  for(i in 1:length(abundant_taxa)) {
    try(space_split_vector <- strsplit(abundant_taxa[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      try(space_split_vector[[1]] <- head(space_split_vector[[1]], -1))#; space_split_vector[[1]]
    })
    try(abundant_taxa[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }
  print(abundant_taxa)

  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames(paste("abundant_taxa"), taxonomic_level)) %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level))

  ## merge the two into a new file
  rta_median_table <- rt2 %>% bind_rows(abun_taxa)	%>% droplevels
  rta_median <- median(table(list(rta_median_table$X.SampleID)))

  rta_time <- data.frame(table(list(rta_median_table$twovar))) %>% mutate(Freq=Freq/rta_median)

  rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=c("twovar" = "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))

  #	rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=structure(names=var1, "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))


  sum(rta$rpa)


  if(length(plot_order) == 1) {
    plot_order <- abundant_taxa
  }

  rta$genus2 <- unlist(rta[,paste0(taxonomic_level)])
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))

  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }

  p <- 	 ggplot(rta, aes(x = get(x_val), y = rpa, fill = genus2)) +
    facet_grid(.~get(var1) * get(var2) , drop = TRUE, space = "free", scales = "free") +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=plot_color) + theme_cowplot()+ 	theme(legend.position=legend_position)

  jpeg(h=800, w=1600, paste("taxon_plot_condensed_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  plot(p)
  dev.off()

  return(p)
}


make_taxon_plots_1var <- function(urlpath = urlpath, file_path, map_path="", var1, mapper_file = mapper_file, taxonomic_level="genus", the_two_var = "onevar", the_id="X.SampleID", rpa_in_chart = 0.05, plot_color = "", plot_order = "",x_val=the_id, read_depth = read_depth, legend_position=legend_position) {

  taxon_url = "https://raw.githubusercontent.com/emarsh25/chaston_lab/master/Files/CURA_project/"

  otu_table <- read.table(url(paste(urlpath,'core-metrics-results-',file_path,'/rarefied_table.txt',sep="")), comment.char="", header=T, sep="\t", fill=T, skip=1)%>%left_join(read.csv(url(paste0(taxon_url,'taxonomy',map_path,'/taxonomy_forR.csv'))), by=c("X.OTU.ID"="Feature.ID"))
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t")

  print(colnames(otu_table))
  ## make melted OTU table
  otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species","OTU")]
  otu2$OTU <- as.character(otu2$OTU)
  melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
  for (i in 2:(dim(otu_table)[2]-7-1)) {
    rm(new_table)
    new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
    melted_table <- rbind(melted_table, new_table)
  }

  str(map2)
  ## sanity check
  dim(melted_table)[1]/(59)==dim(new_table)[1]

  melted_table_no0 <- melted_table %>% filter(Count>-1)
  mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))
  mt0[1,]

  mt0

  mt0$sample2 <- gsub("\\.","", mt0$Sample)
  mt0$sample2 <- gsub("_","", mt0$sample2)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  mt0[1,]
  mt1 <- mt0 %>% inner_join(map2)
  mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
  mt2$sort_order <- c(1:dim(mt2)[1])
  mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
  mt2$order <- paste(mt2$class, mt2$order,sep="_")
  mt2$family <- paste(mt2$order, mt2$family,sep="_")
  mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
  mt2$species <- paste(mt2$genus, mt2$species,sep="_")
  mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")

  mt2$onevar <- mt2[,var1]
  table(list(mt2$onevar))

  mt2[1,]
  ## find the rare taxa
  rare_taxa <- mt2 %>% group_by_(.dots = taxonomic_level) %>% summarize(phylum.abun=sum(Count)) %>% mutate(rpa = phylum.abun / sum(phylum.abun)) %>% filter(rpa<=rpa_in_chart) %>% select_(.dots = taxonomic_level) %>% unlist
  rt2 <- mt2 %>% group_by_(.dots=list(var1, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames("paste('other')", taxonomic_level))

  ## find the abundant taxa
  abun_taxa <- mt2 %>% group_by_(.dots=list(var1, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) %>% data.frame()

  abundant_taxa <- names(table(list(abun_taxa[,paste(taxonomic_level)])))

  print(abundant_taxa)
  for(i in 1:length(abundant_taxa)) {
    try(space_split_vector <- strsplit(abundant_taxa[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      space_split_vector[[1]] <- head(space_split_vector[[1]], -1); space_split_vector[[1]]
    })
    try(abundant_taxa[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }
  print(abundant_taxa)

  abun_taxa <- mt2 %>% group_by_(.dots=list(var1, the_two_var, the_id, taxonomic_level)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames(paste("abundant_taxa"), taxonomic_level)) %>% group_by_(.dots=list(var1, the_two_var, the_id, taxonomic_level))

  ## merge the two into a new file
  rta <- rt2 %>% bind_rows(abun_taxa) %>% group_by(onevar) %>% mutate(rpa = phylum.abun / read_depth) %>% filter(rpa>0) %>% droplevels()%>%dplyr::mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))

  rta

  if(length(plot_order) == 1) {
    plot_order <- abundant_taxa
  }

  rta$genus2 <- unlist(rta[,paste0(taxonomic_level)])
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))

  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }

  rta


  p <- 	 ggplot(rta, aes(x = X.SampleID, y = rpa, fill = genus2)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=plot_color) +
    theme_cowplot() +
    theme(legend.position = legend_position,
          axis.text.x = element_blank()) +
    ylab("fractional abundance")


  jpeg(h=800, w=1600, paste("taxon_plot_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  plot(p)
  dev.off()

  return(p)


}

make_taxon_plots_2var <- function(file_path, map_path="", var1, var2, mapper_file = mapper_file, taxonomic_level="genus", the_two_var = "twovar", the_id="X.SampleID", rpa_in_chart = 0.05, plot_color = "", plot_order = "",x_val=the_id, read_depth = read_depth, legend_position=legend_position) {

  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t")# %>% select_(.dots = list(the_id, var1, var2))

  print(colnames(otu_table))
  ## make melted OTU table
  otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species","OTU")]
  otu2$OTU <- as.character(otu2$OTU)
  melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
  for (i in 2:(dim(otu_table)[2]-7-1)) {
    rm(new_table)
    new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
    melted_table <- rbind(melted_table, new_table)
  }

  str(map2)
  ## sanity check
  dim(melted_table)[1]/(59)==dim(new_table)[1]

  ## start merging
  melted_table_no0 <- melted_table %>% filter(Count>-1)
  mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))
  mt0[1,]
  mt0$sample2 <- gsub("\\.","", mt0$Sample)
  mt0$sample2 <- gsub("_","", mt0$sample2)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  mt0[1,]
  mt1 <- mt0 %>% inner_join(map2)
  mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
  mt2$sort_order <- c(1:dim(mt2)[1])
  mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
  mt2$order <- paste(mt2$class, mt2$order,sep="_")
  mt2$family <- paste(mt2$order, mt2$family,sep="_")
  mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
  mt2$species <- paste(mt2$genus, mt2$species,sep="_")
  mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")

  ## make the interactive variable
  mt2$twovar <- paste(mt2[,var1], mt2[,var2],sep="-")
  table(list(mt2$twovar))

  mt2[1,]
  ## find the rare taxa
  rare_taxa <- mt2 %>% group_by_(.dots = taxonomic_level) %>% summarize(phylum.abun=sum(Count)) %>% mutate(rpa = phylum.abun / sum(phylum.abun)) %>% filter(rpa<=rpa_in_chart) %>% select_(.dots = taxonomic_level) %>% unlist
  rt2 <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames("paste('other')", taxonomic_level))

  ## find the abundant taxa
  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) %>% data.frame()

  abundant_taxa <- names(table(list(abun_taxa[,paste(taxonomic_level)])))

  print(abundant_taxa)
  for(i in 1:length(abundant_taxa)) {
    try(space_split_vector <- strsplit(abundant_taxa[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      space_split_vector[[1]] <- head(space_split_vector[[1]], -1); space_split_vector[[1]]
    })
    try(abundant_taxa[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }
  print(abundant_taxa)
  print(abundant_taxa)

  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames(paste("abundant_taxa"), taxonomic_level)) %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level))

  ## merge the two into a new file
  rta <- rt2 %>% bind_rows(abun_taxa) %>% group_by(twovar) %>% mutate(rpa = phylum.abun / read_depth) %>% filter(rpa>0) %>% droplevels()%>%dplyr::mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))

  if(length(plot_order) == 1) {
    plot_order <- abundant_taxa
  }

  rta$genus2 <- unlist(rta[,paste0(taxonomic_level)])
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))

  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }

  p <- 	 ggplot(rta, aes(x = X.SampleID, y = rpa, fill = genus2)) +
    facet_grid(.~get(var1) * get(var2) , drop = TRUE, space = "free", scales = "free") +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=plot_color) +
    theme_cowplot() +
    theme(legend.position = legend_position,
          axis.text.x = element_blank()) +
    ylab("fractional abundance")

  jpeg(h=800, w=1600, paste("taxon_plot_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  plot(p)
  dev.off()

  return(p)
}

make_taxon_plots_condense_2var_period <- function(file_path, map_path="", var1, var2, taxonomic_level="genus", the_two_var = "twovar", the_id="X.SampleID", rpa_in_chart = 0.05, plot_color = "", plot_order = "", x_val=the_id, read_depth = read_depth, legend_position = legend_position) {

  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t")# %>% select_(.dots = list(the_id, var1, var2))

  print(colnames(otu_table))
  ## make melted OTU table
  otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species","OTU")]
  otu2$OTU <- as.character(otu2$OTU)
  melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
  for (i in 2:(dim(otu_table)[2]-7-1)) {
    rm(new_table)
    new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
    melted_table <- rbind(melted_table, new_table)
  }

  ## sanity check
  dim(melted_table)[1]/(59)==dim(new_table)[1]

  ## start merging
  melted_table_no0 <- melted_table %>% filter(Count>-1)
  mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))

  mt0$sample2 <- gsub("\\.","", mt0$Sample)
  mt0$sample2 <- gsub("_","", mt0$sample2)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  map2$sample2 <- gsub("\\.","",map2$sample2)
  map2$sample2 <- paste0("X",map2$sample2)

  mt1 <- mt0 %>% inner_join(map2)
  mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
  mt2$sort_order <- c(1:dim(mt2)[1])
  mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
  mt2$order <- paste(mt2$class, mt2$order,sep="_")
  mt2$family <- paste(mt2$order, mt2$family,sep="_")
  mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
  mt2$species <- paste(mt2$genus, mt2$species,sep="_")
  mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")

  ## make the interactive variable
  mt2$twovar <- paste(mt2[,var1], mt2[,var2],sep="-")
  table(list(mt2$twovar))

  ## find the rare taxa
  rare_taxa <- mt2 %>% group_by_(.dots = taxonomic_level) %>% summarize(phylum.abun=sum(Count)) %>% mutate(rpa = phylum.abun / sum(phylum.abun)) %>% filter(rpa<=rpa_in_chart) %>% select_(.dots = taxonomic_level) %>% unlist

  rt2 <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames("paste('other')", taxonomic_level))

  ## find the abundant taxa
  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) %>% data.frame()

  abundant_taxa <- names(table(list(abun_taxa[,paste(taxonomic_level)])))

  #abundant_taxa <- abundant_taxa[!abundant_taxa %in% c("__")]
  print(abundant_taxa)
  for(i in 1:length(abundant_taxa)) {
    try(space_split_vector <- strsplit(abundant_taxa[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      try(space_split_vector[[1]] <- head(space_split_vector[[1]], -1))#; space_split_vector[[1]]
    })
    try(abundant_taxa[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }
  print(abundant_taxa)

  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames(paste("abundant_taxa"), taxonomic_level)) %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level))

  ## merge the two into a new file
  rta_median_table <- rt2 %>% bind_rows(abun_taxa)	%>% droplevels
  rta_median <- median(table(list(rta_median_table$X.SampleID)))

  rta_time <- data.frame(table(list(rta_median_table$twovar))) %>% mutate(Freq=Freq/rta_median)

  rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=c("twovar" = "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))

  #	rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=structure(names=var1, "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))


  sum(rta$rpa)


  if(length(plot_order) == 1) {
    plot_order <- abundant_taxa
  }

  rta$genus2 <- unlist(rta[,paste0(taxonomic_level)])
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))

  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }

  p <- 	 ggplot(rta, aes(x = get(x_val), y = rpa, fill = genus2)) +
    facet_grid(.~get(var1) * get(var2) , drop = TRUE, space = "free", scales = "free") +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=plot_color) + 	theme(legend.position=legend_position)

  jpeg(h=800, w=1600, paste("taxon_plot_condensed_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  plot(p)
  dev.off()

  return(p)
}


make_taxon_plots_condense_2var_period_name <- function(file_path, var1, var2, taxonomic_level="genus", the_two_var = "twovar", the_id="X.SampleID", rpa_in_chart = 0.05, plot_color = "", plot_order = "", x_val=the_id, read_depth = read_depth) {

  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv('taxonomy/taxonomy_forR.csv'), by=c("X.OTU.ID"="Feature.ID"))
  map2 <- read.table("metadata.tsv",comment.char = "", header=T) %>% select_(.dots = list(the_id, var1, var2))

  print(colnames(otu_table))
  ## make melted OTU table
  otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species","OTU")]
  otu2$OTU <- as.character(otu2$OTU)
  melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
  for (i in 2:(dim(otu_table)[2]-7-1)) {
    rm(new_table)
    new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
    melted_table <- rbind(melted_table, new_table)
  }

  ## sanity check
  dim(melted_table)[1]/(59)==dim(new_table)[1]

  ## start merging
  melted_table_no0 <- melted_table %>% filter(Count>-1)
  mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))
  mt0[1,]
  #mt0$sample2 <- gsub("\\.","", mt0$Sample)
  mt0$sample2 <- gsub("X","", mt0$Sample)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  mt0[1,]
  mt1 <- mt0 %>% inner_join(map2)
  mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
  mt2$sort_order <- c(1:dim(mt2)[1])
  mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
  mt2$order <- paste(mt2$class, mt2$order,sep="_")
  mt2$family <- paste(mt2$order, mt2$family,sep="_")
  mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
  mt2$species <- paste(mt2$genus, mt2$species,sep="_")
  mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")

  ## make the interactive variable
  mt2$twovar <- paste(mt2[,var1], mt2[,var2],sep="-")
  table(list(mt2$twovar))

  ## find the rare taxa
  rare_taxa <- mt2 %>% group_by_(.dots = taxonomic_level) %>% summarize(phylum.abun=sum(Count)) %>% mutate(rpa = phylum.abun / sum(phylum.abun)) %>% filter(rpa<=rpa_in_chart) %>% select_(.dots = taxonomic_level) %>% unlist

  rt2 <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames("paste('other')", taxonomic_level))

  ## find the abundant taxa
  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) %>% data.frame()

  abundant_taxa <- names(table(list(abun_taxa[,paste(taxonomic_level)])))

  #abundant_taxa <- abundant_taxa[!abundant_taxa %in% c("__")]
  print(abundant_taxa)
  for(i in 1:length(abundant_taxa)) {
    try(space_split_vector <- strsplit(abundant_taxa[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      try(space_split_vector[[1]] <- head(space_split_vector[[1]], -1))#; space_split_vector[[1]]
    })
    try(abundant_taxa[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }
  print(abundant_taxa)

  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames(paste("abundant_taxa"), taxonomic_level)) %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level))

  ## merge the two into a new file
  rta_median_table <- rt2 %>% bind_rows(abun_taxa)	%>% droplevels
  rta_median <- median(table(list(rta_median_table$X.SampleID)))

  rta_time <- data.frame(table(list(rta_median_table$twovar))) %>% mutate(Freq=Freq/rta_median)

  rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=c("twovar" = "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))

  #	rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=structure(names=var1, "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))


  sum(rta$rpa)


  if(length(plot_order) == 1) {
    plot_order <- abundant_taxa
  }

  rta$genus2 <- unlist(rta[,paste0(taxonomic_level)])
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))

  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }

  p <- 	 ggplot(rta, aes(x = get(x_val), y = rpa, fill = genus2)) +
    facet_grid(.~get(var1) * get(var2) , drop = TRUE, space = "free", scales = "free") +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(name = "Legend", values=plot_color)

  jpeg(h=800, w=1600, paste("taxon_plot_condensed_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  plot(p)
  dev.off()

  return(p)
}

make_taxon_plots_condense_2var_period_name_nolegend <- function(file_path, var1, var2, taxonomic_level="genus", the_two_var = "twovar", the_id="X.SampleID", rpa_in_chart = 0.05, plot_color = "", plot_order = "", x_val=the_id, read_depth = read_depth) {

  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv('taxonomy/taxonomy_forR.csv'), by=c("X.OTU.ID"="Feature.ID"))
  map2 <- read.table("metadata.tsv",comment.char = "", header=T) %>% select_(.dots = list(the_id, var1, var2))

  print(colnames(otu_table))
  ## make melted OTU table
  otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species","OTU")]
  otu2$OTU <- as.character(otu2$OTU)
  melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
  for (i in 2:(dim(otu_table)[2]-7-1)) {
    rm(new_table)
    new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
    melted_table <- rbind(melted_table, new_table)
  }

  ## sanity check
  dim(melted_table)[1]/(59)==dim(new_table)[1]

  ## start merging
  melted_table_no0 <- melted_table %>% filter(Count>-1)
  mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))
  mt0[1,]
  #mt0$sample2 <- gsub("\\.","", mt0$Sample)
  mt0$sample2 <- gsub("X","", mt0$Sample)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  mt0[1,]
  mt1 <- mt0 %>% inner_join(map2)
  mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
  mt2$sort_order <- c(1:dim(mt2)[1])
  mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
  mt2$order <- paste(mt2$class, mt2$order,sep="_")
  mt2$family <- paste(mt2$order, mt2$family,sep="_")
  mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
  mt2$species <- paste(mt2$genus, mt2$species,sep="_")
  mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")

  ## make the interactive variable
  mt2$twovar <- paste(mt2[,var1], mt2[,var2],sep="-")
  table(list(mt2$twovar))

  ## find the rare taxa
  rare_taxa <- mt2 %>% group_by_(.dots = taxonomic_level) %>% summarize(phylum.abun=sum(Count)) %>% mutate(rpa = phylum.abun / sum(phylum.abun)) %>% filter(rpa<=rpa_in_chart) %>% select_(.dots = taxonomic_level) %>% unlist

  rt2 <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames("paste('other')", taxonomic_level))

  ## find the abundant taxa
  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) %>% data.frame()

  abundant_taxa <- names(table(list(abun_taxa[,paste(taxonomic_level)])))

  #abundant_taxa <- abundant_taxa[!abundant_taxa %in% c("__")]
  print(abundant_taxa)
  for(i in 1:length(abundant_taxa)) {
    try(space_split_vector <- strsplit(abundant_taxa[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      try(space_split_vector[[1]] <- head(space_split_vector[[1]], -1))#; space_split_vector[[1]]
    })
    try(abundant_taxa[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }
  print(abundant_taxa)

  abun_taxa <- mt2 %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level)) %>% filter(get(taxonomic_level)%in%rare_taxa==F) %>% summarise(phylum.abun = sum(Count)) %>% mutate_(.dots = setNames(paste("abundant_taxa"), taxonomic_level)) %>% group_by_(.dots=list(var1,var2, the_two_var, the_id, taxonomic_level))

  ## merge the two into a new file
  rta_median_table <- rt2 %>% bind_rows(abun_taxa)	%>% droplevels
  rta_median <- median(table(list(rta_median_table$X.SampleID)))

  rta_time <- data.frame(table(list(rta_median_table$twovar))) %>% mutate(Freq=Freq/rta_median)

  rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=c("twovar" = "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))

  #	rta <- rt2 %>% bind_rows(abun_taxa) %>% inner_join(rta_time, by=structure(names=var1, "Var1")) %>% group_by(twovar) %>% mutate(rpa = phylum.abun/(Freq* read_depth)) %>% filter(rpa>0) %>% droplevels()%>%mutate_(.dots=setNames("factor(get(taxonomic_level))", taxonomic_level)) #%>% mutate_(.dots=setNames(factor(rta[,paste(taxonomic_level)], levels=c("other",plot_order)), taxonomic_level))


  sum(rta$rpa)


  if(length(plot_order) == 1) {
    plot_order <- abundant_taxa
  }

  rta$genus2 <- unlist(rta[,paste0(taxonomic_level)])
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))

  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }

  rta$sex

  p <- 	 ggplot(rta, aes(x = get(x_val), y = rpa, fill = genus2)) +
    facet_grid(.~get(var1) * get(var2) , drop = TRUE, space = "free", scales = "free") +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(name = "Legend", values=plot_color) + 	theme(legend.position="none")


  jpeg(h=800, w=1600, paste("taxon_plot_condensed_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  plot(p)
  dev.off()

  return(p)
}
