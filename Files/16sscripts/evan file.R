

### taxon plots funciton - read it in
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
  
  rta <- rta %>%droplevels()
  p <- ggplot(rta, aes(x = X.SampleID, y = rpa, fill = genus2)) + 
    facet_grid(.~get(var1) * get(var2) , drop = TRUE, space = "free", scales = "free") + 
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=plot_color) +
    theme_cowplot() + 
    theme(legend.position = legend_position,
          axis.text.x = element_blank()) +
    ylab("fractional abundance")
  
  
  # ggplot(rta, aes(x = X.SampleID, y = rpa, fill = genus2)) + 
  #   facet_grid(.~get(var1) * get(var2) , drop = TRUE, space = "free", scales = "free") + 
  #   geom_bar(stat = "identity", width = 1) +
  #   scale_fill_manual(values=plot_color) +
  #   theme_cowplot() + 
  #   theme(legend.position = legend_position,
  #         axis.text.x = element_blank()) +
  #   ylab("fractional abundance") + 
  #   facet_wrap(.~genotype_sex , strip.position = "bottom", scales = "free_x", nrow = 1, as.table = F, drop = T) 
  #   
  # 
  jpeg(h=800, w=1600, paste("taxon_plot_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  plot(p)
  dev.off()
  
  return(p)
}

### pcoa with no cirlces - read it in
pcoa_2var <- function(folder_name, mapper_file = mapper_file, var1, var2, source_var1, var1_colors, var2_shape, title_name_add, legend_position = "bottom") {
  
  wfpc <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt',sep=""), sep="\t", fill=T, skip = 9,blank.lines.skip = T, header=F),-2)
  pc2 <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt', sep = ""), sep="\t", fill=T, skip = 3, header=F, colClasses = "character"),2)
  pc_values <- as.numeric(as.character(pc2[2,]))*100
  wfpc <- wfpc[,1:4] 
  fut2 <- read.table(paste('core-metrics-results-',folder_name,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t", stringsAsFactors = T) %>% dplyr::select(X) %>% mutate(X=as.character(X)) %>% left_join(read.table(paste0(mapper_file),comment.char = "", header=T, sep="\t"), by=c("X"="X.SampleID"))
  
  uwmpc_all <- wfpc %>% inner_join(fut2, by=c("V1"="X")) #%>% mutate_(inoculated_with=as.factor(inoculated_with))
  
  
  uwmpc_all[,paste(var1)] <- as.factor(uwmpc_all[,paste(var1)])
  uwmpc_all[,paste(var2)] <- as.factor(uwmpc_all[,paste(var2)])
  
  uwmpc_all[1,]
  var1_table <- table(list(uwmpc_all[,var1]))
  var2_table <- table(list(uwmpc_all[,var2]))
  
  ggplot(uwmpc_all, aes(x = V2, y = V3, colour=get(var1), shape=get(var2), fill=get(var1))) + 
    geom_point(alpha = 1, size=2.5) + #, shape=plot_shapes, size=3, col=plot_colors) + 
    scale_color_manual(name=var1, labels=names(var1_table)[var1_table!=0], values=var1_colors) + 
    scale_fill_manual(name=var1, labels=names(var1_table)[var1_table!=0], values=var1_colors) + 
    scale_shape_manual(name=var2,labels=names(var2_table)[var2_table!=0], values=var2_shape) + 
    #		scale_fill_manual(name="Legend", values=c("red","black","white")) +
    theme(axis.text=element_text(size=12), 
          panel.background = element_blank(), 
          axis.line = element_line(), 
          axis.ticks=element_line(), 
          axis.title=element_text(size=14),	
          #           title=element_text(size=16),
          legend.position=legend_position,
          plot.title = element_text(size=16, hjust=0))+
    labs(x = paste("PCo1 ( ",round(as.numeric(pc_values[1]),1),"% )",sep=""), 
         y = paste("PCo2 ( ",round(as.numeric(pc_values[2]),1),"% )",sep=""), 
         title = title_name_add) + 
    stat_ellipse(aes(x=V2, y=V3,group= get(var1)),            level = .9, show.legend = F, type = "t", geom = "polygon", alpha = 0, inherit.aes=T)
}

### pcoa with circles - read it in
pcoa_2var_circle_2var <- function(folder_name, mapper_file, var1, var2, source_var1, var1_colors, var2_shape, title_name_add, legend_position = "bottom", var2_line = var2_line) {
  
  wfpc <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt',sep=""), sep="\t", fill=T, skip = 9,blank.lines.skip = T, header=F),-2)
  pc2 <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt', sep = ""), sep="\t", fill=T, skip = 3, header=F, colClasses = "character"),2)
  pc_values <- as.numeric(as.character(pc2[2,]))*100
  wfpc <- wfpc[,1:4] 
  fut2 <- read.table(paste('core-metrics-results-',folder_name,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t", stringsAsFactors = T) %>% dplyr::select(X) %>% mutate(X=as.character(X)) %>% left_join(read.table(mapper_file,comment.char = "", header=T, sep = "\t"), by=c("X"="X.SampleID"))
  
  uwmpc_all <- wfpc %>% inner_join(fut2, by=c("V1"="X")) %>% mutate(two_var=as.factor(paste0(get(var1),"_",get(var2))))
  
  uwmpc_all[,paste(var1)] <- as.factor(uwmpc_all[,paste(var1)])
  uwmpc_all[,paste(var2)] <- as.factor(uwmpc_all[,paste(var2)])
  
  uwmpc_all[1,]
  var1_table <- table(list(uwmpc_all[,var1]))
  var2_table <- table(list(uwmpc_all[,var2]))
  twovar_table <- table(list(uwmpc_all$two_var))
  # 
  # 	ggplot(uwmpc_all, aes(x = V2, y = V3, colour=get(var1), shape=get(var2), fill=get(var1)))+
  # 		geom_point(alpha = 1, size=2.5) + #, shape=plot_shapes, size=3, col=plot_colors) + 
  # 		scale_color_manual(name=var1, labels=names(twovar_table)[twovar_table!=0], values=var1_colors) + 
  # 		scale_fill_manual(name=var1, labels=names(twovar_table)[twovar_table!=0], values=var1_colors) + 
  # 		scale_shape_manual(name=var2,labels=names(var2_table)[var2_table!=0], values=var2_shape) + 
  # #		scale_linetype_manual(name=var2, labels=names(var2_table)[var2_table!=0], values=c(1,2)) +
  # 		#		scale_fill_manual(name="Legend", values=c("red","black","white")) +
  # 		theme(axis.text=element_text(size=12), 
  # 					panel.background = element_blank(), 
  # 					axis.line = element_line(), 
  # 					axis.ticks=element_line(), 
  # 					axis.title=element_text(size=14),	
  # 					#           title=element_text(size=16),
  # 					legend.position="right",
  # 					plot.title = element_text(size=16, hjust=0))+
  # 		labs(x = paste("PCo1 ( ",round(as.numeric(pc_values[1]),1),"% )",sep=""), 
  # 				 y = paste("PCo2 ( ",round(as.numeric(pc_values[2]),1),"% )",sep=""), 
  # 				 title = title_name_add) + 
  # 		stat_ellipse(aes(x=V2, y=V3,group=two_var, lty=get(var2), colour=get(var1)),            level = .9, show.legend = T, type = "t", geom = "polygon", alpha = 0, inherit.aes=T) +   scale_linetype_manual(name=var2, labels=names(var2_table)[var2_table!=0], values=var2_line)
  
  ggplot(uwmpc_all, aes(x = V2, y = V3, colour=two_var, shape=two_var, fill=two_var))+
    geom_point(alpha = 1, size=2.5) + #, shape=plot_shapes, size=3, col=plot_colors) + 
    scale_color_manual(name = "", labels=names(twovar_table)[twovar_table!=0], values=var1_colors) + 
    scale_fill_manual(name = "", labels=names(twovar_table)[twovar_table!=0], values=var1_colors) + 
    scale_shape_manual(name = "", labels=names(twovar_table)[twovar_table!=0], values=var2_shape) + 
    #		scale_linetype_manual(name=var2, labels=names(var2_table)[var2_table!=0], values=c(1,2)) +
    #		scale_fill_manual(name="Legend", values=c("red","black","white")) +
    theme(axis.text=element_text(size=12), 
          panel.background = element_blank(), 
          axis.line = element_line(), 
          axis.ticks=element_line(), 
          axis.title=element_text(size=14),	
          #           title=element_text(size=16),
          legend.position=legend_position,
          plot.title = element_text(size=16, hjust=0))+
    labs(x = paste("PCo1 ( ",round(as.numeric(pc_values[1]),1),"% )",sep=""), 
         y = paste("PCo2 ( ",round(as.numeric(pc_values[2]),1),"% )",sep=""), 
         title = title_name_add) + 
    stat_ellipse(aes(x=V2, y=V3,group=two_var, lty=two_var),            level = .9, show.legend = T, type = "t", geom = "polygon", alpha = 0, inherit.aes=T) +   scale_linetype_manual(name = "", labels=names(twovar_table)[twovar_table!=0], values=var2_line)
  
}

### g legend - read it in
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


## work through this text

setwd('~/Dropbox/sequencing/Aug2019/robert/')
file_path <- "wW"
mapper_file <- "robert_mapper.tsv"

## unweighted
flies_unweighted <- read.table(paste('core-metrics-results-',file_path,'/unweighted_unifrac_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
flies_unweighted_dm <- as.dist(flies_unweighted[,2:dim(flies_unweighted)[2]])
flies_unifrac_table <- flies_unweighted %>% select(X) %>% left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% droplevels()

str(flies_unifrac_table)
table(list(flies_unifrac_table$pop)) # no duplicates
table(list(flies_unifrac_table$line)) # 5-6 reps each
table(list(flies_unifrac_table$generation)) # f m
table(list(flies_unifrac_table$popline)) #
table(list(flies_unifrac_table$Site)) #
table(list(flies_unifrac_table$Genotype)) #

## unweighted
adonis(flies_unweighted_dm ~ pop * line, strata = flies_unifrac_table$generation, flies_unifrac_table, permutations=1000)
adonis(flies_unweighted_dm ~ pop * line + generation, flies_unifrac_table, permutations=1000)



fig1_Suw_permanova <- adonis(flies_unweighted_dm ~ pop * line, strata = flies_unifrac_table$generation, flies_unifrac_table, permutations=1000)
fig1_permSU <- tableGrob(round(fig1_Suw_permanova$aov.tab, 2), theme=ttheme_minimal())
plot(fig1_permSU)

## weighted
flies_weighted <- read.table(paste('core-metrics-results-',file_path,'/weighted_unifrac_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t")
flies_weighted_dm <- as.dist(flies_weighted[,2:dim(flies_weighted)[2]])

adonis(flies_weighted_dm ~ pop * line, strata = flies_unifrac_table$generation, flies_unifrac_table, permutations = 1000)

fig1_Sw_permanova <- adonis(flies_weighted_dm ~ pop * line, strata = flies_unifrac_table$generation, flies_unifrac_table, permutations = 1000)
fig1_permSW <- tableGrob(round(fig1_Sw_permanova$aov.tab, 2), theme=ttheme_minimal())


## bray curtis
flies_bc <- read.table(paste('core-metrics-results-',file_path,'/bray_curtis_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t")
flies_bc_dm <- as.dist(flies_bc[,2:dim(flies_bc)[2]])

adonis(flies_bc_dm ~ pop * line, strata = flies_unifrac_table$generation, flies_unifrac_table, permutations = 1000)

fig1_Sb_permanova <- adonis(flies_bc_dm ~ pop * line, strata = flies_unifrac_table$generation, flies_unifrac_table, permutations = 1000)
fig1_permSB <- tableGrob(round(fig1_Sb_permanova$aov.tab, 2), theme=ttheme_minimal())

## taxon plot 
var1 = "popline"
var2 = "popline"
taxonomic_level = "order"
rpa_in_chart = 0.015
the_two_var = "twovar"
the_id = "X.SampleID"
plot_color = c("black","brown","blue","red")
plot_order = c("Rickettsiales", "Lactobacillales", "Rhodospirillales")
x_val = "popline"
read_depth = 8400
map_path = "-robert" # this is the extension of the taxonomy file
legend_position = "none"

plot_2016S <- make_taxon_plots_2var_nolegend(file_path = file_path,	map_path = map_path, var1 = var1,var2 = var2,	taxonomic_level = taxonomic_level,	rpa_in_chart = rpa_in_chart,	the_two_var = "twovar",	the_id = "X.SampleID",	plot_color = plot_color,	plot_order = plot_order,	x_val = x_val,read_depth = read_depth, mapper_file = mapper_file, legend_position = legend_position)
plot_2016S
legend_position = "bottom"
plot_fig1legend <- g_legend(make_taxon_plots_2var_nolegend(file_path = file_path,	map_path = map_path, var1 = var1,var2 = var2,	taxonomic_level = taxonomic_level,	rpa_in_chart = rpa_in_chart,	the_two_var = "twovar",	the_id = "X.SampleID",	plot_color = plot_color,	plot_order = plot_order,	x_val = x_val,read_depth = read_depth, mapper_file = mapper_file, legend_position = legend_position))

folder_name = paste0(file_path,"/unweighted_unifrac")
var1 = "pop"
var2 = "generation"
var1_colors <-c("blue","red")
var2_shape <- c(0,15)
title_name_add <- "unweighted Unifrac"
legend_position = "none"
plot_uw2016 <- pcoa_2var(folder_name=folder_name, mapper_file = mapper_file, var1=var1, var2=var2, var1_colors=var1_colors, var2_shape=var2_shape, title_name_add=title_name_add, legend_position = legend_position)
plot_uw2016

folder_name = paste0(file_path,"/weighted_unifrac")
title_name_add <- "weighted Unifrac"
plot_w2016 <- pcoa_2var(folder_name=folder_name, mapper_file = mapper_file, var1=var1, var2=var2, var1_colors=var1_colors, var2_shape=var2_shape, title_name_add=title_name_add, legend_position = legend_position)
plot_w2016

folder_name = paste0(file_path,"/bray_curtis")
title_name_add <- "Bray Curtis"
plot_bc2016 <- pcoa_2var(folder_name=folder_name, mapper_file = mapper_file, var1=var1, var2=var2, var1_colors=var1_colors, var2_shape=var2_shape, title_name_add=title_name_add, legend_position = legend_position)
plot_bc2016

legend_position = "bottom"
plot_wuwlegend <- g_legend(pcoa_2var(folder_name=folder_name, mapper_file = mapper_file, var1=var1, var2=var2, var1_colors=var1_colors, var2_shape=var2_shape, title_name_add=title_name_add, legend_position = legend_position))
plot(plot_wuwlegend)

jpeg(h=1200, w=1600, "w_pabloFigure_wWolb_nov2019.jpg", units = "px", quality = 0.9, pointsize = 1/300)
grid.arrange(	
  plot_2016S, plot_fig1legend, 
  plot_bc2016,plot_uw2016, plot_w2016,
  plot_wuwlegend,
  fig1_permSB, fig1_permSU, fig1_permSW, 
  widths = c(1,1,1),	
  heights = c(2,.3,2,.3,1.5),	
  layout_matrix=rbind(		
    c(1,1,NA),
    c(2,2,2),
    c(3,4,5),
    c(6,6,6),
    c(7,8,9))
)
dev.off()

getwd()

