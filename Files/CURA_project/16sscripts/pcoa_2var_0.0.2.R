# bc_plot <- pcoa_2var(
# 	folder_name = "Efemale/bray_curtis",
# 	var1 = "time_point",
# 	var2 = "genotype",
# 	source_var1 = c(1,2,3),
# 	var1_colors = c("red","blue","black"),
# 	var2_shape = c(21,22,23,24),
# 	title_name_add = "Bray Curtis",
# 	legend_position = "none")

# 0.0.2 - allows use of variable mapper name
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

pcoa_1var <- function(folder_name, mapper_file = mapper_file, var1, source_var1, var1_colors, title_name_add, legend_position = "bottom") {

  wfpc <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt',sep=""), sep="\t", fill=T, skip = 9,blank.lines.skip = T, header=F),-2)
  pc2 <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt', sep = ""), sep="\t", fill=T, skip = 3, header=F, colClasses = "character"),2)
  pc_values <- as.numeric(as.character(pc2[2,]))*100
  wfpc <- wfpc[,1:4]
  fut2 <- read.table(paste('core-metrics-results-',folder_name,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t", stringsAsFactors = T) %>% dplyr::select(X) %>% mutate(X=as.character(X)) %>% left_join(read.table(paste0(mapper_file),comment.char = "", header=T, sep="\t"), by=c("X"="X.SampleID"))

  uwmpc_all <- wfpc %>% inner_join(fut2, by=c("V1"="X")) #%>% mutate_(inoculated_with=as.factor(inoculated_with))

  uwmpc_all[,paste(var1)] <- as.factor(uwmpc_all[,paste(var1)])

  uwmpc_all[1,]
  var1_table <- table(list(uwmpc_all[,var1]))

  ggplot(uwmpc_all, aes(x = V2, y = V3, colour=get(var1), fill=get(var1))) +
    geom_point(alpha = 1, size=2.5, shape = 15) + #, shape=plot_shapes, size=3, col=plot_colors) +
    # scale_color_manual(name=var1, labels=names(var1_table)[var1_table!=0], values=var1_colors) +
    # scale_fill_manual(name=var1, labels=names(var1_table)[var1_table!=0], values=var1_colors) +
    # scale_shape_manual(name=var2,labels=names(var2_table)[var2_table!=0], values=var2_shape) +
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
    stat_ellipse(aes(x=V2, y=V3,group= get(var1)),            level = .9, show.legend = T, type = "t", geom = "polygon", alpha = 0, inherit.aes=T)
}

pcoa_2var_circle_2var <- function(folder_name, mapper_file, var1, var2, source_var1, var1_colors, var2_shape, title_name_add, legend_position = "bottom", var2_line = var2_line) {

	wfpc <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt',sep=""), sep="\t", fill=T, skip = 9,blank.lines.skip = T, header=F),-2)
	pc2 <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt', sep = ""), sep="\t", fill=T, skip = 3, header=F, colClasses = "character"),2)
	pc_values <- as.numeric(as.character(pc2[2,]))*100
	wfpc <- wfpc[,1:4]
	fut2 <- read.table(paste('core-metrics-results-',folder_name,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t", stringsAsFactors = T) %>% select(X) %>% mutate(X=as.character(X)) %>% left_join(read.table(mapper_file,comment.char = "", header=T), by=c("X"="X.SampleID"))

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
