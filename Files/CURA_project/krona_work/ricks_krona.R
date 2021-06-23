library(ggrepel)

taxonomy <- "https://raw.githubusercontent.com/emarsh25/chaston_lab/master/Park-Proj/Evan_stuff/CURA_project/krona_work/taxonomy.tsv"
taxonomyR <- "https://raw.githubusercontent.com/emarsh25/chaston_lab/master/Park-Proj/Evan_stuff/CURA_project/krona_work/taxonomy_forR.csv"
table <- "https://raw.githubusercontent.com/emarsh25/chaston_lab/master/Park-Proj/Evan_stuff/CURA_project/krona_work/table.from_biom.tsv"

aaa <- read.table(url(taxonomy), header = T, sep = "\t")
aaa[1,]

mt3$phylum

otu_table <- read.table(url(table), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(url(taxonomyR)), by=c("X.OTU.ID"="Feature.ID"))

## make melted OTU table
otu_table$OTU <- as.character(otu_table$X)
otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species","OTU")]
otu2$OTU <- as.character(otu2$OTU)
melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
for (i in 2:(dim(otu_table)[2]-7-1)) {
  rm(new_table)
  new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
  melted_table <- rbind(melted_table, new_table)
}

names(otu_table)

## sanity check
colnames(otu_table)
dim(melted_table)[1]/(59)==dim(new_table)[1]


## start merging
melted_table_no0 <- melted_table %>% filter(Count>-1)
mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))
mt0[1,]
# mt0$sample2 <- gsub("\\.","", mt0$Sample)
# mt0$sample2 <- gsub("_","", mt0$sample2)
# map2$sample2 <- gsub("_","",map2$id)
# map2$sample2 <- gsub("-","",map2$sample2)
# mt0[1,]
#mt1 <- mt0 %>% inner_join(map2)
mt1 <- mt0
mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
mt2$sort_order <- c(1:dim(mt2)[1])
mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
mt2$order <- paste(mt2$class, mt2$order,sep="_")
mt2$family <- paste(mt2$order, mt2$family,sep="_")
mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
mt2$species <- paste(mt2$genus, mt2$species,sep="_")
mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")

subset_column <- ""
subset_category <- ""
perc_rep <- 0.025
in_title <- ""

#make_krona <- function(mt2, subset_column, subset_category, perc_rep, in_title) {
 # mt3 <- mt2 %>% filter(get(subset_column)==subset_category)
mt3 <- mt2
mt3[1,]

mt3
mt4 <- mt3 %>% filter(Count > 0)


  jphy <- mt3 %>% group_by(phylum) %>% filter(Count > 0) %>% summarize(phylum.abun = sum(Count), sort_order=first(sort_order), color="black")
  jphy_order <- order(jphy[,"sort_order"])
  jphy <- jphy[jphy_order, , drop = FALSE]
  jphy <- jphy %>% mutate(ymax=cumsum(phylum.abun), ymin=c(0, head(ymax, n=-1)), xmin=1,xmax=2, level = "phylum", taxon=phylum, perc_phylabun = phylum.abun/sum(phylum.abun)) %>% dplyr::select(-phylum)
  jphy %>% filter(perc_phylabun > perc_rep)
  phy_colors <- jphy %>% filter(perc_phylabun>perc_rep) %>% mutate(color2 = c("blue","red")) %>% dplyr::select(taxon, color2)
  jphy2 <- jphy %>%
    full_join(phy_colors, by = "taxon") %>%
    rowwise() %>%
    mutate(color2 = ifelse(is.na(color2),"black",color2)) %>%
    ungroup()

  jclass <- mt3 %>% group_by(class) %>% filter(Count > 0) %>%  summarize(phylum.abun = sum(Count), sort_order=first(sort_order), color="black")
  jclass_order <- order(jclass[,"sort_order"])
  jclass <- jclass[jclass_order, , drop = FALSE]
  jclass <- jclass %>% mutate(ymax=cumsum(phylum.abun), ymin=c(0, head(ymax, n=-1)), xmin=2,xmax=3, level = "class", taxon=class, perc_phylabun = phylum.abun/sum(phylum.abun)) %>% dplyr::select(-class)
  jclass %>% filter(perc_phylabun > perc_rep)
  class_colors <- jclass %>% filter(perc_phylabun>perc_rep) %>% mutate(color2 = c("blue","red","red")) %>% dplyr::select(taxon, color2)
  jclass2 <- jclass %>%
    full_join(class_colors, by = "taxon") %>%
    rowwise() %>%
    mutate(color2 = ifelse(is.na(color2),"black",color2)) %>%
    ungroup()

  jorder <- mt3 %>% group_by(order) %>% filter(Count > 0) %>% summarize(phylum.abun = sum(Count), sort_order=first(sort_order), color="black")
  jorder_order <- order(jorder[,"sort_order"])
  jorder <- jorder[jorder_order, , drop = FALSE]
  jorder <- jorder %>% mutate(ymax=cumsum(phylum.abun), ymin=c(0, head(ymax, n=-1)), xmin=3,xmax=4, level = "order", taxon=order, perc_phylabun = phylum.abun/sum(phylum.abun)) %>% dplyr::select(-order)
  jorder %>% filter(perc_phylabun > perc_rep)
  order_colors <- jorder %>% filter(perc_phylabun>perc_rep) %>% mutate(color2 = c("blue","red","yellow")) %>% dplyr::select(taxon, color2)
  jorder2 <- jorder %>%
    full_join(order_colors, by = "taxon") %>%
    rowwise() %>%
    mutate(color2 = ifelse(is.na(color2),"black",color2)) %>%
    ungroup()

  jfamily <- mt3 %>% group_by(family) %>% filter(Count > 0) %>% summarize(phylum.abun = sum(Count), sort_order=first(sort_order), color="black")
  jfamily_order <- order(jfamily[,"sort_order"])
  jfamily <- jfamily[jfamily_order, , drop = FALSE]
  jfamily <- jfamily %>% mutate(ymax=cumsum(phylum.abun), ymin=c(0, head(ymax, n=-1)), xmin=4,xmax=5, level = "family", taxon=family, perc_phylabun = phylum.abun/sum(phylum.abun)) %>% dplyr::select(-family)
  jfamily %>% filter(perc_phylabun > perc_rep)
  family_colors <- jfamily %>% filter(perc_phylabun>perc_rep) %>% mutate(color2 = c("green","blue","red","yellow")) %>% dplyr::select(taxon, color2)
  jfamily2 <- jfamily %>%
    full_join(family_colors, by = "taxon") %>%
    rowwise() %>%
    mutate(color2 = ifelse(is.na(color2),"black",color2)) %>%
    ungroup()

  jgenus <- mt3 %>% group_by(genus) %>% filter(Count > 0) %>% summarize(phylum.abun = sum(Count), sort_order=first(sort_order), color="black")
  jgenus_order <- order(jgenus[,"sort_order"])
  jgenus <- jgenus[jgenus_order, , drop = FALSE]
  jgenus <- jgenus %>% mutate(ymax=cumsum(phylum.abun), ymin=c(0, head(ymax, n=-1)), xmin=5,xmax=6, level = "genus", taxon=genus, perc_phylabun = phylum.abun/sum(phylum.abun)) %>% dplyr::select(-genus)
  jgenus %>% filter(perc_phylabun > perc_rep)  %>% dplyr::select(ymax, ymin, taxon)
  genus_colors <- jgenus %>% filter(perc_phylabun>perc_rep) %>% mutate(color2 = c("green","brown","blue","red","yellow")) %>% dplyr::select(taxon, color2)
  jgenus2 <- jgenus %>%
    full_join(genus_colors, by = "taxon") %>%
    rowwise() %>%
    mutate(color2 = ifelse(is.na(color2),"black",color2)) %>%
    ungroup()

  jspecies <- mt3 %>% group_by(species) %>% filter(Count > 0) %>% summarize(sort_order=first(sort_order), phylum.abun = sum(Count), color="black")
  jspecies_order <- order(jspecies[,"sort_order"])
  jspecies <- jspecies[jspecies_order, , drop = FALSE]
  jspecies <- jspecies %>% mutate(ymax=cumsum(phylum.abun), ymin=c(0, head(ymax, n=-1)), xmin=6,xmax=7, level = "species", taxon=species, perc_phylabun = phylum.abun/sum(phylum.abun)) %>% dplyr::select(-species)
  jspecies %>% filter(perc_phylabun > perc_rep) %>% dplyr::select(ymax, ymin, taxon)
  species_colors <- jspecies %>% filter(perc_phylabun>perc_rep) %>% mutate(color2 = c("green","brown","blue","red","yellow")) %>% dplyr::select(taxon, color2)
  jspecies2 <- jspecies %>%
    full_join(species_colors, by = "taxon") %>%
    rowwise() %>%
    mutate(color2 = ifelse(is.na(color2),"black",color2)) %>%
    ungroup()

  jOTU <- mt3 %>% group_by(OTU) %>% filter(Count > 0) %>% summarize(phylum.abun = sum(Count), sort_order=first(sort_order), color="black")
  jOTU_order <- order(jOTU[,"sort_order"])
  jOTU <- jOTU[jOTU_order, , drop = FALSE]
  jOTU <- jOTU %>% mutate(ymax=cumsum(phylum.abun), ymin=c(0, head(ymax, n=-1)), xmin=7,xmax=8, level = "OTU", taxon=OTU, perc_phylabun = phylum.abun/sum(phylum.abun)) %>% dplyr::select(-OTU)
  jOTU %>% filter(perc_phylabun > perc_rep)
  OTU_colors <- jOTU %>% filter(perc_phylabun>perc_rep) %>% mutate(color2 = c("blue","magenta","green","brown","black","grey","orange","red")) %>% dplyr::select(taxon, color2)
  OTU_colors <- jOTU %>% filter(perc_phylabun>perc_rep) %>% mutate(color2 = c("black")) %>% dplyr::select(taxon, color2)
  jOTU2 <- jOTU %>%
    full_join(OTU_colors, by = "taxon") %>%
    rowwise() %>%
    mutate(color2 = ifelse(is.na(color2),"black",color2)) %>%
    ungroup()

  jphy$taxon <- as.character(jphy$taxon); jclass$taxon <- as.character(jclass$taxon); jorder$taxon <- as.character(jorder$taxon); jfamily$taxon <- as.character(jfamily$taxon); jgenus$taxon <- as.character(jgenus$taxon);  jspecies$taxon <- as.character(jspecies$taxon); jOTU$taxon <- as.character(jOTU$taxon)

  krona <- bind_rows(jphy2,jclass2,jorder2,jfamily2,jgenus2,jspecies2,jOTU2)
  str(krona)
  col <- krona$color2
  names(col) <- krona$taxon


  krona_labels <- c(order_colors$taxon, family_colors$taxon,genus_colors$taxon,species_colors$taxon,OTU_colors$taxon)


  labels_df <- krona %>%
    mutate(label_ypos=ifelse(taxon %in% krona_labels,
                             (ymax+ymin)/2,
                             NA),
           label_xpos=ifelse(taxon %in% krona_labels,
                             (xmin+xmax)/2,
                             NA),
           taxon3 = gsub("^.*_","", taxon),
           taxon3 = gsub("^.*_","", taxon)
           ) %>%
    mutate(taxon2=ifelse((level=="OTU"|taxon3==""),
                         NA,
                         taxon3)) %>%
    dplyr::select(taxon, taxon2, label_ypos,label_xpos,phylum.abun)
  labels_df

  label_ypos <-labels_df$label_ypos; names(label_ypos) <- as.character(labels_df$taxon)
  label_xpos <- labels_df$label_xpos; names(label_xpos) <- as.character(labels_df$taxon)
  label_taxa <- as.character(labels_df$taxon2); names(label_taxa) <- as.character(labels_df$taxon)


  in_title <- "krona for PD datasets"

  library(ggrepel)
  ggplot(krona, aes(fill=color2,alpha = level,ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin)) +
    # set xlim based on how many circles i plan on
    geom_rect(colour="grey30") +  coord_polar(theta="y") + xlim(c(0, 8)) +
    theme_bw() + theme(legend.position = "none") +theme(panel.grid=element_blank()) + theme(axis.text=element_blank()) +  theme(axis.ticks=element_blank(),  panel.border = element_blank()) +
    scale_fill_identity() +
    scale_alpha_discrete(range = c(1,1)) +
    labs(title=in_title) +
    #geom_label(aes(x=label_xpos, y=label_ypos), label=label_taxa, na.rm = T, inherit.aes = F, position = position_jitter(height = 1.5, width = 1.5))
    geom_label_repel(aes(x=label_xpos, y=label_ypos), fill = "white", alpha = 0.75, label=label_taxa, na.rm = T, inherit.aes = F, segment.colour = "white", box.padding = 0, label.padding = 0.1, min.segment.length = 0, direction = "y")


table(label_xpos)
table(label_ypos)

