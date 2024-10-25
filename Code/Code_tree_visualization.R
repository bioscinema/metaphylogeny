library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(TDbook)

# Code for generate phylogenetic tree figure
## Phylogenetic tree

# load data from TDbook, including tree_hmptree, 
# df_tippoint (the abundance and types of microbes),
# df_ring_heatmap (the abundance of microbes at different body sites),
# and df_barplot_attr (the abundance of microbes of greatest prevalence)
tree <- tree_hmptree
dat1 <- df_tippoint
dat1 <- dat1 %>% filter(!Phylum=="NA")
dat2 <- df_ring_heatmap
dat3 <- df_barplot_attr

# adjust the order
dat2$Sites <- factor(dat2$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)","Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))
dat3$Sites <- factor(dat3$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)", "Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))

nodeids <- nodeid(tree, tree$node.label[nchar(tree$node.label)>4])
nodedf <- data.frame(node=nodeids)
nodelab <- gsub("[\\.0-9]", "", tree$node.label[nchar(tree$node.label)>4])
# The layers of clade and hightlight
poslist <- c(1.6, 1.4, 1.6, 0.8, 0.1, 0.25, 1.6, 1.6, 1.2, 0.4,
             1.2, 1.8, 0.3, 0.8, 0.4, 0.3, 0.4, 0.4, 0.4, 0.6,
             0.3, 0.4, 0.3)
labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)

# The circular layout tree.
p <- ggtree(tree, layout="fan", size=0.15, open.angle=5) +
  geom_hilight(data=nodedf, mapping=aes(node=node),
               extendto=6.8, alpha=0.3, fill="grey", color="grey50",
               size=0.05) 
# geom_cladelab(data=labdf, 
#               mapping=aes(node=node, 
#                           label=label,
#                           offset.text=pos),
#               hjust=0.5,
#               angle="auto",
#               barsize=NA,
#               horizontal=FALSE, 
#               fontsize=1.4,
#               fontface="italic"
# )

p <- p %<+% dat1 + geom_star(
  mapping=aes(fill=Phylum, starshape=Type, size=Size),
  position="identity",starstroke=0.1) +
  scale_fill_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
                             "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                             "#EE6A50","#8DEEEE", "#006400","#800000",
                             "#B0171F","#191970")
                    ,na.tra)+
  scale_starshape_manual(values=c(15, 1),
                         guide="none",
                         na.translate=FALSE)+
  scale_size_continuous(range = c(1, 2.5),
                        guide = "none")+
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
  )
p
p <- ggtree(tree, layout="fan", size=0.15, open.angle=5) +
  geom_hilight(data=nodedf, mapping=aes(node=node),
               extendto=6.8, alpha=0.3, fill="grey", color="grey50",size=0.05)

p <- p %<+% dat1 + geom_star(
  mapping=aes(fill=Phylum, starshape=Type, size=Size),
  position="identity", starstroke=0.1) +
  scale_fill_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
                             "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                             "#EE6A50","#8DEEEE", "#006400","#800000",
                             "#B0171F","#191970"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(shape=22)),na.translate=FALSE) +
  scale_starshape_manual(values=c(15, 1),
                         guide="none") +
  scale_size_continuous(range = c(1, 2.5),
                        guide="none") +
  guides(
    fill = guide_legend(override.aes = list(shape = 25, size = 5))
  ) +
  theme(legend.position = "right",
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))

p <- p %<+% dat1 + 
  geom_point(mapping = aes(fill = Phylum, size = Size), 
             position = "identity", shape = 21, stroke = 0.1, show.legend = FALSE) +
  scale_fill_manual(values = c("#FFC125","#87CEFA","#7B68EE","#808080",
                               "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                               "#EE6A50","#8DEEEE", "#006400","#800000",
                               "#B0171F","#191970")) +
  scale_size_continuous(range = c(1, 2.5), guide = "none") +
  theme(legend.position = "none")

print(p)


## Taxonomy tree visualization
# Installing BAZE package
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("bioscinema/BAZE")

mytax <- TDbook::df_alltax_info
table(mytax$Phylum)

num_features <- 1351
num_patients <- 50

# Create a matrix with all values set to 1
otu_table <- matrix(1, nrow = num_features, ncol = num_patients)

# Convert the matrix to a data frame
otu_df <- as.data.frame(otu_table)

# Assign feature names (OTU IDs) and patient IDs
rownames(otu_df) <- paste0("OTU", 1:num_features)
colnames(otu_df) <- paste0("Patient", 1:num_patients)
rownames(mytax) <- rownames(otu_df)
physeq <- phyloseq(otu_table(otu_df,taxa_are_rows=TRUE),tax_table(as.matrix(mytax)))
phy1 <- fix_duplicate_tax(physeq)
tree1 <- phy_to_tax(phy1)

taxview <- function(ps, tree, branch_thickness=0.5, layout='circular', level="Phylum") {
  # Extract taxonomy table
  tax <- as.data.frame(ps@tax_table)
  
  # Generate a color palette for the taxonomic levels
  taxon_colors <- c("#FFC125","#87CEFA","#7B68EE","#808080",
                    "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                    "#EE6A50","#8DEEEE", "#006400","#800000",
                    "#B0171F","#191970")
  names(taxon_colors) <- sort(unique(tax[[level]]))
  
  # Create a mapping from tips to their taxonomy level
  tip_labels <- data.frame(label = rownames(tax), taxon = tax[[level]])
  
  # Create the ggtree plot
  p <- ggtree(tree, size=branch_thickness, layout=layout)
  
  # Highlight each taxon clade based on internal node labels and add clade labels
  for (taxon in unique(tip_labels$taxon)) {
    node_label <- paste0(tolower(substr(level, 1, 1)), "__", taxon)
    if (node_label %in% tree@phylo[["node.label"]]) {
      node_index <- which(tree@phylo[["node.label"]] == node_label) + length(tree@phylo[["tip.label"]])
      p <- p + geom_hilight(node = node_index, fill = taxon_colors[taxon], alpha = 0.3)

  
  # Add a dummy dataframe to create the legend
  legend_df <- data.frame(taxon = names(taxon_colors), fill = taxon_colors)
  
  # Add the legend manually
  p <- p + geom_point(data = legend_df, aes(x = 0, y = 0, fill = taxon), size = 0) +
    scale_fill_manual(values = taxon_colors, name = level) +
    guides(fill = guide_legend(override.aes = list(size = 5, shape = 21))) +
    theme(legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.key = element_rect(fill = "white"))
  
  # Return the plot
  return(p)
}

p12 <- taxview(phy1,tree1)
p12 <- rotate_tree(p12,155)

p<- p+theme(plot.title = element_text(size = 30, face = "bold"),
            legend.title=element_text(size=30),
            legend.text=element_text(size=24))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # axis.line.x = element_line(color="black", size = 2),
        # axis.line.y = element_line(color="black", size = 2),
        axis.text=element_text(size=24,face="bold"),
        axis.title=element_text(size=30,face="bold"))
p12 <- p12+theme(plot.title = element_text(size = 30, face = "bold"),
                 legend.title=element_text(size=30),
                 legend.text=element_text(size=24),
                 legend.key = element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # axis.line.x = element_line(color="black", size = 2),
        # axis.line.y = element_line(color="black", size = 2),
        axis.text=element_text(size=24,face="bold"),
        axis.title=element_text(size=30,face="bold"))
p12
library(cowplot)
combined_plot <- plot_grid(p,p12, labels = c("A","B"), ncol = 2, align = "hv",label_size = 10,label_y = 0.8)
print(combined_plot)