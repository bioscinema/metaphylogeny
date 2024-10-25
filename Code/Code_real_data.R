
## Visualization for 16S rRNA data processed by LotuS 2 (/phyloseq_objetcs/lotus_phy.RData)

## Load phyloseq object
load("/phyloseq_objects/lotus_phy.RData")

## Data cleaning
mytax <- data.frame(physeq@tax_table)
mytax[mytax == "?"] = "unknown"
colnames(mytax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_table(physeq) <- tax_table(as.matrix(mytax))

## Plot taxonomy tree
ps <- fix_duplicate_tax(ps1)
ps1 <- subset_taxa(physeq,Kingdom=="Bacteria")
mytax <- data.frame(ps1@tax_table)
tr1 <- phy_to_tax(ps)
p <- taxview(ps,tr1)
p
ggsave(filename = "~/lotus2.pdf", scale = 1.5,width = 10,height = 13)

## Visualziationfor 16S rRNA data processed by QIIME 2 (/phyloseq_objetcs/qiime_phy.RData)
## Load phyloseq object
load("/phyloseq_objects/qiime_phy.RData")

## Data cleaning
tax1 <- data.frame(physeq@tax_table)
colnames(tax1) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax1[] <- lapply(tax1, function(x) gsub("[a-zA-Z]__", "", x))
tax1[is.na(tax1)]="unknown"
rows_to_keep <- apply(tax1[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 
                      1, 
                      function(x) !all(x == "unknown"))

# Filter the tax_table to keep only the rows where not all ranks are "unknown"
tax1 <- tax1[-rows_all_unknown, ]

tax_table(physeq) <- tax_table(as.matrix(tax1))

## Plot taxonomy tree
ps <- fix_duplicate_tax(physeq)
tr1 <- phy_to_tax(ps)
p<- taxview(ps,tr1)
p
ggsave(filename = "~/qiime2.pdf", scale = 1.5,height = 13,width = 10)


## Visualization for WGS data processed by Woltka (/phyloseq_objects/woltka_phy.RData)
## Load phyloseq data
load("/phyloseq_objects/woltka_phy.RData")

## Data cleaning
ps2 <- subset_taxa(physeq, kingdom=="Bacteria")
mytax <- data.frame(ps2@tax_table)
colnames(mytax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_table(ps2) <- tax_table(as.matrix(mytax))

## Plot taxonomy tree
ps3 <- fix_duplicate_tax(ps2)
tr2 <- phy_to_tax(ps3)
p1 <- taxview(ps3,tr2)
p1
ggsave(filename = "~/woltka.pdf",width = 10,height = 8,bg="white")


## Visualization for WGS data processed by MetaPhlAn 4 (/phyloseq_objects/metaphlan_phy.RData)
## Load phyloseq data
load("/phyloseq_objects/woltka_phy.RData")

## Data cleaning
mytax <- data.frame(ps@tax_table)
columns_to_modify <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
mytax[columns_to_modify] <- lapply(mytax[columns_to_modify], function(x) gsub("^[A-Za-z]+__", "", x))
mytax$Phylum <- gsub("^p__", "", mytax$Phylum)
tax_table(ps) <- tax_table(as.matrix(mytax))

## Plot taxonomy tree
ps1 <- fix_duplicate_tax(ps)
tr1 <- phy_to_tax(ps1)
p <- taxview(ps1,tr1)
p
ggsave(filename = "~/metaphlan.pdf",width = 10,height = 8,bg="white")