### Import necessary packages

install.packages("BiocManager")
install.packages("msa")
install.packages("Biostrings")
install.packages("ape")
install.packages("cluster")
install.packages("factoextra")

library(msa)
library(Biostrings)
library(ape)
library(cluster)
library(ggplot2)
library(factoextra)

################################################################################
# I. Read in and fill out accession list
acclist <- NULL
tryCatch( acclist <- as.data.frame(read.table(
  choose.files(filters = Filters[c("zip", "All"),]),header=F,sep="\t")),error = function(e)
      print("Could not read species & accession list!"))

if ((dim(acclist)[1] < 2) | (dim(acclist)[2] != 2)) {
  print("Improper format of accession file!")
  exit()
}

# give acclist data frame column names
colnames(acclist) <- c("species","accession")
# add sequence column to acclist
acclist$sequence <- NA
acclist$sequence_length <- 0

# Exit if there are species with multiple accessions
if (length(as.list(names(which(table(acclist$species)>1)))) > 0) {
  print("Multiple accessions detected for at least one species! Please edit accession file!")
  exit()
}

# number of species
n_species <- dim(acclist)[1]
# list of all products: genes, tRNA, rRNA
products <- c()
# maximum length of DNA
max_sequence_length <- 0
# maximum number of genes
max_n_genes <- 200 # > 150, since plastids have ~150 genes

# genemx is a matrix containing info for n species and all genes
# associated with it are genemx_start and genemx_end
gene_matrix <- matrix(nrow=n_species,ncol=max_n_genes)
gene_start_matrix <- matrix(nrow=n_species,ncol=max_n_genes)
gene_end_matrix <- matrix(nrow=n_species,ncol=max_n_genes)
rownames(gene_matrix) <- acclist$species
rownames(gene_start_matrix) <- acclist$species
rownames(gene_end_matrix) <- acclist$species

# set counter to zero
c <- 0
# go through list, download GenBank annotation, package into variables
for (accession in acclist$accession) {
  c <- c + 1
  species <- acclist$species[c]
  
  print(paste("Downloading ",accession,"...",sep=""))
  tryCatch( annot <- getAnnotationsGenBank(accession),
            error = function(e)
              print("getAnnotationsGenBank unsuccessful!"))
  # sleep so there aren't too many requests all at once
  Sys.sleep(0.5)
  
  # put gene/RNA products into prod
  products <- unique(c(sort(products),sort(annot$product)))
  # get end of last element
  annotation_end <- annot$end[length(annot$end)]
  # if longer than longest position, overwrite it
  if (annotation_end > max_sequence_length) max_sequence_length <- annotation_end
  # get number of genes/RNAs
  num_genes <- length(annot$product)
  # if more than largest number, overwrite
  if (num_genes > max_n_genes) max_n_genes <- num_genes
  
  # go through all of species' genes/RNAs (products)
  for (i in 1:num_genes) {
    if (!is.na(annot$product[i])) {
      # add product elements to genemx, genemx_start, genemx_end
      gene_matrix[species,i] = annot$product[i]
      gene_start_matrix[species,i] = annot$start[i]
      gene_end_matrix[species,i] = annot$end[i]
    }
  }
  # list of genes/RNAs for species
  print(paste("Writing elements for ",species,"...",sep=""))
  species_genes_rnas <- paste(species," |",sep="")
  # go through all products and add them to sp_genes_rnas and write them into a file
  for (pr in annot$product) {
    element <- pr
    if (is.na(pr)) {pr <- 'misc_feature'}
    species_genes_rnas <- paste(species_genes_rnas,pr,sep="\t")
  }
  write(species_genes_rnas, file=paste("species_products.txt"), sep="\t", append = T)

  # get sequence
  GBi <- read.GenBank(accession,as.character=T)
  Sys.sleep(1)
  # the name of the column in GBi is some accession umber, so we rename it
  names(GBi) <- "id"
  sequence_length <- length(GBi$id)
  sequence <- paste(GBi$id,sep="",collapse="")
  acclist[acclist$accession==accession,]$sequence <- sequence
  acclist[acclist$accession==accession,]$sequence_length <- sequence_length
}

################################################################################
# II. Calculate matrixes and write them to file
# Distance matrix and sequence similarity array

# distance matrix
distance_matrix <- matrix(nrow = n_species, ncol = n_species)
# sequence similarity matrix
sequence_similarity_matrix <- matrix(nrow = n_species, ncol = n_species)
# column and row names for both are species, since both are square matrixes
colnames(distance_matrix) <- acclist$species
rownames(distance_matrix) <- acclist$species
colnames(sequence_similarity_matrix) <- acclist$species
rownames(sequence_similarity_matrix) <- acclist$species

# go through all species
for (i in seq(1,n_species)) {
  # get first species (i)
  species_i <- acclist[i,1]
  print(paste("Analyzinging ",species_i,"...",sep=""))
  # get its accession number
  accession_i <- acclist[i,2]
  # get its sequence
  sequence_i <- acclist[acclist$accession == accession_i,]$sequence
  
  for (j in seq(i,n_species)) {
    # get data and sequence for species j
    species_j <- acclist[j,1]
    accession_j <- acclist[j,2]
    # get its sequence
    sequence_j <- acclist[acclist$accession == accession_j,]$sequence
    
    # % similarity with percent identity = pid function
    percent_similarity <- pid(pairwiseAlignment(sequence_i,sequence_j,type="global"),"PID2")/100
    # for pixels in diagonal
    if (i == j) {
      # distance is 0, and similarity is 100%
      distance_matrix[species_i,species_j] = 0
      sequence_similarity_matrix[species_i,species_j] = 1
    }
    # for all other species pairs
    if (i < j) {
      sequence_similarity_matrix[species_i,species_j] = percent_similarity
      sequence_similarity_matrix[species_j,species_i] = percent_similarity
      # initialize distance to 0
      d <- 0
      # get gene products for species i, and exlude NA elements
      x <- gene_matrix[species_i,]
      x <- x[!is.na(x)]
      # same for species j
      y <- gene_matrix[species_j,]
      y <- y[!is.na(y)]
      # get common elements to species i and j
      xy <- intersect(x,y)
      # go through all common elements (products/genes/RNAs)
      for (k in seq(1,length(xy))) {
        # an element might be present more than once (duplicate),
        # thus take the mean index for species i and species j
        d <- d + abs(mean(which(x==xy[k])) - mean(which(y==xy[k])))
      }
      # Add mean distance of common elements and also add number of elements unique to species i
      # and elements unique to species j
      # by subtracting number of common elements from number of elements in species i
      # and in species j
      d <- d/length(xy) + (length(x) - length(xy)) + (length(y) - length(xy))
      distance_matrix[species_i,species_j] <- d
      distance_matrix[species_j,species_i] <- d
    }
  }
}

# Create similarity heatmap and write matrix
# gene order matrix is 1 minus normalized distance matrix
# normalize so you get a scale from 0 to 1,
# but subtract norm.distance from 1 to get gene order similarity
similarity_matrix <- 1 - distance_matrix/max(distance_matrix)
# write matrix output
write.table(distance_matrix, file="distance_matrix.txt", sep="\t", quote = F)
write.table(similarity_matrix, file="similarity_matrix.txt", sep="\t", quote = F)
combined_matrix = (similarity_matrix+sequence_similarity_matrix)/2
write.table(combined_matrix, file="combined_matrix.txt", sep="\t", quote = F)

################################################################################
# III. Create Silhouette plot, heat map and clusters/stats file for
# sequence similarity, gene order and combined matrixes

# General parameters
myBreaks <- c(seq(0,1,by=0.01))
cexx <- length(species)/(6*length(species))
ceyy <- cexx
clr <- colorRampPalette(c("white","yellow","orange","red"))(100)
clusmeth="ward.D2" # ward.D ward.D2 single median average mcquitty complete centroid

# A. Sequence similarity matrix
res_sequence_similarity <- get_clust_tendency(sequence_similarity_matrix, n = nrow(sequence_similarity_matrix)-1, graph = FALSE)
hopkins_sequence_similarity <- res_sequence_similarity$hopkins_stat

# A.1. Silhouette plot
jpeg("Silhouette_sequence_similarity.jpg")
fviz_nbclust(sequence_similarity_matrix, pam, method = "silhouette", k.max=n_species-1) + theme_classic()
dev.off()
silhouette_sequence_similarity <- fviz_nbclust(sequence_similarity_matrix, pam, method = "silhouette", k.max=n_species-1) + theme_classic()

# max clusters gene order similarity
k_sequence_similarity <- which(silhouette_sequence_similarity$data$y == max(silhouette_sequence_similarity$data$y))

# A.2. Heat map
# normalize
sequence_similarity_matrix_normalized <- (sequence_similarity_matrix - min(sequence_similarity_matrix))/(max(sequence_similarity_matrix) - min(sequence_similarity_matrix))

# make heat map
sequence_similarity_heatmap_name="sequence_similarity_heatmap.jpg"
jpeg(filename = sequence_similarity_heatmap_name, height=3000, width=3000, units = "px", res=300)
h <- heatmap(sequence_similarity_matrix_normalized, symkey =F, symbreaks=F, scale="none", dendrogram = F, Rowv=F, Colv=F,col = clr, breaks = myBreaks, border_color=NA, na.color="white", margin = c(10,10),
             cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5,1), hclustfun = function(d) hclust(d, method=clusmeth),
             labCol=as.expression(lapply(colnames(sequence_similarity_matrix_normalized), function(a) bquote(italic(.(a))))),labRow=as.expression(lapply(rownames(sequence_similarity_matrix_normalized), function(a) bquote(italic(.(a))))))
invisible(dev.off())

# A.3. clustering and stats
row.clusters = hclust(dist(sequence_similarity_matrix),method=clusmeth)
ctk_sequence_similarity <- cutree(row.clusters,k=k_sequence_similarity)
filename=paste("clusters_sequence_similarity.txt")
write.table(ctk_sequence_similarity, file=filename, col.names=F, quote=F, sep="\t")

# get cluster stats for sequence similarity
header = "cluster\tno. species\tMean DNA length±sd\tmin\tmean\tmax\tSEM\tp-value\tneglog"
write(header, file=paste("statistics_sequence_similarity.txt"), sep="\t", append=T)

cluster_sizes <- table(ctk_sequence_similarity)
non_groups = c()

for (n_cluster in 1:k_sequence_similarity) {
  csize = cluster_sizes[n_cluster]
  if (csize >= 2) {
    m1 = as.matrix(sequence_similarity_matrix[ctk_sequence_similarity == n_cluster,ctk_sequence_similarity == n_cluster])
    
    x = m1[upper.tri(m1)]
    xx = as.numeric(as.list(x))
    nn = dim(m1)[1]
    
    m2 = as.matrix(cbind(sequence_similarity_matrix[ctk_sequence_similarity != n_cluster,ctk_sequence_similarity == n_cluster],t(sequence_similarity_matrix[ctk_sequence_similarity == n_cluster,ctk_sequence_similarity != n_cluster])))
    m2b = m2[!duplicated(colnames(m2))]
    non_groups = c(non_groups, m2b)
    
    mean_dna_length <- sprintf("%.3f",mean(acclist[ctk_sequence_similarity == n_cluster,]$sequence_length))
    sd_dna_length <- sprintf("%.3f",sd(acclist[ctk_sequence_similarity == n_cluster,]$sequence_length))
    
    if (csize>= 3) {
      t = t.test(x,m2b)
    } else if (csize == 2) {
      t = t.test(m2b,mu=x)
    }
    p_value = t$p.value
    neglog = -log10(p_value)
    min = min(x)
    max = max(x)
    
    mean2 = sprintf("%.3f", mean(x))
    sd2 = sprintf("%.3f", sd(x)/sqrt(csize))
    min2 = sprintf("%.3f", min)
    max2 = sprintf("%.3f", max)
    p_value2 = sprintf("%.3f", p_value)
    neglog2 = sprintf("%.3f", neglog)
    
    statistics = paste(n_cluster, nn, paste(mean_dna_length,"±",sd_dna_length,sep=""), min2, mean2, max2, sd2, p_value, neglog2, sep="\t")
    statistics2 = gsub("\n","\t",statistics)
    write(statistics, file="statistics_sequence_similarity.txt", sep="\t", append = T)
  }
}

# B. Gene order similarity matrix
res_similarity <- get_clust_tendency(similarity_matrix, n = nrow(similarity_matrix)-1, graph = FALSE)
hopkins_similarity <- res_similarity$hopkins_stat

# B.1. Silhouette plot
jpeg("Silhouette_similarity.jpg")
fviz_nbclust(similarity_matrix, pam, method = "silhouette", k.max=n_species-1) + theme_classic()
dev.off()
silhouette_similarity <- fviz_nbclust(similarity_matrix, pam, method = "silhouette", k.max=n_species-1) + theme_classic()

# max clusters gene order similarity
k_similarity <- which(silhouette_similarity$data$y == max(silhouette_similarity$data$y))

# B.2. Heat map
# normalize
similarity_matrix_normalized <- (similarity_matrix - min(similarity_matrix))/(max(similarity_matrix) - min(similarity_matrix))

# make heat map
similarity_heatmap_name="similarity_heatmap.jpg"
jpeg(filename = similarity_heatmap_name, height=3000, width=3000, units = "px", res=300)
h <- heatmap(similarity_matrix_normalized, symkey =F, symbreaks=F, scale="none", dendrogram = F, Rowv=F, Colv=F,col = clr, breaks = myBreaks, border_color=NA, na.color="white", margin = c(10,10),
             cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5,1), hclustfun = function(d) hclust(d, method=clusmeth),
             labCol=as.expression(lapply(colnames(similarity_matrix_normalized), function(a) bquote(italic(.(a))))),labRow=as.expression(lapply(rownames(similarity_matrix_normalized), function(a) bquote(italic(.(a))))))
invisible(dev.off())

# B.3. clustering and stats
row.clusters = hclust(dist(similarity_matrix),method=clusmeth)
ctk_similarity <- cutree(row.clusters,k=k_similarity)
filename=paste("clusters_similarity.txt")
write.table(ctk_similarity, file=filename, col.names=F, quote=F, sep="\t")

# get cluster stats for sequence similarity
header = "cluster\tno. species\tMean DNA length±sd\tmin\tmean\tmax\tSEM\tp-value\tneglog"
write(header, file=paste("statistics_similarity.txt"), sep="\t", append=T)

cluster_sizes <- table(ctk_similarity)
non_groups = c()

for (n_cluster in 1:k_similarity) {
  csize = cluster_sizes[n_cluster]
  if (csize >= 2) {
    m1 = as.matrix(similarity_matrix[ctk_similarity == n_cluster,ctk_similarity == n_cluster])
    
    x = m1[upper.tri(m1)]
    xx = as.numeric(as.list(x))
    nn = dim(m1)[1]
    
    m2 = as.matrix(cbind(similarity_matrix[ctk_similarity != n_cluster,ctk_similarity == n_cluster],t(sequence_similarity_matrix[ctk_similarity == n_cluster,ctk_similarity != n_cluster])))
    m2b = m2[!duplicated(colnames(m2))]
    non_groups = c(non_groups, m2b)
    
    mean_dna_length <- sprintf("%.3f",mean(acclist[ctk_similarity == n_cluster,]$sequence_length))
    sd_dna_length <- sprintf("%.3f",sd(acclist[ctk_similarity == n_cluster,]$sequence_length))
    
    if (csize>= 3) {
      t = t.test(x,m2b)
    } else if (csize == 2) {
      t = t.test(m2b,mu=x)
    }
    p_value = t$p.value
    neglog = -log10(p_value)
    min = min(x)
    max = max(x)
    
    mean2 = sprintf("%.3f", mean(x))
    sd2 = sprintf("%.3f", sd(x)/sqrt(csize))
    min2 = sprintf("%.3f", min)
    max2 = sprintf("%.3f", max)
    p_value2 = sprintf("%.3f", p_value)
    neglog2 = sprintf("%.3f", neglog)
    
    statistics = paste(n_cluster, nn, paste(mean_dna_length,"±",sd_dna_length,sep=""), min2, mean2, max2, sd2, p_value, neglog2, sep="\t")
    statistics2 = gsub("\n","\t",statistics)
    write(statistics, file="statistics_similarity.txt", sep="\t", append = T)
  }
}

# C. Combined matrix
res_combined <- get_clust_tendency(combined_matrix, n = nrow(combined_matrix)-1, graph = FALSE)
hopkins_combined <- res_combined$hopkins_stat

# C.1. Silhouette plot
jpeg("Silhouette_combined.jpg")
fviz_nbclust(combined_matrix, pam, method = "silhouette", k.max=n_species-1) + theme_classic()
dev.off()
silhouette_combined <- fviz_nbclust(combined_matrix, pam, method = "silhouette", k.max=n_species-1) + theme_classic()

# max clusters gene order similarity
k_combined <- which(silhouette_combined$data$y == max(silhouette_combined$data$y))

# C.2. Heat map
# normalize
combined_matrix_normalized <- (combined_matrix - min(combined_matrix))/(max(combined_matrix) - min(combined_matrix))

# make heat map
combined_heatmap_name="combined_heatmap.jpg"
jpeg(filename = combined_heatmap_name, height=3000, width=3000, units = "px", res=300)
h <- heatmap(combined_matrix_normalized, symkey =F, symbreaks=F, scale="none", dendrogram = F, Rowv=F, Colv=F,col = clr, breaks = myBreaks, border_color=NA, na.color="white", margin = c(10,10),
             cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5,1), hclustfun = function(d) hclust(d, method=clusmeth),
             labCol=as.expression(lapply(colnames(combined_matrix_normalized), function(a) bquote(italic(.(a))))),labRow=as.expression(lapply(rownames(combined_matrix_normalized), function(a) bquote(italic(.(a))))))
invisible(dev.off())

# C.3. clustering and stats
row.clusters = hclust(dist(combined_matrix),method=clusmeth)
ctk_combined <- cutree(row.clusters,k=k_combined)
filename=paste("clusters_combined.txt")
write.table(ctk_combined, file=filename, col.names=F, quote=F, sep="\t")

# get cluster stats for sequence similarity
header = "cluster\tno. species\tMean DNA length±sd\tmin\tmean\tmax\tSEM\tp-value\tneglog"
write(header, file=paste("statistics_combined.txt"), sep="\t", append=T)

cluster_sizes <- table(ctk_sequence_similarity)
non_groups = c()

for (n_cluster in 1:k_combined) {
  csize = cluster_sizes[n_cluster]
  if (csize >= 2) {
    m1 = as.matrix(combined_matrix[ctk_combined == n_cluster,ctk_combined == n_cluster])
    
    x = m1[upper.tri(m1)]
    xx = as.numeric(as.list(x))
    nn = dim(m1)[1]
    
    m2 = as.matrix(cbind(combined_matrix[ctk_combined != n_cluster,ctk_combined == n_cluster],t(combined_matrix[ctk_combined == n_cluster,ctk_combined != n_cluster])))
    m2b = m2[!duplicated(colnames(m2))]
    non_groups = c(non_groups, m2b)
    
    mean_dna_length <- sprintf("%.3f",mean(acclist[ctk_combined == n_cluster,]$sequence_length))
    sd_dna_length <- sprintf("%.3f",sd(acclist[ctk_combined == n_cluster,]$sequence_length))
    
    if (csize>= 3) {
      t = t.test(x,m2b)
    } else if (csize == 2) {
      t = t.test(m2b,mu=x)
    }
    p_value = t$p.value
    neglog = -log10(p_value)
    min = min(x)
    max = max(x)
    
    mean2 = sprintf("%.3f", mean(x))
    sd2 = sprintf("%.3f", sd(x)/sqrt(csize))
    min2 = sprintf("%.3f", min)
    max2 = sprintf("%.3f", max)
    p_value2 = sprintf("%.3f", p_value)
    neglog2 = sprintf("%.3f", neglog)
    
    statistics = paste(n_cluster, nn, paste(mean_dna_length,"±",sd_dna_length,sep=""), min2, mean2, max2, sd2, p_value, neglog2, sep="\t")
    statistics2 = gsub("\n","\t",statistics)
    write(statistics, file="statistics_combined.txt", sep="\t", append = T)
  }
}

### Write stats
write(date(), file="results.txt", sep="\t", append=T)
write(paste("Number of species: ",n_species,sep=""), file="results.txt", append=T)
write(paste("Hopkins clustering stat (sequence similarity): ",round(hopkins_sequence_similarity,4),sep=""), file="results.txt", append=T)
write(paste("Hopkins clustering stat (gene order similarity): ",round(hopkins_similarity,4),sep=""), file="results.txt", append=T)
write(paste("Hopkins clustering stat (combined): ",round(hopkins_combined,4),sep=""), file="results.txt", append=T)
write(paste("Number of optimal clusters (sequence similarity): ",k_sequence_similarity,sep=""), file="results.txt", append=T)
write(paste("Number of optimal clusters (gene order similarity): ",k_similarity,sep=""), file="results.txt", append=T)
write(paste("Number of optimal clusters (combined): ",k_combined,sep=""), file="results.txt", append=T)

################################################################################
# IV. Determine parameters for gene map
# height of gene order for each species
offset_gene_height <- 250
# length of stretch of species name
species_name_offset <- 5000
# length of left offset
left_offset <- 250
# total height of figure
figure_height <- n_species * offset_gene_height + k_combined * offset_gene_height
min_figure_height <- 3300
if (figure_height < min_figure_height) {
  figure_height <- min_figure_height
}

# width of figure: species name offset for species names,
# plus longest genome plus 250 padding for sequence length
# if more than 100 products, then make figure taller
n_legend_columns <- 2
if (length(products) > 100) {
  figure_height <- 10000
  n_legend_columns <- 3
}
max_sequence_length <- max(acclist$sequence_length)
figure_width <- max_sequence_length + species_name_offset + left_offset
max_figure_width <- 50000
# in case plastids are being studied
shrink_factor <- 1
if (max_sequence_length > max_figure_width) {
  shrink_factor <- 0.9*max_figure_width/(max_sequence_length+species_name_offset+left_offset)
}
if (figure_width > max_figure_width) {
  figure_width <- max_figure_width
}

### Draw gene order figure ###
acclist$cluster <- ctk_combined[acclist$species]
acclist <- acclist[order(acclist$species),]

c <- 0
jpeg("gene_order_map_per_cluster.jpg",height=figure_height,width=figure_width)
par(mfrow=c(1,2))

plot(c(1, figure_width), c(1, figure_height), type= "n", xlab = "", ylab = "", axes=T)
par(mar=c(0,0,0,0))
# go through each cluster
for (k in 1:k_combined) {
  c <- c + 1
  # make a mini acclist per cluster
  cluster_accession_list <- acclist[acclist$cluster == k,]
  text(1,figure_height-(c-0.5)*offset_gene_height,paste("cluster",k),cex=11,adj=0,font.lab=2)
  # go through each accession belonging to that cluster
  ck <- 0
  for (acc in cluster_accession_list$accession) {
    ck <- ck + 1
    c <- c + 1
    species <- cluster_accession_list$species[ck]
    text(1,figure_height - (c-0.5) * offset_gene_height,species,cex=11,adj=0)
    for (i in 1:length(gene_matrix[species,])) {
      if (!is.na(gene_matrix[species,i])) {
        s <- shrink_factor * gene_start_matrix[species,i]
        e <- shrink_factor * gene_end_matrix[species,i]
        rect(s+species_name_offset,figure_height - c * offset_gene_height + 1,
             e+species_name_offset,
             figure_height-offset_gene_height * (c-1) - 30,
             col=gene_colors[gene_matrix[species,i]])
      }
    }
    # place sequence length at right
    sequence_length <- acclist[acclist$accession==acc,]$sequence_length
    text(25+e+species_name_offset,
         figure_height - (c-0.5) * offset_gene_height,sequence_length,cex=11,adj=0)
  }
}

# draw legend with no plot, just legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
par(mar=c(0,0,0,0))
legend('topleft',legend=names(gene_colors),fill=gene_colors,cex=15,ncol=n_legend_columns,bty='o')
dev.off()

########## END ##########
