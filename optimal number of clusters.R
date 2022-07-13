

##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################


#' Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/beta-diversity/)
#' Note: the path is denoted by forward slash "/"
setwd("C:/....../....../.....")  #<--- CHANGE ACCORDINGLY !!!

#' Please give the name of the meta-file that contains individual sample information
input_meta = "4.mapping_file.tab"

#' Please give the name of the OTUs table
file_name = "1.OTUs-Table.tab"      

#' Please insert "YES" if the OTUs table is normalized otherwise type "NO"
normalized = "NO"

#' Please give the name of the input tree
input_tree = "2.OTUs-NJTree.nwk"      

#' Please provide the name of the column (in the mapping file) based on which the samples are partitioned into individual groups
mapping_column = "Category"

#' Please place in the vector one or more names which will be used to identify the samples that composing the reference group  (e.g reference_name <- c("a","b"))
reference_name = c("OB")

#' De-Novo Clustering will be performed for the number of samples or maximal for the set limit
kmers_limit = 10



######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######

##################################################################################
######                             Main Script                              ######
##################################################################################

###################       Load all required libraries     ########################



# Check if required packages are already installed, and install if missing
packages <-c("GUniFrac","phangorn","factoextra","cluster","fpc","dplyr","graphics","ade4","vegan","stats") 

# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos ="http://cloud.r-project.org/")
  } 
}

# Applying the installation on the list of packages
lapply(packages, InsPack)

# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)

# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))


###################################################################################
#################     Functions to be used  in main Script.       #################
###################################################################################

#------------------------ Function that generates MDS plots-----------------------#

# distance -> the given distance matrix 
# groups -> a vector that contains the indexes of the groups 

sclass <- function(dist,all_groups) {

  
  mds<- cmdscale(dist,eig=T, x.ret=T)
  mds.variation<- round(mds$eig/sum(mds$eig)*100,1) # axes variance calculation
  plot_color<-rainbow(length(levels(all_groups)))[all_groups]
  all_groups_comp <- all_groups[!is.na(all_groups)]
  all_groups_comp<-factor(all_groups_comp,levels(all_groups_comp)[unique(all_groups_comp)])
  s.class(
    mds$points, col = unique(plot_color), cpoint =
      2, fac = all_groups_comp )
  graphics:: title (main=paste0("MDS graph k= ",k ), xlab = paste("MDS1:",mds.variation[1],"%"), ylab=paste("MDS2:",mds.variation[2],"%"))

}

##################  Create all the necessary paths and directories ###############

# Save the current path in one variable (will be used later)
OriginalPath <- getwd()

# Create the directory where all output files will be saved 
results_folder <- "Optimal number of Clusters"
results_path <- paste0(OriginalPath,"/",results_folder)
dir.create(results_folder)


#################################################################
##############          Data Loading          ###################
#################################################################


#------------------- Meta_file--------------------------#
# Load the mapping file containing individual sample information (sample names in the first column)
meta_file <- read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")

# Clean table from empty lines
meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),])

# Order meta file
meta_file <- meta_file[order(meta_file[,mapping_column]),] 

# Convert the selected column (mapping_column) of meta_file to factor
meta_file[,mapping_column] <- as.factor(meta_file[,mapping_column])
#-------------------------------------------------------#


#------------------------ OTUs table--------------------------#
# Load the tab-delimited file containing the values to be be checked (rownames in the first column)
otu_table <-  read.table (file_name,check.names = FALSE,header = TRUE,dec = ".",sep = "\t", row.names = 1,comment.char = "")

# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]

# Delete if exists the column with taxonomy information
taxonomy <- otu_table %>% select_if(is.factor)
if (ncol(taxonomy)==0) {taxonomy <-  otu_table %>% select_if(is.character)}
if (ncol(taxonomy)!=0) {
  otu_table[,colnames(taxonomy)] <- NULL
}

#--Normalizethe OTUs table (if the user ask for)--#
if (normalized=="YES") {
  
  # Check if Sample are in rows or columns
  if (any(colnames(otu_table) %in% rownames(meta_file))==TRUE) {
    
    # keep only those rows that appear in the mapping file
    otu_table <- otu_table[,rownames(meta_file)]
    
    # Transpose OTU-table and convert format to a data frame
    otu_table<- data.frame(t(otu_table))
    
    # Create the otu_file
    otu_file <- otu_table[rownames(meta_file),] } else {
      
      # keep only those rows that appear in the mapping file
      otu_table <- otu_table[rownames(meta_file),]
      
      # Convert format to a data frame
      otu_table<- data.frame(otu_table)
      
      # Create the otu_file
      otu_file <- otu_table[rownames(meta_file),] }} else {
        
        # Calculate the minimum sum of all columns/samples
        min_sum <- min(colSums(otu_table))
        
        # Divide each value by the sum of the sample and multiply by the minimal sample sum
        otu_file <- t(min_sum * t(otu_table) / colSums(otu_table))
        
        # Clean table from empty lines
        otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
        
        # Transpose OTU-table and convert format to a data frame
        otu_file <- data.frame(t(otu_file))
        
        # Create the otu_file
        otu_file <- otu_file[rownames(meta_file),]
        
      } 
#-------------------------------------------------------#


#---------------- Tree -----------------------#
# Load the phylogenetic tree calculated from the OTU sequences 
tree_file <- read.tree(input_tree)
tree_file <- keep.tip(tree_file,colnames(otu_file))


# Root the OTU tree at midpoint 
rooted_tree <- midpoint(tree_file)
#--------------------------------------------#

##############################################################################
############    Unifrac Distance Matrix and PAM clustering        ############
##############################################################################


# OTUs table and mapping file transformation in order to include only the chosen by the user groups
otu_file <- rbind(otu_file[meta_file[,mapping_column] %in% reference_name,])
meta_file <- rbind(meta_file[meta_file[,mapping_column] %in% reference_name,])


# Calculate the UniFrac distance matrix for comparing microbial communities
unifracs <- GUniFrac(otu_file, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
# Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
unifract_dist <- unifracs[, , "d_0.5"]



#############################################################################
##################### Calculation of the Indices  ###########################
#############################################################################

setwd(results_path)

ch_nclusters=NULL
sil_nclusters=NULL

# Create a PDF file where the MDS plot for the different k will be stored
pdf("Reference clusters.pdf")

if (dim(otu_file)[1]-1 <= kmers_limit) {
  kmers_limit=dim(otu_file)[1]-1
}
for (k in 1:kmers_limit) { 
  if (k==1) {
    ch_nclusters[k]=NA 
    sil_nclusters[k]=NA
    data_cluster=as.vector(pam(as.dist(unifract_dist), k, diss=TRUE)$clustering)
    
    #MDS plot for k=1
    sclass(unifract_dist,as.factor(data_cluster))
  } else {
    # Partitioning the data into k clusters (max k is number of samples within the dataset)
    data_cluster=as.vector(pam(as.dist(unifract_dist), k, diss=TRUE)$clustering)
    
    # Calculate Calinski-Harabasz and silhouette Index 
    index=cluster.stats(as.dist(unifract_dist),data_cluster)
    ch_nclusters[k] <- index[["ch"]]
    sil_nclusters[k] <-index[["avg.silwidth"]]
    print(k)
 
    #  MDS plot for k 2-kmers_limit
    sclass(unifract_dist,as.factor(data_cluster))
    
  
  }
}
dev.off()

# Calculate Within sum of squares
wss <- fviz_nbclust(otu_file ,diss=unifract_dist, k.max = kmers_limit  ,cluster::pam, method =  "wss")
# Calculate gap Statistics index
gap <- fviz_nbclust(otu_file ,diss=unifract_dist, k.max = kmers_limit  ,cluster::pam, method =  "gap_stat")
# Calculate prediction strength index
prediction_strength <- prediction.strength(as.dist(unifract_dist), Gmin = 2, Gmax = 10, clustermethod = claraCBI)


#################################################################################
######                        Write Output Files                           ######
#################################################################################

# Generated plots showing the optimal number of clusters
pdf(paste0("optimal number of clusters-",paste(reference_name,collapse="+"),".pdf"))

# Plot Calinski-Harabasz plot
plot(ch_nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters (CH)")
# Plot silhouette Index plot
plot(sil_nclusters, type="h", xlab="k clusters", ylab="Average silhouette width",main="Optimal number of clusters(Silhouette)")
# Plot WSS plot
plot(1:kmers_limit,wss[["data"]][["y"]], type="b", xlab="k clusters", ylab="Within sum of squares",main="Optimal number of clusters(WSS)",pch=19)
# Plot prediction strength plot
plot(2:10, prediction_strength$mean.pred[2:10], ylab="Prediction Strength",ylim=c(0,1.1), type="b", pch=1, xlab="Number of Clusters")
abline(.9,0, lty=5, col="grey70")
abline(0.8,0,lty=8, col="grey70")
# Plot gap Statistics plot
plot(gap)
dev.off()

setwd(OriginalPath)

# Graphical output files are generated in the main part of the script
if(!flag) { stop("
    It was not possible to install all required R libraries properly.
                 Please check the installation of all required libraries manually.\n
                 Required libaries: GUniFrac, phangorn, factoextra, cluster, fpc,dplyr")
}

#################################################################################
######                           End of Script                             ######
#################################################################################
