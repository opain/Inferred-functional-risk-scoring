#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at Cardiff University under the supervision of Richard Anney.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--twas_results", action="store", default=NA, type='character',
		help="File containing TWAS results [required]"),
make_option("--target_gene_exp", action="store", default=NA, type='character',
		help="File containing gene expression values for the target sample [required]"),
make_option("--output", action="store", default=NA, type='character',
		help="Name of output files [required]"),
make_option("--prune_thresh", action="store", default=0.9, type='numeric',
		help="r2 threshold for pruning [optional]"),
make_option("--cor_window", action="store", default=5e6, type='numeric',
		help="Window for deriving pruning blocks in bases[optional]"),
make_option("--pTs", action="store", default='5e-1,1e-1,5e-2,1e-2,1e-3,1e-4,1e-5,1e-6', type='character',
		help="Window for deriving pruning blocks in bases[optional]"),
make_option("--prune_mhc", action="store", default=T, type='logical',
		help="Retain only the most significant gene within the MHC region [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$pTs<- as.numeric(unlist(strsplit(opt$pTs,',')))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# Inferred Risk Scoring (IRIS)
# V1.0 31/07/2018
#################################################################

Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')

if(file.exists(paste(opt$output,'-GeRS.csv',sep=''))){
	cat('A file named',paste(opt$output,'-GeRS.csv',sep=''),'already exists.\n')
	q()
}

sink()
sink("/dev/null")
suppressMessages(library(data.table))

###########
# Read in the input files and find intersect
###########

# Read in the TWAS results
if(substr(opt$twas_results, nchar(opt$twas_results)-2, nchar(opt$twas_results)) == '.gz'){
	TWAS<-data.frame(fread(paste('zcat ',opt$twas_results,sep='')))
} else {
	TWAS<-data.frame(fread(opt$twas_results))
}

sink()
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('twas_results contains', dim(TWAS)[1],'rows.\n')

# Convert TWAS FILE column to match the gene expression column names in expression_ref
TWAS$FILE<-sub(".*/", "", TWAS$FILE)
TWAS$FILE<-sub(".wgt.RDat", "", TWAS$FILE)
TWAS$FILE<-gsub(":", ".", TWAS$FILE)

# Remove rows with containing NAs or duplicate FILE ID.
TWAS<-TWAS[!is.na(TWAS$TWAS.Z),]
TWAS<-TWAS[!duplicated(TWAS$FILE),]

cat('twas_results contains', dim(TWAS)[1],'unique features with non-missing TWAS.Z values.\n')
cat('The minimum TWAS.P value is ', min(TWAS$TWAS.P),'.\n',sep='')

# Add the 0.5 Mb window to P0 and P1 values t account for window size when deriving weights.
TWAS$P0<-TWAS$P0-5e5
TWAS$P0[TWAS$P0 < 0]<-0
TWAS$P1<-TWAS$P1+5e5

# Read in the gene expression values
sink()
sink("/dev/null")

if(substr(opt$target_gene_exp, nchar(opt$target_gene_exp)-2, nchar(opt$target_gene_exp)) == '.gz'){
	AllGene<-data.frame(fread(paste('zcat ',opt$target_gene_exp,sep='')))
} else {
	AllGene<-data.frame(fread(opt$target_gene_exp))
}

sink()
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('target_gene_exp contains', dim(AllGene)[1],'individuals and', dim(AllGene)[2],'features.\n')

# Remove features with zero variance.
var_all<-NULL
for(i in 3:dim(AllGene)[2]){
	var<-data.frame(column=i,
					variance=var(AllGene[,i]))
	var_all<-rbind(var_all,var)
}

var_all$column[var_all$variance == 0]
AllGene<-AllGene[-var_all$column[var_all$variance == 0]]
cat(length(var_all$column[var_all$variance == 0]),'features were removed due to zero variance.\n')

# Find intersecting features.
intersecting_genes<-intersect(TWAS$FILE, names(AllGene))

cat(length(intersecting_genes), 'features are present in both twas_results and target_gene_exp.\n')

AllGene_intersect<-AllGene[c('FID','IID',intersecting_genes)]
TWAS_intersect<-TWAS[(TWAS$FILE %in% intersecting_genes),]

rm(AllGene)
rm(TWAS)

if(opt$prune_mhc == T){
###########
# Prune the MHC region to only contain the most significant gene in that region
###########

TWAS_intersect_notMHC<-TWAS_intersect[!(TWAS_intersect$CHR == 6 & TWAS_intersect$P0 > 26e6 & TWAS_intersect$P1 < 34e6),]
TWAS_intersect_MHC<-TWAS_intersect[TWAS_intersect$CHR == 6 & TWAS_intersect$P0 > 26e6 & TWAS_intersect$P1 < 34e6,]
TWAS_intersect_MHC_retain<-TWAS_intersect_MHC[TWAS_intersect_MHC$TWAS.P == min(TWAS_intersect_MHC$TWAS.P),]
TWAS_intersect<-rbind(TWAS_intersect_notMHC,TWAS_intersect_MHC_retain)

AllGene_intersect<-AllGene_intersect[c('FID','IID',TWAS_intersect$FILE)]

cat(dim(TWAS_intersect)[1], 'after pruning the MHC region to contain only the top gene.\n')

}

###########
# Create a table showing the number of genes in each pT
###########
# The number of features after pruning will be filled in as each pT is processed.

NGenes_table<-NULL
for(i in 1:length(opt$pTs)){
	NGenes<-data.frame(	pT=opt$pTs[i],
						NGenes=sum(TWAS_intersect$TWAS.P < opt$pTs[i]),
						NGenes_post_prune=NA)
	NGenes_table<-rbind(NGenes_table,NGenes)
}

#########
# Determine gene blocks.
#########
TWAS_intersect$Block<-NA
for(i in 1:dim(TWAS_intersect)[1]){
if(i == 1){
		TWAS_intersect$Block[i]<-1
		} else {
		if(i > 1 & TWAS_intersect$CHR[i] == TWAS_intersect$CHR[i-1] & TWAS_intersect$P1[i] > (TWAS_intersect$P0[i-1] - opt$cor_window) & TWAS_intersect$P0[i] < (TWAS_intersect$P1[i-1] + opt$cor_window)){
				TWAS_intersect$Block[i]<-TWAS_intersect$Block[i-1]
				}
		if(!(i > 1 & TWAS_intersect$CHR[i] == TWAS_intersect$CHR[i-1] & TWAS_intersect$P1[i] > (TWAS_intersect$P0[i-1] - opt$cor_window) & TWAS_intersect$P0[i] < (TWAS_intersect$P1[i-1] + opt$cor_window))){
				TWAS_intersect$Block[i]<-TWAS_intersect$Block[i-1]+1
				}
		}
}

cat('The features were be separated into',length(unique(TWAS_intersect$Block)),'blocks.\n')

########
# Weight the gene expression in each individuals by TWAS.Z
########

AllGene_intersect_weighted<-AllGene_intersect

for(j in 3:dim(AllGene_intersect)[2]){
AllGene_intersect_weighted[,j]<-scale(AllGene_intersect[,j])*TWAS_intersect$TWAS.Z[TWAS_intersect$FILE == names(AllGene_intersect[j])]
}

##################
# For each pT, prune genes, and then calculate GeRS.
##################

GeneX_Risk<-data.frame(	FID=AllGene_intersect$FID,
						IID=AllGene_intersect$IID)

for(i in 1:sum(NGenes_table$pT > min(TWAS_intersect$TWAS.P))){
TWAS_intersect_pT<-TWAS_intersect[TWAS_intersect$TWAS.P < opt$pT[i],]

########
# Prune genes
########
# Calculate correlation matrix for each block and remove genes that have any correlation value greater than opt$prune_thresh.
TWAS_intersect_pT_pruned_all<-NULL
for(j in unique(TWAS_intersect_pT$Block)){
	if(sum(TWAS_intersect_pT$Block == j) == 1){
		TWAS_intersect_pT_pruned<-TWAS_intersect_pT[TWAS_intersect_pT$Block == j,]
	} else {
		cor_block<-WGCNA::cor(as.matrix(AllGene_intersect[(names(AllGene_intersect) %in% TWAS_intersect_pT$FILE[TWAS_intersect_pT$Block == j])]), method='pearson')
		cor_block_pruned<-cor_block
		j<-1
		while(TRUE){
			tmp<-cor_block_pruned
			tmp[!lower.tri(tmp)] <- 0
			if(max(tmp) < sqrt(opt$prune_thresh)){
				break()
			}
			
			if(sum(abs(tmp[,j]) > sqrt(opt$prune_thresh)) > 0){
					cor_block_pruned<-cor_block_pruned[-j,-j]
					j<-1
			} else {
			j<-j+1
			}
		}
		TWAS_intersect_pT_pruned<-TWAS_intersect_pT[(TWAS_intersect_pT$FILE %in% colnames(cor_block_pruned)),]
	}
	
	TWAS_intersect_pT_pruned_all<-rbind(TWAS_intersect_pT_pruned_all, TWAS_intersect_pT_pruned)
}

n_prune<-sum(!(TWAS_intersect_pT$FILE %in% TWAS_intersect_pT_pruned_all$FILE))

AllGene_intersect_pT_pruned_weighted<-AllGene_intersect_weighted[c('FID','IID',names(AllGene_intersect)[(names(AllGene_intersect) %in% TWAS_intersect_pT_pruned_all$FILE)])]

NGenes_table$NGenes_post_prune[i]<-dim(TWAS_intersect_pT_pruned_all)[1]

# Sum effect size weighted gene expression values. 1e6 only has one gene in it so don't use row sums.
	SCORE<-data.frame(SCORE=rowSums(AllGene_intersect_pT_pruned_weighted[-1:-2]))
	names(SCORE)<-paste('SCORE_pT_',as.character(opt$pTs[i]),sep='')
	GeneX_Risk<-cbind(GeneX_Risk,SCORE)

cat(i,'of',length(opt$pTs),'pTs complete\n')
rm(AllGene_intersect_pT_pruned_weighted)
}

########
# Save the gene expression risk scores
########

write.csv(GeneX_Risk, paste(opt$output,'-GeRS.csv',sep=''), row.names=F, quote=F)

########
# Save the gene expression risk scores
########

write.csv(NGenes_table, paste(opt$output,'-NGene_Table.csv',sep=''), row.names=F, quote=F)

end.time <- Sys.time()
time.taken <- end.time - start.time

cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,3)),attr(time.taken, 'units'),'\n')
sink()
