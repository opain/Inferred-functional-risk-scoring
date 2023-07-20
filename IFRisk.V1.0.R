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
make_option("--clump_thresh", action="store", default=0.9, type='numeric',
		help="r2 threshold for clumping [optional]"),
make_option("--cor_window", action="store", default=5e6, type='numeric',
		help="Window for deriving pruning blocks in bases[optional]"),
make_option("--pTs", action="store", default='5e-1,1e-1,5e-2,1e-2,1e-3,1e-4,1e-5,1e-6', type='character',
		help="Window for deriving pruning blocks in bases[optional]"),
make_option("--clump_mhc", action="store", default=T, type='logical',
		help="Retain only the most significant gene within the MHC region [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)

if(file.exists(paste(opt$output,'-GeRS.csv',sep=''))){
	cat('A file named',paste(opt$output,'-GeRS.csv',sep=''),'already exists.\n')
	q()
}

system(paste0('mkdir -p ',opt$output_dir))

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

AllGene<-AllGene[c(1,2,var_all$column[which(var_all$variance != 0)])]
cat(length(var_all$column[var_all$variance == 0]),'features were removed due to zero variance.\n')

# Find intersecting features.
intersecting_genes<-intersect(TWAS$FILE, names(AllGene))

cat(length(intersecting_genes), 'features are present in both twas_results and target_gene_exp.\n')

AllGene_intersect<-AllGene[c('FID','IID',intersecting_genes)]
TWAS_intersect<-TWAS[(TWAS$FILE %in% intersecting_genes),]
TWAS_intersect<-TWAS_intersect[order(TWAS_intersect$CHR, TWAS_intersect$P0),]

rm(AllGene)
rm(TWAS)

if(opt$clump_mhc == T){
###########
# Clump the MHC region to only contain the most significant gene in that region
###########

TWAS_intersect_notMHC<-TWAS_intersect[!(TWAS_intersect$CHR == 6 & TWAS_intersect$P0 > 26e6 & TWAS_intersect$P1 < 34e6),]
TWAS_intersect_MHC<-TWAS_intersect[TWAS_intersect$CHR == 6 & TWAS_intersect$P0 > 26e6 & TWAS_intersect$P1 < 34e6,]
TWAS_intersect_MHC_retain<-TWAS_intersect_MHC[TWAS_intersect_MHC$TWAS.P == min(TWAS_intersect_MHC$TWAS.P),]
TWAS_intersect<-rbind(TWAS_intersect_notMHC,TWAS_intersect_MHC_retain)

AllGene_intersect<-AllGene_intersect[c('FID','IID',TWAS_intersect$FILE)]

cat(dim(TWAS_intersect)[1], 'after clumping the MHC region to contain only the top gene.\n')

}

###########
# Create a table showing the number of genes in each pT
###########
# The number of features after pruning will be filled in as each pT is processed.

NGenes_table<-NULL
for(i in 1:length(opt$pTs)){
	NGenes<-data.frame(	pT=opt$pTs[i],
						NGenes=sum(TWAS_intersect$TWAS.P < opt$pTs[i]),
						NGenes_post_clump=NA)
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

###
# Clump genes based on their correlation
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Clumping features...')
sink()

TWAS_clumped<-NULL
for(j in unique(TWAS_intersect$Block)){
	if(sum(TWAS_intersect$Block == j) == 1){
		TWAS_clumped<-rbind(TWAS_clumped,TWAS_intersect[TWAS_intersect$Block == j,])
	} else {
		TWAS_Block<-TWAS_intersect[TWAS_intersect$Block == j,]
		TWAS_Block<-TWAS_Block[order(TWAS_Block$TWAS.P),]
		cor_block<-cor(as.matrix(AllGene_intersect[,TWAS_Block$FILE]), method='pearson')

		i<-1
		while(i){
		  # Subset rows within range of variable i
		  TWAS_Block_i<-TWAS_Block[TWAS_Block$P0 < (TWAS_Block$P1[i] + opt$cor_window) & TWAS_Block$P1 > (TWAS_Block$P0[i] - opt$cor_window),]
		  cor_block_i<-cor_block[(dimnames(cor_block)[[1]] %in% TWAS_Block_i$FILE),(dimnames(cor_block)[[2]] %in% TWAS_Block_i$FILE)]
		  
		  # If no nearby features (after previous pruning) skip
		  if(dim(TWAS_Block_i)[1] == 1){
			i<-i+1
		  	next()
		  }
		  
		  # Skip if no variables in range have correlation greater than threshold
		  if(!(max(abs(cor_block_i[-1,1])) > sqrt(opt$clump_thresh))){
		    i<-i+1
				if(i > dim(TWAS_Block)[1]){
					break()
				}
		    next()
		  }
		  
		  # Create list of variables that correlate to highly with variable i
		  exclude<-dimnames(cor_block_i)[[1]][-1][abs(cor_block_i[-1,1]) > sqrt(opt$clump_thresh)]
		  
		  # Remove these correlated variables
		  TWAS_Block<-TWAS_Block[!(TWAS_Block$FILE %in% exclude),]
		  cor_block<-cor_block[!(dimnames(cor_block)[[1]] %in% exclude),!(dimnames(cor_block)[[2]] %in% exclude)]

			i<-i+1

			# Break if we have gone through all variables.
			if(i > dim(TWAS_Block)[1]){
				break()
			}
		}	
	TWAS_clumped<-rbind(TWAS_clumped,TWAS_Block)
	}
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Subset columns
TWAS_clumped<-TWAS_clumped[, names(TWAS_clumped) %in% c('FILE','PANEL','ID','TWAS.Z','TWAS.P')]

# Save list of clumped TWAS results with TWAS.Z and TWAS.P values
fwrite(TWAS_clumped, paste0(opt$output,'.score'), sep=' ')

# Update NGenes_table after clumping
for(i in 1:length(opt$pTs)){
	NGenes_table$NGenes_post_clump[i]<-sum(TWAS_clumped$TWAS.P <= opt$pTs[i])
}

write.table(NGenes_table, paste(opt$output,'.NFeat',sep=''), row.names=F, quote=F, sep='\t')

########
# Weight the gene expression in each individuals by TWAS.Z, scaling predicted expression in advance
########

IDs<-AllGene_intersect[,1:2]
TWAS_intersect<-TWAS_intersect[match(names(AllGene_intersect)[-1:-2], TWAS_intersect$FILE),]
AllGene_intersect_noID<-data.table(scale(AllGene_intersect[,-1:-2]))
AllGene_intersect_noID<-t(t(AllGene_intersect_noID) * TWAS_intersect$TWAS.Z)
AllGene_intersect<-cbind(IDs,AllGene_intersect_noID)
rm(AllGene_intersect_noID)

##################
# For each pT calculate GeRS.
##################

GeneX_Risk<-data.frame(	FID=AllGene_intersect$FID,
						IID=AllGene_intersect$IID)

for(i in 1:sum(NGenes_table$pT > min(TWAS_clumped$TWAS.P))){
	tmp<-rowSums(AllGene_intersect[which(names(AllGene_intersect) %in% TWAS_clumped$FILE[TWAS_clumped$TWAS.P <= NGenes_table$pT[i]])])
	GeneX_Risk<-cbind(GeneX_Risk,tmp)
	names(GeneX_Risk)[2+i]<-paste0('SCORE_',NGenes_table$pT[i])
}

########
# Save the gene expression risk scores
########

write.csv(GeneX_Risk, paste(opt$output,'-GeRS.csv',sep=''), row.names=F, quote=F)

end.time <- Sys.time()
time.taken <- end.time - start.time

cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,3)),attr(time.taken, 'units'),'\n')
sink()
