suppressMessages(suppressWarnings(library("DESeq2", character.only=T, warn.conflicts = F, quietly = T)))

# args from command line:
args<-commandArgs(TRUE)
RAW_COUNT_MATRIX<-args[1]
ANN_FILE<-args[2]
ORIG_ANN_COL<-args[3]
ORIG_CONDITION_A<-args[4]
ORIG_CONDITION_B<-args[5]
OUTPUT_DESEQ_FILE_BASE <- 'deseq2_results'
OUTPUT_NORMALIZED_COUNTS_BASE <- 'deseq2_normalized_counts'

# change the working directory to co-locate with the counts file:
working_dir <- dirname(RAW_COUNT_MATRIX)
setwd(working_dir)

# convert names in case they are not "proper" for R
CONDITION_A <- make.names(ORIG_CONDITION_A)
CONDITION_B <- make.names(ORIG_CONDITION_B)
ANN_COL <- make.names(ORIG_ANN_COL)

# create a string to identify the contrast:
contrast_str = paste0(ANN_COL, '_', CONDITION_B, '_vs_', CONDITION_A)

# read the annotation file:
annotations <- read.table(ANN_FILE, sep='\t', header = T, row.names = 1)

if (!(ANN_COL %in% colnames(annotations))) {
    message(sprintf('The column "%s" was not found in your annotation file. Please check your inputs. Note that this input is case-sensitive and must match exactly.', ORIG_ANN_COL))
    quit(status=1)
}

# filter to only those samples which match the conditions
base_samples_filter <- annotations[, ANN_COL] == CONDITION_A
exp_samples_filter <- annotations[, ANN_COL] == CONDITION_B
orig_base_samples <- rownames(annotations)[base_samples_filter]
orig_exp_samples <- rownames(annotations)[exp_samples_filter]
base_samples <- make.names(orig_base_samples)
exp_samples <- make.names(orig_exp_samples)

if(
    (length(base_samples) == 0)
    ||
    (length(exp_samples) == 0)
){
    message(sprintf('One or both of your sample sets contained zero samples. Check that both %s and %s are valid values in column "%s".', ORIG_CONDITION_A, ORIG_CONDITION_B, ORIG_ANN_COL))
    quit(status=1)
}

if(
    (length(base_samples) < 2)
    ||
    (length(exp_samples) < 2)
){
    message(sprintf('One or both of your sample sets contained fewer than 2 samples. To perform differential expression analysis, replicates are required. Check that both conditions %s and %s have 2 or more samples each.', ORIG_CONDITION_A, ORIG_CONDITION_B))
    quit(status=1)
}

intersection_list = intersect(base_samples, exp_samples)

if (length(intersection_list) > 0){
    sample_list = paste0(intersection_list, collapse=',')
    message(paste(
       'The following samples were in both contrast groups. Fix this and try again: ',
       sample_list
    ))
    quit(status=1)
}
all_samples <- c(base_samples, exp_samples)
original_sample_names <- c(orig_base_samples, orig_exp_samples)
colname_mapping = data.frame(
    orig_names = original_sample_names,
    row.names=all_samples,
    stringsAsFactors=F)

condition_a_list <- rep(CONDITION_A, length(base_samples))
condition_b_list <- rep(CONDITION_B, length(exp_samples))
all_conditions <- c(condition_a_list, condition_b_list)

# full annotation matrix:
annotations <- as.data.frame(cbind(all_samples, all_conditions), stringsAsFactors = F)
colnames(annotations) <- c('sample', 'condition')

# read the raw count matrix, genes as row names:
count_data <- read.table(RAW_COUNT_MATRIX, sep='\t', header = T, row.names = 1, stringsAsFactors = F)

# subset to keep only the samples in the count table.  This is important if the annotation
# file has more samples than the count matrix:
count_mtx_cols = colnames(count_data)
annotations <- annotations[annotations$sample %in% count_mtx_cols,]

# subset to only keep samples corresponding to the current groups in the count_data dataframe
count_data <- count_data[,annotations$sample]

# DESeq2 expects that the rownames of the annotation data frame are the sample names.  Set the rownames and drop that col
rownames(annotations) <- annotations$sample
annotations <- annotations[-1]

# Need to set the condition as a factor since it's going to be used as a design matrix
annotations$condition <- as.factor(annotations$condition)

fl <- length(levels(as.factor(annotations$condition)))
if(fl < 2){
    message(sprintf('After subsetting the matrix for the samples of interest (%d found), only one cohort of samples was present. Please double-check your inputs or sample names.',  dim(annotations)[1]))
    quit(status=1)
}

if (dim(count_data)[2] == 0){
    message('After subsetting the matrix for the samples of interest, the matrix was empty. Please check the input samples and matrix')
    quit(status=1)
}
# run the actual differential expression:
dds <- DESeqDataSetFromMatrix(countData = count_data,
							  colData = annotations,
							  design = ~condition)


# wraps the typical DESeq call to catch edge cases where the 
# table of expressions is small and we cannot use the typical
# methods to estimate the mean-dispersion relationship.
run_dge_func <- function(dds){
    tryCatch(
        {
            dds <- DESeq(dds)
            return (dds)
        },
        error=function(x){
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersionsGeneEst(dds)
            dispersions(dds) <- mcols(dds)$dispGeneEst
            dds <- nbinomWaldTest(dds)
            return (dds)
        }
    )
}

dds <- run_dge_func(dds)
res <- results(dds, contrast=c("condition", CONDITION_B, CONDITION_A))
original_colnames = colnames(res)
n = length(original_colnames)
baseMeanA = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == CONDITION_A]) 
baseMeanB = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == CONDITION_B]) 
res = cbind(rownames(res), res[,1],baseMeanA, baseMeanB, as.data.frame(res[,2:n])) 
colnames(res) = c('Gene', 'overall_mean', CONDITION_A, CONDITION_B, original_colnames[2:n])
resOrdered <- res[order(res$padj),]

# extract and output the normalized counts:
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)

# map the potentially 'mangled' names back to the original
nc_cols = colnames(nc)
remapped_cols = colname_mapping[nc_cols, 'orig_names']
colnames(nc) = remapped_cols
nc <- cbind(gene=rownames(nc), nc)
fout2 <- paste(OUTPUT_NORMALIZED_COUNTS_BASE, contrast_str, 'tsv', sep='.')
fout2 <- paste(working_dir, fout2, sep='/')
write.table(nc, fout2, sep='\t', quote=F, row.names=F)

# merge to create a single table, which makes frontend work easier
m <- merge(resOrdered, nc, by="row.names")
rownames(m) <- m[,'Row.names']
drops <- c("Row.names", "Gene", "gene")
m <- m[, !(names(m) %in% drops)]

# change column name for the 'stat' column which will match other dge-type analyses
cols <- colnames(m)
cols[which(names(m) == 'stat')] = 'statistic'
colnames(m) <- cols

output_filename <- paste(OUTPUT_DESEQ_FILE_BASE, contrast_str, 'tsv', sep='.')
output_filename <- paste(working_dir, output_filename, sep='/')
write.table(m, output_filename, sep='\t', quote=F)

json_str = paste0(
       '{"dge_results":"', output_filename, '",',
       '"normalized_counts":"', fout2, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)

