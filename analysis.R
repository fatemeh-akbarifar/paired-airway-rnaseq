suppressPackageStartupMessages({library(airway);library(DESeq2);library(ggplot2);library(org.Hs.eg.db);library(AnnotationDbi);library(jsonlite)})
options(stringsAsFactors=FALSE)
set.seed(42)
dir.create('results',showWarnings=FALSE)
data(airway)
cts <- assay(airway)
meta <- as.data.frame(colData(airway))
stopifnot(identical(colnames(cts),rownames(meta)), ncol(cts)==8, length(unique(meta$cell))==4,
          all(table(meta$cell,meta$dex)==1),all(is.finite(cts)),all(cts>=0),all(cts==round(cts)),!anyDuplicated(rownames(cts)))
meta$cell <- factor(meta$cell)
meta$dex <- relevel(factor(meta$dex),ref='untrt')
qc <- data.frame(sample=rownames(meta),cell=meta$cell,treatment=meta$dex,total_counts=colSums(cts),detected_genes=colSums(cts>0))
write.csv(qc,'results/sample-quality.csv',row.names=FALSE)
keep <- rowSums(cts>=10)>=4
filtered <- cts[keep,]
fit <- function(design) DESeq(DESeqDataSetFromMatrix(filtered,meta,design=design),quiet=TRUE,parallel=FALSE)
dds <- fit(~cell+dex)
res <- as.data.frame(results(dds,contrast=c('dex','trt','untrt'),alpha=.05))
res$gene_id <- rownames(res)
symbols <- mapIds(org.Hs.eg.db,keys=res$gene_id,column='SYMBOL',keytype='ENSEMBL',multiVals='first')
res$symbol <- unname(symbols[res$gene_id])
res$significant <- !is.na(res$padj)&res$padj<.05
res$twofold_flag <- abs(res$log2FoldChange)>=1
res <- res[,c('gene_id','symbol','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj','significant','twofold_flag')]
write.csv(res,'results/paired-results.csv',row.names=FALSE,na='NA')
# Save the transformed PCA coordinates so the figure can be audited.
v <- vst(dds,blind=FALSE)
pc <- plotPCA(v,intgroup=c('dex','cell'),returnData=TRUE)
pct <- round(100*attr(pc,'percentVar'),1)
write.csv(pc,'results/pca.csv',row.names=FALSE)
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),legend.position='bottom'))
colors <- c(untrt='#52616b',trt='#b44c74')
p <- ggplot(pc,aes(PC1,PC2,color=dex,label=cell))+geom_point(size=4)+geom_text(nudge_y=1.5,size=3,show.legend=FALSE)+scale_color_manual(values=colors)+labs(title='Eight samples, four paired cell lines',subtitle='Variance-stabilised expression; PCA of the 500 most variable genes',x=paste0('PC1 (',pct[1],'%)'),y=paste0('PC2 (',pct[2],'%)'),color='Dexamethasone')
ggsave('results/sample-variation.png',p,width=9,height=5,dpi=170,bg='white')
p <- ggplot(res,aes(baseMean,log2FoldChange,color=significant))+geom_point(size=.6,alpha=.5,na.rm=TRUE)+scale_x_log10()+scale_color_manual(values=c('FALSE'='#b5b9bf','TRUE'='#b44c74'))+geom_hline(yintercept=0,linewidth=.3)+labs(title='Treatment-associated expression changes',subtitle='Paired DESeq2 model; treated minus untreated; unshrunk estimates',x='Mean normalised count (log scale)',y='Estimated log2 fold change',color='Adjusted p < 0.05')
ggsave('results/expression-response.png',p,width=9,height=5,dpi=170,bg='white')
reference <- res$gene_id[!is.na(res$symbol)&res$symbol=='CRISPLD2']
stopifnot(length(reference)==1)
choose <- function(direction) {
 z <- res[res$significant & !is.na(res$log2FoldChange) & sign(res$log2FoldChange)==direction & !res$gene_id %in% reference,]
 if(!nrow(z)) return(character())
 z$gene_id[order(z$padj,z$gene_id)][1]
}
selected <- c(reference,choose(1),choose(-1))
norm <- counts(dds,normalized=TRUE)
gene_counts <- do.call(rbind,lapply(selected,function(id) data.frame(gene_id=id,symbol=res$symbol[match(id,res$gene_id)],sample=colnames(norm),cell=meta$cell,treatment=meta$dex,normalised_count=as.numeric(norm[id,]))))
write.csv(gene_counts,'results/selected-gene-counts.csv',row.names=FALSE)
p <- ggplot(gene_counts,aes(treatment,normalised_count,group=cell))+geom_line(color='#aeb5bd',linewidth=.7)+geom_point(aes(color=treatment),size=3)+facet_wrap(~symbol,scales='free_y')+scale_color_manual(values=colors)+labs(title='Keep the pairs visible',subtitle='One line per cell line; each panel has its own vertical scale',x=NULL,y='Normalised count',color='Dexamethasone')
ggsave('results/paired-genes.png',p,width=10,height=4.5,dpi=170,bg='white')
unpaired <- fit(~dex)
u <- as.data.frame(results(unpaired,contrast=c('dex','trt','untrt'),alpha=.05))
comparison <- res[match(selected,res$gene_id),c('gene_id','symbol','log2FoldChange','lfcSE','padj')]
for(k in c('log2FoldChange','lfcSE','padj')) comparison[[paste0('unpaired_',k)]] <- u[selected,k]
write.csv(comparison,'results/design-comparison.csv',row.names=FALSE,na='NA')
# Direction audit: compare fitted signs with paired normalised-count ratios.
# Match paired columns explicitly; never rely on an incidental treatment ordering.
ratios <- sapply(selected,function(id) mean(sapply(levels(meta$cell),function(cell){a=norm[id,meta$cell==cell & meta$dex=='trt'];b=norm[id,meta$cell==cell & meta$dex=='untrt'];log2((a+.5)/(b+.5))})))
stopifnot(all(sign(ratios)==sign(comparison$log2FoldChange)))
summary <- list(raw_genes=nrow(cts),retained_genes=sum(keep),available_pvalues=sum(!is.na(res$pvalue)),available_adjusted_pvalues=sum(!is.na(res$padj)),significant_genes=sum(res$significant),significant_and_twofold=sum(res$significant & res$twofold_flag,na.rm=TRUE),missing_pvalues=sum(is.na(res$pvalue)),missing_adjusted_pvalues=sum(is.na(res$padj)),cell_lines=4,samples=8,paired_sign_check=TRUE,selected_genes=as.list(comparison$symbol),filter_threshold=metadata(results(dds,contrast=c('dex','trt','untrt'),alpha=.05))$filterThreshold)
write_json(summary,'results/summary.json',pretty=TRUE,auto_unbox=TRUE)
packages <- c('airway','DESeq2','SummarizedExperiment','ggplot2','org.Hs.eg.db','AnnotationDbi','jsonlite')
write.csv(data.frame(package=packages,version=sapply(packages,function(x)as.character(packageVersion(x)))),'results/package-versions.csv',row.names=FALSE)
# Session information is useful, but installation paths identify a local machine.
session <- capture.output(sessionInfo())
session <- gsub('/(Users|home)/[^ ]+','<local-installation>',session)
writeLines(session,'results/session-info.txt')
writeLines(c('Dataset: airway object, GSE52778','Source: https://bioconductor.org/packages/airway/','Package licence: LGPL; upstream counts obtained through package, not redistributed','Model: ~ cell + dex; contrast trt - untrt',paste('analysis.R MD5:',unname(tools::md5sum('analysis.R')))),'results/provenance.txt')
print(summary);print(comparison)
