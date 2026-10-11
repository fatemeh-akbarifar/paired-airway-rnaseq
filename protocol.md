# Paired airway RNA-seq reanalysis

Data: the `airway` object in Bioconductor airway, derived from GSE52778. Eight samples, four cell lines, one untreated and one dexamethasone-treated sample per line. Counts are analysed without external normalisation. Retain genes with at least 10 counts in at least four samples. Primary DESeq2 design: `~ cell + dex`; contrast: treated minus untreated; alpha = 0.05. Keep independent filtering and Cook's-distance handling at their documented defaults. Report the number with missing p-values separately from those with missing adjusted p-values.

Select CRISPLD2 as the source-paper reference, the most significant other positively associated gene and the most significant negatively associated gene (adjusted p-value, then gene ID). A twofold flag is descriptive, not a formal minimum-effect test. Show VST PCA, an MA plot and paired normalised counts. Compare the selected estimates with an unpaired `~ dex` fit to the identical filtered matrix; never select the model by its discovery count. Four cell lines, not genes, are the biological replication units.

This is a secondary analysis of cell-culture data, not a clinical efficacy study. Current gene-symbol annotation may differ from the historical study. Preserve Ensembl keys and report annotation package versions. Raw reads and packaged counts are not republished.
