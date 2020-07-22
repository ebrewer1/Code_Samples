# Code_Samples
This is a diverse sampling of my coding experience in bioinformatics, and I will give a description of each project below.

The "Coronavirus_Alignment.sh" and "Coronavirus_Taxonomy_PCA.R" files were for a project in my Genomics class, where I aligned SARS-Cov-2, SARS, MERS, and 3 common cold-causing coronaviruses in a multiple alignment, as well as to SARS-Cov-2 individually. In the R script, I made a taxonomy based on the multiple alignment file, and then performed a principle components analysis (PCA) on the samples to see how related the viruses were to eachother on a 2D graph.

The "Nebulin_KO_DEG.Rmd" file was a project for my Statistical Genomics class, in which I took microarray data from this repository: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70213, on the the Gene Expression Omnibus (GEO). I used the limma, qvalue, and statmod packages among others to perform fitting and find significance between various sample contrasts. I then found differentially expressed genes (DEGs) between the contrasts and displayed them using vennCounts. Next, I performed a Gene Ontology using goSTAG for 3 categories: Molecular Function, Biological Process, and Cellular Components. Finally, I performed a PCA to see how similar the samples were to eachother on a 2D graph.

The "Drosophila_RNA_seq_analysis.txt" file was for my Applied Bioinformatics course, in which I downloaded RNA-seq data from this project https://www.ncbi.nlm.nih.gov/bioproject/PRJNA604807. I aligned the sample sequences to the Drosophila genome, used featureCounts to find the transcript counts that overlap annotations, and used DESeq to find the differentially expressed transcripts. 

