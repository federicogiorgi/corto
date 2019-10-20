# corto
corto (Correlation Tool): an R package to generate correlation-based DPI networks
![alt text](https://giorgilaborg.files.wordpress.com/2019/10/cortoicon.png)

# Introduction
The _corto_ ("correlation tool") package provides a pipeline to infer networks between "centroid" and "target" variables in a dataset, using a combination of Spearman correlation and Data Processing Inequality (DPI), first proposed in [1]. The main application of _corto_ is in the field of Bioinformatics and Transcriptomics, where co-occurrence between variables can be used as a mean to infer regulatory mechanisms [2] or gene functions [3]. In this field, usually the tested features are genes (or rather, their expression profile across samples), whereas the centroids are Transcription Factors (TFs) and their targets are Target Genes (TGs). The TF-TG co-expression can hint at a causal regulatory relationship, as proven in many studies [4,5,6]. The _corto_ tool replicates the well-established pipeline of the ARACNe family of tools [7,8,9]

In brief, _corto_ operates using the following steps:

1. Calculate all pairwise Spearman correlation coefficients between centroid and target features. The rank transformation operated by the Spearman correlation coefficient is a common procedure used in co-expression based studies, due to its benefits in reducing the effects of outliers [10].
2. Filter out all centroid-target features whose correlation coefficient _absolute value_ is below the provided p-value threshold _p_.
3. Apply DPI to all centroid-centroid-target triplets, in order to identify which centroid-target correlation is stronger and identify the most likely association (e.g. which TF is the most likely regulator of the TG in the dataset).
4. All edges are tested for Data Processing Inequality in a number of bootstrapped versions of the same input matrix (specified by the _nbootstraps_ parameter, 100 by default, as in ARACNE-AP [8]). This step can be run in parallel by specifiying the number of threads using the _nthreads_ parameter. The number of occurrences of each edge (if they survived DPI) is calculated. This number can range from 0 (the edge is not significant in any bootstrap) to _nbootstraps_+1 (the edge is significant in all bootstraps, and also in the original matrix).
5. A network is generated in the form of a list, where each element is a centroid-based list of targets. In order to follow the structure of the _regulon_ class implemented by downstream analysis packages (e.g. VIPER [11]), each target is characterized by two parameters:
    + The _tfmode_, the mode of action, i.e. the Spearman correlation coefficient between the target and the centroid in the original input matrix
    + The _likelihood_ of the interaction, calculated as the number of bootstraps in which the edge survives DPI, divided by the total number of bootstraps (in this impelmentation, the presence of an edge in the original, non-bootstrapped matrix is considered as an additional evidence)


# Running _corto_
Here is how to run _corto_. First, install the package:
```{r install, eval=FALSE}
install.packages("corto")
```

Then, load it:
```{r load}
library(corto)
```

Then, you can see how the input matrix looks like. For example, this dataset comes from the TCGA mesothelioma project [12] and measures the expression of 10021 genes across 87 samples:
```{r load1}
load(system.file("extdata","inmat.rda",package="corto"))
inmat[1:5,1:5]
```
```{r load2}
dim(inmat)
```
Another input needed by _corto_ is a list of centroid features. In our case, we can specify a list of TFs generated from Gene Ontology with the term "Regulation of Transcription" [13].

```{r load3}
load(system.file("extdata","centroids.rda",package="corto"))
centroids[15]
```
```{r load4}
length(centroids)
```

Finally, we can run _corto_. In this example, we will run it with p-value threshold of 1e-30, 10 bootstraps and 2 threads
```{r runcorto,message=FALSE,results="hide"}
regulon<-corto(inmat,centroids=centroids,nbootstraps=10,p=1e-30,nthreads=2)
# Input Matrix has 87 samples and 10021 features
# Correlation Coefficient Threshold is: 0.889962633618839
# Removed 112 features with zero variance
# Calculating pairwise correlations
# Initial testing of triplets for DPI
# 246 edges passed the initial threshold
# Building DPI network from 37 centroids and 136 targets
# Running 100 bootstraps with 2 thread(s)
# Calculating edge likelihood
# Generating regulon object
```

The regulon object is a list:
```{r prinregulon}
regulon[1:2]
```

The regulon in this dataset is composed of 34 final centroids with at least one target:
```{r prinregulon2}
length(regulon)
```
