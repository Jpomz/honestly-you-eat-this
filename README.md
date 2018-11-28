# Inferring predator-prey interactions in food webs

The code included in this project reproduces the analysis presented in: Pomeranz et al. _Inferring predator-prey interactions in food webs. Methods in Ecology and Evolution_

The full dataset is available on DataDryad [doi:10.5061/dryad.k59m37f]. 

To run the full analysis, download the data and place the files in the Data/Raw_data folder. 

The scripts should be run in order, e.g. 1_trait_matching_models.R, 1b_..., 2_..., 3_... etc. 

If you have any questions, please feel free to contact me. 
jfpomeranz@gmail.com

Below is an example of the novel methods used in this analysis. 

For a full description of the WebBulder method, please see Gray et al. 2015, Joining the dots: An automated method for constructing fod webs from compendia of published interactions. Food Webs. http://dx.doi.org/10.1016/j.fooweb.2015.09.001 
The supplementary material of Gray et al. 2015 contains the code for the `WebBuilder` functions, as well as an example. 

For a full description of the initial trait-matching model, please see Gravel et al 2013, Inferring food web structure from predator-prey body size relationships. Methods in Ecology and Evolution. doi: 10.1111/2041-210X.12103. 
The supplementary information for Gravel et al. 2013 contains functions and an example for inferring interactions using predator-prey body sizes. 

## Trait-matching

Read in the example data.

```{r}
# load example data
example.data <- readRDS("Data/example_data.RDS")
```

`example.data` is a list containing data for the Demspters Creek site. 

```{r}
summary(example.data)
```

`example.data$observed.A` is the empirical adjacency matrix. Note that this matrix has been modified from the original, as described in the main text. 

`example.data$dw.ab` is a data.frame with four columns and 33 rows. 
`taxa` is the taxonomic name
`dw` is the estimated dry weight in grams
`no.m2` is the estimated number of individuals per meter^2^. Note that the abundances of the fish taxa have _not_ already been "corrected" as described in main text. 

`example.data$tm.initial` is an adjacency matrix of inferred interactions using the method described in Gravel et al. 2013, Methods in Ecology and Evolution doi: 10.1111/2041-210X.12103. 
`example.data$niche.forbidden` is a character vector of taxa names which are considered to be niche-forbidden (see main text for details) 

Read in scripts which contain useful functions for this example.

```{r}
# load function scripts
# useful food web functions, modified from Petchey
source("Functions/FoodWeb_Functions.R")
# functions written for this manuscript
source("Functions/Inference_MS_functions.R")
```

Plot the empirical and initial trait-matching inference adjacency matrices, and caluclate the TSS and AUC:
```{r}
Plot.matrix2(example.data$observed.A,
             sp.pt.ch = 18, point.cex = 1)
title(main = "Observed")
Plot.matrix2(example.data$tm.initial,
             sp.pt.ch = 18, point.cex = 1)
title(main = "TM initial")

list(auc = get_auc(observed = example.data$observed.A,
                   inferred = example.data$tm.initial),
     tss = get_tss(observed = example.data$observed.A,
                   inferred = example.data$tm.initial))

```

As we can see from the plot, the initial inference greatly over predicts the number of links. Likewise, the AUC is 0.59 (e.g. only slightly better than a coin toss), and the TSS is low. In order to improve our inference, we remove "niche forbidden" links using the function `rm_niche()` found in `Inference_MS_functions.R`. This function requires a vector of taxa names that we want to restrict.

```{r}
# remove niche forbidden taxa
tm.niche <- rm_niche(inf = example.data$tm.initial,
                     taxa = example.data$niche.forbidden)
# plot new inference
Plot.matrix2(tm.niche,
             sp.pt.ch = 18, point.cex = 1)
title(main = "TM Niche")

list(auc = get_auc(observed = example.data$observed.A,
                   inferred = tm.niche),
     tss = get_tss(observed = example.data$observed.A,
                   inferred = tm.niche))
```

Now our AUC and TSS are both much higher, indicating a better predicition. 

Next, we want to restrict "neutrally forbidden" links based on taxa densities. However, first we need to create a matrix N, where each element Nij is the product of species' _i_ and _j_ relative abundances. 

The `get_rel_ab()` function in the `Inference_MS_functions.R` script calculates relative abundances for each taxa, and then creates a matrix of their crossproducts. The arguments for the function are `vec`, which is a vector of numerical abundances, and `taxa`, which is a vector of taxa names.  
__Make sure that the order of the two vectors matches the order of the adjacency matrices!!!__  

```{r}
# make sure order of taxa names in dw_ab object match order in adjaceny matrix
identical(rownames(example.data$observed.A),
          example.data$dw.ab$taxa)

# calculate relative abundance matrix, N
N <- get_rel_ab(vec = example.data$dw.ab$no.m2,
                taxa = example.data$dw.ab$taxa)

# make sure rownames of N match rownames of observed
identical(rownames(N),
          rownames(example.data$observed.A))
```

Now we need to "correct" fish abundances as described in the main text. We use the `f_ab_corr()` function in the `Inference_MS_functions.R`

```{r}
# only multiply column names that match the "taxa" vector
N.fish <- f_ab_corr(Nij = N,
                    taxa = c("Salmo", "Galaxias",
                             "Anguilla", "Gobiomorphus"),
                    cf = 1000)
```

Next we convert values in N.fish to a binary matrix by setting values < n' to 0, and values > n' to 1 using the `rm_neutral()` function in the `Inference_MS_functions.R` script. n' is an arbitrarily defined threshold, below which species are considered to be too rare to be likely to interact (see main text). Here we use the neutral abundance threshold of 3e-4 as in the main text (see analysis script 1 and 1b for details on how we determined this threshold). 

__Make sure to multiply the right relative abundance matrix!__

```{r}
N.binary <- rm_neutral(Nij = N.fish,
                       threshold = 3.0e-04)
```

Now we multiply `N.binary` by `tm.niche` to remove neutrally forbidden links. e.g. only retaining links which appear in both. First, make sure that both matrices are the same dimensions, and have the same rownames

```{r}
identical(dim(tm.niche), dim(N.binary))
identical(rownames(tm.niche), rownames(N.binary))
# multiply two matrices together
tm.niche.neutral <- tm.niche * N.binary
# plot new inference matrix
Plot.matrix2(tm.niche.neutral,
             sp.pt.ch = 18, point.cex = 1)
title(main = "TM niche + neutral")
list(auc = get_auc(observed = example.data$observed.A,
                   inferred = tm.niche.neutral),
     tss = get_tss(observed = example.data$observed.A,
                   inferred = tm.niche.neutral))
```

In this case the AUC and TSS decrease slightly when restricting links at this neutral abundance threshold. The threshold was selected by averaging the AUC and TSS across all sites. Therefore, if you were trying to maximize predictive power at this site, a slightly different threshold or fish correction factor may need to be used. 
