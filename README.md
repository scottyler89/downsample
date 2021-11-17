# downsample
downsample columns such that the sum is equal in all columns

# Installation

`devtools::install_github('scottyler89/downsample')`

# Examples
## How to downsample a count matrix

```r
library("downsample")
## For illustrative purposes, we'll generate a 
## negative binomial distributed count matrix
in_mat<-matrix(rnbinom(1000*20000,.1,c(0.5,0.5,0.5)),ncol=1000,nrow=20000)

```

In the context of scRNAseq, you can think of this as your count matrix with cells in columns and genes in rows.

Now you'll typically want to look at the distribution of colSums & filter out the garbage with low sums that frequently make it through.

```r
mat_col_sums <- colSums(in_mat)
hist(log10(mat_col_sums),breaks=20)
```

In real scRNAseq data this spread will be much bigger, will likely be multimodal, and will almost always have a long tail of low depth cells. In a more typical case, we'll filter out that long tail of low count cells & use a cutof for the minimum number of UMI that we actually want to include in the final cleaned dataset

```r
umi_cutoff <- 1900
keep_cols <- which(mat_col_sums>=umi_cutoff)
## in this case we're keeping ~95% of cells as noted below
print(length(keep_cols)/dim(in_mat)[2])
## now we'll actually do the downsampling
out_mat <- downsample_mat(in_mat[,keep_cols])
## Now you'll see that all 'cells' have the same total UMI
print(head(colSums(out_mat)))
```

If you're using downsampling as a way to normalize for UMI depth on different datasets to be merged. Then you'll want to use the same cutoff for each dataset!




