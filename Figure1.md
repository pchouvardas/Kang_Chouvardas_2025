# My Project

## Loading libraries

<details>
<summary>Code</summary>

``` r
library(ggplot2)
library(pals)
library(factoextra)
library(reshape)
library(ggpubr)
library(matrixStats)
library(tidyverse)
library(decoupleR)
library(rcartocolor)
`%!in%` = Negate(`%in%`)
library(DESeq2)
library(ggrepel)
library(rcartocolor)
library(colourvalues)
library(readxl)
library(ggh4x)
```

</details>

## Reading master file

<details>
<summary>Code</summary>

``` r
cores <- read_xlsx("masterfile_GS_TC_allcores.xlsx")
cores$Sample <- paste0(cores$`Lab Code`, cores$Core)
cores$condition <- "Tumor"
cores$condition[which(cores$`Pathology evaluation` == "normal")] <- "Benign"
```

</details>

## Figure 1C. Cohort.

<details>
<summary>Code</summary>

``` r
cores_gs <- melt(table(cores$`Lab Code`, cores$`Pathology evaluation`))
colnames(cores_gs) <- c("Sample", "GS", "count")
cores_gs$Sample <- factor(cores_gs$Sample)
cores_gs$GS <- gsub("normal", "Benign", cores_gs$GS)
cores_gs$GS <- factor(cores_gs$GS, levels = c("Benign", "3+3", "3+4", "4+3", "4+4", "4+5", "5+4", "5+5"))
ggplot(cores_gs, aes(x=Sample, y=count, fill=GS)) + geom_bar(stat="identity", position="fill", col="black") +
  theme_classic2() + 
  scale_fill_manual(values = color_values(16:1,"piyg")[c(2,9:16)]) + 
  #geom_text(aes(label=count),position = position_fill(vjust = 0.5), data = cores_gs[which(cores_gs$count > 0),], size=2) +
  theme(legend.position = "top") + ylab("Fraction (n = 4/sample)") + theme(axis.text = element_text(size = 7)) 
```

</details>

![](Figure1_files/figure-commonmark/1C-1.png)
