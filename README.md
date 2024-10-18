# parseggNOG
For parsing the output results of eggNOG, used for clusterprofiler

**Support:**

- KEGG: pathway，EC，Module，reaction

- GO: Cellular Component, Biological Process, Molecular Function
- CAZy

# Install

```r
install.packages("devtools")
devtools::install_github("51cat/parseggNOG")

library(parseggNOG)
```

# load eggNOG output file

```r
df <- load.eggNOG("test_data/eggNOG_output.tsv")

> head(df, 3)
# A tibble: 3 × 21
  Gene         seed_ortholog   evalue score eggNOG_OGs max_annot_lvl COG_category Description Preferred_name GOs   EC    KEGG_ko
  <chr>        <chr>            <dbl> <dbl> <chr>      <chr>         <chr>        <chr>       <chr>          <chr> <chr> <chr>  
1 WP_00315066… 224308.BSU00… 6.10e-40   132 COG1396@1… 91061|Bacilli K            transcript… yazB           -     -     -      
2 WP_00315067… 326423.RBAM_… 1.13e-81   241 COG1539@1… 91061|Bacilli H            Catalyzes … folB           -     1.13… ko:K01…
3 WP_00315068… 326423.RBAM_… 6.44e-50   159 2E9CH@1|r… 91061|Bacilli S            Sigma-K fa… bofA           -     -     ko:K06…
# ℹ 9 more variables: KEGG_Pathway <chr>, KEGG_Module <chr>, KEGG_Reaction <chr>, KEGG_rclass <chr>, BRITE <chr>,
#   KEGG_TC <chr>, CAZy <chr>, BiGG_Reaction <chr>, PFAMs <chr>
```

# kegg

- Pathway

```r
res <- parse.eggNOG.KEGGP("test_data/eggNOG_output.tsv")
# species: A character string specifying the species for pathway filtering
# res <- parse.eggNOG.KEGGP("test_data/eggNOG_output.tsv", species="bac")

# enrich analysis
test.gene <- res$TERM2GENE$Gene |> sample(50) |> unique() # Randomly select 50 unique genes
enrich_res <- clusterProfiler::enricher(
  test.gene,
  TERM2GENE = res$TERM2GENE,
  TERM2NAME = res$TERM2NAME
)

enrichplot::dotplot(enrich_res)
```

- Moduel

```r
res <- parse.eggNOG.KEGGM("test_data/eggNOG_output.tsv")

test.gene <- res$TERM2GENE$Gene |> sample(50) |> unique() # Randomly select 50 unique genes

enrich_res <- clusterProfiler::enricher(
  test.gene,
  pvalueCutoff  = 1,
  TERM2GENE = res$TERM2GENE,
  TERM2NAME = res$TERM2NAME
)

enrichplot::dotplot(enrich_res)
```

- EC

```R
res <- parse.eggNOG.EC("test_data/eggNOG_output.tsv")

test.gene <- res$TERM2GENE$Gene |> sample(50) |> unique() # Randomly select 50 unique genes

enrich_res <- clusterProfiler::enricher(
  test.gene,
  pvalueCutoff  = 1,
  TERM2GENE = res$TERM2GENE,
  TERM2NAME = res$TERM2NAME
)

enrichplot::dotplot(enrich_res)
```

- reaction

```r
res <- parse.eggNOG.EC("test_data/eggNOG_output.tsv")

test.gene <- res$TERM2GENE$Gene |> sample(50) |> unique() # Randomly select 50 unique genes

enrich_res <- clusterProfiler::enricher(
  test.gene,
  pvalueCutoff  = 1,
  TERM2GENE = res$TERM2GENE,
  TERM2NAME = res$TERM2NAME
)

enrichplot::dotplot(enrich_res)
```

# GO

```r
res <- parse.eggNOG.GO("test_data/eggNOG_output.tsv")
test.gene <- res$TERM2GENE$Gene |> sample(50) |> unique() # Randomly select 50 unique genes

enrich_res <- clusterProfiler::enricher(
  test.gene,
  pvalueCutoff  = 1,
  TERM2GENE = res$TERM2GENE.CC., # TERM2GENE.CC or TERM2GENE.BP or TERM2GENE.MF or TERM2GENE.ALL
  TERM2NAME = res$TERM2NAME
)
enrichplot::dotplot(enrich_res)
```

# CAZy

```r
res <- parse.eggNOG.CAZy("test_data/eggNOG_output.tsv")
test.gene <- res$TERM2GENE$Gene |> sample(50) |> unique() # Randomly select 50 unique genes

enrich_res <- clusterProfiler::enricher(
  test.gene,
  pvalueCutoff  = 1,
  TERM2GENE = res$TERM2GENE,
  TERM2NAME = res$TERM2NAME
)
```

