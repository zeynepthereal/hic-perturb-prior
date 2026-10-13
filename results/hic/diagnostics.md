# Hi-C gene graph diagnostics

## Graph structure

| variant         |   nodes_total |   nodes_embedded |   isolated |   edges |   components |   largest_component |   deg_min |   deg_median |   deg_mean |   deg_p95 |   deg_max |
|:----------------|--------------:|-----------------:|-----------:|--------:|-------------:|--------------------:|----------:|-------------:|-----------:|----------:|----------:|
| hic_cis         |         20263 |            20165 |         98 |  135400 |           39 |                1490 |         1 |           12 |      13.43 |        23 |        92 |
| hic_cis_trans   |         20263 |            20256 |          7 |  195476 |            1 |               20256 |         3 |           16 |      19.3  |        36 |       367 |
| linear_distance |         20263 |            20232 |         31 |  119312 |          111 |                 567 |         1 |           11 |      11.79 |        16 |        20 |

## Variance of embedding dimensions explained by chromosome / position (R², mean over dims)

| variant                     | dims     |   R2_chromosome |   R2_chrom+linear_pos |   R2_chrom+poly3_pos |
|:----------------------------|:---------|----------------:|----------------------:|---------------------:|
| hic_cis                     | first 10 |           0.005 |                 0.577 |                0.814 |
| hic_cis                     | first 50 |           0.004 |                 0.24  |                0.63  |
| hic_cis_trans               | first 10 |           0.222 |                 0.361 |                0.507 |
| hic_cis_trans               | first 50 |           0.114 |                 0.192 |                0.321 |
| linear_distance             | first 10 |           0     |                 0.232 |                0.617 |
| linear_distance             | first 50 |           0     |                 0.088 |                0.35  |
| random_gaussian (reference) | first 10 |           0.001 |                 0.002 |                0.005 |

`chromosome` = 23 chromosome indicators; `chrom+linear_pos` adds a per-chromosome linear position term; `chrom+poly3_pos` a per-chromosome cubic. Low-frequency eigenvectors of a path-like graph are cosines of position, so the cubic term is the fair 'linear genomic position' benchmark.

**FLAG hic_cis**: first-10 mean R² of chromosome+cubic position = 0.81 (> 0.5): the embedding is largely chromosome/position.
**FLAG hic_cis_trans**: first-10 mean R² of chromosome+cubic position = 0.51 (> 0.5): the embedding is largely chromosome/position.
**FLAG linear_distance**: first-10 mean R² of chromosome+cubic position = 0.62 (> 0.5): the embedding is largely chromosome/position.

## Similarity of Hi-C and linear-distance constructions

| pair                             |   n_common_nodes |   mean_canonical_corr_10d | top3_canonical        |   min_canonical |
|:---------------------------------|-----------------:|--------------------------:|:----------------------|----------------:|
| hic_cis vs linear_distance       |            20165 |                     0.295 | [0.968, 0.704, 0.656] |           0     |
| hic_cis_trans vs linear_distance |            20230 |                     0.112 | [0.301, 0.277, 0.215] |           0.004 |
| hic_cis_trans vs hic_cis         |            20165 |                     0.141 | [0.422, 0.298, 0.268] |           0.001 |

Graph level (after O/E):

```
{
 "edges_hic_cis": 135400,
 "edges_linear": 119312,
 "edges_shared": 30888,
 "jaccard": 0.13800128672528414,
 "frac_hic_cis_edges_also_in_linear": 0.22812407680945346,
 "spearman_OE_vs_neg_distance_all_candidate_pairs": 0.06363535477075892,
 "spearman_OE_vs_neg_distance_kept_edges": -0.1268602605040259,
 "n_candidate_pairs": 1059092,
 "spearman_weights_on_shared_edges": 0.016064286600997738
}
```

## Predictable GEARS perturbations (eval-set sizes)

| scope                          |   train_pca_expr |   hic_cis |   hic_cis_trans |   linear_distance |   common_pca_and_hic_cis_trans |
|:-------------------------------|-----------------:|----------:|----------------:|------------------:|-------------------------------:|
| all 1092 perturbations         |              412 |      1092 |            1092 |              1092 |                            412 |
| test split gears_seed1 (n=273) |              104 |       273 |             273 |               273 |                            104 |
| test split gears_seed2 (n=273) |              105 |       273 |             273 |               273 |                            105 |
| test split gears_seed3 (n=273) |              105 |       273 |             273 |               273 |                            105 |
| test split gears_seed4 (n=273) |              100 |       273 |             273 |               273 |                            100 |
| test split gears_seed5 (n=273) |              103 |       273 |             273 |               273 |                            103 |

`common` eval set for PCA-vs-Hi-C comparisons = PCA ∩ Hi-C column; `extended:<hic>` = the Hi-C column (vs mean baseline only).
