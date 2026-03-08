comb-p pipeline \
    -c 4 \          # p-values in 4th column
    --seed 1e-3 \   # require a p-value of 1e-3 to start a region
    --dist 200      # extend region if find another p-value within this dist
    -p $OUT_PREFIX \
    --region-filter-p 0.1 \ # post-filter reported regions
    --anno mm9 \            # annotate with genome mm9 from UCSC
    $PVALS                  # sorted BED file with pvals in 4th column

# Limma Meta
### seed 0.001
comb-p pipeline -c 4 -p Meta_Females_Random_seed001 --seed 0.001 --dist 750 --region-filter-n 3 --region-filter-p 0.05 ../meta_analysis/Meta_BIP_Female_Random.bed
comb-p pipeline -c 4 -p Meta_Males_Random_seed001 --seed 0.001 --dist 750 --region-filter-n 3 --region-filter-p 0.05 ../meta_analysis/Meta_BIP_Male_Random.bed
comb-p pipeline -c 4 -p Meta_SexAgnostic_Random_seed001 --seed 0.001 --dist 750 --region-filter-n 3 --region-filter-p 0.05 ../meta_analysis/Meta_BIP_SexAgnostic_Random.bed

comb-p pipeline -c 4 -p Meta_Females_Random_seed001_dist500 --seed 0.001 --dist 500 --region-filter-n 3 --region-filter-p 0.05 ../meta_analysis/Meta_BIP_Female_Random.bed
comb-p pipeline -c 4 -p Meta_Males_Random_seed001_dist500 --seed 0.001 --dist 500 --region-filter-n 3 --region-filter-p 0.05 ../meta_analysis/Meta_BIP_Male_Random.bed
comb-p pipeline -c 4 -p Meta_SexAgnostic_Random_seed001_dist500 --seed 0.001 --dist 500 --region-filter-n 3 --region-filter-p 0.05 ../meta_analysis/Meta_BIP_SexAgnostic_Random.bed


