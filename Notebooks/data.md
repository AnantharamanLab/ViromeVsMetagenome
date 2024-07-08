Gather and Organize Data
================
James C. Kosmopoulos
2024-07-08

# Load packages

``` r
library(tidyverse); packageVersion("tidyverse")
```

    ## [1] '2.0.0'

``` r
library(reshape2); packageVersion("reshape2")
```

    ## [1] '1.4.4'

# Load data

``` r
raw_data <- read.csv("../Tables/stats_and_metadata.csv", header = TRUE)
```

# Reformat dataframe and make desired categories

``` r
raw_data$viral_contig_percent_10000 <- (raw_data$viral.scaffold.no / raw_data$X..contigs.....10000.bp.) * 100
raw_data$viral_reads_percent <- (raw_data$sum_paired_reads_mapped_to_viral_contigs / raw_data$Reads_paired_reads) * 100
raw_data$contig_10000_reads_percent <- (raw_data$sum_paired_reads_mapped_to_contigs_10000 / raw_data$Reads_paired_reads) * 100
raw_data$Reads_paired_reads_million <- raw_data$Reads_paired_reads / 1e6
raw_data$lytic_state_na_no <- raw_data$virus.no - (raw_data$int_prophage_no + raw_data$lysogen_virus_no + raw_data$lysogen_scaff_no + raw_data$lytic_virus_no+ raw_data$lytic_scaff_no) # Remaining lytic state unknown

raw_data$lytic_no <- raw_data$lytic_scaff_no + raw_data$lytic_virus_no # Combine lytic virus and scaffold into one lytic category
raw_data$lysogen_no <- raw_data$lysogen_scaff_no + raw_data$lysogen_virus_no  # Combine lysogenic virus and scaffold into one lysogenic category

data_reform <- raw_data %>% 
  mutate(env2 = case_when(Environment =="human_gut"~"Human gut",
                          Environment =="freshwater"~"Freshwater",
                          Environment =="marine"~"Marine",
                          Environment =="soil"~"Soil"
                      )) %>%
  mutate(method2 = case_when(Method == "metagenome"~"Mixed MG",
                             Method == "virome"~"Virome",
                            ))

data_reform$lytic_scaff_no_norm_VMAGs <- data_reform$lytic_scaff_no / data_reform$virus.no
data_reform$lytic_virus_no_norm_VMAGs <- data_reform$lytic_virus_no / data_reform$virus.no
data_reform$lysogen_scaff_no_norm_VMAGs <- data_reform$lysogen_scaff_no / data_reform$virus.no
data_reform$lysogen_virus_no_norm_VMAGs <- data_reform$lysogen_virus_no / data_reform$virus.no
data_reform$int_prophage_no_norm_VMAGs <- data_reform$int_prophage_no / data_reform$virus.no
data_reform$lytic_state_na_no_norm_VMAGs <- data_reform$lytic_state_na_no / data_reform$virus.no
data_reform$checkv_low_norm_VMAGs <- data_reform$checkv_low / data_reform$virus.no
data_reform$checkv_med_norm_VMAGs <- data_reform$checkv_med / data_reform$virus.no
data_reform$checkv_high_norm_VMAGs <- data_reform$checkv_high / data_reform$virus.no
data_reform$checkv_na_norm_VMAGs <- data_reform$checkv_na / data_reform$virus.no
data_reform$lytic_no_norm_vmags <- (data_reform$lytic_scaff_no + data_reform$lytic_virus_no) / data_reform$virus.no
data_reform$lytic_no_norm_reads <- (data_reform$lytic_scaff_no + data_reform$lytic_virus_no) / data_reform$X100M_paired_reads
data_reform$lysogen_no_norm_vmags <- (data_reform$lysogen_scaff_no + data_reform$lysogen_virus_no) / data_reform$virus.no
data_reform$lysogen_no_norm_reads <- (data_reform$lysogen_scaff_no + data_reform$lysogen_virus_no) / data_reform$X100M_paired_reads
data_reform$int_prophage_no_norm_vmags <- data_reform$int_prophage_no / data_reform$virus.no
data_reform$int_prophage_no_norm_reads <- data_reform$int_prophage_no / data_reform$X100M_paired_reads
data_reform$lytic_state_na_no_norm_vmags <- data_reform$lytic_state_na_no / data_reform$virus.no
data_reform$lytic_state_na_no_norm_reads <- data_reform$lytic_state_na_no / data_reform$X100M_paired_reads

data_reform[is.na(data_reform)] <- 0 # Convert NA values to 0
data_reform$method2 <- factor(data_reform$method2, levels=c("Mixed MG", "Virome"))

head(data_reform)
```

    ##       Sample sample_source Environment     Method Avg_paired_reads
    ## 1 SRR9161505    APC055_916   human_gut     virome            281.4
    ## 2 SRR9162904    APC055_916   human_gut metagenome            149.8
    ## 3 SRR9161504    APC055_917   human_gut     virome            276.0
    ## 4 SRR9162903    APC055_917   human_gut metagenome            149.7
    ## 5 SRR9161507    APC055_918   human_gut     virome            289.2
    ## 6 SRR9162906    APC055_918   human_gut metagenome            149.8
    ##   Bases_paired_reads Max_paired_reads Median_paired_reads Min_paired_reads
    ## 1         4161259362              300                 300               93
    ## 2          735582373              150                 150               51
    ## 3         4720624501              300                 300               94
    ## 4          661566686              150                 150               51
    ## 5         4282250065              300                 300               94
    ## 6          784243722              150                 150               51
    ##   Mode_paired_reads Reads_paired_reads Std_Dev_paired_reads
    ## 1               300           14788544                 33.3
    ## 2               150            4909488                  2.7
    ## 3               300           17102848                 35.9
    ## 4               150            4418516                  3.4
    ## 5               300           14808998                 26.3
    ## 6               150            5236862                  3.3
    ##   X..contigs.....0.bp. X..contigs.....2000.bp. X..contigs.....5000.bp.
    ## 1               195247                   13654                    3534
    ## 2                90373                    6203                    1352
    ## 3               170002                   12035                    3135
    ## 4                84992                    4820                    1545
    ## 5               120975                   11399                    2898
    ## 6               103476                    4419                    1000
    ##   X..contigs.....10000.bp. X..contigs.....25000.bp. X..contigs.....50000.bp.
    ## 1                     1489                      495                      190
    ## 2                      442                      140                       46
    ## 3                     1339                      476                      178
    ## 4                      614                      152                       63
    ## 5                     1114                      486                      265
    ## 6                      535                      223                       96
    ##   Total.length.....0.bp. Total.length.....2000.bp. Total.length.....5000.bp.
    ## 1              217108812                  87843191                  58337795
    ## 2               81458979                  32134071                  17825580
    ## 3              193515957                  77386482                  51631395
    ## 4               75606946                  31939177                  22270022
    ## 5              168572364                  87215534                  62057307
    ## 6               89937920                  31769357                  22124734
    ##   Total.length.....10000.bp. Total.length.....25000.bp.
    ## 1                   44494048                   29374791
    ## 2                   11761781                    7072004
    ## 3                   39420061                   26432307
    ## 4                   15915970                    8951144
    ## 5                   49993244                   40570009
    ## 6                   18891075                   13978638
    ##   Total.length.....50000.bp. X..contigs Largest.contig Total.length   N50  N90
    ## 1                   18608150      13654         481186     87843191 10324 2483
    ## 2                    3866591       6203         237571     32134071  5859 2384
    ## 3                   16107803      12035         496893     77386482 10453 2470
    ## 4                    5959624       4820         312785     31939177  9938 2585
    ## 5                   32648308      11399         881056     87215534 18449 2636
    ## 6                    9552820       4419         498610     31769357 18285 2469
    ##       auN  L50  L90 X..N.s.per.100.kbp viral.scaffold.no virus.no
    ## 1 37597.0 1433 9698                  0               183      164
    ## 2 21229.3 1027 4730                  0                66       65
    ## 3 34647.4 1268 8530                  0               154      126
    ## 4 32025.1  620 3407                  0                72       70
    ## 5 80547.4  630 7579                  0               191      175
    ## 6 54327.8  312 2977                  0                65       60
    ##   species.cluster.no genus.cluster.no no.of.virus.taxonomy.info
    ## 1                147              124                        54
    ## 2                 48               58                        10
    ## 3                107              102                        37
    ## 4                 56               60                        22
    ## 5                146              140                        49
    ## 6                 42               52                         6
    ##   no.of.virus.with.host.prediction lytic_scaff_no lytic_virus_no
    ## 1                                0             77              6
    ## 2                                0             36              0
    ## 3                                0             43              8
    ## 4                                0             35              1
    ## 5                                0             81              2
    ## 6                                0             26              2
    ##   lysogen_scaff_no lysogen_virus_no int_prophage_no checkv_low checkv_med
    ## 1               25                5              51         82         24
    ## 2                7                1              21         41          5
    ## 3               20                2              53         53         21
    ## 4               11                0              23         42          5
    ## 5               33                2              57         75         24
    ## 6                8                0              24         17          6
    ##   checkv_high checkv_na sum_paired_reads_mapped_to_viral_contigs
    ## 1          19        35                                   843299
    ## 2           2        17                                    77219
    ## 3          13        32                                   755959
    ## 4           3        19                                    74390
    ## 5          11        45                                  2359024
    ## 6           3        33                                   127500
    ##   sum_paired_reads_mapped_to_contigs_10000 X100M_paired_reads
    ## 1                                  7473888         0.14788544
    ## 2                                  1195205         0.04909488
    ## 3                                  9611588         0.17102848
    ## 4                                  1428781         0.04418516
    ## 5                                 10379345         0.14808998
    ## 6                                  2008971         0.05236862
    ##   X..N.s.per.100.kbp_norm X..contigs_norm X..contigs.....0.bp._norm
    ## 1                       0        92328.22                 1320258.4
    ## 2                       0       126347.19                 1840782.6
    ## 3                       0        70368.40                  993998.2
    ## 4                       0       109086.40                 1923541.8
    ## 5                       0        76973.47                  816902.0
    ## 6                       0        84382.59                 1975916.1
    ##   X..contigs.....10000.bp._norm X..contigs.....2000.bp._norm
    ## 1                     10068.604                     92328.22
    ## 2                      9002.975                    126347.19
    ## 3                      7829.105                     70368.40
    ## 4                     13896.068                    109086.40
    ## 5                      7522.454                     76973.47
    ## 6                     10216.042                     84382.59
    ##   X..contigs.....25000.bp._norm X..contigs.....5000.bp._norm
    ## 1                      3347.185                     23896.88
    ## 2                      2851.621                     27538.51
    ## 3                      2783.162                     18330.28
    ## 4                      3440.069                     34966.49
    ## 5                      3281.789                     19569.18
    ## 6                      4258.275                     19095.40
    ##   X..contigs.....50000.bp._norm Avg_paired_reads_norm Bases_paired_reads_norm
    ## 1                     1284.7783              1902.824             28138397952
    ## 2                      936.9612              3051.235             14982873428
    ## 3                     1040.7623              1613.766             27601394230
    ## 4                     1425.8181              3388.015             14972599081
    ## 5                     1789.4526              1952.867             28916541585
    ## 6                     1833.1589              2860.492             14975451368
    ##    L50_norm L90_norm Largest.contig_norm Max_paired_reads_norm
    ## 1  9689.933 65577.79             3253775              2028.597
    ## 2 20918.678 96344.06             4839018              3055.308
    ## 3  7413.970 49874.73             2905323              1754.094
    ## 4 14031.860 77107.34             7078960              3394.805
    ## 5  4254.170 51178.34             5949464              2025.795
    ## 6  5957.766 56847.02             9521160              2864.311
    ##   Median_paired_reads_norm Min_paired_reads_norm Mode_paired_reads_norm
    ## 1                 2028.597              628.8652               2028.597
    ## 2                 3055.308             1038.8049               3055.308
    ## 3                 1754.094              549.6161               1754.094
    ## 4                 3394.805             1154.2337               3394.805
    ## 5                 2025.795              634.7492               2025.795
    ## 6                 2864.311              973.8656               2864.311
    ##    N50_norm N90_norm Reads_paired_reads_norm Std_Dev_paired_reads_norm
    ## 1  69810.79 16790.02               100000000                 225.17430
    ## 2 119340.35 48559.04               100000000                  54.99555
    ## 3  61118.48 14442.04               100000000                 209.90656
    ## 4 224917.14 58503.81               100000000                  76.94891
    ## 5 124579.66 17799.99               100000000                 177.59473
    ## 6 349159.48 47146.55               100000000                  63.01484
    ##   Total.length_norm Total.length.....0.bp._norm Total.length.....10000.bp._norm
    ## 1         593994858                  1468087812                       300868348
    ## 2         654529984                  1659215360                       239572456
    ## 3         452477166                  1131483815                       230488285
    ## 4         722848508                  1711138898                       360210758
    ## 5         588936091                  1138310398                       337586946
    ## 6         606648734                  1717400993                       360732725
    ##   Total.length.....2000.bp._norm Total.length.....25000.bp._norm
    ## 1                      593994858                       198632070
    ## 2                      654529984                       144047689
    ## 3                      452477166                       154549155
    ## 4                      722848508                       202582587
    ## 5                      588936091                       273955125
    ## 6                      606648734                       266927752
    ##   Total.length.....5000.bp._norm Total.length.....50000.bp._norm  auN_norm
    ## 1                      394479639                       125828141  254230.6
    ## 2                      363084297                        78757520  432413.7
    ## 3                      301887703                        94181992  202582.6
    ## 4                      504015873                       134878407  724793.1
    ## 5                      419051356                       220462640  543908.5
    ## 6                      422480753                       182414965 1037411.3
    ##   checkv_high_norm checkv_low_norm checkv_med_norm checkv_na_norm
    ## 1        128.47783        554.4833        162.2878       236.6697
    ## 2         40.73745        835.1176        101.8436       346.2683
    ## 3         76.01073        309.8899        122.7866       187.1033
    ## 4         67.89610        950.5454        113.1602       430.0086
    ## 5         74.27916        506.4488        162.0636       303.8693
    ## 6         57.28621        324.6219        114.5724       630.1484
    ##   genus.cluster.no_norm int_prophage_no_norm lysogen_scaff_no_norm
    ## 1              838.4869             344.8615              169.0498
    ## 2             1181.3859             427.7432              142.5811
    ## 3              596.3919             309.8899              116.9396
    ## 4             1357.9220             520.5368              248.9524
    ## 5              945.3712             384.9011              222.8375
    ## 6              992.9611             458.2897              152.7632
    ##   lysogen_virus_no_norm lytic_scaff_no_norm lytic_virus_no_norm
    ## 1              33.80995            520.6733            40.57195
    ## 2              20.36872            733.2740             0.00000
    ## 3              11.69396            251.4201            46.77584
    ## 4               0.00000            792.1212            22.63203
    ## 5              13.50530            546.9648            13.50530
    ## 6               0.00000            496.4805            38.19081
    ##   no.of.virus.taxonomy.info_norm no.of.virus.with.host.prediction_norm
    ## 1                       365.1475                                     0
    ## 2                       203.6872                                     0
    ## 3                       216.3382                                     0
    ## 4                       497.9047                                     0
    ## 5                       330.8799                                     0
    ## 6                       114.5724                                     0
    ##   species.cluster.no_norm viral.scaffold.no_norm virus.no_norm
    ## 1                994.0127              1237.4443     1108.9665
    ## 2                977.6987              1344.3357     1323.9670
    ## 3                625.6268               900.4348      736.7194
    ## 4               1267.3938              1629.5064     1584.2423
    ## 5                985.8871              1289.7564     1181.7140
    ## 6                802.0070              1241.2013     1145.7243
    ##   viral_contig_percent_10000 viral_reads_percent contig_10000_reads_percent
    ## 1                   12.29013            5.702380                   50.53836
    ## 2                   14.93213            1.572852                   24.34480
    ## 3                   11.50112            4.420077                   56.19876
    ## 4                   11.72638            1.683597                   32.33622
    ## 5                   17.14542           15.929667                   70.08810
    ## 6                   12.14953            2.434664                   38.36211
    ##   Reads_paired_reads_million lytic_state_na_no lytic_no lysogen_no      env2
    ## 1                  14.788544                 0       83         30 Human gut
    ## 2                   4.909488                 0       36          8 Human gut
    ## 3                  17.102848                 0       51         22 Human gut
    ## 4                   4.418516                 0       36         11 Human gut
    ## 5                  14.808998                 0       83         35 Human gut
    ## 6                   5.236862                 0       28          8 Human gut
    ##    method2 lytic_scaff_no_norm_VMAGs lytic_virus_no_norm_VMAGs
    ## 1   Virome                 0.4695122                0.03658537
    ## 2 Mixed MG                 0.5538462                0.00000000
    ## 3   Virome                 0.3412698                0.06349206
    ## 4 Mixed MG                 0.5000000                0.01428571
    ## 5   Virome                 0.4628571                0.01142857
    ## 6 Mixed MG                 0.4333333                0.03333333
    ##   lysogen_scaff_no_norm_VMAGs lysogen_virus_no_norm_VMAGs
    ## 1                   0.1524390                  0.03048780
    ## 2                   0.1076923                  0.01538462
    ## 3                   0.1587302                  0.01587302
    ## 4                   0.1571429                  0.00000000
    ## 5                   0.1885714                  0.01142857
    ## 6                   0.1333333                  0.00000000
    ##   int_prophage_no_norm_VMAGs lytic_state_na_no_norm_VMAGs checkv_low_norm_VMAGs
    ## 1                  0.3109756                            0             0.5000000
    ## 2                  0.3230769                            0             0.6307692
    ## 3                  0.4206349                            0             0.4206349
    ## 4                  0.3285714                            0             0.6000000
    ## 5                  0.3257143                            0             0.4285714
    ## 6                  0.4000000                            0             0.2833333
    ##   checkv_med_norm_VMAGs checkv_high_norm_VMAGs checkv_na_norm_VMAGs
    ## 1            0.14634146             0.11585366            0.2134146
    ## 2            0.07692308             0.03076923            0.2615385
    ## 3            0.16666667             0.10317460            0.2539683
    ## 4            0.07142857             0.04285714            0.2714286
    ## 5            0.13714286             0.06285714            0.2571429
    ## 6            0.10000000             0.05000000            0.5500000
    ##   lytic_no_norm_vmags lytic_no_norm_reads lysogen_no_norm_vmags
    ## 1           0.5060976            561.2452             0.1829268
    ## 2           0.5538462            733.2740             0.1230769
    ## 3           0.4047619            298.1959             0.1746032
    ## 4           0.5142857            814.7532             0.1571429
    ## 5           0.4742857            560.4701             0.2000000
    ## 6           0.4666667            534.6713             0.1333333
    ##   lysogen_no_norm_reads int_prophage_no_norm_vmags int_prophage_no_norm_reads
    ## 1              202.8597                  0.3109756                   344.8615
    ## 2              162.9498                  0.3230769                   427.7432
    ## 3              128.6335                  0.4206349                   309.8899
    ## 4              248.9524                  0.3285714                   520.5368
    ## 5              236.3428                  0.3257143                   384.9011
    ## 6              152.7632                  0.4000000                   458.2897
    ##   lytic_state_na_no_norm_vmags lytic_state_na_no_norm_reads
    ## 1                            0                            0
    ## 2                            0                            0
    ## 3                            0                            0
    ## 4                            0                            0
    ## 5                            0                            0
    ## 6                            0                            0

``` r
saveRDS(data_reform, file = "../Data/stats_and_metadata.RDS")
```

## Table S1

``` r
tableS1 <- data_reform[, c("Sample", "sample_source", "Environment", "Method", "Avg_paired_reads", "Bases_paired_reads", "Reads_paired_reads", "X..contigs.....0.bp.", "X..contigs.....2000.bp.", "X..contigs.....5000.bp.", "X..contigs.....10000.bp.", "X..contigs.....25000.bp.", "X..contigs.....50000.bp.", "Total.length.....0.bp.", "Total.length.....2000.bp.", "Total.length.....5000.bp.", "Total.length.....10000.bp.", "Total.length.....25000.bp.", "Total.length.....50000.bp.", "X..contigs", "Total.length", "N50", "N90", "L50", "L90", "viral.scaffold.no", "virus.no", "species.cluster.no", "genus.cluster.no", "no.of.virus.taxonomy.info", "lytic_no", "lysogen_no", "int_prophage_no", "checkv_low", "checkv_med", "checkv_high", "checkv_na")]
colnames(tableS1) <- c("Sample", "Sample source",   "Environment", "Method", "Average read length", "Total bases",  "Total read pairs", "Number of contigs (>= 0 bp)",  "Number of contigs (>= 2000 bp)",   "Number of contigs (>= 5000 bp)",   "Number of contigs (>= 10000 bp)",  "Number of contigs (>= 25000 bp)",  "Number of contigs (>= 50000 bp)", "Total length (>= 0 bp)",    "Total length (>= 2000 bp)",    "Total length (>= 5000 bp)",    "Total length (>= 10000 bp)",   "Total length (>= 25000 bp)",   "Total length (>= 50000 bp)",   "Number of contigs",    "Total length", "N50",  "N90",  "L50",  "L90",  "Number of viral contigs",  "Number of VMAGs",  "Number of species clusters",   "Number of genus clusters", "Number of VMAGs with taxonomy information",    "Number of predicted lytic VMAGs",  "Number of predicted lysogenic VMAGs",  "Number of predicted integrated prophage VMAGs",    "Number of CheckV low-quality VMAGs",   "Number of CheckV medium-quality VMAGs",    "Number of CheckV high-quality VMAGs", "Number of CheckV not-determined")
tableS1 <- tableS1 %>% 
  mutate(Environment = case_when(Environment =="human_gut"~"Human gut",
                          Environment =="freshwater"~"Freshwater",
                          Environment =="marine"~"Marine",
                          Environment =="soil"~"Soil"
                      )) %>%
  mutate(Method = case_when(Method == "metagenome"~"Mixed MG",
                             Method == "virome"~"Virome",
                            ))
write.csv(tableS1, file = "../Tables/TableS1.csv")
saveRDS(tableS1, file = "../Data/TableS1.RDS")
head(tableS1)
```

    ##       Sample Sample source Environment   Method Average read length Total bases
    ## 1 SRR9161505    APC055_916   Human gut   Virome               281.4  4161259362
    ## 2 SRR9162904    APC055_916   Human gut Mixed MG               149.8   735582373
    ## 3 SRR9161504    APC055_917   Human gut   Virome               276.0  4720624501
    ## 4 SRR9162903    APC055_917   Human gut Mixed MG               149.7   661566686
    ## 5 SRR9161507    APC055_918   Human gut   Virome               289.2  4282250065
    ## 6 SRR9162906    APC055_918   Human gut Mixed MG               149.8   784243722
    ##   Total read pairs Number of contigs (>= 0 bp) Number of contigs (>= 2000 bp)
    ## 1         14788544                      195247                          13654
    ## 2          4909488                       90373                           6203
    ## 3         17102848                      170002                          12035
    ## 4          4418516                       84992                           4820
    ## 5         14808998                      120975                          11399
    ## 6          5236862                      103476                           4419
    ##   Number of contigs (>= 5000 bp) Number of contigs (>= 10000 bp)
    ## 1                           3534                            1489
    ## 2                           1352                             442
    ## 3                           3135                            1339
    ## 4                           1545                             614
    ## 5                           2898                            1114
    ## 6                           1000                             535
    ##   Number of contigs (>= 25000 bp) Number of contigs (>= 50000 bp)
    ## 1                             495                             190
    ## 2                             140                              46
    ## 3                             476                             178
    ## 4                             152                              63
    ## 5                             486                             265
    ## 6                             223                              96
    ##   Total length (>= 0 bp) Total length (>= 2000 bp) Total length (>= 5000 bp)
    ## 1              217108812                  87843191                  58337795
    ## 2               81458979                  32134071                  17825580
    ## 3              193515957                  77386482                  51631395
    ## 4               75606946                  31939177                  22270022
    ## 5              168572364                  87215534                  62057307
    ## 6               89937920                  31769357                  22124734
    ##   Total length (>= 10000 bp) Total length (>= 25000 bp)
    ## 1                   44494048                   29374791
    ## 2                   11761781                    7072004
    ## 3                   39420061                   26432307
    ## 4                   15915970                    8951144
    ## 5                   49993244                   40570009
    ## 6                   18891075                   13978638
    ##   Total length (>= 50000 bp) Number of contigs Total length   N50  N90  L50
    ## 1                   18608150             13654     87843191 10324 2483 1433
    ## 2                    3866591              6203     32134071  5859 2384 1027
    ## 3                   16107803             12035     77386482 10453 2470 1268
    ## 4                    5959624              4820     31939177  9938 2585  620
    ## 5                   32648308             11399     87215534 18449 2636  630
    ## 6                    9552820              4419     31769357 18285 2469  312
    ##    L90 Number of viral contigs Number of VMAGs Number of species clusters
    ## 1 9698                     183             164                        147
    ## 2 4730                      66              65                         48
    ## 3 8530                     154             126                        107
    ## 4 3407                      72              70                         56
    ## 5 7579                     191             175                        146
    ## 6 2977                      65              60                         42
    ##   Number of genus clusters Number of VMAGs with taxonomy information
    ## 1                      124                                        54
    ## 2                       58                                        10
    ## 3                      102                                        37
    ## 4                       60                                        22
    ## 5                      140                                        49
    ## 6                       52                                         6
    ##   Number of predicted lytic VMAGs Number of predicted lysogenic VMAGs
    ## 1                              83                                  30
    ## 2                              36                                   8
    ## 3                              51                                  22
    ## 4                              36                                  11
    ## 5                              83                                  35
    ## 6                              28                                   8
    ##   Number of predicted integrated prophage VMAGs
    ## 1                                            51
    ## 2                                            21
    ## 3                                            53
    ## 4                                            23
    ## 5                                            57
    ## 6                                            24
    ##   Number of CheckV low-quality VMAGs Number of CheckV medium-quality VMAGs
    ## 1                                 82                                    24
    ## 2                                 41                                     5
    ## 3                                 53                                    21
    ## 4                                 42                                     5
    ## 5                                 75                                    24
    ## 6                                 17                                     6
    ##   Number of CheckV high-quality VMAGs Number of CheckV not-determined
    ## 1                                  19                              35
    ## 2                                   2                              17
    ## 3                                  13                              32
    ## 4                                   3                              19
    ## 5                                  11                              45
    ## 6                                   3                              33

## Contig stats

``` r
contig_stats <- data_reform[, c("Sample", "sample_source", "Environment", "Method", "env2", "method2")]
contig_stats <- merge(contig_stats, data_reform[,c(1,16:35)], by="Sample")
contig_stats_melt <- melt(contig_stats[, c("Sample", "sample_source", "Environment", "env2", "Method", "method2", "N50", "N90", "L50", "L90")], id = c("Sample", "sample_source", "Environment", "env2", "Method", "method2"))
contig_stats_melt$method2 <- factor(contig_stats_melt$method2, levels=c("Mixed MG", "Virome"))
contig_stats_melt$Method <- factor(contig_stats_melt$Method, levels=c("metagenome", "virome"))

saveRDS(contig_stats_melt, file = "../Data/contig_stats_melt.RDS")
head(contig_stats_melt)
```

    ##      Sample sample_source Environment   env2 Method method2 variable value
    ## 1 ERR594353  TARA_070_SRF      marine Marine virome  Virome      N50  8294
    ## 2 ERR594354  TARA_076_SRF      marine Marine virome  Virome      N50  9912
    ## 3 ERR594355  TARA_076_DCM      marine Marine virome  Virome      N50  9253
    ## 4 ERR594359  TARA_065_SRF      marine Marine virome  Virome      N50  6416
    ## 5 ERR594362  TARA_066_SRF      marine Marine virome  Virome      N50 13218
    ## 6 ERR594364  TARA_072_SRF      marine Marine virome  Virome      N50  6995

## Virus info

``` r
vir_info <- read.csv("../Tables/virus_summary_info_combined.csv", header=TRUE)
# Combine "scaffold" and "virus" distinctions
vir_info <- vir_info %>% mutate(lytic_state = case_when(lytic_state == "lytic_scaffold" ~ "Lytic",
                                                        lytic_state == "lytic_virus" ~ "Lytic",
                                                        lytic_state == "lysogenic_scaffold" ~ "Lysogenic",
                                                        lytic_state == "lysogenic_virus" ~ "Lysogenic",
                                                        lytic_state == "integrated_prophage" ~ "Int. prophage"))
vir_info$sample <- gsub("_contigs", "", vir_info$sample)
vir_info$genome <- gsub("_contigs", "", vir_info$genome)
saveRDS(vir_info, file= "../Data/virus_summary_info_combined.RDS")
head(vir_info)
```

    ##       sample                         genome genome_size scaffold_num
    ## 1 SRR9162906 SRR9162906__vRhyme_unbinned_58        9211            1
    ## 2 SRR9162906 SRR9162906__vRhyme_unbinned_54       24741            1
    ## 3 SRR9162906 SRR9162906__vRhyme_unbinned_24       26135            1
    ## 4 SRR9162906 SRR9162906__vRhyme_unbinned_57       12235            1
    ## 5 SRR9162906 SRR9162906__vRhyme_unbinned_20       14202            1
    ## 6 SRR9162906 SRR9162906__vRhyme_unbinned_12       17764            1
    ##   protein_count   AMG_KOs   lytic_state checkv_quality  miuvig_quality
    ## 1            11           Int. prophage    Low-quality Genome-fragment
    ## 2            28 K20170(1) Int. prophage    Low-quality Genome-fragment
    ## 3            40 K00558(1) Int. prophage Medium-quality Genome-fragment
    ## 4            11                   Lytic Not-determined Genome-fragment
    ## 5            18 K00655(1) Int. prophage   High-quality    High-quality
    ## 6            17                   Lytic Not-determined Genome-fragment
    ##   completeness         completeness_method
    ## 1         3.70     HMM-based (lower-bound)
    ## 2        15.55     HMM-based (lower-bound)
    ## 3        67.69 AAI-based (high-confidence)
    ## 4           NA                            
    ## 5       100.00     HMM-based (lower-bound)
    ## 6           NA

## Virus species cluster representatives

``` r
species_reps <- read.csv("../Tables/drep_cluster_reps.csv", header=TRUE) %>%
  rename(`Species cluster representative` = Species.cluster.representative)
saveRDS(species_reps, file= "../Data/species_reps.RDS")
head(species_reps)
```

    ##      Sample                    Genome Species cluster representative
    ## 1 ERR594353   ERR594353__vRhyme_bin_1        ERR594353__vRhyme_bin_1
    ## 2 ERR594353  ERR594353__vRhyme_bin_10       ERR594353__vRhyme_bin_10
    ## 3 ERR594353 ERR594353__vRhyme_bin_100      ERR594353__vRhyme_bin_100
    ## 4 ERR594353 ERR594353__vRhyme_bin_101      ERR594353__vRhyme_bin_101
    ## 5 ERR594353 ERR594353__vRhyme_bin_102      ERR594353__vRhyme_bin_102
    ## 6 ERR594353 ERR594353__vRhyme_bin_103      ERR594353__vRhyme_bin_103

## Virus genome completeness

Filter the dataframe such that: - There are no NA values for
completeness - Each species cluster representative is has member genomes
in both Viromes and Mixed MGs for each sample source

``` r
completeness <- vir_info %>%
  select(sample, genome, completeness) %>%
  rename(Sample = sample, Genome = genome, Completeness = completeness) %>%
  inner_join(data_reform %>%
               select(Sample, sample_source, Environment, Method)) %>%
  inner_join(species_reps %>%
               select(Genome, `Species cluster representative`)) %>%
  mutate(Method = case_when(Method == "virome" ~ "Virome",
                            Method == "metagenome" ~ "Mixed MG")) %>%
  mutate(Method = factor(Method, levels = c("Mixed MG", "Virome"))) %>%
  mutate(Environment = case_when(Environment == "adult_gut" ~ "Human gut",
                                 Environment == "freshwater" ~ "Freshwater",
                                 Environment == "marine" ~ "Marine",
                                 Environment == "soil" ~ "Soil")) %>%
  filter(is.na(Completeness) == F)

completeness <- completeness %>%
  group_by(`Species cluster representative`, sample_source) %>%
  filter(n_distinct(Method) == 2) %>%
  ungroup()

saveRDS(completeness, file= "../Data/virus_genome_completeness.RDS")
head(completeness)
```

    ## # A tibble: 6 × 7
    ##   Sample    Genome                 Completeness sample_source Environment Method
    ##   <chr>     <chr>                         <dbl> <chr>         <chr>       <fct> 
    ## 1 ERR599148 ERR599148__vRhyme_unb…        100   TARA_076_DCM  Marine      Mixed…
    ## 2 ERR599148 ERR599148__vRhyme_unb…         50.9 TARA_076_DCM  Marine      Mixed…
    ## 3 ERR599148 ERR599148__vRhyme_unb…        100   TARA_076_DCM  Marine      Mixed…
    ## 4 ERR599148 ERR599148__vRhyme_unb…         31.5 TARA_076_DCM  Marine      Mixed…
    ## 5 ERR599148 ERR599148__vRhyme_unb…         18.9 TARA_076_DCM  Marine      Mixed…
    ## 6 ERR599148 ERR599148__vRhyme_unb…         95.0 TARA_076_DCM  Marine      Mixed…
    ## # ℹ 1 more variable: `Species cluster representative` <chr>

## Virus taxonomy

``` r
vir_tax <- read.csv("../Tables/virus_tax_classification_results_combined.csv", header=TRUE)
vir_tax <- separate(data = vir_tax, col = taxonomy, into = c("realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = "\\;")
vir_tax$sample <- gsub("_contigs", "", vir_tax$sample)
vir_tax$genome <- gsub("_contigs", "", vir_tax$genome)
# Remove lingering rank labels
vir_tax <- vir_tax %>% mutate_all(funs(str_replace(., "r__", "")))
vir_tax <- vir_tax %>% mutate_all(funs(str_replace(., "k__", "")))
vir_tax <- vir_tax %>% mutate_all(funs(str_replace(., "p__", "")))
vir_tax <- vir_tax %>% mutate_all(funs(str_replace(., "c__", "")))
vir_tax <- vir_tax %>% mutate_all(funs(str_replace(., "o__", "")))
vir_tax <- vir_tax %>% mutate_all(funs(str_replace(., "f__", "")))
vir_tax <- vir_tax %>% mutate_all(funs(str_replace(., "g__", "")))
vir_tax <- vir_tax %>% mutate_all(funs(str_replace(., "s__", "")))
# Add metadata columns
cols <- data_reform[, c("Sample", "sample_source", "Environment", "Method", "env2", "method2")]
names(cols)[names(cols) == 'Sample'] <- 'sample'
vir_tax <- merge(vir_tax, cols, by="sample", all.x = TRUE)
saveRDS(vir_tax, file= "../Data/virus_tax_classification_results_combined.RDS")
head(vir_tax)
```

    ##      sample                          genome         realm        kingdom
    ## 1 ERR594353       ERR594353__vRhyme_bin_101 Duplodnaviria Heunggongvirae
    ## 2 ERR594353        ERR594353__vRhyme_bin_87 Duplodnaviria Heunggongvirae
    ## 3 ERR594353        ERR594353__vRhyme_bin_13 Duplodnaviria Heunggongvirae
    ## 4 ERR594353 ERR594353__vRhyme_unbinned_1061 Duplodnaviria Heunggongvirae
    ## 5 ERR594353   ERR594353__vRhyme_unbinned_66 Duplodnaviria Heunggongvirae
    ## 6 ERR594353  ERR594353__vRhyme_unbinned_286 Duplodnaviria Heunggongvirae
    ##        phylum          class order       family      genus
    ## 1 Uroviricota Caudoviricetes    NA           NA         NA
    ## 2 Uroviricota Caudoviricetes    NA           NA         NA
    ## 3 Uroviricota Caudoviricetes    NA Kyanoviridae Haifavirus
    ## 4 Uroviricota Caudoviricetes    NA           NA         NA
    ## 5 Uroviricota Caudoviricetes    NA           NA         NA
    ## 6 Uroviricota Caudoviricetes    NA Kyanoviridae Haifavirus
    ##                           species                              source
    ## 1     Pelagibacter phage HTVC008M NCBI RefSeq viral protein searching
    ## 2 Puniceispirillum phage HMO-2011 NCBI RefSeq viral protein searching
    ## 3   Prochlorococcus phage P-TIM68 NCBI RefSeq viral protein searching
    ## 4     Pelagibacter phage HTVC008M NCBI RefSeq viral protein searching
    ## 5     Pelagibacter phage HTVC008M NCBI RefSeq viral protein searching
    ## 6   Prochlorococcus phage P-TIM68 NCBI RefSeq viral protein searching
    ##   sample_source Environment Method   env2 method2
    ## 1  TARA_070_SRF      marine virome Marine  Virome
    ## 2  TARA_070_SRF      marine virome Marine  Virome
    ## 3  TARA_070_SRF      marine virome Marine  Virome
    ## 4  TARA_070_SRF      marine virome Marine  Virome
    ## 5  TARA_070_SRF      marine virome Marine  Virome
    ## 6  TARA_070_SRF      marine virome Marine  Virome

## Genome breadths (covered fractions)

``` r
breadth.gut <- read.csv("../Tables/covered_fraction_human_gut.csv", header=TRUE, row.names = 1)
breadth.fw <- read.csv("../Tables/covered_fraction_freshwater.csv", header=TRUE, row.names = 1)
breadth.mar <- read.csv("../Tables/covered_fraction_marine.csv", header=TRUE, row.names = 1)
breadth.soil <- read.csv("../Tables/covered_fraction_soil.csv", header=TRUE, row.names = 1)
```

### Fix rownames, add an extra “\_” between sample ID and bin name

``` r
fix_rownames <- function(x) {
  gsub("_vRhyme", "__vRhyme", x)
}

rownames(breadth.gut) <- fix_rownames(rownames(breadth.gut))
head(breadth.gut)
```

    ##                                 SRR9161502 SRR9161504 SRR9161503 SRR9162903
    ## SRR9161506__vRhyme_unbinned_143  0.6467704  0.3172089 0.61669020  0.1833098
    ## SRR9161505__vRhyme_unbinned_108  0.3009803  0.4573876 0.91389210  0.8856310
    ## SRR9161504__vRhyme_unbinned_104  0.0000000  1.0000000 0.00000000  0.1834470
    ## SRR9161504__vRhyme_bin_20        0.2757188  1.0000000 0.00000000  0.7504299
    ## SRR9162903__vRhyme_unbinned_28   0.7433444  0.8978298 0.01579839  0.9999453
    ## SRR9161506__vRhyme_unbinned_183  0.0000000  0.0000000 0.00000000  0.0000000
    ##                                  SRR9161510 SRR9161505 SRR9162901 SRR9162902
    ## SRR9161506__vRhyme_unbinned_143 0.028288543  0.5354078 0.33738804 0.32032060
    ## SRR9161505__vRhyme_unbinned_108 0.000000000  0.9773028 0.02649475 0.85825310
    ## SRR9161504__vRhyme_unbinned_104 0.000000000  0.0000000 0.00000000 0.00000000
    ## SRR9161504__vRhyme_bin_20       0.003966565  0.0000000 0.06305647 0.00000000
    ## SRR9162903__vRhyme_unbinned_28  0.000000000  0.0000000 0.06275625 0.02547422
    ## SRR9161506__vRhyme_unbinned_183 0.000000000  0.0000000 0.00000000 0.00000000
    ##                                 SRR9162906  SRR9161500 SRR9162899 SRR9161509
    ## SRR9161506__vRhyme_unbinned_143  0.5335220 0.543045760 0.54446020  0.7647336
    ## SRR9161505__vRhyme_unbinned_108  0.5257441 0.000000000 0.00000000  0.4078425
    ## SRR9161504__vRhyme_unbinned_104  0.0000000 0.000000000 0.00000000  0.0000000
    ## SRR9161504__vRhyme_bin_20        0.1596670 0.004238947 0.00000000  0.0000000
    ## SRR9162903__vRhyme_unbinned_28   0.7231181 0.029847482 0.07970262  0.3933745
    ## SRR9161506__vRhyme_unbinned_183  0.0000000 0.050531300 0.00000000  0.0000000
    ##                                 SRR9161501 SRR9161507 SRR9162900  SRR9162905
    ## SRR9161506__vRhyme_unbinned_143  0.4336634 0.56878830  0.3852900 0.600754400
    ## SRR9161505__vRhyme_unbinned_108  0.8244281 0.17124437  0.1588802 0.089463920
    ## SRR9161504__vRhyme_unbinned_104  0.0000000 0.00000000  0.0000000 0.000000000
    ## SRR9161504__vRhyme_bin_20        0.0000000 0.05597453  0.0000000 0.005107165
    ## SRR9162903__vRhyme_unbinned_28   0.7812278 0.11288471  0.4466189 0.038812660
    ## SRR9161506__vRhyme_unbinned_183  0.0000000 0.00000000  0.0000000 0.672961400
    ##                                 SRR9162908 SRR9161506 SRR9162904 SRR9162907
    ## SRR9161506__vRhyme_unbinned_143 0.64912780 1.00000000  0.3325790 0.36435643
    ## SRR9161505__vRhyme_unbinned_108 0.04071359 0.24207366  0.8095028 0.12690982
    ## SRR9161504__vRhyme_unbinned_104 0.00000000 0.01438314  0.0000000 0.00000000
    ## SRR9161504__vRhyme_bin_20       0.00000000 0.22561754  0.0000000 0.00000000
    ## SRR9162903__vRhyme_unbinned_28  0.07931996 0.81692450  0.0000000 0.04302192
    ## SRR9161506__vRhyme_unbinned_183 0.00000000 0.99670830  0.0000000 0.00000000

``` r
#rownames(breadth.fw) <- fix_rownames(rownames(breadth.fw)) # Don't need this for these data
head(breadth.fw)
```

    ##                                Ga0485186  Ga0485170 Ga0485158  Ga0485173
    ## Ga0485186__vRhyme_unbinned_706 0.9932029 0.04531380 0.0000000 0.03368326
    ## Ga0485178__vRhyme_bin_183      0.4215565 0.05090661 0.4884495 0.17350566
    ## Ga0485175__vRhyme_unbinned_245 0.2289822 0.00000000 0.9950814 0.38592508
    ## Ga0485176__vRhyme_unbinned_646 0.3564419 0.00000000 0.9872688 0.01286387
    ## Ga0485158__vRhyme_bin_8        0.3061446 0.03499513 0.9635441 0.45979014
    ## Ga0485162__vRhyme_bin_87       0.3661096 0.08233111 0.8750479 0.11348867
    ##                                 Ga0485184  Ga0485175 Ga0485162 Ga0485165
    ## Ga0485186__vRhyme_unbinned_706 0.99811190 0.04425648 0.0000000 0.0000000
    ## Ga0485178__vRhyme_bin_183      0.35141742 0.67411304 0.4655677 0.5553525
    ## Ga0485175__vRhyme_unbinned_245 0.13250095 1.00000000 0.8621264 0.2454030
    ## Ga0485176__vRhyme_unbinned_646 0.08437769 0.92195475 0.9156554 0.9296134
    ## Ga0485158__vRhyme_bin_8        0.38109824 0.95226650 0.8271448 0.4755442
    ## Ga0485162__vRhyme_bin_87       0.34773594 0.79757290 0.8163678 0.8528405
    ##                                 Ga0485161  Ga0485159  Ga0485179   Ga0485168
    ## Ga0485186__vRhyme_unbinned_706 0.01465146 0.00000000 0.01208368 0.179971310
    ## Ga0485178__vRhyme_bin_183      0.55181270 0.04481773 0.72443220 0.040665545
    ## Ga0485175__vRhyme_unbinned_245 0.02270148 0.00000000 0.38168746 0.000000000
    ## Ga0485176__vRhyme_unbinned_646 0.90398514 0.00000000 0.67667925 0.006233009
    ## Ga0485158__vRhyme_bin_8        0.42664844 0.07126252 0.31886720 0.065875040
    ## Ga0485162__vRhyme_bin_87       0.93218510 0.02548738 0.52914363 0.091120355
    ##                                Ga0485157  Ga0485174  Ga0485183 Ga0485181
    ## Ga0485186__vRhyme_unbinned_706 0.0000000 0.07658032 0.08141379 0.1131334
    ## Ga0485178__vRhyme_bin_183      0.3390953 0.78666025 0.42967668 0.7331747
    ## Ga0485175__vRhyme_unbinned_245 0.0538025 0.99909190 0.01967461 0.2072645
    ## Ga0485176__vRhyme_unbinned_646 0.9894901 0.98077050 0.66361650 0.9596181
    ## Ga0485158__vRhyme_bin_8        0.3859674 0.98173280 0.13586530 0.4562561
    ## Ga0485162__vRhyme_bin_87       0.8631633 0.89691454 0.26467120 0.6443362
    ##                                 Ga0485167 Ga0485180  Ga0485178  Ga0485172
    ## Ga0485186__vRhyme_unbinned_706 0.09870856 0.1616192 0.08609622 0.13866022
    ## Ga0485178__vRhyme_bin_183      0.49442890 0.6244412 0.97283155 0.07241435
    ## Ga0485175__vRhyme_unbinned_245 0.14264093 0.1240257 0.46205070 0.00000000
    ## Ga0485176__vRhyme_unbinned_646 0.77627480 0.6714078 0.87285330 0.02237915
    ## Ga0485158__vRhyme_bin_8        0.45415136 0.1911853 0.57211070 0.05816291
    ## Ga0485162__vRhyme_bin_87       0.86778740 0.3809761 0.89094620 0.11770107
    ##                                Ga0485164  Ga0485176 Ga0485166 Ga0485182
    ## Ga0485186__vRhyme_unbinned_706 0.0000000 0.02688619 0.0226569 0.9964504
    ## Ga0485178__vRhyme_bin_183      0.4879417 0.80594254 0.3990580 0.4240707
    ## Ga0485175__vRhyme_unbinned_245 0.0000000 0.07816875 0.0000000 0.1909951
    ## Ga0485176__vRhyme_unbinned_646 0.4494397 0.89781845 0.3501426 0.1517141
    ## Ga0485158__vRhyme_bin_8        0.3627996 0.56411590 0.3311658 0.3504853
    ## Ga0485162__vRhyme_bin_87       0.8814281 0.99973970 0.8752467 0.3153762
    ##                                 Ga0485177 Ga0485185  Ga0485169  Ga0485171
    ## Ga0485186__vRhyme_unbinned_706 0.02945397 0.9978854 0.09863304 0.13125896
    ## Ga0485178__vRhyme_bin_183      0.67508390 0.4066953 0.50263370 0.08653875
    ## Ga0485175__vRhyme_unbinned_245 0.99977297 0.1082104 0.06840711 0.00000000
    ## Ga0485176__vRhyme_unbinned_646 0.99323654 0.1179630 0.65874280 0.08172535
    ## Ga0485158__vRhyme_bin_8        0.98019350 0.3800930 0.43335533 0.05564980
    ## Ga0485162__vRhyme_bin_87       0.94112104 0.3833710 0.88699883 0.15443414

``` r
# rownames(breadth.mar) <- fix_rownames(rownames(breadth.mar)) # Don't need this for these data
head(breadth.mar)
```

    ##                                 ERR594377  ERR594415 ERR594364 ERR594392
    ## ERR594392__vRhyme_unbinned_120 0.04932208 0.19825678 0.4203791 0.9930825
    ## ERR599144__vRhyme_unbinned_1   0.00000000 0.00000000 0.0000000 0.0000000
    ## ERR594388__vRhyme_unbinned_985 0.00000000 0.09034749 0.3858628 0.1199881
    ## ERR594382__vRhyme_unbinned_389 0.00000000 0.54447450 0.9505592 0.9982033
    ## ERR594379__vRhyme_unbinned_287 0.00000000 0.01835908 0.0000000 0.0000000
    ## ERR594364__vRhyme_unbinned_706 0.00000000 0.33362000 0.9977885 0.6206155
    ##                                 ERR594353  ERR599118  ERR594385  ERR594389
    ## ERR594392__vRhyme_unbinned_120 0.11434698 0.01189817 0.99398170 0.76957667
    ## ERR599144__vRhyme_unbinned_1   0.00000000 0.00000000 0.00000000 0.00000000
    ## ERR594388__vRhyme_unbinned_985 0.02340362 0.00000000 0.07906148 0.04407485
    ## ERR594382__vRhyme_unbinned_389 0.95096296 0.00000000 0.99895024 0.48950216
    ## ERR594379__vRhyme_unbinned_287 0.00000000 0.00000000 0.00000000 0.03165358
    ## ERR594364__vRhyme_unbinned_706 0.42782113 0.00000000 0.23680815 0.13861416
    ##                                 ERR599121  ERR599148  ERR594412  ERR599174
    ## ERR594392__vRhyme_unbinned_120 0.02040675 0.16401494 0.11012729 0.01383509
    ## ERR599144__vRhyme_unbinned_1   0.00000000 0.00000000 0.00000000 0.00000000
    ## ERR594388__vRhyme_unbinned_985 0.00000000 0.01538462 0.02221562 0.01793882
    ## ERR594382__vRhyme_unbinned_389 0.00000000 0.03494569 0.03421892 0.00000000
    ## ERR594379__vRhyme_unbinned_287 0.00000000 0.00000000 0.00000000 0.01405419
    ## ERR594364__vRhyme_unbinned_706 0.00000000 0.01182505 0.28469193 0.00000000
    ##                                ERR594404  ERR594388  ERR594411   ERR599065
    ## ERR594392__vRhyme_unbinned_120 0.0000000 0.09518539 0.01950747 0.006917543
    ## ERR599144__vRhyme_unbinned_1   0.1447247 0.00000000 0.00000000 0.000000000
    ## ERR594388__vRhyme_unbinned_985 0.0000000 0.99833680 0.01621622 0.012295812
    ## ERR594382__vRhyme_unbinned_389 0.0000000 0.97508780 0.30677918 0.000000000
    ## ERR594379__vRhyme_unbinned_287 0.0000000 0.41884020 0.00000000 0.000000000
    ## ERR594364__vRhyme_unbinned_706 0.0000000 0.99926287 0.12500767 0.000000000
    ##                                  ERR599044 ERR594379 ERR594354  ERR594359
    ## ERR594392__vRhyme_unbinned_120 0.006848367 0.2558799 0.2909519 0.50498060
    ## ERR599144__vRhyme_unbinned_1   0.000000000 0.0000000 0.0000000 0.00000000
    ## ERR594388__vRhyme_unbinned_985 0.000000000 1.0000000 0.0000000 0.02263142
    ## ERR594382__vRhyme_unbinned_389 0.000000000 0.8153592 0.7580248 0.98100290
    ## ERR594379__vRhyme_unbinned_287 0.000000000 1.0000000 0.0000000 0.00000000
    ## ERR594364__vRhyme_unbinned_706 0.000000000 0.0000000 0.2260274 0.63403773
    ##                                  ERR599126  ERR599146  ERR599122   ERR599150
    ## ERR594392__vRhyme_unbinned_120 0.200885440 0.12921970 0.02296624 0.453514100
    ## ERR599144__vRhyme_unbinned_1   0.004188916 0.00000000 0.00000000 0.000000000
    ## ERR594388__vRhyme_unbinned_985 0.000000000 0.00000000 0.00000000 0.012177012
    ## ERR594382__vRhyme_unbinned_389 0.080106590 0.01487867 0.00000000 0.007267735
    ## ERR594379__vRhyme_unbinned_287 0.006330716 0.00000000 0.00000000 0.000000000
    ## ERR594364__vRhyme_unbinned_706 0.017998649 0.00577431 0.00000000 0.000000000
    ##                                  ERR599173  ERR598982   ERR599006  ERR594362
    ## ERR594392__vRhyme_unbinned_120 0.004150526 0.01307416 0.163738240 0.38620642
    ## ERR599144__vRhyme_unbinned_1   0.000000000 0.00000000 0.000000000 0.00000000
    ## ERR594388__vRhyme_unbinned_985 0.016097415 0.00000000 0.016038015 0.02946243
    ## ERR594382__vRhyme_unbinned_389 0.000000000 0.00000000 0.004946098 0.42964430
    ## ERR594379__vRhyme_unbinned_287 0.000000000 0.00000000 0.010888833 0.00000000
    ## ERR594364__vRhyme_unbinned_706 0.000000000 0.00000000 0.018827938 0.13879845
    ##                                ERR599144  ERR599017  ERR594391   ERR599005
    ## ERR594392__vRhyme_unbinned_120 0.0000000 0.03195905 0.25698670 0.089305475
    ## ERR599144__vRhyme_unbinned_1   0.9991944 0.00000000 0.00000000 0.000000000
    ## ERR594388__vRhyme_unbinned_985 0.0000000 0.01823582 0.03932284 0.066943870
    ## ERR594382__vRhyme_unbinned_389 0.0000000 0.00000000 0.38444300 0.005410425
    ## ERR594379__vRhyme_unbinned_287 0.0000000 0.01107875 0.00000000 0.011205369
    ## ERR594364__vRhyme_unbinned_706 0.0000000 0.00000000 0.12448553 0.146108480
    ##                                 ERR594382  ERR598984  ERR594355   ERR599165
    ## ERR594392__vRhyme_unbinned_120 0.29219702 0.55160487 0.29371887 0.000000000
    ## ERR599144__vRhyme_unbinned_1   0.01150801 0.00000000 0.00000000 0.000000000
    ## ERR594388__vRhyme_unbinned_985 0.02381942 0.19501040 0.01829522 0.009088209
    ## ERR594382__vRhyme_unbinned_389 0.98298140 0.52081400 0.74247990 0.000000000
    ## ERR594379__vRhyme_unbinned_287 0.00000000 0.02158774 0.00000000 0.000000000
    ## ERR594364__vRhyme_unbinned_706 0.58707535 0.96756560 0.36762086 0.000000000
    ##                                  ERR599133   ERR594380  ERR594409   ERR599110
    ## ERR594392__vRhyme_unbinned_120 0.006917543 0.000000000 0.00000000 0.109504710
    ## ERR599144__vRhyme_unbinned_1   0.000000000 0.000000000 0.00000000 0.000000000
    ## ERR594388__vRhyme_unbinned_985 0.983962000 0.008375408 0.00000000 0.005940006
    ## ERR594382__vRhyme_unbinned_389 0.000000000 0.001110349 0.00000000 0.002987847
    ## ERR594379__vRhyme_unbinned_287 0.000000000 0.000000000 0.00000000 0.000000000
    ## ERR594364__vRhyme_unbinned_706 0.000000000 0.000000000 0.02150009 0.000000000
    ##                                  ERR599023   ERR594407
    ## ERR594392__vRhyme_unbinned_120 0.267847270 0.000000000
    ## ERR599144__vRhyme_unbinned_1   0.000000000 0.000000000
    ## ERR594388__vRhyme_unbinned_985 0.005999406 0.000000000
    ## ERR594382__vRhyme_unbinned_389 0.018330844 0.006117011
    ## ERR594379__vRhyme_unbinned_287 0.000000000 0.000000000
    ## ERR594364__vRhyme_unbinned_706 0.005866454 0.000000000

``` r
rownames(breadth.soil) <- fix_rownames(rownames(breadth.soil))
head(breadth.soil)
```

    ##                                 SRR8487029 SRR8487034 SRR8487015  SRR8487030
    ## SRR8487014__vRhyme_bin_28       0.00000000  0.6293566  0.5848656 0.008949829
    ## SRR8487018__vRhyme_bin_107      0.02230951  0.9014331  0.2205788 0.003805199
    ## SRR8487020__vRhyme_unbinned_111 0.04263318  0.9253154  0.6870597 0.100999690
    ## SRR8487022__vRhyme_unbinned_642 0.00000000  0.9144880  0.6488250 0.006226809
    ## SRR8487022__vRhyme_unbinned_703 0.01787150  0.6307151  0.1990886 0.000000000
    ## SRR8487017__vRhyme_unbinned_41  0.26196885  0.9989906  0.9983177 0.401076700
    ##                                 SRR8487021 SRR8487010 SRR8487022 SRR8487020
    ## SRR8487014__vRhyme_bin_28        0.8992857  0.7432661  0.5294742  0.5161355
    ## SRR8487018__vRhyme_bin_107       0.8515666  0.7367355  0.9394544  0.8258201
    ## SRR8487020__vRhyme_unbinned_111  0.9350617  0.7833812  0.9898360  0.9959065
    ## SRR8487022__vRhyme_unbinned_642  0.9976969  0.6637310  0.9830682  0.9347251
    ## SRR8487022__vRhyme_unbinned_703  0.6500760  0.8624192  0.9262205  0.9208888
    ## SRR8487017__vRhyme_unbinned_41   0.9997116  0.9977889  0.9997116  0.9990386
    ##                                 SRR8487018 SRR8487035 SRR8487011  SRR8487013
    ## SRR8487014__vRhyme_bin_28        0.6936118  0.9288030  0.4099423 0.019190500
    ## SRR8487018__vRhyme_bin_107       0.9830300  0.8294105  0.8583484 0.021787830
    ## SRR8487020__vRhyme_unbinned_111  0.9474256  0.8705968  0.9962128 0.042020550
    ## SRR8487022__vRhyme_unbinned_642  0.9704440  0.9994029  0.4209281 0.001492728
    ## SRR8487022__vRhyme_unbinned_703  0.7848569  0.5837727  0.7616835 0.000000000
    ## SRR8487017__vRhyme_unbinned_41   0.9995194  0.9990386  0.9990386 0.478657960
    ##                                 SRR8487023  SRR8487039 SRR8487032 SRR8487025
    ## SRR8487014__vRhyme_bin_28        0.8991710 0.002552996 0.00000000 0.00000000
    ## SRR8487018__vRhyme_bin_107       0.9148127 0.012888575 0.02470310 0.01743027
    ## SRR8487020__vRhyme_unbinned_111  0.9762747 0.038985267 0.05095931 0.05104285
    ## SRR8487022__vRhyme_unbinned_642  0.9103937 0.003305327 0.06356890 0.00000000
    ## SRR8487022__vRhyme_unbinned_703  0.7678491 0.001578650 0.00000000 0.00000000
    ## SRR8487017__vRhyme_unbinned_41   1.0000000 0.760382600 0.26874640 0.32070756
    ##                                 SRR8487014 SRR8487036 SRR8487019 SRR8487017
    ## SRR8487014__vRhyme_bin_28        0.9052809  0.5758584  0.7297266  0.7674766
    ## SRR8487018__vRhyme_bin_107       0.8409795  0.7645380  0.8890662  0.7139043
    ## SRR8487020__vRhyme_unbinned_111  0.9233382  0.7467350  0.9877475  0.9973546
    ## SRR8487022__vRhyme_unbinned_642  0.6702777  0.9971212  0.5899049  0.7322472
    ## SRR8487022__vRhyme_unbinned_703  0.5533018  0.4889941  0.8789206  0.9236589
    ## SRR8487017__vRhyme_unbinned_41   0.9991348  0.9991348  0.9986541  0.9990867
    ##                                 SRR8487037  SRR8487024 SRR8487031 SRR8487027
    ## SRR8487014__vRhyme_bin_28       0.00000000 0.006196036 0.01052752 0.02349330
    ## SRR8487018__vRhyme_bin_107      0.02138890 0.000000000 0.00000000 0.01586522
    ## SRR8487020__vRhyme_unbinned_111 0.01559411 0.004176993 0.04163070 0.01356130
    ## SRR8487022__vRhyme_unbinned_642 0.00000000 0.000000000 0.00000000 0.01571630
    ## SRR8487022__vRhyme_unbinned_703 0.01251005 0.004467876 0.00000000 0.00000000
    ## SRR8487017__vRhyme_unbinned_41  0.27014035 0.077437030 0.38136896 0.50230724
    ##                                  SRR8487026  SRR8487028  SRR8487040  SRR8487038
    ## SRR8487014__vRhyme_bin_28       0.022518000 0.000000000 0.019649465 0.000000000
    ## SRR8487018__vRhyme_bin_107      0.015343542 0.000000000 0.009022003 0.006965968
    ## SRR8487020__vRhyme_unbinned_111 0.024282254 0.115090090 0.025841665 0.097741640
    ## SRR8487022__vRhyme_unbinned_642 0.045272317 0.005693692 0.035398986 0.000000000
    ## SRR8487022__vRhyme_unbinned_703 0.004467876 0.017484289 0.003723230 0.006850743
    ## SRR8487017__vRhyme_unbinned_41  0.395212470 0.239665450 0.444337640 0.203134020
    ##                                  SRR8487012 SRR8487016
    ## SRR8487014__vRhyme_bin_28       0.011215972  0.5707237
    ## SRR8487018__vRhyme_bin_107      0.020897904  0.8161537
    ## SRR8487020__vRhyme_unbinned_111 0.052518725  0.9998608
    ## SRR8487022__vRhyme_unbinned_642 0.023670405  0.3914787
    ## SRR8487022__vRhyme_unbinned_703 0.002442439  0.8374885
    ## SRR8487017__vRhyme_unbinned_41  0.562391900  1.0000000

### Split into viromes and metagenomes

``` r
# Human gut
breadth.gut.vir <- breadth.gut[colnames(breadth.gut) %in% subset(data_reform, env2=="Human gut" & method2 == "Virome")$Sample]
breadth.gut.mg <- breadth.gut[colnames(breadth.gut) %in% subset(data_reform, env2=="Human gut" & method2 == "Mixed MG")$Sample]

# Freshwater
breadth.fw.vir <- breadth.fw[colnames(breadth.fw) %in% subset(data_reform, env2=="Freshwater" & method2 == "Virome")$Sample]
breadth.fw.mg <- breadth.fw[colnames(breadth.fw) %in% subset(data_reform, env2=="Freshwater" & method2 == "Mixed MG")$Sample]

# Marine
breadth.mar.vir <- breadth.mar[colnames(breadth.mar) %in% subset(data_reform, env2=="Marine" & method2 == "Virome")$Sample]
breadth.mar.mg <- breadth.mar[colnames(breadth.mar) %in% subset(data_reform, env2=="Marine" & method2 == "Mixed MG")$Sample]

# Soil
breadth.soil.vir <- breadth.soil[colnames(breadth.soil) %in% subset(data_reform, env2=="Soil" & method2 == "Virome")$Sample]
breadth.soil.mg <- breadth.soil[colnames(breadth.soil) %in% subset(data_reform, env2=="Soil" & method2 == "Mixed MG")$Sample]
```

### Get names of present VMAGs in viromes and metagenomes, and both

``` r
# Consider values >= 0.75 to be "present"
## Each environment, separate
present.gut <- rownames(breadth.gut)[apply(breadth.gut, 1, function(row) any(row >= 0.75))]
saveRDS(present.gut, file="../Data/present_gut.RDS")
present.gut.vir <- rownames(breadth.gut.vir)[apply(breadth.gut.vir, 1, function(row) any(row >= 0.75))]
saveRDS(present.gut.vir, file="../Data/present_gut_vir.RDS")
present.gut.mg <- rownames(breadth.gut.mg)[apply(breadth.gut.mg, 1, function(row) any(row >= 0.75))]
saveRDS(present.gut.mg, file="../Data/present_gut_mg.RDS")
present.fw <- rownames(breadth.fw)[apply(breadth.fw, 1, function(row) any(row >= 0.75))]
saveRDS(present.fw, file="../Data/present_fw.RDS")
present.fw.vir <- rownames(breadth.fw.vir)[apply(breadth.fw.vir, 1, function(row) any(row >= 0.75))]
saveRDS(present.fw.vir, file="../Data/present_fw_vir.RDS")
present.fw.mg <- rownames(breadth.fw.mg)[apply(breadth.fw.mg, 1, function(row) any(row >= 0.75))]
saveRDS(present.fw.mg, file="../Data/present_fw_mg.RDS")
present.mar <- rownames(breadth.mar)[apply(breadth.mar, 1, function(row) any(row >= 0.75))]
saveRDS(present.mar, file="../Data/present_mar.RDS")
present.mar.vir <- rownames(breadth.mar.vir)[apply(breadth.mar.vir, 1, function(row) any(row >= 0.75))]
saveRDS(present.mar.vir, file="../Data/present_mar_vir.RDS")
present.mar.mg <- rownames(breadth.mar.mg)[apply(breadth.mar.mg, 1, function(row) any(row >= 0.75))]
saveRDS(present.mar.mg, file="../Data/present_mar_mg.RDS")
present.soil <- rownames(breadth.soil)[apply(breadth.soil, 1, function(row) any(row >= 0.75))]
saveRDS(present.soil, file="../Data/present_soil.RDS")
present.soil.vir <- rownames(breadth.soil.vir)[apply(breadth.soil.vir, 1, function(row) any(row >= 0.75))]
saveRDS(present.soil.vir, file="../Data/present_soil_vir.RDS")
present.soil.mg <- rownames(breadth.soil.mg)[apply(breadth.soil.mg, 1, function(row) any(row >= 0.75))]
saveRDS(present.soil.mg, file="../Data/present_soil_mg.RDS")

## All environments, combined
present.all <- c(present.gut, present.fw, present.mar, present.soil)
saveRDS(present.all, file="../Data/present_all.RDS")
present.all.vir <- c(present.gut.vir, present.fw.vir, present.mar.vir, present.soil.vir)
saveRDS(present.all.vir, file="../Data/present_all_vir.RDS")
present.all.mg <- c(present.gut.mg, present.fw.mg, present.mar.mg, present.soil.mg)
saveRDS(present.all.mg, file="../Data/present_all_mg.RDS")
```

## CoverM trimmed means

``` r
tmeans.gut <- read.csv("../Tables/trimmed_means_human_gut.csv", header=TRUE, row.names = 1)
head(tmeans.gut)
```

    ##                                SRR9162899 SRR9162906 SRR9161507 SRR9162908
    ## SRR9161506_vRhyme_unbinned_143   1.753055  1.2873234  2.3557540   2.912124
    ## SRR9161505_vRhyme_unbinned_108   0.000000  1.0535945  0.1804713   0.000000
    ## SRR9161504_vRhyme_unbinned_104   0.000000  0.0000000  0.0000000   0.000000
    ## SRR9161504_vRhyme_bin_20         0.000000  0.2230545  0.0000000   0.000000
    ## SRR9162903_vRhyme_unbinned_28    0.000000  2.5800722  0.1123767   0.000000
    ## SRR9161506_vRhyme_unbinned_183   0.000000  0.0000000  0.0000000   0.000000
    ##                                SRR9161510 SRR9161509 SRR9162905 SRR9161500
    ## SRR9161506_vRhyme_unbinned_143          0 18.9384770   1.354904   2.307619
    ## SRR9161505_vRhyme_unbinned_108          0  0.8289748   0.000000   0.000000
    ## SRR9161504_vRhyme_unbinned_104          0  0.0000000   0.000000   0.000000
    ## SRR9161504_vRhyme_bin_20                0  0.0000000   0.000000   0.000000
    ## SRR9162903_vRhyme_unbinned_28           0  0.9560904   0.000000   0.000000
    ## SRR9161506_vRhyme_unbinned_183          0  0.0000000   1.912880   0.000000
    ##                                SRR9161503 SRR9162907 SRR9161502 SRR9162900
    ## SRR9161506_vRhyme_unbinned_143   5.584635  0.7404101  8.2197430  0.7956647
    ## SRR9161505_vRhyme_unbinned_108   6.362534  0.1282689  0.5984886  0.1677439
    ## SRR9161504_vRhyme_unbinned_104   0.000000  0.0000000  0.0000000  0.0000000
    ## SRR9161504_vRhyme_bin_20         0.000000  0.0000000  2.7557466  0.0000000
    ## SRR9162903_vRhyme_unbinned_28    0.000000  0.0000000  2.7896993  0.9798518
    ## SRR9161506_vRhyme_unbinned_183   0.000000  0.0000000  0.0000000  0.0000000
    ##                                SRR9161506 SRR9162904 SRR9161505 SRR9161501
    ## SRR9161506_vRhyme_unbinned_143 42.3580930  0.9413452   6.149293   2.250133
    ## SRR9161505_vRhyme_unbinned_108  0.3843094  2.5317688   7.708661   4.257731
    ## SRR9161504_vRhyme_unbinned_104  0.0000000  0.0000000   0.000000   0.000000
    ## SRR9161504_vRhyme_bin_20        0.4139400  0.0000000   0.000000   0.000000
    ## SRR9162903_vRhyme_unbinned_28   3.5868700  0.0000000   0.000000   4.656011
    ## SRR9161506_vRhyme_unbinned_183 20.3844660  0.0000000   0.000000   0.000000
    ##                                 SRR9161504 SRR9162901 SRR9162903 SRR9162902
    ## SRR9161506_vRhyme_unbinned_143   1.8155351  0.5911168   0.467538  0.6798427
    ## SRR9161505_vRhyme_unbinned_108   0.9881675  0.0000000   6.030327  4.3285275
    ## SRR9161504_vRhyme_unbinned_104 110.9030100  0.0000000   0.239466  0.0000000
    ## SRR9161504_vRhyme_bin_20        32.2554500  0.0000000   4.287088  0.0000000
    ## SRR9162903_vRhyme_unbinned_28    6.4746770  0.0000000   8.078755  0.0000000
    ## SRR9161506_vRhyme_unbinned_183   0.0000000  0.0000000   0.000000  0.0000000

``` r
tmeans.fw <- read.csv("../Tables/trimmed_means_freshwater.csv", header=TRUE, row.names = 1)
head(tmeans.fw)
```

    ##                                 Ga0485183 Ga0485176 Ga0485166 Ga0485157
    ## Ga0485186__vRhyme_unbinned_706 0.00000000  0.000000 0.0000000 0.0000000
    ## Ga0485178__vRhyme_bin_183      0.83286107 13.772288 0.5828834 0.5044175
    ## Ga0485175__vRhyme_unbinned_245 0.00000000  0.000000 0.0000000 0.0000000
    ## Ga0485176__vRhyme_unbinned_646 1.10102920  3.034799 0.4081889 7.8085666
    ## Ga0485158__vRhyme_bin_8        0.09614238  1.535473 0.4877258 0.5339996
    ## Ga0485162__vRhyme_bin_87       0.25769553 13.565138 2.4351100 2.4662397
    ##                                Ga0485164 Ga0485170 Ga0485162 Ga0485179
    ## Ga0485186__vRhyme_unbinned_706 0.0000000         0  0.000000 0.0000000
    ## Ga0485178__vRhyme_bin_183      0.9032562         0  1.160608 3.3897026
    ## Ga0485175__vRhyme_unbinned_245 0.0000000         0  2.469683 0.6691045
    ## Ga0485176__vRhyme_unbinned_646 0.6241670         0  2.971161 1.4014142
    ## Ga0485158__vRhyme_bin_8        0.7221287         0  2.233596 0.4374364
    ## Ga0485162__vRhyme_bin_87       2.8559217         0  2.038429 0.8024784
    ##                                Ga0485165  Ga0485184  Ga0485173  Ga0485178
    ## Ga0485186__vRhyme_unbinned_706 0.0000000 12.8446060 0.00000000  0.0000000
    ## Ga0485178__vRhyme_bin_183      1.6622306  1.0470695 0.30552518 14.0164200
    ## Ga0485175__vRhyme_unbinned_245 0.2142189  0.1398078 0.71332600  0.9806106
    ## Ga0485176__vRhyme_unbinned_646 3.0280247  0.0000000 0.00000000  2.5854805
    ## Ga0485158__vRhyme_bin_8        0.8113800  0.6437138 0.81672806  1.3252324
    ## Ga0485162__vRhyme_bin_87       2.1584240  0.4144973 0.07072777  2.4923348
    ##                                Ga0485177 Ga0485167 Ga0485174  Ga0485182
    ## Ga0485186__vRhyme_unbinned_706  0.000000 0.0000000  0.000000 20.8128660
    ## Ga0485178__vRhyme_bin_183       4.601953 1.1345545  6.096799  0.9376590
    ## Ga0485175__vRhyme_unbinned_245 15.876520 0.1048559 16.089634  0.2605664
    ## Ga0485176__vRhyme_unbinned_646  6.587110 1.6229824  4.412113  0.1139123
    ## Ga0485158__vRhyme_bin_8         9.243135 0.7715413  9.832158  0.4987901
    ## Ga0485162__vRhyme_bin_87        3.504211 2.1667080  2.851214  0.3605852
    ##                                 Ga0485171 Ga0485159   Ga0485185 Ga0485161
    ## Ga0485186__vRhyme_unbinned_706 0.09208181         0 19.80879200 0.0000000
    ## Ga0485178__vRhyme_bin_183      0.00000000         0  1.54297000 1.8348746
    ## Ga0485175__vRhyme_unbinned_245 0.00000000         0  0.09915809 0.0000000
    ## Ga0485176__vRhyme_unbinned_646 0.00000000         0  0.07622538 2.7981640
    ## Ga0485158__vRhyme_bin_8        0.00000000         0  0.60317380 0.7748378
    ## Ga0485162__vRhyme_bin_87       0.11601354         0  0.47615638 3.7359164
    ##                                Ga0485175  Ga0485181  Ga0485172 Ga0485186
    ## Ga0485186__vRhyme_unbinned_706  0.000000 0.07171348 0.11304422 9.6803020
    ## Ga0485178__vRhyme_bin_183       3.989475 3.51506520 0.00000000 1.3080574
    ## Ga0485175__vRhyme_unbinned_245 38.130795 0.34042010 0.00000000 0.2842929
    ## Ga0485176__vRhyme_unbinned_646  3.045609 3.85132530 0.00000000 0.4479120
    ## Ga0485158__vRhyme_bin_8        16.333230 1.03591100 0.00000000 0.4589690
    ## Ga0485162__vRhyme_bin_87        1.980477 1.27731810 0.07541486 0.4612472
    ##                                Ga0485158 Ga0485169 Ga0485168 Ga0485180
    ## Ga0485186__vRhyme_unbinned_706  0.000000 0.0000000 0.1468217 0.1261987
    ## Ga0485178__vRhyme_bin_183       1.427775 1.0919192 0.0000000 1.7555721
    ## Ga0485175__vRhyme_unbinned_245 12.599200 0.0000000 0.0000000 0.1110639
    ## Ga0485176__vRhyme_unbinned_646  6.572709 1.0436473 0.0000000 1.1999851
    ## Ga0485158__vRhyme_bin_8         5.059863 0.8346835 0.0000000 0.1838506
    ## Ga0485162__vRhyme_bin_87        2.317922 2.5905268 0.0000000 0.4682620

``` r
tmeans.mar <- read.csv("../Tables/trimmed_means_marine.csv", header=TRUE, row.names = 1)
head(tmeans.mar)
```

    ##                                 ERR594391 ERR594354  ERR594389 ERR599121
    ## ERR594392__vRhyme_unbinned_120 0.36556384 0.3986486 1.96707050         0
    ## ERR599144__vRhyme_unbinned_1   0.00000000 0.0000000 0.00000000         0
    ## ERR594388__vRhyme_unbinned_985 0.00000000 0.0000000 0.00000000         0
    ## ERR594382__vRhyme_unbinned_389 0.55838543 2.5649102 0.71830990         0
    ## ERR594379__vRhyme_unbinned_287 0.00000000 0.0000000 0.00000000         0
    ## ERR594364__vRhyme_unbinned_706 0.08344761 0.2316923 0.09921832         0
    ##                                ERR599118 ERR594392 ERR594380 ERR594404
    ## ERR594392__vRhyme_unbinned_120         0 6.3503420         0 0.0000000
    ## ERR599144__vRhyme_unbinned_1           0 0.0000000         0 0.1207339
    ## ERR594388__vRhyme_unbinned_985         0 0.1035493         0 0.0000000
    ## ERR594382__vRhyme_unbinned_389         0 7.2180624         0 0.0000000
    ## ERR594379__vRhyme_unbinned_287         0 0.0000000         0 0.0000000
    ## ERR594364__vRhyme_unbinned_706         0 1.9586533         0 0.0000000
    ##                                ERR599023 ERR594415 ERR599122  ERR594364
    ## ERR594392__vRhyme_unbinned_120 0.3543803 0.1770736         0  0.5949829
    ## ERR599144__vRhyme_unbinned_1   0.0000000 0.0000000         0  0.0000000
    ## ERR594388__vRhyme_unbinned_985 0.0000000 0.0000000         0  0.5746154
    ## ERR594382__vRhyme_unbinned_389 0.0000000 1.1001664         0  6.8503350
    ## ERR594379__vRhyme_unbinned_287 0.0000000 0.0000000         0  0.0000000
    ## ERR594364__vRhyme_unbinned_706 0.0000000 0.4734983         0 22.1227720
    ##                                ERR594355   ERR594379  ERR599110 ERR599133
    ## ERR594392__vRhyme_unbinned_120 0.3689811   0.3257223 0.06748991  0.000000
    ## ERR599144__vRhyme_unbinned_1   0.0000000   0.0000000 0.00000000  0.000000
    ## ERR594388__vRhyme_unbinned_985 0.0000000 109.9719000 0.00000000  5.610908
    ## ERR594382__vRhyme_unbinned_389 2.1656168   3.9211178 0.00000000  0.000000
    ## ERR594379__vRhyme_unbinned_287 0.0000000  35.0715100 0.00000000  0.000000
    ## ERR594364__vRhyme_unbinned_706 0.5737795   0.0000000 0.00000000  0.000000
    ##                                 ERR599146 ERR599017 ERR599126 ERR594407
    ## ERR594392__vRhyme_unbinned_120 0.09948742         0 0.3348089         0
    ## ERR599144__vRhyme_unbinned_1   0.00000000         0 0.0000000         0
    ## ERR594388__vRhyme_unbinned_985 0.00000000         0 0.0000000         0
    ## ERR594382__vRhyme_unbinned_389 0.00000000         0 0.0000000         0
    ## ERR594379__vRhyme_unbinned_287 0.00000000         0 0.0000000         0
    ## ERR594364__vRhyme_unbinned_706 0.00000000         0 0.0000000         0
    ##                                ERR599005 ERR599165  ERR594353 ERR598984
    ## ERR594392__vRhyme_unbinned_120 0.0000000         0 0.07292637 1.1387854
    ## ERR599144__vRhyme_unbinned_1   0.0000000         0 0.00000000 0.0000000
    ## ERR594388__vRhyme_unbinned_985 0.0000000         0 0.00000000 0.2442565
    ## ERR594382__vRhyme_unbinned_389 0.0000000         0 6.29896970 0.9493768
    ## ERR594379__vRhyme_unbinned_287 0.0000000         0 0.00000000 0.0000000
    ## ERR594364__vRhyme_unbinned_706 0.1142691         0 1.20179650 5.5198507
    ##                                ERR594382 ERR594377  ERR594412  ERR594362
    ## ERR594392__vRhyme_unbinned_120 0.2927928         0 0.06818888 0.49029200
    ## ERR599144__vRhyme_unbinned_1   0.0000000         0 0.00000000 0.00000000
    ## ERR594388__vRhyme_unbinned_985 0.0000000         0 0.00000000 0.00000000
    ## ERR594382__vRhyme_unbinned_389 4.6751113         0 0.00000000 0.55064570
    ## ERR594379__vRhyme_unbinned_287 0.0000000         0 0.00000000 0.00000000
    ## ERR594364__vRhyme_unbinned_706 1.2986492         0 0.40530717 0.09942403
    ##                                ERR594359 ERR599148 ERR594409 ERR599173
    ## ERR594392__vRhyme_unbinned_120 0.7819975 0.1745884         0         0
    ## ERR599144__vRhyme_unbinned_1   0.0000000 0.0000000         0         0
    ## ERR594388__vRhyme_unbinned_985 0.0000000 0.0000000         0         0
    ## ERR594382__vRhyme_unbinned_389 4.3719570 0.0000000         0         0
    ## ERR594379__vRhyme_unbinned_287 0.0000000 0.0000000         0         0
    ## ERR594364__vRhyme_unbinned_706 1.6465304 0.0000000         0         0
    ##                                 ERR594388 ERR599150 ERR598982 ERR599144
    ## ERR594392__vRhyme_unbinned_120  0.0000000 0.7895309         0   0.00000
    ## ERR599144__vRhyme_unbinned_1    0.0000000 0.0000000         0  18.48268
    ## ERR594388__vRhyme_unbinned_985 13.8304590 0.0000000         0   0.00000
    ## ERR594382__vRhyme_unbinned_389 17.2626110 0.0000000         0   0.00000
    ## ERR594379__vRhyme_unbinned_287  0.8211191 0.0000000         0   0.00000
    ## ERR594364__vRhyme_unbinned_706 18.5302730 0.0000000         0   0.00000
    ##                                 ERR594385 ERR599006 ERR599044 ERR599174
    ## ERR594392__vRhyme_unbinned_120 12.3780680 0.2492234         0         0
    ## ERR599144__vRhyme_unbinned_1    0.0000000 0.0000000         0         0
    ## ERR594388__vRhyme_unbinned_985  0.0000000 0.0000000         0         0
    ## ERR594382__vRhyme_unbinned_389 16.8601900 0.0000000         0         0
    ## ERR594379__vRhyme_unbinned_287  0.0000000 0.0000000         0         0
    ## ERR594364__vRhyme_unbinned_706  0.2342979 0.0000000         0         0
    ##                                ERR599065  ERR594411
    ## ERR594392__vRhyme_unbinned_120         0 0.00000000
    ## ERR599144__vRhyme_unbinned_1           0 0.00000000
    ## ERR594388__vRhyme_unbinned_985         0 0.00000000
    ## ERR594382__vRhyme_unbinned_389         0 0.32572109
    ## ERR594379__vRhyme_unbinned_287         0 0.00000000
    ## ERR594364__vRhyme_unbinned_706         0 0.08403044

``` r
tmeans.soil <- read.csv("../Tables/trimmed_means_soil.csv", header=TRUE, row.names = 1)
head(tmeans.soil)
```

    ##                                SRR8487021 SRR8487020 SRR8487015 SRR8487022
    ## SRR8487014_vRhyme_bin_28         3.906170   1.326919  1.1343838   1.943366
    ## SRR8487018_vRhyme_bin_107        3.355174   3.568017  0.2265391   4.452218
    ## SRR8487020_vRhyme_unbinned_111   5.398291   7.504707  1.6354202   6.316918
    ## SRR8487022_vRhyme_unbinned_642   6.961636   2.755764  1.1474684   5.667198
    ## SRR8487022_vRhyme_unbinned_703   1.261062   3.482630  0.2001263   3.596456
    ## SRR8487017_vRhyme_unbinned_41  125.396880  81.483920 23.6239380 111.131360
    ##                                SRR8487040 SRR8487010 SRR8487035 SRR8487024
    ## SRR8487014_vRhyme_bin_28        0.0000000   3.083626  4.8675537          0
    ## SRR8487018_vRhyme_bin_107       0.0000000   1.737534  2.7461715          0
    ## SRR8487020_vRhyme_unbinned_111  0.0000000   3.775703  3.7395060          0
    ## SRR8487022_vRhyme_unbinned_642  0.0000000   1.186641 10.9988360          0
    ## SRR8487022_vRhyme_unbinned_703  0.0000000   2.156611  0.9558193          0
    ## SRR8487017_vRhyme_unbinned_41   0.7174287  87.216354 99.8252260          0
    ##                                SRR8487013 SRR8487037 SRR8487032 SRR8487025
    ## SRR8487014_vRhyme_bin_28        0.0000000  0.0000000  0.0000000  0.0000000
    ## SRR8487018_vRhyme_bin_107       0.0000000  0.0000000  0.0000000  0.0000000
    ## SRR8487020_vRhyme_unbinned_111  0.0000000  0.0000000  0.0000000  0.0000000
    ## SRR8487022_vRhyme_unbinned_642  0.0000000  0.0000000  0.0000000  0.0000000
    ## SRR8487022_vRhyme_unbinned_703  0.0000000  0.0000000  0.0000000  0.0000000
    ## SRR8487017_vRhyme_unbinned_41   0.7308768  0.2745562  0.2838623  0.3750404
    ##                                 SRR8487034 SRR8487027 SRR8487017 SRR8487029
    ## SRR8487014_vRhyme_bin_28         1.8451132  0.0000000   3.898163  0.0000000
    ## SRR8487018_vRhyme_bin_107        3.0258784  0.0000000   1.750370  0.0000000
    ## SRR8487020_vRhyme_unbinned_111   3.0839212  0.0000000   9.986702  0.0000000
    ## SRR8487022_vRhyme_unbinned_642   2.6765392  0.0000000   1.444378  0.0000000
    ## SRR8487022_vRhyme_unbinned_703   0.9800538  0.0000000   2.973771  0.0000000
    ## SRR8487017_vRhyme_unbinned_41  103.2752000  0.8369015  83.088540  0.2916622
    ##                                SRR8487014 SRR8487028 SRR8487023 SRR8487012
    ## SRR8487014_vRhyme_bin_28        3.4940104  0.0000000   3.326531  0.0000000
    ## SRR8487018_vRhyme_bin_107       2.4005299  0.0000000   3.379848  0.0000000
    ## SRR8487020_vRhyme_unbinned_111  3.1733105  0.0728911   6.802610  0.0000000
    ## SRR8487022_vRhyme_unbinned_642  1.3404564  0.0000000   3.103257  0.0000000
    ## SRR8487022_vRhyme_unbinned_703  0.8799907  0.0000000   1.750573  0.0000000
    ## SRR8487017_vRhyme_unbinned_41  99.7961300  0.2162453 123.793920  0.9870361
    ##                                SRR8487011 SRR8487019 SRR8487030  SRR8487016
    ## SRR8487014_vRhyme_bin_28        1.1512060   2.322108 0.00000000   2.7510898
    ## SRR8487018_vRhyme_bin_107       3.2740630   3.324030 0.00000000   4.2726180
    ## SRR8487020_vRhyme_unbinned_111  9.1455340   6.774087 0.05831909  18.9131580
    ## SRR8487022_vRhyme_unbinned_642  0.5676967   1.013169 0.00000000   0.5150939
    ## SRR8487022_vRhyme_unbinned_703  1.7656660   2.770353 0.00000000   2.2078722
    ## SRR8487017_vRhyme_unbinned_41  79.4266300  79.711890 0.56304467 116.3724600
    ##                                SRR8487036 SRR8487026 SRR8487031 SRR8487038
    ## SRR8487014_vRhyme_bin_28        1.3942075  0.0000000  0.0000000  0.0000000
    ## SRR8487018_vRhyme_bin_107       2.5837780  0.0000000  0.0000000  0.0000000
    ## SRR8487020_vRhyme_unbinned_111  1.7520274  0.0000000  0.0000000  0.0000000
    ## SRR8487022_vRhyme_unbinned_642  6.6055620  0.0000000  0.0000000  0.0000000
    ## SRR8487022_vRhyme_unbinned_703  0.6959875  0.0000000  0.0000000  0.0000000
    ## SRR8487017_vRhyme_unbinned_41  83.9803100  0.5764927  0.5683163  0.1812803
    ##                                SRR8487039 SRR8487018
    ## SRR8487014_vRhyme_bin_28         0.000000   4.081496
    ## SRR8487018_vRhyme_bin_107        0.000000   5.669018
    ## SRR8487020_vRhyme_unbinned_111   0.000000   3.933012
    ## SRR8487022_vRhyme_unbinned_642   0.000000   4.322201
    ## SRR8487022_vRhyme_unbinned_703   0.000000   2.200259
    ## SRR8487017_vRhyme_unbinned_41    1.908553 159.976820

## Subset trimmed means tables for each method

``` r
tmeans.gut.vir <- tmeans.gut[colnames(tmeans.gut) %in% subset(data_reform, env2=="Human gut" & method2 == "Virome")$Sample]
saveRDS(tmeans.gut.vir, file="../Data/tmeans_gut_vir.RDS")
tmeans.gut.mg <- tmeans.gut[colnames(tmeans.gut) %in% subset(data_reform, env2=="Human gut" & method2 == "Mixed MG")$Sample]
saveRDS(tmeans.gut.mg, file="../Data/tmeans_gut_mg.RDS")
tmeans.fw.vir <- tmeans.fw[colnames(tmeans.fw) %in% subset(data_reform, env2=="Freshwater" & method2 == "Virome")$Sample]
saveRDS(tmeans.fw.vir, file="../Data/tmeans_fw_vir.RDS")
tmeans.fw.mg <- tmeans.fw[colnames(tmeans.fw) %in% subset(data_reform, env2=="Freshwater" & method2 == "Mixed MG")$Sample]
saveRDS(tmeans.fw.mg, file="../Data/tmeans_fw_mg.RDS")
tmeans.mar.vir <- tmeans.mar[colnames(tmeans.mar) %in% subset(data_reform, env2=="Marine" & method2 == "Virome")$Sample]
saveRDS(tmeans.mar.vir, file="../Data/tmeans_mar_vir.RDS")
tmeans.mar.mg <- tmeans.mar[colnames(tmeans.mar) %in% subset(data_reform, env2=="Marine" & method2 == "Mixed MG")$Sample]
saveRDS(tmeans.mar.mg, file="../Data/tmeans_mar_mg.RDS")
tmeans.soil.vir <- tmeans.soil[colnames(tmeans.soil) %in% subset(data_reform, env2=="Soil" & method2 == "Virome")$Sample]
saveRDS(tmeans.soil.vir, file="../Data/tmeans_soil_vir.RDS")
tmeans.soil.mg <- tmeans.soil[colnames(tmeans.soil) %in% subset(data_reform, env2=="Soil" & method2 == "Mixed MG")$Sample]
saveRDS(tmeans.soil.mg, file="../Data/tmeans_soil_mg.RDS")
```

# Get differences in genome coverage (breadth) between viromes and metagenomes across all samples

## Get vMAGs from every environment that were detected in both viromes and metagenomes

## Format and filter breadths data

``` r
breadths <- rbind(melt(breadth.fw %>%
                         rownames_to_column(var = "Genome")) %>%
                    rename(Sample=variable, Breadth=value),
                  melt(breadth.gut %>%
                         rownames_to_column(var = "Genome")) %>%
                    rename(Sample=variable, Breadth=value),
                  melt(breadth.mar %>%
                         rownames_to_column(var = "Genome")) %>%
                    rename(Sample=variable, Breadth=value),
                  melt(breadth.soil %>%
                         rownames_to_column(var = "Genome")) %>%
                    rename(Sample=variable, Breadth=value)) %>%
  inner_join(data_reform %>%
               select(Sample, sample_source, Environment, Method)) %>%
  inner_join(vir_info %>%
               select(genome, completeness),
             by = join_by(Genome==genome)) %>%
  filter(Breadth >= 0.75, completeness==100) 

breadths <- breadths %>%
  filter(Genome %in% subset(breadths, Method == "virome")$Genome & Genome %in% subset(breadths, Method == "metagenome")$Genome) %>%
  select(-completeness) %>%
  mutate(Method = case_when(Method == "virome" ~ "Virome",
                            Method == "metagenome" ~ "Mixed MG"),
         Environment = case_when(Environment == "human_gut" ~ "Human Gut",
                                 Environment == "freshwater" ~ "Freshwater",
                                 Environment == "marine" ~ "Marine",
                                 Environment == "soil" ~ "Soil"))

saveRDS(breadths, file="../Data/breadths.RDS")
head(breadths)
```

    ##                           Genome    Sample   Breadth       sample_source
    ## 1 Ga0485159__vRhyme_unbinned_300 Ga0485186 0.9257947 ME_2020-10-19_23.5m
    ## 2      Ga0485175__vRhyme_bin_121 Ga0485186 0.8385605 ME_2020-10-19_23.5m
    ## 3 Ga0485173__vRhyme_unbinned_394 Ga0485186 0.8625234 ME_2020-10-19_23.5m
    ## 4 Ga0485178__vRhyme_unbinned_862 Ga0485186 0.7603572 ME_2020-10-19_23.5m
    ## 5      Ga0485175__vRhyme_bin_447 Ga0485186 0.9271312 ME_2020-10-19_23.5m
    ## 6       Ga0485159__vRhyme_bin_64 Ga0485186 0.8584101 ME_2020-10-19_23.5m
    ##   Environment Method
    ## 1  Freshwater Virome
    ## 2  Freshwater Virome
    ## 3  Freshwater Virome
    ## 4  Freshwater Virome
    ## 5  Freshwater Virome
    ## 6  Freshwater Virome

# Depth per position for highlighted vMAGs

## Load and format coverage tables

``` r
Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov <- read.csv("../Tables/Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov", header=FALSE, sep="\t")
colnames(Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov) <- c("Contig", "Position", "Depth")

Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov <- read.csv("../Tables/Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov", header=FALSE, sep="\t")
colnames(Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov) <- c("Contig", "Position", "Depth")

Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov <- read.csv("../Tables/Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov", header=FALSE, sep="\t")
colnames(Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov) <- c("Contig", "Position", "Depth")

Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov <- read.csv("../Tables/Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov", header=FALSE, sep="\t")
colnames(Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov) <- c("Contig", "Position", "Depth")
```

## Add a relative position column to account for the incompletness in the Ga0485172\_\_vRhyme_unbinned_38 assembly

### Define a function to modify the Ga0485172\_\_vRhyme_unbinned_38 dataframes

``` r
# Define the offset for assembly B
offset_B <- 2789

# Define the range that is missing in assembly B
missing_range_start <- 6822 # Reversed
missing_range_end <- 6692 # Reversed
additional_missing_start <- 3905
additional_missing_end <- 3908

# Function to calculate relative position for assembly B
calculate_relative_position <- function(position_B) {
  if (position_B < missing_range_start) {
    return(position_B + offset_B - 1)
  } else if (position_B >= missing_range_end) {
    return(position_B + offset_B - (missing_range_end - missing_range_start + 1))
  } else {
    return(NA)  # Positions in the missing range are not valid in assembly B
  }
}

calculate_relative_position_2 <- function(position_A) {
  if (position_A >= additional_missing_end) {
    return(position_A - (additional_missing_end - additional_missing_start + 1))
  } else {
    return(position_A)  # Positions in the missing range are not valid in assembly A
  }
}
```

### Execute the fucntion for Ga0485184_to_Ga0485172\_\_vRhyme_unbinned_38.depth_per_base.cov and Ga0485172_to_Ga0485172\_\_vRhyme_unbinned_38.depth_per_base.cov

``` r
# Apply the function to the position column in assembly B dataframe
Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Relative.position <- sapply(Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Position, calculate_relative_position)
Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Relative.position <- sapply(Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Position, calculate_relative_position)
```

### Add relative position to Ga0485172_to_Ga0485184\_\_vRhyme_unbinned_566.depth_per_base.cov and Ga0485184_to_Ga0485184\_\_vRhyme_unbinned_566.depth_per_base.cov

``` r
Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Relative.position <- sapply(Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Position, calculate_relative_position_2)
Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Relative.position <- sapply(Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Position, calculate_relative_position_2)
```

### Add categories and normnalized depths to Ga0485184_to_Ga0485172\_\_vRhyme_unbinned_38.depth_per_base.cov and Ga0485172_to_Ga0485172\_\_vRhyme_unbinned_38.depth_per_base.cov

``` r
# Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov
Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth <- rev(Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth) # REVERSE because position plots below include the reverse complement of the reference, reads were mapped to original sequence
Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth.per.100M.reads <- Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth / (data_reform[which(data_reform$Sample == "Ga0485184"),]$Reads_paired_reads / 100000000) # normalize depth by 100M reads in sample
Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth.normalized <- log10(Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth.per.100M.reads) # log10 transform depth axis
Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Genome <- "Ga0485172__vRhyme_unbinned_38"
Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Read.sample <- "Ga0485184"
Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Read.sample.method <- "Virome"
Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Genome.method <- "Mixed MG"
head(Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov)
```

    ##                                  Contig Position Depth Relative.position
    ## 1 vRhyme_unbinned_38__Ga0485172_0000163        1   135              2789
    ## 2 vRhyme_unbinned_38__Ga0485172_0000163        2   211              2790
    ## 3 vRhyme_unbinned_38__Ga0485172_0000163        3   245              2791
    ## 4 vRhyme_unbinned_38__Ga0485172_0000163        4   260              2792
    ## 5 vRhyme_unbinned_38__Ga0485172_0000163        5   288              2793
    ## 6 vRhyme_unbinned_38__Ga0485172_0000163        6   322              2794
    ##   Depth.per.100M.reads Depth.normalized                        Genome
    ## 1             48.68402         1.687386 Ga0485172__vRhyme_unbinned_38
    ## 2             76.09131         1.881335 Ga0485172__vRhyme_unbinned_38
    ## 3             88.35247         1.946219 Ga0485172__vRhyme_unbinned_38
    ## 4             93.76181         1.972026 Ga0485172__vRhyme_unbinned_38
    ## 5            103.85923         2.016445 Ga0485172__vRhyme_unbinned_38
    ## 6            116.12039         2.064908 Ga0485172__vRhyme_unbinned_38
    ##   Read.sample Read.sample.method Genome.method
    ## 1   Ga0485184             Virome      Mixed MG
    ## 2   Ga0485184             Virome      Mixed MG
    ## 3   Ga0485184             Virome      Mixed MG
    ## 4   Ga0485184             Virome      Mixed MG
    ## 5   Ga0485184             Virome      Mixed MG
    ## 6   Ga0485184             Virome      Mixed MG

``` r
# Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov
Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth <- rev(Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth) # REVERSE because position plots below include the reverse complement of the reference, reads were mapped to original sequence
Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth.per.100M.reads <- Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth / (data_reform[which(data_reform$Sample == "Ga0485184"),]$Reads_paired_reads / 100000000) # normalize depth by 100M reads in sample
Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth.normalized <- log10(Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Depth.per.100M.reads) # log10 transform depth axis
Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Genome <- "Ga0485172__vRhyme_unbinned_38"
Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Read.sample <- "Ga0485184"
Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Read.sample.method <- "Mixed MG"
Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov$Genome.method <- "Mixed MG"
head(Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov)
```

    ##                                  Contig Position Depth Relative.position
    ## 1 vRhyme_unbinned_38__Ga0485172_0000163        3     3              2791
    ## 2 vRhyme_unbinned_38__Ga0485172_0000163        4     3              2792
    ## 3 vRhyme_unbinned_38__Ga0485172_0000163        5     3              2793
    ## 4 vRhyme_unbinned_38__Ga0485172_0000163        6     3              2794
    ## 5 vRhyme_unbinned_38__Ga0485172_0000163        7     3              2795
    ## 6 vRhyme_unbinned_38__Ga0485172_0000163        8     4              2796
    ##   Depth.per.100M.reads Depth.normalized                        Genome
    ## 1             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 2             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 3             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 4             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 5             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 6             1.442489       0.15911262 Ga0485172__vRhyme_unbinned_38
    ##   Read.sample Read.sample.method Genome.method
    ## 1   Ga0485184           Mixed MG      Mixed MG
    ## 2   Ga0485184           Mixed MG      Mixed MG
    ## 3   Ga0485184           Mixed MG      Mixed MG
    ## 4   Ga0485184           Mixed MG      Mixed MG
    ## 5   Ga0485184           Mixed MG      Mixed MG
    ## 6   Ga0485184           Mixed MG      Mixed MG

### Add categories and normnalized depths to Ga0485172_to_Ga0485184\_\_vRhyme_unbinned_566.depth_per_base.cov and Ga0485184_to_Ga0485184\_\_vRhyme_unbinned_566.depth_per_base.cov

``` r
# Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov
Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Depth.per.100M.reads <- Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Depth / (data_reform[which(data_reform$Sample == "Ga0485172"),]$Reads_paired_reads / 100000000) # normalize depth by 100M reads in sample
Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Depth.normalized <- log10(Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Depth.per.100M.reads) # log10 transform depth axis
Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Genome <- "Ga0485184__vRhyme_unbinned_566"
Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Read.sample <- "Ga0485172"
Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Read.sample.method <- "Mixed MG"
Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Genome.method <- "Virome"
head(Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov)
```

    ##                                  Contig Position Depth Relative.position
    ## 1 vRhyme_unbinned_566__Ga0485184_000256        6     1                 6
    ## 2 vRhyme_unbinned_566__Ga0485184_000256        7     1                 7
    ## 3 vRhyme_unbinned_566__Ga0485184_000256        8     1                 8
    ## 4 vRhyme_unbinned_566__Ga0485184_000256        9     1                 9
    ## 5 vRhyme_unbinned_566__Ga0485184_000256       10     1                10
    ## 6 vRhyme_unbinned_566__Ga0485184_000256       11     1                11
    ##   Depth.per.100M.reads Depth.normalized                         Genome
    ## 1            0.7677562       -0.1147767 Ga0485184__vRhyme_unbinned_566
    ## 2            0.7677562       -0.1147767 Ga0485184__vRhyme_unbinned_566
    ## 3            0.7677562       -0.1147767 Ga0485184__vRhyme_unbinned_566
    ## 4            0.7677562       -0.1147767 Ga0485184__vRhyme_unbinned_566
    ## 5            0.7677562       -0.1147767 Ga0485184__vRhyme_unbinned_566
    ## 6            0.7677562       -0.1147767 Ga0485184__vRhyme_unbinned_566
    ##   Read.sample Read.sample.method Genome.method
    ## 1   Ga0485172           Mixed MG        Virome
    ## 2   Ga0485172           Mixed MG        Virome
    ## 3   Ga0485172           Mixed MG        Virome
    ## 4   Ga0485172           Mixed MG        Virome
    ## 5   Ga0485172           Mixed MG        Virome
    ## 6   Ga0485172           Mixed MG        Virome

``` r
# Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov
Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Depth.per.100M.reads <- Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Depth / (data_reform[which(data_reform$Sample == "Ga0485172"),]$Reads_paired_reads / 100000000) # normalize depth by 100M reads in sample
Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Depth.normalized <- log10(Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Depth.per.100M.reads) # log10 transform depth axis
Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Genome <- "Ga0485184__vRhyme_unbinned_566"
Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Read.sample <- "Ga0485172"
Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Read.sample.method <- "Virome"
Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov$Genome.method <- "Virome"
head(Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov)
```

    ##                                  Contig Position Depth Relative.position
    ## 1 vRhyme_unbinned_566__Ga0485184_000256        1    56                 1
    ## 2 vRhyme_unbinned_566__Ga0485184_000256        2    81                 2
    ## 3 vRhyme_unbinned_566__Ga0485184_000256        3    88                 3
    ## 4 vRhyme_unbinned_566__Ga0485184_000256        4    98                 4
    ## 5 vRhyme_unbinned_566__Ga0485184_000256        5   107                 5
    ## 6 vRhyme_unbinned_566__Ga0485184_000256        6   120                 6
    ##   Depth.per.100M.reads Depth.normalized                         Genome
    ## 1             42.99435         1.633411 Ga0485184__vRhyme_unbinned_566
    ## 2             62.18825         1.793708 Ga0485184__vRhyme_unbinned_566
    ## 3             67.56254         1.829706 Ga0485184__vRhyme_unbinned_566
    ## 4             75.24011         1.876449 Ga0485184__vRhyme_unbinned_566
    ## 5             82.14991         1.914607 Ga0485184__vRhyme_unbinned_566
    ## 6             92.13074         1.964405 Ga0485184__vRhyme_unbinned_566
    ##   Read.sample Read.sample.method Genome.method
    ## 1   Ga0485172             Virome        Virome
    ## 2   Ga0485172             Virome        Virome
    ## 3   Ga0485172             Virome        Virome
    ## 4   Ga0485172             Virome        Virome
    ## 5   Ga0485172             Virome        Virome
    ## 6   Ga0485172             Virome        Virome

### Combine dataframes for both assemblies

``` r
depth_per_base <- rbind(Ga0485172_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov,
                        Ga0485172_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov,
                        Ga0485184_to_Ga0485184__vRhyme_unbinned_566.depth_per_base.cov,
                        Ga0485184_to_Ga0485172__vRhyme_unbinned_38.depth_per_base.cov)
saveRDS(depth_per_base, file="../Data/depth_per_base.RDS")
head(depth_per_base)
```

    ##                                  Contig Position Depth Relative.position
    ## 1 vRhyme_unbinned_38__Ga0485172_0000163        3     3              2791
    ## 2 vRhyme_unbinned_38__Ga0485172_0000163        4     3              2792
    ## 3 vRhyme_unbinned_38__Ga0485172_0000163        5     3              2793
    ## 4 vRhyme_unbinned_38__Ga0485172_0000163        6     3              2794
    ## 5 vRhyme_unbinned_38__Ga0485172_0000163        7     3              2795
    ## 6 vRhyme_unbinned_38__Ga0485172_0000163        8     4              2796
    ##   Depth.per.100M.reads Depth.normalized                        Genome
    ## 1             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 2             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 3             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 4             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 5             1.081867       0.03417388 Ga0485172__vRhyme_unbinned_38
    ## 6             1.442489       0.15911262 Ga0485172__vRhyme_unbinned_38
    ##   Read.sample Read.sample.method Genome.method
    ## 1   Ga0485184           Mixed MG      Mixed MG
    ## 2   Ga0485184           Mixed MG      Mixed MG
    ## 3   Ga0485184           Mixed MG      Mixed MG
    ## 4   Ga0485184           Mixed MG      Mixed MG
    ## 5   Ga0485184           Mixed MG      Mixed MG
    ## 6   Ga0485184           Mixed MG      Mixed MG

# Gene annotations for highlighted vMAGs

``` r
features.Ga0485184__vRhyme_unbinned_566 <- read.csv("../Tables/Ga0485184__vRhyme_unbinned_566_cds_final_merged_output.tsv", header=TRUE, sep="\t")
features.Ga0485184__vRhyme_unbinned_566$Genome <- "Ga0485184__vRhyme_unbinned_566"
features.Ga0485184__vRhyme_unbinned_566$Genome.method <- "Virome"
features.Ga0485184__vRhyme_unbinned_566$color[features.Ga0485184__vRhyme_unbinned_566$color == 'None'] <- '#c9c9c9' # 'None' is invalid, change to grey
features.Ga0485184__vRhyme_unbinned_566$start.relative <- sapply(features.Ga0485184__vRhyme_unbinned_566$start, calculate_relative_position_2) # Apply the function to get the relative position
features.Ga0485184__vRhyme_unbinned_566$stop.relative <- sapply(features.Ga0485184__vRhyme_unbinned_566$stop, calculate_relative_position_2) # Apply the function to get the relative position
head(features.Ga0485184__vRhyme_unbinned_566)
```

    ##                gene start stop frame                                contig
    ## 1 EGNDLYYE_CDS_0001     1 1491     + vRhyme_unbinned_566__Ga0485184_000256
    ## 2 EGNDLYYE_CDS_0002  1491 1748     + vRhyme_unbinned_566__Ga0485184_000256
    ## 3 EGNDLYYE_CDS_0003  1748 1918     + vRhyme_unbinned_566__Ga0485184_000256
    ## 4 EGNDLYYE_CDS_0004  1921 2118     + vRhyme_unbinned_566__Ga0485184_000256
    ## 5 EGNDLYYE_CDS_0005  2118 2258     + vRhyme_unbinned_566__Ga0485184_000256
    ## 6 EGNDLYYE_CDS_0006  2260 4440     + vRhyme_unbinned_566__Ga0485184_000256
    ##           score mmseqs_phrog mmseqs_alnScore mmseqs_seqIdentity
    ## 1 -1.019086e+06         5902           117.0               0.26
    ## 2 -3.201130e+00        10346            64.0               0.52
    ## 3 -2.131409e+00        16028            84.0              0.857
    ## 4 -6.910225e+00     No_PHROG        No_PHROG           No_PHROG
    ## 5 -9.105793e-01         2423            53.0              0.595
    ## 6 -5.480994e+10          156           302.0              0.329
    ##              mmseqs_eVal mmseqs_top_hit pyhmmer_phrog pyhmmer_bitscore
    ## 1 2.4770000000000003e-30 p7962 VI_07028          5902       161.420181
    ## 2              1.742e-16 NC_008562_p164         10346        87.035629
    ## 3 3.1169999999999998e-24  NC_025461_p45         16028        89.217369
    ## 4               No_PHROG       No_PHROG No_PHROGs_HMM    No_PHROGs_HMM
    ## 5              4.017e-13   KY435490_p57          2423        62.463036
    ## 6 2.2119999999999998e-88   NC_016767_p7           156       446.008759
    ##            pyhmmer_evalue custom_hmm_id custom_hmm_bitscore custom_hmm_evalue
    ## 1   6.999983409973393e-47 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 2  2.4109747012294182e-24 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 3  3.0277880929044873e-25 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 4           No_PHROGs_HMM No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 5   5.732501544624617e-17 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 6 4.6952711300628354e-133 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ##      phrog    Method Region   color                annot           category
    ## 1     5902 PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 2    10346 PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 3    16028 PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 4 No_PHROG PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 5     2423 PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 6      156 PHANOTATE    CDS #3e83f6       portal protein head and packaging
    ##   vfdb_hit vfdb_alnScore vfdb_seqIdentity vfdb_eVal vfdb_short_name
    ## 1     None          None             None      None            None
    ## 2     None          None             None      None            None
    ## 3     None          None             None      None            None
    ## 4     None          None             None      None            None
    ## 5     None          None             None      None            None
    ## 6     None          None             None      None            None
    ##   vfdb_description vfdb_species CARD_hit CARD_alnScore CARD_seqIdentity
    ## 1             None         None     None          None             None
    ## 2             None         None     None          None             None
    ## 3             None         None     None          None             None
    ## 4             None         None     None          None             None
    ## 5             None         None     None          None             None
    ## 6             None         None     None          None             None
    ##   CARD_eVal CARD_species ARO_Accession CARD_short_name Protein_Accession
    ## 1      None         None          None            None              None
    ## 2      None         None          None            None              None
    ## 3      None         None          None            None              None
    ## 4      None         None          None            None              None
    ## 5      None         None          None            None              None
    ## 6      None         None          None            None              None
    ##   DNA_Accession AMR_Gene_Family Drug_Class Resistance_Mechanism
    ## 1          None            None       None                 None
    ## 2          None            None       None                 None
    ## 3          None            None       None                 None
    ## 4          None            None       None                 None
    ## 5          None            None       None                 None
    ## 6          None            None       None                 None
    ##                           Genome Genome.method start.relative stop.relative
    ## 1 Ga0485184__vRhyme_unbinned_566        Virome              1          1491
    ## 2 Ga0485184__vRhyme_unbinned_566        Virome           1491          1748
    ## 3 Ga0485184__vRhyme_unbinned_566        Virome           1748          1918
    ## 4 Ga0485184__vRhyme_unbinned_566        Virome           1921          2118
    ## 5 Ga0485184__vRhyme_unbinned_566        Virome           2118          2258
    ## 6 Ga0485184__vRhyme_unbinned_566        Virome           2260          4436

``` r
features.Ga0485172__vRhyme_unbinned_38 <- read.csv("../Tables/Ga0485172__vRhyme_unbinned_38_reversed_cds_final_merged_output.tsv", header=TRUE, sep="\t")
features.Ga0485172__vRhyme_unbinned_38$Genome <- "Ga0485172__vRhyme_unbinned_38"
features.Ga0485172__vRhyme_unbinned_38$Genome.method <- "Mixed MG"
features.Ga0485172__vRhyme_unbinned_38$color[features.Ga0485172__vRhyme_unbinned_38$color == 'None'] <- '#c9c9c9' # 'None' is invalid, change to grey
features.Ga0485172__vRhyme_unbinned_38$start.relative <- sapply(features.Ga0485172__vRhyme_unbinned_38$start, calculate_relative_position) # Apply the function to get the relative position
features.Ga0485172__vRhyme_unbinned_38$stop.relative <- sapply(features.Ga0485172__vRhyme_unbinned_38$stop, calculate_relative_position) # Apply the function to get the relative position
head(features.Ga0485172__vRhyme_unbinned_38)
```

    ##                gene start stop frame
    ## 1 VXJRLFAH_CDS_0001     1 1653     +
    ## 2 VXJRLFAH_CDS_0002  1650 2012     +
    ## 3 VXJRLFAH_CDS_0003  2016 2318     +
    ## 4 VXJRLFAH_CDS_0004  2318 2551     +
    ## 5 VXJRLFAH_CDS_0005  2563 4527     +
    ## 6 VXJRLFAH_CDS_0006  4524 4667     +
    ##                                                     contig         score
    ## 1 vRhyme_unbinned_38__Ga0485172_0000163_reverse_complement -1.988937e+08
    ## 2 vRhyme_unbinned_38__Ga0485172_0000163_reverse_complement -3.223830e+02
    ## 3 vRhyme_unbinned_38__Ga0485172_0000163_reverse_complement -7.543475e+00
    ## 4 vRhyme_unbinned_38__Ga0485172_0000163_reverse_complement -6.291279e+00
    ## 5 vRhyme_unbinned_38__Ga0485172_0000163_reverse_complement -1.045630e+06
    ## 6 vRhyme_unbinned_38__Ga0485172_0000163_reverse_complement -8.628040e+00
    ##   mmseqs_phrog mmseqs_alnScore mmseqs_seqIdentity            mmseqs_eVal
    ## 1          156           223.0              0.332              1.755e-62
    ## 2          858            85.0              0.395 1.1939999999999999e-22
    ## 3     No_PHROG        No_PHROG           No_PHROG               No_PHROG
    ## 4     No_PHROG        No_PHROG           No_PHROG               No_PHROG
    ## 5        22896            54.0              0.297              4.073e-10
    ## 6     No_PHROG        No_PHROG           No_PHROG               No_PHROG
    ##   mmseqs_top_hit pyhmmer_phrog pyhmmer_bitscore         pyhmmer_evalue
    ## 1   NC_016767_p7           156       339.562775 7.598456368943292e-101
    ## 2  NC_027364_p55           858        94.858177 1.0638259486086216e-26
    ## 3       No_PHROG No_PHROGs_HMM    No_PHROGs_HMM          No_PHROGs_HMM
    ## 4       No_PHROG No_PHROGs_HMM    No_PHROGs_HMM          No_PHROGs_HMM
    ## 5  NC_021540_p65          9942       115.893623  3.093625409242715e-33
    ## 6       No_PHROG No_PHROGs_HMM    No_PHROGs_HMM          No_PHROGs_HMM
    ##   custom_hmm_id custom_hmm_bitscore custom_hmm_evalue    phrog    Method Region
    ## 1 No_custom_HMM       No_custom_HMM     No_custom_HMM      156 PHANOTATE    CDS
    ## 2 No_custom_HMM       No_custom_HMM     No_custom_HMM      858 PHANOTATE    CDS
    ## 3 No_custom_HMM       No_custom_HMM     No_custom_HMM No_PHROG PHANOTATE    CDS
    ## 4 No_custom_HMM       No_custom_HMM     No_custom_HMM No_PHROG PHANOTATE    CDS
    ## 5 No_custom_HMM       No_custom_HMM     No_custom_HMM    22896 PHANOTATE    CDS
    ## 6 No_custom_HMM       No_custom_HMM     No_custom_HMM No_PHROG PHANOTATE    CDS
    ##     color                                          annot
    ## 1 #3e83f6                                 portal protein
    ## 2 #a861e3 starvation-inducible transcriptional regulator
    ## 3 #c9c9c9                           hypothetical protein
    ## 4 #c9c9c9                           hypothetical protein
    ## 5 #07e9a2                                   tail protein
    ## 6 #c9c9c9                           hypothetical protein
    ##                   category vfdb_hit vfdb_alnScore vfdb_seqIdentity vfdb_eVal
    ## 1       head and packaging     None          None             None      None
    ## 2 transcription regulation     None          None             None      None
    ## 3         unknown function     None          None             None      None
    ## 4         unknown function     None          None             None      None
    ## 5                     tail     None          None             None      None
    ## 6         unknown function     None          None             None      None
    ##   vfdb_short_name vfdb_description vfdb_species CARD_hit CARD_alnScore
    ## 1            None             None         None     None          None
    ## 2            None             None         None     None          None
    ## 3            None             None         None     None          None
    ## 4            None             None         None     None          None
    ## 5            None             None         None     None          None
    ## 6            None             None         None     None          None
    ##   CARD_seqIdentity CARD_eVal CARD_species ARO_Accession CARD_short_name
    ## 1             None      None         None          None            None
    ## 2             None      None         None          None            None
    ## 3             None      None         None          None            None
    ## 4             None      None         None          None            None
    ## 5             None      None         None          None            None
    ## 6             None      None         None          None            None
    ##   Protein_Accession DNA_Accession AMR_Gene_Family Drug_Class
    ## 1              None          None            None       None
    ## 2              None          None            None       None
    ## 3              None          None            None       None
    ## 4              None          None            None       None
    ## 5              None          None            None       None
    ## 6              None          None            None       None
    ##   Resistance_Mechanism                        Genome Genome.method
    ## 1                 None Ga0485172__vRhyme_unbinned_38      Mixed MG
    ## 2                 None Ga0485172__vRhyme_unbinned_38      Mixed MG
    ## 3                 None Ga0485172__vRhyme_unbinned_38      Mixed MG
    ## 4                 None Ga0485172__vRhyme_unbinned_38      Mixed MG
    ## 5                 None Ga0485172__vRhyme_unbinned_38      Mixed MG
    ## 6                 None Ga0485172__vRhyme_unbinned_38      Mixed MG
    ##   start.relative stop.relative
    ## 1           2789          4441
    ## 2           4438          4800
    ## 3           4804          5106
    ## 4           5106          5339
    ## 5           5351          7315
    ## 6           7312          7455

``` r
features <- rbind(features.Ga0485184__vRhyme_unbinned_566, features.Ga0485172__vRhyme_unbinned_38)
saveRDS(features, file="../Data/features_combined.RDS")
head(features)
```

    ##                gene start stop frame                                contig
    ## 1 EGNDLYYE_CDS_0001     1 1491     + vRhyme_unbinned_566__Ga0485184_000256
    ## 2 EGNDLYYE_CDS_0002  1491 1748     + vRhyme_unbinned_566__Ga0485184_000256
    ## 3 EGNDLYYE_CDS_0003  1748 1918     + vRhyme_unbinned_566__Ga0485184_000256
    ## 4 EGNDLYYE_CDS_0004  1921 2118     + vRhyme_unbinned_566__Ga0485184_000256
    ## 5 EGNDLYYE_CDS_0005  2118 2258     + vRhyme_unbinned_566__Ga0485184_000256
    ## 6 EGNDLYYE_CDS_0006  2260 4440     + vRhyme_unbinned_566__Ga0485184_000256
    ##           score mmseqs_phrog mmseqs_alnScore mmseqs_seqIdentity
    ## 1 -1.019086e+06         5902           117.0               0.26
    ## 2 -3.201130e+00        10346            64.0               0.52
    ## 3 -2.131409e+00        16028            84.0              0.857
    ## 4 -6.910225e+00     No_PHROG        No_PHROG           No_PHROG
    ## 5 -9.105793e-01         2423            53.0              0.595
    ## 6 -5.480994e+10          156           302.0              0.329
    ##              mmseqs_eVal mmseqs_top_hit pyhmmer_phrog pyhmmer_bitscore
    ## 1 2.4770000000000003e-30 p7962 VI_07028          5902       161.420181
    ## 2              1.742e-16 NC_008562_p164         10346        87.035629
    ## 3 3.1169999999999998e-24  NC_025461_p45         16028        89.217369
    ## 4               No_PHROG       No_PHROG No_PHROGs_HMM    No_PHROGs_HMM
    ## 5              4.017e-13   KY435490_p57          2423        62.463036
    ## 6 2.2119999999999998e-88   NC_016767_p7           156       446.008759
    ##            pyhmmer_evalue custom_hmm_id custom_hmm_bitscore custom_hmm_evalue
    ## 1   6.999983409973393e-47 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 2  2.4109747012294182e-24 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 3  3.0277880929044873e-25 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 4           No_PHROGs_HMM No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 5   5.732501544624617e-17 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ## 6 4.6952711300628354e-133 No_custom_HMM       No_custom_HMM     No_custom_HMM
    ##      phrog    Method Region   color                annot           category
    ## 1     5902 PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 2    10346 PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 3    16028 PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 4 No_PHROG PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 5     2423 PHANOTATE    CDS #c9c9c9 hypothetical protein   unknown function
    ## 6      156 PHANOTATE    CDS #3e83f6       portal protein head and packaging
    ##   vfdb_hit vfdb_alnScore vfdb_seqIdentity vfdb_eVal vfdb_short_name
    ## 1     None          None             None      None            None
    ## 2     None          None             None      None            None
    ## 3     None          None             None      None            None
    ## 4     None          None             None      None            None
    ## 5     None          None             None      None            None
    ## 6     None          None             None      None            None
    ##   vfdb_description vfdb_species CARD_hit CARD_alnScore CARD_seqIdentity
    ## 1             None         None     None          None             None
    ## 2             None         None     None          None             None
    ## 3             None         None     None          None             None
    ## 4             None         None     None          None             None
    ## 5             None         None     None          None             None
    ## 6             None         None     None          None             None
    ##   CARD_eVal CARD_species ARO_Accession CARD_short_name Protein_Accession
    ## 1      None         None          None            None              None
    ## 2      None         None          None            None              None
    ## 3      None         None          None            None              None
    ## 4      None         None          None            None              None
    ## 5      None         None          None            None              None
    ## 6      None         None          None            None              None
    ##   DNA_Accession AMR_Gene_Family Drug_Class Resistance_Mechanism
    ## 1          None            None       None                 None
    ## 2          None            None       None                 None
    ## 3          None            None       None                 None
    ## 4          None            None       None                 None
    ## 5          None            None       None                 None
    ## 6          None            None       None                 None
    ##                           Genome Genome.method start.relative stop.relative
    ## 1 Ga0485184__vRhyme_unbinned_566        Virome              1          1491
    ## 2 Ga0485184__vRhyme_unbinned_566        Virome           1491          1748
    ## 3 Ga0485184__vRhyme_unbinned_566        Virome           1748          1918
    ## 4 Ga0485184__vRhyme_unbinned_566        Virome           1921          2118
    ## 5 Ga0485184__vRhyme_unbinned_566        Virome           2118          2258
    ## 6 Ga0485184__vRhyme_unbinned_566        Virome           2260          4436

## Save dataframe representing shaded absent regions

``` r
rect <- data.frame(xmin = c(missing_range_start, 0), xmax = c(missing_range_end, 0), 
                 ymin = c(-0.75, 0), ymax = c(max(depth_per_base$Depth.normalized), 0), 
                 alpha = c(1, 1),
                 fill = c("white", "white"))
head(rect)
```

    ##   xmin xmax  ymin     ymax alpha  fill
    ## 1 6822 6692 -0.75 3.401494     1 white
    ## 2    0    0  0.00 0.000000     1 white

``` r
saveRDS(rect, file="../Data/features_absent.RDS")
```

# Gene counts

``` r
counts.gut <- read.csv('../Tables/gene_counts_human_gut.csv', header = TRUE, row.names = 1)
saveRDS(counts.gut, file = "../Data/gene_counts_human_gut.RDS")
head(counts.gut)
```

    ##                                       SRR9162900 SRR9161506 SRR9161509
    ## vRhyme_unbinned_19__SRR9162906_435_1           0          0          0
    ## vRhyme_unbinned_19__SRR9162906_435_3           2        309        330
    ## vRhyme_unbinned_19__SRR9162906_435_13          0         10          0
    ## vRhyme_unbinned_53__SRR9162906_141_2           0          4          0
    ## vRhyme_unbinned_53__SRR9162906_141_3           0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_14          0          0          0
    ##                                       SRR9161503 SRR9162908 SRR9162903
    ## vRhyme_unbinned_19__SRR9162906_435_1           0          0          0
    ## vRhyme_unbinned_19__SRR9162906_435_3          81        160         63
    ## vRhyme_unbinned_19__SRR9162906_435_13          0         46          0
    ## vRhyme_unbinned_53__SRR9162906_141_2           0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_3           0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_14          0          0          0
    ##                                       SRR9162902 SRR9161501 SRR9162899
    ## vRhyme_unbinned_19__SRR9162906_435_1           0          0          0
    ## vRhyme_unbinned_19__SRR9162906_435_3          95         34         59
    ## vRhyme_unbinned_19__SRR9162906_435_13          0          0          9
    ## vRhyme_unbinned_53__SRR9162906_141_2           0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_3           0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_14          0          0          0
    ##                                       SRR9161504 SRR9162904 SRR9161505
    ## vRhyme_unbinned_19__SRR9162906_435_1           0          0          0
    ## vRhyme_unbinned_19__SRR9162906_435_3          38          4          3
    ## vRhyme_unbinned_19__SRR9162906_435_13          0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_2           0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_3           0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_14          0          0          0
    ##                                       SRR9162907 SRR9161507 SRR9162906
    ## vRhyme_unbinned_19__SRR9162906_435_1           0          0          1
    ## vRhyme_unbinned_19__SRR9162906_435_3         121         49        158
    ## vRhyme_unbinned_19__SRR9162906_435_13          2          0         22
    ## vRhyme_unbinned_53__SRR9162906_141_2           0          0          6
    ## vRhyme_unbinned_53__SRR9162906_141_3           0          0          3
    ## vRhyme_unbinned_53__SRR9162906_141_14          0          0          0
    ##                                       SRR9162901 SRR9161510 SRR9161502
    ## vRhyme_unbinned_19__SRR9162906_435_1           0          0          0
    ## vRhyme_unbinned_19__SRR9162906_435_3          57          0        400
    ## vRhyme_unbinned_19__SRR9162906_435_13          4          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_2           0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_3           0          0          0
    ## vRhyme_unbinned_53__SRR9162906_141_14          0          0          0
    ##                                       SRR9162905 SRR9161500
    ## vRhyme_unbinned_19__SRR9162906_435_1           0          0
    ## vRhyme_unbinned_19__SRR9162906_435_3          61          0
    ## vRhyme_unbinned_19__SRR9162906_435_13          5          0
    ## vRhyme_unbinned_53__SRR9162906_141_2           4          0
    ## vRhyme_unbinned_53__SRR9162906_141_3           0          0
    ## vRhyme_unbinned_53__SRR9162906_141_14          1          0

``` r
counts.fw <- read.csv('../Tables/gene_counts_freshwater.csv', header = TRUE, row.names = 1)
saveRDS(counts.fw, file = "../Data/gene_counts_freshwater.RDS")
head(counts.fw)
```

    ##                                           Ga0485173 Ga0485185 Ga0485166
    ## vRhyme_unbinned_369__Ga0485159_0000068_5       1092         2         1
    ## vRhyme_unbinned_369__Ga0485159_0000068_9       2774        11         1
    ## vRhyme_unbinned_369__Ga0485159_0000068_11       464         2         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_14       707         1         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_22       156         3         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_29       596         1         0
    ##                                           Ga0485170 Ga0485162 Ga0485181
    ## vRhyme_unbinned_369__Ga0485159_0000068_5          0         4         5
    ## vRhyme_unbinned_369__Ga0485159_0000068_9          0         6         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_11         0         1         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_14         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_22         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_29         0         0         0
    ##                                           Ga0485161 Ga0485169 Ga0485171
    ## vRhyme_unbinned_369__Ga0485159_0000068_5          0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_9          0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_11         0         1         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_14         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_22         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_29         0         0         0
    ##                                           Ga0485178 Ga0485158 Ga0485176
    ## vRhyme_unbinned_369__Ga0485159_0000068_5          0         1         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_9          0         2         2
    ## vRhyme_unbinned_369__Ga0485159_0000068_11         0         0         1
    ## vRhyme_unbinned_369__Ga0485159_0000068_14         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_22         0         0         1
    ## vRhyme_unbinned_369__Ga0485159_0000068_29         0         0         5
    ##                                           Ga0485167 Ga0485182 Ga0485165
    ## vRhyme_unbinned_369__Ga0485159_0000068_5          0         7         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_9          0        18         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_11         0         1         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_14         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_22         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_29         0         1         0
    ##                                           Ga0485180 Ga0485168 Ga0485174
    ## vRhyme_unbinned_369__Ga0485159_0000068_5          0         0         1
    ## vRhyme_unbinned_369__Ga0485159_0000068_9          1         0         6
    ## vRhyme_unbinned_369__Ga0485159_0000068_11         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_14         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_22         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_29         0         1         3
    ##                                           Ga0485184 Ga0485177 Ga0485175
    ## vRhyme_unbinned_369__Ga0485159_0000068_5          2         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_9          3         7         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_11         3         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_14         0         5         1
    ## vRhyme_unbinned_369__Ga0485159_0000068_22         0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_29         0         1         0
    ##                                           Ga0485157 Ga0485172 Ga0485186
    ## vRhyme_unbinned_369__Ga0485159_0000068_5          0         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_9          2         0         4
    ## vRhyme_unbinned_369__Ga0485159_0000068_11         1         2         8
    ## vRhyme_unbinned_369__Ga0485159_0000068_14         0         0         8
    ## vRhyme_unbinned_369__Ga0485159_0000068_22         0         0         2
    ## vRhyme_unbinned_369__Ga0485159_0000068_29         0         0         2
    ##                                           Ga0485159 Ga0485164 Ga0485183
    ## vRhyme_unbinned_369__Ga0485159_0000068_5        129         1         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_9        260         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_11        42         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_14        27         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_22        20         0         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_29        39         0         0
    ##                                           Ga0485179
    ## vRhyme_unbinned_369__Ga0485159_0000068_5          0
    ## vRhyme_unbinned_369__Ga0485159_0000068_9          0
    ## vRhyme_unbinned_369__Ga0485159_0000068_11         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_14         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_22         0
    ## vRhyme_unbinned_369__Ga0485159_0000068_29         0

``` r
counts.mar <- read.csv('../Tables/gene_counts_marine.csv', header = TRUE, row.names = 1)
saveRDS(counts.mar, file = "../Data/gene_counts_marine.RDS")
head(counts.mar)
```

    ##                                                                     ERR594382
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599146
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594411
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594362
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR598984
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         6
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         5
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33        22
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         4
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     2
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                   897
    ##                                                                     ERR599006
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599173
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594377
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599144
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594392
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594391
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                    21
    ##                                                                     ERR599165
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30       163
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32        88
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33       168
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34        32
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                   227
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                   413
    ##                                                                     ERR594354
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     4
    ##                                                                     ERR599110
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599174
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30       107
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32        39
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33       159
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34        15
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     7
    ##                                                                     ERR594355
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594388
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594353
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                  1259
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                  2110
    ##                                                                     ERR599065
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         2
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                   660
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                  1085
    ##                                                                     ERR599122
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599044
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30      2237
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32       997
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33      2551
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34       453
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599133
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30        81
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32        14
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33        80
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34        20
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     4
    ##                                                                     ERR594409
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599148
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30       165
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32        68
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33       205
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34        40
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                    14
    ##                                                                     ERR599017
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30       141
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32        63
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33       149
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34        25
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     3
    ##                                                                     ERR599023
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594379
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                    12
    ##                                                                     ERR599126
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         8
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         2
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33        46
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34        17
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     1
    ##                                                                     ERR599118
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                   308
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                   339
    ##                                                                     ERR594364
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                    98
    ##                                                                     ERR594389
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594415
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594380
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     2
    ##                                                                     ERR594407
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     2
    ##                                                                     ERR594404
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599150
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR598982
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599121
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594385
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR594412
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                    40
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                    70
    ##                                                                     ERR594359
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33         0
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34         0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                     0
    ##                                                                     ERR599005
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30      1082
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_32       484
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_33      1245
    ## vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_34       236
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_1                     0
    ## vRhyme_unbinned_8__NODE_256_length_20384_cov_0.136554_2                   257

``` r
counts.soil <- read.csv('../Tables/gene_counts_soil.csv', header = TRUE, row.names = 1)
saveRDS(counts.soil, file = "../Data/gene_counts_soil.RDS")
head(counts.soil)
```

    ##                                     SRR8487032 SRR8487034 SRR8487037 SRR8487027
    ## vRhyme_bin_2__SRR8487012_24_5                0         38          0          0
    ## vRhyme_bin_2__SRR8487012_24_6                0          4          0          0
    ## vRhyme_unbinned_1__SRR8487012_17_4           0        166         13         24
    ## vRhyme_unbinned_1__SRR8487012_17_5           0         97         13          1
    ## vRhyme_unbinned_1__SRR8487012_17_7           0         41          7          3
    ## vRhyme_unbinned_1__SRR8487012_17_11          0          4          0          0
    ##                                     SRR8487038 SRR8487024 SRR8487036 SRR8487030
    ## vRhyme_bin_2__SRR8487012_24_5                0          0         39          0
    ## vRhyme_bin_2__SRR8487012_24_6                0          0          3          0
    ## vRhyme_unbinned_1__SRR8487012_17_4           0          0         73          6
    ## vRhyme_unbinned_1__SRR8487012_17_5           0          0         82          8
    ## vRhyme_unbinned_1__SRR8487012_17_7           0          0         19          0
    ## vRhyme_unbinned_1__SRR8487012_17_11          0          0          8          0
    ##                                     SRR8487040 SRR8487025 SRR8487011 SRR8487010
    ## vRhyme_bin_2__SRR8487012_24_5                0          0         23          5
    ## vRhyme_bin_2__SRR8487012_24_6                0          0          0          2
    ## vRhyme_unbinned_1__SRR8487012_17_4           6          0         22        126
    ## vRhyme_unbinned_1__SRR8487012_17_5           0          0         18        108
    ## vRhyme_unbinned_1__SRR8487012_17_7           0          0          8         29
    ## vRhyme_unbinned_1__SRR8487012_17_11          0          0          3          9
    ##                                     SRR8487018 SRR8487016 SRR8487035 SRR8487023
    ## vRhyme_bin_2__SRR8487012_24_5               73        111         64        174
    ## vRhyme_bin_2__SRR8487012_24_6               10         26          0         33
    ## vRhyme_unbinned_1__SRR8487012_17_4         172         48         40         32
    ## vRhyme_unbinned_1__SRR8487012_17_5         120         48         25         31
    ## vRhyme_unbinned_1__SRR8487012_17_7          48         24         10         12
    ## vRhyme_unbinned_1__SRR8487012_17_11         18         10          4          3
    ##                                     SRR8487019 SRR8487026 SRR8487028 SRR8487021
    ## vRhyme_bin_2__SRR8487012_24_5              168          0          0         12
    ## vRhyme_bin_2__SRR8487012_24_6               33          0          0          2
    ## vRhyme_unbinned_1__SRR8487012_17_4          93          0          0        336
    ## vRhyme_unbinned_1__SRR8487012_17_5          86          0          0        195
    ## vRhyme_unbinned_1__SRR8487012_17_7          26          0          0         67
    ## vRhyme_unbinned_1__SRR8487012_17_11          5          0          0         20
    ##                                     SRR8487020 SRR8487039 SRR8487031 SRR8487015
    ## vRhyme_bin_2__SRR8487012_24_5              105          0          0         14
    ## vRhyme_bin_2__SRR8487012_24_6               11          0          0          0
    ## vRhyme_unbinned_1__SRR8487012_17_4         177          0        194         12
    ## vRhyme_unbinned_1__SRR8487012_17_5         164          0        221         19
    ## vRhyme_unbinned_1__SRR8487012_17_7          65          0         71          8
    ## vRhyme_unbinned_1__SRR8487012_17_11          8          0          8          0
    ##                                     SRR8487029 SRR8487017 SRR8487013 SRR8487012
    ## vRhyme_bin_2__SRR8487012_24_5                0        106          0        316
    ## vRhyme_bin_2__SRR8487012_24_6                0          9          0         39
    ## vRhyme_unbinned_1__SRR8487012_17_4           0        102         15        371
    ## vRhyme_unbinned_1__SRR8487012_17_5           0        102          3        248
    ## vRhyme_unbinned_1__SRR8487012_17_7           0         27          3        116
    ## vRhyme_unbinned_1__SRR8487012_17_11          0          7          2         30
    ##                                     SRR8487014 SRR8487022
    ## vRhyme_bin_2__SRR8487012_24_5              639        654
    ## vRhyme_bin_2__SRR8487012_24_6              112         81
    ## vRhyme_unbinned_1__SRR8487012_17_4          50        230
    ## vRhyme_unbinned_1__SRR8487012_17_5          31        143
    ## vRhyme_unbinned_1__SRR8487012_17_7          17         85
    ## vRhyme_unbinned_1__SRR8487012_17_11          6         11

# Pharokka annotations for genes

``` r
annot.pharokka.gut <- read.csv("../Tables/pharokka_human_gut_virus_proteins.csv")
annot.pharokka.gut$gene <- sub("_CDS_0001", "", annot.pharokka.gut$gene)
saveRDS(annot.pharokka.gut, file="../Data/pharokka_human_gut_genes.RDS")
head(annot.pharokka.gut)
```

    ##                                   gene   color                annot
    ## 1 vRhyme_unbinned_19__SRR9162906_435_1    None hypothetical protein
    ## 2 vRhyme_unbinned_19__SRR9162906_435_2    None hypothetical protein
    ## 3 vRhyme_unbinned_19__SRR9162906_435_3 #c9c9c9 hypothetical protein
    ## 4 vRhyme_unbinned_19__SRR9162906_435_4    None hypothetical protein
    ## 5 vRhyme_unbinned_19__SRR9162906_435_5    None hypothetical protein
    ## 6 vRhyme_unbinned_19__SRR9162906_435_6    None hypothetical protein
    ##           category
    ## 1 unknown function
    ## 2 unknown function
    ## 3 unknown function
    ## 4 unknown function
    ## 5 unknown function
    ## 6 unknown function

``` r
annot.pharokka.fw <- read.csv("../Tables/pharokka_freshwater_virus_proteins.csv")
annot.pharokka.fw$gene <- sub("_CDS_0001", "", annot.pharokka.fw$gene)
saveRDS(annot.pharokka.fw, file="../Data/pharokka_freshwater_genes.RDS")
head(annot.pharokka.fw)
```

    ##                                       gene   color
    ## 1 vRhyme_unbinned_369__Ga0485159_0000068_1 #ffdf59
    ## 2 vRhyme_unbinned_369__Ga0485159_0000068_2 #ffdf59
    ## 3 vRhyme_unbinned_369__Ga0485159_0000068_3 #ffdf59
    ## 4 vRhyme_unbinned_369__Ga0485159_0000068_4 #ffdf59
    ## 5 vRhyme_unbinned_369__Ga0485159_0000068_5 #c9c9c9
    ## 6 vRhyme_unbinned_369__Ga0485159_0000068_6 #ff59f5
    ##                                annot
    ## 1     clamp loader of DNA polymerase
    ## 2                DNA binding protein
    ## 3                 DNA endonuclease V
    ## 4 DNA polymerase processivity factor
    ## 5               hypothetical protein
    ## 6             porphyrin biosynthesis
    ##                                            category
    ## 1                DNA, RNA and nucleotide metabolism
    ## 2                DNA, RNA and nucleotide metabolism
    ## 3                DNA, RNA and nucleotide metabolism
    ## 4                DNA, RNA and nucleotide metabolism
    ## 5                                  unknown function
    ## 6 moron, auxiliary metabolic gene and host takeover

``` r
annot.pharokka.mar <- read.csv("../Tables/pharokka_marine_virus_proteins.csv")
annot.pharokka.mar$gene <- sub("_CDS_0001", "", annot.pharokka.mar$gene)
saveRDS(annot.pharokka.mar, file="../Data/pharokka_marine_genes.RDS")
head(annot.pharokka.mar)
```

    ##                                                                  gene   color
    ## 1 vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_26 #ff59f5
    ## 2 vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_27 #fea328
    ## 3 vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_28    None
    ## 4 vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_29 #c9c9c9
    ## 5 vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_30    None
    ## 6 vRhyme_unbinned_19__NODE_79_length_37623_cov_0.132795_fragment_1_31 #c9c9c9
    ##                      annot                                          category
    ## 1 anti-restriction protein moron, auxiliary metabolic gene and host takeover
    ## 2              transposase                          integration and excision
    ## 3     hypothetical protein                                  unknown function
    ## 4     hypothetical protein                                  unknown function
    ## 5     hypothetical protein                                  unknown function
    ## 6     hypothetical protein                                  unknown function

``` r
annot.pharokka.soil <- read.csv("../Tables/pharokka_soil_virus_proteins.csv")
annot.pharokka.soil$gene <- sub("_CDS_0001", "", annot.pharokka.soil$gene)
saveRDS(annot.pharokka.soil, file="../Data/pharokka_soil_genes.RDS")
head(annot.pharokka.soil)
```

    ##                            gene color                annot         category
    ## 1 vRhyme_bin_2__SRR8487012_29_2  None hypothetical protein unknown function
    ## 2 vRhyme_bin_2__SRR8487012_29_3  None hypothetical protein unknown function
    ## 3 vRhyme_bin_2__SRR8487012_29_4  None hypothetical protein unknown function
    ## 4 vRhyme_bin_2__SRR8487012_29_5  None hypothetical protein unknown function
    ## 5 vRhyme_bin_2__SRR8487012_29_7  None hypothetical protein unknown function
    ## 6 vRhyme_bin_2__SRR8487012_29_8  None hypothetical protein unknown function
