# Lab 19: Cancer Mutation Mini-Project


``` r
library(bio3d)
sequence <- read.fasta("A17502778_mutant_seq.fa")
```

``` r
sequence
```

                   1        .         .         .         .         .         60 
    wt_healthy     MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEH
    mutant_tumor   MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEH
                   ************************************************************ 
                   1        .         .         .         .         .         60 

                  61        .         .         .         .         .         120 
    wt_healthy     IEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTV
    mutant_tumor   IEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTV
                   ************************************************************ 
                  61        .         .         .         .         .         120 

                 121        .         .         .         .         .         180 
    wt_healthy     TSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDS
    mutant_tumor   TSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDS
                   ************************************************************ 
                 121        .         .         .         .         .         180 

                 181        .         .         .         .         .         240 
    wt_healthy     LKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRK
    mutant_tumor   LKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRK
                   ************************************************************ 
                 181        .         .         .         .         .         240 

                 241        .         .         .         .         .         300 
    wt_healthy     TFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPI
    mutant_tumor   TFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPI
                   ************************************************************ 
                 241        .         .         .         .         .         300 

                 301        .         .         .         .         .         360 
    wt_healthy     PQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQR
    mutant_tumor   PQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQR
                   ************************************************************ 
                 301        .         .         .         .         .         360 

                 361        .         .         .         .         .         420 
    wt_healthy     DRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSP
    mutant_tumor   DRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSP
                   ************************************************************ 
                 361        .         .         .         .         .         420 

                 421        .         .         .         .         .         480 
    wt_healthy     GPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDV
    mutant_tumor   GPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDV
                   ************************************************************ 
                 421        .         .         .         .         .         480 

                 481        .         .         .         .         .         540 
    wt_healthy     AVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHH
    mutant_tumor   AVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPVLAIVTQWCEGSSLYHH
                   ******************************************* **************** 
                 481        .         .         .         .         .         540 

                 541        .         .         .         .         .         600 
    wt_healthy     LHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATV
    mutant_tumor   LHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATV
                   ************************************************************ 
                 541        .         .         .         .         .         600 

                 601        .         .         .         .         .         660 
    wt_healthy     KSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNIN
    mutant_tumor   KSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAEGIVRYELMTGQLPYSNIN
                   ***************************************** *** ************** 
                 601        .         .         .         .         .         660 

                 661        .         .         .         .         .         720 
    wt_healthy     NRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARS
    mutant_tumor   NRDQIIYMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARS
                   ******^***************************************************** 
                 661        .         .         .         .         .         720 

                 721        .         .         .         .     766 
    wt_healthy     LPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH
    mutant_tumor   LPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH
                   ********************************************** 
                 721        .         .         .         .     766 

    Call:
      read.fasta(file = "A17502778_mutant_seq.fa")

    Class:
      fasta

    Alignment dimensions:
      2 sequence rows; 766 position columns (766 non-gap, 0 gap) 

    + attr: id, ali, call

``` r
wt_healthy <- "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEHIEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTVTSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDSLKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRKTFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPIPQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQRDRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSPGPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH"
mutant_tumor <- "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEHIEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTVTSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDSLKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRKTFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPIPQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQRDRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSPGPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPVLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAEGIVRYELMTGQLPYSNINNRDQIIYMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH"
```

``` r
conserv(sequence)
```

      [1]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [16]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [31]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [46]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [61]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [76]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [91]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [106]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [121]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [136]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [151]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [166]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [181]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [196]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [211]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [226]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [241]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [256]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [271]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [286]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [301]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [316]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [331]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [346]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [361]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [376]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [391]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [406]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [421]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [436]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [451]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [466]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [481]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [496]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [511]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0 -0.3  1.0
    [526]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [541]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [556]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [571]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [586]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [601]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [616]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [631]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0 -0.5  1.0  1.0  1.0
    [646] -0.3  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [661]  1.0  1.0  1.0  1.0  1.0  1.0  0.5  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [676]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [691]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [706]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [721]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [736]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [751]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [766]  1.0

``` r
blast.pdb(sequence)
```

    Warning in blast.pdb(sequence): Multiple sequences detected - using only the
    first sequence in input object

     Searching ... please wait (updates every 5 seconds) RID = V61TE0AB016 
     ...............................
     Reporting 2417 hits

    $hit.tbl
               queryid subjectids identity alignmentlength mismatches gapopens
    1    Query_1615341     6Q0K_A   99.869             766          1        0
    2    Query_1615341     7MFD_A   99.739             766          2        0
    3    Query_1615341     6Q0J_A   99.739             766          2        0
    4    Query_1615341     6NYB_A   99.739             766          2        0
    5    Query_1615341     8VYO_A   99.739             766          2        0
    6    Query_1615341     7ZR0_K   99.869             766          1        0
    7    Query_1615341     6UAN_B   99.608             766          3        0
    8    Query_1615341     8VYV_C   99.608             766          3        0
    9    Query_1615341     8VYU_A   99.608             766          3        0
    10   Query_1615341     8VYP_C   99.608             766          3        0
    11   Query_1615341     8VYS_A   99.478             766          4        0
    12   Query_1615341     8VYW_C   99.200             375          3        0
    13   Query_1615341     8VYQ_C   99.200             375          3        0
    14   Query_1615341     9MMS_A   57.058             673        245       13
    15   Query_1615341     9MMR_A   56.909             673        246       13
    16   Query_1615341     9MMQ_A   56.761             673        247       13
    17   Query_1615341     8CHF_A   60.359             613        206       12
    18   Query_1615341     8U1L_C   56.464             673        249       13
    19   Query_1615341     9MMP_A   56.464             673        249       13
    20   Query_1615341   7Z37_CP1   59.547             618        213       12
    21   Query_1615341     4MNF_A   99.672             305          1        0
    22   Query_1615341     4MNE_B  100.000             295          0        0
    23   Query_1615341     3II5_A  100.000             295          0        0
    24   Query_1615341     3D4Q_A  100.000             295          0        0
    25   Query_1615341     7SHV_A  100.000             295          0        0
    26   Query_1615341     5FD2_A  100.000             294          0        0
    27   Query_1615341     4EHG_A   99.661             295          1        0
    28   Query_1615341     3IDP_A   99.322             295          2        0
    29   Query_1615341     4MBJ_A  100.000             292          0        0
    30   Query_1615341     5HIE_A   98.305             295          0        1
    31   Query_1615341     4DBN_A  100.000             282          0        0
    32   Query_1615341     3Q96_A  100.000             282          0        0
    33   Query_1615341     6PP9_A   99.645             282          1        0
    34   Query_1615341     9ECU_A  100.000             279          0        0
    35   Query_1615341     2FB8_A  100.000             279          0        0
    36   Query_1615341     9AXX_B  100.000             279          0        0
    37   Query_1615341     4H58_A  100.000             275          0        0
    38   Query_1615341     1UWH_A   99.638             276          1        0
    39   Query_1615341     1UWJ_A   99.638             276          1        0
    40   Query_1615341     9BFB_A   98.925             279          3        0
    41   Query_1615341     6N0P_A  100.000             273          0        0
    42   Query_1615341     6U2H_C   94.810             289         15        0
    43   Query_1615341     4YHT_A   98.162             272          5        0
    44   Query_1615341     4XV9_A   94.643             280         15        0
    45   Query_1615341     4JVG_A   94.964             278         14        0
    46   Query_1615341     5CSW_A   94.643             280         15        0
    47   Query_1615341     3C4C_A   94.964             278         14        0
    48   Query_1615341     4RZV_A   94.624             279         15        0
    49   Query_1615341     4WO5_A   94.964             278         14        0
    50   Query_1615341     4FK3_A   94.286             280         16        0
    51   Query_1615341     8F7O_A   93.929             280         17        0
    52   Query_1615341     6XFP_A   94.286             280         16        0
    53   Query_1615341     8C7Y_A   94.604             278         15        0
    54   Query_1615341     4XV1_A   93.571             280         18        0
    55   Query_1615341     6V34_A   93.571             280         18        0
    56   Query_1615341     3OG7_A   93.525             278         18        0
    57   Query_1615341     5ITA_A   93.548             279         18        0
    58   Query_1615341     5JRQ_A   94.526             274         15        0
    59   Query_1615341     6P3D_A   93.190             279         19        0
    60   Query_1615341     4CQE_A   94.526             274         15        0
    61   Query_1615341     7P3V_A   94.853             272         14        0
    62   Query_1615341     5HI2_A   92.857             280         15        1
    63   Query_1615341     9AXC_A   74.854             342         79        4
    64   Query_1615341     9AXA_A   74.561             342         80        4
    65   Query_1615341     3OMV_A   77.627             295         66        0
    66   Query_1615341     9O0U_A   81.004             279         53        0
    67   Query_1615341     9AY7_A   81.004             279         53        0
    68   Query_1615341     8GFT_D   80.000             285         57        0
    69   Query_1615341     9AXM_B   76.703             279         65        0
    70   Query_1615341     8GAE_D   83.495             206         34        0
    71   Query_1615341     6PTS_D   62.308             130         46        1
    72   Query_1615341     7JHP_C   62.308             130         46        1
    73   Query_1615341     6XI7_B   62.308             130         46        1
    74   Query_1615341     2Y4I_B   36.519             293        167        8
    75   Query_1615341     5KKR_B   36.519             293        167        8
    76   Query_1615341     8BW9_D   38.644             295        165        9
    77   Query_1615341     2L05_A  100.000              84          0        0
    78   Query_1615341     7JUQ_B   36.519             293        167        8
    79   Query_1615341     3NY5_A   97.647              85          2        0
    80   Query_1615341     6XGU_B   61.538             130         47        1
    81   Query_1615341     5J17_A   88.542              96          2        1
    82   Query_1615341     7JUW_B   35.254             295        170        8
    83   Query_1615341     5VR3_A   88.298              94          5        1
    84   Query_1615341     9AXH_C   34.982             283        167        7
    85   Query_1615341     3PPZ_A   36.934             287        165       10
    86   Query_1615341     4CSV_A   36.559             279        153       10
    87   Query_1615341     3P86_A   36.585             287        166       10
    88   Query_1615341     5VYK_A  100.000              73          0        0
    89   Query_1615341     8DEG_A   37.729             273        148        8
    90   Query_1615341     5CEN_A   37.729             273        148        8
    91   Query_1615341     9NS1_A   33.574             277        162        9
    92   Query_1615341     8JF3_A   33.574             277        162        9
    93   Query_1615341     2BDF_A   33.574             277        162        9
    94   Query_1615341     7NG7_A   33.574             277        162        9
    95   Query_1615341     7OTE_A   33.574             277        162        9
    96   Query_1615341     4MXO_A   33.574             277        162        9
    97   Query_1615341     6E6E_A   33.818             275        160        9
    98   Query_1615341     8HAQ_A   33.818             275        160        9
    99   Query_1615341     3A4O_X   35.075             268        151       10
    100  Query_1615341     4MXX_A   33.574             277        162        9
    101  Query_1615341     1YOL_A   33.213             277        163        9
    102  Query_1615341     9NS0_A   33.574             277        162        9
    103  Query_1615341     3OEZ_A   33.574             277        162        9
    104  Query_1615341     2OIQ_A   33.574             277        162        9
    105  Query_1615341     1YI6_A   33.818             275        160        9
    106  Query_1615341     2OG8_A   32.852             277        160       10
    107  Query_1615341     1YOJ_A   32.971             276        165        9
    108  Query_1615341     3D7U_B   33.696             276        161        9
    109  Query_1615341     4MXY_A   33.574             277        162        9
    110  Query_1615341     2PL0_A   32.734             278        161       10
    111  Query_1615341     2OFV_A   32.734             278        161       10
    112  Query_1615341     3BYS_A   32.734             278        161       10
    113  Query_1615341     5XY1_A   34.962             266        150       10
    114  Query_1615341     2OF2_A   32.734             278        161       10
    115  Query_1615341     2DQ7_X   34.058             276        160        9
    116  Query_1615341     3GEQ_A   33.574             277        162        9
    117  Query_1615341     3KXZ_A   32.727             275        165        9
    118  Query_1615341     5T0P_A   33.213             277        163        9
    119  Query_1615341     5SWH_A   33.213             277        163        9
    120  Query_1615341     2HWO_A   33.213             277        163        9
    121  Query_1615341     2ZM1_A   32.727             275        165        9
    122  Query_1615341     1QPC_A   32.727             275        165        9
    123  Query_1615341     3SVV_A   33.213             277        163        9
    124  Query_1615341     2OFU_A   32.727             275        165        9
    125  Query_1615341     3KMM_A   32.727             275        165        9
    126  Query_1615341     6PDJ_A   32.727             275        165        9
    127  Query_1615341     3U4W_A   33.700             273        159        9
    128  Query_1615341     4MCV_A   33.213             277        163        9
    129  Query_1615341     3G6H_A   33.213             277        163        9
    130  Query_1615341     3BYO_A   32.727             275        165        9
    131  Query_1615341     3BYM_A   32.727             275        165        9
    132  Query_1615341     4LGH_A   33.333             276        162        9
    133  Query_1615341     2H8H_A   33.574             277        162        9
    134  Query_1615341     8JN8_A   33.574             277        162        9
    135  Query_1615341     9IRL_A   33.574             277        162        9
    136  Query_1615341     2ZV7_A   34.586             266        151       10
    137  Query_1615341     1FMK_A   33.574             277        162        9
    138  Query_1615341     2HK5_A   33.333             270        158        9
    139  Query_1615341     5ZJ6_A   33.333             270        158        9
    140  Query_1615341     1Y57_A   33.574             277        162        9
    141  Query_1615341     2QI8_A   32.852             277        164        9
    142  Query_1615341     9BT8_C   33.935             277        161        9
    143  Query_1615341     3DQW_A   33.091             275        166        9
    144  Query_1615341     3MPM_A   32.246             276        161       10
    145  Query_1615341     8XN8_A   33.213             277        163        9
    146  Query_1615341     6F3F_A   33.213             277        163        9
    147  Query_1615341     4K11_A   33.213             277        163        9
    148  Query_1615341     1KSW_A   33.213             277        163        9
    149  Query_1615341     1QPD_A   32.364             275        166        9
    150  Query_1615341     2PTK_A   33.213             277        163        9
    151  Query_1615341     6IN0_A   33.559             295        173       10
    152  Query_1615341     2QOK_A   33.898             295        172       10
    153  Query_1615341     2QOC_A   33.784             296        171       11
    154  Query_1615341     3DZQ_A   33.559             295        173       10
    155  Query_1615341     1AD5_A   32.721             272        157        8
    156  Query_1615341     2QOO_A   33.559             295        173       10
    157  Query_1615341     2QOI_A   33.559             295        173       10
    158  Query_1615341     2QOF_A   33.559             295        173       10
    159  Query_1615341     2QOD_A   33.559             295        173       10
    160  Query_1615341     3FXX_A   33.559             295        173       10
    161  Query_1615341     2GSF_A   33.559             295        173       10
    162  Query_1615341     1QCF_A   32.727             275        163        9
    163  Query_1615341     2QOL_A   33.559             295        173       10
    164  Query_1615341     9BYJ_A   32.727             275        163        9
    165  Query_1615341     2YN8_A   30.928             291        181        9
    166  Query_1615341     2QOB_A   33.784             296        171       11
    167  Query_1615341     6FNI_A   30.928             291        181        9
    168  Query_1615341     2QON_A   33.559             295        173       10
    169  Query_1615341     2QO7_A   33.559             295        173       10
    170  Query_1615341     6CZ2_A   33.582             268        163        8
    171  Query_1615341     2XYU_A   32.534             292        177        9
    172  Query_1615341     2Y6M_A   32.534             292        177        9
    173  Query_1615341     6CZ3_A   33.582             268        163        8
    174  Query_1615341     2HEL_A   32.534             292        177        9
    175  Query_1615341     3KUL_A   32.997             297        181        9
    176  Query_1615341     2VWU_A   30.605             281        175        9
    177  Query_1615341     8JNA_B   71.605              81         23        0
    178  Query_1615341     3KUL_B   32.886             298        180        9
    179  Query_1615341     3ZEW_A   30.241             291        183        9
    180  Query_1615341     4AW5_A   29.897             291        184        8
    181  Query_1615341     5MJA_A   31.673             281        172        8
    182  Query_1615341     5MJB_A   31.673             281        172        8
    183  Query_1615341     7UY0_B   33.803             284        164       10
    184  Query_1615341     7UY0_A   33.803             284        164       10
    185  Query_1615341     7KPL_A   31.673             281        172        8
    186  Query_1615341     4LGG_A   33.333             261        152        9
    187  Query_1615341     3ZFX_A   31.673             281        172        8
    188  Query_1615341     6UMW_A   31.915             282        170       10
    189  Query_1615341     2R2P_A   32.975             279        169        8
    190  Query_1615341     3ZFY_A   31.973             294        178       10
    191  Query_1615341     3ZFM_A   31.317             281        173        7
    192  Query_1615341     1JPA_A   31.317             281        173        7
    193  Query_1615341     5D7V_A   32.955             264        162        8
    194  Query_1615341     5DA3_A   32.955             264        162        8
    195  Query_1615341     1K3A_A   31.488             289        166        9
    196  Query_1615341     5H2U_A   32.955             264        162        8
    197  Query_1615341     2ZM3_A   31.488             289        166        9
    198  Query_1615341     5I9U_A   31.429             280        173        8
    199  Query_1615341     8XPV_A   31.429             280        173        8
    200  Query_1615341     5FXS_A   31.741             293        160       11
    201  Query_1615341     3O23_A   31.741             293        160       11
    202  Query_1615341     1JQH_A   31.741             293        160       11
    203  Query_1615341     3QQU_A   31.741             293        160       11
    204  Query_1615341     1MQB_A   31.429             280        173        8
    205  Query_1615341     1M7N_A   31.741             293        160       11
    206  Query_1615341     2OJ9_A   31.741             293        160       11
    207  Query_1615341     4TRL_A   32.028             281        168        9
    208  Query_1615341     4D2R_A   31.741             293        160       11
    209  Query_1615341     5FXQ_A   31.741             293        160       11
    210  Query_1615341     3I81_A   31.741             293        160       11
    211  Query_1615341   8PYI_AAA   31.741             293        160       11
    212  Query_1615341     5FXR_A   31.741             293        160       11
    213  Query_1615341     3D94_A   31.741             293        160       11
    214  Query_1615341     4P2K_A   31.541             279        172        8
    215  Query_1615341     5EK7_A   31.541             279        172        8
    216  Query_1615341     8XKP_A   31.525             295        178       10
    217  Query_1615341     7KJA_A   31.429             280        173        8
    218  Query_1615341     6JMF_A   31.525             295        178       10
    219  Query_1615341     7KJC_A   31.429             280        173        8
    220  Query_1615341     7KJB_A   31.429             280        173        8
    221  Query_1615341     3LVP_A   31.741             293        160       11
    222  Query_1615341     3BKB_A   31.525             295        178       10
    223  Query_1615341     7EEF_A   31.752             274        167        9
    224  Query_1615341     5X5O_A   31.618             272        155        9
    225  Query_1615341     3LW0_A   31.399             293        161       11
    226  Query_1615341     5HES_A   31.618             272        155        9
    227  Query_1615341     1P4O_A   31.399             293        161       11
    228  Query_1615341     7EEC_A   31.429             280        172        9
    229  Query_1615341     2REI_A   31.752             274        167        9
    230  Query_1615341     7EED_A   31.429             280        172        9
    231  Query_1615341     4XLV_A   30.492             305        179       10
    232  Query_1615341     1RQQ_A   30.492             305        179       10
    233  Query_1615341     3CD3_A   31.058             293        182        8
    234  Query_1615341   8ATL_BBB   30.712             267        164        7
    235  Query_1615341     4UYA_A   32.632             285        142       14
    236  Query_1615341   8ATL_AAA   30.712             267        164        7
    237  Query_1615341     2Z8C_A   30.877             285        165        9
    238  Query_1615341   8ATB_AAA   30.712             267        164        7
    239  Query_1615341     8FLN_A   28.938             273        170        9
    240  Query_1615341     1GAG_A   30.877             285        165        9
    241  Query_1615341     7TYJ_A   31.724             290        164       10
    242  Query_1615341     6J6M_A   28.623             276        173        9
    243  Query_1615341     6MNY_A   28.938             273        170        9
    244  Query_1615341     3DK3_A   30.605             281        171       10
    245  Query_1615341     1I44_A   30.968             310        171       12
    246  Query_1615341     6VXQ_A   28.713             303        189       11
    247  Query_1615341     4IBM_A   30.744             309        173       12
    248  Query_1615341     6BL8_A   30.466             279        170       10
    249  Query_1615341     8DWN_A   30.492             305        179       10
    250  Query_1615341     5HU9_A   29.720             286        177       10
    251  Query_1615341     4Z3V_A   28.571             273        171        9
    252  Query_1615341     6NZM_A   28.571             273        171        9
    253  Query_1615341     4XEY_A   29.932             294        178       11
    254  Query_1615341     6W7O_A   28.571             273        171        9
    255  Query_1615341     4YHF_A   28.571             273        171        9
    256  Query_1615341     3GEN_A   28.571             273        171        9
    257  Query_1615341     8YVV_A   28.571             273        171        9
    258  Query_1615341     5J87_A   28.571             273        171        9
    259  Query_1615341     3P08_A   28.571             273        171        9
    260  Query_1615341     6TFP_A   28.571             273        171        9
    261  Query_1615341     5P9F_A   28.571             273        171        9
    262  Query_1615341     4OT5_A   28.767             292        182       10
    263  Query_1615341     3OCT_A   28.571             273        171        9
    264  Query_1615341     6S90_A   28.571             273        171        9
    265  Query_1615341     3PIX_A   28.571             273        171        9
    266  Query_1615341     4WA9_A   29.825             285        176       10
    267  Query_1615341     1IRK_A   30.744             309        173       12
    268  Query_1615341     6AUB_A   28.571             273        171        9
    269  Query_1615341     6O8I_A   28.571             273        171        9
    270  Query_1615341     5XYZ_A   28.571             273        171        9
    271  Query_1615341     6NFI_A   28.571             273        171        9
    272  Query_1615341     7KXQ_A   28.571             273        171        9
    273  Query_1615341     8GC8_A   28.571             273        171        9
    274  Query_1615341     5ZZ4_A   28.571             273        171        9
    275  Query_1615341     2XYN_A   30.961             281        170       10
    276  Query_1615341     1FPU_A   30.108             279        171       10
    277  Query_1615341     5FBN_C   28.571             273        171        9
    278  Query_1615341     4XLI_A   30.824             279        169       10
    279  Query_1615341     4ZLY_A   28.571             273        171        9
    280  Query_1615341     8FLL_A   28.571             273        171        9
    281  Query_1615341     9ZLJ_A   28.938             273        170        9
    282  Query_1615341     6XRG_A   30.466             279        170       10
    283  Query_1615341     2F4J_A   29.825             285        176       10
    284  Query_1615341     3OCS_A   28.571             273        171        9
    285  Query_1615341     3QRI_A   30.108             279        171       10
    286  Query_1615341     8SSN_A   29.932             294        178       11
    287  Query_1615341     2HZI_A   30.108             279        171       10
    288  Query_1615341     2NRY_A   30.712             267        164        7
    289  Query_1615341     4RX5_A   28.571             273        171        9
    290  Query_1615341     2E2B_A   30.108             279        171       10
    291  Query_1615341     4OTF_A   28.571             273        171        9
    292  Query_1615341     3OXZ_A   30.108             279        171       10
    293  Query_1615341     6O8U_A   30.712             267        164        7
    294  Query_1615341     2QOH_A   30.108             279        171       10
    295  Query_1615341     5K72_A   30.712             267        164        7
    296  Query_1615341     6EG9_A   30.712             267        164        7
    297  Query_1615341     6XE4_A   28.571             273        171        9
    298  Query_1615341     6XR6_A   30.108             279        171       10
    299  Query_1615341     4ZOG_A   30.108             279        171       10
    300  Query_1615341     4Y95_A   27.839             273        173        9
    301  Query_1615341     2NRU_A   30.712             267        164        7
    302  Query_1615341     6THW_A   30.712             267        164        7
    303  Query_1615341     2HIW_A   30.108             279        171       10
    304  Query_1615341     7N9G_A   30.108             279        171       10
    305  Query_1615341     6THX_A   30.712             267        164        7
    306  Query_1615341     6E4F_A   28.571             273        171        9
    307  Query_1615341     2HZ0_A   30.108             279        171       10
    308  Query_1615341     9PSU_A   30.712             267        164        7
    309  Query_1615341     6LXY_A   30.712             267        164        7
    310  Query_1615341     2HYY_A   30.108             279        171       10
    311  Query_1615341     2G1T_A   30.108             279        171       10
    312  Query_1615341     6BIK_A   28.571             273        171        9
    313  Query_1615341     1P14_A   30.421             309        174       12
    314  Query_1615341     2G2F_A   30.108             279        171       10
    315  Query_1615341     3EKK_A   30.293             307        177       12
    316  Query_1615341     8SCV_A   30.712             267        164        7
    317  Query_1615341     7C2W_A   31.298             262        159        7
    318  Query_1615341     8X2A_A   31.387             274        162       10
    319  Query_1615341     8H7F_A   29.720             286        177       10
    320  Query_1615341     3DK6_A   30.605             281        171       10
    321  Query_1615341     6F3D_A   31.298             262        159        7
    322  Query_1615341    9RPV_N1   31.618             272        155        9
    323  Query_1615341     3SXR_A   31.022             274        163       10
    324  Query_1615341     6F3I_A   30.712             267        164        7
    325  Query_1615341     6MOM_A   30.712             267        164        7
    326  Query_1615341     4RMZ_A   30.712             267        164        7
    327  Query_1615341     2OIB_A   30.712             267        164        7
    328  Query_1615341     5UIT_A   30.712             267        164        7
    329  Query_1615341     3PYY_A   30.108             279        171       10
    330  Query_1615341     4Y73_A   30.712             267        164        7
    331  Query_1615341     7CC2_A   30.108             279        171       10
    332  Query_1615341     7C2V_A   30.712             267        164        7
    333  Query_1615341     7QG3_A   30.712             267        164        7
    334  Query_1615341     9NA2_A   30.712             267        164        7
    335  Query_1615341     5W84_A   30.712             267        164        7
    336  Query_1615341   8BR5_AAA   30.712             267        164        7
    337  Query_1615341     6N8G_A   31.298             262        159        7
    338  Query_1615341     6I99_A   31.022             274        163       10
    339  Query_1615341     6F3E_A   31.298             262        159        7
    340  Query_1615341     6JK8_A   31.741             293        160       11
    341  Query_1615341     5E1S_A   30.421             309        174       12
    342  Query_1615341     1OPK_A   29.452             292        182       10
    343  Query_1615341     5HHW_A   30.421             309        174       12
    344  Query_1615341     6F3G_A   31.298             262        159        7
    345  Query_1615341     4Y93_A   28.205             273        172        9
    346  Query_1615341     3K54_A   28.571             273        171        9
    347  Query_1615341     5UIS_A   30.712             267        164        7
    348  Query_1615341     5UIQ_A   30.712             267        164        7
    349  Query_1615341     2GQG_A   30.108             279        171       10
    350  Query_1615341     6O94_A   30.712             267        164        7
    351  Query_1615341     8V1O_A   31.298             262        159        7
    352  Query_1615341     7SL1_A   30.392             306        178       11
    353  Query_1615341     7W7Y_A   30.108             279        171       10
    354  Query_1615341     1OPL_A   29.452             292        182       10
    355  Query_1615341     7W7X_A   30.108             279        171       10
    356  Query_1615341     8EYR_A   31.507             292        163       11
    357  Query_1615341     8DTL_A   30.392             306        178       11
    358  Query_1615341     3QRJ_A   29.749             279        172       10
    359  Query_1615341     2FO0_A   29.452             292        182       10
    360  Query_1615341     4YFF_A   31.227             269        166        8
    361  Query_1615341     4U97_A   30.337             267        165        7
    362  Query_1615341     4YFI_A   31.227             269        166        8
    363  Query_1615341     6AUA_A   28.205             273        172        9
    364  Query_1615341     7MGJ_A   31.227             269        166        8
    365  Query_1615341     6EGD_A   30.337             267        165        7
    366  Query_1615341     8S9F_A   28.938             273        170        9
    367  Query_1615341     8WTF_A   30.337             267        165        7
    368  Query_1615341     8EYX_A   30.392             306        178       11
    369  Query_1615341     4XI2_A   28.938             273        170        9
    370  Query_1615341     2Z60_A   29.749             279        172       10
    371  Query_1615341     3OY3_A   29.749             279        172       10
    372  Query_1615341     8DKS_A   31.298             262        159        7
    373  Query_1615341     4TWP_A   29.893             281        169       11
    374  Query_1615341     5MO4_A   29.110             292        183       10
    375  Query_1615341     6PYH_A   31.164             292        164       11
    376  Query_1615341     5BPY_A   28.309             272        171        9
    377  Query_1615341     3DK7_A   29.893             281        173       10
    378  Query_1615341     4NWM_A   28.309             272        171        9
    379  Query_1615341     2V7A_A   29.749             279        172       10
    380  Query_1615341     5AMN_A   31.469             286        166        8
    381  Query_1615341     7KHJ_A   31.419             296        163       10
    382  Query_1615341     7KHK_A   31.419             296        163       10
    383  Query_1615341     8FD9_A   27.574             272        173        9
    384  Query_1615341     1K2P_A   28.195             266        167        9
    385  Query_1615341     8GMB_A   28.571             273        171       10
    386  Query_1615341     6FEK_A   31.579             285        166        8
    387  Query_1615341     9KS5_A   28.623             276        179        8
    388  Query_1615341     7BW7_A   30.619             307        176       12
    389  Query_1615341     8U4B_A   30.619             307        176       12
    390  Query_1615341     7PG0_A   30.619             307        176       12
    391  Query_1615341     8VJB_A   30.619             307        176       12
    392  Query_1615341     2IVV_A   30.667             300        164        8
    393  Query_1615341     3ETA_A   30.769             286        158       11
    394  Query_1615341     8S93_A   28.309             272        171       10
    395  Query_1615341     6NE7_A   30.333             300        165        8
    396  Query_1615341     4HVS_A   31.081             296        164       10
    397  Query_1615341     6PXV_A   30.293             307        177       12
    398  Query_1615341     6NJA_A   30.333             300        165        8
    399  Query_1615341     4GU9_A   29.286             280        173        9
    400  Query_1615341     9SDI_A   29.286             280        173        9
    401  Query_1615341     2IVT_A   30.333             300        165        8
    402  Query_1615341     3PXK_A   29.286             280        173        9
    403  Query_1615341     6I8Z_A   29.286             280        173        9
    404  Query_1615341     2ETM_A   29.286             280        173        9
    405  Query_1615341     2IVS_A   30.333             300        165        8
    406  Query_1615341     2O8Y_A   29.963             267        166        7
    407  Query_1615341     6I83_A   30.333             300        165        8
    408  Query_1615341     1MP8_A   29.286             280        173        9
    409  Query_1615341     2R4B_A   30.943             265        166        7
    410  Query_1615341     2J0M_B   29.286             280        173        9
    411  Query_1615341     2JKM_A   29.286             280        173        9
    412  Query_1615341     3BBT_B   30.827             266        167        7
    413  Query_1615341     4EBV_A   29.286             280        173        9
    414  Query_1615341     6GQJ_A   31.186             295        163       10
    415  Query_1615341     8PQ9_A   31.186             295        163       10
    416  Query_1615341     1U54_A   30.483             269        163        9
    417  Query_1615341     3ZZW_A   28.720             289        164       12
    418  Query_1615341     4EWH_A   30.483             269        163        9
    419  Query_1615341     6GQK_A   31.186             295        163       10
    420  Query_1615341     6VQM_A   30.483             269        163        9
    421  Query_1615341     1U46_A   30.483             269        163        9
    422  Query_1615341     4XCU_A   30.201             298        166       12
    423  Query_1615341     7KP6_A   30.483             269        163        9
    424  Query_1615341     3BZ3_A   29.242             277        171        9
    425  Query_1615341     4GT4_A   28.720             289        164       12
    426  Query_1615341     5ZXB_A   30.483             269        163        9
    427  Query_1615341     2JKK_A   29.286             280        173        9
    428  Query_1615341     8PQG_A   31.525             295        162       10
    429  Query_1615341     8FE9_A   30.483             269        163        9
    430  Query_1615341     3EQP_A   30.483             269        163        9
    431  Query_1615341     2J0L_A   29.286             280        173        9
    432  Query_1615341     4ID7_A   30.483             269        163        9
    433  Query_1615341     4HZR_A   30.483             269        163        9
    434  Query_1615341     5FM2_A   30.132             302        163       10
    435  Query_1615341     4CKI_A   30.000             300        166        8
    436  Query_1615341     3QUP_A   28.912             294        176       10
    437  Query_1615341     6SDC_A   30.996             271        163       10
    438  Query_1615341     8ANS_A   30.996             271        163       10
    439  Query_1615341     9EF1_A   29.568             301        168       10
    440  Query_1615341     4HZS_A   30.483             269        163        9
    441  Query_1615341     7RUN_A   30.132             302        163       10
    442  Query_1615341     8AU5_A   30.996             271        163       10
    443  Query_1615341     6VHG_A   30.000             300        166        8
    444  Query_1615341     2J0J_A   27.479             353        225       11
    445  Query_1615341     6JPE_A   29.316             307        166       12
    446  Query_1615341     4UXQ_A   29.316             307        166       12
    447  Query_1615341     8KH9_A   29.316             307        166       12
    448  Query_1615341     7DTZ_A   29.316             307        166       12
    449  Query_1615341     5JKG_A   29.316             307        166       12
    450  Query_1615341     8W5C_A   29.316             307        166       12
    451  Query_1615341     7F3M_A   29.316             307        166       12
    452  Query_1615341     6YI8_A   29.316             307        166       12
    453  Query_1615341     4QQT_A   29.316             307        166       12
    454  Query_1615341     5NUD_A   29.316             307        166       12
    455  Query_1615341     5XFJ_A   29.316             307        166       12
    456  Query_1615341     3D7T_A   30.403             273        164       10
    457  Query_1615341     7VJL_A   29.316             307        166       12
    458  Query_1615341     4TYE_A   29.316             307        166       12
    459  Query_1615341     6I82_A   30.000             300        166        8
    460  Query_1615341     3T9T_A   27.143             280        178       10
    461  Query_1615341     1BYG_A   30.403             273        164       10
    462  Query_1615341     7WCW_A   29.316             307        166       12
    463  Query_1615341     5XFF_A   29.316             307        166       12
    464  Query_1615341     7WCX_A   29.316             307        166       12
    465  Query_1615341     6NFY_A   29.123             285        170        9
    466  Query_1615341     3D7U_A   30.403             273        164       10
    467  Query_1615341     4KNB_A   30.627             271        164       10
    468  Query_1615341     8AW1_A   30.627             271        164       10
    469  Query_1615341     4R1V_A   30.627             271        164       10
    470  Query_1615341   7PI4_DDD   28.986             276        171        9
    471  Query_1615341     4QQC_A   29.316             307        166       12
    472  Query_1615341     6SD9_A   30.627             271        164       10
    473  Query_1615341     8AN8_A   30.627             271        164       10
    474  Query_1615341     8VI1_A   30.627             271        164       10
    475  Query_1615341     6TY3_A   27.479             353        225       11
    476  Query_1615341     4QQJ_A   29.316             307        166       12
    477  Query_1615341     2WD1_A   30.627             271        164       10
    478  Query_1615341     2J0K_A   27.479             353        225       11
    479  Query_1615341     3GQI_A   28.231             294        174       10
    480  Query_1615341     7B3Q_A   30.627             271        164       10
    481  Query_1615341     6IUO_A   29.316             307        166       12
    482  Query_1615341     3F66_A   30.627             271        164       10
    483  Query_1615341     2G15_A   30.627             271        164       10
    484  Query_1615341     4EEV_A   30.627             271        164       10
    485  Query_1615341     3Q6U_A   30.627             271        164       10
    486  Query_1615341     8PAS_A   29.123             285        170        9
    487  Query_1615341     3LQ8_A   30.627             271        164       10
    488  Query_1615341     8FH4_A   28.826             281        176        9
    489  Query_1615341     6CQD_A   28.826             281        176        9
    490  Query_1615341     9C1R_A   30.627             271        164       10
    491  Query_1615341     2WGJ_A   30.627             271        164       10
    492  Query_1615341     8CDW_A   28.826             281        176        9
    493  Query_1615341     4QQ5_A   29.316             307        166       12
    494  Query_1615341     8PAU_A   28.669             293        177        9
    495  Query_1615341     6NG0_A   28.826             281        176        9
    496  Query_1615341     3VW8_A   30.627             271        164       10
    497  Query_1615341     8PAR_A   28.826             281        176        9
    498  Query_1615341     5A46_A   28.231             294        174       10
    499  Query_1615341     7M0L_A   28.374             289        183        9
    500  Query_1615341     9SZJ_A   30.627             271        164       10
    501  Query_1615341     4GG5_A   30.627             271        164       10
    502  Query_1615341     5UAB_A   30.627             271        164       10
    503  Query_1615341     5ZV2_A   28.231             294        174       10
    504  Query_1615341     3I5N_A   30.627             271        164       10
    505  Query_1615341     7M0M_A   28.374             289        183        9
    506  Query_1615341     8YKI_A   28.231             294        174       10
    507  Query_1615341     5A4C_A   28.231             294        174       10
    508  Query_1615341     3C4F_A   28.231             294        174       10
    509  Query_1615341     7M0K_A   28.374             289        183        9
    510  Query_1615341     4ZSA_A   28.231             294        174       10
    511  Query_1615341     3RHX_A   28.231             294        174       10
    512  Query_1615341     2RFN_A   30.627             271        164       10
    513  Query_1615341     9H8D_A   28.826             281        176        9
    514  Query_1615341     8XN7_A   28.826             281        176        9
    515  Query_1615341     4F63_A   28.231             294        174       10
    516  Query_1615341     4WUN_A   28.231             294        174       10
    517  Query_1615341     5FLF_A   28.231             294        174       10
    518  Query_1615341     5HLW_A   30.627             271        164       10
    519  Query_1615341     7L25_A   28.328             293        178        9
    520  Query_1615341     1AGW_A   28.231             294        174       10
    521  Query_1615341     6LVM_A   28.571             294        173       11
    522  Query_1615341     4F0F_A   28.269             283        169       11
    523  Query_1615341     3KXX_A   28.231             294        174       10
    524  Query_1615341     4PWN_A   33.460             263        147       10
    525  Query_1615341     5TQ5_A   29.787             282        161       11
    526  Query_1615341     3JS2_A   28.231             294        174       10
    527  Query_1615341     4HCT_A   26.855             283        181       10
    528  Query_1615341     7WCL_A   28.231             294        174       10
    529  Query_1615341     6CQF_A   28.772             285        171        9
    530  Query_1615341     5AM7_A   28.231             294        174       10
    531  Query_1615341     8Y22_A   28.231             294        174       10
    532  Query_1615341     9BI8_A   28.826             281        176        9
    533  Query_1615341     3GQL_A   28.231             294        174       10
    534  Query_1615341     6MZW_A   28.231             294        174       10
    535  Query_1615341     7L24_A   28.328             293        178        9
    536  Query_1615341     6HH1_A   31.399             293        163       10
    537  Query_1615341     3DKG_A   30.258             271        165       10
    538  Query_1615341     6NVL_A   28.231             294        174       10
    539  Query_1615341     8G6Z_A   29.537             281        163       10
    540  Query_1615341     3Q6W_A   30.627             271        164       10
    541  Query_1615341     6CQE_A   28.772             285        171        9
    542  Query_1615341     4RWI_A   28.231             294        174       10
    543  Query_1615341     7L26_A   28.328             293        178        9
    544  Query_1615341     6NFZ_A   28.826             281        176        9
    545  Query_1615341     7LL5_A   31.095             283        156       13
    546  Query_1615341     5J5T_A   29.151             271        173        7
    547  Query_1615341     7KAC_A   29.123             285        170        9
    548  Query_1615341     7LL4_A   31.095             283        156       13
    549  Query_1615341     3KRR_A   30.303             297        168       13
    550  Query_1615341     7RN6_A   31.095             283        156       13
    551  Query_1615341     5TQ4_A   29.787             282        161       11
    552  Query_1615341     8AU3_A   30.627             271        164       10
    553  Query_1615341     3QTI_A   30.258             271        165       10
    554  Query_1615341     3CE3_A   30.258             271        165       10
    555  Query_1615341     3C1X_A   30.258             271        165       10
    556  Query_1615341     4IWD_A   30.627             271        164       10
    557  Query_1615341     5GRN_A   28.527             319        175       10
    558  Query_1615341     3A4P_A   30.258             271        165       10
    559  Query_1615341     2P2H_A   28.912             294        170        9
    560  Query_1615341     8UDT_A   28.723             282        176        9
    561  Query_1615341     5WEV_A   29.874             318        172       14
    562  Query_1615341     3DKC_A   30.258             271        165       10
    563  Query_1615341     8PQJ_A   28.527             319        175       10
    564  Query_1615341     1R0P_A   30.258             271        165       10
    565  Query_1615341     4E4M_A   29.801             302        173       13
    566  Query_1615341     5TQ6_A   29.787             282        161       11
    567  Query_1615341     3MIY_A   27.338             278        176       10
    568  Query_1615341     5VND_A   28.231             294        174       10
    569  Query_1615341     5TQ3_A   29.787             282        161       11
    570  Query_1615341     5USY_A   31.095             283        156       13
    571  Query_1615341     9GB9_A   30.986             213        134        5
    572  Query_1615341     6VGL_A   31.095             283        156       13
    573  Query_1615341     6TU9_A   29.861             288        162       12
    574  Query_1615341     8UDV_A   28.723             282        176        9
    575  Query_1615341     4HGE_A   29.801             302        173       13
    576  Query_1615341     3UGC_A   30.303             297        168       13
    577  Query_1615341     8G8O_A   29.787             282        161       11
    578  Query_1615341     1K9A_A   30.403             273        164       10
    579  Query_1615341     3TJC_A   31.095             283        156       13
    580  Query_1615341     6AAJ_A   29.801             302        173       13
    581  Query_1615341     3C7Q_A   28.767             292        171        9
    582  Query_1615341     3TT0_A   28.231             294        174       10
    583  Query_1615341     7MO7_B   30.627             271        164       10
    584  Query_1615341     2OH4_A   28.669             293        171        9
    585  Query_1615341     9V70_A   32.178             202        123        6
    586  Query_1615341     8BX6_A   30.201             298        169       13
    587  Query_1615341     3E62_A   31.095             283        156       13
    588  Query_1615341     2XA4_A   31.095             283        156       13
    589  Query_1615341     4D0W_A   31.095             283        156       13
    590  Query_1615341     7Q7I_A   29.787             282        161       11
    591  Query_1615341     2B7A_A   31.095             283        156       13
    592  Query_1615341     4AQC_A   31.095             283        156       13
    593  Query_1615341     8BM2_A   30.201             298        169       13
    594  Query_1615341     3EWH_A   28.912             294        170        9
    595  Query_1615341     3U6J_A   28.912             294        170        9
    596  Query_1615341     3ZMM_A   31.095             283        156       13
    597  Query_1615341     9V71_A   32.178             202        123        6
    598  Query_1615341     4BBE_A   29.866             298        170       13
    599  Query_1615341     5HEZ_A   29.470             302        174       12
    600  Query_1615341     3CJF_A   28.966             290        171        9
    601  Query_1615341     1T46_A   30.976             297        163       10
    602  Query_1615341     4F1O_A   27.972             286        172       11
    603  Query_1615341     3RVG_A   29.766             299        171       13
    604  Query_1615341     6WTN_A   31.095             283        156       13
    605  Query_1615341     2W1I_A   31.095             283        156       13
    606  Query_1615341     7UYW_A   31.095             283        156       13
    607  Query_1615341     2OGV_A   30.241             291        171        9
    608  Query_1615341     2I0V_A   29.966             297        170        9
    609  Query_1615341     1YWN_A   28.571             294        171        9
    610  Query_1615341     4U0I_A   30.976             297        163       10
    611  Query_1615341     6MOB_A   30.976             297        163       10
    612  Query_1615341     6WXJ_A   29.966             297        170        9
    613  Query_1615341     6TPD_A   29.787             282        161       11
    614  Query_1615341     3Q32_A   31.095             283        156       13
    615  Query_1615341     4YTC_A   31.095             283        156       13
    616  Query_1615341     4F1M_A   27.915             283        170       11
    617  Query_1615341     9D3F_A   33.840             263        145       11
    618  Query_1615341     2X7F_A   29.964             277        163       13
    619  Query_1615341     8W1L_A   29.966             297        170        9
    620  Query_1615341     7TNH_A   29.966             297        170        9
    621  Query_1615341     3JY9_A   30.000             300        171       13
    622  Query_1615341     5HOA_A   29.889             271        166       10
    623  Query_1615341     3CJG_A   28.966             290        171        9
    624  Query_1615341     7SIU_A   28.772             285        171        9
    625  Query_1615341     5UGL_A   27.211             294        177       10
    626  Query_1615341     3IO7_A   30.000             300        171       13
    627  Query_1615341     7XZQ_A   29.964             277        163       13
    628  Query_1615341     7ZVS_A   29.787             282        163       14
    629  Query_1615341     7FCZ_A   28.253             269        164        9
    630  Query_1615341     5AX9_A   29.964             277        163       13
    631  Query_1615341     6RA5_A   29.964             277        163       13
    632  Query_1615341     8ZML_A   29.964             277        163       13
    633  Query_1615341     3V5J_A   26.978             278        177       10
    634  Query_1615341     6C4D_A   28.253             269        164        9
    635  Query_1615341     9HY8_A   28.253             269        164        9
    636  Query_1615341     1VR2_A   28.571             294        171        9
    637  Query_1615341     7UOS_A   33.840             263        145       11
    638  Query_1615341     1SM2_A   26.978             278        177       10
    639  Query_1615341     4ZIM_A   31.095             283        156       13
    640  Query_1615341     6NYH_A   28.253             269        164        9
    641  Query_1615341     9WYR_A   31.683             202        124        6
    642  Query_1615341     9KFU_A   28.231             294        174       11
    643  Query_1615341     6OL2_A   33.840             263        145       11
    644  Query_1615341     8X88_A   29.964             277        163       13
    645  Query_1615341     4ITH_A   28.253             269        164        9
    646  Query_1615341     5DRB_A   33.840             263        145       11
    647  Query_1615341     2P2I_A   28.571             294        171        9
    648  Query_1615341     3VNT_A   28.571             294        171        9
    649  Query_1615341     3PLS_A   29.044             272        169        9
    650  Query_1615341     9CD7_A   28.716             296        170       13
    651  Query_1615341     8W3D_A   27.211             294        177       10
    652  Query_1615341     9MZZ_A   28.358             268        163        9
    653  Query_1615341     8W3B_A   27.211             294        177       10
    654  Query_1615341     9MZY_A   28.358             268        163        9
    655  Query_1615341     9MZX_A   28.358             268        163        9
    656  Query_1615341     8W2X_A   27.211             294        177       10
    657  Query_1615341     8JQI_B   28.231             294        174       10
    658  Query_1615341     4K33_A   28.716             296        170       13
    659  Query_1615341     1PKG_A   30.976             297        163       10
    660  Query_1615341     8W38_A   27.211             294        177       10
    661  Query_1615341     4Q2A_A   33.840             263        145       11
    662  Query_1615341     7TEU_A   29.766             299        171       13
    663  Query_1615341     2XIR_A   28.571             294        171        9
    664  Query_1615341     3LCD_A   30.241             291        171        9
    665  Query_1615341     9GZH_A   28.213             319        176       10
    666  Query_1615341     8PQH_A   28.213             319        176       10
    667  Query_1615341     5W7T_A   33.840             263        145       11
    668  Query_1615341     5HX6_A   28.253             269        164        9
    669  Query_1615341     4APC_A   28.148             270        171        7
    670  Query_1615341     7ZW8_A   31.313             297        162       11
    671  Query_1615341     5HOR_A   29.889             271        166       10
    672  Query_1615341     6ITV_A   30.976             297        163       10
    673  Query_1615341     6A32_A   28.213             319        176       10
    674  Query_1615341     6C3E_A   28.253             269        164        9
    675  Query_1615341     6GQO_A   28.571             294        171        9
    676  Query_1615341     4E6D_A   29.431             299        172       13
    677  Query_1615341     3WZD_A   28.571             294        171        9
    678  Query_1615341     6PNX_A   28.231             294        174       11
    679  Query_1615341     6NW2_A   28.253             269        164        9
    680  Query_1615341     5UGX_A   27.211             294        177       10
    681  Query_1615341     3G0E_A   30.976             297        163       10
    682  Query_1615341     9GTG_A   28.253             269        164        9
    683  Query_1615341     1T45_A   30.976             297        163       10
    684  Query_1615341     6JOI_A   28.213             319        176       10
    685  Query_1615341     3G0F_A   30.976             297        163       10
    686  Query_1615341     4J98_A   27.211             294        177       10
    687  Query_1615341     1FAQ_A   73.077              52         14        0
    688  Query_1615341     2PVF_A   26.871             294        178       10
    689  Query_1615341     6FD3_A   28.794             257        165        9
    690  Query_1615341     6AGX_A   26.871             294        178       10
    691  Query_1615341     8XRR_A   28.213             319        176       10
    692  Query_1615341     3CLY_A   26.871             294        178       10
    693  Query_1615341     5TF9_A   33.080             263        147       10
    694  Query_1615341     4USF_A   27.239             268        179        9
    695  Query_1615341     4AGC_A   28.571             294        171        9
    696  Query_1615341     3BEA_A   29.392             296        171       10
    697  Query_1615341     8BEM_A   27.239             268        179        9
    698  Query_1615341     4GL9_A   28.053             303        183       10
    699  Query_1615341     2I1M_A   30.068             296        169       12
    700  Query_1615341     2PZ5_A   26.871             294        178       10
    701  Query_1615341     8E1X_A   26.871             294        178       10
    702  Query_1615341     6BDN_A   30.078             256        163        7
    703  Query_1615341     2J51_A   27.239             268        179        9
    704  Query_1615341     2JFM_A   27.239             268        179        9
    705  Query_1615341     2JFL_A   27.239             268        179        9
    706  Query_1615341     2PZP_A   26.871             294        178       10
    707  Query_1615341     4J97_A   26.871             294        178       10
    708  Query_1615341     4KIO_A   26.978             278        177       10
    709  Query_1615341     8H75_A   26.531             294        179        9
    710  Query_1615341     4NEU_A   28.253             269        164        9
    711  Query_1615341     2PWL_A   26.871             294        178       10
    712  Query_1615341     6LVK_A   26.531             294        179        9
    713  Query_1615341     9U7E_A   26.531             294        179        9
    714  Query_1615341     8SWE_A   26.871             294        178       10
    715  Query_1615341     1GJO_A   26.871             294        178       10
    716  Query_1615341     3RI1_A   26.871             294        178       10
    717  Query_1615341     5EG3_A   26.871             294        178       10
    718  Query_1615341     3B2T_A   26.871             294        178       10
    719  Query_1615341     8STG_A   26.871             294        178       10
    720  Query_1615341     2PVY_A   26.871             294        178       10
    721  Query_1615341     6LVL_A   26.531             294        179        9
    722  Query_1615341     9GFZ_A   32.836             201        121        6
    723  Query_1615341     3QGW_A   26.619             278        178       10
    724  Query_1615341     9LBG_A   26.045             311        206       11
    725  Query_1615341   7OZY_AAA   26.531             294        179        9
    726  Query_1615341     8E4T_A   28.417             278        173       10
    727  Query_1615341     4J96_A   26.871             294        178       10
    728  Query_1615341     4J99_A   26.871             294        178       10
    729  Query_1615341     2PSQ_A   26.871             294        178       10
    730  Query_1615341     9U3N_A   26.871             294        178       10
    731  Query_1615341     5UHN_A   26.871             294        178       10
    732  Query_1615341     7KIA_A   26.531             294        179        9
    733  Query_1615341     2Q0B_A   26.871             294        178       10
    734  Query_1615341     8JOT_A   28.025             314        172       10
    735  Query_1615341     2PZR_A   26.531             294        179       10
    736  Query_1615341     9LBF_A   26.045             311        206       11
    737  Query_1615341     2PY3_A   26.871             294        178       10
    738  Query_1615341     3LCO_A   28.947             304        171        9
    739  Query_1615341     5UI0_A   26.871             294        178       10
    740  Query_1615341     6T2W_A   28.947             304        171        9
    741  Query_1615341     9D51_A   27.237             257        169        9
    742  Query_1615341     4HW7_A   28.383             303        174        9
    743  Query_1615341     7AAY_A   25.000             296        190       10
    744  Query_1615341     8CGC_A   28.947             304        171        9
    745  Query_1615341     4O27_B   27.652             264        172        9
    746  Query_1615341     9BHI_A   24.671             304        193       10
    747  Query_1615341     5U6B_A   25.784             287        183        8
    748  Query_1615341     9D7Q_A   32.558             258        141       12
    749  Query_1615341     9M41_A   25.723             311        207       11
    750  Query_1615341     7DXL_A   24.667             300        190       10
    751  Query_1615341     6V6Q_A   26.531             294        179       10
    752  Query_1615341     2GCD_A   29.070             258        167        7
    753  Query_1615341     3S95_A   28.520             277        176        7
    754  Query_1615341     1U5Q_A   29.070             258        167        7
    755  Query_1615341     4QML_A   27.652             264        172        9
    756  Query_1615341     7AAZ_A   24.667             300        190       10
    757  Query_1615341     7OAM_A   24.667             300        190       10
    758  Query_1615341     3A7F_A   27.652             264        172        9
    759  Query_1615341     7CQE_A   24.667             300        190       10
    760  Query_1615341     4W8E_A   27.652             264        172        9
    761  Query_1615341     5TC0_A   24.667             300        190       10
    762  Query_1615341     4U8Z_A   27.652             264        172        9
    763  Query_1615341     7B30_A   27.652             264        172        9
    764  Query_1615341     2P0C_A   24.667             300        190       10
    765  Query_1615341     5O1V_A   32.558             258        141       12
    766  Query_1615341     7Z5W_A   30.124             322        177       15
    767  Query_1615341     5U6C_A   24.667             300        190       10
    768  Query_1615341     3ZHP_C   27.652             264        172        9
    769  Query_1615341     5O2B_A   32.558             258        141       12
    770  Query_1615341     8SE1_A   25.541             462        290       18
    772  Query_1615341     8SE2_A   25.541             462        290       18
    774  Query_1615341     2XIK_A   29.658             263        168        9
    775  Query_1615341     3CKW_A   27.652             264        172        9
    776  Query_1615341     1LUF_A   29.630             297        159       11
    777  Query_1615341     3CKX_A   27.652             264        172        9
    778  Query_1615341     7Z4V_A   29.658             263        168        9
    779  Query_1615341     4V0G_B   28.470             281        167        9
    780  Query_1615341     3PFQ_A   25.054             463        291       18
    782  Query_1615341     4RIO_A   28.470             281        167        9
    783  Query_1615341     5CNN_A   29.057             265        171        7
    784  Query_1615341     8HV2_A   29.057             265        171        7
    785  Query_1615341     5HVJ_A   28.159             277        177        7
    786  Query_1615341     3LXN_A   29.348             276        156       10
    787  Query_1615341     4ZJV_A   29.057             265        171        7
    788  Query_1615341     5HVK_A   28.159             277        177        7
    789  Query_1615341     4GVJ_A   29.348             276        156       10
    790  Query_1615341     8EDH_A   33.074             257        141       12
    791  Query_1615341     7SYD_A   28.947             266        170        7
    792  Query_1615341     5O26_A   32.296             257        143       11
    793  Query_1615341     4V0G_A   28.470             281        167        9
    794  Query_1615341     5CAV_A   29.057             265        171        7
    795  Query_1615341     1M14_A   29.057             265        171        7
    796  Query_1615341     7UKV_A   29.057             265        171        7
    797  Query_1615341     1YVJ_A   28.114             281        168        9
    798  Query_1615341     3D5V_A   30.667             225        139        7
    799  Query_1615341     4HVD_A   28.315             279        170        8
    800  Query_1615341     5O2C_A   32.558             258        141       12
    801  Query_1615341     4TKS_A   29.057             265        171        7
    802  Query_1615341     4Z16_A   28.114             281        168        9
    803  Query_1615341     4I23_A   29.057             265        171        7
    804  Query_1615341     3D5U_A   30.667             225        139        7
    805  Query_1615341     3VJO_A   29.057             265        171        7
    806  Query_1615341     2RFD_A   29.057             265        171        7
    807  Query_1615341     2ITW_A   29.057             265        171        7
    808  Query_1615341     9FZR_A   29.057             265        171        7
    809  Query_1615341     3D5W_A   30.667             225        139        7
    810  Query_1615341     4G5J_A   29.057             265        171        7
    811  Query_1615341     4JQ7_A   29.057             265        171        7
    812  Query_1615341     4LI5_A   29.057             265        171        7
    813  Query_1615341     4RIW_B   29.057             265        171        7
    814  Query_1615341     9H42_A   29.057             265        171        7
    815  Query_1615341     2GS2_A   29.057             265        171        7
    816  Query_1615341     5FED_A   29.057             265        171        7
    817  Query_1615341     9UBI_A   26.568             271        179        8
    818  Query_1615341     4WKQ_A   29.057             265        171        7
    819  Query_1615341     3D5X_A   30.667             225        139        7
    820  Query_1615341     6JZ0_A   29.057             265        171        7
    821  Query_1615341     6CN9_A   33.716             261        148       11
    822  Query_1615341     7UYV_A   28.470             281        167        9
    823  Query_1615341   7Q6H_AAA   28.470             281        167        9
    824  Query_1615341     2J5E_A   29.057             265        171        7
    825  Query_1615341     8A27_A   29.057             265        171        7
    826  Query_1615341     7AEM_A   29.057             265        171        7
    827  Query_1615341     8PO3_A   29.057             265        171        7
    828  Query_1615341     2RGP_A   29.057             265        171        7
    829  Query_1615341     8PO4_A   29.057             265        171        7
    830  Query_1615341     1XKK_A   29.057             265        171        7
    831  Query_1615341     3ZC6_A   28.470             281        167        9
    832  Query_1615341     4PY1_A   28.986             276        157       10
    833  Query_1615341     5W86_A   28.470             281        167        9
    834  Query_1615341     8H7X_A   29.057             265        171        7
    835  Query_1615341     9N6G_A   29.057             265        171        7
    836  Query_1615341     5NG0_A   27.241             290        192       10
    837  Query_1615341     5XGN_A   29.057             265        171        7
    838  Query_1615341     9BY4_A   29.057             265        171        7
    839  Query_1615341     7XDY_A   26.568             271        179        8
    840  Query_1615341     7XDV_A   26.568             271        179        8
    841  Query_1615341     7XDX_A   26.568             271        179        8
    842  Query_1615341     2GS7_A   29.057             265        171        7
    843  Query_1615341     5TOZ_A   28.315             279        170        8
    844  Query_1615341     6AAM_A   28.986             276        157       10
    845  Query_1615341     5W5J_A   27.241             290        192       10
    846  Query_1615341     3LXK_A   28.315             279        170        8
    847  Query_1615341     9QBG_A   29.434             265        170        6
    848  Query_1615341     9QBF_A   29.434             265        170        6
    849  Query_1615341     6HMX_A   27.437             277        185        9
    850  Query_1615341     4HJO_A   29.057             265        171        7
    851  Query_1615341     4TPT_A   30.986             284        165       11
    852  Query_1615341     8X2O_A   27.241             290        192       10
    853  Query_1615341     9F3V_A   27.241             290        192       10
    854  Query_1615341     6ES0_A   27.891             294        185       11
    855  Query_1615341     6UL8_A   27.241             290        192       10
    856  Query_1615341     5NXD_A   30.986             284        165       11
    857  Query_1615341   9QXN_AAA   29.057             265        171        7
    858  Query_1615341     5AR2_A   27.437             277        185        9
    859  Query_1615341     8AZA_A   27.891             294        185       11
    860  Query_1615341     5NG3_B   27.076             277        186        9
    861  Query_1615341     6HZV_A   28.470             281        167        9
    862  Query_1615341     3LZB_A   29.057             265        171        7
    863  Query_1615341     4C8B_A   28.114             281        178       10
    864  Query_1615341     5ZWJ_A   29.057             265        171        7
    865  Query_1615341     6LUB_A   28.679             265        172        7
    866  Query_1615341     8HV4_A   28.679             265        172        7
    867  Query_1615341     2V5Q_A   32.673             202        119        7
    868  Query_1615341     4E4L_A   25.623             281        178        9
    869  Query_1615341     6N7A_A   25.623             281        178        9
    870  Query_1615341     3IKA_A   28.679             265        172        7
    871  Query_1615341     4GIH_A   28.986             276        157       10
    872  Query_1615341     6GTT_A   26.354             277        176       10
    873  Query_1615341     8WD4_A   28.679             265        172        7
    874  Query_1615341     9D3V_A   28.679             265        172        7
    875  Query_1615341     5TD2_A   24.579             297        188       10
    876  Query_1615341     6ELR_A   25.623             281        178        9
    877  Query_1615341     7TVD_A   28.352             261        175        6
    878  Query_1615341     6GGH_A   25.623             281        178        9
    879  Query_1615341     3PJC_A   28.114             281        168        9
    880  Query_1615341     2YAC_A   32.673             202        119        7
    881  Query_1615341     9S3X_A   28.679             265        172        7
    882  Query_1615341     3PP0_A   28.679             265        172        6
    883  Query_1615341     7PCD_A   28.679             265        172        6
    884  Query_1615341     2QKW_B   27.372             274        177        7
    885  Query_1615341     9QEK_A   30.097             309        170       14
    886  Query_1615341     3ZBF_A   30.097             309        170       14
    887  Query_1615341     7JXH_A   28.679             265        172        6
    888  Query_1615341     3KB7_A   32.673             202        119        7
    889  Query_1615341     2JIU_A   28.679             265        172        7
    890  Query_1615341     4I24_A   28.679             265        172        7
    891  Query_1615341     5KHW_A   25.623             281        178        9
    892  Query_1615341     7SZ0_A   28.571             266        171        7
    893  Query_1615341     5GMP_A   28.679             265        172        7
    894  Query_1615341     5Y9T_A   28.679             265        172        7
    895  Query_1615341     5FEE_A   28.679             265        172        7
    896  Query_1615341     3EYG_A   25.623             281        178        9
    897  Query_1615341     4LQM_A   28.679             265        172        7
    898  Query_1615341     3HGK_A   27.372             274        177        7
    899  Query_1615341     6TPE_A   25.623             281        178        9
    900  Query_1615341     5GNK_A   28.679             265        172        7
    901  Query_1615341     5AJQ_A   28.319             226        137        8
    902  Query_1615341     4G5P_A   28.679             265        172        7
    903  Query_1615341     2J7T_A   28.319             226        137        8
    904  Query_1615341     2JIT_A   28.679             265        172        7
    905  Query_1615341     6EIM_A   28.319             226        137        8
    906  Query_1615341     4BC6_A   28.319             226        137        8
    907  Query_1615341     5XDL_A   28.679             265        172        7
    908  Query_1615341     4YNE_A   27.619             315        173       12
    909  Query_1615341     2OU7_A   32.673             202        119        7
    910  Query_1615341     4NZW_B   29.278             263        169        9
    911  Query_1615341     5J9Z_A   28.679             265        172        7
    912  Query_1615341     5J9Y_A   28.679             265        172        7
    913  Query_1615341     4ZSE_A   28.679             265        172        7
    914  Query_1615341     5NG3_A   27.240             279        183        9
    915  Query_1615341     2EB2_A   28.679             265        172        7
    916  Query_1615341     4QPS_A   28.114             281        168        9
    917  Query_1615341     2ITN_A   28.679             265        172        7
    918  Query_1615341     6JWL_A   28.679             265        172        7
    919  Query_1615341     3THB_A   32.673             202        119        7
    920  Query_1615341     3GGF_A   27.757             263        173        9
    921  Query_1615341     6HXF_A   27.602             221        145        7
    922  Query_1615341     7C3N_A   27.957             279        171        8
    923  Query_1615341     4J52_A   32.353             204        121        7
    924  Query_1615341     2EB3_A   28.679             265        172        7
    925  Query_1615341     6P1D_A   28.679             265        172        7
    926  Query_1615341     2ITT_A   28.679             265        172        7
    927  Query_1615341     4H1J_A   26.022             269        174        8
    928  Query_1615341     9JQ1_A   28.679             265        172        7
    929  Query_1615341     6V5N_A   28.679             265        172        7
    930  Query_1615341     2RKU_A   32.353             204        121        7
    931  Query_1615341     6C9D_A   29.524             210        138        5
    932  Query_1615341     6TFU_A   28.679             265        172        7
    933  Query_1615341     8PEH_A   31.683             202        124        6
    934  Query_1615341     3ET7_A   26.022             269        174        8
    935  Query_1615341     4R3P_A   28.679             265        172        7
    936  Query_1615341     5TO8_A   26.022             269        174        8
    937  Query_1615341     5LWM_A   28.114             281        168        9
    938  Query_1615341     8S9P_C   29.153             295        163       11
    939  Query_1615341     7K1H_A   28.679             265        172        7
    940  Query_1615341     9P9U_A   29.057             265        171        7
    941  Query_1615341     9P9U_A   29.057             265        171        7
    942  Query_1615341     9XU9_A   28.679             265        172        7
    943  Query_1615341     3CC6_A   26.022             269        174        8
    944  Query_1615341     4I20_A   28.679             265        172        7
    945  Query_1615341     4RJ4_A   28.302             265        173        7
    946  Query_1615341     5JFS_A   27.445             317        175       12
    947  Query_1615341     7UYR_A   28.986             276        157       10
    948  Query_1615341     7VKO_A   27.619             315        173       12
    949  Query_1615341     6D22_A   27.445             317        175       12
    950  Query_1615341     8HY7_A   28.302             265        173        7
    951  Query_1615341     4RIY_A   27.407             270        177        7
    952  Query_1615341     7VKM_A   27.619             315        173       12
    953  Query_1615341     5TA6_A   32.673             202        119        7
    954  Query_1615341     3GOP_A   28.679             265        172        7
    955  Query_1615341     9DF3_A   29.104             268        170        8
    956  Query_1615341     9DF2_A   29.104             268        170        8
    957  Query_1615341     4KS7_A   31.658             199        123        7
    958  Query_1615341     6C7Y_A   25.806             279        176        9
    959  Query_1615341     6S9B_A   28.302             265        173        7
    960  Query_1615341     6S9C_A   28.302             265        173        7
    961  Query_1615341     4J7B_A   31.401             207        125        7
    962  Query_1615341     5H3Q_A   27.832             309        168       12
    963  Query_1615341     7MN5_B   28.679             265        172        6
    964  Query_1615341     4AOJ_A   27.832             309        168       12
    965  Query_1615341     7MN6_B   28.679             265        172        6
    966  Query_1615341     8A2B_A   28.679             265        172        7
    967  Query_1615341     5KMI_A   27.832             309        168       12
    968  Query_1615341     3DAK_A   26.950             282        181        9
    969  Query_1615341     5KML_A   27.832             309        168       12
    970  Query_1615341     2HAK_A   29.524             210        138        5
    971  Query_1615341     3NZ0_A   28.623             276        158       10
    972  Query_1615341     2VWI_A   26.950             282        181        9
    973  Query_1615341     2C30_A   31.658             199        123        7
    974  Query_1615341     4ZP5_A   28.571             273        172        9
    975  Query_1615341     8DSW_A   29.104             268        171        8
    976  Query_1615341     4OBO_A   28.571             273        172        9
    977  Query_1615341     8J5W_A   27.832             309        168       12
    978  Query_1615341     5DI1_A   28.571             273        172        9
    979  Query_1615341     5J95_A   28.571             273        172        9
    980  Query_1615341     4F0I_A   27.832             309        168       12
    981  Query_1615341     2QNJ_A   27.600             250        157        8
    982  Query_1615341     4U3Z_A   28.571             273        172        9
    983  Query_1615341     4LRM_A   30.370             270        164        9
    984  Query_1615341     6NSP_A   27.832             309        168       12
    985  Query_1615341     6D1Y_A   27.832             309        168       12
    986  Query_1615341     7XAF_A   27.832             309        168       12
    987  Query_1615341     6IQN_A   27.832             309        168       12
    988  Query_1615341     2A19_B   28.421             285        167        8
    989  Query_1615341     4RVT_A   28.571             273        172        9
    990  Query_1615341     4E1Z_A   27.857             280        163       10
    991  Query_1615341     4E20_A   27.857             280        163       10
    992  Query_1615341     4GT5_A   27.832             309        168       12
    993  Query_1615341     7MN5_A   27.407             270        177        7
    994  Query_1615341     8J5X_A   27.653             311        170       12
    995  Query_1615341     4LL0_A   28.302             265        173        7
    996  Query_1615341     8PO0_A   30.370             270        164        9
    997  Query_1615341     8UOI_A   27.200             250        158        8
    998  Query_1615341     7AAX_A   24.667             300        190       10
    999  Query_1615341     4I5M_A   26.459             257        172        6
    1000 Query_1615341     9U8C_A   30.370             270        164        9
    1001 Query_1615341     7M5Z_A   24.667             300        190       10
    1002 Query_1615341     6NSS_A   27.653             311        170       12
    1003 Query_1615341     4I21_A   28.302             265        173        7
    1004 Query_1615341     8V5I_A   28.571             273        172        9
    1005 Query_1615341     3LMG_A   27.037             270        178        7
    1006 Query_1615341     3UG1_A   28.302             265        173        7
    1007 Query_1615341     7P1L_A   28.455             246        160        9
    1008 Query_1615341     4RIX_A   27.037             270        178        7
    1009 Query_1615341     4RIW_A   27.037             270        178        7
    1010 Query_1615341     7FEH_A   25.697             323        180       12
    1011 Query_1615341     8PO1_A   30.370             270        164        9
    1012 Query_1615341     5F1Z_A   28.261             276        159       10
    1013 Query_1615341     6YAT_A   26.724             232        162        5
    1014 Query_1615341     3NYX_A   28.261             276        159       10
    1015 Query_1615341     9KLW_A   28.302             265        173        7
    1016 Query_1615341     6OP9_A   27.037             270        178        7
    1017 Query_1615341     5SAU_A   25.697             323        180       12
    1018 Query_1615341     7MX3_A   28.519             270        156       10
    1019 Query_1615341   8C12_AAA   30.392             204        123        7
    1020 Query_1615341     8UOH_A   27.200             250        158        8
    1021 Query_1615341     3FE3_A   27.200             250        158        8
    1022 Query_1615341     3COM_A   26.724             232        162        5
    1023 Query_1615341     6CTH_A   28.788             198        132        5
    1024 Query_1615341     7OXB_A   28.302             265        173        7
    1025 Query_1615341     8UOJ_A   27.200             250        158        8
    1026 Query_1615341   9R1W_AAA   32.178             202        120        7
    1027 Query_1615341     3W2O_A   28.302             265        173        7
    1028 Query_1615341     6M0U_A   31.000             200        128        4
    1029 Query_1615341     3KEX_A   27.037             270        178        7
    1030 Query_1615341     8PAV_A   26.724             232        162        5
    1031 Query_1615341     7LGS_A   30.370             270        164        9
    1032 Query_1615341     4I1Z_A   28.302             265        173        7
    1033 Query_1615341     5Y25_A   28.302             265        173        7
    1034 Query_1615341     5WNI_A   28.205             273        182        8
    1035 Query_1615341     8KFQ_A   28.302             265        173        7
    1036 Query_1615341     9FQP_A   30.370             270        164        9
    1037 Query_1615341     7MON_B   28.519             270        156       10
    1038 Query_1615341     2F57_A   30.392             204        123        7
    1039 Query_1615341     8U8X_A   29.520             271        166        8
    1040 Query_1615341     4FZA_B   27.376             263        174        9
    1041 Query_1615341     8VB5_A   29.520             271        166        8
    1042 Query_1615341     9IIC_A   27.004             237        155        6
    1043 Query_1615341     5WR7_A   27.331             311        171       12
    1044 Query_1615341     8D73_A   28.302             265        173        7
    1045 Query_1615341     6NPT_A   27.331             311        171       12
    1046 Query_1615341     4PMM_A   28.155             309        167       13
    1047 Query_1615341     9LFU_A   28.309             272        158       10
    1048 Query_1615341     6D3K_A   27.797             295        163        8
    1049 Query_1615341     2WZJ_A   28.634             227        150        6
    1050 Query_1615341     6PL1_A   27.508             309        169       12
    1051 Query_1615341     4B6L_A   27.907             258        167        8
    1052 Query_1615341     2MSE_D   60.563              71         28        0
    1053 Query_1615341     2R0I_A   29.075             227        149        6
    1054 Query_1615341     4OLI_A   27.737             274        163        9
    1055 Query_1615341     1ZMU_A   29.075             227        149        6
    1056 Query_1615341     5EAK_A   29.439             214        139        6
    1057 Query_1615341     1WXM_A   60.563              71         28        0
    1058 Query_1615341     1RFA_A   54.667              75         31        1
    1059 Query_1615341     6JRK_A   29.057             265        171        7
    1060 Query_1615341     3UIU_A   27.797             295        163        8
    1061 Query_1615341     5DBX_A   26.596             282        182        8
    1062 Query_1615341     8A5J_A   27.830             212        145        5
    1063 Query_1615341     8QEL_B   28.889             270        155        8
    1064 Query_1615341     5BVK_A   25.228             329        180       12
    1065 Query_1615341     3ZOS_A   25.228             329        180       12
    1066 Query_1615341     5D9H_A   26.596             282        182        8
    1067 Query_1615341     1ZMW_A   28.634             227        150        6
    1068 Query_1615341     6S89_A   29.057             265        171        7
    1069 Query_1615341     1ZMV_A   28.634             227        150        6
    1070 Query_1615341     1C1Y_B   55.405              74         30        1
    1071 Query_1615341     4G0N_B   55.405              74         30        1
    1072 Query_1615341     5KZ7_A   29.439             214        139        6
    1073 Query_1615341     6JRJ_A   29.057             265        171        7
    1074 Query_1615341     6VJJ_B   55.405              74         30        1
    1075 Query_1615341     6Y23_A   25.228             329        180       12
    1076 Query_1615341     6BRJ_A   25.228             329        180       12
    1077 Query_1615341     1GUA_B   55.405              74         30        1
    1078 Query_1615341     5E8U_A   30.216             278        148       13
    1079 Query_1615341     6VC0_A   25.758             264        173       10
    1080 Query_1615341     5FDP_A   25.228             329        180       12
    1081 Query_1615341     8JOF_A   53.659              82         33        2
    1082 Query_1615341     4ASZ_A   27.609             297        175       12
    1083 Query_1615341     8TXY_A   28.505             214        141        6
    1084 Query_1615341     1RW8_A   30.686             277        148       13
    1085 Query_1615341     1B6C_B   30.686             277        148       13
    1086 Query_1615341     5USQ_A   30.686             277        148       13
    1087 Query_1615341     9F6X_A   30.686             277        148       13
    1088 Query_1615341     8YHF_A   30.686             277        148       13
    1089 Query_1615341     9J9D_A   30.686             277        148       13
    1090 Query_1615341     1VJY_A   30.686             277        148       13
    1091 Query_1615341     5DH3_A   27.897             233        154        7
    1092 Query_1615341     8A66_B   28.507             221        144        7
    1093 Query_1615341     3TZM_A   30.686             277        148       13
    1094 Query_1615341     1PY5_A   30.686             277        148       13
    1095 Query_1615341     4X0M_A   30.686             277        148       13
    1096 Query_1615341     5FRI_A   30.686             277        148       13
    1097 Query_1615341     8A66_A   28.761             226        137        8
    1098 Query_1615341     2WOT_A   30.686             277        148       13
    1099 Query_1615341     5I8A_A   27.986             293        156       12
    1100 Query_1615341     5KVT_A   27.986             293        156       12
    1101 Query_1615341     5E8T_A   30.686             277        148       13
    1102 Query_1615341     4YNZ_A   28.295             258        164        9
    1103 Query_1615341     4YOM_B   28.295             258        164        9
    1104 Query_1615341     8XFL_A   26.531             245        166        8
    1105 Query_1615341     1RRB_A   54.667              75         31        1
    1106 Query_1615341     5ES1_A   26.531             245        166        8
    1107 Query_1615341     3IEC_A   29.954             217        134        7
    1108 Query_1615341     5LPZ_A   28.239             301        184        9
    1109 Query_1615341     4M68_A   25.177             282        175       11
    1110 Query_1615341     6FER_A   26.174             298        173       11
    1111 Query_1615341     3KUD_B   54.054              74         31        1
    1112 Query_1615341     1Y8G_A   28.634             227        150        6
    1113 Query_1615341     3KUC_B   54.054              74         31        1
    1114 Query_1615341     3V5Q_A   26.689             296        172       10
    1115 Query_1615341     4BFM_A   31.220             205        130        8
    1116 Query_1615341     5LPV_A   28.309             272        166        7
    1117 Query_1615341     5LPB_A   28.309             272        166        7
    1118 Query_1615341     7AYM_A   25.378             331        177       13
    1119 Query_1615341     6N3N_A   30.667             225        115        8
    1120 Query_1615341     4OH4_A   28.309             272        166        7
    1121 Query_1615341     4CQG_A   31.220             205        130        8
    1122 Query_1615341     4LG4_A   28.319             226        138        8
    1123 Query_1615341     5XD6_A   30.556             216        130        6
    1124 Query_1615341     4YUR_A   26.357             258        174        7
    1125 Query_1615341     3MDY_A   30.357             280        147       14
    1126 Query_1615341     5WNO_A   27.027             259        174        6
    1127 Query_1615341     4XUF_A   28.980             245        133        9
    1128 Query_1615341     9Y9N_A   26.357             258        174        7
    1129 Query_1615341     6AO5_A   28.054             221        145        7
    1130 Query_1615341     3COK_A   26.357             258        174        7
    1131 Query_1615341     4JXF_A   26.357             258        174        7
    1132 Query_1615341     4LGD_A   28.054             221        145        7
    1133 Query_1615341     2OO8_X   29.352             293        160       15
    1134 Query_1615341     3TL8_A   31.797             217        132        8
    1135 Query_1615341     6JQR_A   28.980             245        133        9
    1136 Query_1615341     6OYW_A   25.670             261        171       10
    1137 Query_1615341     1FVR_A   29.110             292        162       15
    1138 Query_1615341     4BIB_A   25.692             253        165       10
    1139 Query_1615341     6VRE_A   26.087             253        164       10
    1140 Query_1615341     2BMC_A   24.573             293        182       10
    1141 Query_1615341     3UIM_A   31.797             217        132        8
    1142 Query_1615341     4YMJ_A   26.351             296        173       10
    1143 Query_1615341     1U59_A   25.856             263        175        8
    1144 Query_1615341     3COH_A   24.915             293        181       10
    1145 Query_1615341     3QBN_A   25.338             296        176       11
    1146 Query_1615341     4PRJ_A   24.915             293        181       10
    1147 Query_1615341     8XB1_A   28.980             245        133        9
    1148 Query_1615341     5X02_A   28.980             245        133        9
    1149 Query_1615341     6BFN_A   32.857             210        118        9
    1150 Query_1615341     1RJB_A   28.980             245        133        9
    1151 Query_1615341     3R21_A   24.662             296        184       10
    1152 Query_1615341     4BF2_A   25.670             261        171       10
    1153 Query_1615341     1GZK_A   28.400             250        156       10
    1154 Query_1615341     7T6F_A   25.267             281        179        9
    1155 Query_1615341     3D0E_A   28.400             250        156       10
    1156 Query_1615341     6N3L_A   30.222             225        116        8
    1157 Query_1615341     2XRU_A   24.573             293        182       10
    1158 Query_1615341     1O6K_A   28.400             250        156       10
    1159 Query_1615341     1GZN_A   28.400             250        156       10
    1160 Query_1615341     1MRV_A   28.400             250        156       10
    1161 Query_1615341     6E2M_A   25.670             261        171       10
    1162 Query_1615341     2JED_A   28.155             206        136        5
    1163 Query_1615341     5VIL_A   25.670             261        171       10
    1164 Query_1615341     6IL3_A   28.980             245        133        9
    1165 Query_1615341     2J4Z_A   24.573             293        182       10
    1166 Query_1615341     2JDO_A   28.400             250        156       10
    1167 Query_1615341     5F9E_A   28.155             206        136        5
    1168 Query_1615341     7ZTL_A   24.573             293        182       10
    1169 Query_1615341     3VW6_A   26.087             253        164       10
    1170 Query_1615341     2CLQ_A   25.670             261        171       10
    1171 Query_1615341     4Q9Z_A   28.155             206        136        5
    1172 Query_1615341     6HJJ_A   25.000             296        177       11
    1173 Query_1615341     6SEQ_A   28.175             252        150        9
    1174 Query_1615341     3UNZ_A   24.573             293        182       10
    1175 Query_1615341     3NRM_A   24.662             296        184       10
    1176 Query_1615341     1O6L_A   28.400             250        156       10
    1177 Query_1615341     6XIH_A   26.087             253        164       10
    1178 Query_1615341     5UOR_A   26.087             253        164       10
    1179 Query_1615341     8C1M_A   24.585             301        188       10
    1180 Query_1615341     4PX6_A   28.832             274        161       11
    1181 Query_1615341     5OBJ_A   24.585             301        188       10
    1182 Query_1615341     8Q61_A   28.400             250        156       10
    1183 Query_1615341     8GUW_A   25.000             296        177       11
    1184 Query_1615341     3EMG_A   28.832             274        161       11
    1185 Query_1615341     5DN3_A   24.585             301        188       10
    1186 Query_1615341     4J8M_A   24.573             293        182       10
    1187 Query_1615341     1MQ4_A   24.585             301        188       10
    1188 Query_1615341     2WQB_A   29.010             293        161       15
    1189 Query_1615341     5AAD_A   25.000             296        177       11
    1190 Query_1615341     4RX7_A   28.832             274        161       11
    1191 Query_1615341     2OZO_A   25.475             263        176        8
    1192 Query_1615341     5UOX_A   26.087             253        164       10
    1193 Query_1615341     4RX9_A   28.832             274        161       11
    1194 Query_1615341     6OYT_A   26.087             253        164       10
    1195 Query_1615341     3TUB_A   28.832             274        161       11
    1196 Query_1615341     9BZG_A   25.000             296        177       11
    1197 Query_1615341     1XJD_A   28.155             206        136        5
    1198 Query_1615341     4CEG_A   25.000             296        177       11
    1199 Query_1615341     4Q5J_A   32.178             202        123        5
    1200 Query_1615341     3TUC_A   28.832             274        161       11
    1201 Query_1615341     4YJO_A   28.832             274        161       11
    1202 Query_1615341     6J5T_A   27.018             285        171       11
    1203 Query_1615341     9C1W_A   28.400             250        156       10
    1204 Query_1615341     8OF5_A   25.000             296        177       11
    1205 Query_1615341     5LXM_A   25.000             296        177       11
    1206 Query_1615341     4X3J_A   29.010             293        161       15
    1207 Query_1615341     4K2R_A   25.475             263        176        8
    1208 Query_1615341     6VPJ_A   25.091             275        169       10
    1209 Query_1615341     8X5K_A   28.832             274        161       11
    1210 Query_1615341     2J50_A   24.573             293        182       10
    1211 Query_1615341     4F4P_A   28.832             274        161       11
    1212 Query_1615341     1MUO_A   24.573             293        182       10
    1213 Query_1615341     5ORL_A   24.573             293        182       10
    1214 Query_1615341     3SRV_B   28.832             274        161       11
    1215 Query_1615341     3FDN_A   24.573             293        182       10
    1216 Query_1615341     4BTF_A   25.000             280        178       11
    1217 Query_1615341     3W16_A   24.573             293        182       10
    1218 Query_1615341     4C3P_A   24.573             293        182       10
    1219 Query_1615341     5OS2_A   25.000             296        177       11
    1220 Query_1615341     7QQ6_A   30.804             224        113        9
    1221 Query_1615341     4RSS_A   28.832             274        161       11
    1222 Query_1615341     8SSP_A   24.573             293        182       10
    1223 Query_1615341     8FO7_C   30.435             276        165       11
    1224 Query_1615341     3W10_A   24.573             293        182       10
    1225 Query_1615341     8SMC_C   30.435             276        165       11
    1226 Query_1615341     5OSD_A   24.573             293        182       10
    1227 Query_1615341     7QWK_A   30.804             224        113        9
    1228 Query_1615341     8C1K_A   25.000             284        171       10
    1229 Query_1615341     5TR6_A   28.832             274        161       11
    1230 Query_1615341     6C83_A   24.573             293        182       10
    1231 Query_1615341     1XBA_A   28.832             274        161       11
    1232 Query_1615341     5ZAN_A   24.573             293        182       10
    1233 Query_1615341     6I2U_A   24.573             293        182       10
    1234 Query_1615341     5OS5_A   24.573             293        182       10
    1235 Query_1615341     4PV0_A   28.832             274        161       11
    1236 Query_1615341     2X6D_A   24.573             293        182       10
    1237 Query_1615341     8C15_A   25.000             284        171       10
    1238 Query_1615341     6VNO_A   30.576             278        162       13
    1239 Query_1615341     5C26_A   28.832             274        161       11
    1240 Query_1615341     8TXZ_A   30.576             278        162       13
    1241 Query_1615341     9DMI_A   30.576             278        162       13
    1242 Query_1615341     1OL5_A   24.573             293        182       10
    1243 Query_1615341     2WTW_A   24.573             293        182       10
    1244 Query_1615341     6VOV_A   28.832             274        161       11
    1245 Query_1615341     3LAU_A   24.573             293        182       10
    1246 Query_1615341     2WTV_A   24.573             293        182       10
    1247 Query_1615341     4YJQ_A   28.832             274        161       11
    1248 Query_1615341     2C6E_A   24.573             293        182       10
    1249 Query_1615341     6VP8_A   30.576             278        162       13
    1250 Query_1615341     2C6D_A   24.573             293        182       10
    1251 Query_1615341     2XNG_A   24.573             293        182       10
    1252 Query_1615341     5Y5T_A   28.832             274        161       11
    1253 Query_1615341     4FL3_A   28.832             274        161       11
    1254 Query_1615341     9KDS_A   25.000             296        177       11
    1255 Query_1615341     4DFL_A   28.832             274        161       11
    1256 Query_1615341     8RRQ_A   28.832             274        161       11
    1257 Query_1615341     8PS7_A   27.018             285        146       12
    1258 Query_1615341     7O2V_A   24.573             293        182       10
    1259 Query_1615341     4FL2_A   28.832             274        161       11
    1260 Query_1615341     4BN1_A   24.573             293        182       10
    1261 Query_1615341     8HOD_A   30.396             227        139        6
    1262 Query_1615341     4UZD_A   25.000             296        177       11
    1263 Query_1615341     3E5A_A   24.573             293        182       10
    1264 Query_1615341     7LHT_A   30.576             278        162       13
    1265 Query_1615341     8TZB_A   30.576             278        162       13
    1266 Query_1615341     3HA6_A   24.573             293        182       10
    1267 Query_1615341     4JAI_A   24.252             301        189       10
    1268 Query_1615341     8OKU_A   26.263             297        177       13
    1269 Query_1615341     5DOS_A   24.555             281        176        9
    1270 Query_1615341     3EFW_A   24.573             293        182       10
    1271 Query_1615341     5EW9_A   24.632             272        174        9
    1272 Query_1615341     3H0Y_A   24.573             293        182       10
    1273 Query_1615341     5DT4_A   24.555             281        176        9
    1274 Query_1615341     9CI3_A   30.576             278        162       13
    1275 Query_1615341     5DT3_A   24.555             281        176        9
    1276 Query_1615341     9D8Z_A   30.556             216        120        9
    1277 Query_1615341     8R4O_A   25.728             206        143        5
    1278 Query_1615341     6HM6_A   28.727             275        160       12
    1279 Query_1615341     3SRV_A   28.727             275        160       12
    1280 Query_1615341     6VPG_A   24.632             272        174        9
    1281 Query_1615341     3H9R_A   30.556             216        120        9
    1282 Query_1615341     9RDA_A   30.556             216        120        9
    1283 Query_1615341     6VPH_A   24.632             272        174        9
    1284 Query_1615341     6XKA_A   24.632             272        174        9
    1285 Query_1615341     9KS6_A   24.573             293        182       10
    1286 Query_1615341     9CHO_A   30.576             278        162       13
    1287 Query_1615341     6VPL_A   24.632             272        174        9
    1288 Query_1615341     9D8F_A   30.556             216        120        9
    1289 Query_1615341     6JUX_A   30.556             216        120        9
    1290 Query_1615341     8HO6_A   30.396             227        139        6
    1291 Query_1615341     4DYM_A   30.556             216        120        9
    1292 Query_1615341     5TOS_A   25.566             309        193       10
    1293 Query_1615341   9ESA_AAA   25.912             274        177        9
    1294 Query_1615341     6VPI_A   24.727             275        170       10
    1295 Query_1615341     8JF4_A   25.091             275        169       10
    1296 Query_1615341     9P6A_A   26.996             263        163       11
    1297 Query_1615341     8TZG_A   30.435             276        165       11
    1298 Query_1615341     3W2C_A   25.091             275        169       10
    1299 Query_1615341     2W1C_A   25.091             275        169       10
    1300 Query_1615341     8HOA_A   29.956             227        140        6
    1301 Query_1615341     2H6D_A   29.224             219        142        7
    1302 Query_1615341     2W1D_A   25.091             275        169       10
    1303 Query_1615341     6UNQ_A   30.556             216        120        9
    1304 Query_1615341     6UNR_A   30.556             216        120        9
    1305 Query_1615341     6GVX_A   30.244             205        132        8
    1306 Query_1615341     9F31_A   30.244             205        132        8
    1307 Query_1615341     4BKY_A   30.244             205        132        8
    1308 Query_1615341     5TWU_A   30.244             205        132        8
    1309 Query_1615341     7CTX_A   28.431             204        131        6
    1310 Query_1615341     6XR4_A   30.576             278        162       13
    1311 Query_1615341     6EIX_A   30.093             216        121        9
    1312 Query_1615341     5TWL_A   30.244             205        132        8
    1313 Query_1615341     2WQE_A   25.091             275        169       10
    1314 Query_1615341     3MY0_A   31.019             216        119        9
    1315 Query_1615341     2YZA_A   28.440             218        145        6
    1316 Query_1615341     6VXR_A   30.244             205        132        8
    1317 Query_1615341     3DAJ_A   24.080             299        184       11
    1318 Query_1615341     6KZC_A   26.129             310        170       11
    1319 Query_1615341     4ZHX_A   29.167             216        140        7
    1320 Query_1615341     5TVT_A   30.244             205        132        8
    1321 Query_1615341     5ISO_A   29.167             216        140        7
    1322 Query_1615341     4CFE_A   29.167             216        140        7
    1323 Query_1615341     2Y7J_A   25.357             280        176       13
    1324 Query_1615341     8TZC_A   30.216             278        163       13
    1325 Query_1615341     8BIK_A   29.167             216        140        7
    1326 Query_1615341     7KX8_A   26.471             238        152        9
    1327 Query_1615341     5IH8_A   30.244             205        132        8
    1328 Query_1615341     7F3G_A   26.471             238        152        9
    1329 Query_1615341     7KXW_A   26.471             238        152        9
    1330 Query_1615341     5JZN_A   26.471             238        152        9
    1331 Query_1615341     7LI3_A   30.216             278        163       13
    1332 Query_1615341     4M69_A   27.925             265        160        9
    1333 Query_1615341     4UMQ_A   30.049             203        131        8
    1334 Query_1615341     8TYQ_A   30.216             278        163       13
    1335 Query_1615341     3MTF_A   30.233             215        120        9
    1336 Query_1615341     6VBZ_A   24.275             276        185       11
    1337 Query_1615341     9IC2_A   29.167             216        140        7
    1338 Query_1615341     4UMT_A   30.049             203        131        8
    1339 Query_1615341     7OPO_A   28.270             237        146        9
    1340 Query_1615341     5EZV_A   29.167             216        140        7
    1341 Query_1615341     5D9K_A   28.270             237        146        9
    1342 Query_1615341     3G51_A   28.270             237        146        9
    1343 Query_1615341     5JZJ_A   26.471             238        152        9
    1344 Query_1615341     4QFG_A   29.167             216        140        7
    1345 Query_1615341     4NUS_A   28.270             237        146        9
    1346 Query_1615341     3D14_A   24.080             299        184       11
    1347 Query_1615341     6NPZ_A   28.287             251        156       11
    1348 Query_1615341     8XFY_A   28.270             237        146        9
    1349 Query_1615341     3CQU_A   28.016             257        149       12
    1350 Query_1615341     9S9T_A   29.851             201        129        6
    1351 Query_1615341     6BX6_A   28.372             215        143        6
    1352 Query_1615341     8UWR_A   30.233             215        120        9
    1353 Query_1615341     6KYQ_A   26.267             217        142        7
    1354 Query_1615341     4GV1_A   28.287             251        156       11
    1355 Query_1615341     4RER_A   30.189             212        131        9
    1356 Query_1615341     4REW_A   30.189             212        131        9
    1357 Query_1615341     3OCB_A   28.016             257        149       12
    1358 Query_1615341     6BUU_A   28.016             257        149       12
    1359 Query_1615341     2I0E_A   29.851             201        129        6
    1360 Query_1615341     4EL9_A   28.270             237        146        9
    1361 Query_1615341     6KYR_A   26.267             217        142        7
    1362 Query_1615341     6CCY_A   28.287             251        156       11
    1363 Query_1615341     3UBD_A   28.270             237        146        9
    1364 Query_1615341     4D2P_A   29.557             203        132        8
    1365 Query_1615341     2Z7Q_A   27.397             219        136        7
    1366 Query_1615341     3OHT_A   28.033             239        138        9
    1367 Query_1615341     4CFH_A   29.817             218        136        9
    1368 Query_1615341     4CZU_A   25.118             211        150        6
    1369 Query_1615341     6GR8_A   25.912             274        177        9
    1370 Query_1615341     4RED_A   30.189             212        131        9
    1371 Query_1615341     4CZT_A   25.118             211        150        6
    1372 Query_1615341     7JIJ_A   30.189             212        131        9
    1373 Query_1615341     5DE2_A   26.087             207        138        7
    1374 Query_1615341     8DFP_A   33.696             184        109        5
    1376 Query_1615341     7JHG_A   30.189             212        131        9
    1377 Query_1615341     3ZDU_A   25.839             298        197       12
    1378 Query_1615341     6C9F_A   30.189             212        131        9
    1379 Query_1615341     8DFQ_A   33.696             184        109        5
    1381 Query_1615341     8DFM_A   33.696             184        109        5
    1383 Query_1615341     6NPY_B   26.087             207        138        7
    1384 Query_1615341     6C9H_A   30.189             212        131        9
    1385 Query_1615341     2DWB_A   24.232             293        183       10
    1386 Query_1615341     1ZYC_A   25.468             267        166        5
    1387 Query_1615341     8SXN_A   25.604             207        139        7
    1388 Query_1615341     1ZY4_A   25.468             267        166        5
    1389 Query_1615341     8E05_A   27.891             294        157       12
    1390 Query_1615341     2FH9_A   29.048             210        135        8
    1391 Query_1615341     3MN3_A   29.048             210        135        8
    1392 Query_1615341     3HYH_A   29.048             210        135        8
    1393 Query_1615341     3O96_A   29.435             248        157       11
    1394 Query_1615341     3MTL_A   28.972             214        140        7
    1395 Query_1615341     4EJN_A   29.435             248        157       11
    1396 Query_1615341     8E04_A   27.891             294        157       12
    1397 Query_1615341     3DAE_A   29.048             210        135        8
    1398 Query_1615341     8FAC_A   27.891             294        157       12
    1399 Query_1615341     6S73_A   26.087             207        138        7
    1400 Query_1615341     4M66_A   27.925             265        160        9
    1401 Query_1615341     9H59_C   25.604             207        139        7
    1402 Query_1615341     9NFQ_C   25.604             207        139        7
    1403 Query_1615341     8WS0_A   25.604             207        139        7
    1404 Query_1615341     4FSN_A   27.854             219        133        9
    1405 Query_1615341     4A4X_A   27.273             220        131        8
    1406 Query_1615341     1IA8_A   27.854             219        133        9
    1407 Query_1615341     4QYE_A   27.854             219        133        9
    1408 Query_1615341     5KCV_A   29.435             248        157       11
    1409 Query_1615341     7AKO_A   28.641             206        129        8
    1410 Query_1615341     2WQM_A   25.604             207        139        7
    1411 Query_1615341     7AKM_A   28.641             206        129        8
    1412 Query_1615341     1ZLT_A   27.854             219        133        9
    1413 Query_1615341     4FSM_A   27.854             219        133        9
    1414 Query_1615341     4FSY_A   27.854             219        133        9
    1415 Query_1615341     2X8E_A   27.854             219        133        9
    1416 Query_1615341     2BR1_A   27.854             219        133        9
    1417 Query_1615341     8QGY_A   25.692             253        165       10
    1418 Query_1615341     2HOG_A   27.982             218        132        9
    1419 Query_1615341     5OQ5_A   28.311             219        132        9
    1420 Query_1615341     2QHM_A   27.982             218        132        9
    1421 Query_1615341     7APJ_A   29.435             248        157       11
    1422 Query_1615341     3ORM_A   28.251             223        135        8
    1423 Query_1615341     4FSW_A   27.854             219        133        9
    1424 Query_1615341     5F4N_A   27.854             219        133        9
    1425 Query_1615341     4D28_A   26.923             208        142        7
    1426 Query_1615341     3JVR_A   27.854             219        133        9
    1427 Query_1615341     6TM5_S   26.941             219        133        7
    1428 Query_1615341     2ACX_A   27.586             232        134        7
    1429 Query_1615341     2E9V_A   27.854             219        133        9
    1430 Query_1615341     4FSZ_A   27.854             219        133        9
    1431 Query_1615341     3OT3_A   27.854             219        133        9
    1432 Query_1615341     8UW7_A   29.435             248        157       11
    1433 Query_1615341     3NYN_A   27.586             232        134        7
    1434 Query_1615341     3IW4_A   28.358             201        132        5
    1435 Query_1615341     2AYP_A   27.854             219        133        9
    1436 Query_1615341     8UVY_A   29.435             248        157       11
    1437 Query_1615341     6TM5_Q   26.941             219        133        7
    1438 Query_1615341     2E9N_A   27.854             219        133        9
    1439 Query_1615341     2YDJ_A   27.854             219        133        9
    1440 Query_1615341     4FT3_A   27.854             219        133        9
    1441 Query_1615341     6TUA_A   26.073             303        189       10
    1442 Query_1615341     8EJ4_K   25.118             211        143        7
    1443 Query_1615341     4FST_A   27.854             219        133        9
    1444 Query_1615341     1MRU_A   28.251             223        135        8
    1445 Query_1615341     2W5A_A   26.941             219        133        7
    1446 Query_1615341     4FR4_A   29.500             200        130        6
    1447 Query_1615341     2JAV_A   26.941             219        133        7
    1448 Query_1615341     8A3T_S   30.286             175        111        7
    1449 Query_1615341     5U94_A   28.251             223        135        8
    1450 Query_1615341     4AF3_A   25.589             297        186       11
    1451 Query_1615341     3F69_A   28.251             223        135        8
    1452 Query_1615341     6B2P_A   28.251             223        135        8
    1453 Query_1615341     4RA4_A   28.358             201        132        5
    1454 Query_1615341     6I2P_A   28.251             223        135        8
    1455 Query_1615341     7PUE_A   26.210             248        165        7
    1456 Query_1615341     3ORI_A   28.251             223        135        8
    1457 Query_1615341     1ZYS_A   27.854             219        133        9
    1458 Query_1615341     3F61_A   28.251             223        135        8
    1459 Query_1615341     4G31_A   29.224             219        121        7
    1460 Query_1615341     1ZXE_A   25.094             267        167        5
    1461 Query_1615341     1O6Y_A   28.251             223        135        8
    1462 Query_1615341     8WS1_A   25.121             207        140        7
    1463 Query_1615341     4FSU_A   27.854             219        133        9
    1464 Query_1615341     2R5T_A   26.210             248        165        7
    1465 Query_1615341     9IWX_A   27.547             265        161        9
    1466 Query_1615341     6G78_A   26.432             227        153        7
    1467 Query_1615341     5I3O_A   24.219             256        172        9
    1468 Query_1615341     4PNI_A   25.616             203        136        5
    1469 Query_1615341     4W9W_A   24.219             256        172        9
    1470 Query_1615341     8VSU_C   23.894             226        162        4
    1471 Query_1615341     4W9X_A   24.219             256        172        9
    1472 Query_1615341     3QC9_A   25.616             203        136        5
    1473 Query_1615341     3C4X_A   25.616             203        136        5
    1474 Query_1615341     3C4W_A   25.616             203        136        5
    1475 Query_1615341     3T8O_A   25.616             203        136        5
    1476 Query_1615341     6G76_A   26.432             227        153        7
    1477 Query_1615341     6G77_A   26.432             227        153        7
    1478 Query_1615341     7MT8_G   25.616             203        136        5
    1479 Query_1615341     4WBO_A   25.616             203        136        5
    1480 Query_1615341     6NTD_B   48.649              74         35        1
    1481 Query_1615341     4E5A_X   27.197             239        140        9
    1482 Query_1615341     3H4J_A   29.717             212        131        8
    1483 Query_1615341     4L9I_A   25.616             203        136        5
    1484 Query_1615341     2PUU_A   27.197             239        140        9
    1485 Query_1615341     3HNG_A   31.176             170        101        5
    1486 Query_1615341     2FSL_X   27.197             239        140        9
    1487 Query_1615341     3P4K_A   27.197             239        140        9
    1488 Query_1615341     8EFJ_A   27.197             239        140        9
    1489 Query_1615341     1LEW_A   27.197             239        140        9
    1490 Query_1615341     2FST_X   27.197             239        140        9
    1491 Query_1615341     4EQM_A   25.391             256        170       11
    1492 Query_1615341     4LOO_A   27.197             239        140        9
    1493 Query_1615341     3VHE_A   30.909             165        102        4
    1494 Query_1615341     3VHK_A   30.909             165        102        4
    1495 Query_1615341     2GTM_A   27.197             239        140        9
    1496 Query_1615341     3NNX_A   27.197             239        140        9
    1497 Query_1615341     1BMK_A   27.197             239        140        9
    1498 Query_1615341     5O90_A   27.197             239        140        9
    1499 Query_1615341     5MRD_A   26.087             207        144        5
    1500 Query_1615341     3TG1_A   27.197             239        140        9
    1501 Query_1615341     3VID_A   30.909             165        102        4
    1502 Query_1615341     2FSO_X   27.197             239        140        9
    1503 Query_1615341     1Y6A_A   30.909             165        102        4
    1504 Query_1615341     2BAQ_A   27.615             239        139        9
    1505 Query_1615341     4TYH_B   27.197             239        140        9
    1506 Query_1615341     8YPE_A   27.197             239        140        9
    1507 Query_1615341     5O8U_A   27.197             239        140        9
    1508 Query_1615341     5O8V_A   27.197             239        140        9
    1509 Query_1615341     3ODZ_X   27.197             239        140        9
    1510 Query_1615341   8ACM_AAA   27.197             239        140        9
    1511 Query_1615341     3PWY_A   27.230             213        144        6
    1512 Query_1615341     3HVC_A   27.197             239        140        9
    1513 Query_1615341     3S3I_A   27.197             239        140        9
    1514 Query_1615341     3OD6_X   27.197             239        140        9
    1515 Query_1615341     3NNU_A   27.197             239        140        9
    1516 Query_1615341     6SO1_A   27.197             239        140        9
    1517 Query_1615341     3ODY_X   27.197             239        140        9
    1518 Query_1615341     5NZZ_E   27.197             239        140        9
    1519 Query_1615341     3FI4_A   27.197             239        140        9
    1520 Query_1615341     3HEC_A   27.197             239        140        9
    1521 Query_1615341     6SOI_A   27.197             239        140        9
    1522 Query_1615341     1A9U_A   27.197             239        140        9
    1523 Query_1615341     9MHB_A   27.197             239        140        9
    1524 Query_1615341     1YW2_A   27.197             239        140        9
    1525 Query_1615341     2YIS_A   27.197             239        140        9
    1526 Query_1615341     1OZ1_A   27.197             239        140        9
    1527 Query_1615341     2GFS_A   27.197             239        140        9
    1528 Query_1615341     3K3I_A   27.197             239        140        9
    1529 Query_1615341     8X3M_A   27.197             239        140        9
    1530 Query_1615341     3GCU_A   27.197             239        140        9
    1531 Query_1615341     3OEF_X   27.197             239        140        9
    1532 Query_1615341     1YWR_A   27.197             239        140        9
    1533 Query_1615341     2GHL_A   27.197             239        140        9
    1534 Query_1615341     3KL8_A   27.317             205        135        7
    1535 Query_1615341     2BAJ_A   27.197             239        140        9
    1536 Query_1615341     3K3J_A   27.197             239        140        9
    1537 Query_1615341     3ZS5_A   27.197             239        140        9
    1538 Query_1615341     1DI9_A   27.197             239        140        9
    1539 Query_1615341     1ZZL_A   27.197             239        140        9
    1540 Query_1615341     2NPQ_A   27.197             239        140        9
    1541 Query_1615341     3KQ7_A   27.197             239        140        9
    1542 Query_1615341     3MPT_A   27.197             239        140        9
    1543 Query_1615341     8VXE_A   27.197             239        140        9
    1544 Query_1615341     5ETC_A   27.197             239        140        9
    1545 Query_1615341     6KA4_A   30.435             184        109        4
    1546 Query_1615341     3PY3_A   27.197             239        140        9
    1547 Query_1615341     2BAL_A   27.197             239        140        9
    1548 Query_1615341     6TCA_B   27.197             239        140        9
    1549 Query_1615341     2OZA_B   27.197             239        140        9
    1550 Query_1615341     3KK9_A   27.317             205        135        7
    1551 Query_1615341     1M7Q_A   27.197             239        140        9
    1552 Query_1615341     3E92_A   27.197             239        140        9
    1553 Query_1615341     3KK8_A   27.317             205        135        7
    1554 Query_1615341     2BDW_A   27.317             205        135        7
    1555 Query_1615341     6ZQS_A   27.197             239        140        9
    1556 Query_1615341     3D7Z_A   27.197             239        140        9
    1557 Query_1615341     9CJ1_A   27.197             239        140        9
    1558 Query_1615341     4F9W_A   27.197             239        140        9
    1559 Query_1615341     4X7H_A   28.378             222        122        7
    1560 Query_1615341     2LGC_A   27.197             239        140        9
    1561 Query_1615341     6Y7W_A   27.197             239        140        9
    1562 Query_1615341     1H1W_A   26.761             213        145        6
    1563 Query_1615341     6SOV_A   27.197             239        140        9
    1564 Query_1615341     4R3C_A   27.197             239        140        9
    1565 Query_1615341     3DT1_A   27.197             239        140        9
    1566 Query_1615341     3D83_A   27.197             239        140        9
    1567 Query_1615341     7PVU_A   27.197             239        140        9
    1568 Query_1615341   7Z6I_AAA   27.197             239        140        9
    1569 Query_1615341     9NYT_A   27.197             239        140        9
    1570 Query_1615341     4F9Y_A   27.197             239        140        9
    1571 Query_1615341     1IAN_A   27.197             239        140        9
    1572 Query_1615341     1UU9_A   26.761             213        145        6
    1573 Query_1615341     1V0O_A   26.190             210        145        6
    1574 Query_1615341     3ION_A   26.761             213        145        6
    1575 Query_1615341     1V0B_A   26.190             210        145        6
    1576 Query_1615341     2Y8O_A   27.197             239        140        9
    1577 Query_1615341     4CT1_A   26.087             207        144        5
    1578 Query_1615341     5ETF_A   26.778             239        141        9
    1579 Query_1615341     6QYX_A   27.848             237        147        8
    1580 Query_1615341     1FOT_A   27.982             218        137        7
    1581 Query_1615341     2R7B_A   26.761             213        145        6
    1582 Query_1615341     5WJJ_A   27.197             239        140        9
    1583 Query_1615341     3H9O_A   26.761             213        145        6
    1584 Query_1615341     4A07_A   26.570             207        143        5
    1585 Query_1615341     3NUS_A   26.761             213        145        6
    1586 Query_1615341     5L4Q_A   24.561             285        188       10
    1587 Query_1615341     3ORX_A   26.761             213        145        6
    1588 Query_1615341     1Z5M_A   26.761             213        145        6
    1589 Query_1615341     5MZ3_A   27.197             239        140        9
    1590 Query_1615341     9YLH_A   27.197             239        140        9
    1591 Query_1615341     5TE0_A   24.561             285        188       10
    1592 Query_1615341     3QC4_A   26.761             213        145        6
    1593 Query_1615341     3HRC_A   26.570             207        143        5
    1594 Query_1615341     5ETI_A   26.778             239        141        9
    1595 Query_1615341     2XCH_A   26.761             213        145        6
    1596 Query_1615341     1OKY_A   26.761             213        145        6
    1597 Query_1615341     1OVE_A   27.197             239        140        9
    1598 Query_1615341     7LVH_A   24.912             285        187       10
    1599 Query_1615341     4XX9_A   26.570             207        143        5
    1600 Query_1615341     4EWQ_A   27.197             239        140        9
    1601 Query_1615341     1OB3_A   26.190             210        145        6
    1602 Query_1615341     3Q4T_A   32.353             204        111        7
    1603 Query_1615341     6M9L_A   27.197             239        140        9
    1604 Query_1615341     3RWQ_A   26.761             213        145        6
    1605 Query_1615341     2BIY_A   26.761             213        145        6
    1606 Query_1615341     6M95_A   27.197             239        140        9
    1607 Query_1615341     3RWP_A   26.761             213        145        6
    1608 Query_1615341     4WSQ_A   24.561             285        188       10
    1609 Query_1615341     7XBR_F   26.733             202        134        7
    1610 Query_1615341     3MH2_A   26.778             239        141        9
    1611 Query_1615341     2WTK_C   23.451             226        163        4
    1612 Query_1615341     4IC7_A   23.881             268        178        8
    1613 Query_1615341     3HKO_A   23.810             294        167       12
    1614 Query_1615341     8A8M_A   27.197             239        140        9
    1615 Query_1615341     9QB5_A   24.561             285        188       10
    1616 Query_1615341     4KIK_A   26.697             221        132        9
    1617 Query_1615341     3NUN_A   26.761             213        145        6
    1618 Query_1615341     4B99_A   23.881             268        178        8
    1619 Query_1615341     2PHK_A   29.412             204        121       10
    1620 Query_1615341     2QLU_A   28.700             223        127        8
    1621 Query_1615341     4GEO_A   26.778             239        141        9
    1622 Query_1615341     9F58_A   28.351             194        120        4
    1623 Query_1615341     3GCP_A   26.778             239        141        9
    1624 Query_1615341     1PHK_A   29.412             204        121       10
    1625 Query_1615341     9CMZ_A   27.119             236        124       12
    1626 Query_1615341     8GMC_A   24.561             285        188       10
    1627 Query_1615341     8OMV_A   26.244             221        133        9
    1628 Query_1615341     4E3C_A   26.244             221        133        9
    1629 Query_1615341     8U2O_A   27.119             236        124       12
    1630 Query_1615341     3NAX_A   26.761             213        145        6
    1631 Query_1615341     2XCK_A   26.291             213        146        6
    1632 Query_1615341     9BF3_A   30.070             143         86        2
    1633 Query_1615341     3NAY_A   26.761             213        145        6
    1634 Query_1615341     2W96_B   26.818             220        133        9
    1635 Query_1615341     3V3V_A   25.676             296        155       12
    1636 Query_1615341     3WE4_A   25.806             248        155       10
    1637 Query_1615341     3A60_A   25.806             248        155       10
    1638 Query_1615341     3GP0_A   26.522             230        140        9
    1639 Query_1615341     3A62_A   25.806             248        155       10
    1640 Query_1615341     4KIK_B   26.244             221        133        9
    1641 Query_1615341     4L45_A   26.210             248        154       10
    1642 Query_1615341     4L46_A   26.210             248        154       10
    1643 Query_1615341     4L43_A   26.210             248        154       10
    1644 Query_1615341     7N91_A   26.210             248        154       10
    1645 Query_1615341     8YGW_A   26.522             230        140        9
    1646 Query_1615341     7XBR_A   26.238             202        135        7
    1647 Query_1615341     1QL6_A   29.412             204        121       10
    1648 Query_1615341     4ZSJ_A   24.793             242        151        7
    1649 Query_1615341     1CM8_A   27.928             222        126        9
    1650 Query_1615341     4L42_A   25.806             248        155       10
    1651 Query_1615341     5BYY_A   24.793             242        151        7
    1652 Query_1615341     6P8E_B   26.818             220        133        9
    1653 Query_1615341     2W99_B   26.818             220        133        9
    1654 Query_1615341     5BYZ_A   24.793             242        151        7
    1655 Query_1615341     6HKM_A   24.793             242        151        7
    1656 Query_1615341     4ZSG_A   24.793             242        151        7
    1657 Query_1615341     8T7T_A   34.043             141         75        4
    1659 Query_1615341     4L3J_A   27.570             214        136        9
    1660 Query_1615341     5O7I_A   24.793             242        151        7
    1661 Query_1615341     4Y83_A   25.403             248        166       10
    1662 Query_1615341     4AGU_A   25.592             211        148        5
    1663 Query_1615341     6UNA_A   27.727             220        125        9
    1664 Query_1615341     4C57_A   26.016             246        148       11
    1665 Query_1615341     5LOH_A   25.472             212        143        6
    1666 Query_1615341     9CSK_B   26.818             220        133        9
    1667 Query_1615341     7CGA_A   27.928             222        126        9
    1668 Query_1615341     6V6A_A   24.354             271        172        9
    1669 Query_1615341     4RLO_A   27.570             214        136        9
    1670 Query_1615341     3GC8_A   26.522             230        140        9
    1671 Query_1615341     2CN5_A   25.221             226        149        6
    1672 Query_1615341     9LTA_A   25.207             242        150        8
    1673 Query_1615341     2W9F_B   26.818             220        133        9
    1674 Query_1615341     2XK9_A   25.221             226        149        6
    1675 Query_1615341     6LBA_A   29.843             191        109        5
    1676 Query_1615341     5Y8U_A   31.818             154         97        4
    1677 Query_1615341     2YCF_A   25.221             226        149        6
    1678 Query_1615341     2YCR_A   25.221             226        149        6
    1679 Query_1615341     2BFX_A   22.172             221        152        6
    1680 Query_1615341     5Z1D_A   31.818             154         97        4
    1681 Query_1615341     4C2W_A   22.172             221        152        6
    1682 Query_1615341     2W0J_A   25.221             226        149        6
    1683 Query_1615341     4B8M_A   22.172             221        152        6
    1684 Query_1615341     2VRX_A   22.172             221        152        6
    1685 Query_1615341     3TXO_A   28.019             207        137        7
    1686 Query_1615341     4C2V_A   22.172             221        152        6
    1687 Query_1615341     8FP1_A   28.019             207        137        7
    1688 Query_1615341     3COI_A   25.847             236        133       10
    1689 Query_1615341     3GC9_A   26.522             230        140        9
    1690 Query_1615341     4B8L_A   22.172             221        152        6
    1691 Query_1615341     6IB0_A   31.818             154         97        4
    1692 Query_1615341     2DYL_A   31.169             154         98        4
    1693 Query_1615341     5Y90_A   31.818             154         97        4
    1694 Query_1615341     7OVJ_A   31.818             154         97        4
    1695 Query_1615341     2XRW_A   25.161             310        163       13
    1696 Query_1615341     3DLS_A   26.244             221        124        7
    1697 Query_1615341     3NIZ_A   27.273             209        142        6
    1698 Query_1615341     4ZSL_A   24.686             239        149        7
    1699 Query_1615341     6YG4_A   31.818             154         97        4
    1700 Query_1615341     6YFZ_A   31.818             154         97        4
    1701 Query_1615341     6HKN_A   24.561             228        148        6
    1702 Query_1615341     5B2L_A   31.818             154         97        4
    1703 Query_1615341     3ALN_A   26.882             279        176       11
    1704 Query_1615341     3WZU_A   31.818             154         97        4
    1705 Query_1615341     3J4Q_D   27.572             243        152       11
    1706 Query_1615341     3FHI_A   27.572             243        152       11
    1707 Query_1615341     8V5H_A   25.463             216        141        7
    1708 Query_1615341     6YG1_A   31.169             154         98        4
    1709 Query_1615341     4EYJ_A   25.847             236        133       10
    1710 Query_1615341     7DV6_A   27.803             223        126       10
    1711 Query_1615341     6ZR5_A   25.338             296        156       12
    1712 Query_1615341     3G33_A   26.457             223        133        9
    1713 Query_1615341     2QUR_A   27.572             243        152       11
    1714 Query_1615341     2QKR_A   27.273             209        142        6
    1715 Query_1615341     8X23_A   25.847             236        133       10
    1716 Query_1615341     2XS0_A   26.897             290        159       13
    1717 Query_1615341     8YGZ_A   27.803             223        126       10
    1718 Query_1615341     7SJ3_A   26.457             223        133        9
    1719 Query_1615341     5Y7Z_A   25.532             235        145       10
    1720 Query_1615341     5FWK_K   26.457             223        133        9
    1721 Query_1615341     1UKH_A   25.338             296        156       12
    1722 Query_1615341     6YG0_A   31.169             154         98        4
    1723 Query_1615341     8X5M_A   25.338             296        156       12
    1724 Query_1615341     4YR8_A   25.338             296        156       12
    1725 Query_1615341     3QAL_E   27.572             243        152       11
    1726 Query_1615341     1ATP_E   27.572             243        152       11
    1727 Query_1615341     5E8V_A   27.803             223        126       10
    1728 Query_1615341     2ERZ_E   27.572             243        152       11
    1729 Query_1615341     7E0Z_A   27.572             243        152       11
    1730 Query_1615341     1J3H_A   27.572             243        152       11
    1731 Query_1615341     3QA8_A   27.149             221        131       10
    1732 Query_1615341     4YHJ_A   27.273             165        112        4
    1733 Query_1615341     3O17_A   25.753             299        151       13
    1734 Query_1615341     2CPK_E   27.572             243        152       11
    1735 Query_1615341     3RZF_A   27.149             221        131       10
    1736 Query_1615341     4IAC_A   27.572             243        152       11
    1737 Query_1615341     8R99_A   29.557             203        115        9
    1738 Query_1615341     4DG3_E   27.572             243        152       11
    1739 Query_1615341     3X2U_A   27.572             243        152       11
    1740 Query_1615341     4MYG_A   25.847             236        133       10
    1741 Query_1615341     9GLA_A   26.364             220        151        7
    1742 Query_1615341     2G01_A   25.753             299        151       13
    1743 Query_1615341     6CCF_A   26.364             220        149        5
    1744 Query_1615341     8X5L_A   26.939             245        151       10
    1745 Query_1615341     2GNJ_A   27.049             244        152       10
    1746 Query_1615341     3I6U_A   23.667             300        201        9
    1747 Query_1615341     6CMJ_A   25.503             298        177       10
    1748 Query_1615341     1VZO_A   25.481             208        136        7
    1749 Query_1615341     3O7L_B   27.572             243        152       11
    1750 Query_1615341     2GU8_A   28.509             228        139       11
    1751 Query_1615341     1ZRZ_A   25.263             190        129        5
    1752 Query_1615341     3QD2_B   34.862             109         54        3
    1753 Query_1615341     3A8W_A   26.207             145        101        3
    1754 Query_1615341     8R3X_A   26.207             145        101        3
    1755 Query_1615341     3I6W_A   23.667             300        201        9
    1756 Query_1615341     3L9M_A   27.160             243        153       10
    1757 Query_1615341     5LI1_A   26.207             145        101        3
    1758 Query_1615341     2GNF_A   27.049             244        152       10
    1759 Query_1615341     4DC2_A   26.207             145        101        3
    1760 Query_1615341     2PK9_A   25.822             213        145        8
    1761 Query_1615341     8RU8_A   25.551             227        148        9
    1762 Query_1615341     5LI9_A   26.207             145        101        3
    1763 Query_1615341     1CMK_E   27.869             244        150       11
    1764 Query_1615341     1CTP_E   27.869             244        150       11
    1765 Query_1615341     5LIH_A   26.207             145        101        3
    1766 Query_1615341     6FRX_A   26.939             245        151       10
    1767 Query_1615341     1CDK_A   27.869             244        150       11
    1768 Query_1615341     3ZH8_A   26.207             145        101        3
    1769 Query_1615341     3RP9_A   27.459             244        124       11
    1770 Query_1615341     8R9B_A   29.557             203        115        9
    1771 Query_1615341     4NTS_A   27.500             240        150       11
    1772 Query_1615341     9EJK_B   26.207             145        101        3
    1774 Query_1615341     2V7O_A   26.087             207        132        7
    1775 Query_1615341     4DFX_E   27.500             240        150       11
    1776 Query_1615341     4WNK_A   29.560             159        104        4
    1777 Query_1615341     8R3Y_I   26.207             145        101        3
    1778 Query_1615341     4AE9_A   26.971             241        152       10
    1779 Query_1615341     1SYK_A   27.160             243        153       11
    1780 Query_1615341     10BL_A   26.339             224        147        7
    1781 Query_1615341     4AE6_A   26.971             241        152       10
    1782 Query_1615341     10SL_A   26.339             224        147        7
    1783 Query_1615341     6IB2_A   31.169             154         98        4
    1784 Query_1615341     9HIX_J   29.557             203        115        9
    1785 Query_1615341     5UV4_A   24.413             213        150        6
    1786 Query_1615341     2BFY_A   21.560             218        153        6
    1787 Query_1615341     5OO1_A   25.225             222        155        7
    1788 Query_1615341     4DFY_A   27.160             243        153       11
    1789 Query_1615341     8P7L_J   29.557             203        115        9
    1790 Query_1615341     1JBP_E   28.070             228        140       11
    1791 Query_1615341     8JFK_C   26.961             204        126       10
    1792 Query_1615341     1BKX_A   28.070             228        140       11
    1793 Query_1615341     1L3R_E   28.070             228        140       11
    1794 Query_1615341     5EYK_A   21.560             218        153        6
    1795 Query_1615341     3FI3_A   25.830             271        140       12
    1796 Query_1615341     1APM_E   28.070             228        140       11
    1797 Query_1615341     5X3F_B   28.070             228        140       11
    1798 Query_1615341     8ORM_J   29.557             203        115        9
    1799 Query_1615341     5K3Y_A   21.560             218        153        6
    1800 Query_1615341     2WEL_A   26.087             207        132        7
    1801 Query_1615341     6MM5_E   28.070             228        140       11
    1802 Query_1615341     7UJR_A   25.000             212        132        8
    1803 Query_1615341     8UKP_E   28.070             228        140       11
    1804 Query_1615341     3PVB_A   28.070             228        140       11
    1805 Query_1615341     7E11_A   28.070             228        140       11
    1806 Query_1615341     8UKN_C   28.070             228        140       11
    1807 Query_1615341     1UA2_A   29.557             203        115        9
    1808 Query_1615341     6O9L_8   29.557             203        115        9
    1809 Query_1615341     5VLO_A   25.714             210        129        8
    1810 Query_1615341     1RDQ_E   27.572             243        152       11
    1811 Query_1615341     1XH9_A   26.829             246        150       10
    1812 Query_1615341     3NX8_A   26.749             243        154       10
    1813 Query_1615341     4WB5_A   26.749             243        154       10
    1814 Query_1615341     3TTI_A   25.830             271        140       12
    1815 Query_1615341     3MVJ_A   26.749             243        154       10
    1816 Query_1615341     6C0U_A   26.749             243        154       10
    1817 Query_1615341     8P4Z_A   29.557             203        115        9
    1818 Query_1615341     3AGM_A   26.749             243        154       10
    1819 Query_1615341     3QAM_E   27.160             243        153       11
    1820 Query_1615341     4ERW_A   25.225             222        155        7
    1821 Query_1615341     1Q61_A   26.639             244        153       10
    1822 Query_1615341     4WHZ_A   25.830             271        140       12
    1823 Query_1615341     3KVX_A   25.830             271        140       12
    1824 Query_1615341     7B55_B   25.604             207        133        7
    1825 Query_1615341     4O21_A   28.070             228        140       11
    1826 Query_1615341     9BLH_A   26.087             207        132        7
    1827 Query_1615341     2VN9_A   26.087             207        132        7
    1828 Query_1615341     1SZM_A   26.639             244        153       10
    1829 Query_1615341     1H4L_A   24.091             220        158        5
    1830 Query_1615341     4AU8_A   24.186             215        154        5
    1831 Query_1615341     2GNG_A   26.639             244        153       10
    1832 Query_1615341     2R9S_A   25.830             271        140       12
    1833 Query_1615341     2IW6_A   25.446             224        152        9
    1834 Query_1615341     4Z84_A   26.639             244        153       10
    1835 Query_1615341     3VUL_A   26.667             240        130        9
    1836 Query_1615341     2VZ6_A   25.604             207        133        7
    1837 Query_1615341     7VDP_A   24.186             215        154        5
    1838 Query_1615341     8JPB_G   26.147             218        142        6
    1839 Query_1615341     7B5O_J   29.557             203        115        9
    1840 Query_1615341     3VUM_A   26.667             240        130        9
    1841 Query_1615341     4CFU_A   25.225             222        155        7
    1842 Query_1615341     3AMA_A   26.531             245        152       10
    1843 Query_1615341     3ELJ_A   25.418             299        152       13
    1844 Query_1615341     3PTG_A   25.830             271        140       12
    1845 Query_1615341     6Q4G_A   25.225             222        155        7
    1846 Query_1615341     1GZ8_A   25.225             222        155        7
    1847 Query_1615341     3FV8_A   25.830             271        140       12
    1848 Query_1615341     7ORE_A   25.830             271        140       12
    1849 Query_1615341     1PMN_A   25.830             271        140       12
    1850 Query_1615341     3PZE_A   25.418             299        152       13
    1851 Query_1615341     5N23_A   26.423             246        152       10
    1852 Query_1615341     3OXI_A   25.830             271        140       12
    1853 Query_1615341     1GII_A   25.974             231        142       10
    1854 Query_1615341     1JNK_A   25.830             271        140       12
    1855 Query_1615341     4X21_A   25.830             271        140       12
    1856 Query_1615341     4EOS_A   25.225             222        155        7
    1857 Query_1615341     4UX9_A   25.418             299        152       13
    1858 Query_1615341     1OGU_A   25.225             222        155        7
    1859 Query_1615341     3FI2_A   25.830             271        140       12
    1860 Query_1615341     4W4V_A   25.830             271        140       12
    1861 Query_1615341     2O0U_A   25.830             271        140       12
    1862 Query_1615341     2OK1_A   25.830             271        140       12
    1863 Query_1615341     4QTD_A   25.418             299        152       13
    1864 Query_1615341     9FT9_A   25.418             299        152       13
    1865 Query_1615341     4EOP_A   25.225             222        155        7
    1866 Query_1615341     4BCM_A   25.225             222        155        7
    1867 Query_1615341     4EOQ_A   25.225             222        155        7
    1868 Query_1615341     1VYW_A   25.225             222        155        7
    1869 Query_1615341     2B1P_A   25.830             271        140       12
    1870 Query_1615341     4EOM_A   25.446             224        152        9
    1871 Query_1615341     3PXF_A   25.225             222        155        7
    1872 Query_1615341     1OIT_A   25.225             222        155        7
    1873 Query_1615341     4X3F_A   25.333             225        131        8
    1874 Query_1615341     9CKO_A   28.931             159        105        4
    1875 Query_1615341     5OO0_A   25.225             222        155        7
    1876 Query_1615341     2F7E_E   26.230             244        154       10
    1877 Query_1615341     3VUD_A   26.667             240        130        9
    1878 Query_1615341     2EXC_X   25.830             271        140       12
    1879 Query_1615341     6W4O_A   25.238             210        130        8
    1880 Query_1615341     2IW8_A   25.217             230        145        9
    1881 Query_1615341     4EOO_A   25.225             222        155        7
    1882 Query_1615341     6WJF_A   25.806             248        160       10
    1883 Query_1615341     5LW1_B   25.418             299        152       13
    1884 Query_1615341     4O38_A   26.016             246        148       11
    1885 Query_1615341     4EOJ_A   25.446             224        152        9
    1886 Query_1615341     4X3F_C   25.806             217        124        8
    1887 Query_1615341     8H6P_A   25.225             222        155        7
    1888 Query_1615341     4EON_A   25.446             224        152        9
    1889 Query_1615341     9I9J_K   25.225             222        155        7
    1890 Query_1615341     7E34_A   25.225             222        155        7
    1891 Query_1615341     1H1P_A   25.225             222        155        7
    1892 Query_1615341     1W98_A   25.225             222        155        7
    1893 Query_1615341     4EOK_A   25.446             224        152        9
    1894 Query_1615341     8USO_A   25.714             210        129        8
    1895 Query_1615341     8UV0_A   25.225             222        155        7
    1896 Query_1615341     7NVQ_A   25.225             222        155        7
    1897 Query_1615341     9DC6_A   25.806             248        160       10
    1898 Query_1615341     9NFS_A   25.806             248        160       10
    1899 Query_1615341     3EZR_A   25.225             222        155        7
    1900 Query_1615341     6GUE_A   25.225             222        155        7
    1901 Query_1615341     4OW8_A   25.333             225        131        8
    1902 Query_1615341     8FEC_B   25.806             248        160       10
    1903 Query_1615341     3BHT_A   25.225             222        155        7
    1904 Query_1615341     4TNB_A   28.931             159        105        4
    1905 Query_1615341     6PJX_A   28.931             159        105        4
    1906 Query_1615341     4WB7_A   25.806             248        160       10
    1907 Query_1615341     9HIU_B   25.225             222        155        7
    1908 Query_1615341     4EOI_A   25.217             230        145        9
    1909 Query_1615341     6INL_A   25.225             222        155        7
    1910 Query_1615341     2JGZ_A   25.225             222        155        7
    1911 Query_1615341     5N3N_A   27.160             243        153       11
    1912 Query_1615341     4I3Z_A   25.225             222        155        7
    1913 Query_1615341     1GY3_A   25.225             222        155        7
    1914 Query_1615341     5UQ1_A   25.225             222        155        7
    1915 Query_1615341     4L7F_A   25.418             299        152       13
    1916 Query_1615341     8BZO_A   25.225             222        155        7
    1917 Query_1615341     1E9H_A   25.225             222        155        7
    1918 Query_1615341     9NYQ_A   25.225             222        155        7
    1919 Query_1615341     2JDT_A   26.639             244        153       10
    1920 Query_1615341     6XD3_J   29.064             203        116        9
    1921 Query_1615341     5U6Y_A   25.238             210        130        8
    1922 Query_1615341     5K4J_A   25.225             222        155        7
    1923 Query_1615341     3PJ8_A   25.225             222        155        7
    1924 Query_1615341     1YDR_E   26.230             244        154       10
    1925 Query_1615341     6F14_A   27.160             243        153       11
    1926 Query_1615341     2UVY_A   26.639             244        153       10
    1927 Query_1615341     6Q4I_A   25.225             222        155        7
    1928 Query_1615341     1AQ1_A   25.225             222        155        7
    1929 Query_1615341     1FQ1_B   25.225             222        155        7
    1930 Query_1615341     1B38_A   25.225             222        155        7
    1931 Query_1615341     6OQI_A   25.225             222        155        7
    1932 Query_1615341     1SMH_A   26.423             246        151       10
    1933 Query_1615341     4RT7_A   29.697             165         96        6
    1935 Query_1615341     6Y05_A   27.160             243        153       11
    1936 Query_1615341     9OB2_A   25.225             222        155        7
    1937 Query_1615341     9PE7_A   24.000             225        152        5
    1938 Query_1615341     8PYR_A   29.064             203        116        9
    1939 Query_1615341     6EM7_A   27.160             243        153       11
    1940 Query_1615341     1SVH_A   26.230             244        154       10
    1941 Query_1615341     5VI9_A   26.230             244        154       10
    1942 Query_1615341     3QHR_A   25.225             222        155        7
    1943 Query_1615341     3G2F_A   27.523             218        117        9
    1944 Query_1615341     4X3F_B   31.250             112         73        2
    1945 Query_1615341     3DND_A   26.230             244        154       10
    1946 Query_1615341     5UUU_A   25.688             218        143        6
    1947 Query_1615341     1Q24_A   26.230             244        154       10
    1948 Query_1615341     1Q8W_A   26.230             244        154       10
    1949 Query_1615341     3AGL_A   26.337             243        155       10
    1950 Query_1615341     2C1A_A   26.230             244        154       10
    1951 Query_1615341     6VZK_A   24.528             212        133        8
    1952 Query_1615341     3SV0_A   26.389             216        140        6
    1953 Query_1615341     3VUG_A   26.667             240        130        9
    1954 Query_1615341     2JDS_A   26.230             244        154       10
    1955 Query_1615341     3VUK_A   26.667             240        130        9
    1956 Query_1615341     6UNP_A   27.523             218        117        9
    1957 Query_1615341     8SF8_A   26.230             244        154       10
    1958 Query_1615341    8GXQ_HI   29.064             203        116        9
    1959 Query_1615341     8YNG_A   23.982             221        159        6
    1960 Query_1615341     6XBZ_J   29.064             203        116        9
    1961 Query_1615341     1BI7_A   23.556             225        153        5
    1962 Query_1615341     3VUI_A   26.667             240        130        9
    1963 Query_1615341     3VUH_A   26.667             240        130        9
    1964 Query_1615341     4C33_A   26.230             244        154       10
    1965 Query_1615341     3KRW_A   25.688             218        143        6
    1966 Query_1615341     4WB8_A   27.193             228        142       10
    1967 Query_1615341     1JOW_B   23.556             225        153        5
    1968 Query_1615341     7PWD_A   25.688             218        143        6
    1969 Query_1615341     6TD3_B   25.991             227        149        8
    1970 Query_1615341     3CIK_A   25.688             218        143        6
    1971 Query_1615341     5NW8_A   28.636             220        133       11
    1972 Query_1615341     6C2Y_A   25.688             218        143        6
    1973 Query_1615341     8BYA_A   25.225             222        155        7
    1974 Query_1615341     5DYK_A   28.571             273        154       14
    1975 Query_1615341     5N1F_A   28.636             220        133       11
    1976 Query_1615341     1OMW_A   25.688             218        143        6
    1977 Query_1615341     3PSC_A   25.688             218        143        6
    1978 Query_1615341     4C0T_A   24.597             248        137       11
    1979 Query_1615341     4WIH_A   28.636             220        133       11
    1980 Query_1615341     6B2Q_A   25.333             225        131        8
    1981 Query_1615341     1UNG_A   23.636             220        159        5
    1982 Query_1615341     8I0M_A   23.556             225        153        5
    1983 Query_1615341     5UKK_A   25.688             218        143        6
    1984 Query_1615341     5HE1_A   25.688             218        143        6
    1985 Query_1615341     3NUP_A   23.556             225        153        5
    1986 Query_1615341     4MK0_A   25.688             218        143        6
    1987 Query_1615341     9I9I_K   26.540             211        142        6
    1988 Query_1615341     5MHQ_A   24.775             222        156        7
    1989 Query_1615341     2UZT_A   26.638             229        142       10
    1990 Query_1615341     5HE3_A   25.688             218        143        6
    1991 Query_1615341     5HE0_A   25.688             218        143        6
    1992 Query_1615341     4IZ5_A   26.471             238        150        9
    1993 Query_1615341     2VO0_A   27.074             229        141       10
    1994 Query_1615341     4WB6_B   26.337             243        155       10
    1995 Query_1615341     4WB6_A   26.337             243        155       10
    1996 Query_1615341     5EFQ_A   26.471             204        134        6
    1997 Query_1615341     2CJM_A   25.225             222        155        7
    1998 Query_1615341     5YV8_A   31.126             151         94        4
    1999 Query_1615341     2ZV2_A   31.126             151         94        4
    2000 Query_1615341     1Q8T_A   25.820             244        155       10
    2001 Query_1615341     6RFP_A   26.471             238        150        9
    2002 Query_1615341     7VDU_A   24.775             222        156        7
    2003 Query_1615341     1XH7_A   25.820             244        155       10
    2004 Query_1615341     3ZO2_A   25.820             244        155       10
    2005 Query_1615341     1H01_A   24.775             222        156        7
    2006 Query_1615341     1OIR_A   24.775             222        156        7
    2007 Query_1615341     3VN9_A   29.747             158        103        5
    2008 Query_1615341     8S79_A   24.380             242        130        7
    2009 Query_1615341     8A8M_B   23.839             323        219       10
    2010 Query_1615341     8XEY_A   24.576             236        139       11
    2011 Query_1615341     5UY6_A   31.126             151         94        4
    2012 Query_1615341     4C34_A   26.638             229        142       10
    2013 Query_1615341     2AC5_A   30.769             143         88        5
    2014 Query_1615341     4XRL_A   26.471             238        150        9
    2015 Query_1615341     2QR8_A   24.576             236        139       11
    2016 Query_1615341     6OQL_A   23.611             216        146        5
    2017 Query_1615341     4JG6_A   24.576             236        139       11
    2018 Query_1615341     3RNY_A   26.033             242        151       12
    2019 Query_1615341     4BCF_A   23.529             204        139        6
    2020 Query_1615341     8WF4_A   26.293             232        143       12
    2021 Query_1615341     4NIF_A   25.941             239        135       13
    2022 Query_1615341     9HW6_B   26.180             233        134        9
    2023 Query_1615341    9IJJ_4Z   26.147             218        142        6
    2024 Query_1615341     8I0L_A   23.529             204        139        6
    2025 Query_1615341     4IZA_A   26.471             238        150        9
    2026 Query_1615341     4S2Z_A   26.471             238        150        9
    2027 Query_1615341     5O1S_A   24.576             236        139       11
    2028 Query_1615341     3KN5_A   25.604             207        125        9
    2029 Query_1615341     6OPG_A   26.471             238        150        9
    2030 Query_1615341     5WP1_A   26.471             238        150        9
    2031 Query_1615341     7UKZ_A   26.415             212        143        6
    2032 Query_1615341     2BHH_A   24.775             222        156        7
    2033 Query_1615341     3NPC_A   27.468             233        137        9
    2034 Query_1615341     4EC8_A   23.529             204        139        6
    2035 Query_1615341     7N8T_A   27.468             233        137        9
    2036 Query_1615341     9HVX_A   26.180             233        134        9
    2037 Query_1615341     2WNT_A   25.941             239        135       13
    2038 Query_1615341     5V62_A   26.471             238        150        9
    2039 Query_1615341     3E7O_A   27.468             233        137        9
    2040 Query_1615341     7OPM_A   26.471             238        150        9
    2041 Query_1615341     3BLH_A   23.529             204        139        6
    2042 Query_1615341     6U2G_A   26.394             269        172       11
    2043 Query_1615341     6OPK_A   26.471             238        150        9
    2044 Query_1615341     7MFD_B   26.394             269        172       11
    2045 Query_1615341     4RZ7_A   37.500             104         54        4
    2046 Query_1615341     2QR7_A   24.891             229        147       10
    2047 Query_1615341     3ZUV_A   26.471             238        150        9
    2048 Query_1615341     5KKR_C   26.394             269        172       11
    2049 Query_1615341     2ERK_A   26.471             238        150        9
    2050 Query_1615341     4QP1_A   26.471             238        150        9
    2051 Query_1615341     5EZR_A   37.500             104         54        4
    2052 Query_1615341     4IZ7_A   26.471             238        150        9
    2053 Query_1615341     7CML_A   27.468             233        137        9
    2054 Query_1615341     7W5O_A   26.471             238        150        9
    2055 Query_1615341     4XJ0_A   26.471             238        150        9
    2056 Query_1615341     9HW6_A   26.180             233        134        9
    2057 Query_1615341     6G9J_A   26.471             238        150        9
    2058 Query_1615341     1PME_A   27.197             239        147       10
    2059 Query_1615341     6GZH_A   23.529             204        139        6
    2060 Query_1615341     6G9K_A   26.471             238        150        9
    2061 Query_1615341     2GPH_A   26.471             238        150        9
    2062 Query_1615341     6W4O_O   24.286             210        132        8
    2063 Query_1615341     5LCK_A   26.471             238        150        9
    2064 Query_1615341     9Z8K_B   25.862             232        134        8
    2065 Query_1615341     5K4I_A   26.471             238        150        9
    2066 Query_1615341     8ZJV_A   26.471             238        150        9
    2067 Query_1615341     4QTA_A   26.471             238        150        9
    2068 Query_1615341     4XOY_A   26.471             238        150        9
    2069 Query_1615341     3FME_A   29.114             158        104        5
    2070 Query_1615341     8RMB_A   26.471             238        150        9
    2071 Query_1615341     8PSR_A   26.471             238        150        9
    2072 Query_1615341     2Y9Q_A   26.471             238        150        9
    2073 Query_1615341     4QP4_A   26.471             238        150        9
    2074 Query_1615341     3O71_A   26.471             238        150        9
    2075 Query_1615341     4FV6_A   26.471             238        150        9
    2076 Query_1615341     5AWM_A   24.370             238        138        7
    2077 Query_1615341     8CHF_E   26.022             269        173       10
    2078 Query_1615341     6W9E_A   23.529             204        139        6
    2079 Query_1615341     7XQK_A   25.328             229        136        9
    2080 Query_1615341     3MI9_A   23.529             204        139        6
    2081 Query_1615341     2OJG_A   26.471             238        150        9
    2082 Query_1615341     3C9W_A   26.471             238        150        9
    2083 Query_1615341     5L1Z_A   23.529             204        139        6
    2084 Query_1615341     4OR5_A   23.529             204        139        6
    2085 Query_1615341     6OTS_A   26.471             238        150        9
    2086 Query_1615341     9AXA_B   26.022             269        173       10
    2087 Query_1615341     1WZY_A   26.471             238        150        9
    2088 Query_1615341     1TVO_A   26.471             238        150        9
    2089 Query_1615341     6RFO_A   25.833             240        143       10
    2090 Query_1615341     10JU_B   25.862             232        134        8
    2091 Query_1615341     2FYS_A   26.471             238        150        9
    2092 Query_1615341     6OT6_A   26.471             238        150        9
    2093 Query_1615341     3ZU7_A   26.471             238        150        9
    2094 Query_1615341     4QYY_A   26.471             238        150        9
    2095 Query_1615341     2Z7L_A   26.471             238        150        9
    2096 Query_1615341     6NBS_A   26.471             238        150        9
    2097 Query_1615341     4IMY_A   23.529             204        139        6
    2098 Query_1615341     4ZZM_A   26.471             238        150        9
    2099 Query_1615341     3R63_A   26.471             238        150        9
    2100 Query_1615341     9SKQ_A   24.268             239        137        8
    2101 Query_1615341     4XOZ_A   26.471             238        150        9
    2102 Query_1615341     9O0U_B   26.022             269        173       10
    2103 Query_1615341     8BW9_C   29.048             210        135        9
    2104 Query_1615341     10JU_A   25.862             232        134        8
    2105 Query_1615341     6PP9_B   26.022             269        173       10
    2106 Query_1615341     4FUX_A   26.471             238        150        9
    2107 Query_1615341     7UGB_A   26.471             238        150        9
    2108 Query_1615341     9TYG_A   26.022             269        173       10
    2109 Query_1615341     8ELC_A   27.468             233        137        9
    2110 Query_1615341     5BUE_A   26.471             238        150        9
    2111 Query_1615341     9QQJ_A   26.471             238        150        9
    2112 Query_1615341     3SOA_A   24.286             210        132        8
    2113 Query_1615341     3QYW_A   26.471             238        150        9
    2114 Query_1615341     4Y72_A   24.268             239        137        8
    2115 Query_1615341     6DMG_A   26.471             238        150        9
    2116 Query_1615341     9Z8K_A   25.862             232        134        8
    2117 Query_1615341     7NJ0_B   24.268             239        137        8
    2118 Query_1615341     6NYB_B   26.022             269        173       10
    2119 Query_1615341     5LCJ_A   26.471             238        150        9
    2120 Query_1615341     6Q0T_C   26.022             269        173       10
    2121 Query_1615341     8RLX_A   26.471             238        150        9
    2122 Query_1615341     9TYG_B   26.471             238        150        9
    2123 Query_1615341     4AN2_A   24.922             321        201       12
    2124 Query_1615341     4YC6_A   24.268             239        137        8
    2125 Query_1615341     6FI6_A   26.471             238        150        9
    2126 Query_1615341     4XNE_A   26.471             238        150        9
    2127 Query_1615341     4XP2_A   26.471             238        150        9
    2128 Query_1615341     8AOC_A   26.471             238        150        9
    2129 Query_1615341     6Z45_A   23.039             204        140        6
    2130 Query_1615341     4NST_A   25.490             204        136        6
    2131 Query_1615341     4O6E_A   26.027             219        144        7
    2132 Query_1615341     4FV7_A   26.471             238        150        9
    2133 Query_1615341     4GSB_A   26.471             238        150        9
    2134 Query_1615341     6FXV_A   26.471             238        150        9
    2135 Query_1615341     3SA0_A   26.471             238        150        9
    2136 Query_1615341     2ZOQ_A   23.853             218        154        5
    2137 Query_1615341     5NHH_A   26.471             238        150        9
    2138 Query_1615341     1GOL_A   26.050             238        151        9
    2139 Query_1615341     6GDM_A   26.471             238        150        9
    2140 Query_1615341     2B9H_A   24.199             281        182       11
    2141 Query_1615341     8RM2_A   26.471             238        150        9
    2142 Query_1615341     5NGU_A   26.471             238        150        9
    2143 Query_1615341     8K5R_A   23.039             204        140        6
    2144 Query_1615341     4N0S_A   26.471             238        150        9
    2145 Query_1615341     4JG8_A   24.153             236        140       11
    2146 Query_1615341     4S30_A   26.471             238        150        9
    2147 Query_1615341     3TEI_A   26.471             238        150        9
    2148 Query_1615341     9TU0_A   23.853             218        154        5
    2149 Query_1615341     4QTB_A   23.853             218        154        5
    2150 Query_1615341     4I5H_A   27.197             239        147       10
    2151 Query_1615341     9JK1_A   25.490             204        136        6
    2152 Query_1615341     2Y4I_C   28.230             209        138        8
    2153 Query_1615341     4CXA_A   25.490             204        136        6
    2154 Query_1615341     4UN0_C   25.490             204        136        6
    2155 Query_1615341     2AC3_A   30.070             143         89        5
    2156 Query_1615341     8XFM_A   30.070             143         89        5
    2157 Query_1615341     6B3E_A   25.490             204        136        6
    2158 Query_1615341     7F2X_A   24.062             320        209       10
    2159 Query_1615341     2B9F_A   24.199             281        182       11
    2160 Query_1615341     7E73_A   26.471             238        150        9
    2161 Query_1615341     6XI8_A   25.676             222        140       10
    2162 Query_1615341     2F9G_A   24.199             281        182       11
    2163 Query_1615341     4H3Q_A   26.471             238        150        9
    2164 Query_1615341     7KUE_A   25.676             222        140       10
    2165 Query_1615341     4CRS_A   27.136             199        130        6
    2166 Query_1615341     3ORN_A   28.230             209        138        8
    2167 Query_1615341     4XHL_A   24.528             265        177        8
    2168 Query_1615341     3N9X_A   26.695             236        125       11
    2169 Query_1615341     7E75_A   26.050             238        151        9
    2170 Query_1615341     2P55_A   28.230             209        138        8
    2171 Query_1615341     5WVD_A   30.128             156         97        7
    2172 Query_1615341     7W5C_A   24.638             207        126        7
    2173 Query_1615341     4U7Z_A   28.230             209        138        8
    2174 Query_1615341     1S9J_A   28.230             209        138        8
    2175 Query_1615341     5EYM_A   28.230             209        138        8
    2176 Query_1615341     2HW6_A   30.128             156         97        7
    2177 Query_1615341     3EQC_A   28.230             209        138        8
    2178 Query_1615341     5YXI_A   40.845              71         39        1
    2179 Query_1615341     7JUQ_C   28.230             209        138        8
    2180 Query_1615341     8YP4_A   28.230             209        138        8
    2181 Query_1615341     3ZLY_A   28.230             209        138        8
    2182 Query_1615341     3ZLS_A   28.230             209        138        8
    2183 Query_1615341     4AW2_A   23.016             252        175       10
    2184 Query_1615341     3DV3_A   28.230             209        138        8
    2185 Query_1615341     5HZE_A   28.230             209        138        8
    2186 Query_1615341     3MBL_A   28.230             209        138        8
    2187 Query_1615341     3W8Q_A   28.230             209        138        8
    2188 Query_1615341     7PQV_A   28.230             209        138        8
    2189 Query_1615341     6PXN_A   26.852             216        133        7
    2190 Query_1615341     5YT3_A   27.751             209        139        7
    2191 Query_1615341     6Z3U_B   25.743             202        133        9
    2192 Query_1615341     6BDL_A   27.184             206        127        8
    2193 Query_1615341     1S9I_A   28.230             209        138        8
    2194 Query_1615341     5V5Y_A   25.701             214        122        9
    2195 Query_1615341     7N3U_A   25.701             214        122        9
    2196 Query_1615341     8WDK_W   25.701             214        122        9
    2197 Query_1615341     2IN6_A   25.701             214        122        9
    2198 Query_1615341     7LV3_A   27.184             206        127        8
    2199 Query_1615341     1X8B_A   25.701             214        122        9
    2200 Query_1615341     3BI6_A   25.701             214        122        9
    2201 Query_1615341     7NAA_A   23.394             218        146        8
    2202 Query_1615341     2Z2W_A   25.701             214        122        9
    2203 Query_1615341     4BGQ_A   23.973             292        199       11
    2204 Query_1615341     8H59_A   28.736             174        103        9
    2205 Query_1615341     7T4T_A   27.184             206        127        8
    2206 Query_1615341     6DTL_A   24.510             204        124        6
    2207 Query_1615341     3ENM_A   28.662             157        104        5
    2208 Query_1615341     5Z33_A   28.736             174        103        9
    2209 Query_1615341     4OTD_A   25.180             278        180       11
    2210 Query_1615341     9AXH_A   28.230             209        138        8
    2211 Query_1615341     5CZO_A   23.755             261        178        7
    2212 Query_1615341     4UAK_A   24.506             253        170       10
    2213 Query_1615341     8ZTC_A   28.324             173        105        8
    2214 Query_1615341     3QFV_A   24.506             253        170       10
    2215 Query_1615341     7JV7_A   25.714             210        143        6
    2216 Query_1615341     5CI6_A   24.510             204        124        6
    2217 Query_1615341     3SLS_A   28.230             209        138        8
    2218 Query_1615341     6P5M_A   23.415             205        146        6
    2219 Query_1615341     3TKU_A   24.506             253        170       10
    2220 Query_1615341     7JNT_A   23.415             205        146        6
    2221 Query_1615341     7B3M_A   28.230             209        138        8
    2222 Query_1615341     5U7Q_A   23.415             205        146        6
    2223 Query_1615341     4WOT_A   23.415             205        146        6
    2224 Query_1615341     4AAA_A   22.936             218        160        4
    2225 Query_1615341     5U7R_A   23.415             205        146        6
    2226 Query_1615341     3NIE_A   26.800             250        130       12
    2227 Query_1615341     8X8X_A   23.415             205        146        6
    2228 Query_1615341     9AXM_A   27.751             209        139        7
    2229 Query_1615341     6ED6_A   23.415             205        146        6
    2230 Query_1615341     3ZLW_A   27.751             209        139        8
    2231 Query_1615341     4L6Q_A   23.415             205        146        6
    2232 Query_1615341     6PXP_A   27.315             216        132        7
    2233 Query_1615341     2F2U_A   22.927             205        147        6
    2234 Query_1615341     4KB8_A   27.074             229        136        8
    2235 Query_1615341     4XH0_A   23.596             267        177        8
    2236 Query_1615341     6RUU_A   25.271             277        180       11
    2237 Query_1615341     2H34_A   29.487             156         97        6
    2238 Query_1615341     6ZIW_I   27.586             203        131        8
    2239 Query_1615341     6E9W_A   22.927             205        147        6
    2240 Query_1615341     9JCU_A   22.927             205        147        6
    2241 Query_1615341     3TV7_A   23.448             145        107        3
    2242 Query_1615341     4W7P_A   23.448             145        107        3
    2243 Query_1615341     4QNY_A   29.139             151         89        8
    2244 Query_1615341     7JOU_A   23.448             145        107        3
    2245 Query_1615341     2ESM_A   23.448             145        107        3
    2246 Query_1615341     8ZH5_A   23.448             145        107        3
    2247 Query_1615341     2V55_A   23.448             145        107        3
    2248 Query_1615341     7S26_A   23.448             145        107        3
    2249 Query_1615341     9B3S_A   26.852             216        133        7
    2250 Query_1615341     7S25_A   23.448             145        107        3
    2251 Query_1615341     4TN6_A   26.852             216        133        7
    2252 Query_1615341     8VXF_A   26.852             216        133        7
    2253 Query_1615341     8D7M_A   26.852             216        133        7
    2254 Query_1615341     5MQV_A   26.852             216        133        7
    2255 Query_1615341     3UYS_A   26.852             216        133        7
    2256 Query_1615341     4JJR_A   26.852             216        133        7
    2257 Query_1615341     8VXD_A   26.852             216        133        7
    2258 Query_1615341     5OKT_A   26.852             216        133        7
    2259 Query_1615341     5IH4_A   26.852             216        133        7
    2260 Query_1615341     5X17_A   26.852             216        133        7
    2261 Query_1615341     8IZC_A   26.852             216        133        7
    2262 Query_1615341     7P7F_A   26.852             216        133        7
    2263 Query_1615341     6RCG_A   26.852             216        133        7
    2264 Query_1615341     1CKI_A   26.852             216        133        7
    2265 Query_1615341     4TW9_A   26.852             216        133        7
    2266 Query_1615341     4HNI_A   26.291             213        138        7
    2267 Query_1615341     5CYZ_A   23.372             261        179        7
    1658 Query_1615341     8T7T_A   28.906             128         79        3
    2268 Query_1615341    9FQR_Xr   31.915              94         55        3
    2269 Query_1615341     7UP4_A   25.121             207        106        9
    2270 Query_1615341     2ELI_A   37.037              54         30        1
    2271 Query_1615341     5L2Q_A   25.478             157        110        4
    2272 Query_1615341     2VD5_A   26.180             233        145       10
    2273 Query_1615341     3OZ6_A   23.239             284        172       11
    2274 Query_1615341     7WTT_a   29.936             157         93        6
    2275 Query_1615341     6GZD_A   29.936             157         93        6
    2276 Query_1615341     4FI1_A   24.832             149         95        6
    2277 Query_1615341     8TQ2_A   24.537             216        133        9
    2278 Query_1615341     4O2Z_A   27.523             109         71        3
    2279 Query_1615341     9H8C_A   24.537             216        133        9
    2280 Query_1615341     5X18_A   22.400             250        173        8
    2281 Query_1615341     4CRL_A   24.537             216        133        9
    2282 Query_1615341     3RGF_A   24.537             216        133        9
    2283 Query_1615341     6T41_A   24.537             216        133        9
    2284 Query_1615341     6TPA_A   24.537             216        133        9
    2285 Query_1615341     6QTG_A   24.537             216        133        9
    2286 Query_1615341     9R59_A   23.721             215        148        8
    2287 Query_1615341     8XU4_A   23.721             215        148        8
    2288 Query_1615341     6BXI_A   24.569             232        133        9
    2289 Query_1615341     5IDN_A   24.537             216        133        9
    2290 Query_1615341     5HNB_A   24.537             216        133        9
    2291 Query_1615341     5M4U_A   23.973             146         94        6
    2292 Query_1615341     5XQX_A   24.537             216        133        9
    2293 Query_1615341     5FQD_C   30.189             159         90        6
    2294 Query_1615341     9IHG_A   23.973             146         94        6
    2295 Query_1615341     3OFM_A   23.973             146         94        6
    2296 Query_1615341     5OOI_A   23.973             146         94        6
    2297 Query_1615341     5FGK_A   24.537             216        133        9
    2298 Query_1615341     7KPV_A   32.000             100         57        3
    2299 Query_1615341     6HMD_A   23.973             146         94        6
    2300 Query_1615341     9HK5_A   23.973             146         94        6
    2301 Query_1615341     3KA0_A   24.074             216        146        9
    2302 Query_1615341     6L20_A   23.973             146         94        6
    2303 Query_1615341     6QY8_A   23.973             146         94        6
    2304 Query_1615341     9OTY_C   30.189             159         90        6
    2305 Query_1615341     3E3B_X   23.973             146         94        6
    2306 Query_1615341     6T8X_A   22.120             217        152        8
    2307 Query_1615341     3UFF_A   33.333              51         30        1
    2308 Query_1615341     7PSU_A   24.242             132         86        5
    2309 Query_1615341     3UGD_A   33.333              51         30        1
    2310 Query_1615341     2JBO_A   24.074             216        146        9
    2311 Query_1615341     3GOK_A   24.074             216        146        9
    2312 Query_1615341     6TLL_A   24.242             132         86        5
    2313 Query_1615341     2PZY_A   24.074             216        146        9
    2314 Query_1615341     4MD8_E   24.242             132         86        5
    2315 Query_1615341     6TCA_A   24.074             216        146        9
    2316 Query_1615341     3R2B_A   24.074             216        146        9
    2317 Query_1615341     4MD7_E   24.242             132         86        5
    2318 Query_1615341     5OMY_A   24.242             132         86        5
    2319 Query_1615341     4TYH_A   24.074             216        146        9
    2320 Query_1615341     2OZA_A   24.074             216        146        9
    2321 Query_1615341     3FPM_A   24.074             216        146        9
    2322 Query_1615341     4FKD_A   35.294              51         29        1
    2323 Query_1615341     3R2Y_A   24.074             216        146        9
    2324 Query_1615341     2P3G_X   24.074             216        146        9
    2325 Query_1615341     1KWP_A   24.074             216        146        9
    2326 Query_1615341     3UEJ_A   33.333              51         30        1
    2327 Query_1615341     2ONL_C   24.074             216        146        9
    2328 Query_1615341     1PTQ_A   34.000              50         29        1
    2329 Query_1615341     7KND_A   34.000              50         29        1
    2330 Query_1615341     6O6Q_A   21.888             233        150        9
    2331 Query_1615341     2ENZ_A   36.000              50         28        1
    2332 Query_1615341     3UEY_A   33.333              51         30        1
    2333 Query_1615341     6L22_A   23.288             146         95        6
    2334 Query_1615341     3U87_A   24.242             132         86        5
    2335 Query_1615341     6UWA_A   36.538              52         29        1
    2336 Query_1615341     1NA7_A   22.477             218        140       10
    2337 Query_1615341     6L24_A   23.288             146         95        6
    2338 Query_1615341     5XVU_A   25.455             110         75        4
    2339 Query_1615341     2CSN_A   24.299             214        143        5
    2340 Query_1615341     1CSN_A   24.299             214        143        5
    2341 Query_1615341     2ZJW_A   24.242             132         86        5
    2342 Query_1615341     3JUH_A   24.242             132         86        5
    2343 Query_1615341     1NXK_A   26.027             146         95        6
    2344 Query_1615341     6L23_A   24.242             132         86        5
    2345 Query_1615341     5CSP_A   22.477             218        140       10
    2346 Query_1615341     5OSL_A   22.477             218        140       10
    2347 Query_1615341     6QY7_A   24.242             132         86        5
    2348 Query_1615341     5MOV_A   22.477             218        140       10
    2349 Query_1615341     9HXU_A   24.242             132         86        5
    2350 Query_1615341     4MD9_E   23.358             137         91        5
    2351 Query_1615341     1JWH_A   24.242             132         86        5
    2352 Query_1615341     3Q9W_A   24.242             132         86        5
    2353 Query_1615341     6HME_A   24.242             132         86        5
    2354 Query_1615341     5N1V_A   24.242             132         86        5
    2355 Query_1615341     7A4Q_A   22.477             218        140       10
    2356 Query_1615341     1XA6_A   36.364              66         34        2
    2357 Query_1615341     3W8L_A   24.242             132         86        5
    2358 Query_1615341     5KU8_A   24.242             132         86        5
    2359 Query_1615341     1PJK_A   24.242             132         86        5
    2360 Query_1615341     2R7I_A   24.242             132         86        5
    2361 Query_1615341     3BQC_A   24.242             132         86        5
    2362 Query_1615341   6Z83_AAA   24.242             132         86        5
    2363 Query_1615341     3NSZ_A   24.242             132         86        5
    2364 Query_1615341     3H30_A   24.242             132         86        5
    2365 Query_1615341     3NGA_A   24.242             132         86        5
    2366 Query_1615341     6Z19_B   24.242             132         86        5
    2367 Query_1615341     3MB6_A   24.242             132         86        5
    2368 Query_1615341     5ZN5_A   24.242             132         86        5
    2369 Query_1615341     3Q04_A   24.242             132         86        5
    2370 Query_1615341     6Q38_A   24.242             132         86        5
    2371 Query_1615341     6YPH_A   24.242             132         86        5
    2372 Query_1615341     1TBN_A   35.185              54         31        1
    2373 Query_1615341     5CVG_A   24.242             132         86        5
    2374 Query_1615341     5CLP_A   24.242             132         86        5
    2375 Query_1615341     7QUX_A   24.242             132         86        5
    2376 Query_1615341     6L21_A   24.242             132         86        5
    2377 Query_1615341     6YZH_A   24.242             132         86        5
    2378 Query_1615341     5ZN0_A   23.358             137         91        5
    2379 Query_1615341     3CXL_A   26.000             150         89        5
    2380 Query_1615341     9EDY_A   21.834             229        130        9
    2381 Query_1615341     6U69_A   21.586             227        133        9
    2382 Query_1615341     4DGL_C   24.545             110         76        4
    2383 Query_1615341     6JKK_A   19.910             221        152        5
    2384 Query_1615341     1KBE_A   34.091              44         28        1
    1375 Query_1615341     8DFP_A   29.545             132         79        6
    1380 Query_1615341     8DFQ_A   29.545             132         79        6
    1382 Query_1615341     8DFM_A   29.545             132         79        6
    2385 Query_1615341     5M07_A   22.749             211        128        7
    2386 Query_1615341     5M06_A   22.749             211        128        7
    2387 Query_1615341     2PVH_A   22.603             146        103        5
    2388 Query_1615341     4JRN_A   22.086             163        110        3
    2389 Query_1615341     5XKA_A   22.596             208        132        6
    2390 Query_1615341     4ANM_A   22.603             146        103        5
    2391 Query_1615341     1DS5_A   23.636             110         77        4
    2392 Query_1615341     3PVG_A   22.603             146        103        5
    2393 Query_1615341     1M2P_A   23.636             110         77        4
    2394 Query_1615341     1DAW_A   23.636             110         77        4
    2395 Query_1615341     2QC6_A   22.603             146        103        5
    2396 Query_1615341     4DGN_A   22.603             146        103        5
    2397 Query_1615341     3KXG_A   22.603             146        103        5
    2398 Query_1615341     5TS8_A   22.603             146        103        5
    2399 Query_1615341     4DGM_A   22.603             146        103        5
    781  Query_1615341     3PFQ_A   31.034              58         36        1
    2400 Query_1615341     1Y8F_A   42.857              42         20        1
    771  Query_1615341     8SE1_A   31.034              58         36        1
    773  Query_1615341     8SE2_A   31.034              58         36        1
    2401 Query_1615341     6RA0_A   39.474              38         19        1
    2402 Query_1615341     3KGA_A   21.860             215        131        8
    1934 Query_1615341     4RT7_A   27.523             109         68        4
    2403 Query_1615341     7P5Z_1   26.115             157         88        6
    2404 Query_1615341     2YUU_A   22.388              67         43        2
    2405 Query_1615341     7T7C_A   37.736              53         29        1
    2406 Query_1615341     5UE8_A   37.736              53         29        1
    2407 Query_1615341     6VP8_B   28.358              67         43        2
    2408 Query_1615341     7DG2_A   27.119              59         42        1
    2409 Query_1615341     3UIB_A   23.200             250        132       12
    2410 Query_1615341     3PG1_A   23.200             250        132       12
    2411 Query_1615341     2ENN_A   28.000              50         32        1
    2412 Query_1615341     2E73_A   32.000              50         30        1
    2413 Query_1615341     2DB6_A   29.167              48         30        1
    2414 Query_1615341     9C5F_A   45.833              24         13        0
    2415 Query_1615341     5HNV_A   19.915             236        138        6
    2416 Query_1615341     4IX3_A   23.214             168         98        7
    1773 Query_1615341     9EJK_B   30.000              50         31        1
    2417 Query_1615341     4B6D_A   32.000              50         31        1
         q.start q.end s.start s.end    evalue bitscore positives mlog.evalue
    1          1   766      28   793  0.00e+00   1591.0     99.87 709.1962086
    2          1   766       1   766  0.00e+00   1590.0     99.74 709.1962086
    3          1   766      28   793  0.00e+00   1590.0     99.87 709.1962086
    4          1   766      28   793  0.00e+00   1589.0     99.74 709.1962086
    5          1   766       2   767  0.00e+00   1589.0     99.74 709.1962086
    6          1   766      63   828  0.00e+00   1588.0     99.87 709.1962086
    7          1   766       3   768  0.00e+00   1587.0     99.74 709.1962086
    8          1   766       2   767  0.00e+00   1587.0     99.74 709.1962086
    9          1   766       2   767  0.00e+00   1586.0     99.61 709.1962086
    10         1   766       2   767  0.00e+00   1586.0     99.61 709.1962086
    11         1   766       2   767  0.00e+00   1583.0     99.48 709.1962086
    12       360   734       1   375  0.00e+00    774.0     99.20 709.1962086
    13       360   734       1   375  0.00e+00    774.0     99.20 709.1962086
    14        97   755      27   669  0.00e+00    711.0     68.65 709.1962086
    15        97   755      27   669  0.00e+00    709.0     68.50 709.1962086
    16        97   755      27   669  0.00e+00    705.0     68.35 709.1962086
    17       157   755      58   647  0.00e+00    703.0     72.27 709.1962086
    18        97   755       6   648  0.00e+00    701.0     68.35 709.1962086
    19        97   755      27   669  0.00e+00    701.0     68.05 709.1962086
    20       157   760      58   652  0.00e+00    700.0     71.84 709.1962086
    21       432   736      25   329  0.00e+00    637.0     99.67 709.1962086
    22       432   726      13   307  0.00e+00    621.0    100.00 709.1962086
    23       432   726      12   306  0.00e+00    621.0    100.00 709.1962086
    24       432   726      13   307  0.00e+00    621.0    100.00 709.1962086
    25       432   726      13   307  0.00e+00    621.0    100.00 709.1962086
    26       433   726       7   300  0.00e+00    620.0    100.00 709.1962086
    27       432   726      13   307  0.00e+00    618.0     99.66 709.1962086
    28       433   727       6   300  0.00e+00    616.0     99.32 709.1962086
    29       432   723      13   304  0.00e+00    614.0    100.00 709.1962086
    30       432   726      13   302  0.00e+00    604.0     98.31 709.1962086
    31       445   726       3   284  0.00e+00    595.0    100.00 709.1962086
    32       446   727       1   282  0.00e+00    594.0    100.00 709.1962086
    33       442   723       2   283  0.00e+00    590.0     99.65 709.1962086
    34       445   723       6   284  0.00e+00    588.0    100.00 709.1962086
    35       445   723       3   281  0.00e+00    588.0    100.00 709.1962086
    36       445   723       2   280  0.00e+00    587.0    100.00 709.1962086
    37       448   722       1   275  0.00e+00    580.0    100.00 709.1962086
    38       448   723       1   276  0.00e+00    580.0     99.64 709.1962086
    39       448   723       1   276  0.00e+00    579.0     99.64 709.1962086
    40       445   723      16   294  0.00e+00    578.0     98.92 709.1962086
    41       449   721       1   273  0.00e+00    576.0    100.00 709.1962086
    42       447   735       2   290  0.00e+00    572.0     95.50 709.1962086
    43       449   720       1   272  0.00e+00    563.0     98.16 709.1962086
    44       442   721      11   290  0.00e+00    554.0     95.36 709.1962086
    45       444   721       5   282  0.00e+00    552.0     95.68 709.1962086
    46       442   721       1   280  0.00e+00    552.0     95.36 709.1962086
    47       444   721       1   278  0.00e+00    552.0     95.68 709.1962086
    48       443   721       1   279  0.00e+00    552.0     95.34 709.1962086
    49       444   721      21   298  0.00e+00    551.0     95.68 709.1962086
    50       442   721      11   290  0.00e+00    551.0     95.00 709.1962086
    51       442   721       2   281  0.00e+00    551.0     95.36 709.1962086
    52       442   721       7   286  0.00e+00    550.0     95.00 709.1962086
    53       444   721       3   280  0.00e+00    550.0     95.32 709.1962086
    54       442   721      11   290  0.00e+00    550.0     95.00 709.1962086
    55       442   721       2   281  0.00e+00    548.0     95.00 709.1962086
    56       442   719      11   288  0.00e+00    546.0     94.96 709.1962086
    57       443   721      17   295  0.00e+00    544.0     94.27 709.1962086
    58       448   721       5   278  0.00e+00    542.0     95.26 709.1962086
    59       443   721      17   295  0.00e+00    541.0     93.91 709.1962086
    60       448   721       3   276  0.00e+00    541.0     95.26 709.1962086
    61       448   719       2   273  0.00e+00    540.0     95.59 709.1962086
    62       442   721       7   281  0.00e+00    536.0     93.57 709.1962086
    63       418   755     224   562  0.00e+00    532.0     86.26 709.1962086
    64       418   755     224   562  0.00e+00    531.0     85.96 709.1962086
    65       432   726      13   307 1.70e-171    496.0     89.15 393.2114227
    66       445   723       2   280 4.39e-169    489.0     89.96 387.6575515
    67       445   723       2   280 6.14e-169    489.0     89.96 387.3220560
    68       442   726       1   285 9.84e-169    489.0     89.12 386.8504250
    69       445   723       2   280 6.04e-163    473.0     90.32 373.5229661
    70       512   717       1   206 1.40e-124    372.0     92.72 285.1840793
    71       157   283       3   132  1.19e-54    186.0     76.92 124.1656417
    72       157   283       4   133  1.23e-54    186.0     76.92 124.1325809
    73       157   283       7   136  1.40e-54    186.0     76.92 124.0031228
    74       449   731      27   310  2.81e-54    191.0     59.04 123.3064105
    75       449   731      27   310  2.95e-54    191.0     59.04 123.2577899
    76       446   731      18   305  3.60e-54    191.0     59.66 123.0586612
    77       149   232      12    95  4.11e-54    182.0    100.00 122.9261720
    78       449   731      50   333  5.21e-54    191.0     59.04 122.6890152
    79       153   237      12    96  5.75e-54    182.0     97.65 122.5903952
    80       157   283       7   136  4.92e-53    181.0     76.15 120.4437014
    81       137   232       6    92  1.58e-52    178.0     89.58 119.2770000
    82       449   731      40   325  3.19e-48    175.0     55.93 109.3640635
    83        21   114       1    88  5.10e-48    166.0     89.36 108.8948439
    84       449   719       6   283  2.53e-46    168.0     56.18 104.9906950
    85       448   727      30   307  1.11e-45    167.0     59.93 103.5119692
    86       448   717       4   267  1.57e-45    166.0     55.20 103.1652536
    87       448   727      30   307  3.58e-45    166.0     59.93 102.3409664
    88        38   110     160   232  2.82e-42    155.0    100.00  95.6718370
    89       446   716      25   277  3.74e-42    157.0     54.95  95.3894883
    90       446   716       3   255  4.05e-42    157.0     54.95  95.3098570
    91       446   715       9   270  4.63e-40    150.0     54.51  90.5708469
    92       446   715       6   267  1.25e-39    149.0     54.15  89.5776751
    93       446   715       2   263  1.30e-39    149.0     54.15  89.5384544
    94       446   715      10   271  1.53e-39    149.0     54.15  89.3755509
    95       446   715      16   277  1.78e-39    149.0     54.15  89.2242053
    96       446   715       9   270  1.85e-39    149.0     54.15  89.1856330
    97       448   715       1   260  2.19e-39    148.0     54.18  89.0169171
    98       448   715       2   261  2.30e-39    148.0     54.18  88.9679095
    99       448   707       6   258  3.07e-39    148.0     54.48  88.6791411
    100      446   715       9   270  3.19e-39    148.0     54.15  88.6407977
    101      446   715       6   267  3.77e-39    148.0     54.15  88.4737436
    102      446   715       9   270  4.40e-39    148.0     53.79  88.3192141
    103      446   715       9   270  4.67e-39    147.0     53.79  88.2596596
    104      446   715       9   270  5.05e-39    147.0     53.79  88.1814304
    105      448   715       1   260  5.07e-39    147.0     53.82  88.1774778
    106      448   715       1   260  5.19e-39    147.0     57.04  88.1540849
    107      446   715       6   267  5.25e-39    147.0     53.99  88.1425906
    108      447   715       1   261  5.36e-39    147.0     53.62  88.1218547
    109      446   715       9   270  5.95e-39    147.0     53.79  88.0174274
    110      447   715      15   275  7.29e-39    147.0     56.83  87.8143151
    111      447   715      10   270  7.39e-39    147.0     56.83  87.8006909
    112      447   715      11   271  7.61e-39    147.0     56.83  87.7713555
    113      450   707       7   257  7.64e-39    147.0     54.51  87.7674210
    114      447   715       5   265  9.43e-39    146.0     56.83  87.5569225
    115      447   715       1   261  1.11e-38    146.0     53.99  87.3938735
    116      446   715       9   270  1.34e-38    146.0     53.43  87.2055639
    117      447   715      13   273  1.76e-38    146.0     56.73  86.9329197
    118      446   715       9   270  1.80e-38    146.0     53.43  86.9104469
    119      446   715       9   270  1.80e-38    146.0     53.43  86.9104469
    120      446   715       9   270  1.91e-38    146.0     53.43  86.8511303
    121      447   715      11   271  2.08e-38    146.0     56.73  86.7658656
    122      447   715       5   265  2.41e-38    145.0     56.73  86.6186068
    123      446   715       9   270  2.46e-38    145.0     53.43  86.5980722
    124      447   715       7   267  2.51e-38    145.0     56.73  86.5779508
    125      447   715      14   274  2.52e-38    145.0     56.73  86.5739746
    126      447   715      14   274  2.59e-38    145.0     56.73  86.5465757
    127      450   715       2   259  2.66e-38    145.0     53.85  86.5199074
    128      446   715       1   262  2.72e-38    145.0     53.43  86.4976017
    129      446   715       9   270  2.91e-38    145.0     53.43  86.4300805
    130      447   715       5   265  2.95e-38    145.0     56.73  86.4164284
    131      447   715       6   266  3.02e-38    145.0     56.73  86.3929767
    132      447   715       1   261  4.02e-38    145.0     53.26  86.1069516
    133      446   715     258   519  4.35e-38    151.0     54.15  86.0280577
    134      446   715     259   520  4.66e-38    151.0     54.15  85.9592181
    135      446   715     174   435  5.44e-38    149.0     54.15  85.8044545
    136      450   707       7   257  5.45e-38    144.0     54.51  85.8026179
    137      446   715     175   436  5.75e-38    149.0     54.15  85.7490337
    138      448   710       8   262  5.90e-38    144.0     53.70  85.7232812
    139      448   710      19   273  6.48e-38    144.0     53.70  85.6295130
    140      446   715     175   436  6.56e-38    149.0     54.15  85.6172429
    141      446   715       9   270  1.12e-37    144.0     53.07  85.0823198
    142      446   715     200   461  1.12e-37    149.0     54.15  85.0823198
    143      446   715       9   270  1.40e-37    143.0     53.45  84.8591762
    144      449   715       3   261  1.40e-37    143.0     56.88  84.8591762
    145      446   715     174   435  1.47e-37    148.0     54.15  84.8103860
    146      446   715     178   439  1.49e-37    148.0     54.15  84.7968723
    147      446   715     173   434  3.83e-37    146.0     53.79  83.8527836
    148      446   715     175   436  4.26e-37    146.0     53.79  83.7463793
    149      447   715       5   265  6.11e-37    141.0     56.00  83.3857217
    150      446   715     176   437  8.24e-37    145.0     53.43  83.0866481
    151      451   732       4   288  9.36e-37    142.0     54.24  82.9592032
    152      451   732      41   325  9.47e-37    144.0     54.24  82.9475195
    153      451   732      12   296  1.09e-36    142.0     53.72  82.8068857
    154      451   732      29   313  1.78e-36    142.0     54.24  82.3164500
    155      448   715     175   424  1.85e-36    144.0     53.68  82.2778777
    156      451   732      41   325  1.85e-36    143.0     54.24  82.2778777
    157      451   732      41   325  2.61e-36    142.0     54.24  81.9337131
    158      451   732      41   325  2.77e-36    142.0     54.24  81.8742160
    159      451   732      41   325  2.77e-36    142.0     54.24  81.8742160
    160      451   732      39   323  2.84e-36    142.0     54.24  81.8492593
    161      451   732      41   325  3.02e-36    142.0     54.24  81.7878065
    162      448   715     181   440  3.06e-36    144.0     53.82  81.7746484
    163      451   732      41   325  3.13e-36    142.0     53.90  81.7520303
    164      448   715     181   440  3.33e-36    144.0     53.82  81.6900910
    165      451   730      15   296  4.41e-36    139.0     54.64  81.4091887
    166      451   732      12   296  4.43e-36    141.0     53.38  81.4046638
    167      451   730      14   295  6.91e-36    139.0     54.64  80.9600937
    168      451   732      41   325  8.88e-36    140.0     53.90  80.7092618
    169      451   732      41   325  1.16e-35    140.0     53.90  80.4420582
    170      447   710       2   258  1.74e-35    137.0     54.10  80.0365931
    171      451   732       4   285  1.99e-35    137.0     54.11  79.9023436
    172      451   732      10   291  2.08e-35    137.0     54.11  79.8581104
    173      447   710       2   258  2.11e-35    136.0     54.10  79.8437903
    174      451   732      25   306  2.18e-35    138.0     54.11  79.8111534
    175      442   729      36   323  2.30e-35    138.0     51.52  79.7575691
    176      451   720      12   283  3.44e-35    137.0     55.16  79.3550068
    177      157   237      14    94  4.52e-35    129.0     83.95  79.0819663
    178      442   729      36   323  4.90e-35    137.0     51.34  79.0012430
    179      451   730      15   296  5.21e-35    136.0     53.95  78.9398984
    180      451   730      10   291  5.39e-35    136.0     53.26  78.9059329
    181      451   720      22   293  5.79e-35    136.0     54.80  78.8343460
    182      451   720      22   293  8.27e-35    136.0     54.80  78.4778437
    183      439   715     176   442  9.64e-35    139.0     52.82  78.3245571
    184      439   715     176   442  1.01e-34    139.0     52.82  78.2779428
    185      451   720       3   274  1.34e-34    135.0     54.80  77.9952235
    186      462   715       9   254  1.38e-34    134.0     53.64  77.9658097
    187      451   720      15   286  1.42e-34    135.0     54.80  77.9372363
    188      451   720      18   289  1.46e-34    135.0     55.67  77.9094567
    189      451   720      18   287  2.66e-34    134.0     53.76  77.3095670
    190      451   732      15   298  4.58e-34    134.0     54.76  76.7661942
    191      451   720      15   286  6.78e-34    133.0     53.38  76.3739161
    192      451   720      29   300  1.79e-33    132.0     53.38  75.4030924
    193      451   710       2   254  1.98e-33    131.0     53.79  75.3022112
    194      451   710       2   254  2.08e-33    131.0     53.79  75.2529402
    195      448   715       3   280  2.33e-33    132.0     50.87  75.1394398
    196      451   710       2   254  2.55e-33    130.0     53.79  75.0492147
    197      448   715      12   289  2.83e-33    132.0     50.87  74.9450314
    198      451   720      13   283  3.68e-33    131.0     53.57  74.6823953
    199      451   720      12   282  3.92e-33    131.0     53.57  74.6192164
    200      448   715      12   289  5.33e-33    131.0     51.19  74.3119568
    201      448   715       9   286  5.86e-33    130.0     51.19  74.2171585
    202      448   715      12   289  6.57e-33    130.0     51.19  74.1027942
    203      448   715       5   282  6.60e-33    130.0     51.19  74.0982384
    204      451   720      40   310  7.23e-33    131.0     53.57  74.0070690
    205      448   715      18   295  8.14e-33    131.0     51.19  73.8885179
    206      448   715      11   288  8.48e-33    130.0     51.19  73.8475976
    207      451   719      30   299  8.65e-33    130.0     52.67  73.8277487
    208      448   715       6   283  9.00e-33    130.0     51.19  73.7880835
    209      448   715      12   289  9.19e-33    130.0     51.19  73.7671921
    210      448   715      11   288  9.48e-33    130.0     51.19  73.7361238
    211      448   715      18   295  9.64e-33    130.0     51.19  73.7193870
    212      448   715      12   289  9.73e-33    130.0     51.19  73.7100942
    213      448   715       5   282  1.09e-32    130.0     51.19  73.5965453
    214      451   719      30   299  1.22e-32    130.0     53.41  73.4838721
    215      451   719      39   308  1.40e-32    130.0     53.41  73.3462507
    216      440   726      97   375  1.42e-32    131.0     50.17  73.3320661
    217      451   720      21   291  2.71e-32    131.0     53.57  72.6857743
    218      440   726      97   375  2.86e-32    130.0     49.83  72.6319014
    219      451   720      21   291  2.89e-32    131.0     53.57  72.6214665
    220      451   720      21   291  3.00e-32    131.0     53.57  72.5841107
    221      448   715      40   317  3.01e-32    129.0     51.19  72.5807829
    222      440   726      99   377  3.06e-32    130.0     49.83  72.5643081
    223      457   720      35   298  6.32e-32    128.0     53.65  71.8390038
    224      461   722      21   271  6.69e-32    128.0     51.84  71.7821091
    225      448   715       8   285  6.82e-32    127.0     50.51  71.7628635
    226      461   722      18   268  7.15e-32    127.0     51.84  71.7156106
    227      448   715      18   295  7.44e-32    128.0     50.51  71.6758521
    228      451   720      29   298  7.67e-32    127.0     53.21  71.6454064
    229      457   720      45   308  7.73e-32    128.0     53.65  71.6376141
    230      451   720      29   298  8.02e-32    127.0     53.21  71.6007846
    231      448   730      32   325  8.45e-32    128.0     50.82  71.5485565
    232      448   730      10   303  1.04e-31    127.0     50.82  71.3409172
    233      440   726      99   377  1.06e-31    129.0     48.81  71.3218690
    234      460   710      27   288  1.40e-31    127.0     51.31  71.0436656
    235      456   714      27   287  1.45e-31    127.0     51.58  71.0085743
    236      460   710      27   288  1.58e-31    126.0     51.31  70.9227130
    237      448   711       7   280  1.91e-31    126.0     50.88  70.7330346
    238      460   710      27   288  1.99e-31    126.0     51.31  70.6920032
    239      450   714       9   265  2.38e-31    125.0     52.38  70.5130374
    240      448   711      10   283  2.42e-31    126.0     50.88  70.4963703
    241      448   715     970  1247  2.76e-31    133.0     51.03  70.3649072
    242      447   714       2   261  2.83e-31    125.0     51.81  70.3398612
    243      450   714      12   268  3.24e-31    125.0     52.01  70.2045646
    244      446   717       2   267  3.31e-31    125.0     50.53  70.1831897
    245      448   730      10   303  3.32e-31    125.0     50.00  70.1801731
    246      421   714       1   285  3.43e-31    125.0     50.83  70.1475776
    247      448   730      10   303  3.48e-31    125.0     51.13  70.1331056
    248      448   717       1   264  3.61e-31    124.0     50.54  70.0964301
    249      448   730      10   303  3.65e-31    125.0     50.49  70.0854107
    250      441   717       3   273  3.72e-31    125.0     50.35  70.0664142
    251      450   714      19   275  3.75e-31    125.0     52.01  70.0583820
    252      450   714      16   272  3.77e-31    125.0     52.01  70.0530629
    253      435   717     105   381  3.94e-31    128.0     50.00  70.0089572
    254      450   714      13   269  4.10e-31    124.0     52.01  69.9691509
    255      450   714      14   270  4.16e-31    124.0     52.01  69.9546228
    256      450   714      19   275  4.54e-31    125.0     52.01  69.8672109
    257      450   714       2   258  4.54e-31    124.0     52.01  69.8672109
    258      450   714      11   267  4.63e-31    124.0     52.01  69.8475810
    259      450   714       3   259  4.65e-31    124.0     52.01  69.8432707
    260      450   714      12   268  4.76e-31    124.0     52.01  69.8198902
    261      450   714      15   271  4.86e-31    124.0     52.01  69.7990994
    262      431   714       2   275  5.09e-31    124.0     50.68  69.7528601
    263      450   714       4   260  5.42e-31    124.0     52.01  69.6900421
    264      450   714       3   259  5.66e-31    124.0     52.01  69.6467140
    265      450   714      10   266  5.71e-31    124.0     52.01  69.6379189
    266      442   717       1   270  5.94e-31    124.0     50.18  69.5984287
    267      448   730      10   303  5.99e-31    125.0     50.81  69.5900465
    268      450   714       1   257  6.17e-31    124.0     52.01  69.5604390
    269      450   714       5   261  6.32e-31    124.0     52.01  69.5364187
    270      450   714       4   260  6.32e-31    124.0     52.01  69.5364187
    271      450   714       4   260  6.42e-31    124.0     52.01  69.5207198
    272      450   714       6   262  6.60e-31    124.0     52.01  69.4930682
    273      450   714       7   263  6.93e-31    124.0     52.38  69.4442781
    274      450   714       7   263  7.13e-31    124.0     52.01  69.4158266
    275      446   717      23   288  7.21e-31    124.0     49.47  69.4046689
    276      448   717      11   274  7.24e-31    124.0     50.54  69.4005167
    277      450   714       7   263  7.37e-31    124.0     52.01  69.3827202
    278      448   717       1   264  7.41e-31    124.0     49.46  69.3773074
    279      450   714       7   263  7.48e-31    124.0     52.01  69.3679051
    280      450   714       9   265  7.80e-31    124.0     52.01  69.3260141
    281      450   714      14   270  8.20e-31    124.0     52.01  69.2760037
    282      448   717       5   268  8.26e-31    124.0     50.54  69.2687133
    283      442   717       1   270  8.33e-31    124.0     50.18  69.2602744
    284      450   714       4   260  8.51e-31    124.0     52.01  69.2388959
    285      448   717      11   274  8.65e-31    124.0     50.54  69.2225786
    286      435   717     158   434  8.88e-31    127.0     50.00  69.1963363
    287      448   717      10   273  8.90e-31    124.0     50.54  69.1940866
    288      460   710      36   297  9.05e-31    124.0     50.94  69.1773731
    289      450   714       4   260  9.27e-31    123.0     52.01  69.1533545
    290      448   717      11   274  9.28e-31    124.0     50.54  69.1522763
    291      450   714       4   260  9.45e-31    123.0     52.01  69.1341231
    292      448   717       6   269  9.98e-31    124.0     50.54  69.0795548
    293      460   710      31   292  1.00e-30    124.0     50.94  69.0775528
    294      448   717       6   269  1.02e-30    124.0     50.54  69.0577502
    295      460   710      30   291  1.04e-30    124.0     50.94  69.0383321
    296      460   710      41   302  1.04e-30    124.0     50.94  69.0383321
    297      450   714       4   260  1.04e-30    124.0     52.01  69.0383321
    298      448   717       5   268  1.05e-30    124.0     50.54  69.0287626
    299      448   717       5   268  1.05e-30    124.0     50.54  69.0287626
    300      450   714       2   258  1.06e-30    123.0     52.75  69.0192839
    301      460   710      36   297  1.06e-30    124.0     50.94  69.0192839
    302      460   710      37   298  1.10e-30    124.0     50.94  68.9822426
    303      448   717       8   271  1.10e-30    124.0     50.54  68.9822426
    304      448   717       5   268  1.10e-30    123.0     50.54  68.9822426
    305      460   710      37   298  1.14e-30    124.0     50.94  68.9465245
    306      450   714      26   282  1.17e-30    124.0     52.01  68.9205490
    307      448   717       6   269  1.17e-30    123.0     50.54  68.9205490
    308      460   710      33   294  1.21e-30    124.0     50.94  68.8869324
    309      460   710      34   295  1.21e-30    124.0     50.94  68.8869324
    310      448   717       6   269  1.25e-30    123.0     50.54  68.8544092
    311      448   717       8   271  1.28e-30    123.0     50.54  68.8306927
    312      450   714       4   260  1.32e-30    123.0     52.01  68.7999211
    313      448   730      10   303  1.32e-30    124.0     51.13  68.7999211
    314      448   717       8   271  1.33e-30    123.0     50.54  68.7923738
    315      448   730      11   304  1.34e-30    124.0     51.47  68.7848832
    316      460   710      34   295  1.36e-30    124.0     50.94  68.7700681
    317      460   705      27   283  1.40e-30    124.0     51.15  68.7410806
    318      451   715       2   258  1.41e-30    123.0     51.82  68.7339631
    319      441   717       3   273  1.44e-30    123.0     50.00  68.7129097
    320      446   717       2   267  1.47e-30    123.0     50.18  68.6922904
    321      460   705      26   282  1.49e-30    124.0     51.15  68.6787767
    322      461   722      20   270  1.62e-30    130.0     51.84  68.5951266
    323      451   715       4   260  1.72e-30    122.0     51.82  68.5352285
    324      460   710      51   312  1.79e-30    124.0     50.94  68.4953372
    325      460   710      32   293  1.80e-30    124.0     50.56  68.4897661
    326      460   710      36   297  1.80e-30    124.0     50.56  68.4897661
    327      460   710      30   291  1.84e-30    123.0     50.56  68.4677872
    328      460   710      52   313  1.89e-30    124.0     50.94  68.4409760
    329      448   717      19   282  1.91e-30    123.0     50.54  68.4304495
    330      460   710      31   292  1.91e-30    123.0     50.56  68.4304495
    331      448   717      33   296  1.95e-30    124.0     50.54  68.4097234
    332      460   710      29   290  2.00e-30    123.0     50.56  68.3844056
    333      460   710      37   298  2.02e-30    124.0     50.56  68.3744553
    334      460   710      33   294  2.06e-30    123.0     50.56  68.3548468
    335      460   710      34   295  2.06e-30    123.0     50.56  68.3548468
    336      460   710      28   289  2.21e-30    123.0     50.56  68.2845603
    337      460   705      26   282  2.26e-30    123.0     50.76  68.2621880
    338      451   715      14   270  2.33e-30    122.0     51.82  68.2316845
    339      460   705      26   282  2.47e-30    123.0     50.76  68.1733346
    340      448   715     990  1267  2.73e-30    130.0     51.19  68.0732512
    341      448   730      12   305  2.74e-30    123.0     50.81  68.0695949
    342      435   717     200   476  2.75e-30    127.0     50.00  68.0659519
    343      448   730      11   304  2.79e-30    123.0     50.81  68.0515112
    344      460   705      26   282  2.83e-30    122.0     50.76  68.0372761
    345      450   714     183   439  2.92e-30    126.0     52.75  68.0059692
    346      450   714      19   275  3.07e-30    122.0     51.65  67.9558752
    347      460   710      52   313  3.25e-30    123.0     50.56  67.8988978
    348      460   710      52   313  3.32e-30    123.0     50.56  67.8775880
    349      448   717      11   274  3.38e-30    122.0     50.18  67.8596771
    350      460   710      46   307  3.43e-30    123.0     50.56  67.8449925
    351      460   705      56   312  3.43e-30    123.0     50.76  67.8449925
    352      448   730    1004  1297  3.74e-30    130.0     50.65  67.7584672
    353      448   717       5   268  3.87e-30    122.0     50.18  67.7242983
    354      435   717     239   515  4.08e-30    127.0     50.00  67.6714558
    355      448   717       5   268  4.35e-30    121.0     50.18  67.6073769
    356      448   715     965  1243  4.41e-30    129.0     51.03  67.5936781
    357      448   730     977  1270  4.47e-30    129.0     50.65  67.5801644
    358      448   717      11   274  4.70e-30    121.0     50.18  67.5299903
    359      435   717     197   473  4.72e-30    126.0     50.00  67.5257440
    360      456   715      43   301  4.84e-30    122.0     52.04  67.5006381
    361      460   710      41   302  4.89e-30    122.0     50.94  67.4903605
    362      456   715      62   320  4.90e-30    123.0     52.42  67.4883176
    363      450   714       2   258  5.11e-30    121.0     51.65  67.4463534
    364      456   715      62   320  5.18e-30    123.0     52.04  67.4327477
    365      460   710      31   292  5.21e-30    122.0     50.94  67.4269729
    366      450   714      15   271  5.37e-30    121.0     51.65  67.3967249
    367      460   710      26   287  5.38e-30    122.0     50.94  67.3948644
    368      448   730     977  1270  5.41e-30    129.0     50.65  67.3893037
    369      450   714     182   438  5.48e-30    125.0     52.01  67.3764477
    370      448   717       6   269  5.54e-30    122.0     50.18  67.3655583
    371      448   717       6   269  6.25e-30    121.0     50.18  67.2449713
    372      460   705     189   445  6.74e-30    125.0     51.53  67.1694929
    373      448   717       1   264  7.06e-30    121.0     50.18  67.1231077
    374      435   717     200   476  1.10e-29    125.0     49.66  66.6796575
    375      448   715     961  1239  1.39e-29    128.0     51.03  66.4456639
    376      451   714       4   259  1.46e-29    120.0     51.84  66.3965313
    377      446   717       2   267  1.71e-29    120.0     49.82  66.2384743
    378      451   714       2   257  1.85e-29    119.0     51.84  66.1597821
    379      448   717       7   270  2.18e-29    120.0     49.82  65.9956428
    380      450   715      18   293  2.56e-29    120.0     48.25  65.8349604
    381      450   715      41   326  3.06e-29    120.0     46.96  65.6565528
    382      450   715      41   326  3.21e-29    120.0     46.96  65.6086968
    383      451   714       2   257  3.23e-29    119.0     52.57  65.6024856
    384      457   714       6   255  4.20e-29    118.0     52.26  65.3398832
    385      450   714     395   651  4.57e-29    125.0     53.11  65.2554545
    386      450   715      18   292  6.20e-29    119.0     48.07  64.9504184
    387      448   717       5   268  8.59e-29    118.0     48.91  64.6243690
    388      448   730     987  1280  1.01e-28    125.0     51.14  64.4624323
    389      448   730    1014  1307  1.02e-28    125.0     51.14  64.4525800
    390      448   730    1002  1295  1.07e-28    125.0     51.14  64.4047240
    391      448   730    1002  1295  1.09e-28    125.0     51.14  64.3862049
    392      450   715      18   307  2.03e-28    118.0     46.33  63.7643468
    393      451   711      12   282  2.04e-28    118.0     50.70  63.7594328
    394      451   714     194   449  2.51e-28    120.0     52.94  63.5520999
    395      450   715      18   307  2.55e-28    117.0     46.33  63.5362892
    396      450   715      41   326  2.74e-28    118.0     46.62  63.4644247
    397      448   730     975  1268  3.05e-28    124.0     51.14  63.3572410
    398      450   715      18   307  5.14e-28    117.0     46.00  62.8353295
    399      446   714       4   269  6.85e-28    115.0     49.29  62.5481340
    400      446   714       2   267  7.00e-28    115.0     49.29  62.5264725
    401      450   715      18   307  7.16e-28    116.0     46.00  62.5038726
    402      446   714       4   269  7.67e-28    115.0     49.29  62.4350660
    403      446   714       3   268  7.72e-28    115.0     49.29  62.4285682
    404      446   714       3   268  7.80e-28    115.0     49.29  62.4182589
    405      450   715      18   307  9.60e-28    116.0     46.00  62.2106195
    406      460   710      27   288  9.83e-28    115.0     49.06  62.1869437
    407      450   715      18   307  9.97e-28    116.0     46.00  62.1728020
    408      446   714       6   271  1.03e-27    115.0     49.29  62.1402387
    409      463   719      46   301  1.06e-27    116.0     49.43  62.1115286
    410      446   714       1   266  1.08e-27    115.0     49.29  62.0928365
    411      446   714       1   266  1.11e-27    115.0     49.29  62.0654375
    412      463   720      23   279  1.46e-27    115.0     49.25  61.7913611
    413      446   714      29   294  1.47e-27    115.0     49.29  61.7845351
    414      450   715      36   319  1.86e-27    115.0     47.46  61.5492210
    415      450   715      35   318  1.86e-27    115.0     47.46  61.5492210
    416      452   707      15   272  1.87e-27    114.0     47.58  61.5438591
    417      451   710       5   280  1.91e-27    114.0     48.44  61.5226943
    418      452   707       5   262  2.10e-27    114.0     47.58  61.4278602
    419      450   715      36   319  2.12e-27    115.0     47.46  61.4183814
    420      452   707      15   272  2.25e-27    114.0     47.58  61.3588673
    421      452   707      15   272  2.43e-27    114.0     47.58  61.2819063
    422      450   720      12   294  2.60e-27    114.0     47.65  61.2142861
    423      452   707      17   274  2.97e-27    114.0     47.58  61.0812356
    424      449   714       1   263  3.03e-27    113.0     49.10  61.0612349
    425      451   710      22   297  3.06e-27    114.0     48.44  61.0513826
    426      452   707       5   262  3.08e-27    113.0     47.58  61.0448679
    427      446   714       1   266  3.12e-27    113.0     48.93  61.0319645
    428      450   715      35   318  3.16e-27    115.0     46.78  61.0192255
    429      452   707      17   274  3.18e-27    114.0     47.58  61.0129163
    430      452   707       5   262  3.24e-27    113.0     47.58  60.9942242
    431      446   714       1   266  3.63e-27    113.0     48.93  60.8805649
    432      452   707       5   262  3.63e-27    113.0     47.58  60.8805649
    433      452   707       9   266  3.64e-27    113.0     47.58  60.8778138
    434      450   715      59   348  3.96e-27    115.0     46.36  60.7935535
    435      450   715      18   307  4.01e-27    114.0     45.67  60.7810063
    436      444   715      12   294  5.99e-27    114.0     48.64  60.3797061
    437      457   714      41   300  7.89e-27    113.0     49.45  60.1042014
    438      457   714      28   287  7.90e-27    112.0     49.45  60.1029348
    439      448   715    1362  1651  8.19e-27    119.0     46.84  60.0668836
    440      452   707       9   266  8.44e-27    114.0     47.58  60.0368152
    441      450   715      38   327  8.64e-27    114.0     46.36  60.0133949
    442      457   714      28   287  9.33e-27    112.0     49.45  59.9365625
    443      450   715      15   304  1.00e-26    113.0     45.67  59.8672124
    444      378   714     309   646  1.02e-26    117.0     46.18  59.8474098
    445      450   720      18   309  1.05e-26    113.0     46.58  59.8184223
    446      450   720      16   307  1.15e-26    112.0     46.58  59.7274505
    447      450   720       7   298  1.20e-26    112.0     46.58  59.6848909
    448      450   720      16   307  1.34e-26    112.0     46.58  59.5745428
    449      450   720      18   309  1.36e-26    112.0     46.58  59.5597277
    450      450   720       8   299  1.40e-26    112.0     46.25  59.5307402
    451      450   720      18   309  1.59e-26    112.0     46.25  59.4034784
    452      450   720      14   305  1.60e-26    112.0     46.25  59.3972088
    453      450   720      30   321  1.61e-26    112.0     46.58  59.3909782
    454      450   720      14   305  1.65e-26    112.0     46.25  59.3664371
    455      450   720      18   309  1.65e-26    112.0     46.58  59.3664371
    456      450   715       7   260  1.66e-26    111.0     48.72  59.3603948
    457      450   720       8   299  1.67e-26    112.0     46.25  59.3543888
    458      450   720      18   309  1.70e-26    112.0     46.25  59.3365842
    459      450   715      18   307  1.74e-26    112.0     45.67  59.3133273
    460      450   720       3   265  1.78e-26    111.0     50.00  59.2905991
    461      450   715      16   269  1.81e-26    111.0     48.72  59.2738856
    462      450   720      16   307  1.83e-26    112.0     46.58  59.2628965
    463      450   720      18   309  1.83e-26    112.0     46.58  59.2628965
    464      450   720      16   307  1.88e-26    112.0     46.58  59.2359406
    465      461   736      14   275  1.92e-26    112.0     47.37  59.2148872
    466      450   715       1   254  1.96e-26    110.0     48.72  59.1942679
    467      457   714      19   278  2.09e-26    111.0     49.08  59.1300484
    468      457   714      28   287  2.23e-26    111.0     49.08  59.0652108
    469      457   714      24   283  2.25e-26    111.0     49.08  59.0562822
    470      450   714       1   262  2.26e-26    110.0     48.91  59.0518476
    471      450   720      30   321  2.27e-26    112.0     46.25  59.0474326
    472      457   714      41   300  2.28e-26    112.0     49.08  59.0430370
    473      457   714      28   287  2.31e-26    111.0     49.08  59.0299649
    474      457   714      41   300  2.36e-26    112.0     49.08  59.0085508
    475      378   714     311   648  2.49e-26    116.0     45.89  58.9549297
    476      450   720      30   321  2.49e-26    112.0     46.58  58.9549297
    477      457   714      24   283  2.52e-26    111.0     49.08  58.9429535
    478      378   714     309   646  2.60e-26    116.0     45.89  58.9117010
    479      450   715      23   307  2.74e-26    112.0     48.30  58.8592545
    480      457   714      31   290  2.86e-26    111.0     49.08  58.8163908
    481      450   720      29   320  2.87e-26    112.0     46.25  58.8129004
    482      457   714      27   286  2.96e-26    111.0     49.08  58.7820231
    483      457   714      50   309  2.98e-26    112.0     49.08  58.7752891
    484      457   714      49   308  3.09e-26    111.0     49.08  58.7390413
    485      457   714      32   291  3.26e-26    111.0     49.08  58.6854852
    486      461   736      15   276  3.37e-26    110.0     47.37  58.6522997
    487      457   714      30   289  3.43e-26    111.0     49.08  58.6346522
    488      461   736      45   306  3.59e-26    112.0     48.40  58.5890602
    489      461   736      22   283  3.68e-26    110.0     48.40  58.5642997
    490      457   714      32   291  3.75e-26    111.0     49.08  58.5454566
    491      457   714      29   288  3.79e-26    111.0     49.08  58.5348464
    492      461   736      21   282  3.79e-26    111.0     48.40  58.5348464
    493      450   720      30   321  3.81e-26    111.0     46.25  58.5295832
    494      453   736       7   276  3.84e-26    110.0     46.76  58.5217401
    495      461   736      23   284  3.87e-26    111.0     48.40  58.5139579
    496      457   714      54   313  3.87e-26    111.0     49.08  58.5139579
    497      461   736      22   283  3.90e-26    111.0     48.40  58.5062359
    498      450   715      35   319  3.93e-26    112.0     47.96  58.4985730
    499      453   736       9   278  4.00e-26    110.0     47.75  58.4809181
    500      457   714      41   300  4.17e-26    111.0     49.08  58.4392964
    501      457   714      51   310  4.18e-26    111.0     49.08  58.4369012
    502      457   714      61   320  4.28e-26    112.0     49.08  58.4132594
    503      450   715      12   296  4.32e-26    110.0     48.30  58.4039570
    504      457   714      31   290  4.41e-26    111.0     49.08  58.3833377
    505      453   736      12   281  4.42e-26    110.0     47.75  58.3810727
    506      450   715      14   298  4.55e-26    111.0     48.30  58.3520852
    507      450   715      13   297  4.58e-26    110.0     48.30  58.3455134
    508      450   715       8   292  4.65e-26    110.0     48.30  58.3303452
    509      453   736       8   277  4.69e-26    110.0     47.75  58.3217798
    510      450   715      14   298  4.80e-26    110.0     48.30  58.2985965
    511      450   715      12   296  4.93e-26    110.0     48.30  58.2718734
    512      457   714      32   291  4.99e-26    110.0     49.08  58.2597765
    513      461   736      22   283  5.01e-26    110.0     48.40  58.2557765
    514      461   736      24   285  5.02e-26    110.0     48.40  58.2537825
    515      450   715      15   299  5.38e-26    110.0     48.30  58.1845240
    516      450   715      17   301  5.54e-26    110.0     48.30  58.1552179
    517      450   715      16   300  5.65e-26    110.0     48.30  58.1355569
    518      457   714      22   281  5.98e-26    110.0     49.08  58.0787918
    519      453   736       8   277  6.01e-26    110.0     46.76  58.0737877
    520      450   715      16   300  6.03e-26    110.0     48.30  58.0704654
    521      450   715      19   303  6.04e-26    110.0     48.30  58.0688084
    522      452   714      16   284  6.04e-26    110.0     50.18  58.0688084
    523      450   715      23   307  6.10e-26    110.0     48.30  58.0589236
    524      463   711      18   266  6.22e-26    109.0     50.95  58.0394425
    525      461   715      19   290  6.32e-26    110.0     51.06  58.0234932
    526      450   715      23   307  6.52e-26    110.0     48.30  57.9923380
    527      447   720       2   267  6.75e-26    109.0     49.47  57.9576699
    528      450   715      16   300  6.89e-26    110.0     47.96  57.9371413
    529      461   736      22   283  7.21e-26    110.0     47.37  57.8917435
    530      450   715      16   300  7.29e-26    110.0     48.30  57.8807089
    531      450   715      14   298  7.35e-26    110.0     47.96  57.8725121
    532      461   736      36   297  7.50e-26    110.0     48.40  57.8523094
    533      450   715      23   307  7.64e-26    110.0     48.30  57.8338148
    534      450   715      17   301  7.64e-26    110.0     47.96  57.8338148
    535      453   736       7   276  7.75e-26    109.0     46.76  57.8195196
    536      450   715      18   299  7.97e-26    110.0     47.44  57.7915279
    537      457   714      33   292  8.02e-26    110.0     49.45  57.7852740
    538      450   715      15   299  8.08e-26    110.0     47.96  57.7778205
    539      461   715      39   310  8.10e-26    110.0     50.53  57.7753484
    540      457   714      31   290  8.15e-26    110.0     48.71  57.7691945
    541      461   736      22   283  8.29e-26    110.0     47.37  57.7521624
    542      450   715      23   307  8.33e-26    110.0     48.30  57.7473490
    543      453   736       6   275  8.67e-26    109.0     46.76  57.7073436
    544      461   736      23   284  9.76e-26    110.0     48.04  57.5889200
    545      461   715      14   285  9.90e-26    109.0     50.53  57.5746777
    546      461   721      10   271  1.00e-25    111.0     48.71  57.5646273
    547      461   736      42   303  1.02e-25    110.0     47.37  57.5448247
    548      461   715      18   289  1.03e-25    109.0     50.53  57.5350685
    549      447   715       2   287  1.07e-25    109.0     49.49  57.4969687
    550      461   715      19   290  1.09e-25    109.0     50.53  57.4784496
    551      461   715      19   290  1.14e-25    109.0     50.71  57.4335991
    552      457   714      28   287  1.14e-25    109.0     48.71  57.4335991
    553      457   714      32   291  1.29e-25    109.0     49.08  57.3099851
    554      457   714      32   291  1.30e-25    109.0     49.08  57.3022631
    555      457   714      91   350  1.36e-25    111.0     49.08  57.2571426
    556      457   714      32   291  1.44e-25    109.0     48.71  57.1999842
    557      450   725      41   349  1.44e-25    110.0     46.08  57.1999842
    558      457   714      37   296  1.45e-25    109.0     49.08  57.1930638
    559      450   714      13   296  1.47e-25    109.0     47.96  57.1793649
    560      450   715      12   284  1.48e-25    109.0     48.58  57.1725852
    561      426   715      17   311  1.48e-25    109.0     48.11  57.1725852
    562      457   714      33   292  1.48e-25    109.0     49.08  57.1725852
    563      450   725      37   345  1.48e-25    110.0     46.08  57.1725852
    564      457   714      30   289  1.55e-25    109.0     49.08  57.1263724
    565      442   715       4   294  1.57e-25    109.0     49.01  57.1135517
    566      461   715      25   296  1.58e-25    109.0     50.71  57.1072025
    567      452   720       4   264  1.58e-25    108.0     50.00  57.1072025
    568      450   715      15   299  1.61e-25    109.0     47.96  57.0883931
    569      461   715      25   296  1.64e-25    109.0     50.71  57.0699311
    570      461   715      37   308  1.64e-25    109.0     50.53  57.0699311
    571      446   653      10   214  1.71e-25    108.0     53.52  57.0281340
    572      461   715      29   300  1.82e-25    109.0     50.53  56.9657908
    573      451   710      15   290  1.84e-25    109.0     48.26  56.9548618
    574      450   715      12   284  1.84e-25    108.0     48.58  56.9548618
    575      442   715       2   292  1.86e-25    108.0     49.01  56.9440508
    576      447   715       2   287  1.94e-25    108.0     49.49  56.9019394
    577      461   715      39   310  1.97e-25    109.0     50.71  56.8865938
    578      450   715     188   441  1.97e-25    112.0     48.72  56.8865938
    579      461   715      19   290  2.08e-25    108.0     50.18  56.8322594
    580      442   715       1   291  2.21e-25    108.0     48.68  56.7716348
    581      450   714      24   305  2.24e-25    109.0     47.95  56.7581515
    582      450   715      64   348  2.25e-25    110.0     48.30  56.7536971
    583      457   714    1078  1337  2.25e-25    114.0     49.08  56.7536971
    584      450   714      23   305  2.35e-25    109.0     47.78  56.7102120
    585      458   653      50   243  2.36e-25    109.0     53.96  56.7059657
    586      446   715      22   308  2.51e-25    108.0     49.33  56.6443446
    587      461   715      15   286  2.58e-25    108.0     50.18  56.6168379
    588      461   715      19   290  2.59e-25    108.0     49.82  56.6129694
    589      461   715      19   290  2.64e-25    108.0     50.18  56.5938484
    590      461   715      37   308  2.66e-25    108.0     50.71  56.5863012
    591      461   715      14   285  2.68e-25    108.0     50.18  56.5788105
    592      461   715      22   293  2.70e-25    108.0     50.18  56.5713756
    593      446   715      23   309  2.73e-25    108.0     49.33  56.5603257
    594      450   714      13   296  2.76e-25    108.0     47.96  56.5493966
    595      450   714      13   296  2.79e-25    108.0     47.96  56.5385857
    596      461   715      19   290  2.79e-25    108.0     49.82  56.5385857
    597      458   653      23   216  2.86e-25    108.0     53.96  56.5138057
    598      446   715       2   288  2.96e-25    108.0     49.66  56.4794381
    599      442   715       4   294  2.97e-25    108.0     48.34  56.4760654
    600      450   714      22   301  2.98e-25    108.0     48.28  56.4727040
    601      450   715      18   303  3.01e-25    108.0     46.80  56.4626872
    602      452   717      16   287  3.06e-25    108.0     49.65  56.4462124
    603      445   715       4   291  3.33e-25    108.0     49.16  56.3616550
    604      461   715      20   291  3.33e-25    108.0     50.18  56.3616550
    605      461   715      47   318  3.35e-25    108.0     50.18  56.3556670
    606      461   715      12   283  3.38e-25    108.0     50.18  56.3467516
    607      447   715      30   310  3.49e-25    108.0     45.70  56.3147256
    608      447   715      38   324  3.51e-25    108.0     45.45  56.3090113
    609      450   714      22   305  3.53e-25    108.0     47.62  56.3033295
    610      450   715      20   305  3.53e-25    108.0     46.80  56.3033295
    611      450   715      17   302  3.53e-25    108.0     46.80  56.3033295
    612      447   715      35   321  3.56e-25    108.0     45.45  56.2948668
    613      461   715      12   283  3.64e-25    107.0     50.71  56.2726436
    614      461   715      16   287  3.76e-25    108.0     50.18  56.2402084
    615      461   715      17   288  3.78e-25    108.0     50.18  56.2349033
    616      452   714      16   284  3.81e-25    107.0     49.82  56.2269981
    617      463   711      34   281  3.82e-25    107.0     50.95  56.2243769
    618      453   711      22   285  3.86e-25    108.0     50.90  56.2139601
    619      447   715      36   322  3.90e-25    108.0     45.45  56.2036508
    620      447   715      41   327  4.03e-25    108.0     45.45  56.1708609
    621      444   715      17   305  4.03e-25    108.0     49.33  56.1708609
    622      457   714      30   289  4.10e-25    108.0     49.08  56.1536404
    623      450   714      22   301  4.14e-25    108.0     48.28  56.1439315
    624      461   736      36   297  4.14e-25    108.0     46.67  56.1439315
    625      450   715      30   314  4.16e-25    108.0     48.64  56.1391123
    626      444   715      17   305  4.38e-25    108.0     49.33  56.0875786
    627      453   711      13   276  4.64e-25    108.0     50.90  56.0299130
    628      449   705      22   293  4.66e-25    108.0     50.35  56.0256119
    629      462   710      25   284  4.66e-25    107.0     46.84  56.0256119
    630      453   711      15   278  4.86e-25    108.0     50.90  55.9835889
    631      453   711      13   276  4.95e-25    107.0     50.90  55.9652397
    632      453   711      12   275  5.09e-25    107.0     50.90  55.9373495
    633      452   720       4   264  5.14e-25    107.0     49.64  55.9275742
    634      462   710      25   284  5.21e-25    107.0     46.84  55.9140475
    635      462   710      25   284  5.21e-25    107.0     46.84  55.9140475
    636      450   714      22   305  5.29e-25    108.0     47.62  55.8988091
    637      463   711      34   281  5.33e-25    107.0     50.95  55.8912761
    638      452   720       2   262  5.36e-25    106.0     49.64  55.8856633
    639      461   715      36   307  5.53e-25    108.0     50.18  55.8544395
    640      462   710      24   283  5.62e-25    107.0     46.84  55.8382957
    641      458   653      23   216  5.70e-25    107.0     53.96  55.8241612
    642      450   715      18   302  5.71e-25    107.0     48.64  55.8224083
    643      463   711      34   281  5.91e-25    107.0     50.95  55.7879815
    644      453   711      11   274  5.92e-25    107.0     50.90  55.7862909
    645      462   710      22   281  5.94e-25    107.0     46.84  55.7829182
    646      463   711      35   282  6.09e-25    107.0     50.95  55.7579792
    647      450   714      13   296  6.10e-25    107.0     47.62  55.7563386
    648      450   714      24   307  6.21e-25    107.0     47.62  55.7384664
    649      463   721      29   289  6.31e-25    107.0     47.79  55.7224916
    650      450   715      18   302  6.33e-25    107.0     50.00  55.7193271
    651      450   715      30   314  6.36e-25    108.0     48.30  55.7145989
    652      463   710      17   275  6.36e-25    107.0     47.01  55.7145989
    653      450   715      30   314  6.48e-25    108.0     48.30  55.6959068
    654      463   710      16   274  6.60e-25    107.0     47.01  55.6775577
    655      463   710      18   276  6.62e-25    107.0     47.01  55.6745320
    656      450   715      30   314  6.86e-25    108.0     48.30  55.6389199
    657      450   715     471   755  6.90e-25    112.0     47.96  55.6331059
    658      450   715      31   315  7.18e-25    107.0     50.00  55.5933279
    659      450   715      34   319  7.36e-25    108.0     46.80  55.5685674
    660      450   715      30   314  7.39e-25    107.0     48.30  55.5644996
    661      463   711      34   281  7.39e-25    107.0     50.95  55.5644996
    662      445   715      22   309  7.42e-25    107.0     49.16  55.5604483
    663      450   714      22   305  7.43e-25    107.0     47.62  55.5591015
    664      447   715      38   318  7.57e-25    107.0     45.70  55.5404343
    665      450   725      40   348  7.78e-25    108.0     45.77  55.5130710
    666      450   725      37   345  7.91e-25    108.0     45.77  55.4964995
    667      463   711      22   269  7.92e-25    106.0     50.95  55.4952361
    668      462   710      31   290  8.08e-25    107.0     46.84  55.4752355
    669      461   712      30   294  8.25e-25    108.0     51.48  55.4544141
    670      450   715      32   317  8.25e-25    107.0     47.81  55.4544141
    671      457   714      30   289  8.81e-25    107.0     48.71  55.3887399
    672      450   715      36   321  8.94e-25    107.0     46.80  55.3740917
    673      450   725      41   349  9.01e-25    108.0     45.77  55.3662923
    674      462   710      25   284  9.74e-25    107.0     46.84  55.2883862
    675      450   714      22   305  9.91e-25    107.0     47.62  55.2710830
    676      445   715       3   290  1.00e-24    107.0     49.16  55.2620422
    677      450   714      14   297  1.01e-24    107.0     47.62  55.2520919
    678      450   715      30   314  1.03e-24    107.0     48.64  55.2324834
    679      462   710      38   297  1.03e-24    107.0     46.84  55.2324834
    680      450   715      30   314  1.04e-24    107.0     48.30  55.2228215
    681      450   715      41   326  1.04e-24    107.0     46.80  55.2228215
    682      462   710      41   300  1.06e-24    107.0     46.84  55.2037733
    683      450   715      36   321  1.08e-24    107.0     46.80  55.1850812
    684      450   725      60   368  1.09e-24    108.0     45.77  55.1758645
    685      450   715      41   326  1.27e-24    107.0     46.80  55.0230253
    686      450   715      30   314  1.32e-24    107.0     48.30  54.9844105
    687      232   283       1    52  1.33e-24     98.6     82.69  54.9768633
    688      450   715      30   314  1.48e-24    107.0     48.30  54.8700001
    689      453   704      20   263  1.54e-24    106.0     50.97  54.8302598
    690      450   715       8   292  1.58e-24    106.0     47.96  54.8046174
    691      450   725      41   349  1.59e-24    107.0     45.45  54.7983082
    692      450   715      30   314  1.65e-24    107.0     48.30  54.7612669
    693      463   711      23   270  1.68e-24    105.0     49.81  54.7432484
    694      448   711      16   271  1.75e-24    106.0     49.63  54.7024264
    695      450   714      59   342  1.81e-24    107.0     47.62  54.6687154
    696      447   715      38   322  1.96e-24    107.0     46.28  54.5890978
    697      448   711      16   271  1.96e-24    106.0     49.63  54.5890978
    698      445   721       2   295  2.00e-24    105.0     48.51  54.5688951
    699      447   715      38   322  2.08e-24    106.0     47.64  54.5296743
    700      450   715      30   314  2.14e-24    106.0     48.30  54.5012364
    701      450   715      36   320  2.76e-24    106.0     48.30  54.2468116
    702      461   711      32   276  2.98e-24    105.0     48.44  54.1701189
    703      448   711      37   292  3.00e-24    106.0     49.63  54.1634299
    704      448   711      37   292  3.05e-24    106.0     49.63  54.1469006
    705      448   711      37   292  3.48e-24    105.0     49.63  54.0150099
    706      450   715      30   314  3.59e-24    105.0     48.30  53.9838900
    707      450   715      30   314  3.59e-24    105.0     48.30  53.9838900
    708      452   720       4   264  3.61e-24    104.0     49.28  53.9783345
    709      450   715      18   302  3.76e-24    105.0     47.28  53.9376233
    710      462   710      31   290  3.85e-24    105.0     46.84  53.9139691
    711      450   715      30   314  4.13e-24    105.0     48.30  53.8437648
    712      450   715      19   303  4.29e-24    105.0     47.28  53.8057555
    713      450   715      10   294  4.55e-24    105.0     47.28  53.7469150
    714      450   715      22   306  4.66e-24    105.0     47.96  53.7230268
    715      450   715      22   306  4.66e-24    105.0     47.96  53.7230268
    716      450   715      19   303  4.89e-24    105.0     47.96  53.6748499
    717      450   715      30   314  5.14e-24    105.0     47.96  53.6249892
    718      450   715      17   301  5.23e-24    105.0     47.96  53.6076310
    719      450   715      22   306  5.67e-24    105.0     47.96  53.5268531
    720      450   715      30   314  5.74e-24    105.0     47.96  53.5145830
    721      450   715      19   303  6.31e-24    104.0     47.28  53.4199066
    722      458   653      22   213  6.39e-24    104.0     55.22  53.4073080
    723      452   720      24   284  6.59e-24    104.0     49.28  53.3764889
    724      399   704     139   430  6.77e-24    107.0     49.52  53.3495411
    725      450   715      14   298  6.83e-24    104.0     47.28  53.3407176
    726      451   717       9   271  6.89e-24    103.0     48.56  53.3319711
    727      450   715      30   314  7.05e-24    105.0     47.96  53.3090146
    728      450   715      30   314  7.60e-24    105.0     47.96  53.2338940
    729      450   715      76   360  7.82e-24    105.0     47.96  53.2053577
    730      450   715      17   301  8.70e-24    104.0     47.96  53.0987192
    731      450   715      30   314  8.92e-24    104.0     47.96  53.0737463
    732      450   715      14   298  8.96e-24    104.0     47.28  53.0692720
    733      450   715      30   314  9.44e-24    104.0     47.96  53.0170863
    734      447   715      42   346  9.51e-24    105.0     45.22  53.0096984
    735      450   715      30   314  9.62e-24    104.0     48.30  52.9981980
    736      399   704     152   443  1.10e-23    107.0     49.52  52.8641470
    737      450   715      30   314  1.32e-23    104.0     47.96  52.6818254
    738      447   715      23   316  1.40e-23    103.0     43.75  52.6229849
    739      450   715      30   314  1.60e-23    103.0     47.62  52.4894535
    740      447   715      31   324  1.85e-23    103.0     43.75  52.3442715
    741      453   704      18   261  1.86e-23    103.0     51.75  52.3388807
    742      447   715      42   335  2.01e-23    103.0     45.21  52.2613224
    743      448   724      12   294  2.08e-23    102.0     48.31  52.2270892
    744      447   715      42   335  2.62e-23    103.0     43.75  51.9962828
    745      453   711       4   253  2.79e-23    102.0     51.89  51.9334155
    746      448   728      29   319  2.81e-23    103.0     47.37  51.9262727
    747      456   722      24   300  3.48e-23    102.0     47.39  51.7124248
    748      463   704      24   264  3.53e-23    102.0     51.16  51.6981593
    749      399   704     220   511  3.84e-23    105.0     49.52  51.6139848
    750      448   724      10   296  3.88e-23    102.0     47.67  51.6036220
    751      450   715      64   348  3.90e-23    104.0     47.62  51.5984806
    752      453   705      13   259  4.32e-23    102.0     49.22  51.4962017
    753      460   721      15   284  4.35e-23    102.0     48.38  51.4892813
    754      453   705      52   298  4.47e-23    103.0     49.22  51.4620687
    755      453   711      27   276  4.61e-23    102.0     52.27  51.4312293
    756      448   724      10   296  4.91e-23    102.0     47.67  51.3681832
    757      448   724      10   296  4.96e-23    102.0     47.67  51.3580514
    758      453   711      20   269  5.14e-23    102.0     52.27  51.3224041
    759      448   724       8   294  5.82e-23    101.0     47.67  51.1981569
    760      453   711      12   261  5.87e-23    101.0     52.27  51.1896025
    761      448   724      28   314  6.03e-23    102.0     47.67  51.1627101
    762      453   711      12   261  6.04e-23    101.0     52.27  51.1610531
    763      453   711      30   279  6.49e-23    101.0     52.27  51.0891946
    764      448   724      27   313  7.56e-23    101.0     47.67  50.9365859
    765      463   704      24   264  7.93e-23    100.0     51.16  50.8888041
    766      452   743      13   316  8.20e-23    102.0     46.89  50.8553230
    767      448   724      29   315  8.60e-23    101.0     47.67  50.8076949
    768      453   711      25   274  8.65e-23    101.0     51.52  50.8018978
    769      463   704      24   264  8.80e-23    100.0     51.16  50.7847054
    770      231   656      99   542  9.33e-23    105.0     41.99  50.7262221
    772      231   656      99   542  1.00e-22    105.0     41.99  50.6568720
    774      453   711      17   266  1.01e-22    100.0     49.05  50.6469217
    775      453   711       5   254  1.04e-22    101.0     51.52  50.6176513
    776      451   710      43   326  1.05e-22    102.0     44.11  50.6080819
    777      453   711       5   254  1.08e-22    101.0     51.52  50.5799110
    778      453   711      18   267  1.09e-22    100.0     49.05  50.5706943
    779      462   717      12   283  1.10e-22    100.0     46.62  50.5615619
    780      231   656      99   542  1.20e-22    105.0     42.33  50.4745505
    782      462   717      20   291  1.36e-22    100.0     46.62  50.3493873
    783      463   719      26   281  1.41e-22    101.0     47.55  50.3132823
    784      463   719      27   282  1.42e-22    101.0     47.55  50.3062152
    785      460   721      19   288  1.46e-22    100.0     48.38  50.2784356
    786      463   711      39   302  1.50e-22    100.0     45.65  50.2514069
    787      463   719      27   282  1.58e-22    101.0     47.55  50.1994472
    788      460   721      19   288  1.61e-22    100.0     48.38  50.1806379
    789      463   711      22   285  1.64e-22    100.0     45.65  50.1621758
    790      463   704      36   276  1.64e-22    100.0     50.97  50.1621758
    791      463   719     718   973  1.69e-22    105.0     46.99  50.1321435
    792      463   704      23   263  1.77e-22    100.0     50.19  50.0858925
    793      462   717      12   283  1.79e-22    100.0     46.62  50.0746564
    794      463   719      25   280  1.82e-22    100.0     47.55  50.0580355
    795      463   719      29   284  1.88e-22    100.0     47.55  50.0256003
    796      463   719      24   279  1.93e-22    100.0     47.55  49.9993520
    797      462   717      14   285  1.95e-22    100.0     46.62  49.9890427
    798      460   676      47   262  1.96e-22    100.0     49.33  49.9839276
    799      462   717      17   288  1.97e-22    100.0     46.59  49.9788385
    800      463   704      33   273  1.99e-22    101.0     51.16  49.9687374
    801      463   719      28   283  2.05e-22    100.0     47.55  49.9390323
    802      462   717      19   290  2.06e-22    100.0     46.62  49.9341661
    803      463   719      25   280  2.06e-22    100.0     47.55  49.9341661
    804      460   676      47   262  2.07e-22    100.0     49.33  49.9293234
    805      463   719      30   285  2.07e-22    100.0     47.55  49.9293234
    806      463   719      20   275  2.14e-22    100.0     47.55  49.8960662
    807      463   719      23   278  2.15e-22    100.0     47.55  49.8914042
    808      463   719      25   280  2.16e-22    100.0     47.55  49.8867638
    809      460   676      47   262  2.19e-22    100.0     49.33  49.8729705
    810      463   719      26   281  2.22e-22    100.0     47.55  49.8593648
    811      463   719      25   280  2.29e-22    100.0     47.55  49.8283202
    812      463   719      26   281  2.29e-22    100.0     47.55  49.8283202
    813      463   719      41   296  2.35e-22    100.0     47.55  49.8024567
    814      463   719      24   279  2.35e-22    100.0     47.55  49.8024567
    815      463   719      26   281  2.37e-22    100.0     47.55  49.7939821
    816      463   719      24   279  2.40e-22    100.0     47.55  49.7814033
    817      463   716      49   316  2.43e-22    100.0     46.86  49.7689808
    818      463   719      26   281  2.43e-22    100.0     47.55  49.7689808
    819      460   676      31   246  2.48e-22    100.0     49.33  49.7486135
    820      463   719      23   278  2.55e-22    100.0     47.55  49.7207787
    821      463   711      34   281  2.56e-22     99.8     48.28  49.7168648
    822      462   717      18   289  2.61e-22     99.8     46.26  49.6975218
    823      462   717      30   301  2.69e-22    100.0     46.62  49.6673309
    824      463   719      23   278  2.74e-22    100.0     47.55  49.6489141
    825      463   719      22   277  2.75e-22    100.0     47.55  49.6452711
    826      463   719      42   297  2.80e-22    100.0     47.92  49.6272526
    827      463   719      26   281  2.85e-22    100.0     47.55  49.6095531
    828      463   719      17   272  2.94e-22    100.0     47.55  49.5784625
    829      463   719      24   279  2.95e-22    100.0     47.55  49.5750669
    830      463   719      48   303  3.12e-22    100.0     47.55  49.5190390
    831      462   717      15   286  3.56e-22     99.4     46.26  49.3871115
    832      463   711      39   302  3.63e-22     99.8     45.65  49.3676394
    833      462   717      14   285  3.63e-22     99.0     46.26  49.3676394
    834      463   719      26   281  3.64e-22    100.0     47.55  49.3648884
    835      463   719      23   278  3.70e-22    100.0     47.92  49.3485392
    836      452   724      17   304  3.70e-22     99.4     50.34  49.3485392
    837      463   719      27   282  3.73e-22    100.0     47.55  49.3404638
    838      463   719      34   289  3.83e-22    100.0     47.55  49.3140072
    839      463   716      39   306  3.98e-22     99.8     46.86  49.2755902
    840      463   716      39   306  4.02e-22     99.8     46.86  49.2655901
    841      463   716      39   306  4.10e-22     99.8     46.86  49.2458851
    842      463   719      26   281  4.22e-22     99.8     47.55  49.2170369
    843      462   717      24   295  4.31e-22     99.8     46.59  49.1959341
    844      463   711      19   282  4.37e-22     99.4     45.65  49.1821090
    845      452   724      12   299  4.48e-22     99.4     50.34  49.1572490
    846      462   717      30   301  4.50e-22     99.8     46.59  49.1527946
    847      463   719     760  1015  4.54e-22    104.0     47.55  49.1439450
    848      463   719     760  1015  4.62e-22    104.0     47.55  49.1264773
    849      462   724      27   301  5.05e-22     99.4     50.18  49.0374838
    850      463   719      33   288  5.08e-22     99.8     47.55  49.0315608
    851      460   721       5   279  5.20e-22     99.0     48.24  49.0082134
    852      452   724      13   300  5.30e-22     99.4     50.34  48.9891652
    853      452   724      17   304  5.41e-22     99.4     50.34  48.9686230
    854      452   724      13   300  5.44e-22     99.4     49.66  48.9630930
    855      452   724       9   296  5.59e-22     99.0     50.00  48.9358928
    856      460   721       7   281  5.60e-22     99.0     48.24  48.9341054
    857      463   719      42   297  5.73e-22     99.8     47.55  48.9111565
    858      462   724      42   316  5.81e-22     99.4     50.18  48.8972915
    859      452   724      13   300  5.98e-22     99.4     49.66  48.8684515
    860      462   724      30   304  6.21e-22     99.0     50.18  48.8307111
    861      462   717      13   284  6.22e-22     98.6     46.26  48.8291021
    862      463   719      23   278  6.48e-22     99.4     47.55  48.7881515
    863      462   724      43   317  6.72e-22     99.4     49.47  48.7517839
    864      463   719      45   300  7.31e-22     99.8     47.55  48.6676288
    865      463   719      25   280  7.35e-22     99.4     47.55  48.6621717
    866      463   719      25   280  7.38e-22     99.4     47.55  48.6580984
    867      463   656      29   221  7.66e-22     99.0     49.50  48.6208601
    868      463   719      29   302  7.66e-22     98.6     47.69  48.6208601
    869      463   719      31   304  7.86e-22     98.6     47.69  48.5950854
    870      463   719      27   282  7.88e-22     99.0     47.17  48.5925441
    871      463   711      22   285  8.03e-22     98.6     45.65  48.5736875
    872      445   711      16   274  8.26e-22     98.6     49.46  48.5454475
    873      463   719      25   280  8.30e-22     99.0     47.55  48.5406165
    874      463   719      26   281  8.35e-22     99.0     47.55  48.5346105
    875      448   721       2   285  8.58e-22     97.8     47.47  48.5074381
    876      463   719      28   301  8.58e-22     98.2     47.69  48.5074381
    877      463   719      34   286  8.71e-22     99.0     47.51  48.4924003
    878      463   719      19   292  8.72e-22     98.2     47.69  48.4912528
    879      462   717      18   289  8.81e-22     98.6     46.26  48.4809846
    880      463   656      25   217  8.96e-22     98.6     49.50  48.4641018
    881      463   719      23   278  9.06e-22     99.0     47.55  48.4530029
    882      463   719      25   280  9.38e-22     99.0     46.79  48.4182923
    883      463   719      24   279  9.51e-22     99.0     46.79  48.4045282
    884      463   719      47   315  9.56e-22     98.6     49.27  48.3992843
    885      452   730       9   301  9.70e-22     98.2     47.25  48.3847462
    886      452   730      30   322  9.72e-22     98.6     47.25  48.3826864
    887      463   719      28   283  9.80e-22     98.6     46.79  48.3744897
    888      463   656      25   217  9.84e-22     98.6     49.50  48.3704163
    889      463   719      24   279  9.94e-22     98.6     47.17  48.3603050
    890      463   719      25   280  1.00e-21     99.0     47.17  48.3542870
    891      463   719      43   316  1.00e-21     98.6     47.69  48.3542870
    892      463   719     718   973  1.03e-21    102.0     46.62  48.3247282
    893      463   719      27   282  1.05e-21     98.6     47.17  48.3054968
    894      463   719      29   284  1.06e-21     99.0     47.17  48.2960180
    895      463   719      24   279  1.08e-21     98.6     47.17  48.2773259
    896      463   719      17   290  1.08e-21     97.8     47.69  48.2773259
    897      463   719      27   282  1.09e-21     98.6     47.55  48.2681093
    898      463   719      47   315  1.09e-21     98.6     49.27  48.2681093
    899      463   719      18   291  1.10e-21     97.8     47.69  48.2589768
    900      463   719      27   282  1.12e-21     97.8     47.17  48.2409583
    901      445   660      13   223  1.14e-21     97.8     52.65  48.2232587
    902      463   719      26   281  1.15e-21     98.6     47.17  48.2145250
    903      445   660      16   226  1.15e-21     98.2     52.65  48.2145250
    904      463   719      23   278  1.17e-21     98.6     47.17  48.1972832
    905      445   660      16   226  1.17e-21     97.8     52.65  48.1972832
    906      445   660       8   218  1.20e-21     97.8     52.65  48.1719654
    907      463   719      27   282  1.24e-21     98.6     47.17  48.1391756
    908      437   722       8   296  1.25e-21     97.8     46.67  48.1311434
    909      463   656      49   241  1.26e-21     98.6     49.50  48.1231752
    910      453   711      18   267  1.28e-21     97.8     48.67  48.1074269
    911      463   719      23   278  1.31e-21     98.2     47.17  48.0842598
    912      463   719      22   277  1.33e-21     98.2     47.17  48.0691080
    913      463   719      27   282  1.37e-21     98.6     47.17  48.0394762
    914      462   724      30   304  1.41e-21     97.8     49.10  48.0106972
    915      463   719      30   285  1.45e-21     98.6     47.17  47.9827234
    916      462   717      17   288  1.49e-21     97.4     46.26  47.9555108
    917      463   719      23   278  1.50e-21     98.2     47.17  47.9488218
    918      463   719      27   282  1.54e-21     98.2     47.17  47.9225045
    919      463   656      47   239  1.55e-21     98.2     49.50  47.9160320
    920      453   711      21   270  1.58e-21     97.4     50.57  47.8968621
    921      445   660      16   226  1.61e-21     97.4     52.94  47.8780528
    922      462   717      18   289  1.62e-21     97.8     45.88  47.8718608
    923      463   658      22   216  1.62e-21     97.4     49.02  47.8718608
    924      463   719      30   285  1.63e-21     98.2     47.17  47.8657069
    925      463   719      23   278  1.65e-21     98.2     47.17  47.8535117
    926      463   719      23   278  1.67e-21     98.2     47.17  47.8414633
    927      457   714      26   280  1.69e-21     97.4     47.58  47.8295584
    928      463   719      27   282  1.69e-21     98.2     47.17  47.8295584
    929      463   719      24   279  1.69e-21     98.2     47.17  47.8295584
    930      463   658      23   217  1.70e-21     97.4     49.02  47.8236587
    931      455   660      21   224  1.75e-21    100.0     51.90  47.7946712
    932      463   719      29   284  1.75e-21     98.2     47.17  47.7946712
    933      463   655      23   219  1.81e-21     97.4     55.45  47.7609601
    934      457   714      10   264  1.84e-21     96.7     47.58  47.7445214
    935      463   719      23   278  1.85e-21     97.8     47.55  47.7391013
    936      457   714      15   269  1.90e-21     97.1     47.58  47.7124331
    937      462   717      18   289  1.90e-21     97.1     46.26  47.7124331
    938      451   710     569   852  1.94e-21    101.0     45.08  47.6915990
    939      463   719      27   282  1.95e-21     97.8     47.17  47.6864576
    940      463   719      23   278  1.95e-21    101.0     47.55  47.6864576
    941      463   719     374   629  1.95e-21    101.0     47.55  47.6864576
    942      463   719      27   282  1.97e-21     97.8     47.17  47.6762534
    943      457   714      14   268  2.03e-21     96.7     47.58  47.6462512
    944      463   719      25   280  2.03e-21     97.8     47.17  47.6462512
    945      463   719      25   280  2.04e-21     97.8     47.17  47.6413371
    946      435   722      11   301  2.12e-21     97.4     46.37  47.6028709
    947      463   711      15   278  2.14e-21     97.1     44.93  47.5934811
    948      437   722      31   319  2.15e-21     97.8     46.67  47.5888191
    949      435   722      10   300  2.22e-21     97.4     46.37  47.5567798
    950      463   719      24   279  2.29e-21     97.8     47.17  47.5257351
    951      463   723      21   280  2.50e-21     97.4     46.30  47.4379962
    952      437   722      31   319  2.54e-21     97.4     46.67  47.4221229
    953      463   656      72   264  2.55e-21     98.2     49.50  47.4181936
    954      463   719      57   312  2.70e-21     98.2     47.17  47.3610352
    955      463   719      26   284  2.71e-21     97.4     47.01  47.3573383
    956      463   719      26   284  2.71e-21     97.4     47.01  47.3573383
    957      462   656      30   219  2.73e-21     96.7     55.28  47.3499853
    958      463   717      14   285  2.74e-21     96.7     47.67  47.3463290
    959      463   719      22   277  2.74e-21     97.4     47.17  47.3463290
    960      463   719      25   280  2.77e-21     97.4     47.17  47.3354396
    961      460   658      27   224  2.79e-21     96.7     50.72  47.3282454
    962      443   722      35   317  2.89e-21     97.4     46.60  47.2930305
    963      463   719     726   981  2.91e-21    101.0     46.79  47.2861339
    964      443   722      40   322  2.98e-21     97.4     46.60  47.2623637
    965      463   719     726   981  3.04e-21    101.0     46.79  47.2424294
    966      463   719      22   277  3.07e-21     97.4     47.17  47.2326094
    967      443   722      34   316  3.27e-21     97.4     46.60  47.1694970
    968      450   711       5   281  3.30e-21     96.3     47.16  47.1603645
    969      443   722      35   317  3.33e-21     97.4     46.60  47.1513146
    970      455   660      15   218  3.35e-21     97.1     51.43  47.1453266
    971      463   711      22   285  3.35e-21     96.7     45.29  47.1453266
    972      450   711      10   286  3.54e-21     96.7     47.16  47.0901602
    973      462   656      52   241  3.75e-21     97.1     55.28  47.0325311
    974      453   711      21   284  3.92e-21     96.7     48.35  46.9881953
    975      463   719      27   286  3.98e-21     97.1     46.64  46.9730051
    976      453   711      22   285  4.02e-21     97.1     48.35  46.9630051
    977      443   722      24   306  4.27e-21     96.7     46.60  46.9026731
    978      453   711      21   284  4.30e-21     96.7     48.35  46.8956719
    979      453   711      20   283  4.36e-21     96.7     48.35  46.8818149
    980      443   722      11   293  4.41e-21     96.3     46.60  46.8704123
    981      455   694      12   247  4.43e-21     96.7     48.40  46.8658874
    982      453   711      22   285  4.45e-21     97.1     48.35  46.8613829
    983      463   719      27   285  4.47e-21     97.1     47.78  46.8568985
    984      443   722       8   290  4.51e-21     96.3     46.60  46.8479898
    985      443   722      31   313  4.53e-21     96.7     46.60  46.8435650
    986      443   722      10   292  4.59e-21     96.3     46.60  46.8304069
    987      443   722       7   289  4.69e-21     96.3     46.60  46.8088544
    988      463   729      19   284  4.92e-21     95.9     46.32  46.7609784
    989      453   711      24   287  5.00e-21     96.7     48.35  46.7448490
    990      463   715      17   284  5.18e-21     95.9     45.00  46.7094819
    991      463   715      16   283  5.23e-21     95.9     45.00  46.6998757
    992      443   722      17   299  5.27e-21     96.3     46.60  46.6922566
    993      463   723     715   974  5.52e-21    100.0     46.30  46.6459091
    994      441   722       8   292  5.97e-21     95.9     46.62  46.5675400
    995      463   719      27   282  6.02e-21     96.7     46.79  46.5591997
    996      463   719      25   283  6.11e-21     96.7     47.78  46.5443602
    997      455   694      12   247  6.19e-21     96.3     48.40  46.5313519
    998      448   724      12   298  6.21e-21     95.9     47.67  46.5281261
    999      460   711      33   277  6.40e-21     95.9     46.69  46.4979890
    1000     463   719      31   289  6.45e-21     96.7     47.41  46.4902068
    1001     448   724      27   313  6.51e-21     95.9     47.67  46.4809475
    1002     441   722      21   305  6.85e-21     95.9     46.62  46.4300383
    1003     463   719      25   280  6.97e-21     96.3     46.79  46.4126717
    1004     453   711      21   284  7.09e-21     95.9     48.35  46.3956016
    1005     463   723      39   298  7.65e-21     96.3     45.56  46.3195813
    1006     463   719      30   285  7.88e-21     96.3     46.79  46.2899590
    1007     455   694       7   242  8.09e-21     95.9     50.00  46.2636582
    1008     463   723      21   280  8.64e-21     95.9     45.56  46.1978844
    1009     463   723      21   280  8.80e-21     95.9     45.56  46.1795352
    1010     451   734      17   318  8.88e-21     95.9     43.34  46.1704854
    1011     463   719      25   283  8.96e-21     95.9     47.78  46.1615167
    1012     463   711      22   285  9.00e-21     95.5     45.29  46.1570624
    1013     439   669      13   237  9.10e-21     95.5     50.00  46.1460125
    1014     463   711      22   285  9.15e-21     95.5     45.29  46.1405331
    1015     463   719      26   281  9.28e-21     95.9     46.79  46.1264254
    1016     463   723      42   301  9.34e-21     95.9     45.56  46.1199807
    1017     451   734      12   313  9.45e-21     95.9     43.34  46.1082722
    1018     463   711      31   284  9.68e-21     95.5     48.89  46.0842251
    1019     462   658      30   221  9.69e-21     95.1     52.45  46.0831925
    1020     455   694      12   247  9.78e-21     95.9     48.40  46.0739475
    1021     455   694      15   250  9.88e-21     95.9     48.40  46.0637744
    1022     439   669      13   237  9.99e-21     95.5     50.00  46.0527024
    1023     461   653      21   214  1.02e-20     95.5     53.54  46.0318992
    1024     463   719      24   279  1.02e-20     95.9     46.79  46.0318992
    1025     455   694      12   247  1.03e-20     95.9     48.40  46.0221431
    1026     463   656      44   236  1.03e-20     95.5     49.01  46.0221431
    1027     463   719      27   282  1.07e-20     95.9     46.79  45.9840432
    1028     463   655      24   220  1.11e-20     95.1     53.00  45.9473418
    1029     463   723      21   280  1.14e-20     95.5     45.56  45.9206736
    1030     439   669      13   237  1.14e-20     95.1     50.00  45.9206736
    1031     463   719      24   282  1.15e-20     95.5     47.78  45.9119399
    1032     463   719      25   280  1.15e-20     95.5     46.79  45.9119399
    1033     463   719      28   283  1.18e-20     95.5     46.79  45.8861874
    1034     455   718      20   287  1.26e-20     95.9     46.89  45.8205901
    1035     463   719      22   277  1.27e-20     95.5     46.79  45.8126850
    1036     463   719      25   283  1.32e-20     95.5     47.41  45.7740701
    1037     463   711      27   280  1.35e-20     95.1     48.89  45.7515973
    1038     462   658      52   243  1.36e-20     95.1     52.45  45.7442172
    1039     463   719      35   294  1.44e-20     95.9     47.60  45.6870587
    1040     453   711       6   255  1.60e-20     94.4     50.19  45.5816982
    1041     463   719      25   284  1.61e-20     95.5     47.60  45.5754677
    1042     439   669       3   227  1.82e-20     94.4     49.79  45.4528654
    1043     441   722      17   301  1.83e-20     94.4     46.30  45.4473859
    1044     463   719      44   299  1.85e-20     95.5     46.79  45.4365162
    1045     441   722      15   299  1.86e-20     94.4     46.30  45.4311254
    1046     441   720       5   287  1.99e-20     94.0     47.25  45.3635672
    1047     461   711      56   311  2.13e-20     95.1     48.90  45.2955799
    1048     463   726      45   320  2.37e-20     94.7     44.41  45.1888119
    1049     442   663       1   220  2.42e-20     94.7     51.10  45.1679343
    1050     443   722      23   305  2.45e-20     94.4     46.28  45.1556138
    1051     460   711      14   258  2.52e-20     93.6     47.67  45.1274430
    1052     157   227       3    73  2.56e-20     87.0     71.83  45.1116946
    1053     442   663       1   220  2.64e-20     94.4     50.66  45.0809229
    1054     463   711     361   624  2.67e-20     97.8     44.89  45.0696234
    1055     442   663       1   220  2.81e-20     94.4     50.66  45.0185174
    1056     455   663      15   221  3.35e-20     94.4     51.40  44.8427415
    1057     157   227      10    80  3.70e-20     87.0     71.83  44.7433690
    1058     157   228       5    79  3.72e-20     87.0     72.00  44.7379782
    1059     463   719      27   282  3.77e-20     94.0     47.92  44.7246269
    1060     463   726      20   295  3.88e-20     93.6     44.41  44.6958667
    1061     450   711       6   282  4.03e-20     94.0     45.74  44.6579355
    1062     459   669       5   209  4.60e-20     92.8     51.42  44.5256456
    1063     463   714      18   268  4.71e-20     92.4     46.67  44.5020140
    1064     451   734      15   322  5.12e-20     93.6     42.55  44.4185474
    1065     451   734       6   313  5.36e-20     93.2     42.55  44.3727379
    1066     450   711      13   289  5.37e-20     94.0     45.74  44.3708740
    1067     442   663       1   220  5.56e-20     93.6     50.66  44.3361038
    1068     463   719      29   284  5.74e-20     93.6     47.92  44.3042426
    1069     442   663       1   220  5.82e-20     93.6     50.66  44.2904016
    1070     157   227       4    77  5.90e-20     86.3     71.62  44.2767495
    1071     157   227       5    78  6.06e-20     86.3     71.62  44.2499921
    1072     455   663      33   239  6.35e-20     93.6     51.40  44.2032470
    1073     463   719      27   282  6.35e-20     93.6     47.92  44.2032470
    1074     157   227       7    80  6.42e-20     86.3     71.62  44.1922837
    1075     451   734      41   348  7.28e-20     93.6     42.55  44.0665710
    1076     451   734      42   349  7.31e-20     93.6     42.55  44.0624586
    1077     157   227       8    81  8.02e-20     85.9     71.62  43.9697634
    1078     457   696       9   278  9.52e-20     92.4     49.28  43.7983070
    1079     466   714      27   282  9.82e-20     92.0     50.76  43.7672807
    1080     451   734      28   335  9.86e-20     93.2     42.55  43.7632157
    1081     151   227       5    86  1.21e-19     85.5     68.29  43.5584964
    1082     452   722      10   292  1.23e-19     92.0     49.49  43.5421026
    1083     455   663       5   211  1.32e-19     92.4     51.40  43.4714850
    1084     457   696       6   275  1.40e-19     91.7     48.74  43.4126445
    1085     457   696      44   313  1.45e-19     92.4     48.74  43.3775532
    1086     457   696       6   275  1.46e-19     91.7     48.74  43.3706803
    1087     457   696       8   277  1.51e-19     91.7     48.74  43.3370071
    1088     457   696      45   314  1.52e-19     92.4     48.74  43.3304064
    1089     457   696       8   277  1.55e-19     91.7     48.74  43.3108618
    1090     457   696       5   274  1.56e-19     91.7     48.74  43.3044309
    1091     441   669      13   235  1.58e-19     92.0     50.64  43.2916919
    1092     453   669       9   219  1.60e-19     91.7     51.58  43.2791131
    1093     457   696      11   280  1.61e-19     91.7     48.74  43.2728826
    1094     457   696      31   300  1.67e-19     92.0     48.74  43.2362931
    1095     457   696       7   276  1.69e-19     91.7     48.74  43.2243882
    1096     457   696       8   277  1.71e-19     91.7     48.74  43.2126234
    1097     453   669       9   219  1.72e-19     91.7     51.33  43.2067925
    1098     457   696       8   277  1.79e-19     91.7     48.74  43.1669011
    1099     441   704       8   274  1.92e-19     91.3     46.76  43.0967916
    1100     441   704       5   271  2.06e-19     91.3     46.42  43.0264108
    1101     457   696       9   278  2.07e-19     91.7     48.74  43.0215682
    1102     455   699      10   259  2.22e-19     92.0     48.84  42.9516096
    1103     455   699      18   267  2.40e-19     92.0     48.84  42.8736480
    1104     455   694      14   249  2.42e-19     91.7     49.80  42.8653492
    1105     157   228      21    95  2.48e-19     85.5     72.00  42.8408582
    1106     455   694      15   250  2.85e-19     91.7     49.80  42.7017978
    1107     455   663       7   213  2.95e-19     91.3     50.69  42.6673116
    1108     463   740      25   316  3.10e-19     91.3     47.18  42.6177147
    1109     462   723      16   281  3.29e-19     90.5     48.23  42.5582292
    1110     451   710       5   293  3.47e-19     90.9     44.97  42.5049622
    1111     157   227       8    81  3.92e-19     84.0     70.27  42.3830251
    1112     442   663       1   220  4.26e-19     90.9     49.78  42.2998476
    1113     157   227       8    81  4.45e-19     84.0     70.27  42.2562127
    1114     457   722      17   297  4.88e-19     90.1     46.62  42.1639715
    1115     461   660      15   213  6.33e-19     90.5     54.15  41.9038165
    1116     463   712      25   289  6.48e-19     89.7     47.79  41.8803963
    1117     463   712      25   289  6.61e-19     89.7     47.79  41.8605331
    1118     451   734      28   335  6.78e-19     90.5     43.20  41.8351397
    1119     461   649      21   240  7.19e-19     90.1     49.78  41.7764256
    1120     463   712      50   314  7.24e-19     90.5     47.79  41.7694956
    1121     461   660      15   213  8.48e-19     90.1     54.15  41.6114063
    1122     453   669       9   219  8.76e-19     89.4     51.33  41.5789209
    1123     463   664      30   239  9.20e-19     89.4     49.07  41.5299133
    1124     459   711      17   263  1.14e-18     88.6     47.67  41.3155034
    1125     456   696      38   308  1.17e-18     89.7     48.93  41.2895279
    1126     456   707      23   273  1.19e-18     89.7     47.10  41.2725784
    1127     450   659       4   242  1.26e-18     89.0     48.16  41.2154200
    1128     459   711      11   257  1.55e-18     87.8     47.67  41.0082767
    1129     453   669       8   218  1.59e-18     89.7     51.58  40.9827977
    1130     459   711      15   261  1.65e-18     88.2     47.67  40.9457564
    1131     459   711      11   257  1.74e-18     87.8     47.67  40.8926466
    1132     453   669      24   234  1.79e-18     89.7     51.58  40.8643161
    1133     449   714      13   285  2.04e-18     88.6     46.08  40.7335819
    1134     463   668      46   257  2.06e-18     89.4     49.77  40.7238257
    1135     450   659      33   271  2.15e-18     89.0     48.16  40.6810638
    1136     463   712      29   277  2.37e-18     88.2     49.04  40.5836417
    1137     449   714      23   295  2.40e-18     88.6     46.58  40.5710629
    1138     463   704      27   267  2.45e-18     88.6     49.80  40.5504436
    1139     463   704      26   266  2.55e-18     87.8     49.80  40.5104383
    1140     450   732      31   294  2.68e-18     88.2     45.39  40.4607149
    1141     463   668      38   249  2.69e-18     88.6     49.31  40.4569905
    1142     457   722      17   297  2.91e-18     88.2     46.28  40.3783786
    1143     454   707       9   260  3.35e-18     87.4     46.01  40.2375713
    1144     450   732       5   268  3.36e-18     87.0     45.05  40.2345907
    1145     450   732       6   269  3.49e-18     87.4     44.93  40.1966299
    1146     450   732       6   269  3.62e-18     87.0     45.05  40.1600576
    1147     450   659      42   280  3.77e-18     88.6     48.16  40.1194567
    1148     450   659      43   281  3.82e-18     88.6     48.16  40.1062813
    1149     462   653      29   233  4.05e-18     88.2     50.95  40.0478148
    1150     450   659      40   278  4.19e-18     88.2     48.16  40.0138309
    1151     447   732       5   271  4.19e-18     86.7     45.27  40.0138309
    1152     463   712      43   291  4.23e-18     88.2     49.04  40.0043297
    1153     463   704      13   247  4.41e-18     87.8     49.20  39.9626570
    1154     463   719     880  1153  4.60e-18     91.3     46.62  39.9204754
    1155     463   704      13   247  5.20e-18     87.8     49.20  39.7978730
    1156     461   649      34   253  5.24e-18     87.8     49.78  39.7902102
    1157     450   732       5   268  5.33e-18     86.7     45.39  39.7731804
    1158     463   704      13   247  5.47e-18     87.8     49.20  39.7472531
    1159     463   704      13   247  5.55e-18     87.8     49.20  39.7327337
    1160     463   704      16   250  5.68e-18     87.8     49.20  39.7095804
    1161     463   712      30   278  5.70e-18     87.0     48.66  39.7060655
    1162     457   656      20   219  5.82e-18     87.8     49.03  39.6852314
    1163     463   712      28   276  5.84e-18     87.0     48.66  39.6818009
    1164     450   659      66   304  5.90e-18     88.2     48.16  39.6715793
    1165     450   732      31   294  5.93e-18     87.0     45.05  39.6665075
    1166     463   704      18   252  6.36e-18     87.8     49.20  39.5965033
    1167     457   656      21   220  6.40e-18     87.8     49.03  39.5902337
    1168     450   732      34   297  6.49e-18     87.0     45.39  39.5762691
    1169     463   704      16   256  6.55e-18     86.3     49.41  39.5670666
    1170     463   712      30   278  6.68e-18     86.7     48.66  39.5474137
    1171     457   656       7   206  6.77e-18     87.4     49.03  39.5340306
    1172     450   732      10   273  7.26e-18     86.3     44.93  39.4641518
    1173     461   688     137   381  7.26e-18     90.5     48.81  39.4641518
    1174     450   732       6   269  7.42e-18     86.3     45.39  39.4423526
    1175     447   732       5   271  7.64e-18     86.3     44.93  39.4131341
    1176     463   704      13   247  7.65e-18     87.4     49.20  39.4118260
    1177     463   704      21   261  7.83e-18     86.7     49.41  39.3885692
    1178     463   704      17   257  8.12e-18     85.9     49.41  39.3522015
    1179     442   732       2   273  8.24e-18     85.9     44.52  39.3375313
    1180     461   720      23   276  8.30e-18     86.3     46.35  39.3302762
    1181     442   732       1   272  8.33e-18     85.9     44.52  39.3266682
    1182     463   704     160   394  8.42e-18     89.0     49.20  39.3159218
    1183     450   732      29   292  8.43e-18     86.7     44.93  39.3147349
    1184     461   720      31   284  8.46e-18     86.3     46.35  39.3111825
    1185     442   732       1   272  8.57e-18     85.9     44.52  39.2982639
    1186     450   732       6   269  8.62e-18     86.3     45.39  39.2924466
    1187     442   732       1   272  8.98e-18     85.9     44.52  39.2515318
    1188     449   714      20   292  9.33e-18     87.0     46.08  39.2132967
    1189     450   732      10   273  9.35e-18     86.3     44.93  39.2111553
    1190     461   720      23   276  9.43e-18     86.3     46.35  39.2026356
    1191     454   707     335   586  9.51e-18     89.4     46.01  39.1941878
    1192     463   704      17   257  9.51e-18     85.9     49.41  39.1941878
    1193     461   720      23   276  9.58e-18     86.3     46.35  39.1868541
    1194     463   704      21   261  9.63e-18     85.9     49.41  39.1816484
    1195     461   720      33   286  9.68e-18     86.3     46.35  39.1764698
    1196     450   732      10   273  9.89e-18     85.9     44.93  39.1550075
    1197     457   656      19   218  1.01e-17     87.0     48.54  39.1339963
    1198     450   732      10   273  1.02e-17     85.9     44.93  39.1241440
    1199     463   655      50   246  1.02e-17     87.0     50.99  39.1241440
    1200     461   720      33   286  1.06e-17     86.3     46.35  39.0856777
    1201     461   720      21   274  1.07e-17     85.9     46.35  39.0762879
    1202     463   717      92   369  1.08e-17     88.2     45.61  39.0669855
    1203     463   704     157   391  1.10e-17     88.2     49.20  39.0486364
    1204     450   732      10   273  1.12e-17     85.9     44.59  39.0306179
    1205     450   732       8   271  1.12e-17     85.9     44.93  39.0306179
    1206     449   714      19   291  1.13e-17     86.7     46.08  39.0217289
    1207     454   707     335   586  1.15e-17     89.4     46.01  39.0041846
    1208     450   711      31   281  1.17e-17     85.9     45.82  38.9869428
    1209     461   720      20   273  1.17e-17     85.9     46.35  38.9869428
    1210     450   732       5   268  1.22e-17     85.5     45.05  38.9450957
    1211     461   720      11   264  1.28e-17     85.5     46.35  38.8970865
    1212     450   732      22   285  1.31e-17     85.9     45.05  38.8739194
    1213     450   732       2   265  1.32e-17     85.1     45.39  38.8663148
    1214     461   720      17   270  1.36e-17     85.5     46.35  38.8364619
    1215     450   732       6   269  1.41e-17     85.5     45.05  38.8003569
    1216     462   723     202   467  1.42e-17     88.2     48.93  38.7932897
    1217     450   732       3   266  1.45e-17     85.5     45.05  38.7723830
    1218     450   732       7   270  1.45e-17     85.5     45.05  38.7723830
    1219     450   732       2   265  1.46e-17     85.1     44.93  38.7655101
    1220     461   649      20   236  1.47e-17     86.3     49.55  38.7586842
    1221     461   720      20   273  1.47e-17     85.5     46.35  38.7586842
    1222     450   732       9   272  1.50e-17     85.5     45.05  38.7384815
    1223     463   723     559   822  1.50e-17     89.7     48.91  38.7384815
    1224     450   732       3   266  1.53e-17     85.5     45.05  38.7186788
    1225     463   723     559   822  1.54e-17     89.4     48.91  38.7121642
    1226     450   732       4   267  1.55e-17     85.1     45.39  38.7056916
    1227     461   649      18   234  1.58e-17     86.3     49.55  38.6865217
    1228     441   711       1   255  1.59e-17     85.1     45.77  38.6802126
    1229     461   720      23   276  1.59e-17     85.5     46.35  38.6802126
    1230     450   732      10   273  1.61e-17     85.5     45.05  38.6677124
    1231     461   720      23   276  1.61e-17     85.5     46.35  38.6677124
    1232     450   732       6   269  1.63e-17     85.5     45.05  38.6553666
    1233     450   732      10   273  1.64e-17     85.5     45.05  38.6492503
    1234     450   732       4   267  1.65e-17     85.1     45.39  38.6431713
    1235     461   720      13   266  1.66e-17     85.1     46.35  38.6371290
    1236     450   732      10   273  1.70e-17     85.5     45.05  38.6133183
    1237     441   711       1   255  1.75e-17     84.7     45.77  38.5843308
    1238     463   723     559   822  1.77e-17     89.4     50.00  38.5729670
    1239     461   720      33   286  1.78e-17     85.5     45.99  38.5673332
    1240     463   723     552   815  1.78e-17     89.4     50.00  38.5673332
    1241     463   723     553   816  1.81e-17     89.4     50.00  38.5506197
    1242     450   732       7   270  1.81e-17     85.1     45.05  38.5506197
    1243     450   732      10   273  1.82e-17     85.1     44.71  38.5451101
    1244     461   720      13   266  1.84e-17     85.1     46.35  38.5341810
    1245     450   732       5   268  1.84e-17     85.1     45.05  38.5341810
    1246     450   732      10   273  1.85e-17     85.1     44.71  38.5287609
    1247     461   720      21   274  1.85e-17     85.1     45.99  38.5287609
    1248     450   732      10   273  1.88e-17     85.1     45.05  38.5126748
    1249     463   723     556   819  1.89e-17     89.4     50.00  38.5073698
    1250     450   732       5   268  1.91e-17     85.1     45.05  38.4968433
    1251     450   732       8   271  1.92e-17     85.1     45.05  38.4916214
    1252     461   720      23   276  1.93e-17     85.5     46.35  38.4864266
    1253     461   720     375   628  1.94e-17     88.6     46.35  38.4812586
    1254     450   732       3   266  1.96e-17     84.7     44.93  38.4710021
    1255     461   720      13   266  1.97e-17     85.1     46.35  38.4659130
    1256     461   720      22   275  1.97e-17     85.1     45.99  38.4659130
    1257     469   711      25   289  1.98e-17     85.5     44.21  38.4608497
    1258     450   732      24   287  1.99e-17     85.5     45.05  38.4558119
    1259     461   720     376   629  1.99e-17     88.6     46.35  38.4558119
    1260     450   732      10   273  2.03e-17     85.1     45.05  38.4359108
    1261     463   675      43   264  2.05e-17     85.9     46.70  38.4261068
    1262     450   732       5   268  2.18e-17     85.1     44.59  38.3646217
    1263     450   732       5   268  2.25e-17     84.7     45.05  38.3330164
    1264     463   723    1885  2148  2.39e-17     89.4     50.00  38.2726532
    1265     463   723    1885  2148  2.43e-17     89.4     50.00  38.2560553
    1266     450   732       5   268  2.50e-17     84.3     45.05  38.2276558
    1267     442   732       8   279  2.51e-17     84.7     44.19  38.2236638
    1268     455   715      18   308  2.51e-17     85.9     47.81  38.2236638
    1269     441   711       1   255  2.55e-17     84.7     45.91  38.2078532
    1270     450   732       4   267  2.55e-17     84.3     45.05  38.2078532
    1271     450   711       9   259  2.58e-17     84.3     45.59  38.1961572
    1272     450   732       5   268  2.64e-17     84.3     45.05  38.1731677
    1273     441   711       1   255  2.72e-17     84.3     45.91  38.1433147
    1274     463   723    1920  2183  2.73e-17     89.0     50.00  38.1396450
    1275     441   711       1   255  2.77e-17     84.3     45.91  38.1250993
    1276     456   649      30   237  2.92e-17     85.5     52.31  38.0723630
    1277     455   656       7   206  2.93e-17     85.5     50.97  38.0689442
    1278     461   720      16   269  3.11e-17     84.3     47.27  38.0093239
    1279     461   720      17   270  3.13e-17     84.3     47.27  38.0029136
    1280     450   711      31   281  3.18e-17     84.7     45.59  37.9870654
    1281     456   649      38   245  3.21e-17     85.5     52.31  37.9776756
    1282     456   649      39   246  3.24e-17     85.5     52.31  37.9683733
    1283     450   711      31   281  3.30e-17     84.7     45.59  37.9500241
    1284     450   711      31   281  3.40e-17     84.7     45.59  37.9201711
    1285     450   732       2   265  3.42e-17     84.0     45.05  37.9143060
    1286     463   723    1343  1606  3.53e-17     88.6     49.64  37.8826487
    1287     450   711      31   281  3.53e-17     84.7     45.59  37.8826487
    1288     456   649      15   222  3.57e-17     85.1     52.31  37.8713810
    1289     456   649       7   214  3.61e-17     84.7     52.31  37.8602388
    1290     463   675      43   264  3.70e-17     85.1     46.26  37.8356138
    1291     456   649       9   216  3.75e-17     84.7     52.31  37.8221907
    1292     463   740      73   375  4.15e-17     85.9     46.93  37.7208382
    1293     443   711      30   282  4.16e-17     84.7     47.08  37.7184315
    1294     450   711      31   281  4.54e-17     84.3     45.82  37.6310196
    1295     450   711       2   252  4.59e-17     83.6     45.09  37.6200666
    1296     459   703      38   289  4.64e-17     84.3     47.91  37.6092322
    1297     463   723     553   816  4.68e-17     87.8     48.55  37.6006485
    1298     450   711       1   251  4.87e-17     83.6     45.09  37.5608526
    1299     450   711       8   258  5.22e-17     83.6     45.09  37.4914492
    1300     463   675      43   264  5.36e-17     84.7     46.70  37.4649826
    1301     446   658       2   213  5.36e-17     83.6     49.32  37.4649826
    1302     450   711       8   258  5.47e-17     83.6     45.09  37.4446680
    1303     456   649      38   245  5.80e-17     84.7     52.31  37.3860887
    1304     456   649      38   245  5.85e-17     84.7     52.31  37.3775049
    1305     461   660      19   217  5.89e-17     84.7     52.20  37.3706906
    1306     461   660      16   214  6.24e-17     84.7     52.20  37.3129664
    1307     461   660      16   214  6.44e-17     84.7     52.20  37.2814180
    1308     461   660      19   217  6.45e-17     84.7     52.20  37.2798665
    1309     463   655      49   248  6.80e-17     84.3     50.98  37.2270240
    1310     463   723    1885  2148  6.83e-17     87.8     49.64  37.2226219
    1311     456   649      38   245  6.94e-17     84.3     52.31  37.2066448
    1312     461   660      16   214  7.22e-17     84.3     52.20  37.1670916
    1313     450   711       2   252  7.22e-17     83.2     45.09  37.1670916
    1314     456   649       9   216  7.70e-17     84.0     50.46  37.1027263
    1315     446   658       2   213  8.41e-17     83.2     47.71  37.0145251
    1316     461   660      16   214  8.55e-17     84.3     52.20  36.9980153
    1317     446   732       5   272  8.57e-17     83.2     44.15  36.9956788
    1318     457   722      11   305  9.15e-17     83.6     45.16  36.9301927
    1319     449   658      21   229  9.46e-17     86.3     49.54  36.8968742
    1320     461   660      16   214  9.85e-17     84.0     52.20  36.8564751
    1321     449   658       8   216  1.00e-16     85.9     49.54  36.8413615
    1322     449   658      27   235  1.05e-16     85.9     49.54  36.7925713
    1323     438   693      74   344  1.09e-16     84.3     47.86  36.7551838
    1324     463   723     553   816  1.11e-16     86.7     49.64  36.7370015
    1325     449   658      15   223  1.12e-16     85.9     49.54  36.7280328
    1326     442   667       3   229  1.14e-16     83.6     50.00  36.7103332
    1327     461   660      14   212  1.18e-16     83.6     52.20  36.6758470
    1328     442   667       5   231  1.25e-16     82.8     50.00  36.6182179
    1329     442   667       3   229  1.28e-16     82.8     50.00  36.5945014
    1330     442   667       3   229  1.37e-16     82.4     50.00  36.5265507
    1331     463   723    1885  2148  1.41e-16     86.7     49.64  36.4977718
    1332     463   707      28   281  1.42e-16     83.2     46.42  36.4907046
    1333     461   658      35   231  1.46e-16     84.0     52.22  36.4629251
    1334     463   723    1885  2148  1.46e-16     86.7     49.64  36.4629251
    1335     457   649      10   216  1.47e-16     82.8     52.09  36.4560991
    1336     462   723      28   293  1.57e-16     82.8     49.28  36.3902859
    1337     449   658       4   212  1.58e-16     82.4     49.54  36.3839366
    1338     461   658      35   231  1.60e-16     83.6     52.22  36.3713579
    1339     456   677      33   260  1.62e-16     83.2     51.05  36.3589353
    1340     449   658      22   230  1.62e-16     85.5     49.54  36.3589353
    1341     456   677      31   258  1.69e-16     83.2     51.05  36.3166330
    1342     456   677      25   252  1.73e-16     83.2     51.05  36.2932401
    1343     442   667      20   246  1.84e-16     82.4     50.00  36.2315959
    1344     449   658      10   218  1.85e-16     84.7     50.00  36.2261758
    1345     456   677      31   258  1.90e-16     82.8     51.05  36.1995076
    1346     446   732       5   272  1.94e-16     82.0     43.81  36.1786735
    1347     463   704      34   269  1.96e-16     83.6     49.00  36.1684170
    1348     456   677      28   255  2.10e-16     82.8     51.05  36.0994241
    1349     463   704      18   253  2.20e-16     83.2     49.42  36.0529041
    1350     463   657      28   222  2.21e-16     83.2     49.25  36.0483690
    1351     449   658       5   213  2.23e-16     82.0     47.91  36.0393599
    1352     457   649      32   238  2.26e-16     82.8     52.09  36.0259967
    1353     459   667      16   222  2.32e-16     82.8     50.69  35.9997943
    1354     463   704      16   251  2.42e-16     82.8     49.00  35.9575939
    1355     455   658       6   208  2.44e-16     84.7     50.94  35.9493634
    1356     455   658       6   208  2.46e-16     84.7     50.94  35.9412001
    1357     463   704      17   252  2.50e-16     82.8     49.42  35.9250708
    1358     463   704      13   248  2.66e-16     82.8     49.42  35.8630354
    1359     463   657      28   222  2.83e-16     82.8     48.76  35.8010848
    1360     456   677      26   253  2.88e-16     82.0     51.05  35.7835712
    1361     459   667      16   222  3.01e-16     82.4     50.69  35.7394214
    1362     463   704      16   251  3.21e-16     82.4     49.00  35.6750906
    1363     456   677      25   252  3.24e-16     82.0     51.05  35.6657882
    1364     461   658      35   231  3.32e-16     82.8     52.22  35.6413967
    1365     463   669      36   243  3.40e-16     82.0     50.68  35.6175861
    1366     450   670      47   269  3.72e-16     82.8     49.79  35.5276378
    1367     449   658      27   235  4.05e-16     83.6     50.46  35.4426446
    1368     455   660      11   218  4.10e-16     83.6     52.13  35.4303745
    1369     443   711       5   257  4.28e-16     80.9     47.08  35.3874085
    1370     455   658      14   216  4.29e-16     82.4     50.94  35.3850748
    1371     455   660      11   218  4.81e-16     83.2     52.13  35.2706644
    1372     455   658       4   206  4.90e-16     83.6     50.94  35.2521263
    1373     459   656      36   236  4.91e-16     81.6     51.21  35.2500875
    1374     556   735     770   944  5.20e-16     84.3     48.91  35.1927029
    1376     455   658       4   206  5.21e-16     83.2     50.94  35.1907816
    1377     458   743       7   292  5.23e-16     81.6     47.99  35.1869502
    1378     455   658      14   216  5.35e-16     83.2     50.94  35.1642649
    1379     556   735     766   940  5.60e-16     84.3     48.91  35.1185949
    1381     556   735     768   942  5.71e-16     84.3     48.91  35.0991425
    1383     459   656      13   213  5.79e-16     80.9     51.69  35.0852292
    1384     455   658      14   216  5.90e-16     83.2     50.94  35.0664091
    1385     450   732      10   273  5.95e-16     80.9     44.03  35.0579703
    1386     463   701      14   275  7.23e-16     80.9     42.32  34.8631225
    1387     459   656      36   236  7.26e-16     80.9     51.21  34.8589817
    1388     463   701      14   275  7.29e-16     80.9     42.32  34.8548579
    1389     460   715    1246  1522  7.76e-16     84.3     47.28  34.7923792
    1390     455   658       8   209  7.79e-16     80.1     50.95  34.7885206
    1391     455   658       4   205  7.86e-16     80.1     50.95  34.7795749
    1392     455   658      13   214  8.06e-16     80.1     50.95  34.7544479
    1393     463   704     159   394  8.49e-16     82.4     49.19  34.7024725
    1394     462   669       9   216  8.65e-16     80.9     49.53  34.6838022
    1395     463   704     156   391  8.88e-16     82.4     49.19  34.6575599
    1396     460   715    1226  1502  8.95e-16     84.0     47.28  34.6497080
    1397     455   658      14   215  8.95e-16     80.1     50.95  34.6497080
    1398     460   715    1278  1554  9.15e-16     84.0     47.28  34.6276076
    1399     459   656      36   236  9.47e-16     80.5     51.21  34.5932326
    1400     463   707      28   281  9.87e-16     80.9     46.42  34.5518616
    1401     459   656      38   238  1.01e-15     80.5     51.21  34.5288261
    1402     459   656      36   236  1.05e-15     80.5     51.21  34.4899862
    1403     459   656      36   236  1.06e-15     80.5     51.21  34.4805075
    1404     448   657       5   207  1.08e-15     79.7     50.23  34.4618154
    1405     463   663      14   223  1.09e-15     79.7     49.55  34.4525987
    1406     448   657       7   209  1.10e-15     80.1     50.23  34.4434662
    1407     448   657       7   209  1.11e-15     80.1     50.23  34.4344164
    1408     463   704     170   405  1.11e-15     82.0     49.19  34.4344164
    1409     461   657      16   212  1.13e-15     80.1     50.97  34.4165588
    1410     459   656      36   236  1.14e-15     80.5     51.21  34.4077481
    1411     461   657      16   212  1.19e-15     80.1     50.97  34.3648231
    1412     448   657       7   209  1.20e-15     80.1     50.23  34.3564548
    1413     448   657       6   208  1.21e-15     79.7     50.23  34.3481560
    1414     448   657       6   208  1.22e-15     79.7     50.23  34.3399255
    1415     448   657       7   209  1.22e-15     79.7     50.23  34.3399255
    1416     448   657       7   209  1.24e-15     80.1     50.23  34.3236650
    1417     463   704     600   840  1.26e-15     83.2     49.41  34.3076647
    1418     448   656       6   207  1.31e-15     80.5     50.46  34.2687493
    1419     448   657       7   209  1.35e-15     80.1     50.68  34.2386718
    1420     448   656       7   208  1.39e-15     80.5     50.46  34.2094726
    1421     463   704     151   386  1.43e-15     81.6     49.19  34.1811020
    1422     447   657      11   220  1.47e-15     80.1     46.19  34.1535140
    1423     448   657       6   208  1.48e-15     79.3     50.23  34.1467343
    1424     448   657       7   209  1.48e-15     79.3     50.23  34.1467343
    1425     455   657       9   211  1.52e-15     81.6     52.40  34.1200661
    1426     448   657       6   208  1.53e-15     79.3     50.23  34.1135087
    1427     463   663      14   223  1.56e-15     81.6     48.86  34.0940906
    1428     500   710     233   451  1.56e-15     82.4     46.55  34.0940906
    1429     448   657       6   208  1.57e-15     79.3     50.23  34.0877008
    1430     448   657       6   208  1.58e-15     79.3     50.23  34.0813515
    1431     448   657       6   208  1.62e-15     79.3     50.23  34.0563502
    1432     463   704     148   383  1.64e-15     81.6     49.19  34.0440802
    1433     500   710     233   451  1.68e-15     82.0     46.55  34.0199826
    1434     463   657      27   221  1.69e-15     80.5     47.26  34.0140479
    1435     448   657       7   209  1.70e-15     79.0     50.23  34.0081481
    1436     463   704     148   383  1.71e-15     81.3     49.19  34.0022830
    1437     463   663      14   223  1.72e-15     81.3     48.86  33.9964521
    1438     448   657       6   208  1.74e-15     79.0     50.23  33.9848913
    1439     448   657       7   209  1.78e-15     79.3     50.23  33.9621630
    1440     448   657       6   208  1.90e-15     79.3     50.23  33.8969225
    1441     434   714      17   306  1.96e-15     79.7     41.91  33.8658319
    1442     455   656      13   217  2.03e-15     79.0     50.24  33.8307406
    1443     448   657       6   208  2.13e-15     79.0     50.23  33.7826544
    1444     447   657      11   220  2.49e-15     79.3     45.74  33.6264937
    1445     463   663      14   223  2.66e-15     78.6     49.32  33.5604503
    1446     463   656      23   217  2.70e-15     80.1     49.00  33.5455246
    1447     463   663      14   223  2.81e-15     78.6     49.32  33.5055919
    1448     501   672     164   330  2.91e-15     82.4     56.00  33.4706233
    1449     447   657      13   222  3.07e-15     78.6     45.74  33.4170988
    1450     423   712       3   271  3.27e-15     78.6     48.15  33.3539864
    1451     447   657      11   220  3.60e-15     79.0     45.74  33.2578425
    1452     447   657       8   217  3.82e-15     78.2     45.74  33.1985260
    1453     463   657      28   222  3.86e-15     79.3     46.77  33.1881092
    1454     447   657       9   218  3.96e-15     78.2     45.74  33.1625324
    1455     463   704      48   283  4.27e-15     79.7     44.76  33.0871626
    1456     447   657      11   220  4.41e-15     78.6     45.74  33.0549017
    1457     448   657       7   209  4.50e-15     77.8     49.77  33.0346990
    1458     447   657      11   220  4.53e-15     78.6     45.74  33.0280545
    1459     461   650      11   224  4.80e-15     78.2     46.58  32.9701605
    1460     463   701      14   275  5.02e-15     78.2     41.20  32.9253465
    1461     447   657      28   237  5.03e-15     78.2     45.74  32.9233564
    1462     459   656      36   236  5.19e-15     78.6     51.21  32.8920427
    1463     448   657       6   208  5.25e-15     77.8     49.77  32.8805483
    1464     463   704      46   281  5.46e-15     79.3     44.76  32.8413276
    1465     463   707      28   281  6.33e-15     78.6     46.04  32.6934762
    1466     455   674      29   248  9.63e-15     77.4     50.66  32.2738932
    1467     456   695       9   258  1.42e-14     77.0     46.88  31.8855344
    1468     463   656     193   389  1.52e-14     79.0     48.28  31.8174810
    1469     456   695      12   261  1.66e-14     77.0     46.88  31.7293737
    1470     455   672      63   286  1.71e-14     78.2     44.69  31.6996979
    1471     456   695      15   264  1.89e-14     76.6     46.88  31.5996145
    1472     463   656     193   389  2.05e-14     78.6     48.28  31.5183515
    1473     463   656     193   389  2.07e-14     78.6     48.28  31.5086427
    1474     463   656     193   389  2.09e-14     78.6     48.28  31.4990272
    1475     463   656     193   389  2.09e-14     78.6     48.28  31.4990272
    1476     455   674      29   248  2.09e-14     76.6     50.22  31.4990272
    1477     455   674      29   248  2.11e-14     76.6     50.22  31.4895034
    1478     463   656     193   389  2.22e-14     78.6     48.28  31.4386841
    1479     463   656     193   389  2.23e-14     78.6     48.28  31.4341897
    1480     157   227       8    81  2.53e-14     70.5     64.86  31.3079720
    1481     450   670      18   240  2.63e-14     77.0     49.79  31.2692075
    1482     455   658       9   210  2.67e-14     76.6     48.11  31.2541128
    1483     463   656     164   360  3.06e-14     77.8     48.28  31.1177764
    1484     450   670      14   236  3.07e-14     76.6     49.79  31.1145137
    1485     550   714     198   356  3.21e-14     76.6     51.18  31.0699204
    1486     450   670      25   247  3.23e-14     76.6     49.79  31.0637092
    1487     450   670      28   250  3.31e-14     76.6     49.79  31.0392431
    1488     450   670      25   247  3.57e-14     76.6     49.79  30.9636257
    1489     450   670      18   240  3.65e-14     76.6     49.79  30.9414641
    1490     450   670      25   247  3.67e-14     76.6     49.79  30.9359996
    1491     459   704      15   259  3.68e-14     75.5     48.83  30.9332785
    1492     450   670      19   241  3.74e-14     76.6     49.79  30.9171057
    1493     553   714     195   350  3.75e-14     76.6     50.30  30.9144355
    1494     553   714     202   357  3.83e-14     76.6     50.30  30.8933265
    1495     450   670      14   236  3.88e-14     76.3     49.79  30.8803561
    1496     450   670      18   240  3.89e-14     76.3     49.79  30.8777821
    1497     450   670      37   259  3.89e-14     76.6     49.79  30.8777821
    1498     450   670      18   240  3.96e-14     76.3     49.79  30.8599473
    1499     460   661      36   238  3.99e-14     75.9     47.34  30.8524001
    1500     450   670      38   260  3.99e-14     76.6     49.79  30.8524001
    1501     553   714     193   348  4.07e-14     76.3     50.30  30.8325483
    1502     450   670      25   247  4.16e-14     76.6     49.79  30.8106762
    1503     553   714     200   355  4.16e-14     76.3     50.30  30.8106762
    1504     450   670      23   245  4.17e-14     76.3     49.79  30.8082753
    1505     450   670      13   235  4.25e-14     76.3     49.79  30.7892723
    1506     450   670      20   242  4.49e-14     76.3     49.79  30.7343386
    1507     450   670      18   240  4.57e-14     76.3     49.79  30.7166781
    1508     450   670      18   240  4.57e-14     76.3     49.79  30.7166781
    1509     450   670      18   240  4.61e-14     76.3     49.79  30.7079634
    1510     450   670      18   240  4.62e-14     76.3     49.79  30.7057966
    1511     460   667      37   243  4.63e-14     75.5     47.89  30.7036344
    1512     450   670      20   242  4.91e-14     76.3     49.79  30.6449174
    1513     450   670      15   237  4.97e-14     75.9     49.79  30.6327715
    1514     450   670      18   240  5.01e-14     76.3     49.79  30.6247554
    1515     450   670      18   240  5.01e-14     75.9     49.79  30.6247554
    1516     450   670      39   261  5.13e-14     76.3     49.79  30.6010856
    1517     450   670      18   240  5.19e-14     76.3     49.79  30.5894576
    1518     450   670      18   240  5.19e-14     76.3     49.37  30.5894576
    1519     450   670      30   252  5.24e-14     76.3     49.79  30.5798698
    1520     450   670      14   236  5.29e-14     75.9     49.79  30.5703731
    1521     450   670      39   261  5.33e-14     76.3     49.79  30.5628401
    1522     450   670      37   259  5.34e-14     76.3     49.79  30.5609656
    1523     450   670      20   242  5.38e-14     75.9     49.79  30.5535029
    1524     450   670      18   240  5.38e-14     75.9     49.79  30.5535029
    1525     450   670      17   239  5.38e-14     75.9     49.79  30.5535029
    1526     450   670      30   252  5.43e-14     76.3     49.79  30.5442522
    1527     450   670      30   252  5.43e-14     76.3     49.79  30.5442522
    1528     450   670      16   238  5.44e-14     75.9     49.79  30.5424122
    1529     450   670      23   245  5.47e-14     76.3     49.79  30.5369127
    1530     450   670      18   240  5.48e-14     75.9     49.79  30.5350862
    1531     450   670      18   240  5.48e-14     75.9     49.79  30.5350862
    1532     450   670      18   240  5.53e-14     75.9     49.79  30.5260035
    1533     450   670      14   236  5.59e-14     75.9     49.79  30.5152120
    1534     459   656      10   207  5.63e-14     74.7     48.29  30.5080819
    1535     450   670      23   245  5.67e-14     75.9     49.79  30.5010022
    1536     450   670      20   242  5.68e-14     75.9     49.79  30.4992401
    1537     450   670      20   242  5.78e-14     75.9     49.79  30.4817876
    1538     450   670      18   240  5.84e-14     75.9     49.79  30.4714605
    1539     450   670      15   237  5.85e-14     75.9     49.79  30.4697496
    1540     450   670      25   247  5.87e-14     75.9     49.79  30.4663367
    1541     450   670      38   260  5.88e-14     76.3     49.79  30.4646345
    1542     450   670      29   251  5.95e-14     75.9     49.79  30.4528001
    1543     450   670      45   267  6.03e-14     76.3     49.79  30.4394443
    1544     450   670      21   243  6.10e-14     75.9     49.79  30.4279025
    1545     554   721     334   514  6.11e-14     77.4     49.46  30.4262645
    1546     450   670      38   260  6.21e-14     75.9     49.37  30.4100304
    1547     450   670      23   245  6.37e-14     75.9     49.79  30.3845918
    1548     450   670      22   244  6.44e-14     75.9     49.37  30.3736628
    1549     450   670      24   246  6.54e-14     75.9     49.79  30.3582541
    1550     459   656       9   206  6.59e-14     74.7     48.29  30.3506380
    1551     450   670      24   246  6.60e-14     75.9     49.79  30.3491217
    1552     450   670      29   251  6.63e-14     75.9     49.79  30.3445865
    1553     459   656      10   207  6.64e-14     74.7     48.29  30.3430793
    1554     459   656      33   230  6.80e-14     75.9     48.29  30.3192687
    1555     450   670      20   242  7.06e-14     75.9     49.37  30.2817463
    1556     450   670      18   240  7.19e-14     75.5     49.79  30.2635001
    1557     450   670      18   240  7.19e-14     75.5     49.37  30.2635001
    1558     450   670      41   263  7.20e-14     75.9     49.79  30.2621103
    1559     461   650      25   241  7.23e-14     75.1     45.95  30.2579523
    1560     450   670      23   245  7.33e-14     75.5     49.79  30.2442158
    1561     450   670      18   240  7.39e-14     75.5     49.79  30.2360636
    1562     460   667      15   221  7.44e-14     74.7     47.42  30.2293205
    1563     450   670      19   241  7.45e-14     75.5     49.79  30.2279773
    1564     450   670      41   263  7.53e-14     75.9     49.79  30.2172963
    1565     450   670      41   263  7.67e-14     75.9     49.79  30.1988747
    1566     450   670      18   240  7.73e-14     75.5     49.79  30.1910824
    1567     450   670      18   240  7.81e-14     75.5     49.79  30.1807863
    1568     450   670      18   240  8.02e-14     75.5     49.79  30.1542529
    1569     450   670      42   264  8.38e-14     75.9     49.79  30.1103434
    1570     450   670      41   263  8.39e-14     75.9     49.79  30.1091508
    1571     450   670      24   246  9.31e-14     75.5     49.79  30.0051022
    1572     460   667      14   220  9.44e-14     74.3     47.42  29.9912353
    1573     461   666       8   211  9.52e-14     74.3     48.57  29.9827965
    1574     460   667      38   244  9.87e-14     74.7     47.42  29.9466914
    1575     461   666       8   211  9.88e-14     74.3     48.57  29.9456788
    1576     450   670      20   242  9.95e-14     75.1     49.79  29.9386188
    1577     460   661      37   239  1.02e-13     74.7     47.34  29.9138036
    1578     450   670      18   240  1.05e-13     75.1     49.79  29.8848160
    1579     450   670      18   246  1.07e-13     75.1     48.52  29.8659476
    1580     450   660       3   207  1.08e-13     74.7     49.54  29.8566452
    1581     460   667      38   244  1.08e-13     74.3     47.42  29.8566452
    1582     450   670      23   245  1.09e-13     75.1     49.79  29.8474285
    1583     460   667      37   243  1.09e-13     74.3     47.42  29.8474285
    1584     460   661      37   239  1.09e-13     74.3     46.86  29.8474285
    1585     460   667      13   219  1.10e-13     73.9     47.42  29.8382960
    1586     456   722      19   294  1.10e-13     75.1     45.26  29.8382960
    1587     460   667      42   248  1.11e-13     74.3     47.42  29.8292462
    1588     460   667      12   218  1.12e-13     73.9     47.42  29.8202775
    1589     450   670      27   249  1.12e-13     75.1     49.79  29.8202775
    1590     450   670      18   240  1.14e-13     75.1     49.37  29.8025779
    1591     456   722      20   295  1.19e-13     74.7     45.26  29.7596529
    1592     460   667      40   246  1.22e-13     74.3     47.42  29.7347554
    1593     460   661      37   239  1.24e-13     74.3     46.86  29.7184948
    1594     450   670      24   246  1.25e-13     75.1     49.79  29.7104627
    1595     460   667      35   241  1.26e-13     74.3     47.42  29.7024945
    1596     460   667      35   241  1.27e-13     74.3     47.42  29.6945893
    1597     450   670      24   246  1.28e-13     75.1     49.79  29.6867461
    1598     456   722      33   308  1.29e-13     74.3     45.26  29.6789640
    1599     460   661      36   238  1.31e-13     74.3     46.86  29.6635791
    1600     450   670      41   263  1.32e-13     75.1     49.37  29.6559745
    1601     461   666       8   211  1.34e-13     73.9     48.57  29.6409366
    1602     466   649      35   231  1.35e-13     74.3     47.55  29.6335016
    1603     450   670      46   268  1.35e-13     75.1     49.79  29.6335016
    1604     460   667      37   243  1.40e-13     73.9     47.42  29.5971340
    1605     460   667      35   241  1.41e-13     73.9     47.42  29.5900165
    1606     450   670      46   268  1.46e-13     75.1     49.79  29.5551698
    1607     460   667      37   243  1.49e-13     73.9     47.42  29.5348301
    1608     456   722      17   292  1.49e-13     73.9     45.26  29.5348301
    1609     462   657      35   228  1.55e-13     73.6     50.50  29.4953513
    1610     450   670      18   240  1.57e-13     74.7     49.79  29.4825306
    1611     455   672       5   228  1.60e-13     73.9     44.25  29.4636026
    1612     420   669      20   279  1.64e-13     75.1     45.90  29.4389100
    1613     463   704      34   322  1.66e-13     74.3     45.24  29.4267886
    1614     450   670      39   261  1.67e-13     74.7     49.37  29.4207826
    1615     456   722      18   293  1.83e-13     73.9     45.26  29.3292902
    1616     455   656      19   228  2.22e-13     75.5     47.51  29.1360990
    1617     460   667      19   225  2.53e-13     73.2     47.42  29.0053869
    1618     420   669      19   278  2.56e-13     74.3     45.90  28.9935990
    1619     479   668      31   225  2.60e-13     72.8     51.96  28.9780948
    1620     452   651      12   225  2.77e-13     73.2     47.98  28.9147589
    1621     450   670      25   247  2.82e-13     73.9     48.95  28.8968693
    1622     522   701     781   969  3.18e-13     75.5     45.88  28.7767250
    1623     450   670      18   240  3.29e-13     73.6     49.37  28.7427186
    1624     479   668      44   238  3.83e-13     72.8     51.96  28.5907414
    1625     463   669      17   233  3.94e-13     72.8     50.42  28.5624255
    1626     456   722      33   308  4.10e-13     72.8     45.26  28.5226192
    1627     455   656      15   224  4.21e-13     74.7     47.51  28.4961436
    1628     455   656      13   222  4.31e-13     74.7     47.51  28.4726683
    1629     463   669      17   233  4.49e-13     72.4     50.42  28.4317535
    1630     460   667      34   240  4.53e-13     72.4     47.42  28.4228843
    1631     460   667      35   241  4.55e-13     72.4     46.95  28.4184790
    1632     520   650     792   932  4.63e-13     75.1     49.65  28.4010493
    1633     460   667      34   240  4.78e-13     72.4     47.42  28.3691657
    1634     463   665      12   220  5.19e-13     72.4     47.73  28.2868725
    1635     451   711      14   279  6.10e-13     73.2     44.59  28.1253174
    1636     463   699      27   256  6.66e-13     72.4     47.18  28.0374867
    1637     463   699      25   254  6.84e-13     72.4     47.18  28.0108185
    1638     450   665      16   230  6.84e-13     72.4     49.57  28.0108185
    1639     463   699      25   254  6.96e-13     72.4     47.18  27.9934267
    1640     455   656      19   228  7.18e-13     73.9     47.06  27.9623068
    1641     463   699      26   255  7.30e-13     72.4     45.97  27.9457319
    1642     463   699      26   255  7.36e-13     72.4     45.97  27.9375463
    1643     463   699      26   255  7.43e-13     72.4     45.97  27.9280804
    1644     463   699      16   245  7.44e-13     72.4     45.97  27.9267354
    1645     450   665      18   232  7.50e-13     72.4     49.57  27.9187032
    1646     462   657      35   228  7.51e-13     71.6     50.00  27.9173707
    1647     479   668      44   238  7.86e-13     71.6     51.47  27.8718196
    1648     445   669       1   228  8.37e-13     72.0     45.87  27.8089523
    1649     461   665      31   235  8.62e-13     72.4     49.10  27.7795211
    1650     463   699      26   255  9.23e-13     72.0     47.18  27.7111472
    1651     445   669       2   229  9.40e-13     72.0     45.87  27.6928965
    1652     463   665      14   222  9.89e-13     71.2     47.27  27.6420821
    1653     463   665      12   220  1.01e-12     71.6     47.27  27.6210708
    1654     445   669       3   230  1.03e-12     72.0     45.87  27.6014623
    1655     445   669       2   229  1.03e-12     72.0     45.87  27.6014623
    1656     445   669       4   231  1.05e-12     72.0     45.87  27.5822310
    1657     525   649     799   937  1.07e-12     73.9     55.32  27.5633625
    1659     463   665      26   231  1.08e-12     71.2     48.13  27.5540601
    1660     445   669       6   233  1.08e-12     72.0     45.87  27.5540601
    1661     463   705      81   314  1.17e-12     71.6     46.77  27.4740174
    1662     462   668      10   215  1.23e-12     71.2     48.82  27.4240069
    1663     463   665      27   229  1.27e-12     72.0     49.55  27.3920042
    1664     444   658      12   254  1.35e-12     71.6     47.97  27.3309165
    1665     456   658      37   242  1.43e-12     71.6     48.58  27.2733467
    1666     463   665      12   220  1.47e-12     71.2     47.27  27.2457587
    1667     461   665      26   230  1.48e-12     71.6     49.55  27.2389790
    1668     459   708      14   272  1.54e-12     71.6     44.65  27.1992387
    1669     463   665      13   218  1.56e-12     70.5     48.13  27.1863353
    1670     450   665      24   238  1.58e-12     71.6     49.57  27.1735963
    1671     456   665      17   238  1.88e-12     70.9     45.13  26.9997493
    1672     445   669       3   230  1.89e-12     71.2     46.69  26.9944443
    1673     463   665      12   220  1.93e-12     70.5     47.27  26.9735011
    1674     456   665      10   231  2.04e-12     70.9     45.13  26.9180713
    1675     548   721     333   515  2.06e-12     72.4     47.12  26.9083151
    1676     561   711     125   273  2.08e-12     70.9     49.35  26.8986532
    1677     456   665      11   232  2.17e-12     70.9     45.13  26.8562939
    1678     456   665      11   232  2.19e-12     70.9     45.13  26.8471196
    1679     443   658       9   214  2.25e-12     70.1     45.70  26.8200909
    1680     561   711     126   274  2.27e-12     70.5     49.35  26.8112413
    1681     443   658       9   214  2.31e-12     70.1     45.70  26.7937736
    1682     456   665      11   232  2.34e-12     70.5     45.13  26.7808702
    1683     443   658      10   215  2.34e-12     70.1     45.70  26.7808702
    1684     443   658      10   215  2.36e-12     70.1     45.70  26.7723595
    1685     463   663      31   231  2.44e-12     70.9     45.89  26.7390231
    1686     443   658      11   216  2.45e-12     70.1     45.70  26.7349331
    1687     463   663      31   231  2.59e-12     70.9     45.89  26.6793632
    1688     450   665      20   233  2.69e-12     70.9     49.15  26.6414799
    1689     450   665      24   238  2.73e-12     70.9     49.57  26.6267195
    1690     443   658      10   215  3.01e-12     69.7     45.70  26.5290810
    1691     561   711     137   285  3.06e-12     70.1     49.35  26.5126062
    1692     561   711     134   282  3.09e-12     70.1     49.35  26.5028500
    1693     561   711     125   273  3.10e-12     70.1     49.35  26.4996190
    1694     561   711     137   285  3.15e-12     70.1     49.35  26.4836187
    1695     441   711       1   280  3.15e-12     70.9     43.87  26.4836187
    1696     463   661      32   235  3.18e-12     70.5     41.63  26.4741399
    1697     461   665      27   229  3.27e-12     70.1     49.76  26.4462311
    1698     448   669       1   225  3.29e-12     70.5     46.03  26.4401336
    1699     561   711     129   277  3.32e-12     70.1     49.35  26.4310563
    1700     561   711     129   277  3.35e-12     70.1     49.35  26.4220608
    1701     459   669       4   224  3.35e-12     70.5     46.49  26.4220608
    1702     561   711     126   274  3.40e-12     70.1     49.35  26.4072457
    1703     450   711      17   284  3.46e-12     70.1     44.44  26.3897525
    1704     561   711     126   274  3.49e-12     70.1     49.35  26.3811194
    1705     427   656       8   239  3.50e-12     70.5     44.86  26.3782581
    1706     427   656       7   238  3.51e-12     70.5     44.86  26.3754051
    1707     456   658      36   244  3.66e-12     70.1     49.07  26.3335580
    1708     561   711     170   318  3.92e-12     70.1     49.35  26.2649295
    1709     450   665      38   251  4.12e-12     70.5     49.15  26.2151680
    1710     457   651      10   225  4.15e-12     69.7     47.09  26.2079128
    1711     451   711      16   281  4.17e-12     70.5     44.26  26.2031051
    1712     463   665      17   228  4.25e-12     69.7     46.19  26.1841021
    1713     427   656       7   238  4.28e-12     70.1     44.86  26.1770681
    1714     461   665      27   229  4.40e-12     69.7     49.76  26.1494166
    1715     450   665      19   232  4.50e-12     70.1     49.15  26.1269437
    1716     451   711      16   281  4.61e-12     70.5     45.52  26.1027933
    1717     457   651       9   224  4.69e-12     69.7     47.09  26.0855885
    1718     463   665      12   223  4.76e-12     69.3     46.19  26.0707734
    1719     451   658      12   243  4.77e-12     69.7     48.51  26.0686748
    1720     463   665      19   230  4.85e-12     69.3     46.19  26.0520424
    1721     451   711      14   279  4.86e-12     70.1     44.26  26.0499827
    1722     561   711     129   277  5.15e-12     69.3     49.35  25.9920244
    1723     451   711      13   278  5.20e-12     70.1     44.26  25.9823625
    1724     451   711      14   279  5.25e-12     70.1     44.26  25.9727930
    1725     427   656       7   238  5.31e-12     69.7     44.44  25.9614293
    1726     427   656       7   238  5.41e-12     69.7     44.44  25.9427720
    1727     457   651      11   226  5.53e-12     69.3     47.09  25.9208333
    1728     427   656       8   239  5.81e-12     69.7     44.44  25.8714405
    1729     427   656      11   242  5.95e-12     69.7     44.44  25.8476299
    1730     427   656       7   238  6.14e-12     69.7     44.44  25.8161964
    1731     455   656      14   223  6.18e-12     70.9     47.96  25.8097028
    1732     500   662     234   392  6.38e-12     70.9     47.88  25.7778530
    1733     451   711      14   279  6.58e-12     69.7     43.48  25.7469864
    1734     427   656       7   238  6.60e-12     69.7     44.44  25.7439515
    1735     455   656      15   224  6.69e-12     70.9     47.96  25.7304072
    1736     427   656      13   244  6.90e-12     69.7     44.44  25.6994997
    1737     463   650      43   232  6.93e-12     69.7     47.78  25.6951613
    1738     427   656      28   259  6.93e-12     69.7     44.86  25.6951613
    1739     427   656      14   245  6.95e-12     69.7     44.44  25.6922795
    1740     450   665      39   252  7.04e-12     69.7     48.73  25.6794129
    1741     463   677      10   223  7.27e-12     68.9     47.73  25.6472648
    1742     451   711      14   279  8.23e-12     69.3     43.81  25.5232351
    1743     495   710      72   282  8.37e-12     68.6     45.45  25.5063672
    1744     427   656       9   240  8.39e-12     69.3     44.90  25.5039806
    1745     427   656       7   238  9.13e-12     68.9     43.85  25.4194554
    1746     388   665      64   357  9.74e-12     69.7     42.67  25.3547800
    1747     456   731      35   309  9.84e-12     68.6     42.28  25.3445654
    1748     463   657      62   263  9.90e-12     68.9     48.56  25.3384864
    1749     427   656       7   238  9.99e-12     68.9     44.44  25.3294365
    1750     442   656       9   225  1.01e-11     68.9     46.05  25.3184857
    1751     519   705      90   269  1.07e-11     68.9     48.95  25.2607774
    1752     554   650     167   270  1.08e-11     68.6     56.88  25.2514750
    1753     519   663      79   217  1.10e-11     68.9     52.41  25.2331258
    1754     519   663      82   220  1.12e-11     68.9     52.41  25.2151073
    1755     388   665      78   371  1.13e-11     69.3     42.67  25.2062184
    1756     427   656       8   239  1.17e-11     68.9     44.44  25.1714323
    1757     519   663      80   218  1.21e-11     68.9     52.41  25.1378157
    1758     427   656       7   238  1.22e-11     68.6     43.85  25.1295852
    1759     519   663     122   260  1.24e-11     68.9     52.41  25.1133246
    1760     461   665      11   218  1.33e-11     68.2     52.11  25.0432571
    1761     461   677      16   231  1.33e-11     68.2     48.46  25.0432571
    1762     519   663      76   214  1.34e-11     68.6     52.41  25.0357664
    1763     427   656       7   238  1.36e-11     68.6     45.90  25.0209513
    1764     427   656       7   238  1.36e-11     68.6     45.90  25.0209513
    1765     519   663      75   213  1.37e-11     68.6     52.41  25.0136253
    1766     427   656       8   239  1.38e-11     68.6     44.49  25.0063525
    1767     427   656       7   238  1.39e-11     68.6     45.90  24.9991323
    1768     519   663      75   213  1.41e-11     68.6     52.41  24.9848463
    1769     449   651      48   279  1.45e-11     69.3     46.31  24.9568725
    1770     463   650      44   233  1.47e-11     68.6     47.78  24.9431736
    1771     430   656      10   238  1.47e-11     68.6     44.58  24.9431736
    1772     519   663     322   460  1.48e-11     69.7     52.41  24.9363939
    1774     461   656      37   233  1.51e-11     68.2     45.41  24.9163264
    1775     430   656      10   238  1.54e-11     68.6     44.58  24.8966536
    1776     500   656     233   385  1.55e-11     69.7     48.43  24.8901811
    1777     519   663      75   213  1.55e-11     68.2     52.41  24.8901811
    1778     429   656       2   231  1.55e-11     68.2     45.23  24.8901811
    1779     427   656       7   238  1.57e-11     68.6     44.44  24.8773604
    1780     501   714      57   272  1.58e-11     68.6     47.77  24.8710112
    1781     429   656       2   231  1.60e-11     68.2     45.23  24.8584324
    1782     501   714      57   272  1.68e-11     68.6     47.77  24.8096422
    1783     561   711     137   285  1.79e-11     67.8     48.70  24.7462204
    1784     463   650      21   210  1.83e-11     68.2     48.28  24.7241201
    1785     463   666      33   243  1.84e-11     67.8     46.95  24.7186705
    1786     451   658       5   214  1.87e-11     67.4     44.95  24.7024976
    1787     461   677      13   228  1.88e-11     67.4     47.75  24.6971642
    1788     427   656      28   259  1.89e-11     68.2     45.27  24.6918592
    1789     463   650      51   240  1.96e-11     68.2     47.78  24.6554915
    1790     442   656      22   238  2.05e-11     68.2     46.05  24.6105962
    1791     479   668      45   239  2.05e-11     68.2     52.94  24.6105962
    1792     442   656      22   238  2.15e-11     67.8     46.05  24.5629682
    1793     442   656      22   238  2.21e-11     67.8     46.05  24.5354435
    1794     451   658       2   211  2.21e-11     67.0     44.95  24.5354435
    1795     451   692      14   252  2.22e-11     68.2     45.39  24.5309288
    1796     442   656      22   238  2.25e-11     67.8     46.05  24.5175058
    1797     442   656      65   281  2.28e-11     68.2     46.49  24.5042606
    1798     463   650      21   210  2.33e-11     67.8     47.78  24.4825678
    1799     451   658       6   215  2.49e-11     67.0     44.95  24.4161533
    1800     461   656      10   206  2.57e-11     67.4     45.89  24.3845301
    1801     442   656      11   227  2.66e-11     67.4     46.05  24.3501099
    1802     461   658      11   209  2.69e-11     66.6     45.75  24.3388948
    1803     442   656      11   227  2.71e-11     67.4     46.05  24.3314874
    1804     442   656      17   233  2.76e-11     67.4     46.05  24.3132053
    1805     442   656      10   226  2.79e-11     67.4     46.05  24.3023944
    1806     442   656      11   227  2.81e-11     67.4     46.05  24.2952515
    1807     463   650      18   207  2.83e-11     67.4     47.78  24.2881593
    1808     463   650      18   207  2.91e-11     67.4     47.78  24.2602829
    1809     461   656      19   215  2.91e-11     67.0     46.19  24.2602829
    1810     427   656       7   238  2.91e-11     67.4     44.44  24.2602829
    1811     427   656       7   238  2.91e-11     67.4     44.31  24.2602829
    1812     427   656       8   239  2.94e-11     67.4     44.03  24.2500264
    1813     427   656       7   238  2.97e-11     67.4     44.44  24.2398741
    1814     451   692      52   290  3.19e-11     68.2     45.39  24.1684151
    1815     427   656      28   259  3.20e-11     67.8     44.86  24.1652852
    1816     427   656      11   242  3.34e-11     67.4     44.44  24.1224652
    1817     463   650      11   200  3.41e-11     66.6     47.78  24.1017237
    1818     427   656       8   239  3.42e-11     67.4     44.44  24.0987955
    1819     427   656       7   238  3.58e-11     67.4     44.03  24.0530732
    1820     461   677      16   231  3.61e-11     66.6     47.75  24.0447283
    1821     427   656       7   238  3.62e-11     67.4     43.85  24.0419620
    1822     451   692      14   252  3.79e-11     67.4     45.39  23.9960700
    1823     451   692      14   252  3.80e-11     67.4     45.39  23.9934350
    1824     461   656      19   215  3.82e-11     67.0     45.41  23.9881856
    1825     442   656      14   230  3.87e-11     67.0     46.05  23.9751815
    1826     461   656       9   205  3.90e-11     66.6     45.89  23.9674595
    1827     461   656      10   206  3.94e-11     66.6     45.89  23.9572553
    1828     427   656       7   238  3.96e-11     67.0     43.85  23.9521920
    1829     461   676       8   222  4.00e-11     66.6     46.82  23.9421417
    1830     461   671      12   221  4.06e-11     66.6     46.98  23.9272530
    1831     427   656       7   238  4.10e-11     67.0     43.44  23.9174490
    1832     451   692       7   245  4.20e-11     67.0     45.39  23.8933515
    1833     461   677      12   227  4.20e-11     66.6     48.66  23.8933515
    1834     427   656       8   239  4.21e-11     67.0     44.26  23.8909734
    1835     451   664      14   233  4.23e-11     67.4     45.83  23.8862340
    1836     461   656      28   224  4.26e-11     66.6     45.41  23.8791669
    1837     461   671       8   217  4.43e-11     66.2     46.98  23.8400364
    1838     457   664     190   398  4.49e-11     68.2     45.41  23.8265833
    1839     463   650      63   252  4.51e-11     67.4     47.78  23.8221389
    1840     451   664      14   233  4.59e-11     67.0     45.83  23.8045560
    1841     461   677      13   228  4.65e-11     66.2     47.75  23.7915688
    1842     427   656       8   239  4.68e-11     67.0     44.08  23.7851379
    1843     451   711      19   284  4.81e-11     67.0     43.48  23.7577389
    1844     451   692      14   252  4.85e-11     67.0     45.39  23.7494573
    1845     461   677      16   231  4.92e-11     66.2     47.75  23.7351275
    1846     461   677       9   224  4.98e-11     66.2     47.75  23.7230061
    1847     451   692      14   252  5.03e-11     67.0     45.39  23.7130160
    1848     451   692      15   253  5.04e-11     67.0     45.39  23.7110299
    1849     451   692      15   253  5.06e-11     67.0     45.39  23.7070695
    1850     451   711       8   273  5.09e-11     67.0     43.48  23.7011582
    1851     427   656       8   240  5.11e-11     66.6     43.90  23.6972366
    1852     451   692      13   251  5.13e-11     67.0     45.39  23.6933304
    1853     461   677       8   223  5.16e-11     66.2     48.05  23.6874994
    1854     451   692      52   290  5.21e-11     67.4     45.39  23.6778562
    1855     451   692      17   255  5.29e-11     67.0     45.39  23.6626178
    1856     461   677      11   226  5.37e-11     66.2     47.75  23.6476081
    1857     451   711      14   279  5.38e-11     67.0     43.48  23.6457476
    1858     461   677      12   227  5.42e-11     66.2     47.75  23.6383402
    1859     451   692      14   252  5.43e-11     66.6     45.39  23.6364969
    1860     451   692      16   254  5.45e-11     67.0     45.39  23.6328204
    1861     451   692      14   252  5.48e-11     66.6     45.39  23.6273309
    1862     451   692      15   253  5.51e-11     66.6     45.39  23.6218714
    1863     451   711      15   280  5.53e-11     66.6     43.48  23.6182482
    1864     451   711      10   275  5.59e-11     66.6     43.48  23.6074567
    1865     461   677      11   226  5.62e-11     66.2     47.75  23.6021044
    1866     461   677      11   226  5.73e-11     66.2     47.75  23.5827205
    1867     461   677      12   227  5.73e-11     66.2     47.75  23.5827205
    1868     461   677      13   228  5.75e-11     66.2     47.75  23.5792362
    1869     451   692       7   245  5.76e-11     66.6     45.39  23.5774985
    1870     461   677      12   227  5.78e-11     66.2     48.21  23.5740323
    1871     461   677      16   231  5.79e-11     66.2     47.75  23.5723037
    1872     461   677       9   224  5.82e-11     65.9     47.75  23.5671358
    1873     455   658      12   220  5.91e-11     66.6     44.89  23.5517902
    1874     500   656     233   385  6.01e-11     67.8     47.80  23.5350113
    1875     461   677      13   228  6.06e-11     65.9     47.75  23.5267262
    1876     427   656       8   239  6.08e-11     66.6     44.26  23.5234313
    1877     451   664      14   233  6.10e-11     66.6     45.42  23.5201473
    1878     451   692       8   246  6.17e-11     66.6     45.39  23.5087372
    1879     461   656      12   208  6.21e-11     67.4     45.71  23.5022751
    1880     461   677      12   227  6.28e-11     65.9     47.83  23.4910660
    1881     461   677      10   225  6.50e-11     65.9     47.75  23.4566338
    1882     422   656      57   293  6.52e-11     67.0     44.76  23.4535616
    1883     451   711      24   289  6.55e-11     66.6     43.48  23.4489710
    1884     444   658      14   256  6.63e-11     66.2     46.34  23.4368312
    1885     461   677      12   227  6.63e-11     65.9     48.21  23.4368312
    1886     463   658      20   220  6.64e-11     66.2     44.70  23.4353241
    1887     461   677       9   224  6.73e-11     65.9     47.75  23.4218609
    1888     461   677      10   225  6.81e-11     65.9     48.21  23.4100439
    1889     461   677      11   226  6.87e-11     65.9     47.75  23.4012719
    1890     461   677       8   223  6.92e-11     65.9     47.75  23.3940203
    1891     461   677      13   228  7.07e-11     65.9     47.75  23.3725755
    1892     461   677       9   224  7.11e-11     65.9     47.75  23.3669338
    1893     461   677      11   226  7.19e-11     65.9     48.21  23.3557449
    1894     461   656      19   215  7.21e-11     67.0     46.19  23.3529671
    1895     461   677       7   222  7.24e-11     65.9     47.75  23.3488148
    1896     461   677      10   225  7.25e-11     65.9     47.75  23.3474346
    1897     422   656      57   293  7.31e-11     66.6     44.35  23.3391927
    1898     422   656      57   293  7.31e-11     66.6     44.35  23.3391927
    1899     461   677      10   225  7.39e-11     65.9     47.75  23.3283083
    1900     461   677      12   227  7.40e-11     65.9     47.30  23.3269560
    1901     455   658      11   219  7.43e-11     65.5     44.89  23.3229102
    1902     422   656      57   293  7.44e-11     66.6     44.35  23.3215652
    1903     461   677      10   225  7.46e-11     65.9     47.75  23.3188806
    1904     500   656     233   385  7.64e-11     67.4     47.80  23.2950384
    1905     500   656     234   386  7.65e-11     67.4     47.80  23.2937304
    1906     422   656      57   293  7.71e-11     66.6     44.35  23.2859178
    1907     461   677       8   223  7.79e-11     65.5     47.75  23.2755952
    1908     461   677       9   224  7.80e-11     65.5     47.83  23.2743123
    1909     461   677       9   224  7.80e-11     65.5     47.75  23.2743123
    1910     461   677       9   224  7.93e-11     65.5     47.75  23.2577830
    1911     427   656       8   239  7.96e-11     66.2     44.86  23.2540070
    1912     461   677       8   223  8.00e-11     65.5     47.75  23.2489945
    1913     461   677       9   224  8.02e-11     65.5     47.75  23.2464976
    1914     461   677      11   226  8.03e-11     65.5     47.75  23.2452515
    1915     451   711      10   275  8.05e-11     66.2     43.48  23.2427639
    1916     461   677       8   223  8.08e-11     65.5     47.75  23.2390442
    1917     461   677       9   224  8.15e-11     65.5     47.75  23.2304181
    1918     461   677       9   224  8.24e-11     65.5     47.75  23.2194357
    1919     427   656       8   239  8.25e-11     66.2     43.85  23.2182228
    1920     463   650      21   210  8.28e-11     66.2     47.29  23.2145931
    1921     461   656       6   202  8.29e-11     67.0     45.71  23.2133861
    1922     461   677       9   224  8.31e-11     65.5     47.75  23.2109764
    1923     461   677       9   224  8.31e-11     65.5     47.75  23.2109764
    1924     427   656       7   238  8.49e-11     66.2     43.44  23.1895470
    1925     427   656       8   239  8.55e-11     66.2     44.86  23.1825047
    1926     427   656       8   239  8.55e-11     66.2     43.85  23.1825047
    1927     461   677      16   231  8.65e-11     65.5     47.75  23.1708767
    1928     461   677       8   223  8.70e-11     65.5     47.75  23.1651130
    1929     461   677       8   223  8.70e-11     65.5     47.75  23.1651130
    1930     461   677       9   224  8.78e-11     65.5     47.75  23.1559596
    1931     461   677       9   224  8.86e-11     65.5     47.75  23.1468893
    1932     427   656       7   238  8.96e-11     65.9     44.31  23.1356658
    1933     553   710     245   396  9.11e-11     66.6     50.30  23.1190633
    1935     427   656      10   241  9.22e-11     65.9     44.86  23.1070610
    1936     461   677      11   226  9.29e-11     65.5     47.75  23.0994975
    1937     454   665      10   228  9.31e-11     65.5     45.33  23.0973469
    1938     463   650      18   207  9.62e-11     65.9     47.29  23.0645918
    1939     427   656      10   241  9.65e-11     65.9     44.86  23.0614781
    1940     427   656       7   238  9.71e-11     65.9     43.85  23.0552797
    1941     427   656       8   239  9.78e-11     65.9     43.85  23.0480965
    1942     461   677      10   225  9.79e-11     65.5     47.75  23.0470746
    1943     463   649      21   228  1.01e-10     65.9     43.58  23.0159006
    1944     547   658     113   220  1.02e-10     65.9     52.68  23.0060483
    1945     427   656       7   238  1.06e-10     65.9     43.85  22.9675820
    1946     457   664     200   408  1.07e-10     66.6     45.87  22.9581923
    1947     427   656       7   238  1.07e-10     65.9     43.85  22.9581923
    1948     427   656       7   238  1.07e-10     65.9     43.85  22.9581923
    1949     427   656       8   239  1.10e-10     65.9     44.03  22.9305408
    1950     427   656       8   239  1.12e-10     65.9     43.85  22.9125222
    1951     461   658      11   209  1.13e-10     64.7     45.75  22.9036333
    1952     456   659       8   216  1.13e-10     66.6     46.30  22.9036333
    1953     451   664      14   233  1.14e-10     65.9     45.42  22.8948227
    1954     427   656       8   239  1.20e-10     65.5     43.85  22.8435294
    1955     451   664      14   233  1.20e-10     65.9     45.42  22.8435294
    1956     463   649      43   250  1.28e-10     65.5     43.58  22.7789909
    1957     427   656      10   241  1.33e-10     65.5     43.85  22.7406720
    1958     463   650      18   207  1.35e-10     65.5     47.29  22.7257463
    1959     461   677      11   226  1.36e-10     65.1     47.96  22.7183662
    1960     463   650      63   252  1.45e-10     65.9     47.29  22.6542874
    1961     454   665      10   228  1.46e-10     65.1     44.89  22.6474145
    1962     451   664      14   233  1.54e-10     65.5     45.42  22.5940685
    1963     451   664      14   233  1.59e-10     65.5     45.42  22.5621169
    1964     427   656       8   239  1.60e-10     65.1     43.85  22.5558473
    1965     457   664     190   398  1.62e-10     66.2     44.95  22.5434248
    1966     442   656       9   225  1.65e-10     65.1     46.05  22.5250756
    1967     454   665      10   228  1.68e-10     64.7     44.89  22.5070571
    1968     457   664     191   399  1.70e-10     66.2     44.95  22.4952227
    1969     442   653       3   225  1.72e-10     65.1     47.58  22.4835266
    1970     457   664     191   399  1.72e-10     66.2     44.95  22.4835266
    1971     450   656      31   239  1.97e-10     65.1     47.27  22.3478174
    1972     457   664     191   399  1.97e-10     66.2     44.95  22.3478174
    1973     461   677       8   223  1.98e-10     64.3     47.30  22.3427541
    1974     441   694     525   775  2.10e-10     66.2     45.05  22.2839136
    1975     450   656      33   241  2.12e-10     65.1     47.27  22.2744348
    1976     457   664     191   399  2.13e-10     65.9     44.95  22.2697290
    1977     457   664     191   399  2.14e-10     65.9     44.95  22.2650451
    1978     438   660     258   480  2.18e-10     66.2     42.34  22.2465261
    1979     450   656      33   241  2.18e-10     64.7     47.27  22.2465261
    1980     455   658      32   240  2.24e-10     64.3     44.89  22.2193751
    1981     461   676       8   222  2.26e-10     64.3     46.82  22.2104861
    1982     454   665      10   228  2.27e-10     64.3     44.89  22.2060711
    1983     457   664     162   370  2.27e-10     65.9     44.95  22.2060711
    1984     457   664     163   371  2.27e-10     65.9     44.95  22.2060711
    1985     454   665      10   228  2.40e-10     64.3     44.89  22.1503822
    1986     457   664     162   370  2.43e-10     65.9     44.95  22.1379597
    1987     462   665      90   294  2.55e-10     65.1     47.87  22.0897576
    1988     461   677      16   231  2.56e-10     64.3     47.30  22.0858437
    1989     442   656       8   224  2.64e-10     64.3     45.85  22.0550720
    1990     457   664     162   370  2.68e-10     65.5     44.95  22.0400341
    1991     457   664     162   370  2.79e-10     65.5     44.95  21.9998093
    1992     451   670       8   238  2.83e-10     64.7     45.38  21.9855742
    1993     442   656      23   239  2.84e-10     64.3     45.41  21.9820469
    1994     427   656       7   238  2.98e-10     64.3     43.62  21.9339276
    1995     427   656       7   238  3.01e-10     64.3     43.62  21.9239109
    1996     463   653      19   219  3.20e-10     64.3     48.53  21.8627001
    1997     461   677       8   223  3.24e-10     63.9     47.30  21.8502776
    1998     562   710     148   290  3.36e-10     63.5     49.67  21.8139100
    1999     562   710     148   290  3.55e-10     63.5     49.67  21.7589033
    2000     427   656       7   238  3.66e-10     64.3     43.44  21.7283878
    2001     451   670      30   260  3.85e-10     64.3     44.96  21.6777778
    2002     461   677      10   225  3.96e-10     63.5     47.30  21.6496069
    2003     427   656       7   238  4.00e-10     63.9     43.44  21.6395566
    2004     427   656       8   239  4.03e-10     63.9     43.44  21.6320846
    2005     461   677       8   223  4.14e-10     63.5     47.30  21.6051551
    2006     461   677       9   224  4.18e-10     63.5     47.30  21.5955397
    2007     556   711     158   309  5.00e-10     63.5     51.90  21.4164130
    2008     501   717      74   287  5.42e-10     62.8     40.91  21.3357551
    2009     397   711      10   313  5.83e-10     63.5     43.34  21.2628339
    2010     459   675      25   240  5.86e-10     63.5     48.31  21.2577013
    2011     562   710     140   282  5.90e-10     62.8     49.67  21.2508986
    2012     442   656      23   239  6.03e-10     63.5     45.41  21.2291039
    2013     524   656      85   226  6.29e-10     63.2     51.75  21.1868899
    2014     451   670       2   232  6.52e-10     63.2     44.96  21.1509766
    2015     459   675      26   241  6.57e-10     63.2     48.31  21.1433371
    2016     463   665       9   218  6.64e-10     62.8     45.83  21.1327390
    2017     459   675      39   254  6.75e-10     63.5     48.31  21.1163084
    2018     444   672      17   243  6.86e-10     63.2     47.93  21.1001435
    2019     462   651      25   225  7.17e-10     63.2     45.10  21.0559453
    2020     454   672       6   222  7.18e-10     63.2     48.71  21.0545515
    2021     454   672       8   224  7.27e-10     63.2     47.70  21.0420946
    2022     459   669      30   246  7.64e-10     63.2     45.49  20.9924533
    2023     457   664     191   399  7.66e-10     64.3     44.95  20.9898389
    2024     462   651      26   226  7.69e-10     63.2     45.10  20.9859301
    2025     451   670       8   238  8.11e-10     63.2     44.54  20.9327531
    2026     451   670      10   240  8.50e-10     63.2     44.96  20.8857848
    2027     459   675      48   263  8.56e-10     63.2     48.31  20.8787507
    2028     463   656      19   209  9.08e-10     62.8     48.79  20.8197767
    2029     451   670       6   236  9.17e-10     62.8     44.54  20.8099136
    2030     451   670       9   239  9.25e-10     62.8     44.54  20.8012274
    2031     461   665      19   224  9.50e-10     62.4     47.64  20.7745591
    2032     461   677       8   223  9.66e-10     62.4     46.85  20.7578573
    2033     451   664      14   233  9.78e-10     62.8     45.92  20.7455114
    2034     462   651      25   225  9.80e-10     63.2     45.10  20.7434685
    2035     451   664       5   224  9.90e-10     62.8     45.92  20.7333162
    2036     459   669       7   223  1.02e-09     62.8     45.49  20.7034632
    2037     454   672      27   243  1.04e-09     62.4     47.70  20.6840451
    2038     451   670       3   233  1.04e-09     62.8     44.54  20.6840451
    2039     451   664      12   231  1.07e-09     62.8     45.92  20.6556072
    2040     451   670      16   246  1.08e-09     62.8     44.54  20.6463048
    2041     462   651      25   225  1.08e-09     62.4     45.10  20.6463048
    2042     391   655      17   263  1.09e-09     63.2     45.72  20.6370881
    2043     451   670       4   234  1.10e-09     62.8     44.54  20.6279557
    2044     391   655      16   262  1.10e-09     62.8     45.72  20.6279557
    2045     563   665     645   738  1.12e-09     63.9     56.73  20.6099372
    2046     459   675      26   241  1.13e-09     62.4     48.91  20.6010482
    2047     451   670      16   246  1.14e-09     62.8     44.54  20.5922376
    2048     391   655      18   264  1.14e-09     62.8     45.72  20.5922376
    2049     451   670      17   247  1.14e-09     62.8     44.54  20.5922376
    2050     451   670      21   251  1.15e-09     62.8     44.54  20.5835039
    2051     563   665     644   737  1.15e-09     63.9     56.73  20.5835039
    2052     451   670       8   238  1.17e-09     62.8     44.54  20.5662621
    2053     451   664      31   250  1.17e-09     62.8     45.92  20.5662621
    2054     451   670      20   250  1.18e-09     62.8     44.54  20.5577514
    2055     451   670       1   231  1.21e-09     62.4     44.54  20.5326455
    2056     459   669      30   246  1.21e-09     62.8     45.06  20.5326455
    2057     451   670      20   250  1.23e-09     62.8     44.54  20.5162517
    2058     451   670      32   262  1.23e-09     62.8     45.19  20.5162517
    2059     462   651      24   224  1.28e-09     62.4     45.10  20.4764058
    2060     451   670      20   250  1.30e-09     62.4     44.54  20.4609016
    2061     451   670      16   246  1.32e-09     62.4     44.54  20.4456341
    2062     461   656      12   208  1.37e-09     63.2     45.24  20.4084551
    2063     451   670      12   242  1.38e-09     62.4     44.54  20.4011823
    2064     501   712      57   270  1.39e-09     62.4     45.26  20.3939621
    2065     451   670       4   234  1.45e-09     62.4     44.54  20.3517023
    2066     451   670       3   233  1.46e-09     62.4     44.54  20.3448294
    2067     451   670      13   243  1.46e-09     62.4     44.54  20.3448294
    2068     451   670       3   233  1.49e-09     62.4     44.54  20.3244897
    2069     556   711     114   265  1.51e-09     61.6     51.90  20.3111562
    2070     451   670      17   247  1.52e-09     62.4     44.96  20.3045555
    2071     451   670      16   246  1.53e-09     62.4     44.54  20.2979981
    2072     451   670      14   244  1.56e-09     62.4     44.54  20.2785800
    2073     451   670      21   251  1.58e-09     62.4     44.54  20.2658410
    2074     451   670      10   240  1.60e-09     62.4     44.54  20.2532622
    2075     451   670      12   242  1.60e-09     62.4     44.54  20.2532622
    2076     451   664      34   253  1.61e-09     62.4     43.70  20.2470317
    2077     391   655      16   262  1.61e-09     62.4     43.87  20.2470317
    2078     462   651      24   224  1.63e-09     62.0     45.10  20.2346858
    2079     461   670     236   448  1.63e-09     62.8     46.72  20.2346858
    2080     462   651      24   224  1.64e-09     62.0     45.10  20.2285696
    2081     451   670      32   262  1.65e-09     62.4     44.54  20.2224905
    2082     451   670       9   239  1.65e-09     62.0     44.54  20.2224905
    2083     462   651      24   224  1.66e-09     62.0     45.10  20.2164482
    2084     462   651      18   218  1.69e-09     62.0     45.10  20.1985373
    2085     451   670      10   240  1.69e-09     62.0     44.54  20.1985373
    2086     391   655      17   263  1.69e-09     62.4     43.87  20.1985373
    2087     451   670      20   250  1.70e-09     62.4     44.54  20.1926376
    2088     451   670      12   242  1.70e-09     62.4     44.54  20.1926376
    2089     451   670      30   254  1.71e-09     62.4     45.83  20.1867725
    2090     501   712      57   270  1.71e-09     62.0     45.26  20.1867725
    2091     451   670      16   246  1.71e-09     62.0     44.54  20.1867725
    2092     451   670      10   240  1.72e-09     62.0     44.54  20.1809415
    2093     451   670      17   247  1.72e-09     62.0     44.54  20.1809415
    2094     451   670      17   247  1.72e-09     62.0     44.54  20.1809415
    2095     451   670      18   248  1.73e-09     62.0     44.54  20.1751444
    2096     451   670      19   249  1.74e-09     62.0     44.54  20.1693807
    2097     462   651      26   226  1.74e-09     62.0     45.10  20.1693807
    2098     451   670       2   232  1.75e-09     62.0     44.54  20.1636500
    2099     451   670      10   240  1.75e-09     62.0     44.54  20.1636500
    2100     461   677      13   229  1.77e-09     61.6     44.77  20.1522863
    2101     451   670       3   233  1.78e-09     62.0     44.54  20.1466525
    2102     391   655      18   264  1.78e-09     62.4     43.87  20.1466525
    2103     451   655      81   281  1.78e-09     62.4     51.43  20.1466525
    2104     501   712      57   270  1.79e-09     62.0     45.26  20.1410502
    2105     391   655      20   266  1.81e-09     62.4     43.87  20.1299390
    2106     451   670      12   242  1.82e-09     62.0     44.54  20.1244293
    2107     451   670      32   262  1.82e-09     62.0     44.54  20.1244293
    2108     391   655      22   268  1.83e-09     62.4     43.87  20.1189499
    2109     451   664      14   233  1.84e-09     62.4     45.06  20.1135003
    2110     451   670      13   243  1.84e-09     62.0     44.54  20.1135003
    2111     451   670      10   240  1.85e-09     62.0     44.54  20.1080802
    2112     461   656      17   213  1.90e-09     62.4     45.24  20.0814120
    2113     451   670      16   246  1.91e-09     62.0     44.54  20.0761626
    2114     461   677      13   229  1.92e-09     61.6     44.77  20.0709407
    2115     451   670       2   232  1.93e-09     62.0     44.54  20.0657458
    2116     501   712      57   270  1.94e-09     62.0     45.26  20.0605779
    2117     461   677       8   224  1.94e-09     61.6     44.77  20.0605779
    2118     391   655      38   284  2.00e-09     62.4     43.87  20.0301187
    2119     451   670      20   250  2.03e-09     62.0     44.54  20.0152300
    2120     391   655      38   284  2.05e-09     62.4     43.87  20.0054260
    2121     451   670      17   247  2.08e-09     62.0     44.54  19.9908979
    2122     451   670      15   245  2.21e-09     62.0     44.54  19.9302733
    2123     451   761       5   295  2.23e-09     61.2     42.68  19.9212643
    2124     461   677       8   224  2.25e-09     61.2     44.77  19.9123356
    2125     451   670      16   246  2.28e-09     61.6     44.54  19.8990904
    2126     451   670       2   232  2.40e-09     61.6     44.54  19.8477971
    2127     451   670       3   233  2.47e-09     61.6     44.54  19.8190477
    2128     451   670      20   250  2.47e-09     61.6     44.54  19.8190477
    2129     462   651      24   224  2.53e-09     61.2     45.59  19.7950465
    2130     463   653      21   221  2.54e-09     61.6     48.53  19.7911018
    2131     463   670      39   250  2.56e-09     61.6     44.75  19.7832586
    2132     451   670      12   242  2.62e-09     61.6     44.54  19.7600915
    2133     451   670      16   246  2.62e-09     61.6     44.54  19.7600915
    2134     451   670      16   246  2.65e-09     61.6     44.54  19.7487062
    2135     451   670      12   242  2.71e-09     61.6     44.54  19.7263172
    2136     461   670      49   262  2.83e-09     61.6     43.58  19.6829891
    2137     451   670      31   261  2.92e-09     61.6     44.54  19.6516822
    2138     451   670      16   246  3.05e-09     61.2     44.54  19.6081242
    2139     451   670      20   250  3.09e-09     61.2     44.54  19.5950947
    2140     463   721      19   290  3.21e-09     61.2     44.48  19.5569949
    2141     451   670      17   247  3.24e-09     61.2     44.54  19.5476925
    2142     451   670      31   261  3.31e-09     61.2     44.54  19.5263176
    2143     462   651      24   224  3.34e-09     60.8     45.10  19.5172950
    2144     451   670      12   242  3.45e-09     61.2     44.54  19.4848916
    2145     459   675      39   254  3.58e-09     61.2     47.88  19.4479030
    2146     451   670      10   240  3.60e-09     61.2     44.54  19.4423320
    2147     451   670      14   244  3.67e-09     61.2     44.12  19.4230742
    2148     461   670      46   259  4.16e-09     61.2     43.58  19.2977508
    2149     461   670      47   260  4.25e-09     60.8     43.58  19.2763469
    2150     451   670      11   241  4.73e-09     60.8     44.77  19.1693406
    2151     463   653      21   221  4.78e-09     60.5     48.53  19.1588253
    2152     451   655      64   264  4.98e-09     60.8     48.80  19.1178359
    2153     463   653      21   221  5.04e-09     60.5     48.53  19.1058598
    2154     463   653      21   221  5.10e-09     60.5     48.53  19.0940253
    2155     524   656      85   226  5.20e-09     60.5     51.05  19.0746072
    2156     524   656      83   224  5.37e-09     60.1     51.05  19.0424379
    2157     463   653      17   217  5.40e-09     60.1     48.53  19.0368669
    2158     451   761      21   315  6.22e-09     60.1     42.81  18.8954959
    2159     463   721      19   290  6.37e-09     60.5     44.48  18.8716664
    2160     451   670      15   245  6.94e-09     60.1     44.12  18.7859641
    2161     456   665       4   212  7.10e-09     59.7     45.50  18.7631711
    2162     463   721      19   290  7.54e-09     60.1     44.48  18.7030437
    2163     451   670      14   244  7.75e-09     60.1     44.12  18.6755730
    2164     456   665       6   214  8.01e-09     59.7     45.50  18.6425751
    2165     466   656      23   214  9.14e-09     59.7     44.22  18.5106055
    2166     451   655       2   202  1.01e-08     59.3     48.80  18.4107304
    2167     456   705       8   264  1.21e-08     59.7     43.02  18.2300604
    2168     452   651      24   247  1.25e-08     59.7     43.64  18.1975372
    2169     451   670      15   245  1.26e-08     59.3     44.12  18.1895690
    2170     451   655       2   202  1.27e-08     59.3     48.80  18.1816638
    2171     512   656      75   229  1.27e-08     58.9     49.36  18.1816638
    2172     463   653      32   224  1.33e-08     59.3     45.89  18.1355018
    2173     451   655       2   202  1.35e-08     59.3     48.80  18.1205762
    2174     451   655       2   202  1.37e-08     59.3     48.80  18.1058700
    2175     451   655      28   228  1.42e-08     59.3     48.80  18.0700239
    2176     512   656      72   226  1.44e-08     58.9     49.36  18.0560376
    2177     451   655      29   229  1.44e-08     59.3     48.80  18.0560376
    2178     157   224       6    76  1.49e-08     53.9     63.38  18.0219046
    2179     451   655      53   253  1.50e-08     59.3     48.80  18.0152156
    2180     451   655      30   230  1.58e-08     59.3     48.80  17.9632559
    2181     451   655      27   227  1.65e-08     58.9     48.80  17.9199055
    2182     451   655      27   227  1.67e-08     58.9     48.80  17.9078571
    2183     463   701      82   327  1.67e-08     59.3     45.24  17.9078571
    2184     451   655       2   202  1.70e-08     58.9     48.80  17.8900525
    2185     451   655      27   227  1.71e-08     58.9     48.80  17.8841874
    2186     451   655       2   202  1.76e-08     58.9     48.80  17.8553669
    2187     451   655      24   224  1.81e-08     58.9     48.80  17.8273539
    2188     451   655      24   224  1.86e-08     58.9     48.80  17.8001043
    2189     459   659      11   216  2.19e-08     58.9     44.91  17.6367792
    2190     451   655      24   224  2.20e-08     58.5     46.41  17.6322234
    2191     460   650      95   290  2.31e-08     58.9     47.52  17.5834332
    2192     463   657      42   235  2.39e-08     58.5     44.17  17.5493874
    2193     451   655      12   212  2.45e-08     58.5     48.80  17.5245927
    2194     461   646      17   221  2.51e-08     57.8     46.73  17.5003980
    2195     461   646      13   217  2.81e-08     57.8     46.73  17.3874963
    2196     461   646      14   218  2.83e-08     57.8     46.73  17.3804040
    2197     461   646      15   219  2.91e-08     57.8     46.73  17.3525277
    2198     463   657     312   505  3.07e-08     58.9     43.69  17.2990032
    2199     461   646      17   221  3.09e-08     57.8     46.73  17.2925097
    2200     461   646      15   219  3.10e-08     57.8     46.73  17.2892786
    2201     457   660      12   222  3.23e-08     57.4     43.12  17.2481986
    2202     461   646      13   217  3.37e-08     57.4     46.73  17.2057680
    2203     463   741      20   301  3.44e-08     57.8     46.23  17.1852093
    2204     559   721     128   291  3.96e-08     58.2     52.87  17.0444367
    2205     463   657     164   357  3.98e-08     58.2     43.69  17.0393989
    2206     463   650      38   227  4.03e-08     57.8     45.10  17.0269144
    2207     556   710     141   291  4.15e-08     57.4     49.68  16.9975724
    2208     559   721     142   305  4.31e-08     58.2     52.87  16.9597428
    2209     463   722      20   287  4.31e-08     57.8     43.53  16.9597428
    2210     451   655      27   227  4.41e-08     57.4     48.80  16.9368061
    2211     459   705      12   265  4.57e-08     57.8     41.76  16.9011675
    2212     463   701      84   329  4.59e-08     58.2     45.85  16.8968007
    2213     559   721     143   306  4.76e-08     57.8     52.60  16.8604331
    2214     463   701      82   327  5.04e-08     57.8     45.85  16.8032747
    2215     462   665      35   237  5.09e-08     57.4     44.76  16.7934029
    2216     463   650      44   233  5.15e-08     57.8     45.10  16.7816840
    2217     451   655      21   221  5.30e-08     57.0     48.80  16.7529739
    2218     459   656      82   282  5.53e-08     57.8     46.83  16.7104929
    2219     463   701      98   343  5.62e-08     57.8     45.85  16.6943491
    2220     459   656      76   276  6.24e-08     57.4     46.83  16.5897006
    2221     451   655      43   243  6.32e-08     57.0     48.80  16.5769615
    2222     459   656      72   272  7.03e-08     57.4     46.83  16.4704940
    2223     459   656      73   273  7.05e-08     57.4     46.83  16.4676531
    2224     463   676      33   246  7.56e-08     57.0     44.50  16.3978096
    2225     459   656      72   272  8.02e-08     57.0     46.83  16.3387423
    2226     451   661      25   260  8.24e-08     57.4     43.60  16.3116804
    2227     459   656      68   268  8.95e-08     57.0     46.83  16.2290272
    2228     451   655      27   227  8.96e-08     56.2     46.41  16.2279105
    2229     459   656      92   292  9.05e-08     57.0     46.83  16.2179160
    2230     451   655      27   227  9.17e-08     56.6     48.33  16.2047435
    2231     459   656      79   279  9.35e-08     57.0     46.83  16.1853044
    2232     459   659      11   216  1.05e-07     56.2     44.91  16.0693055
    2233     459   656      79   279  1.08e-07     56.6     46.83  16.0411346
    2234     452   659      12   230  1.17e-07     56.2     43.67  15.9610919
    2235     456   705       8   264  1.21e-07     56.6     41.95  15.9274753
    2236     463   715      32   305  1.25e-07     55.8     41.52  15.8949521
    2237     552   704     135   280  1.38e-07     55.8     48.08  15.7960122
    2238     463   652      32   231  1.46e-07     55.8     43.84  15.7396592
    2239     459   656     106   306  1.50e-07     56.2     46.34  15.7126305
    2240     459   656      68   268  1.64e-07     56.2     46.83  15.6233994
    2241     513   656     132   273  2.07e-07     55.8     48.97  15.3905470
    2242     513   656     138   279  2.07e-07     55.8     48.97  15.3905470
    2243     559   705     133   269  2.17e-07     55.5     52.32  15.3433685
    2244     513   656     135   276  2.25e-07     55.8     48.97  15.3071654
    2245     513   656     137   278  2.27e-07     55.8     48.97  15.2983158
    2246     513   656     136   277  2.40e-07     55.8     48.97  15.2426269
    2247     513   656     137   278  2.61e-07     55.5     48.97  15.1587454
    2248     513   656     133   274  2.64e-07     55.5     48.97  15.1473167
    2249     459   659      11   216  2.89e-07     55.1     44.91  15.0568391
    2250     513   656     131   272  2.92e-07     55.5     48.97  15.0465120
    2251     459   659      11   216  3.15e-07     54.7     44.91  14.9706932
    2252     459   659      12   217  3.29e-07     54.7     44.91  14.9272081
    2253     459   659      18   223  3.40e-07     54.7     44.91  14.8943202
    2254     459   659      31   236  3.75e-07     54.7     44.91  14.7963398
    2255     459   659      13   218  3.77e-07     54.3     44.91  14.7910206
    2256     459   659      11   216  3.81e-07     54.3     44.91  14.7804665
    2257     459   659      12   217  3.83e-07     54.3     44.91  14.7752308
    2258     459   659      31   236  3.89e-07     54.3     44.91  14.7596865
    2259     459   659      11   216  4.03e-07     54.3     44.91  14.7243293
    2260     459   659      13   218  4.16e-07     54.3     44.91  14.6925806
    2261     459   659       8   213  4.24e-07     54.3     44.91  14.6735324
    2262     459   659      13   218  4.51e-07     54.3     44.91  14.6117985
    2263     459   659      13   218  4.72e-07     53.9     44.91  14.5662869
    2264     459   659      11   216  6.07e-07     53.9     44.91  14.3147370
    2265     459   659      11   216  6.54e-07     53.5     44.44  14.2401585
    2266     459   659      13   218  6.70e-07     53.5     46.01  14.2159881
    2267     459   705      12   265  7.30e-07     53.9     40.61  14.1302213
    1658     493   612     339   462  9.30e-07     54.7     50.78  13.8880813
    2268     559   649     166   253  1.06e-06     54.3     54.26  13.7572416
    2269     463   656      20   190  1.21e-06     52.8     43.96  13.6248902
    2270     231   280      14    67  1.85e-06     48.1     51.85  13.2003249
    2271     540   694     142   293  2.07e-06     52.4     49.04  13.0879620
    2272     437   656      52   270  2.08e-06     52.8     45.92  13.0831427
    2273     461   709      15   287  2.12e-06     52.8     44.01  13.0640945
    2274     452   598      12   161  2.22e-06     52.4     47.13  13.0180034
    2275     452   598      49   198  2.49e-06     52.4     47.13  12.9032278
    2276     526   669     151   287  3.99e-06     51.6     53.69  12.4317193
    2277     460   653      24   231  4.23e-06     52.0     45.37  12.3733086
    2278     559   665     158   260  4.93e-06     51.2     51.38  12.2201716
    2279     460   653      35   242  6.21e-06     51.2     45.37  11.9893497
    2280     459   693      11   254  7.10e-06     50.4     42.00  11.8554158
    2281     460   653      27   234  7.20e-06     50.8     45.37  11.8414295
    2282     460   653      26   233  7.64e-06     50.8     45.37  11.7821130
    2283     460   653      27   234  7.87e-06     50.8     45.37  11.7524525
    2284     460   653      26   233  7.98e-06     50.8     45.37  11.7385721
    2285     460   653      26   233  8.26e-06     50.8     45.37  11.7040860
    2286     452   658      16   222  9.04e-06     50.4     44.65  11.6138514
    2287     452   658      14   220  9.28e-06     50.1     44.65  11.5876490
    2288     463   656      14   241  1.04e-05     50.1     42.67  11.4737048
    2289     460   653      26   233  1.04e-05     50.4     45.37  11.4737048
    2290     460   653      27   234  1.08e-05     50.4     45.37  11.4359644
    2291     526   666     125   258  1.20e-05     50.1     52.05  11.3306039
    2292     460   653      33   240  1.20e-05     50.1     45.37  11.3306039
    2293     452   598      16   165  1.25e-05     50.1     45.91  11.2897819
    2294     526   666     131   264  1.35e-05     50.1     52.05  11.2128209
    2295     526   666     111   244  1.36e-05     49.7     52.05  11.2054408
    2296     526   666     125   258  1.38e-05     50.1     52.05  11.1908420
    2297     460   653      26   233  1.44e-05     49.7     45.37  11.1482824
    2298     556   650     266   359  1.47e-05     50.4     54.00  11.1276631
    2299     526   666     126   259  1.50e-05     49.7     52.05  11.1074604
    2300     526   666     125   258  1.63e-05     49.7     52.05  11.0243455
    2301     452   658      14   220  1.79e-05     49.3     45.37  10.9307098
    2302     526   666     116   249  1.82e-05     49.3     52.05  10.9140890
    2303     526   666     112   245  1.82e-05     49.3     52.05  10.9140890
    2304     452   598       3   152  1.91e-05     48.9     45.91  10.8658222
    2305     526   666     116   249  1.91e-05     49.3     52.05  10.8658222
    2306     448   658      42   247  2.02e-05     49.3     45.16  10.8098280
    2307     234   280       9    59  3.62e-05     43.9     50.98  10.2264514
    2308     537   666     144   263  3.67e-05     48.5     50.76  10.2127338
    2309     234   280       9    59  3.70e-05     43.9     50.98  10.2045926
    2310     452   658      22   228  3.82e-05     48.5     44.91  10.1726750
    2311     452   658      30   236  3.97e-05     48.1     44.91  10.1341594
    2312     537   666     124   243  3.98e-05     48.5     51.52  10.1316436
    2313     452   658      20   226  4.00e-05     48.1     44.91  10.1266311
    2314     537   666     124   243  4.05e-05     48.5     51.52  10.1142086
    2315     452   658      22   228  4.07e-05     48.5     44.91  10.1092825
    2316     452   658      14   220  4.19e-05     48.1     44.91  10.0802247
    2317     537   666     124   243  4.19e-05     48.5     51.52  10.0802247
    2318     537   666     124   243  4.19e-05     48.5     51.52  10.0802247
    2319     452   658      10   216  4.29e-05     48.1     44.91  10.0566387
    2320     452   658      16   222  4.37e-05     48.1     44.91  10.0381625
    2321     452   658      21   227  4.50e-05     48.1     44.91  10.0088481
    2322     234   280       9    59  4.68e-05     43.5     50.98   9.9696274
    2323     452   658      15   221  4.89e-05     48.1     44.91   9.9257332
    2324     452   658      16   222  4.99e-05     48.1     44.91   9.9054896
    2325     452   658      60   266  5.18e-05     48.1     44.91   9.8681204
    2326     234   280       9    59  5.42e-05     43.5     50.98   9.8228296
    2327     452   658      66   272  5.43e-05     48.1     44.91   9.8209863
    2328     235   280       1    50  5.95e-05     42.7     50.00   9.7295342
    2329     235   280       3    52  6.41e-05     42.7     50.00   9.6550662
    2330     453   670      39   254  6.46e-05     47.8     45.49   9.6472961
    2331     235   280      13    62  7.50e-05     43.1     50.00   9.4980224
    2332     234   280       9    59  8.35e-05     42.7     50.98   9.3906639
    2333     526   666     115   248  1.02e-04     47.0     50.68   9.1905377
    2334     537   666     124   243  1.03e-04     47.0     51.52   9.1807816
    2335     233   280       1    52  1.05e-04     45.8     51.92   9.1615502
    2336     461   666      43   243  1.06e-04     47.0     45.87   9.1520715
    2337     526   666     115   248  1.16e-04     47.0     50.68   9.0619204
    2338     559   666     133   237  1.40e-04     46.6     53.64   8.8738681
    2339     459   659      13   220  1.42e-04     46.2     41.59   8.8596835
    2340     459   659      14   221  1.47e-04     46.2     41.59   8.8250780
    2341     537   666     129   248  1.61e-04     46.6     51.52   8.7341062
    2342     537   666     124   243  1.66e-04     46.2     51.52   8.7035228
    2343     519   658     128   266  1.70e-04     46.6     47.26   8.6797121
    2344     537   666     129   248  1.71e-04     46.2     51.52   8.6738470
    2345     461   666      42   242  1.72e-04     46.2     45.41   8.6680161
    2346     461   666      42   242  1.78e-04     46.2     45.41   8.6337270
    2347     537   666     125   244  1.84e-04     46.2     51.52   8.6005748
    2348     461   666      41   241  1.86e-04     46.2     45.41   8.5897639
    2349     537   666     144   263  1.87e-04     46.2     51.52   8.5844019
    2350     532   666     119   243  1.90e-04     46.2     50.36   8.5684865
    2351     537   666     124   243  1.94e-04     46.2     51.52   8.5476524
    2352     537   666     124   243  1.95e-04     46.2     51.52   8.5425110
    2353     537   666     138   257  1.95e-04     46.2     51.52   8.5425110
    2354     537   666     124   243  1.98e-04     46.2     51.52   8.5272435
    2355     461   666      41   241  1.99e-04     46.2     45.41   8.5222057
    2356     234   295     212   273  2.03e-04     46.6     50.00   8.5023046
    2357     537   666     124   243  2.06e-04     46.2     51.52   8.4876344
    2358     537   666     123   242  2.15e-04     45.8     51.52   8.4448725
    2359     537   666     123   242  2.15e-04     45.8     51.52   8.4448725
    2360     537   666     124   243  2.16e-04     45.8     51.52   8.4402322
    2361     537   666     124   243  2.17e-04     45.8     51.52   8.4356132
    2362     537   666     146   265  2.23e-04     46.2     51.52   8.4083388
    2363     537   666     123   242  2.30e-04     45.8     51.52   8.3774312
    2364     537   666     124   243  2.31e-04     45.8     51.52   8.3730928
    2365     537   666     124   243  2.34e-04     45.8     51.52   8.3601894
    2366     537   666     123   242  2.38e-04     45.8     51.52   8.3432399
    2367     537   666     124   243  2.45e-04     45.8     51.52   8.3142523
    2368     537   666     124   243  2.48e-04     45.8     51.52   8.3020818
    2369     537   666     122   241  2.49e-04     45.8     51.52   8.2980577
    2370     537   666     122   241  2.57e-04     45.8     51.52   8.2664345
    2371     537   666     137   256  2.72e-04     45.8     51.52   8.2097085
    2372     231   280       7    60  2.89e-04     42.0     48.15   8.1490839
    2373     537   666     146   265  3.01e-04     45.4     51.52   8.1084003
    2374     537   666     147   266  3.02e-04     45.4     51.52   8.1050835
    2375     537   666     147   266  3.05e-04     45.4     51.52   8.0951988
    2376     537   666     129   248  3.19e-04     45.4     50.76   8.0503195
    2377     537   666     171   290  3.29e-04     45.8     51.52   8.0194528
    2378     532   666     119   243  3.34e-04     45.4     50.36   8.0043696
    2379     235   379     210   342  3.36e-04     45.8     39.33   7.9983994
    2380     459   659      12   219  3.37e-04     45.4     41.48   7.9954276
    2381     459   659      10   217  3.69e-04     45.1     42.73   7.9047139
    2382     559   666     139   243  5.10e-04     44.7     53.64   7.5810998
    2383     451   653      50   263  7.05e-04     44.3     45.25   7.2573128
    2384     234   277       4    46  7.49e-04     39.7     56.82   7.1967716
    1375     450   570     580   708  1.00e-03     44.3     49.24   6.9077553
    1380     450   570     576   704  1.00e-03     43.9     49.24   6.9077553
    1382     450   570     578   706  1.00e-03     43.9     49.24   6.9077553
    2385     458   655      33   221  1.00e-03     43.5     39.81   6.9077553
    2386     458   655      33   221  2.00e-03     42.7     39.34   6.2146081
    2387     559   699     154   294  2.00e-03     43.1     50.00   6.2146081
    2388     496   655     152   300  2.00e-03     43.1     42.33   6.2146081
    2389     458   655      21   209  2.00e-03     42.4     38.94   6.2146081
    2390     559   699     135   275  3.00e-03     42.4     50.00   5.8091430
    2391     559   666     134   238  3.00e-03     42.4     52.73   5.8091430
    2392     559   699     133   273  3.00e-03     42.4     50.00   5.8091430
    2393     559   666     133   237  3.00e-03     42.4     52.73   5.8091430
    2394     559   666     133   237  3.00e-03     42.4     52.73   5.8091430
    2395     559   699     134   274  4.00e-03     42.0     50.00   5.5214609
    2396     559   699     133   273  4.00e-03     42.0     50.00   5.5214609
    2397     559   699     133   273  4.00e-03     42.0     50.00   5.5214609
    2398     559   699     134   274  4.00e-03     42.0     50.00   5.5214609
    2399     559   699     133   273  4.00e-03     42.0     50.00   5.5214609
    781      232   285      35    92  6.00e-03     42.0     44.83   5.1159958
    2400     235   272      14    55  6.00e-03     37.7     52.38   5.1159958
    771      232   285      35    92  1.20e-02     40.8     44.83   4.4228486
    773      232   285      35    92  1.20e-02     40.8     44.83   4.4228486
    2401     240   273     106   143  1.40e-02     38.9     60.53   4.2686979
    2402     452   658      16   201  1.40e-02     40.0     41.40   4.2686979
    1934     450   548      60   167  1.50e-02     40.4     48.62   4.1997051
    2403     462   598      38   186  2.50e-02     39.7     45.86   3.6888795
    2404     218   280       6    67  4.80e-02     35.4     47.76   3.0365543
    2405     233   281      45    97  4.80e-02     39.3     49.06   3.0365543
    2406     233   281      37    89  6.40e-02     38.5     49.06   2.7488722
    2407     463   528     216   278  7.30e-02     37.7     55.22   2.6172958
    2408     216   273     151   209  8.40e-02     37.4     47.46   2.4769385
    2409     458   670      25   251  8.70e-02     37.7     41.60   2.4418472
    2410     458   670      25   251  9.30e-02     37.7     41.60   2.3751558
    2411     235   280      24    73  1.10e-01     34.3     46.00   2.2072749
    2412     235   280       8    57  5.60e-01     32.3     44.00   0.5798185
    2413     235   278      18    65  5.80e-01     32.3     43.75   0.5447272
    2414     257   280       4    27  8.50e-01     33.9     58.33   0.1625189
    2415     456   655      22   242  8.50e-01     34.7     38.56   0.1625189
    2416     463   600      28   194  3.60e+00     32.7     37.50  -1.2809338
    1773     235   280     141   190  4.00e+00     32.7     44.00  -1.3862944
    2417     235   281       9    58  5.90e+00     29.3     42.00  -1.7749524
         pdb.id      acc
    1    6Q0K_A   6Q0K_A
    2    7MFD_A   7MFD_A
    3    6Q0J_A   6Q0J_A
    4    6NYB_A   6NYB_A
    5    8VYO_A   8VYO_A
    6    7ZR0_K   7ZR0_K
    7    6UAN_B   6UAN_B
    8    8VYV_C   8VYV_C
    9    8VYU_A   8VYU_A
    10   8VYP_C   8VYP_C
    11   8VYS_A   8VYS_A
    12   8VYW_C   8VYW_C
    13   8VYQ_C   8VYQ_C
    14   9MMS_A   9MMS_A
    15   9MMR_A   9MMR_A
    16   9MMQ_A   9MMQ_A
    17   8CHF_A   8CHF_A
    18   8U1L_C   8U1L_C
    19   9MMP_A   9MMP_A
    20   7Z37_c 7Z37_CP1
    21   4MNF_A   4MNF_A
    22   4MNE_B   4MNE_B
    23   3II5_A   3II5_A
    24   3D4Q_A   3D4Q_A
    25   7SHV_A   7SHV_A
    26   5FD2_A   5FD2_A
    27   4EHG_A   4EHG_A
    28   3IDP_A   3IDP_A
    29   4MBJ_A   4MBJ_A
    30   5HIE_A   5HIE_A
    31   4DBN_A   4DBN_A
    32   3Q96_A   3Q96_A
    33   6PP9_A   6PP9_A
    34   9ECU_A   9ECU_A
    35   2FB8_A   2FB8_A
    36   9AXX_B   9AXX_B
    37   4H58_A   4H58_A
    38   1UWH_A   1UWH_A
    39   1UWJ_A   1UWJ_A
    40   9BFB_A   9BFB_A
    41   6N0P_A   6N0P_A
    42   6U2H_C   6U2H_C
    43   4YHT_A   4YHT_A
    44   4XV9_A   4XV9_A
    45   4JVG_A   4JVG_A
    46   5CSW_A   5CSW_A
    47   3C4C_A   3C4C_A
    48   4RZV_A   4RZV_A
    49   4WO5_A   4WO5_A
    50   4FK3_A   4FK3_A
    51   8F7O_A   8F7O_A
    52   6XFP_A   6XFP_A
    53   8C7Y_A   8C7Y_A
    54   4XV1_A   4XV1_A
    55   6V34_A   6V34_A
    56   3OG7_A   3OG7_A
    57   5ITA_A   5ITA_A
    58   5JRQ_A   5JRQ_A
    59   6P3D_A   6P3D_A
    60   4CQE_A   4CQE_A
    61   7P3V_A   7P3V_A
    62   5HI2_A   5HI2_A
    63   9AXC_A   9AXC_A
    64   9AXA_A   9AXA_A
    65   3OMV_A   3OMV_A
    66   9O0U_A   9O0U_A
    67   9AY7_A   9AY7_A
    68   8GFT_D   8GFT_D
    69   9AXM_B   9AXM_B
    70   8GAE_D   8GAE_D
    71   6PTS_D   6PTS_D
    72   7JHP_C   7JHP_C
    73   6XI7_B   6XI7_B
    74   2Y4I_B   2Y4I_B
    75   5KKR_B   5KKR_B
    76   8BW9_D   8BW9_D
    77   2L05_A   2L05_A
    78   7JUQ_B   7JUQ_B
    79   3NY5_A   3NY5_A
    80   6XGU_B   6XGU_B
    81   5J17_A   5J17_A
    82   7JUW_B   7JUW_B
    83   5VR3_A   5VR3_A
    84   9AXH_C   9AXH_C
    85   3PPZ_A   3PPZ_A
    86   4CSV_A   4CSV_A
    87   3P86_A   3P86_A
    88   5VYK_A   5VYK_A
    89   8DEG_A   8DEG_A
    90   5CEN_A   5CEN_A
    91   9NS1_A   9NS1_A
    92   8JF3_A   8JF3_A
    93   2BDF_A   2BDF_A
    94   7NG7_A   7NG7_A
    95   7OTE_A   7OTE_A
    96   4MXO_A   4MXO_A
    97   6E6E_A   6E6E_A
    98   8HAQ_A   8HAQ_A
    99   3A4O_X   3A4O_X
    100  4MXX_A   4MXX_A
    101  1YOL_A   1YOL_A
    102  9NS0_A   9NS0_A
    103  3OEZ_A   3OEZ_A
    104  2OIQ_A   2OIQ_A
    105  1YI6_A   1YI6_A
    106  2OG8_A   2OG8_A
    107  1YOJ_A   1YOJ_A
    108  3D7U_B   3D7U_B
    109  4MXY_A   4MXY_A
    110  2PL0_A   2PL0_A
    111  2OFV_A   2OFV_A
    112  3BYS_A   3BYS_A
    113  5XY1_A   5XY1_A
    114  2OF2_A   2OF2_A
    115  2DQ7_X   2DQ7_X
    116  3GEQ_A   3GEQ_A
    117  3KXZ_A   3KXZ_A
    118  5T0P_A   5T0P_A
    119  5SWH_A   5SWH_A
    120  2HWO_A   2HWO_A
    121  2ZM1_A   2ZM1_A
    122  1QPC_A   1QPC_A
    123  3SVV_A   3SVV_A
    124  2OFU_A   2OFU_A
    125  3KMM_A   3KMM_A
    126  6PDJ_A   6PDJ_A
    127  3U4W_A   3U4W_A
    128  4MCV_A   4MCV_A
    129  3G6H_A   3G6H_A
    130  3BYO_A   3BYO_A
    131  3BYM_A   3BYM_A
    132  4LGH_A   4LGH_A
    133  2H8H_A   2H8H_A
    134  8JN8_A   8JN8_A
    135  9IRL_A   9IRL_A
    136  2ZV7_A   2ZV7_A
    137  1FMK_A   1FMK_A
    138  2HK5_A   2HK5_A
    139  5ZJ6_A   5ZJ6_A
    140  1Y57_A   1Y57_A
    141  2QI8_A   2QI8_A
    142  9BT8_C   9BT8_C
    143  3DQW_A   3DQW_A
    144  3MPM_A   3MPM_A
    145  8XN8_A   8XN8_A
    146  6F3F_A   6F3F_A
    147  4K11_A   4K11_A
    148  1KSW_A   1KSW_A
    149  1QPD_A   1QPD_A
    150  2PTK_A   2PTK_A
    151  6IN0_A   6IN0_A
    152  2QOK_A   2QOK_A
    153  2QOC_A   2QOC_A
    154  3DZQ_A   3DZQ_A
    155  1AD5_A   1AD5_A
    156  2QOO_A   2QOO_A
    157  2QOI_A   2QOI_A
    158  2QOF_A   2QOF_A
    159  2QOD_A   2QOD_A
    160  3FXX_A   3FXX_A
    161  2GSF_A   2GSF_A
    162  1QCF_A   1QCF_A
    163  2QOL_A   2QOL_A
    164  9BYJ_A   9BYJ_A
    165  2YN8_A   2YN8_A
    166  2QOB_A   2QOB_A
    167  6FNI_A   6FNI_A
    168  2QON_A   2QON_A
    169  2QO7_A   2QO7_A
    170  6CZ2_A   6CZ2_A
    171  2XYU_A   2XYU_A
    172  2Y6M_A   2Y6M_A
    173  6CZ3_A   6CZ3_A
    174  2HEL_A   2HEL_A
    175  3KUL_A   3KUL_A
    176  2VWU_A   2VWU_A
    177  8JNA_B   8JNA_B
    178  3KUL_B   3KUL_B
    179  3ZEW_A   3ZEW_A
    180  4AW5_A   4AW5_A
    181  5MJA_A   5MJA_A
    182  5MJB_A   5MJB_A
    183  7UY0_B   7UY0_B
    184  7UY0_A   7UY0_A
    185  7KPL_A   7KPL_A
    186  4LGG_A   4LGG_A
    187  3ZFX_A   3ZFX_A
    188  6UMW_A   6UMW_A
    189  2R2P_A   2R2P_A
    190  3ZFY_A   3ZFY_A
    191  3ZFM_A   3ZFM_A
    192  1JPA_A   1JPA_A
    193  5D7V_A   5D7V_A
    194  5DA3_A   5DA3_A
    195  1K3A_A   1K3A_A
    196  5H2U_A   5H2U_A
    197  2ZM3_A   2ZM3_A
    198  5I9U_A   5I9U_A
    199  8XPV_A   8XPV_A
    200  5FXS_A   5FXS_A
    201  3O23_A   3O23_A
    202  1JQH_A   1JQH_A
    203  3QQU_A   3QQU_A
    204  1MQB_A   1MQB_A
    205  1M7N_A   1M7N_A
    206  2OJ9_A   2OJ9_A
    207  4TRL_A   4TRL_A
    208  4D2R_A   4D2R_A
    209  5FXQ_A   5FXQ_A
    210  3I81_A   3I81_A
    211  8PYI_a 8PYI_AAA
    212  5FXR_A   5FXR_A
    213  3D94_A   3D94_A
    214  4P2K_A   4P2K_A
    215  5EK7_A   5EK7_A
    216  8XKP_A   8XKP_A
    217  7KJA_A   7KJA_A
    218  6JMF_A   6JMF_A
    219  7KJC_A   7KJC_A
    220  7KJB_A   7KJB_A
    221  3LVP_A   3LVP_A
    222  3BKB_A   3BKB_A
    223  7EEF_A   7EEF_A
    224  5X5O_A   5X5O_A
    225  3LW0_A   3LW0_A
    226  5HES_A   5HES_A
    227  1P4O_A   1P4O_A
    228  7EEC_A   7EEC_A
    229  2REI_A   2REI_A
    230  7EED_A   7EED_A
    231  4XLV_A   4XLV_A
    232  1RQQ_A   1RQQ_A
    233  3CD3_A   3CD3_A
    234  8ATL_b 8ATL_BBB
    235  4UYA_A   4UYA_A
    236  8ATL_a 8ATL_AAA
    237  2Z8C_A   2Z8C_A
    238  8ATB_a 8ATB_AAA
    239  8FLN_A   8FLN_A
    240  1GAG_A   1GAG_A
    241  7TYJ_A   7TYJ_A
    242  6J6M_A   6J6M_A
    243  6MNY_A   6MNY_A
    244  3DK3_A   3DK3_A
    245  1I44_A   1I44_A
    246  6VXQ_A   6VXQ_A
    247  4IBM_A   4IBM_A
    248  6BL8_A   6BL8_A
    249  8DWN_A   8DWN_A
    250  5HU9_A   5HU9_A
    251  4Z3V_A   4Z3V_A
    252  6NZM_A   6NZM_A
    253  4XEY_A   4XEY_A
    254  6W7O_A   6W7O_A
    255  4YHF_A   4YHF_A
    256  3GEN_A   3GEN_A
    257  8YVV_A   8YVV_A
    258  5J87_A   5J87_A
    259  3P08_A   3P08_A
    260  6TFP_A   6TFP_A
    261  5P9F_A   5P9F_A
    262  4OT5_A   4OT5_A
    263  3OCT_A   3OCT_A
    264  6S90_A   6S90_A
    265  3PIX_A   3PIX_A
    266  4WA9_A   4WA9_A
    267  1IRK_A   1IRK_A
    268  6AUB_A   6AUB_A
    269  6O8I_A   6O8I_A
    270  5XYZ_A   5XYZ_A
    271  6NFI_A   6NFI_A
    272  7KXQ_A   7KXQ_A
    273  8GC8_A   8GC8_A
    274  5ZZ4_A   5ZZ4_A
    275  2XYN_A   2XYN_A
    276  1FPU_A   1FPU_A
    277  5FBN_C   5FBN_C
    278  4XLI_A   4XLI_A
    279  4ZLY_A   4ZLY_A
    280  8FLL_A   8FLL_A
    281  9ZLJ_A   9ZLJ_A
    282  6XRG_A   6XRG_A
    283  2F4J_A   2F4J_A
    284  3OCS_A   3OCS_A
    285  3QRI_A   3QRI_A
    286  8SSN_A   8SSN_A
    287  2HZI_A   2HZI_A
    288  2NRY_A   2NRY_A
    289  4RX5_A   4RX5_A
    290  2E2B_A   2E2B_A
    291  4OTF_A   4OTF_A
    292  3OXZ_A   3OXZ_A
    293  6O8U_A   6O8U_A
    294  2QOH_A   2QOH_A
    295  5K72_A   5K72_A
    296  6EG9_A   6EG9_A
    297  6XE4_A   6XE4_A
    298  6XR6_A   6XR6_A
    299  4ZOG_A   4ZOG_A
    300  4Y95_A   4Y95_A
    301  2NRU_A   2NRU_A
    302  6THW_A   6THW_A
    303  2HIW_A   2HIW_A
    304  7N9G_A   7N9G_A
    305  6THX_A   6THX_A
    306  6E4F_A   6E4F_A
    307  2HZ0_A   2HZ0_A
    308  9PSU_A   9PSU_A
    309  6LXY_A   6LXY_A
    310  2HYY_A   2HYY_A
    311  2G1T_A   2G1T_A
    312  6BIK_A   6BIK_A
    313  1P14_A   1P14_A
    314  2G2F_A   2G2F_A
    315  3EKK_A   3EKK_A
    316  8SCV_A   8SCV_A
    317  7C2W_A   7C2W_A
    318  8X2A_A   8X2A_A
    319  8H7F_A   8H7F_A
    320  3DK6_A   3DK6_A
    321  6F3D_A   6F3D_A
    322  9RPV_n  9RPV_N1
    323  3SXR_A   3SXR_A
    324  6F3I_A   6F3I_A
    325  6MOM_A   6MOM_A
    326  4RMZ_A   4RMZ_A
    327  2OIB_A   2OIB_A
    328  5UIT_A   5UIT_A
    329  3PYY_A   3PYY_A
    330  4Y73_A   4Y73_A
    331  7CC2_A   7CC2_A
    332  7C2V_A   7C2V_A
    333  7QG3_A   7QG3_A
    334  9NA2_A   9NA2_A
    335  5W84_A   5W84_A
    336  8BR5_a 8BR5_AAA
    337  6N8G_A   6N8G_A
    338  6I99_A   6I99_A
    339  6F3E_A   6F3E_A
    340  6JK8_A   6JK8_A
    341  5E1S_A   5E1S_A
    342  1OPK_A   1OPK_A
    343  5HHW_A   5HHW_A
    344  6F3G_A   6F3G_A
    345  4Y93_A   4Y93_A
    346  3K54_A   3K54_A
    347  5UIS_A   5UIS_A
    348  5UIQ_A   5UIQ_A
    349  2GQG_A   2GQG_A
    350  6O94_A   6O94_A
    351  8V1O_A   8V1O_A
    352  7SL1_A   7SL1_A
    353  7W7Y_A   7W7Y_A
    354  1OPL_A   1OPL_A
    355  7W7X_A   7W7X_A
    356  8EYR_A   8EYR_A
    357  8DTL_A   8DTL_A
    358  3QRJ_A   3QRJ_A
    359  2FO0_A   2FO0_A
    360  4YFF_A   4YFF_A
    361  4U97_A   4U97_A
    362  4YFI_A   4YFI_A
    363  6AUA_A   6AUA_A
    364  7MGJ_A   7MGJ_A
    365  6EGD_A   6EGD_A
    366  8S9F_A   8S9F_A
    367  8WTF_A   8WTF_A
    368  8EYX_A   8EYX_A
    369  4XI2_A   4XI2_A
    370  2Z60_A   2Z60_A
    371  3OY3_A   3OY3_A
    372  8DKS_A   8DKS_A
    373  4TWP_A   4TWP_A
    374  5MO4_A   5MO4_A
    375  6PYH_A   6PYH_A
    376  5BPY_A   5BPY_A
    377  3DK7_A   3DK7_A
    378  4NWM_A   4NWM_A
    379  2V7A_A   2V7A_A
    380  5AMN_A   5AMN_A
    381  7KHJ_A   7KHJ_A
    382  7KHK_A   7KHK_A
    383  8FD9_A   8FD9_A
    384  1K2P_A   1K2P_A
    385  8GMB_A   8GMB_A
    386  6FEK_A   6FEK_A
    387  9KS5_A   9KS5_A
    388  7BW7_A   7BW7_A
    389  8U4B_A   8U4B_A
    390  7PG0_A   7PG0_A
    391  8VJB_A   8VJB_A
    392  2IVV_A   2IVV_A
    393  3ETA_A   3ETA_A
    394  8S93_A   8S93_A
    395  6NE7_A   6NE7_A
    396  4HVS_A   4HVS_A
    397  6PXV_A   6PXV_A
    398  6NJA_A   6NJA_A
    399  4GU9_A   4GU9_A
    400  9SDI_A   9SDI_A
    401  2IVT_A   2IVT_A
    402  3PXK_A   3PXK_A
    403  6I8Z_A   6I8Z_A
    404  2ETM_A   2ETM_A
    405  2IVS_A   2IVS_A
    406  2O8Y_A   2O8Y_A
    407  6I83_A   6I83_A
    408  1MP8_A   1MP8_A
    409  2R4B_A   2R4B_A
    410  2J0M_B   2J0M_B
    411  2JKM_A   2JKM_A
    412  3BBT_B   3BBT_B
    413  4EBV_A   4EBV_A
    414  6GQJ_A   6GQJ_A
    415  8PQ9_A   8PQ9_A
    416  1U54_A   1U54_A
    417  3ZZW_A   3ZZW_A
    418  4EWH_A   4EWH_A
    419  6GQK_A   6GQK_A
    420  6VQM_A   6VQM_A
    421  1U46_A   1U46_A
    422  4XCU_A   4XCU_A
    423  7KP6_A   7KP6_A
    424  3BZ3_A   3BZ3_A
    425  4GT4_A   4GT4_A
    426  5ZXB_A   5ZXB_A
    427  2JKK_A   2JKK_A
    428  8PQG_A   8PQG_A
    429  8FE9_A   8FE9_A
    430  3EQP_A   3EQP_A
    431  2J0L_A   2J0L_A
    432  4ID7_A   4ID7_A
    433  4HZR_A   4HZR_A
    434  5FM2_A   5FM2_A
    435  4CKI_A   4CKI_A
    436  3QUP_A   3QUP_A
    437  6SDC_A   6SDC_A
    438  8ANS_A   8ANS_A
    439  9EF1_A   9EF1_A
    440  4HZS_A   4HZS_A
    441  7RUN_A   7RUN_A
    442  8AU5_A   8AU5_A
    443  6VHG_A   6VHG_A
    444  2J0J_A   2J0J_A
    445  6JPE_A   6JPE_A
    446  4UXQ_A   4UXQ_A
    447  8KH9_A   8KH9_A
    448  7DTZ_A   7DTZ_A
    449  5JKG_A   5JKG_A
    450  8W5C_A   8W5C_A
    451  7F3M_A   7F3M_A
    452  6YI8_A   6YI8_A
    453  4QQT_A   4QQT_A
    454  5NUD_A   5NUD_A
    455  5XFJ_A   5XFJ_A
    456  3D7T_A   3D7T_A
    457  7VJL_A   7VJL_A
    458  4TYE_A   4TYE_A
    459  6I82_A   6I82_A
    460  3T9T_A   3T9T_A
    461  1BYG_A   1BYG_A
    462  7WCW_A   7WCW_A
    463  5XFF_A   5XFF_A
    464  7WCX_A   7WCX_A
    465  6NFY_A   6NFY_A
    466  3D7U_A   3D7U_A
    467  4KNB_A   4KNB_A
    468  8AW1_A   8AW1_A
    469  4R1V_A   4R1V_A
    470  7PI4_d 7PI4_DDD
    471  4QQC_A   4QQC_A
    472  6SD9_A   6SD9_A
    473  8AN8_A   8AN8_A
    474  8VI1_A   8VI1_A
    475  6TY3_A   6TY3_A
    476  4QQJ_A   4QQJ_A
    477  2WD1_A   2WD1_A
    478  2J0K_A   2J0K_A
    479  3GQI_A   3GQI_A
    480  7B3Q_A   7B3Q_A
    481  6IUO_A   6IUO_A
    482  3F66_A   3F66_A
    483  2G15_A   2G15_A
    484  4EEV_A   4EEV_A
    485  3Q6U_A   3Q6U_A
    486  8PAS_A   8PAS_A
    487  3LQ8_A   3LQ8_A
    488  8FH4_A   8FH4_A
    489  6CQD_A   6CQD_A
    490  9C1R_A   9C1R_A
    491  2WGJ_A   2WGJ_A
    492  8CDW_A   8CDW_A
    493  4QQ5_A   4QQ5_A
    494  8PAU_A   8PAU_A
    495  6NG0_A   6NG0_A
    496  3VW8_A   3VW8_A
    497  8PAR_A   8PAR_A
    498  5A46_A   5A46_A
    499  7M0L_A   7M0L_A
    500  9SZJ_A   9SZJ_A
    501  4GG5_A   4GG5_A
    502  5UAB_A   5UAB_A
    503  5ZV2_A   5ZV2_A
    504  3I5N_A   3I5N_A
    505  7M0M_A   7M0M_A
    506  8YKI_A   8YKI_A
    507  5A4C_A   5A4C_A
    508  3C4F_A   3C4F_A
    509  7M0K_A   7M0K_A
    510  4ZSA_A   4ZSA_A
    511  3RHX_A   3RHX_A
    512  2RFN_A   2RFN_A
    513  9H8D_A   9H8D_A
    514  8XN7_A   8XN7_A
    515  4F63_A   4F63_A
    516  4WUN_A   4WUN_A
    517  5FLF_A   5FLF_A
    518  5HLW_A   5HLW_A
    519  7L25_A   7L25_A
    520  1AGW_A   1AGW_A
    521  6LVM_A   6LVM_A
    522  4F0F_A   4F0F_A
    523  3KXX_A   3KXX_A
    524  4PWN_A   4PWN_A
    525  5TQ5_A   5TQ5_A
    526  3JS2_A   3JS2_A
    527  4HCT_A   4HCT_A
    528  7WCL_A   7WCL_A
    529  6CQF_A   6CQF_A
    530  5AM7_A   5AM7_A
    531  8Y22_A   8Y22_A
    532  9BI8_A   9BI8_A
    533  3GQL_A   3GQL_A
    534  6MZW_A   6MZW_A
    535  7L24_A   7L24_A
    536  6HH1_A   6HH1_A
    537  3DKG_A   3DKG_A
    538  6NVL_A   6NVL_A
    539  8G6Z_A   8G6Z_A
    540  3Q6W_A   3Q6W_A
    541  6CQE_A   6CQE_A
    542  4RWI_A   4RWI_A
    543  7L26_A   7L26_A
    544  6NFZ_A   6NFZ_A
    545  7LL5_A   7LL5_A
    546  5J5T_A   5J5T_A
    547  7KAC_A   7KAC_A
    548  7LL4_A   7LL4_A
    549  3KRR_A   3KRR_A
    550  7RN6_A   7RN6_A
    551  5TQ4_A   5TQ4_A
    552  8AU3_A   8AU3_A
    553  3QTI_A   3QTI_A
    554  3CE3_A   3CE3_A
    555  3C1X_A   3C1X_A
    556  4IWD_A   4IWD_A
    557  5GRN_A   5GRN_A
    558  3A4P_A   3A4P_A
    559  2P2H_A   2P2H_A
    560  8UDT_A   8UDT_A
    561  5WEV_A   5WEV_A
    562  3DKC_A   3DKC_A
    563  8PQJ_A   8PQJ_A
    564  1R0P_A   1R0P_A
    565  4E4M_A   4E4M_A
    566  5TQ6_A   5TQ6_A
    567  3MIY_A   3MIY_A
    568  5VND_A   5VND_A
    569  5TQ3_A   5TQ3_A
    570  5USY_A   5USY_A
    571  9GB9_A   9GB9_A
    572  6VGL_A   6VGL_A
    573  6TU9_A   6TU9_A
    574  8UDV_A   8UDV_A
    575  4HGE_A   4HGE_A
    576  3UGC_A   3UGC_A
    577  8G8O_A   8G8O_A
    578  1K9A_A   1K9A_A
    579  3TJC_A   3TJC_A
    580  6AAJ_A   6AAJ_A
    581  3C7Q_A   3C7Q_A
    582  3TT0_A   3TT0_A
    583  7MO7_B   7MO7_B
    584  2OH4_A   2OH4_A
    585  9V70_A   9V70_A
    586  8BX6_A   8BX6_A
    587  3E62_A   3E62_A
    588  2XA4_A   2XA4_A
    589  4D0W_A   4D0W_A
    590  7Q7I_A   7Q7I_A
    591  2B7A_A   2B7A_A
    592  4AQC_A   4AQC_A
    593  8BM2_A   8BM2_A
    594  3EWH_A   3EWH_A
    595  3U6J_A   3U6J_A
    596  3ZMM_A   3ZMM_A
    597  9V71_A   9V71_A
    598  4BBE_A   4BBE_A
    599  5HEZ_A   5HEZ_A
    600  3CJF_A   3CJF_A
    601  1T46_A   1T46_A
    602  4F1O_A   4F1O_A
    603  3RVG_A   3RVG_A
    604  6WTN_A   6WTN_A
    605  2W1I_A   2W1I_A
    606  7UYW_A   7UYW_A
    607  2OGV_A   2OGV_A
    608  2I0V_A   2I0V_A
    609  1YWN_A   1YWN_A
    610  4U0I_A   4U0I_A
    611  6MOB_A   6MOB_A
    612  6WXJ_A   6WXJ_A
    613  6TPD_A   6TPD_A
    614  3Q32_A   3Q32_A
    615  4YTC_A   4YTC_A
    616  4F1M_A   4F1M_A
    617  9D3F_A   9D3F_A
    618  2X7F_A   2X7F_A
    619  8W1L_A   8W1L_A
    620  7TNH_A   7TNH_A
    621  3JY9_A   3JY9_A
    622  5HOA_A   5HOA_A
    623  3CJG_A   3CJG_A
    624  7SIU_A   7SIU_A
    625  5UGL_A   5UGL_A
    626  3IO7_A   3IO7_A
    627  7XZQ_A   7XZQ_A
    628  7ZVS_A   7ZVS_A
    629  7FCZ_A   7FCZ_A
    630  5AX9_A   5AX9_A
    631  6RA5_A   6RA5_A
    632  8ZML_A   8ZML_A
    633  3V5J_A   3V5J_A
    634  6C4D_A   6C4D_A
    635  9HY8_A   9HY8_A
    636  1VR2_A   1VR2_A
    637  7UOS_A   7UOS_A
    638  1SM2_A   1SM2_A
    639  4ZIM_A   4ZIM_A
    640  6NYH_A   6NYH_A
    641  9WYR_A   9WYR_A
    642  9KFU_A   9KFU_A
    643  6OL2_A   6OL2_A
    644  8X88_A   8X88_A
    645  4ITH_A   4ITH_A
    646  5DRB_A   5DRB_A
    647  2P2I_A   2P2I_A
    648  3VNT_A   3VNT_A
    649  3PLS_A   3PLS_A
    650  9CD7_A   9CD7_A
    651  8W3D_A   8W3D_A
    652  9MZZ_A   9MZZ_A
    653  8W3B_A   8W3B_A
    654  9MZY_A   9MZY_A
    655  9MZX_A   9MZX_A
    656  8W2X_A   8W2X_A
    657  8JQI_B   8JQI_B
    658  4K33_A   4K33_A
    659  1PKG_A   1PKG_A
    660  8W38_A   8W38_A
    661  4Q2A_A   4Q2A_A
    662  7TEU_A   7TEU_A
    663  2XIR_A   2XIR_A
    664  3LCD_A   3LCD_A
    665  9GZH_A   9GZH_A
    666  8PQH_A   8PQH_A
    667  5W7T_A   5W7T_A
    668  5HX6_A   5HX6_A
    669  4APC_A   4APC_A
    670  7ZW8_A   7ZW8_A
    671  5HOR_A   5HOR_A
    672  6ITV_A   6ITV_A
    673  6A32_A   6A32_A
    674  6C3E_A   6C3E_A
    675  6GQO_A   6GQO_A
    676  4E6D_A   4E6D_A
    677  3WZD_A   3WZD_A
    678  6PNX_A   6PNX_A
    679  6NW2_A   6NW2_A
    680  5UGX_A   5UGX_A
    681  3G0E_A   3G0E_A
    682  9GTG_A   9GTG_A
    683  1T45_A   1T45_A
    684  6JOI_A   6JOI_A
    685  3G0F_A   3G0F_A
    686  4J98_A   4J98_A
    687  1FAQ_A   1FAQ_A
    688  2PVF_A   2PVF_A
    689  6FD3_A   6FD3_A
    690  6AGX_A   6AGX_A
    691  8XRR_A   8XRR_A
    692  3CLY_A   3CLY_A
    693  5TF9_A   5TF9_A
    694  4USF_A   4USF_A
    695  4AGC_A   4AGC_A
    696  3BEA_A   3BEA_A
    697  8BEM_A   8BEM_A
    698  4GL9_A   4GL9_A
    699  2I1M_A   2I1M_A
    700  2PZ5_A   2PZ5_A
    701  8E1X_A   8E1X_A
    702  6BDN_A   6BDN_A
    703  2J51_A   2J51_A
    704  2JFM_A   2JFM_A
    705  2JFL_A   2JFL_A
    706  2PZP_A   2PZP_A
    707  4J97_A   4J97_A
    708  4KIO_A   4KIO_A
    709  8H75_A   8H75_A
    710  4NEU_A   4NEU_A
    711  2PWL_A   2PWL_A
    712  6LVK_A   6LVK_A
    713  9U7E_A   9U7E_A
    714  8SWE_A   8SWE_A
    715  1GJO_A   1GJO_A
    716  3RI1_A   3RI1_A
    717  5EG3_A   5EG3_A
    718  3B2T_A   3B2T_A
    719  8STG_A   8STG_A
    720  2PVY_A   2PVY_A
    721  6LVL_A   6LVL_A
    722  9GFZ_A   9GFZ_A
    723  3QGW_A   3QGW_A
    724  9LBG_A   9LBG_A
    725  7OZY_a 7OZY_AAA
    726  8E4T_A   8E4T_A
    727  4J96_A   4J96_A
    728  4J99_A   4J99_A
    729  2PSQ_A   2PSQ_A
    730  9U3N_A   9U3N_A
    731  5UHN_A   5UHN_A
    732  7KIA_A   7KIA_A
    733  2Q0B_A   2Q0B_A
    734  8JOT_A   8JOT_A
    735  2PZR_A   2PZR_A
    736  9LBF_A   9LBF_A
    737  2PY3_A   2PY3_A
    738  3LCO_A   3LCO_A
    739  5UI0_A   5UI0_A
    740  6T2W_A   6T2W_A
    741  9D51_A   9D51_A
    742  4HW7_A   4HW7_A
    743  7AAY_A   7AAY_A
    744  8CGC_A   8CGC_A
    745  4O27_B   4O27_B
    746  9BHI_A   9BHI_A
    747  5U6B_A   5U6B_A
    748  9D7Q_A   9D7Q_A
    749  9M41_A   9M41_A
    750  7DXL_A   7DXL_A
    751  6V6Q_A   6V6Q_A
    752  2GCD_A   2GCD_A
    753  3S95_A   3S95_A
    754  1U5Q_A   1U5Q_A
    755  4QML_A   4QML_A
    756  7AAZ_A   7AAZ_A
    757  7OAM_A   7OAM_A
    758  3A7F_A   3A7F_A
    759  7CQE_A   7CQE_A
    760  4W8E_A   4W8E_A
    761  5TC0_A   5TC0_A
    762  4U8Z_A   4U8Z_A
    763  7B30_A   7B30_A
    764  2P0C_A   2P0C_A
    765  5O1V_A   5O1V_A
    766  7Z5W_A   7Z5W_A
    767  5U6C_A   5U6C_A
    768  3ZHP_C   3ZHP_C
    769  5O2B_A   5O2B_A
    770  8SE1_A   8SE1_A
    772  8SE2_A   8SE2_A
    774  2XIK_A   2XIK_A
    775  3CKW_A   3CKW_A
    776  1LUF_A   1LUF_A
    777  3CKX_A   3CKX_A
    778  7Z4V_A   7Z4V_A
    779  4V0G_B   4V0G_B
    780  3PFQ_A   3PFQ_A
    782  4RIO_A   4RIO_A
    783  5CNN_A   5CNN_A
    784  8HV2_A   8HV2_A
    785  5HVJ_A   5HVJ_A
    786  3LXN_A   3LXN_A
    787  4ZJV_A   4ZJV_A
    788  5HVK_A   5HVK_A
    789  4GVJ_A   4GVJ_A
    790  8EDH_A   8EDH_A
    791  7SYD_A   7SYD_A
    792  5O26_A   5O26_A
    793  4V0G_A   4V0G_A
    794  5CAV_A   5CAV_A
    795  1M14_A   1M14_A
    796  7UKV_A   7UKV_A
    797  1YVJ_A   1YVJ_A
    798  3D5V_A   3D5V_A
    799  4HVD_A   4HVD_A
    800  5O2C_A   5O2C_A
    801  4TKS_A   4TKS_A
    802  4Z16_A   4Z16_A
    803  4I23_A   4I23_A
    804  3D5U_A   3D5U_A
    805  3VJO_A   3VJO_A
    806  2RFD_A   2RFD_A
    807  2ITW_A   2ITW_A
    808  9FZR_A   9FZR_A
    809  3D5W_A   3D5W_A
    810  4G5J_A   4G5J_A
    811  4JQ7_A   4JQ7_A
    812  4LI5_A   4LI5_A
    813  4RIW_B   4RIW_B
    814  9H42_A   9H42_A
    815  2GS2_A   2GS2_A
    816  5FED_A   5FED_A
    817  9UBI_A   9UBI_A
    818  4WKQ_A   4WKQ_A
    819  3D5X_A   3D5X_A
    820  6JZ0_A   6JZ0_A
    821  6CN9_A   6CN9_A
    822  7UYV_A   7UYV_A
    823  7Q6H_a 7Q6H_AAA
    824  2J5E_A   2J5E_A
    825  8A27_A   8A27_A
    826  7AEM_A   7AEM_A
    827  8PO3_A   8PO3_A
    828  2RGP_A   2RGP_A
    829  8PO4_A   8PO4_A
    830  1XKK_A   1XKK_A
    831  3ZC6_A   3ZC6_A
    832  4PY1_A   4PY1_A
    833  5W86_A   5W86_A
    834  8H7X_A   8H7X_A
    835  9N6G_A   9N6G_A
    836  5NG0_A   5NG0_A
    837  5XGN_A   5XGN_A
    838  9BY4_A   9BY4_A
    839  7XDY_A   7XDY_A
    840  7XDV_A   7XDV_A
    841  7XDX_A   7XDX_A
    842  2GS7_A   2GS7_A
    843  5TOZ_A   5TOZ_A
    844  6AAM_A   6AAM_A
    845  5W5J_A   5W5J_A
    846  3LXK_A   3LXK_A
    847  9QBG_A   9QBG_A
    848  9QBF_A   9QBF_A
    849  6HMX_A   6HMX_A
    850  4HJO_A   4HJO_A
    851  4TPT_A   4TPT_A
    852  8X2O_A   8X2O_A
    853  9F3V_A   9F3V_A
    854  6ES0_A   6ES0_A
    855  6UL8_A   6UL8_A
    856  5NXD_A   5NXD_A
    857  9QXN_a 9QXN_AAA
    858  5AR2_A   5AR2_A
    859  8AZA_A   8AZA_A
    860  5NG3_B   5NG3_B
    861  6HZV_A   6HZV_A
    862  3LZB_A   3LZB_A
    863  4C8B_A   4C8B_A
    864  5ZWJ_A   5ZWJ_A
    865  6LUB_A   6LUB_A
    866  8HV4_A   8HV4_A
    867  2V5Q_A   2V5Q_A
    868  4E4L_A   4E4L_A
    869  6N7A_A   6N7A_A
    870  3IKA_A   3IKA_A
    871  4GIH_A   4GIH_A
    872  6GTT_A   6GTT_A
    873  8WD4_A   8WD4_A
    874  9D3V_A   9D3V_A
    875  5TD2_A   5TD2_A
    876  6ELR_A   6ELR_A
    877  7TVD_A   7TVD_A
    878  6GGH_A   6GGH_A
    879  3PJC_A   3PJC_A
    880  2YAC_A   2YAC_A
    881  9S3X_A   9S3X_A
    882  3PP0_A   3PP0_A
    883  7PCD_A   7PCD_A
    884  2QKW_B   2QKW_B
    885  9QEK_A   9QEK_A
    886  3ZBF_A   3ZBF_A
    887  7JXH_A   7JXH_A
    888  3KB7_A   3KB7_A
    889  2JIU_A   2JIU_A
    890  4I24_A   4I24_A
    891  5KHW_A   5KHW_A
    892  7SZ0_A   7SZ0_A
    893  5GMP_A   5GMP_A
    894  5Y9T_A   5Y9T_A
    895  5FEE_A   5FEE_A
    896  3EYG_A   3EYG_A
    897  4LQM_A   4LQM_A
    898  3HGK_A   3HGK_A
    899  6TPE_A   6TPE_A
    900  5GNK_A   5GNK_A
    901  5AJQ_A   5AJQ_A
    902  4G5P_A   4G5P_A
    903  2J7T_A   2J7T_A
    904  2JIT_A   2JIT_A
    905  6EIM_A   6EIM_A
    906  4BC6_A   4BC6_A
    907  5XDL_A   5XDL_A
    908  4YNE_A   4YNE_A
    909  2OU7_A   2OU7_A
    910  4NZW_B   4NZW_B
    911  5J9Z_A   5J9Z_A
    912  5J9Y_A   5J9Y_A
    913  4ZSE_A   4ZSE_A
    914  5NG3_A   5NG3_A
    915  2EB2_A   2EB2_A
    916  4QPS_A   4QPS_A
    917  2ITN_A   2ITN_A
    918  6JWL_A   6JWL_A
    919  3THB_A   3THB_A
    920  3GGF_A   3GGF_A
    921  6HXF_A   6HXF_A
    922  7C3N_A   7C3N_A
    923  4J52_A   4J52_A
    924  2EB3_A   2EB3_A
    925  6P1D_A   6P1D_A
    926  2ITT_A   2ITT_A
    927  4H1J_A   4H1J_A
    928  9JQ1_A   9JQ1_A
    929  6V5N_A   6V5N_A
    930  2RKU_A   2RKU_A
    931  6C9D_A   6C9D_A
    932  6TFU_A   6TFU_A
    933  8PEH_A   8PEH_A
    934  3ET7_A   3ET7_A
    935  4R3P_A   4R3P_A
    936  5TO8_A   5TO8_A
    937  5LWM_A   5LWM_A
    938  8S9P_C   8S9P_C
    939  7K1H_A   7K1H_A
    940  9P9U_A   9P9U_A
    941  9P9U_A   9P9U_A
    942  9XU9_A   9XU9_A
    943  3CC6_A   3CC6_A
    944  4I20_A   4I20_A
    945  4RJ4_A   4RJ4_A
    946  5JFS_A   5JFS_A
    947  7UYR_A   7UYR_A
    948  7VKO_A   7VKO_A
    949  6D22_A   6D22_A
    950  8HY7_A   8HY7_A
    951  4RIY_A   4RIY_A
    952  7VKM_A   7VKM_A
    953  5TA6_A   5TA6_A
    954  3GOP_A   3GOP_A
    955  9DF3_A   9DF3_A
    956  9DF2_A   9DF2_A
    957  4KS7_A   4KS7_A
    958  6C7Y_A   6C7Y_A
    959  6S9B_A   6S9B_A
    960  6S9C_A   6S9C_A
    961  4J7B_A   4J7B_A
    962  5H3Q_A   5H3Q_A
    963  7MN5_B   7MN5_B
    964  4AOJ_A   4AOJ_A
    965  7MN6_B   7MN6_B
    966  8A2B_A   8A2B_A
    967  5KMI_A   5KMI_A
    968  3DAK_A   3DAK_A
    969  5KML_A   5KML_A
    970  2HAK_A   2HAK_A
    971  3NZ0_A   3NZ0_A
    972  2VWI_A   2VWI_A
    973  2C30_A   2C30_A
    974  4ZP5_A   4ZP5_A
    975  8DSW_A   8DSW_A
    976  4OBO_A   4OBO_A
    977  8J5W_A   8J5W_A
    978  5DI1_A   5DI1_A
    979  5J95_A   5J95_A
    980  4F0I_A   4F0I_A
    981  2QNJ_A   2QNJ_A
    982  4U3Z_A   4U3Z_A
    983  4LRM_A   4LRM_A
    984  6NSP_A   6NSP_A
    985  6D1Y_A   6D1Y_A
    986  7XAF_A   7XAF_A
    987  6IQN_A   6IQN_A
    988  2A19_B   2A19_B
    989  4RVT_A   4RVT_A
    990  4E1Z_A   4E1Z_A
    991  4E20_A   4E20_A
    992  4GT5_A   4GT5_A
    993  7MN5_A   7MN5_A
    994  8J5X_A   8J5X_A
    995  4LL0_A   4LL0_A
    996  8PO0_A   8PO0_A
    997  8UOI_A   8UOI_A
    998  7AAX_A   7AAX_A
    999  4I5M_A   4I5M_A
    1000 9U8C_A   9U8C_A
    1001 7M5Z_A   7M5Z_A
    1002 6NSS_A   6NSS_A
    1003 4I21_A   4I21_A
    1004 8V5I_A   8V5I_A
    1005 3LMG_A   3LMG_A
    1006 3UG1_A   3UG1_A
    1007 7P1L_A   7P1L_A
    1008 4RIX_A   4RIX_A
    1009 4RIW_A   4RIW_A
    1010 7FEH_A   7FEH_A
    1011 8PO1_A   8PO1_A
    1012 5F1Z_A   5F1Z_A
    1013 6YAT_A   6YAT_A
    1014 3NYX_A   3NYX_A
    1015 9KLW_A   9KLW_A
    1016 6OP9_A   6OP9_A
    1017 5SAU_A   5SAU_A
    1018 7MX3_A   7MX3_A
    1019 8C12_a 8C12_AAA
    1020 8UOH_A   8UOH_A
    1021 3FE3_A   3FE3_A
    1022 3COM_A   3COM_A
    1023 6CTH_A   6CTH_A
    1024 7OXB_A   7OXB_A
    1025 8UOJ_A   8UOJ_A
    1026 9R1W_a 9R1W_AAA
    1027 3W2O_A   3W2O_A
    1028 6M0U_A   6M0U_A
    1029 3KEX_A   3KEX_A
    1030 8PAV_A   8PAV_A
    1031 7LGS_A   7LGS_A
    1032 4I1Z_A   4I1Z_A
    1033 5Y25_A   5Y25_A
    1034 5WNI_A   5WNI_A
    1035 8KFQ_A   8KFQ_A
    1036 9FQP_A   9FQP_A
    1037 7MON_B   7MON_B
    1038 2F57_A   2F57_A
    1039 8U8X_A   8U8X_A
    1040 4FZA_B   4FZA_B
    1041 8VB5_A   8VB5_A
    1042 9IIC_A   9IIC_A
    1043 5WR7_A   5WR7_A
    1044 8D73_A   8D73_A
    1045 6NPT_A   6NPT_A
    1046 4PMM_A   4PMM_A
    1047 9LFU_A   9LFU_A
    1048 6D3K_A   6D3K_A
    1049 2WZJ_A   2WZJ_A
    1050 6PL1_A   6PL1_A
    1051 4B6L_A   4B6L_A
    1052 2MSE_D   2MSE_D
    1053 2R0I_A   2R0I_A
    1054 4OLI_A   4OLI_A
    1055 1ZMU_A   1ZMU_A
    1056 5EAK_A   5EAK_A
    1057 1WXM_A   1WXM_A
    1058 1RFA_A   1RFA_A
    1059 6JRK_A   6JRK_A
    1060 3UIU_A   3UIU_A
    1061 5DBX_A   5DBX_A
    1062 8A5J_A   8A5J_A
    1063 8QEL_B   8QEL_B
    1064 5BVK_A   5BVK_A
    1065 3ZOS_A   3ZOS_A
    1066 5D9H_A   5D9H_A
    1067 1ZMW_A   1ZMW_A
    1068 6S89_A   6S89_A
    1069 1ZMV_A   1ZMV_A
    1070 1C1Y_B   1C1Y_B
    1071 4G0N_B   4G0N_B
    1072 5KZ7_A   5KZ7_A
    1073 6JRJ_A   6JRJ_A
    1074 6VJJ_B   6VJJ_B
    1075 6Y23_A   6Y23_A
    1076 6BRJ_A   6BRJ_A
    1077 1GUA_B   1GUA_B
    1078 5E8U_A   5E8U_A
    1079 6VC0_A   6VC0_A
    1080 5FDP_A   5FDP_A
    1081 8JOF_A   8JOF_A
    1082 4ASZ_A   4ASZ_A
    1083 8TXY_A   8TXY_A
    1084 1RW8_A   1RW8_A
    1085 1B6C_B   1B6C_B
    1086 5USQ_A   5USQ_A
    1087 9F6X_A   9F6X_A
    1088 8YHF_A   8YHF_A
    1089 9J9D_A   9J9D_A
    1090 1VJY_A   1VJY_A
    1091 5DH3_A   5DH3_A
    1092 8A66_B   8A66_B
    1093 3TZM_A   3TZM_A
    1094 1PY5_A   1PY5_A
    1095 4X0M_A   4X0M_A
    1096 5FRI_A   5FRI_A
    1097 8A66_A   8A66_A
    1098 2WOT_A   2WOT_A
    1099 5I8A_A   5I8A_A
    1100 5KVT_A   5KVT_A
    1101 5E8T_A   5E8T_A
    1102 4YNZ_A   4YNZ_A
    1103 4YOM_B   4YOM_B
    1104 8XFL_A   8XFL_A
    1105 1RRB_A   1RRB_A
    1106 5ES1_A   5ES1_A
    1107 3IEC_A   3IEC_A
    1108 5LPZ_A   5LPZ_A
    1109 4M68_A   4M68_A
    1110 6FER_A   6FER_A
    1111 3KUD_B   3KUD_B
    1112 1Y8G_A   1Y8G_A
    1113 3KUC_B   3KUC_B
    1114 3V5Q_A   3V5Q_A
    1115 4BFM_A   4BFM_A
    1116 5LPV_A   5LPV_A
    1117 5LPB_A   5LPB_A
    1118 7AYM_A   7AYM_A
    1119 6N3N_A   6N3N_A
    1120 4OH4_A   4OH4_A
    1121 4CQG_A   4CQG_A
    1122 4LG4_A   4LG4_A
    1123 5XD6_A   5XD6_A
    1124 4YUR_A   4YUR_A
    1125 3MDY_A   3MDY_A
    1126 5WNO_A   5WNO_A
    1127 4XUF_A   4XUF_A
    1128 9Y9N_A   9Y9N_A
    1129 6AO5_A   6AO5_A
    1130 3COK_A   3COK_A
    1131 4JXF_A   4JXF_A
    1132 4LGD_A   4LGD_A
    1133 2OO8_X   2OO8_X
    1134 3TL8_A   3TL8_A
    1135 6JQR_A   6JQR_A
    1136 6OYW_A   6OYW_A
    1137 1FVR_A   1FVR_A
    1138 4BIB_A   4BIB_A
    1139 6VRE_A   6VRE_A
    1140 2BMC_A   2BMC_A
    1141 3UIM_A   3UIM_A
    1142 4YMJ_A   4YMJ_A
    1143 1U59_A   1U59_A
    1144 3COH_A   3COH_A
    1145 3QBN_A   3QBN_A
    1146 4PRJ_A   4PRJ_A
    1147 8XB1_A   8XB1_A
    1148 5X02_A   5X02_A
    1149 6BFN_A   6BFN_A
    1150 1RJB_A   1RJB_A
    1151 3R21_A   3R21_A
    1152 4BF2_A   4BF2_A
    1153 1GZK_A   1GZK_A
    1154 7T6F_A   7T6F_A
    1155 3D0E_A   3D0E_A
    1156 6N3L_A   6N3L_A
    1157 2XRU_A   2XRU_A
    1158 1O6K_A   1O6K_A
    1159 1GZN_A   1GZN_A
    1160 1MRV_A   1MRV_A
    1161 6E2M_A   6E2M_A
    1162 2JED_A   2JED_A
    1163 5VIL_A   5VIL_A
    1164 6IL3_A   6IL3_A
    1165 2J4Z_A   2J4Z_A
    1166 2JDO_A   2JDO_A
    1167 5F9E_A   5F9E_A
    1168 7ZTL_A   7ZTL_A
    1169 3VW6_A   3VW6_A
    1170 2CLQ_A   2CLQ_A
    1171 4Q9Z_A   4Q9Z_A
    1172 6HJJ_A   6HJJ_A
    1173 6SEQ_A   6SEQ_A
    1174 3UNZ_A   3UNZ_A
    1175 3NRM_A   3NRM_A
    1176 1O6L_A   1O6L_A
    1177 6XIH_A   6XIH_A
    1178 5UOR_A   5UOR_A
    1179 8C1M_A   8C1M_A
    1180 4PX6_A   4PX6_A
    1181 5OBJ_A   5OBJ_A
    1182 8Q61_A   8Q61_A
    1183 8GUW_A   8GUW_A
    1184 3EMG_A   3EMG_A
    1185 5DN3_A   5DN3_A
    1186 4J8M_A   4J8M_A
    1187 1MQ4_A   1MQ4_A
    1188 2WQB_A   2WQB_A
    1189 5AAD_A   5AAD_A
    1190 4RX7_A   4RX7_A
    1191 2OZO_A   2OZO_A
    1192 5UOX_A   5UOX_A
    1193 4RX9_A   4RX9_A
    1194 6OYT_A   6OYT_A
    1195 3TUB_A   3TUB_A
    1196 9BZG_A   9BZG_A
    1197 1XJD_A   1XJD_A
    1198 4CEG_A   4CEG_A
    1199 4Q5J_A   4Q5J_A
    1200 3TUC_A   3TUC_A
    1201 4YJO_A   4YJO_A
    1202 6J5T_A   6J5T_A
    1203 9C1W_A   9C1W_A
    1204 8OF5_A   8OF5_A
    1205 5LXM_A   5LXM_A
    1206 4X3J_A   4X3J_A
    1207 4K2R_A   4K2R_A
    1208 6VPJ_A   6VPJ_A
    1209 8X5K_A   8X5K_A
    1210 2J50_A   2J50_A
    1211 4F4P_A   4F4P_A
    1212 1MUO_A   1MUO_A
    1213 5ORL_A   5ORL_A
    1214 3SRV_B   3SRV_B
    1215 3FDN_A   3FDN_A
    1216 4BTF_A   4BTF_A
    1217 3W16_A   3W16_A
    1218 4C3P_A   4C3P_A
    1219 5OS2_A   5OS2_A
    1220 7QQ6_A   7QQ6_A
    1221 4RSS_A   4RSS_A
    1222 8SSP_A   8SSP_A
    1223 8FO7_C   8FO7_C
    1224 3W10_A   3W10_A
    1225 8SMC_C   8SMC_C
    1226 5OSD_A   5OSD_A
    1227 7QWK_A   7QWK_A
    1228 8C1K_A   8C1K_A
    1229 5TR6_A   5TR6_A
    1230 6C83_A   6C83_A
    1231 1XBA_A   1XBA_A
    1232 5ZAN_A   5ZAN_A
    1233 6I2U_A   6I2U_A
    1234 5OS5_A   5OS5_A
    1235 4PV0_A   4PV0_A
    1236 2X6D_A   2X6D_A
    1237 8C15_A   8C15_A
    1238 6VNO_A   6VNO_A
    1239 5C26_A   5C26_A
    1240 8TXZ_A   8TXZ_A
    1241 9DMI_A   9DMI_A
    1242 1OL5_A   1OL5_A
    1243 2WTW_A   2WTW_A
    1244 6VOV_A   6VOV_A
    1245 3LAU_A   3LAU_A
    1246 2WTV_A   2WTV_A
    1247 4YJQ_A   4YJQ_A
    1248 2C6E_A   2C6E_A
    1249 6VP8_A   6VP8_A
    1250 2C6D_A   2C6D_A
    1251 2XNG_A   2XNG_A
    1252 5Y5T_A   5Y5T_A
    1253 4FL3_A   4FL3_A
    1254 9KDS_A   9KDS_A
    1255 4DFL_A   4DFL_A
    1256 8RRQ_A   8RRQ_A
    1257 8PS7_A   8PS7_A
    1258 7O2V_A   7O2V_A
    1259 4FL2_A   4FL2_A
    1260 4BN1_A   4BN1_A
    1261 8HOD_A   8HOD_A
    1262 4UZD_A   4UZD_A
    1263 3E5A_A   3E5A_A
    1264 7LHT_A   7LHT_A
    1265 8TZB_A   8TZB_A
    1266 3HA6_A   3HA6_A
    1267 4JAI_A   4JAI_A
    1268 8OKU_A   8OKU_A
    1269 5DOS_A   5DOS_A
    1270 3EFW_A   3EFW_A
    1271 5EW9_A   5EW9_A
    1272 3H0Y_A   3H0Y_A
    1273 5DT4_A   5DT4_A
    1274 9CI3_A   9CI3_A
    1275 5DT3_A   5DT3_A
    1276 9D8Z_A   9D8Z_A
    1277 8R4O_A   8R4O_A
    1278 6HM6_A   6HM6_A
    1279 3SRV_A   3SRV_A
    1280 6VPG_A   6VPG_A
    1281 3H9R_A   3H9R_A
    1282 9RDA_A   9RDA_A
    1283 6VPH_A   6VPH_A
    1284 6XKA_A   6XKA_A
    1285 9KS6_A   9KS6_A
    1286 9CHO_A   9CHO_A
    1287 6VPL_A   6VPL_A
    1288 9D8F_A   9D8F_A
    1289 6JUX_A   6JUX_A
    1290 8HO6_A   8HO6_A
    1291 4DYM_A   4DYM_A
    1292 5TOS_A   5TOS_A
    1293 9ESA_a 9ESA_AAA
    1294 6VPI_A   6VPI_A
    1295 8JF4_A   8JF4_A
    1296 9P6A_A   9P6A_A
    1297 8TZG_A   8TZG_A
    1298 3W2C_A   3W2C_A
    1299 2W1C_A   2W1C_A
    1300 8HOA_A   8HOA_A
    1301 2H6D_A   2H6D_A
    1302 2W1D_A   2W1D_A
    1303 6UNQ_A   6UNQ_A
    1304 6UNR_A   6UNR_A
    1305 6GVX_A   6GVX_A
    1306 9F31_A   9F31_A
    1307 4BKY_A   4BKY_A
    1308 5TWU_A   5TWU_A
    1309 7CTX_A   7CTX_A
    1310 6XR4_A   6XR4_A
    1311 6EIX_A   6EIX_A
    1312 5TWL_A   5TWL_A
    1313 2WQE_A   2WQE_A
    1314 3MY0_A   3MY0_A
    1315 2YZA_A   2YZA_A
    1316 6VXR_A   6VXR_A
    1317 3DAJ_A   3DAJ_A
    1318 6KZC_A   6KZC_A
    1319 4ZHX_A   4ZHX_A
    1320 5TVT_A   5TVT_A
    1321 5ISO_A   5ISO_A
    1322 4CFE_A   4CFE_A
    1323 2Y7J_A   2Y7J_A
    1324 8TZC_A   8TZC_A
    1325 8BIK_A   8BIK_A
    1326 7KX8_A   7KX8_A
    1327 5IH8_A   5IH8_A
    1328 7F3G_A   7F3G_A
    1329 7KXW_A   7KXW_A
    1330 5JZN_A   5JZN_A
    1331 7LI3_A   7LI3_A
    1332 4M69_A   4M69_A
    1333 4UMQ_A   4UMQ_A
    1334 8TYQ_A   8TYQ_A
    1335 3MTF_A   3MTF_A
    1336 6VBZ_A   6VBZ_A
    1337 9IC2_A   9IC2_A
    1338 4UMT_A   4UMT_A
    1339 7OPO_A   7OPO_A
    1340 5EZV_A   5EZV_A
    1341 5D9K_A   5D9K_A
    1342 3G51_A   3G51_A
    1343 5JZJ_A   5JZJ_A
    1344 4QFG_A   4QFG_A
    1345 4NUS_A   4NUS_A
    1346 3D14_A   3D14_A
    1347 6NPZ_A   6NPZ_A
    1348 8XFY_A   8XFY_A
    1349 3CQU_A   3CQU_A
    1350 9S9T_A   9S9T_A
    1351 6BX6_A   6BX6_A
    1352 8UWR_A   8UWR_A
    1353 6KYQ_A   6KYQ_A
    1354 4GV1_A   4GV1_A
    1355 4RER_A   4RER_A
    1356 4REW_A   4REW_A
    1357 3OCB_A   3OCB_A
    1358 6BUU_A   6BUU_A
    1359 2I0E_A   2I0E_A
    1360 4EL9_A   4EL9_A
    1361 6KYR_A   6KYR_A
    1362 6CCY_A   6CCY_A
    1363 3UBD_A   3UBD_A
    1364 4D2P_A   4D2P_A
    1365 2Z7Q_A   2Z7Q_A
    1366 3OHT_A   3OHT_A
    1367 4CFH_A   4CFH_A
    1368 4CZU_A   4CZU_A
    1369 6GR8_A   6GR8_A
    1370 4RED_A   4RED_A
    1371 4CZT_A   4CZT_A
    1372 7JIJ_A   7JIJ_A
    1373 5DE2_A   5DE2_A
    1374 8DFP_A   8DFP_A
    1376 7JHG_A   7JHG_A
    1377 3ZDU_A   3ZDU_A
    1378 6C9F_A   6C9F_A
    1379 8DFQ_A   8DFQ_A
    1381 8DFM_A   8DFM_A
    1383 6NPY_B   6NPY_B
    1384 6C9H_A   6C9H_A
    1385 2DWB_A   2DWB_A
    1386 1ZYC_A   1ZYC_A
    1387 8SXN_A   8SXN_A
    1388 1ZY4_A   1ZY4_A
    1389 8E05_A   8E05_A
    1390 2FH9_A   2FH9_A
    1391 3MN3_A   3MN3_A
    1392 3HYH_A   3HYH_A
    1393 3O96_A   3O96_A
    1394 3MTL_A   3MTL_A
    1395 4EJN_A   4EJN_A
    1396 8E04_A   8E04_A
    1397 3DAE_A   3DAE_A
    1398 8FAC_A   8FAC_A
    1399 6S73_A   6S73_A
    1400 4M66_A   4M66_A
    1401 9H59_C   9H59_C
    1402 9NFQ_C   9NFQ_C
    1403 8WS0_A   8WS0_A
    1404 4FSN_A   4FSN_A
    1405 4A4X_A   4A4X_A
    1406 1IA8_A   1IA8_A
    1407 4QYE_A   4QYE_A
    1408 5KCV_A   5KCV_A
    1409 7AKO_A   7AKO_A
    1410 2WQM_A   2WQM_A
    1411 7AKM_A   7AKM_A
    1412 1ZLT_A   1ZLT_A
    1413 4FSM_A   4FSM_A
    1414 4FSY_A   4FSY_A
    1415 2X8E_A   2X8E_A
    1416 2BR1_A   2BR1_A
    1417 8QGY_A   8QGY_A
    1418 2HOG_A   2HOG_A
    1419 5OQ5_A   5OQ5_A
    1420 2QHM_A   2QHM_A
    1421 7APJ_A   7APJ_A
    1422 3ORM_A   3ORM_A
    1423 4FSW_A   4FSW_A
    1424 5F4N_A   5F4N_A
    1425 4D28_A   4D28_A
    1426 3JVR_A   3JVR_A
    1427 6TM5_S   6TM5_S
    1428 2ACX_A   2ACX_A
    1429 2E9V_A   2E9V_A
    1430 4FSZ_A   4FSZ_A
    1431 3OT3_A   3OT3_A
    1432 8UW7_A   8UW7_A
    1433 3NYN_A   3NYN_A
    1434 3IW4_A   3IW4_A
    1435 2AYP_A   2AYP_A
    1436 8UVY_A   8UVY_A
    1437 6TM5_Q   6TM5_Q
    1438 2E9N_A   2E9N_A
    1439 2YDJ_A   2YDJ_A
    1440 4FT3_A   4FT3_A
    1441 6TUA_A   6TUA_A
    1442 8EJ4_K   8EJ4_K
    1443 4FST_A   4FST_A
    1444 1MRU_A   1MRU_A
    1445 2W5A_A   2W5A_A
    1446 4FR4_A   4FR4_A
    1447 2JAV_A   2JAV_A
    1448 8A3T_S   8A3T_S
    1449 5U94_A   5U94_A
    1450 4AF3_A   4AF3_A
    1451 3F69_A   3F69_A
    1452 6B2P_A   6B2P_A
    1453 4RA4_A   4RA4_A
    1454 6I2P_A   6I2P_A
    1455 7PUE_A   7PUE_A
    1456 3ORI_A   3ORI_A
    1457 1ZYS_A   1ZYS_A
    1458 3F61_A   3F61_A
    1459 4G31_A   4G31_A
    1460 1ZXE_A   1ZXE_A
    1461 1O6Y_A   1O6Y_A
    1462 8WS1_A   8WS1_A
    1463 4FSU_A   4FSU_A
    1464 2R5T_A   2R5T_A
    1465 9IWX_A   9IWX_A
    1466 6G78_A   6G78_A
    1467 5I3O_A   5I3O_A
    1468 4PNI_A   4PNI_A
    1469 4W9W_A   4W9W_A
    1470 8VSU_C   8VSU_C
    1471 4W9X_A   4W9X_A
    1472 3QC9_A   3QC9_A
    1473 3C4X_A   3C4X_A
    1474 3C4W_A   3C4W_A
    1475 3T8O_A   3T8O_A
    1476 6G76_A   6G76_A
    1477 6G77_A   6G77_A
    1478 7MT8_G   7MT8_G
    1479 4WBO_A   4WBO_A
    1480 6NTD_B   6NTD_B
    1481 4E5A_X   4E5A_X
    1482 3H4J_A   3H4J_A
    1483 4L9I_A   4L9I_A
    1484 2PUU_A   2PUU_A
    1485 3HNG_A   3HNG_A
    1486 2FSL_X   2FSL_X
    1487 3P4K_A   3P4K_A
    1488 8EFJ_A   8EFJ_A
    1489 1LEW_A   1LEW_A
    1490 2FST_X   2FST_X
    1491 4EQM_A   4EQM_A
    1492 4LOO_A   4LOO_A
    1493 3VHE_A   3VHE_A
    1494 3VHK_A   3VHK_A
    1495 2GTM_A   2GTM_A
    1496 3NNX_A   3NNX_A
    1497 1BMK_A   1BMK_A
    1498 5O90_A   5O90_A
    1499 5MRD_A   5MRD_A
    1500 3TG1_A   3TG1_A
    1501 3VID_A   3VID_A
    1502 2FSO_X   2FSO_X
    1503 1Y6A_A   1Y6A_A
    1504 2BAQ_A   2BAQ_A
    1505 4TYH_B   4TYH_B
    1506 8YPE_A   8YPE_A
    1507 5O8U_A   5O8U_A
    1508 5O8V_A   5O8V_A
    1509 3ODZ_X   3ODZ_X
    1510 8ACM_a 8ACM_AAA
    1511 3PWY_A   3PWY_A
    1512 3HVC_A   3HVC_A
    1513 3S3I_A   3S3I_A
    1514 3OD6_X   3OD6_X
    1515 3NNU_A   3NNU_A
    1516 6SO1_A   6SO1_A
    1517 3ODY_X   3ODY_X
    1518 5NZZ_E   5NZZ_E
    1519 3FI4_A   3FI4_A
    1520 3HEC_A   3HEC_A
    1521 6SOI_A   6SOI_A
    1522 1A9U_A   1A9U_A
    1523 9MHB_A   9MHB_A
    1524 1YW2_A   1YW2_A
    1525 2YIS_A   2YIS_A
    1526 1OZ1_A   1OZ1_A
    1527 2GFS_A   2GFS_A
    1528 3K3I_A   3K3I_A
    1529 8X3M_A   8X3M_A
    1530 3GCU_A   3GCU_A
    1531 3OEF_X   3OEF_X
    1532 1YWR_A   1YWR_A
    1533 2GHL_A   2GHL_A
    1534 3KL8_A   3KL8_A
    1535 2BAJ_A   2BAJ_A
    1536 3K3J_A   3K3J_A
    1537 3ZS5_A   3ZS5_A
    1538 1DI9_A   1DI9_A
    1539 1ZZL_A   1ZZL_A
    1540 2NPQ_A   2NPQ_A
    1541 3KQ7_A   3KQ7_A
    1542 3MPT_A   3MPT_A
    1543 8VXE_A   8VXE_A
    1544 5ETC_A   5ETC_A
    1545 6KA4_A   6KA4_A
    1546 3PY3_A   3PY3_A
    1547 2BAL_A   2BAL_A
    1548 6TCA_B   6TCA_B
    1549 2OZA_B   2OZA_B
    1550 3KK9_A   3KK9_A
    1551 1M7Q_A   1M7Q_A
    1552 3E92_A   3E92_A
    1553 3KK8_A   3KK8_A
    1554 2BDW_A   2BDW_A
    1555 6ZQS_A   6ZQS_A
    1556 3D7Z_A   3D7Z_A
    1557 9CJ1_A   9CJ1_A
    1558 4F9W_A   4F9W_A
    1559 4X7H_A   4X7H_A
    1560 2LGC_A   2LGC_A
    1561 6Y7W_A   6Y7W_A
    1562 1H1W_A   1H1W_A
    1563 6SOV_A   6SOV_A
    1564 4R3C_A   4R3C_A
    1565 3DT1_A   3DT1_A
    1566 3D83_A   3D83_A
    1567 7PVU_A   7PVU_A
    1568 7Z6I_a 7Z6I_AAA
    1569 9NYT_A   9NYT_A
    1570 4F9Y_A   4F9Y_A
    1571 1IAN_A   1IAN_A
    1572 1UU9_A   1UU9_A
    1573 1V0O_A   1V0O_A
    1574 3ION_A   3ION_A
    1575 1V0B_A   1V0B_A
    1576 2Y8O_A   2Y8O_A
    1577 4CT1_A   4CT1_A
    1578 5ETF_A   5ETF_A
    1579 6QYX_A   6QYX_A
    1580 1FOT_A   1FOT_A
    1581 2R7B_A   2R7B_A
    1582 5WJJ_A   5WJJ_A
    1583 3H9O_A   3H9O_A
    1584 4A07_A   4A07_A
    1585 3NUS_A   3NUS_A
    1586 5L4Q_A   5L4Q_A
    1587 3ORX_A   3ORX_A
    1588 1Z5M_A   1Z5M_A
    1589 5MZ3_A   5MZ3_A
    1590 9YLH_A   9YLH_A
    1591 5TE0_A   5TE0_A
    1592 3QC4_A   3QC4_A
    1593 3HRC_A   3HRC_A
    1594 5ETI_A   5ETI_A
    1595 2XCH_A   2XCH_A
    1596 1OKY_A   1OKY_A
    1597 1OVE_A   1OVE_A
    1598 7LVH_A   7LVH_A
    1599 4XX9_A   4XX9_A
    1600 4EWQ_A   4EWQ_A
    1601 1OB3_A   1OB3_A
    1602 3Q4T_A   3Q4T_A
    1603 6M9L_A   6M9L_A
    1604 3RWQ_A   3RWQ_A
    1605 2BIY_A   2BIY_A
    1606 6M95_A   6M95_A
    1607 3RWP_A   3RWP_A
    1608 4WSQ_A   4WSQ_A
    1609 7XBR_F   7XBR_F
    1610 3MH2_A   3MH2_A
    1611 2WTK_C   2WTK_C
    1612 4IC7_A   4IC7_A
    1613 3HKO_A   3HKO_A
    1614 8A8M_A   8A8M_A
    1615 9QB5_A   9QB5_A
    1616 4KIK_A   4KIK_A
    1617 3NUN_A   3NUN_A
    1618 4B99_A   4B99_A
    1619 2PHK_A   2PHK_A
    1620 2QLU_A   2QLU_A
    1621 4GEO_A   4GEO_A
    1622 9F58_A   9F58_A
    1623 3GCP_A   3GCP_A
    1624 1PHK_A   1PHK_A
    1625 9CMZ_A   9CMZ_A
    1626 8GMC_A   8GMC_A
    1627 8OMV_A   8OMV_A
    1628 4E3C_A   4E3C_A
    1629 8U2O_A   8U2O_A
    1630 3NAX_A   3NAX_A
    1631 2XCK_A   2XCK_A
    1632 9BF3_A   9BF3_A
    1633 3NAY_A   3NAY_A
    1634 2W96_B   2W96_B
    1635 3V3V_A   3V3V_A
    1636 3WE4_A   3WE4_A
    1637 3A60_A   3A60_A
    1638 3GP0_A   3GP0_A
    1639 3A62_A   3A62_A
    1640 4KIK_B   4KIK_B
    1641 4L45_A   4L45_A
    1642 4L46_A   4L46_A
    1643 4L43_A   4L43_A
    1644 7N91_A   7N91_A
    1645 8YGW_A   8YGW_A
    1646 7XBR_A   7XBR_A
    1647 1QL6_A   1QL6_A
    1648 4ZSJ_A   4ZSJ_A
    1649 1CM8_A   1CM8_A
    1650 4L42_A   4L42_A
    1651 5BYY_A   5BYY_A
    1652 6P8E_B   6P8E_B
    1653 2W99_B   2W99_B
    1654 5BYZ_A   5BYZ_A
    1655 6HKM_A   6HKM_A
    1656 4ZSG_A   4ZSG_A
    1657 8T7T_A   8T7T_A
    1659 4L3J_A   4L3J_A
    1660 5O7I_A   5O7I_A
    1661 4Y83_A   4Y83_A
    1662 4AGU_A   4AGU_A
    1663 6UNA_A   6UNA_A
    1664 4C57_A   4C57_A
    1665 5LOH_A   5LOH_A
    1666 9CSK_B   9CSK_B
    1667 7CGA_A   7CGA_A
    1668 6V6A_A   6V6A_A
    1669 4RLO_A   4RLO_A
    1670 3GC8_A   3GC8_A
    1671 2CN5_A   2CN5_A
    1672 9LTA_A   9LTA_A
    1673 2W9F_B   2W9F_B
    1674 2XK9_A   2XK9_A
    1675 6LBA_A   6LBA_A
    1676 5Y8U_A   5Y8U_A
    1677 2YCF_A   2YCF_A
    1678 2YCR_A   2YCR_A
    1679 2BFX_A   2BFX_A
    1680 5Z1D_A   5Z1D_A
    1681 4C2W_A   4C2W_A
    1682 2W0J_A   2W0J_A
    1683 4B8M_A   4B8M_A
    1684 2VRX_A   2VRX_A
    1685 3TXO_A   3TXO_A
    1686 4C2V_A   4C2V_A
    1687 8FP1_A   8FP1_A
    1688 3COI_A   3COI_A
    1689 3GC9_A   3GC9_A
    1690 4B8L_A   4B8L_A
    1691 6IB0_A   6IB0_A
    1692 2DYL_A   2DYL_A
    1693 5Y90_A   5Y90_A
    1694 7OVJ_A   7OVJ_A
    1695 2XRW_A   2XRW_A
    1696 3DLS_A   3DLS_A
    1697 3NIZ_A   3NIZ_A
    1698 4ZSL_A   4ZSL_A
    1699 6YG4_A   6YG4_A
    1700 6YFZ_A   6YFZ_A
    1701 6HKN_A   6HKN_A
    1702 5B2L_A   5B2L_A
    1703 3ALN_A   3ALN_A
    1704 3WZU_A   3WZU_A
    1705 3J4Q_D   3J4Q_D
    1706 3FHI_A   3FHI_A
    1707 8V5H_A   8V5H_A
    1708 6YG1_A   6YG1_A
    1709 4EYJ_A   4EYJ_A
    1710 7DV6_A   7DV6_A
    1711 6ZR5_A   6ZR5_A
    1712 3G33_A   3G33_A
    1713 2QUR_A   2QUR_A
    1714 2QKR_A   2QKR_A
    1715 8X23_A   8X23_A
    1716 2XS0_A   2XS0_A
    1717 8YGZ_A   8YGZ_A
    1718 7SJ3_A   7SJ3_A
    1719 5Y7Z_A   5Y7Z_A
    1720 5FWK_K   5FWK_K
    1721 1UKH_A   1UKH_A
    1722 6YG0_A   6YG0_A
    1723 8X5M_A   8X5M_A
    1724 4YR8_A   4YR8_A
    1725 3QAL_E   3QAL_E
    1726 1ATP_E   1ATP_E
    1727 5E8V_A   5E8V_A
    1728 2ERZ_E   2ERZ_E
    1729 7E0Z_A   7E0Z_A
    1730 1J3H_A   1J3H_A
    1731 3QA8_A   3QA8_A
    1732 4YHJ_A   4YHJ_A
    1733 3O17_A   3O17_A
    1734 2CPK_E   2CPK_E
    1735 3RZF_A   3RZF_A
    1736 4IAC_A   4IAC_A
    1737 8R99_A   8R99_A
    1738 4DG3_E   4DG3_E
    1739 3X2U_A   3X2U_A
    1740 4MYG_A   4MYG_A
    1741 9GLA_A   9GLA_A
    1742 2G01_A   2G01_A
    1743 6CCF_A   6CCF_A
    1744 8X5L_A   8X5L_A
    1745 2GNJ_A   2GNJ_A
    1746 3I6U_A   3I6U_A
    1747 6CMJ_A   6CMJ_A
    1748 1VZO_A   1VZO_A
    1749 3O7L_B   3O7L_B
    1750 2GU8_A   2GU8_A
    1751 1ZRZ_A   1ZRZ_A
    1752 3QD2_B   3QD2_B
    1753 3A8W_A   3A8W_A
    1754 8R3X_A   8R3X_A
    1755 3I6W_A   3I6W_A
    1756 3L9M_A   3L9M_A
    1757 5LI1_A   5LI1_A
    1758 2GNF_A   2GNF_A
    1759 4DC2_A   4DC2_A
    1760 2PK9_A   2PK9_A
    1761 8RU8_A   8RU8_A
    1762 5LI9_A   5LI9_A
    1763 1CMK_E   1CMK_E
    1764 1CTP_E   1CTP_E
    1765 5LIH_A   5LIH_A
    1766 6FRX_A   6FRX_A
    1767 1CDK_A   1CDK_A
    1768 3ZH8_A   3ZH8_A
    1769 3RP9_A   3RP9_A
    1770 8R9B_A   8R9B_A
    1771 4NTS_A   4NTS_A
    1772 9EJK_B   9EJK_B
    1774 2V7O_A   2V7O_A
    1775 4DFX_E   4DFX_E
    1776 4WNK_A   4WNK_A
    1777 8R3Y_I   8R3Y_I
    1778 4AE9_A   4AE9_A
    1779 1SYK_A   1SYK_A
    1780 10BL_A   10BL_A
    1781 4AE6_A   4AE6_A
    1782 10SL_A   10SL_A
    1783 6IB2_A   6IB2_A
    1784 9HIX_J   9HIX_J
    1785 5UV4_A   5UV4_A
    1786 2BFY_A   2BFY_A
    1787 5OO1_A   5OO1_A
    1788 4DFY_A   4DFY_A
    1789 8P7L_J   8P7L_J
    1790 1JBP_E   1JBP_E
    1791 8JFK_C   8JFK_C
    1792 1BKX_A   1BKX_A
    1793 1L3R_E   1L3R_E
    1794 5EYK_A   5EYK_A
    1795 3FI3_A   3FI3_A
    1796 1APM_E   1APM_E
    1797 5X3F_B   5X3F_B
    1798 8ORM_J   8ORM_J
    1799 5K3Y_A   5K3Y_A
    1800 2WEL_A   2WEL_A
    1801 6MM5_E   6MM5_E
    1802 7UJR_A   7UJR_A
    1803 8UKP_E   8UKP_E
    1804 3PVB_A   3PVB_A
    1805 7E11_A   7E11_A
    1806 8UKN_C   8UKN_C
    1807 1UA2_A   1UA2_A
    1808 6O9L_8   6O9L_8
    1809 5VLO_A   5VLO_A
    1810 1RDQ_E   1RDQ_E
    1811 1XH9_A   1XH9_A
    1812 3NX8_A   3NX8_A
    1813 4WB5_A   4WB5_A
    1814 3TTI_A   3TTI_A
    1815 3MVJ_A   3MVJ_A
    1816 6C0U_A   6C0U_A
    1817 8P4Z_A   8P4Z_A
    1818 3AGM_A   3AGM_A
    1819 3QAM_E   3QAM_E
    1820 4ERW_A   4ERW_A
    1821 1Q61_A   1Q61_A
    1822 4WHZ_A   4WHZ_A
    1823 3KVX_A   3KVX_A
    1824 7B55_B   7B55_B
    1825 4O21_A   4O21_A
    1826 9BLH_A   9BLH_A
    1827 2VN9_A   2VN9_A
    1828 1SZM_A   1SZM_A
    1829 1H4L_A   1H4L_A
    1830 4AU8_A   4AU8_A
    1831 2GNG_A   2GNG_A
    1832 2R9S_A   2R9S_A
    1833 2IW6_A   2IW6_A
    1834 4Z84_A   4Z84_A
    1835 3VUL_A   3VUL_A
    1836 2VZ6_A   2VZ6_A
    1837 7VDP_A   7VDP_A
    1838 8JPB_G   8JPB_G
    1839 7B5O_J   7B5O_J
    1840 3VUM_A   3VUM_A
    1841 4CFU_A   4CFU_A
    1842 3AMA_A   3AMA_A
    1843 3ELJ_A   3ELJ_A
    1844 3PTG_A   3PTG_A
    1845 6Q4G_A   6Q4G_A
    1846 1GZ8_A   1GZ8_A
    1847 3FV8_A   3FV8_A
    1848 7ORE_A   7ORE_A
    1849 1PMN_A   1PMN_A
    1850 3PZE_A   3PZE_A
    1851 5N23_A   5N23_A
    1852 3OXI_A   3OXI_A
    1853 1GII_A   1GII_A
    1854 1JNK_A   1JNK_A
    1855 4X21_A   4X21_A
    1856 4EOS_A   4EOS_A
    1857 4UX9_A   4UX9_A
    1858 1OGU_A   1OGU_A
    1859 3FI2_A   3FI2_A
    1860 4W4V_A   4W4V_A
    1861 2O0U_A   2O0U_A
    1862 2OK1_A   2OK1_A
    1863 4QTD_A   4QTD_A
    1864 9FT9_A   9FT9_A
    1865 4EOP_A   4EOP_A
    1866 4BCM_A   4BCM_A
    1867 4EOQ_A   4EOQ_A
    1868 1VYW_A   1VYW_A
    1869 2B1P_A   2B1P_A
    1870 4EOM_A   4EOM_A
    1871 3PXF_A   3PXF_A
    1872 1OIT_A   1OIT_A
    1873 4X3F_A   4X3F_A
    1874 9CKO_A   9CKO_A
    1875 5OO0_A   5OO0_A
    1876 2F7E_E   2F7E_E
    1877 3VUD_A   3VUD_A
    1878 2EXC_X   2EXC_X
    1879 6W4O_A   6W4O_A
    1880 2IW8_A   2IW8_A
    1881 4EOO_A   4EOO_A
    1882 6WJF_A   6WJF_A
    1883 5LW1_B   5LW1_B
    1884 4O38_A   4O38_A
    1885 4EOJ_A   4EOJ_A
    1886 4X3F_C   4X3F_C
    1887 8H6P_A   8H6P_A
    1888 4EON_A   4EON_A
    1889 9I9J_K   9I9J_K
    1890 7E34_A   7E34_A
    1891 1H1P_A   1H1P_A
    1892 1W98_A   1W98_A
    1893 4EOK_A   4EOK_A
    1894 8USO_A   8USO_A
    1895 8UV0_A   8UV0_A
    1896 7NVQ_A   7NVQ_A
    1897 9DC6_A   9DC6_A
    1898 9NFS_A   9NFS_A
    1899 3EZR_A   3EZR_A
    1900 6GUE_A   6GUE_A
    1901 4OW8_A   4OW8_A
    1902 8FEC_B   8FEC_B
    1903 3BHT_A   3BHT_A
    1904 4TNB_A   4TNB_A
    1905 6PJX_A   6PJX_A
    1906 4WB7_A   4WB7_A
    1907 9HIU_B   9HIU_B
    1908 4EOI_A   4EOI_A
    1909 6INL_A   6INL_A
    1910 2JGZ_A   2JGZ_A
    1911 5N3N_A   5N3N_A
    1912 4I3Z_A   4I3Z_A
    1913 1GY3_A   1GY3_A
    1914 5UQ1_A   5UQ1_A
    1915 4L7F_A   4L7F_A
    1916 8BZO_A   8BZO_A
    1917 1E9H_A   1E9H_A
    1918 9NYQ_A   9NYQ_A
    1919 2JDT_A   2JDT_A
    1920 6XD3_J   6XD3_J
    1921 5U6Y_A   5U6Y_A
    1922 5K4J_A   5K4J_A
    1923 3PJ8_A   3PJ8_A
    1924 1YDR_E   1YDR_E
    1925 6F14_A   6F14_A
    1926 2UVY_A   2UVY_A
    1927 6Q4I_A   6Q4I_A
    1928 1AQ1_A   1AQ1_A
    1929 1FQ1_B   1FQ1_B
    1930 1B38_A   1B38_A
    1931 6OQI_A   6OQI_A
    1932 1SMH_A   1SMH_A
    1933 4RT7_A   4RT7_A
    1935 6Y05_A   6Y05_A
    1936 9OB2_A   9OB2_A
    1937 9PE7_A   9PE7_A
    1938 8PYR_A   8PYR_A
    1939 6EM7_A   6EM7_A
    1940 1SVH_A   1SVH_A
    1941 5VI9_A   5VI9_A
    1942 3QHR_A   3QHR_A
    1943 3G2F_A   3G2F_A
    1944 4X3F_B   4X3F_B
    1945 3DND_A   3DND_A
    1946 5UUU_A   5UUU_A
    1947 1Q24_A   1Q24_A
    1948 1Q8W_A   1Q8W_A
    1949 3AGL_A   3AGL_A
    1950 2C1A_A   2C1A_A
    1951 6VZK_A   6VZK_A
    1952 3SV0_A   3SV0_A
    1953 3VUG_A   3VUG_A
    1954 2JDS_A   2JDS_A
    1955 3VUK_A   3VUK_A
    1956 6UNP_A   6UNP_A
    1957 8SF8_A   8SF8_A
    1958 8GXQ_h  8GXQ_HI
    1959 8YNG_A   8YNG_A
    1960 6XBZ_J   6XBZ_J
    1961 1BI7_A   1BI7_A
    1962 3VUI_A   3VUI_A
    1963 3VUH_A   3VUH_A
    1964 4C33_A   4C33_A
    1965 3KRW_A   3KRW_A
    1966 4WB8_A   4WB8_A
    1967 1JOW_B   1JOW_B
    1968 7PWD_A   7PWD_A
    1969 6TD3_B   6TD3_B
    1970 3CIK_A   3CIK_A
    1971 5NW8_A   5NW8_A
    1972 6C2Y_A   6C2Y_A
    1973 8BYA_A   8BYA_A
    1974 5DYK_A   5DYK_A
    1975 5N1F_A   5N1F_A
    1976 1OMW_A   1OMW_A
    1977 3PSC_A   3PSC_A
    1978 4C0T_A   4C0T_A
    1979 4WIH_A   4WIH_A
    1980 6B2Q_A   6B2Q_A
    1981 1UNG_A   1UNG_A
    1982 8I0M_A   8I0M_A
    1983 5UKK_A   5UKK_A
    1984 5HE1_A   5HE1_A
    1985 3NUP_A   3NUP_A
    1986 4MK0_A   4MK0_A
    1987 9I9I_K   9I9I_K
    1988 5MHQ_A   5MHQ_A
    1989 2UZT_A   2UZT_A
    1990 5HE3_A   5HE3_A
    1991 5HE0_A   5HE0_A
    1992 4IZ5_A   4IZ5_A
    1993 2VO0_A   2VO0_A
    1994 4WB6_B   4WB6_B
    1995 4WB6_A   4WB6_A
    1996 5EFQ_A   5EFQ_A
    1997 2CJM_A   2CJM_A
    1998 5YV8_A   5YV8_A
    1999 2ZV2_A   2ZV2_A
    2000 1Q8T_A   1Q8T_A
    2001 6RFP_A   6RFP_A
    2002 7VDU_A   7VDU_A
    2003 1XH7_A   1XH7_A
    2004 3ZO2_A   3ZO2_A
    2005 1H01_A   1H01_A
    2006 1OIR_A   1OIR_A
    2007 3VN9_A   3VN9_A
    2008 8S79_A   8S79_A
    2009 8A8M_B   8A8M_B
    2010 8XEY_A   8XEY_A
    2011 5UY6_A   5UY6_A
    2012 4C34_A   4C34_A
    2013 2AC5_A   2AC5_A
    2014 4XRL_A   4XRL_A
    2015 2QR8_A   2QR8_A
    2016 6OQL_A   6OQL_A
    2017 4JG6_A   4JG6_A
    2018 3RNY_A   3RNY_A
    2019 4BCF_A   4BCF_A
    2020 8WF4_A   8WF4_A
    2021 4NIF_A   4NIF_A
    2022 9HW6_B   9HW6_B
    2023 9IJJ_4  9IJJ_4Z
    2024 8I0L_A   8I0L_A
    2025 4IZA_A   4IZA_A
    2026 4S2Z_A   4S2Z_A
    2027 5O1S_A   5O1S_A
    2028 3KN5_A   3KN5_A
    2029 6OPG_A   6OPG_A
    2030 5WP1_A   5WP1_A
    2031 7UKZ_A   7UKZ_A
    2032 2BHH_A   2BHH_A
    2033 3NPC_A   3NPC_A
    2034 4EC8_A   4EC8_A
    2035 7N8T_A   7N8T_A
    2036 9HVX_A   9HVX_A
    2037 2WNT_A   2WNT_A
    2038 5V62_A   5V62_A
    2039 3E7O_A   3E7O_A
    2040 7OPM_A   7OPM_A
    2041 3BLH_A   3BLH_A
    2042 6U2G_A   6U2G_A
    2043 6OPK_A   6OPK_A
    2044 7MFD_B   7MFD_B
    2045 4RZ7_A   4RZ7_A
    2046 2QR7_A   2QR7_A
    2047 3ZUV_A   3ZUV_A
    2048 5KKR_C   5KKR_C
    2049 2ERK_A   2ERK_A
    2050 4QP1_A   4QP1_A
    2051 5EZR_A   5EZR_A
    2052 4IZ7_A   4IZ7_A
    2053 7CML_A   7CML_A
    2054 7W5O_A   7W5O_A
    2055 4XJ0_A   4XJ0_A
    2056 9HW6_A   9HW6_A
    2057 6G9J_A   6G9J_A
    2058 1PME_A   1PME_A
    2059 6GZH_A   6GZH_A
    2060 6G9K_A   6G9K_A
    2061 2GPH_A   2GPH_A
    2062 6W4O_O   6W4O_O
    2063 5LCK_A   5LCK_A
    2064 9Z8K_B   9Z8K_B
    2065 5K4I_A   5K4I_A
    2066 8ZJV_A   8ZJV_A
    2067 4QTA_A   4QTA_A
    2068 4XOY_A   4XOY_A
    2069 3FME_A   3FME_A
    2070 8RMB_A   8RMB_A
    2071 8PSR_A   8PSR_A
    2072 2Y9Q_A   2Y9Q_A
    2073 4QP4_A   4QP4_A
    2074 3O71_A   3O71_A
    2075 4FV6_A   4FV6_A
    2076 5AWM_A   5AWM_A
    2077 8CHF_E   8CHF_E
    2078 6W9E_A   6W9E_A
    2079 7XQK_A   7XQK_A
    2080 3MI9_A   3MI9_A
    2081 2OJG_A   2OJG_A
    2082 3C9W_A   3C9W_A
    2083 5L1Z_A   5L1Z_A
    2084 4OR5_A   4OR5_A
    2085 6OTS_A   6OTS_A
    2086 9AXA_B   9AXA_B
    2087 1WZY_A   1WZY_A
    2088 1TVO_A   1TVO_A
    2089 6RFO_A   6RFO_A
    2090 10JU_B   10JU_B
    2091 2FYS_A   2FYS_A
    2092 6OT6_A   6OT6_A
    2093 3ZU7_A   3ZU7_A
    2094 4QYY_A   4QYY_A
    2095 2Z7L_A   2Z7L_A
    2096 6NBS_A   6NBS_A
    2097 4IMY_A   4IMY_A
    2098 4ZZM_A   4ZZM_A
    2099 3R63_A   3R63_A
    2100 9SKQ_A   9SKQ_A
    2101 4XOZ_A   4XOZ_A
    2102 9O0U_B   9O0U_B
    2103 8BW9_C   8BW9_C
    2104 10JU_A   10JU_A
    2105 6PP9_B   6PP9_B
    2106 4FUX_A   4FUX_A
    2107 7UGB_A   7UGB_A
    2108 9TYG_A   9TYG_A
    2109 8ELC_A   8ELC_A
    2110 5BUE_A   5BUE_A
    2111 9QQJ_A   9QQJ_A
    2112 3SOA_A   3SOA_A
    2113 3QYW_A   3QYW_A
    2114 4Y72_A   4Y72_A
    2115 6DMG_A   6DMG_A
    2116 9Z8K_A   9Z8K_A
    2117 7NJ0_B   7NJ0_B
    2118 6NYB_B   6NYB_B
    2119 5LCJ_A   5LCJ_A
    2120 6Q0T_C   6Q0T_C
    2121 8RLX_A   8RLX_A
    2122 9TYG_B   9TYG_B
    2123 4AN2_A   4AN2_A
    2124 4YC6_A   4YC6_A
    2125 6FI6_A   6FI6_A
    2126 4XNE_A   4XNE_A
    2127 4XP2_A   4XP2_A
    2128 8AOC_A   8AOC_A
    2129 6Z45_A   6Z45_A
    2130 4NST_A   4NST_A
    2131 4O6E_A   4O6E_A
    2132 4FV7_A   4FV7_A
    2133 4GSB_A   4GSB_A
    2134 6FXV_A   6FXV_A
    2135 3SA0_A   3SA0_A
    2136 2ZOQ_A   2ZOQ_A
    2137 5NHH_A   5NHH_A
    2138 1GOL_A   1GOL_A
    2139 6GDM_A   6GDM_A
    2140 2B9H_A   2B9H_A
    2141 8RM2_A   8RM2_A
    2142 5NGU_A   5NGU_A
    2143 8K5R_A   8K5R_A
    2144 4N0S_A   4N0S_A
    2145 4JG8_A   4JG8_A
    2146 4S30_A   4S30_A
    2147 3TEI_A   3TEI_A
    2148 9TU0_A   9TU0_A
    2149 4QTB_A   4QTB_A
    2150 4I5H_A   4I5H_A
    2151 9JK1_A   9JK1_A
    2152 2Y4I_C   2Y4I_C
    2153 4CXA_A   4CXA_A
    2154 4UN0_C   4UN0_C
    2155 2AC3_A   2AC3_A
    2156 8XFM_A   8XFM_A
    2157 6B3E_A   6B3E_A
    2158 7F2X_A   7F2X_A
    2159 2B9F_A   2B9F_A
    2160 7E73_A   7E73_A
    2161 6XI8_A   6XI8_A
    2162 2F9G_A   2F9G_A
    2163 4H3Q_A   4H3Q_A
    2164 7KUE_A   7KUE_A
    2165 4CRS_A   4CRS_A
    2166 3ORN_A   3ORN_A
    2167 4XHL_A   4XHL_A
    2168 3N9X_A   3N9X_A
    2169 7E75_A   7E75_A
    2170 2P55_A   2P55_A
    2171 5WVD_A   5WVD_A
    2172 7W5C_A   7W5C_A
    2173 4U7Z_A   4U7Z_A
    2174 1S9J_A   1S9J_A
    2175 5EYM_A   5EYM_A
    2176 2HW6_A   2HW6_A
    2177 3EQC_A   3EQC_A
    2178 5YXI_A   5YXI_A
    2179 7JUQ_C   7JUQ_C
    2180 8YP4_A   8YP4_A
    2181 3ZLY_A   3ZLY_A
    2182 3ZLS_A   3ZLS_A
    2183 4AW2_A   4AW2_A
    2184 3DV3_A   3DV3_A
    2185 5HZE_A   5HZE_A
    2186 3MBL_A   3MBL_A
    2187 3W8Q_A   3W8Q_A
    2188 7PQV_A   7PQV_A
    2189 6PXN_A   6PXN_A
    2190 5YT3_A   5YT3_A
    2191 6Z3U_B   6Z3U_B
    2192 6BDL_A   6BDL_A
    2193 1S9I_A   1S9I_A
    2194 5V5Y_A   5V5Y_A
    2195 7N3U_A   7N3U_A
    2196 8WDK_W   8WDK_W
    2197 2IN6_A   2IN6_A
    2198 7LV3_A   7LV3_A
    2199 1X8B_A   1X8B_A
    2200 3BI6_A   3BI6_A
    2201 7NAA_A   7NAA_A
    2202 2Z2W_A   2Z2W_A
    2203 4BGQ_A   4BGQ_A
    2204 8H59_A   8H59_A
    2205 7T4T_A   7T4T_A
    2206 6DTL_A   6DTL_A
    2207 3ENM_A   3ENM_A
    2208 5Z33_A   5Z33_A
    2209 4OTD_A   4OTD_A
    2210 9AXH_A   9AXH_A
    2211 5CZO_A   5CZO_A
    2212 4UAK_A   4UAK_A
    2213 8ZTC_A   8ZTC_A
    2214 3QFV_A   3QFV_A
    2215 7JV7_A   7JV7_A
    2216 5CI6_A   5CI6_A
    2217 3SLS_A   3SLS_A
    2218 6P5M_A   6P5M_A
    2219 3TKU_A   3TKU_A
    2220 7JNT_A   7JNT_A
    2221 7B3M_A   7B3M_A
    2222 5U7Q_A   5U7Q_A
    2223 4WOT_A   4WOT_A
    2224 4AAA_A   4AAA_A
    2225 5U7R_A   5U7R_A
    2226 3NIE_A   3NIE_A
    2227 8X8X_A   8X8X_A
    2228 9AXM_A   9AXM_A
    2229 6ED6_A   6ED6_A
    2230 3ZLW_A   3ZLW_A
    2231 4L6Q_A   4L6Q_A
    2232 6PXP_A   6PXP_A
    2233 2F2U_A   2F2U_A
    2234 4KB8_A   4KB8_A
    2235 4XH0_A   4XH0_A
    2236 6RUU_A   6RUU_A
    2237 2H34_A   2H34_A
    2238 6ZIW_I   6ZIW_I
    2239 6E9W_A   6E9W_A
    2240 9JCU_A   9JCU_A
    2241 3TV7_A   3TV7_A
    2242 4W7P_A   4W7P_A
    2243 4QNY_A   4QNY_A
    2244 7JOU_A   7JOU_A
    2245 2ESM_A   2ESM_A
    2246 8ZH5_A   8ZH5_A
    2247 2V55_A   2V55_A
    2248 7S26_A   7S26_A
    2249 9B3S_A   9B3S_A
    2250 7S25_A   7S25_A
    2251 4TN6_A   4TN6_A
    2252 8VXF_A   8VXF_A
    2253 8D7M_A   8D7M_A
    2254 5MQV_A   5MQV_A
    2255 3UYS_A   3UYS_A
    2256 4JJR_A   4JJR_A
    2257 8VXD_A   8VXD_A
    2258 5OKT_A   5OKT_A
    2259 5IH4_A   5IH4_A
    2260 5X17_A   5X17_A
    2261 8IZC_A   8IZC_A
    2262 7P7F_A   7P7F_A
    2263 6RCG_A   6RCG_A
    2264 1CKI_A   1CKI_A
    2265 4TW9_A   4TW9_A
    2266 4HNI_A   4HNI_A
    2267 5CYZ_A   5CYZ_A
    1658 8T7T_A   8T7T_A
    2268 9FQR_x  9FQR_Xr
    2269 7UP4_A   7UP4_A
    2270 2ELI_A   2ELI_A
    2271 5L2Q_A   5L2Q_A
    2272 2VD5_A   2VD5_A
    2273 3OZ6_A   3OZ6_A
    2274 7WTT_a   7WTT_a
    2275 6GZD_A   6GZD_A
    2276 4FI1_A   4FI1_A
    2277 8TQ2_A   8TQ2_A
    2278 4O2Z_A   4O2Z_A
    2279 9H8C_A   9H8C_A
    2280 5X18_A   5X18_A
    2281 4CRL_A   4CRL_A
    2282 3RGF_A   3RGF_A
    2283 6T41_A   6T41_A
    2284 6TPA_A   6TPA_A
    2285 6QTG_A   6QTG_A
    2286 9R59_A   9R59_A
    2287 8XU4_A   8XU4_A
    2288 6BXI_A   6BXI_A
    2289 5IDN_A   5IDN_A
    2290 5HNB_A   5HNB_A
    2291 5M4U_A   5M4U_A
    2292 5XQX_A   5XQX_A
    2293 5FQD_C   5FQD_C
    2294 9IHG_A   9IHG_A
    2295 3OFM_A   3OFM_A
    2296 5OOI_A   5OOI_A
    2297 5FGK_A   5FGK_A
    2298 7KPV_A   7KPV_A
    2299 6HMD_A   6HMD_A
    2300 9HK5_A   9HK5_A
    2301 3KA0_A   3KA0_A
    2302 6L20_A   6L20_A
    2303 6QY8_A   6QY8_A
    2304 9OTY_C   9OTY_C
    2305 3E3B_X   3E3B_X
    2306 6T8X_A   6T8X_A
    2307 3UFF_A   3UFF_A
    2308 7PSU_A   7PSU_A
    2309 3UGD_A   3UGD_A
    2310 2JBO_A   2JBO_A
    2311 3GOK_A   3GOK_A
    2312 6TLL_A   6TLL_A
    2313 2PZY_A   2PZY_A
    2314 4MD8_E   4MD8_E
    2315 6TCA_A   6TCA_A
    2316 3R2B_A   3R2B_A
    2317 4MD7_E   4MD7_E
    2318 5OMY_A   5OMY_A
    2319 4TYH_A   4TYH_A
    2320 2OZA_A   2OZA_A
    2321 3FPM_A   3FPM_A
    2322 4FKD_A   4FKD_A
    2323 3R2Y_A   3R2Y_A
    2324 2P3G_X   2P3G_X
    2325 1KWP_A   1KWP_A
    2326 3UEJ_A   3UEJ_A
    2327 2ONL_C   2ONL_C
    2328 1PTQ_A   1PTQ_A
    2329 7KND_A   7KND_A
    2330 6O6Q_A   6O6Q_A
    2331 2ENZ_A   2ENZ_A
    2332 3UEY_A   3UEY_A
    2333 6L22_A   6L22_A
    2334 3U87_A   3U87_A
    2335 6UWA_A   6UWA_A
    2336 1NA7_A   1NA7_A
    2337 6L24_A   6L24_A
    2338 5XVU_A   5XVU_A
    2339 2CSN_A   2CSN_A
    2340 1CSN_A   1CSN_A
    2341 2ZJW_A   2ZJW_A
    2342 3JUH_A   3JUH_A
    2343 1NXK_A   1NXK_A
    2344 6L23_A   6L23_A
    2345 5CSP_A   5CSP_A
    2346 5OSL_A   5OSL_A
    2347 6QY7_A   6QY7_A
    2348 5MOV_A   5MOV_A
    2349 9HXU_A   9HXU_A
    2350 4MD9_E   4MD9_E
    2351 1JWH_A   1JWH_A
    2352 3Q9W_A   3Q9W_A
    2353 6HME_A   6HME_A
    2354 5N1V_A   5N1V_A
    2355 7A4Q_A   7A4Q_A
    2356 1XA6_A   1XA6_A
    2357 3W8L_A   3W8L_A
    2358 5KU8_A   5KU8_A
    2359 1PJK_A   1PJK_A
    2360 2R7I_A   2R7I_A
    2361 3BQC_A   3BQC_A
    2362 6Z83_a 6Z83_AAA
    2363 3NSZ_A   3NSZ_A
    2364 3H30_A   3H30_A
    2365 3NGA_A   3NGA_A
    2366 6Z19_B   6Z19_B
    2367 3MB6_A   3MB6_A
    2368 5ZN5_A   5ZN5_A
    2369 3Q04_A   3Q04_A
    2370 6Q38_A   6Q38_A
    2371 6YPH_A   6YPH_A
    2372 1TBN_A   1TBN_A
    2373 5CVG_A   5CVG_A
    2374 5CLP_A   5CLP_A
    2375 7QUX_A   7QUX_A
    2376 6L21_A   6L21_A
    2377 6YZH_A   6YZH_A
    2378 5ZN0_A   5ZN0_A
    2379 3CXL_A   3CXL_A
    2380 9EDY_A   9EDY_A
    2381 6U69_A   6U69_A
    2382 4DGL_C   4DGL_C
    2383 6JKK_A   6JKK_A
    2384 1KBE_A   1KBE_A
    1375 8DFP_A   8DFP_A
    1380 8DFQ_A   8DFQ_A
    1382 8DFM_A   8DFM_A
    2385 5M07_A   5M07_A
    2386 5M06_A   5M06_A
    2387 2PVH_A   2PVH_A
    2388 4JRN_A   4JRN_A
    2389 5XKA_A   5XKA_A
    2390 4ANM_A   4ANM_A
    2391 1DS5_A   1DS5_A
    2392 3PVG_A   3PVG_A
    2393 1M2P_A   1M2P_A
    2394 1DAW_A   1DAW_A
    2395 2QC6_A   2QC6_A
    2396 4DGN_A   4DGN_A
    2397 3KXG_A   3KXG_A
    2398 5TS8_A   5TS8_A
    2399 4DGM_A   4DGM_A
    781  3PFQ_A   3PFQ_A
    2400 1Y8F_A   1Y8F_A
    771  8SE1_A   8SE1_A
    773  8SE2_A   8SE2_A
    2401 6RA0_A   6RA0_A
    2402 3KGA_A   3KGA_A
    1934 4RT7_A   4RT7_A
    2403 7P5Z_1   7P5Z_1
    2404 2YUU_A   2YUU_A
    2405 7T7C_A   7T7C_A
    2406 5UE8_A   5UE8_A
    2407 6VP8_B   6VP8_B
    2408 7DG2_A   7DG2_A
    2409 3UIB_A   3UIB_A
    2410 3PG1_A   3PG1_A
    2411 2ENN_A   2ENN_A
    2412 2E73_A   2E73_A
    2413 2DB6_A   2DB6_A
    2414 9C5F_A   9C5F_A
    2415 5HNV_A   5HNV_A
    2416 4IX3_A   4IX3_A
    1773 9EJK_B   9EJK_B
    2417 4B6D_A   4B6D_A

    $raw
               queryid subjectids identity alignmentlength mismatches gapopens
    1    Query_1615341     6Q0K_A   99.869             766          1        0
    2    Query_1615341     7MFD_A   99.739             766          2        0
    3    Query_1615341     6Q0J_A   99.739             766          2        0
    4    Query_1615341     6NYB_A   99.739             766          2        0
    5    Query_1615341     8VYO_A   99.739             766          2        0
    6    Query_1615341     7ZR0_K   99.869             766          1        0
    7    Query_1615341     6UAN_B   99.608             766          3        0
    8    Query_1615341     8VYV_C   99.608             766          3        0
    9    Query_1615341     8VYU_A   99.608             766          3        0
    10   Query_1615341     8VYP_C   99.608             766          3        0
    11   Query_1615341     8VYS_A   99.478             766          4        0
    12   Query_1615341     8VYW_C   99.200             375          3        0
    13   Query_1615341     8VYQ_C   99.200             375          3        0
    14   Query_1615341     9MMS_A   57.058             673        245       13
    15   Query_1615341     9MMR_A   56.909             673        246       13
    16   Query_1615341     9MMQ_A   56.761             673        247       13
    17   Query_1615341     8CHF_A   60.359             613        206       12
    18   Query_1615341     8U1L_C   56.464             673        249       13
    19   Query_1615341     9MMP_A   56.464             673        249       13
    20   Query_1615341   7Z37_CP1   59.547             618        213       12
    21   Query_1615341     4MNF_A   99.672             305          1        0
    22   Query_1615341     4MNE_B  100.000             295          0        0
    23   Query_1615341     3II5_A  100.000             295          0        0
    24   Query_1615341     3D4Q_A  100.000             295          0        0
    25   Query_1615341     7SHV_A  100.000             295          0        0
    26   Query_1615341     5FD2_A  100.000             294          0        0
    27   Query_1615341     4EHG_A   99.661             295          1        0
    28   Query_1615341     3IDP_A   99.322             295          2        0
    29   Query_1615341     4MBJ_A  100.000             292          0        0
    30   Query_1615341     5HIE_A   98.305             295          0        1
    31   Query_1615341     4DBN_A  100.000             282          0        0
    32   Query_1615341     3Q96_A  100.000             282          0        0
    33   Query_1615341     6PP9_A   99.645             282          1        0
    34   Query_1615341     9ECU_A  100.000             279          0        0
    35   Query_1615341     2FB8_A  100.000             279          0        0
    36   Query_1615341     9AXX_B  100.000             279          0        0
    37   Query_1615341     4H58_A  100.000             275          0        0
    38   Query_1615341     1UWH_A   99.638             276          1        0
    39   Query_1615341     1UWJ_A   99.638             276          1        0
    40   Query_1615341     9BFB_A   98.925             279          3        0
    41   Query_1615341     6N0P_A  100.000             273          0        0
    42   Query_1615341     6U2H_C   94.810             289         15        0
    43   Query_1615341     4YHT_A   98.162             272          5        0
    44   Query_1615341     4XV9_A   94.643             280         15        0
    45   Query_1615341     4JVG_A   94.964             278         14        0
    46   Query_1615341     5CSW_A   94.643             280         15        0
    47   Query_1615341     3C4C_A   94.964             278         14        0
    48   Query_1615341     4RZV_A   94.624             279         15        0
    49   Query_1615341     4WO5_A   94.964             278         14        0
    50   Query_1615341     4FK3_A   94.286             280         16        0
    51   Query_1615341     8F7O_A   93.929             280         17        0
    52   Query_1615341     6XFP_A   94.286             280         16        0
    53   Query_1615341     8C7Y_A   94.604             278         15        0
    54   Query_1615341     4XV1_A   93.571             280         18        0
    55   Query_1615341     6V34_A   93.571             280         18        0
    56   Query_1615341     3OG7_A   93.525             278         18        0
    57   Query_1615341     5ITA_A   93.548             279         18        0
    58   Query_1615341     5JRQ_A   94.526             274         15        0
    59   Query_1615341     6P3D_A   93.190             279         19        0
    60   Query_1615341     4CQE_A   94.526             274         15        0
    61   Query_1615341     7P3V_A   94.853             272         14        0
    62   Query_1615341     5HI2_A   92.857             280         15        1
    63   Query_1615341     9AXC_A   74.854             342         79        4
    64   Query_1615341     9AXA_A   74.561             342         80        4
    65   Query_1615341     3OMV_A   77.627             295         66        0
    66   Query_1615341     9O0U_A   81.004             279         53        0
    67   Query_1615341     9AY7_A   81.004             279         53        0
    68   Query_1615341     8GFT_D   80.000             285         57        0
    69   Query_1615341     9AXM_B   76.703             279         65        0
    70   Query_1615341     8GAE_D   83.495             206         34        0
    71   Query_1615341     6PTS_D   62.308             130         46        1
    72   Query_1615341     7JHP_C   62.308             130         46        1
    73   Query_1615341     6XI7_B   62.308             130         46        1
    74   Query_1615341     2Y4I_B   36.519             293        167        8
    75   Query_1615341     5KKR_B   36.519             293        167        8
    76   Query_1615341     8BW9_D   38.644             295        165        9
    77   Query_1615341     2L05_A  100.000              84          0        0
    78   Query_1615341     7JUQ_B   36.519             293        167        8
    79   Query_1615341     3NY5_A   97.647              85          2        0
    80   Query_1615341     6XGU_B   61.538             130         47        1
    81   Query_1615341     5J17_A   88.542              96          2        1
    82   Query_1615341     7JUW_B   35.254             295        170        8
    83   Query_1615341     5VR3_A   88.298              94          5        1
    84   Query_1615341     9AXH_C   34.982             283        167        7
    85   Query_1615341     3PPZ_A   36.934             287        165       10
    86   Query_1615341     4CSV_A   36.559             279        153       10
    87   Query_1615341     3P86_A   36.585             287        166       10
    88   Query_1615341     5VYK_A  100.000              73          0        0
    89   Query_1615341     8DEG_A   37.729             273        148        8
    90   Query_1615341     5CEN_A   37.729             273        148        8
    91   Query_1615341     9NS1_A   33.574             277        162        9
    92   Query_1615341     8JF3_A   33.574             277        162        9
    93   Query_1615341     2BDF_A   33.574             277        162        9
    94   Query_1615341     7NG7_A   33.574             277        162        9
    95   Query_1615341     7OTE_A   33.574             277        162        9
    96   Query_1615341     4MXO_A   33.574             277        162        9
    97   Query_1615341     6E6E_A   33.818             275        160        9
    98   Query_1615341     8HAQ_A   33.818             275        160        9
    99   Query_1615341     3A4O_X   35.075             268        151       10
    100  Query_1615341     4MXX_A   33.574             277        162        9
    101  Query_1615341     1YOL_A   33.213             277        163        9
    102  Query_1615341     9NS0_A   33.574             277        162        9
    103  Query_1615341     3OEZ_A   33.574             277        162        9
    104  Query_1615341     2OIQ_A   33.574             277        162        9
    105  Query_1615341     1YI6_A   33.818             275        160        9
    106  Query_1615341     2OG8_A   32.852             277        160       10
    107  Query_1615341     1YOJ_A   32.971             276        165        9
    108  Query_1615341     3D7U_B   33.696             276        161        9
    109  Query_1615341     4MXY_A   33.574             277        162        9
    110  Query_1615341     2PL0_A   32.734             278        161       10
    111  Query_1615341     2OFV_A   32.734             278        161       10
    112  Query_1615341     3BYS_A   32.734             278        161       10
    113  Query_1615341     5XY1_A   34.962             266        150       10
    114  Query_1615341     2OF2_A   32.734             278        161       10
    115  Query_1615341     2DQ7_X   34.058             276        160        9
    116  Query_1615341     3GEQ_A   33.574             277        162        9
    117  Query_1615341     3KXZ_A   32.727             275        165        9
    118  Query_1615341     5T0P_A   33.213             277        163        9
    119  Query_1615341     5SWH_A   33.213             277        163        9
    120  Query_1615341     2HWO_A   33.213             277        163        9
    121  Query_1615341     2ZM1_A   32.727             275        165        9
    122  Query_1615341     1QPC_A   32.727             275        165        9
    123  Query_1615341     3SVV_A   33.213             277        163        9
    124  Query_1615341     2OFU_A   32.727             275        165        9
    125  Query_1615341     3KMM_A   32.727             275        165        9
    126  Query_1615341     6PDJ_A   32.727             275        165        9
    127  Query_1615341     3U4W_A   33.700             273        159        9
    128  Query_1615341     4MCV_A   33.213             277        163        9
    129  Query_1615341     3G6H_A   33.213             277        163        9
    130  Query_1615341     3BYO_A   32.727             275        165        9
    131  Query_1615341     3BYM_A   32.727             275        165        9
    132  Query_1615341     4LGH_A   33.333             276        162        9
    133  Query_1615341     2H8H_A   33.574             277        162        9
    134  Query_1615341     8JN8_A   33.574             277        162        9
    135  Query_1615341     9IRL_A   33.574             277        162        9
    136  Query_1615341     2ZV7_A   34.586             266        151       10
    137  Query_1615341     1FMK_A   33.574             277        162        9
    138  Query_1615341     2HK5_A   33.333             270        158        9
    139  Query_1615341     5ZJ6_A   33.333             270        158        9
    140  Query_1615341     1Y57_A   33.574             277        162        9
    141  Query_1615341     2QI8_A   32.852             277        164        9
    142  Query_1615341     9BT8_C   33.935             277        161        9
    143  Query_1615341     3DQW_A   33.091             275        166        9
    144  Query_1615341     3MPM_A   32.246             276        161       10
    145  Query_1615341     8XN8_A   33.213             277        163        9
    146  Query_1615341     6F3F_A   33.213             277        163        9
    147  Query_1615341     4K11_A   33.213             277        163        9
    148  Query_1615341     1KSW_A   33.213             277        163        9
    149  Query_1615341     1QPD_A   32.364             275        166        9
    150  Query_1615341     2PTK_A   33.213             277        163        9
    151  Query_1615341     6IN0_A   33.559             295        173       10
    152  Query_1615341     2QOK_A   33.898             295        172       10
    153  Query_1615341     2QOC_A   33.784             296        171       11
    154  Query_1615341     3DZQ_A   33.559             295        173       10
    155  Query_1615341     1AD5_A   32.721             272        157        8
    156  Query_1615341     2QOO_A   33.559             295        173       10
    157  Query_1615341     2QOI_A   33.559             295        173       10
    158  Query_1615341     2QOF_A   33.559             295        173       10
    159  Query_1615341     2QOD_A   33.559             295        173       10
    160  Query_1615341     3FXX_A   33.559             295        173       10
    161  Query_1615341     2GSF_A   33.559             295        173       10
    162  Query_1615341     1QCF_A   32.727             275        163        9
    163  Query_1615341     2QOL_A   33.559             295        173       10
    164  Query_1615341     9BYJ_A   32.727             275        163        9
    165  Query_1615341     2YN8_A   30.928             291        181        9
    166  Query_1615341     2QOB_A   33.784             296        171       11
    167  Query_1615341     6FNI_A   30.928             291        181        9
    168  Query_1615341     2QON_A   33.559             295        173       10
    169  Query_1615341     2QO7_A   33.559             295        173       10
    170  Query_1615341     6CZ2_A   33.582             268        163        8
    171  Query_1615341     2XYU_A   32.534             292        177        9
    172  Query_1615341     2Y6M_A   32.534             292        177        9
    173  Query_1615341     6CZ3_A   33.582             268        163        8
    174  Query_1615341     2HEL_A   32.534             292        177        9
    175  Query_1615341     3KUL_A   32.997             297        181        9
    176  Query_1615341     2VWU_A   30.605             281        175        9
    177  Query_1615341     8JNA_B   71.605              81         23        0
    178  Query_1615341     3KUL_B   32.886             298        180        9
    179  Query_1615341     3ZEW_A   30.241             291        183        9
    180  Query_1615341     4AW5_A   29.897             291        184        8
    181  Query_1615341     5MJA_A   31.673             281        172        8
    182  Query_1615341     5MJB_A   31.673             281        172        8
    183  Query_1615341     7UY0_B   33.803             284        164       10
    184  Query_1615341     7UY0_A   33.803             284        164       10
    185  Query_1615341     7KPL_A   31.673             281        172        8
    186  Query_1615341     4LGG_A   33.333             261        152        9
    187  Query_1615341     3ZFX_A   31.673             281        172        8
    188  Query_1615341     6UMW_A   31.915             282        170       10
    189  Query_1615341     2R2P_A   32.975             279        169        8
    190  Query_1615341     3ZFY_A   31.973             294        178       10
    191  Query_1615341     3ZFM_A   31.317             281        173        7
    192  Query_1615341     1JPA_A   31.317             281        173        7
    193  Query_1615341     5D7V_A   32.955             264        162        8
    194  Query_1615341     5DA3_A   32.955             264        162        8
    195  Query_1615341     1K3A_A   31.488             289        166        9
    196  Query_1615341     5H2U_A   32.955             264        162        8
    197  Query_1615341     2ZM3_A   31.488             289        166        9
    198  Query_1615341     5I9U_A   31.429             280        173        8
    199  Query_1615341     8XPV_A   31.429             280        173        8
    200  Query_1615341     5FXS_A   31.741             293        160       11
    201  Query_1615341     3O23_A   31.741             293        160       11
    202  Query_1615341     1JQH_A   31.741             293        160       11
    203  Query_1615341     3QQU_A   31.741             293        160       11
    204  Query_1615341     1MQB_A   31.429             280        173        8
    205  Query_1615341     1M7N_A   31.741             293        160       11
    206  Query_1615341     2OJ9_A   31.741             293        160       11
    207  Query_1615341     4TRL_A   32.028             281        168        9
    208  Query_1615341     4D2R_A   31.741             293        160       11
    209  Query_1615341     5FXQ_A   31.741             293        160       11
    210  Query_1615341     3I81_A   31.741             293        160       11
    211  Query_1615341   8PYI_AAA   31.741             293        160       11
    212  Query_1615341     5FXR_A   31.741             293        160       11
    213  Query_1615341     3D94_A   31.741             293        160       11
    214  Query_1615341     4P2K_A   31.541             279        172        8
    215  Query_1615341     5EK7_A   31.541             279        172        8
    216  Query_1615341     8XKP_A   31.525             295        178       10
    217  Query_1615341     7KJA_A   31.429             280        173        8
    218  Query_1615341     6JMF_A   31.525             295        178       10
    219  Query_1615341     7KJC_A   31.429             280        173        8
    220  Query_1615341     7KJB_A   31.429             280        173        8
    221  Query_1615341     3LVP_A   31.741             293        160       11
    222  Query_1615341     3BKB_A   31.525             295        178       10
    223  Query_1615341     7EEF_A   31.752             274        167        9
    224  Query_1615341     5X5O_A   31.618             272        155        9
    225  Query_1615341     3LW0_A   31.399             293        161       11
    226  Query_1615341     5HES_A   31.618             272        155        9
    227  Query_1615341     1P4O_A   31.399             293        161       11
    228  Query_1615341     7EEC_A   31.429             280        172        9
    229  Query_1615341     2REI_A   31.752             274        167        9
    230  Query_1615341     7EED_A   31.429             280        172        9
    231  Query_1615341     4XLV_A   30.492             305        179       10
    232  Query_1615341     1RQQ_A   30.492             305        179       10
    233  Query_1615341     3CD3_A   31.058             293        182        8
    234  Query_1615341   8ATL_BBB   30.712             267        164        7
    235  Query_1615341     4UYA_A   32.632             285        142       14
    236  Query_1615341   8ATL_AAA   30.712             267        164        7
    237  Query_1615341     2Z8C_A   30.877             285        165        9
    238  Query_1615341   8ATB_AAA   30.712             267        164        7
    239  Query_1615341     8FLN_A   28.938             273        170        9
    240  Query_1615341     1GAG_A   30.877             285        165        9
    241  Query_1615341     7TYJ_A   31.724             290        164       10
    242  Query_1615341     6J6M_A   28.623             276        173        9
    243  Query_1615341     6MNY_A   28.938             273        170        9
    244  Query_1615341     3DK3_A   30.605             281        171       10
    245  Query_1615341     1I44_A   30.968             310        171       12
    246  Query_1615341     6VXQ_A   28.713             303        189       11
    247  Query_1615341     4IBM_A   30.744             309        173       12
    248  Query_1615341     6BL8_A   30.466             279        170       10
    249  Query_1615341     8DWN_A   30.492             305        179       10
    250  Query_1615341     5HU9_A   29.720             286        177       10
    251  Query_1615341     4Z3V_A   28.571             273        171        9
    252  Query_1615341     6NZM_A   28.571             273        171        9
    253  Query_1615341     4XEY_A   29.932             294        178       11
    254  Query_1615341     6W7O_A   28.571             273        171        9
    255  Query_1615341     4YHF_A   28.571             273        171        9
    256  Query_1615341     3GEN_A   28.571             273        171        9
    257  Query_1615341     8YVV_A   28.571             273        171        9
    258  Query_1615341     5J87_A   28.571             273        171        9
    259  Query_1615341     3P08_A   28.571             273        171        9
    260  Query_1615341     6TFP_A   28.571             273        171        9
    261  Query_1615341     5P9F_A   28.571             273        171        9
    262  Query_1615341     4OT5_A   28.767             292        182       10
    263  Query_1615341     3OCT_A   28.571             273        171        9
    264  Query_1615341     6S90_A   28.571             273        171        9
    265  Query_1615341     3PIX_A   28.571             273        171        9
    266  Query_1615341     4WA9_A   29.825             285        176       10
    267  Query_1615341     1IRK_A   30.744             309        173       12
    268  Query_1615341     6AUB_A   28.571             273        171        9
    269  Query_1615341     6O8I_A   28.571             273        171        9
    270  Query_1615341     5XYZ_A   28.571             273        171        9
    271  Query_1615341     6NFI_A   28.571             273        171        9
    272  Query_1615341     7KXQ_A   28.571             273        171        9
    273  Query_1615341     8GC8_A   28.571             273        171        9
    274  Query_1615341     5ZZ4_A   28.571             273        171        9
    275  Query_1615341     2XYN_A   30.961             281        170       10
    276  Query_1615341     1FPU_A   30.108             279        171       10
    277  Query_1615341     5FBN_C   28.571             273        171        9
    278  Query_1615341     4XLI_A   30.824             279        169       10
    279  Query_1615341     4ZLY_A   28.571             273        171        9
    280  Query_1615341     8FLL_A   28.571             273        171        9
    281  Query_1615341     9ZLJ_A   28.938             273        170        9
    282  Query_1615341     6XRG_A   30.466             279        170       10
    283  Query_1615341     2F4J_A   29.825             285        176       10
    284  Query_1615341     3OCS_A   28.571             273        171        9
    285  Query_1615341     3QRI_A   30.108             279        171       10
    286  Query_1615341     8SSN_A   29.932             294        178       11
    287  Query_1615341     2HZI_A   30.108             279        171       10
    288  Query_1615341     2NRY_A   30.712             267        164        7
    289  Query_1615341     4RX5_A   28.571             273        171        9
    290  Query_1615341     2E2B_A   30.108             279        171       10
    291  Query_1615341     4OTF_A   28.571             273        171        9
    292  Query_1615341     3OXZ_A   30.108             279        171       10
    293  Query_1615341     6O8U_A   30.712             267        164        7
    294  Query_1615341     2QOH_A   30.108             279        171       10
    295  Query_1615341     5K72_A   30.712             267        164        7
    296  Query_1615341     6EG9_A   30.712             267        164        7
    297  Query_1615341     6XE4_A   28.571             273        171        9
    298  Query_1615341     6XR6_A   30.108             279        171       10
    299  Query_1615341     4ZOG_A   30.108             279        171       10
    300  Query_1615341     4Y95_A   27.839             273        173        9
    301  Query_1615341     2NRU_A   30.712             267        164        7
    302  Query_1615341     6THW_A   30.712             267        164        7
    303  Query_1615341     2HIW_A   30.108             279        171       10
    304  Query_1615341     7N9G_A   30.108             279        171       10
    305  Query_1615341     6THX_A   30.712             267        164        7
    306  Query_1615341     6E4F_A   28.571             273        171        9
    307  Query_1615341     2HZ0_A   30.108             279        171       10
    308  Query_1615341     9PSU_A   30.712             267        164        7
    309  Query_1615341     6LXY_A   30.712             267        164        7
    310  Query_1615341     2HYY_A   30.108             279        171       10
    311  Query_1615341     2G1T_A   30.108             279        171       10
    312  Query_1615341     6BIK_A   28.571             273        171        9
    313  Query_1615341     1P14_A   30.421             309        174       12
    314  Query_1615341     2G2F_A   30.108             279        171       10
    315  Query_1615341     3EKK_A   30.293             307        177       12
    316  Query_1615341     8SCV_A   30.712             267        164        7
    317  Query_1615341     7C2W_A   31.298             262        159        7
    318  Query_1615341     8X2A_A   31.387             274        162       10
    319  Query_1615341     8H7F_A   29.720             286        177       10
    320  Query_1615341     3DK6_A   30.605             281        171       10
    321  Query_1615341     6F3D_A   31.298             262        159        7
    322  Query_1615341    9RPV_N1   31.618             272        155        9
    323  Query_1615341     3SXR_A   31.022             274        163       10
    324  Query_1615341     6F3I_A   30.712             267        164        7
    325  Query_1615341     6MOM_A   30.712             267        164        7
    326  Query_1615341     4RMZ_A   30.712             267        164        7
    327  Query_1615341     2OIB_A   30.712             267        164        7
    328  Query_1615341     5UIT_A   30.712             267        164        7
    329  Query_1615341     3PYY_A   30.108             279        171       10
    330  Query_1615341     4Y73_A   30.712             267        164        7
    331  Query_1615341     7CC2_A   30.108             279        171       10
    332  Query_1615341     7C2V_A   30.712             267        164        7
    333  Query_1615341     7QG3_A   30.712             267        164        7
    334  Query_1615341     9NA2_A   30.712             267        164        7
    335  Query_1615341     5W84_A   30.712             267        164        7
    336  Query_1615341   8BR5_AAA   30.712             267        164        7
    337  Query_1615341     6N8G_A   31.298             262        159        7
    338  Query_1615341     6I99_A   31.022             274        163       10
    339  Query_1615341     6F3E_A   31.298             262        159        7
    340  Query_1615341     6JK8_A   31.741             293        160       11
    341  Query_1615341     5E1S_A   30.421             309        174       12
    342  Query_1615341     1OPK_A   29.452             292        182       10
    343  Query_1615341     5HHW_A   30.421             309        174       12
    344  Query_1615341     6F3G_A   31.298             262        159        7
    345  Query_1615341     4Y93_A   28.205             273        172        9
    346  Query_1615341     3K54_A   28.571             273        171        9
    347  Query_1615341     5UIS_A   30.712             267        164        7
    348  Query_1615341     5UIQ_A   30.712             267        164        7
    349  Query_1615341     2GQG_A   30.108             279        171       10
    350  Query_1615341     6O94_A   30.712             267        164        7
    351  Query_1615341     8V1O_A   31.298             262        159        7
    352  Query_1615341     7SL1_A   30.392             306        178       11
    353  Query_1615341     7W7Y_A   30.108             279        171       10
    354  Query_1615341     1OPL_A   29.452             292        182       10
    355  Query_1615341     7W7X_A   30.108             279        171       10
    356  Query_1615341     8EYR_A   31.507             292        163       11
    357  Query_1615341     8DTL_A   30.392             306        178       11
    358  Query_1615341     3QRJ_A   29.749             279        172       10
    359  Query_1615341     2FO0_A   29.452             292        182       10
    360  Query_1615341     4YFF_A   31.227             269        166        8
    361  Query_1615341     4U97_A   30.337             267        165        7
    362  Query_1615341     4YFI_A   31.227             269        166        8
    363  Query_1615341     6AUA_A   28.205             273        172        9
    364  Query_1615341     7MGJ_A   31.227             269        166        8
    365  Query_1615341     6EGD_A   30.337             267        165        7
    366  Query_1615341     8S9F_A   28.938             273        170        9
    367  Query_1615341     8WTF_A   30.337             267        165        7
    368  Query_1615341     8EYX_A   30.392             306        178       11
    369  Query_1615341     4XI2_A   28.938             273        170        9
    370  Query_1615341     2Z60_A   29.749             279        172       10
    371  Query_1615341     3OY3_A   29.749             279        172       10
    372  Query_1615341     8DKS_A   31.298             262        159        7
    373  Query_1615341     4TWP_A   29.893             281        169       11
    374  Query_1615341     5MO4_A   29.110             292        183       10
    375  Query_1615341     6PYH_A   31.164             292        164       11
    376  Query_1615341     5BPY_A   28.309             272        171        9
    377  Query_1615341     3DK7_A   29.893             281        173       10
    378  Query_1615341     4NWM_A   28.309             272        171        9
    379  Query_1615341     2V7A_A   29.749             279        172       10
    380  Query_1615341     5AMN_A   31.469             286        166        8
    381  Query_1615341     7KHJ_A   31.419             296        163       10
    382  Query_1615341     7KHK_A   31.419             296        163       10
    383  Query_1615341     8FD9_A   27.574             272        173        9
    384  Query_1615341     1K2P_A   28.195             266        167        9
    385  Query_1615341     8GMB_A   28.571             273        171       10
    386  Query_1615341     6FEK_A   31.579             285        166        8
    387  Query_1615341     9KS5_A   28.623             276        179        8
    388  Query_1615341     7BW7_A   30.619             307        176       12
    389  Query_1615341     8U4B_A   30.619             307        176       12
    390  Query_1615341     7PG0_A   30.619             307        176       12
    391  Query_1615341     8VJB_A   30.619             307        176       12
    392  Query_1615341     2IVV_A   30.667             300        164        8
    393  Query_1615341     3ETA_A   30.769             286        158       11
    394  Query_1615341     8S93_A   28.309             272        171       10
    395  Query_1615341     6NE7_A   30.333             300        165        8
    396  Query_1615341     4HVS_A   31.081             296        164       10
    397  Query_1615341     6PXV_A   30.293             307        177       12
    398  Query_1615341     6NJA_A   30.333             300        165        8
    399  Query_1615341     4GU9_A   29.286             280        173        9
    400  Query_1615341     9SDI_A   29.286             280        173        9
    401  Query_1615341     2IVT_A   30.333             300        165        8
    402  Query_1615341     3PXK_A   29.286             280        173        9
    403  Query_1615341     6I8Z_A   29.286             280        173        9
    404  Query_1615341     2ETM_A   29.286             280        173        9
    405  Query_1615341     2IVS_A   30.333             300        165        8
    406  Query_1615341     2O8Y_A   29.963             267        166        7
    407  Query_1615341     6I83_A   30.333             300        165        8
    408  Query_1615341     1MP8_A   29.286             280        173        9
    409  Query_1615341     2R4B_A   30.943             265        166        7
    410  Query_1615341     2J0M_B   29.286             280        173        9
    411  Query_1615341     2JKM_A   29.286             280        173        9
    412  Query_1615341     3BBT_B   30.827             266        167        7
    413  Query_1615341     4EBV_A   29.286             280        173        9
    414  Query_1615341     6GQJ_A   31.186             295        163       10
    415  Query_1615341     8PQ9_A   31.186             295        163       10
    416  Query_1615341     1U54_A   30.483             269        163        9
    417  Query_1615341     3ZZW_A   28.720             289        164       12
    418  Query_1615341     4EWH_A   30.483             269        163        9
    419  Query_1615341     6GQK_A   31.186             295        163       10
    420  Query_1615341     6VQM_A   30.483             269        163        9
    421  Query_1615341     1U46_A   30.483             269        163        9
    422  Query_1615341     4XCU_A   30.201             298        166       12
    423  Query_1615341     7KP6_A   30.483             269        163        9
    424  Query_1615341     3BZ3_A   29.242             277        171        9
    425  Query_1615341     4GT4_A   28.720             289        164       12
    426  Query_1615341     5ZXB_A   30.483             269        163        9
    427  Query_1615341     2JKK_A   29.286             280        173        9
    428  Query_1615341     8PQG_A   31.525             295        162       10
    429  Query_1615341     8FE9_A   30.483             269        163        9
    430  Query_1615341     3EQP_A   30.483             269        163        9
    431  Query_1615341     2J0L_A   29.286             280        173        9
    432  Query_1615341     4ID7_A   30.483             269        163        9
    433  Query_1615341     4HZR_A   30.483             269        163        9
    434  Query_1615341     5FM2_A   30.132             302        163       10
    435  Query_1615341     4CKI_A   30.000             300        166        8
    436  Query_1615341     3QUP_A   28.912             294        176       10
    437  Query_1615341     6SDC_A   30.996             271        163       10
    438  Query_1615341     8ANS_A   30.996             271        163       10
    439  Query_1615341     9EF1_A   29.568             301        168       10
    440  Query_1615341     4HZS_A   30.483             269        163        9
    441  Query_1615341     7RUN_A   30.132             302        163       10
    442  Query_1615341     8AU5_A   30.996             271        163       10
    443  Query_1615341     6VHG_A   30.000             300        166        8
    444  Query_1615341     2J0J_A   27.479             353        225       11
    445  Query_1615341     6JPE_A   29.316             307        166       12
    446  Query_1615341     4UXQ_A   29.316             307        166       12
    447  Query_1615341     8KH9_A   29.316             307        166       12
    448  Query_1615341     7DTZ_A   29.316             307        166       12
    449  Query_1615341     5JKG_A   29.316             307        166       12
    450  Query_1615341     8W5C_A   29.316             307        166       12
    451  Query_1615341     7F3M_A   29.316             307        166       12
    452  Query_1615341     6YI8_A   29.316             307        166       12
    453  Query_1615341     4QQT_A   29.316             307        166       12
    454  Query_1615341     5NUD_A   29.316             307        166       12
    455  Query_1615341     5XFJ_A   29.316             307        166       12
    456  Query_1615341     3D7T_A   30.403             273        164       10
    457  Query_1615341     7VJL_A   29.316             307        166       12
    458  Query_1615341     4TYE_A   29.316             307        166       12
    459  Query_1615341     6I82_A   30.000             300        166        8
    460  Query_1615341     3T9T_A   27.143             280        178       10
    461  Query_1615341     1BYG_A   30.403             273        164       10
    462  Query_1615341     7WCW_A   29.316             307        166       12
    463  Query_1615341     5XFF_A   29.316             307        166       12
    464  Query_1615341     7WCX_A   29.316             307        166       12
    465  Query_1615341     6NFY_A   29.123             285        170        9
    466  Query_1615341     3D7U_A   30.403             273        164       10
    467  Query_1615341     4KNB_A   30.627             271        164       10
    468  Query_1615341     8AW1_A   30.627             271        164       10
    469  Query_1615341     4R1V_A   30.627             271        164       10
    470  Query_1615341   7PI4_DDD   28.986             276        171        9
    471  Query_1615341     4QQC_A   29.316             307        166       12
    472  Query_1615341     6SD9_A   30.627             271        164       10
    473  Query_1615341     8AN8_A   30.627             271        164       10
    474  Query_1615341     8VI1_A   30.627             271        164       10
    475  Query_1615341     6TY3_A   27.479             353        225       11
    476  Query_1615341     4QQJ_A   29.316             307        166       12
    477  Query_1615341     2WD1_A   30.627             271        164       10
    478  Query_1615341     2J0K_A   27.479             353        225       11
    479  Query_1615341     3GQI_A   28.231             294        174       10
    480  Query_1615341     7B3Q_A   30.627             271        164       10
    481  Query_1615341     6IUO_A   29.316             307        166       12
    482  Query_1615341     3F66_A   30.627             271        164       10
    483  Query_1615341     2G15_A   30.627             271        164       10
    484  Query_1615341     4EEV_A   30.627             271        164       10
    485  Query_1615341     3Q6U_A   30.627             271        164       10
    486  Query_1615341     8PAS_A   29.123             285        170        9
    487  Query_1615341     3LQ8_A   30.627             271        164       10
    488  Query_1615341     8FH4_A   28.826             281        176        9
    489  Query_1615341     6CQD_A   28.826             281        176        9
    490  Query_1615341     9C1R_A   30.627             271        164       10
    491  Query_1615341     2WGJ_A   30.627             271        164       10
    492  Query_1615341     8CDW_A   28.826             281        176        9
    493  Query_1615341     4QQ5_A   29.316             307        166       12
    494  Query_1615341     8PAU_A   28.669             293        177        9
    495  Query_1615341     6NG0_A   28.826             281        176        9
    496  Query_1615341     3VW8_A   30.627             271        164       10
    497  Query_1615341     8PAR_A   28.826             281        176        9
    498  Query_1615341     5A46_A   28.231             294        174       10
    499  Query_1615341     7M0L_A   28.374             289        183        9
    500  Query_1615341     9SZJ_A   30.627             271        164       10
    501  Query_1615341     4GG5_A   30.627             271        164       10
    502  Query_1615341     5UAB_A   30.627             271        164       10
    503  Query_1615341     5ZV2_A   28.231             294        174       10
    504  Query_1615341     3I5N_A   30.627             271        164       10
    505  Query_1615341     7M0M_A   28.374             289        183        9
    506  Query_1615341     8YKI_A   28.231             294        174       10
    507  Query_1615341     5A4C_A   28.231             294        174       10
    508  Query_1615341     3C4F_A   28.231             294        174       10
    509  Query_1615341     7M0K_A   28.374             289        183        9
    510  Query_1615341     4ZSA_A   28.231             294        174       10
    511  Query_1615341     3RHX_A   28.231             294        174       10
    512  Query_1615341     2RFN_A   30.627             271        164       10
    513  Query_1615341     9H8D_A   28.826             281        176        9
    514  Query_1615341     8XN7_A   28.826             281        176        9
    515  Query_1615341     4F63_A   28.231             294        174       10
    516  Query_1615341     4WUN_A   28.231             294        174       10
    517  Query_1615341     5FLF_A   28.231             294        174       10
    518  Query_1615341     5HLW_A   30.627             271        164       10
    519  Query_1615341     7L25_A   28.328             293        178        9
    520  Query_1615341     1AGW_A   28.231             294        174       10
    521  Query_1615341     6LVM_A   28.571             294        173       11
    522  Query_1615341     4F0F_A   28.269             283        169       11
    523  Query_1615341     3KXX_A   28.231             294        174       10
    524  Query_1615341     4PWN_A   33.460             263        147       10
    525  Query_1615341     5TQ5_A   29.787             282        161       11
    526  Query_1615341     3JS2_A   28.231             294        174       10
    527  Query_1615341     4HCT_A   26.855             283        181       10
    528  Query_1615341     7WCL_A   28.231             294        174       10
    529  Query_1615341     6CQF_A   28.772             285        171        9
    530  Query_1615341     5AM7_A   28.231             294        174       10
    531  Query_1615341     8Y22_A   28.231             294        174       10
    532  Query_1615341     9BI8_A   28.826             281        176        9
    533  Query_1615341     3GQL_A   28.231             294        174       10
    534  Query_1615341     6MZW_A   28.231             294        174       10
    535  Query_1615341     7L24_A   28.328             293        178        9
    536  Query_1615341     6HH1_A   31.399             293        163       10
    537  Query_1615341     3DKG_A   30.258             271        165       10
    538  Query_1615341     6NVL_A   28.231             294        174       10
    539  Query_1615341     8G6Z_A   29.537             281        163       10
    540  Query_1615341     3Q6W_A   30.627             271        164       10
    541  Query_1615341     6CQE_A   28.772             285        171        9
    542  Query_1615341     4RWI_A   28.231             294        174       10
    543  Query_1615341     7L26_A   28.328             293        178        9
    544  Query_1615341     6NFZ_A   28.826             281        176        9
    545  Query_1615341     7LL5_A   31.095             283        156       13
    546  Query_1615341     5J5T_A   29.151             271        173        7
    547  Query_1615341     7KAC_A   29.123             285        170        9
    548  Query_1615341     7LL4_A   31.095             283        156       13
    549  Query_1615341     3KRR_A   30.303             297        168       13
    550  Query_1615341     7RN6_A   31.095             283        156       13
    551  Query_1615341     5TQ4_A   29.787             282        161       11
    552  Query_1615341     8AU3_A   30.627             271        164       10
    553  Query_1615341     3QTI_A   30.258             271        165       10
    554  Query_1615341     3CE3_A   30.258             271        165       10
    555  Query_1615341     3C1X_A   30.258             271        165       10
    556  Query_1615341     4IWD_A   30.627             271        164       10
    557  Query_1615341     5GRN_A   28.527             319        175       10
    558  Query_1615341     3A4P_A   30.258             271        165       10
    559  Query_1615341     2P2H_A   28.912             294        170        9
    560  Query_1615341     8UDT_A   28.723             282        176        9
    561  Query_1615341     5WEV_A   29.874             318        172       14
    562  Query_1615341     3DKC_A   30.258             271        165       10
    563  Query_1615341     8PQJ_A   28.527             319        175       10
    564  Query_1615341     1R0P_A   30.258             271        165       10
    565  Query_1615341     4E4M_A   29.801             302        173       13
    566  Query_1615341     5TQ6_A   29.787             282        161       11
    567  Query_1615341     3MIY_A   27.338             278        176       10
    568  Query_1615341     5VND_A   28.231             294        174       10
    569  Query_1615341     5TQ3_A   29.787             282        161       11
    570  Query_1615341     5USY_A   31.095             283        156       13
    571  Query_1615341     9GB9_A   30.986             213        134        5
    572  Query_1615341     6VGL_A   31.095             283        156       13
    573  Query_1615341     6TU9_A   29.861             288        162       12
    574  Query_1615341     8UDV_A   28.723             282        176        9
    575  Query_1615341     4HGE_A   29.801             302        173       13
    576  Query_1615341     3UGC_A   30.303             297        168       13
    577  Query_1615341     8G8O_A   29.787             282        161       11
    578  Query_1615341     1K9A_A   30.403             273        164       10
    579  Query_1615341     3TJC_A   31.095             283        156       13
    580  Query_1615341     6AAJ_A   29.801             302        173       13
    581  Query_1615341     3C7Q_A   28.767             292        171        9
    582  Query_1615341     3TT0_A   28.231             294        174       10
    583  Query_1615341     7MO7_B   30.627             271        164       10
    584  Query_1615341     2OH4_A   28.669             293        171        9
    585  Query_1615341     9V70_A   32.178             202        123        6
    586  Query_1615341     8BX6_A   30.201             298        169       13
    587  Query_1615341     3E62_A   31.095             283        156       13
    588  Query_1615341     2XA4_A   31.095             283        156       13
    589  Query_1615341     4D0W_A   31.095             283        156       13
    590  Query_1615341     7Q7I_A   29.787             282        161       11
    591  Query_1615341     2B7A_A   31.095             283        156       13
    592  Query_1615341     4AQC_A   31.095             283        156       13
    593  Query_1615341     8BM2_A   30.201             298        169       13
    594  Query_1615341     3EWH_A   28.912             294        170        9
    595  Query_1615341     3U6J_A   28.912             294        170        9
    596  Query_1615341     3ZMM_A   31.095             283        156       13
    597  Query_1615341     9V71_A   32.178             202        123        6
    598  Query_1615341     4BBE_A   29.866             298        170       13
    599  Query_1615341     5HEZ_A   29.470             302        174       12
    600  Query_1615341     3CJF_A   28.966             290        171        9
    601  Query_1615341     1T46_A   30.976             297        163       10
    602  Query_1615341     4F1O_A   27.972             286        172       11
    603  Query_1615341     3RVG_A   29.766             299        171       13
    604  Query_1615341     6WTN_A   31.095             283        156       13
    605  Query_1615341     2W1I_A   31.095             283        156       13
    606  Query_1615341     7UYW_A   31.095             283        156       13
    607  Query_1615341     2OGV_A   30.241             291        171        9
    608  Query_1615341     2I0V_A   29.966             297        170        9
    609  Query_1615341     1YWN_A   28.571             294        171        9
    610  Query_1615341     4U0I_A   30.976             297        163       10
    611  Query_1615341     6MOB_A   30.976             297        163       10
    612  Query_1615341     6WXJ_A   29.966             297        170        9
    613  Query_1615341     6TPD_A   29.787             282        161       11
    614  Query_1615341     3Q32_A   31.095             283        156       13
    615  Query_1615341     4YTC_A   31.095             283        156       13
    616  Query_1615341     4F1M_A   27.915             283        170       11
    617  Query_1615341     9D3F_A   33.840             263        145       11
    618  Query_1615341     2X7F_A   29.964             277        163       13
    619  Query_1615341     8W1L_A   29.966             297        170        9
    620  Query_1615341     7TNH_A   29.966             297        170        9
    621  Query_1615341     3JY9_A   30.000             300        171       13
    622  Query_1615341     5HOA_A   29.889             271        166       10
    623  Query_1615341     3CJG_A   28.966             290        171        9
    624  Query_1615341     7SIU_A   28.772             285        171        9
    625  Query_1615341     5UGL_A   27.211             294        177       10
    626  Query_1615341     3IO7_A   30.000             300        171       13
    627  Query_1615341     7XZQ_A   29.964             277        163       13
    628  Query_1615341     7ZVS_A   29.787             282        163       14
    629  Query_1615341     7FCZ_A   28.253             269        164        9
    630  Query_1615341     5AX9_A   29.964             277        163       13
    631  Query_1615341     6RA5_A   29.964             277        163       13
    632  Query_1615341     8ZML_A   29.964             277        163       13
    633  Query_1615341     3V5J_A   26.978             278        177       10
    634  Query_1615341     6C4D_A   28.253             269        164        9
    635  Query_1615341     9HY8_A   28.253             269        164        9
    636  Query_1615341     1VR2_A   28.571             294        171        9
    637  Query_1615341     7UOS_A   33.840             263        145       11
    638  Query_1615341     1SM2_A   26.978             278        177       10
    639  Query_1615341     4ZIM_A   31.095             283        156       13
    640  Query_1615341     6NYH_A   28.253             269        164        9
    641  Query_1615341     9WYR_A   31.683             202        124        6
    642  Query_1615341     9KFU_A   28.231             294        174       11
    643  Query_1615341     6OL2_A   33.840             263        145       11
    644  Query_1615341     8X88_A   29.964             277        163       13
    645  Query_1615341     4ITH_A   28.253             269        164        9
    646  Query_1615341     5DRB_A   33.840             263        145       11
    647  Query_1615341     2P2I_A   28.571             294        171        9
    648  Query_1615341     3VNT_A   28.571             294        171        9
    649  Query_1615341     3PLS_A   29.044             272        169        9
    650  Query_1615341     9CD7_A   28.716             296        170       13
    651  Query_1615341     8W3D_A   27.211             294        177       10
    652  Query_1615341     9MZZ_A   28.358             268        163        9
    653  Query_1615341     8W3B_A   27.211             294        177       10
    654  Query_1615341     9MZY_A   28.358             268        163        9
    655  Query_1615341     9MZX_A   28.358             268        163        9
    656  Query_1615341     8W2X_A   27.211             294        177       10
    657  Query_1615341     8JQI_B   28.231             294        174       10
    658  Query_1615341     4K33_A   28.716             296        170       13
    659  Query_1615341     1PKG_A   30.976             297        163       10
    660  Query_1615341     8W38_A   27.211             294        177       10
    661  Query_1615341     4Q2A_A   33.840             263        145       11
    662  Query_1615341     7TEU_A   29.766             299        171       13
    663  Query_1615341     2XIR_A   28.571             294        171        9
    664  Query_1615341     3LCD_A   30.241             291        171        9
    665  Query_1615341     9GZH_A   28.213             319        176       10
    666  Query_1615341     8PQH_A   28.213             319        176       10
    667  Query_1615341     5W7T_A   33.840             263        145       11
    668  Query_1615341     5HX6_A   28.253             269        164        9
    669  Query_1615341     4APC_A   28.148             270        171        7
    670  Query_1615341     7ZW8_A   31.313             297        162       11
    671  Query_1615341     5HOR_A   29.889             271        166       10
    672  Query_1615341     6ITV_A   30.976             297        163       10
    673  Query_1615341     6A32_A   28.213             319        176       10
    674  Query_1615341     6C3E_A   28.253             269        164        9
    675  Query_1615341     6GQO_A   28.571             294        171        9
    676  Query_1615341     4E6D_A   29.431             299        172       13
    677  Query_1615341     3WZD_A   28.571             294        171        9
    678  Query_1615341     6PNX_A   28.231             294        174       11
    679  Query_1615341     6NW2_A   28.253             269        164        9
    680  Query_1615341     5UGX_A   27.211             294        177       10
    681  Query_1615341     3G0E_A   30.976             297        163       10
    682  Query_1615341     9GTG_A   28.253             269        164        9
    683  Query_1615341     1T45_A   30.976             297        163       10
    684  Query_1615341     6JOI_A   28.213             319        176       10
    685  Query_1615341     3G0F_A   30.976             297        163       10
    686  Query_1615341     4J98_A   27.211             294        177       10
    687  Query_1615341     1FAQ_A   73.077              52         14        0
    688  Query_1615341     2PVF_A   26.871             294        178       10
    689  Query_1615341     6FD3_A   28.794             257        165        9
    690  Query_1615341     6AGX_A   26.871             294        178       10
    691  Query_1615341     8XRR_A   28.213             319        176       10
    692  Query_1615341     3CLY_A   26.871             294        178       10
    693  Query_1615341     5TF9_A   33.080             263        147       10
    694  Query_1615341     4USF_A   27.239             268        179        9
    695  Query_1615341     4AGC_A   28.571             294        171        9
    696  Query_1615341     3BEA_A   29.392             296        171       10
    697  Query_1615341     8BEM_A   27.239             268        179        9
    698  Query_1615341     4GL9_A   28.053             303        183       10
    699  Query_1615341     2I1M_A   30.068             296        169       12
    700  Query_1615341     2PZ5_A   26.871             294        178       10
    701  Query_1615341     8E1X_A   26.871             294        178       10
    702  Query_1615341     6BDN_A   30.078             256        163        7
    703  Query_1615341     2J51_A   27.239             268        179        9
    704  Query_1615341     2JFM_A   27.239             268        179        9
    705  Query_1615341     2JFL_A   27.239             268        179        9
    706  Query_1615341     2PZP_A   26.871             294        178       10
    707  Query_1615341     4J97_A   26.871             294        178       10
    708  Query_1615341     4KIO_A   26.978             278        177       10
    709  Query_1615341     8H75_A   26.531             294        179        9
    710  Query_1615341     4NEU_A   28.253             269        164        9
    711  Query_1615341     2PWL_A   26.871             294        178       10
    712  Query_1615341     6LVK_A   26.531             294        179        9
    713  Query_1615341     9U7E_A   26.531             294        179        9
    714  Query_1615341     8SWE_A   26.871             294        178       10
    715  Query_1615341     1GJO_A   26.871             294        178       10
    716  Query_1615341     3RI1_A   26.871             294        178       10
    717  Query_1615341     5EG3_A   26.871             294        178       10
    718  Query_1615341     3B2T_A   26.871             294        178       10
    719  Query_1615341     8STG_A   26.871             294        178       10
    720  Query_1615341     2PVY_A   26.871             294        178       10
    721  Query_1615341     6LVL_A   26.531             294        179        9
    722  Query_1615341     9GFZ_A   32.836             201        121        6
    723  Query_1615341     3QGW_A   26.619             278        178       10
    724  Query_1615341     9LBG_A   26.045             311        206       11
    725  Query_1615341   7OZY_AAA   26.531             294        179        9
    726  Query_1615341     8E4T_A   28.417             278        173       10
    727  Query_1615341     4J96_A   26.871             294        178       10
    728  Query_1615341     4J99_A   26.871             294        178       10
    729  Query_1615341     2PSQ_A   26.871             294        178       10
    730  Query_1615341     9U3N_A   26.871             294        178       10
    731  Query_1615341     5UHN_A   26.871             294        178       10
    732  Query_1615341     7KIA_A   26.531             294        179        9
    733  Query_1615341     2Q0B_A   26.871             294        178       10
    734  Query_1615341     8JOT_A   28.025             314        172       10
    735  Query_1615341     2PZR_A   26.531             294        179       10
    736  Query_1615341     9LBF_A   26.045             311        206       11
    737  Query_1615341     2PY3_A   26.871             294        178       10
    738  Query_1615341     3LCO_A   28.947             304        171        9
    739  Query_1615341     5UI0_A   26.871             294        178       10
    740  Query_1615341     6T2W_A   28.947             304        171        9
    741  Query_1615341     9D51_A   27.237             257        169        9
    742  Query_1615341     4HW7_A   28.383             303        174        9
    743  Query_1615341     7AAY_A   25.000             296        190       10
    744  Query_1615341     8CGC_A   28.947             304        171        9
    745  Query_1615341     4O27_B   27.652             264        172        9
    746  Query_1615341     9BHI_A   24.671             304        193       10
    747  Query_1615341     5U6B_A   25.784             287        183        8
    748  Query_1615341     9D7Q_A   32.558             258        141       12
    749  Query_1615341     9M41_A   25.723             311        207       11
    750  Query_1615341     7DXL_A   24.667             300        190       10
    751  Query_1615341     6V6Q_A   26.531             294        179       10
    752  Query_1615341     2GCD_A   29.070             258        167        7
    753  Query_1615341     3S95_A   28.520             277        176        7
    754  Query_1615341     1U5Q_A   29.070             258        167        7
    755  Query_1615341     4QML_A   27.652             264        172        9
    756  Query_1615341     7AAZ_A   24.667             300        190       10
    757  Query_1615341     7OAM_A   24.667             300        190       10
    758  Query_1615341     3A7F_A   27.652             264        172        9
    759  Query_1615341     7CQE_A   24.667             300        190       10
    760  Query_1615341     4W8E_A   27.652             264        172        9
    761  Query_1615341     5TC0_A   24.667             300        190       10
    762  Query_1615341     4U8Z_A   27.652             264        172        9
    763  Query_1615341     7B30_A   27.652             264        172        9
    764  Query_1615341     2P0C_A   24.667             300        190       10
    765  Query_1615341     5O1V_A   32.558             258        141       12
    766  Query_1615341     7Z5W_A   30.124             322        177       15
    767  Query_1615341     5U6C_A   24.667             300        190       10
    768  Query_1615341     3ZHP_C   27.652             264        172        9
    769  Query_1615341     5O2B_A   32.558             258        141       12
    770  Query_1615341     8SE1_A   25.541             462        290       18
    771  Query_1615341     8SE1_A   31.034              58         36        1
    772  Query_1615341     8SE2_A   25.541             462        290       18
    773  Query_1615341     8SE2_A   31.034              58         36        1
    774  Query_1615341     2XIK_A   29.658             263        168        9
    775  Query_1615341     3CKW_A   27.652             264        172        9
    776  Query_1615341     1LUF_A   29.630             297        159       11
    777  Query_1615341     3CKX_A   27.652             264        172        9
    778  Query_1615341     7Z4V_A   29.658             263        168        9
    779  Query_1615341     4V0G_B   28.470             281        167        9
    780  Query_1615341     3PFQ_A   25.054             463        291       18
    781  Query_1615341     3PFQ_A   31.034              58         36        1
    782  Query_1615341     4RIO_A   28.470             281        167        9
    783  Query_1615341     5CNN_A   29.057             265        171        7
    784  Query_1615341     8HV2_A   29.057             265        171        7
    785  Query_1615341     5HVJ_A   28.159             277        177        7
    786  Query_1615341     3LXN_A   29.348             276        156       10
    787  Query_1615341     4ZJV_A   29.057             265        171        7
    788  Query_1615341     5HVK_A   28.159             277        177        7
    789  Query_1615341     4GVJ_A   29.348             276        156       10
    790  Query_1615341     8EDH_A   33.074             257        141       12
    791  Query_1615341     7SYD_A   28.947             266        170        7
    792  Query_1615341     5O26_A   32.296             257        143       11
    793  Query_1615341     4V0G_A   28.470             281        167        9
    794  Query_1615341     5CAV_A   29.057             265        171        7
    795  Query_1615341     1M14_A   29.057             265        171        7
    796  Query_1615341     7UKV_A   29.057             265        171        7
    797  Query_1615341     1YVJ_A   28.114             281        168        9
    798  Query_1615341     3D5V_A   30.667             225        139        7
    799  Query_1615341     4HVD_A   28.315             279        170        8
    800  Query_1615341     5O2C_A   32.558             258        141       12
    801  Query_1615341     4TKS_A   29.057             265        171        7
    802  Query_1615341     4Z16_A   28.114             281        168        9
    803  Query_1615341     4I23_A   29.057             265        171        7
    804  Query_1615341     3D5U_A   30.667             225        139        7
    805  Query_1615341     3VJO_A   29.057             265        171        7
    806  Query_1615341     2RFD_A   29.057             265        171        7
    807  Query_1615341     2ITW_A   29.057             265        171        7
    808  Query_1615341     9FZR_A   29.057             265        171        7
    809  Query_1615341     3D5W_A   30.667             225        139        7
    810  Query_1615341     4G5J_A   29.057             265        171        7
    811  Query_1615341     4JQ7_A   29.057             265        171        7
    812  Query_1615341     4LI5_A   29.057             265        171        7
    813  Query_1615341     4RIW_B   29.057             265        171        7
    814  Query_1615341     9H42_A   29.057             265        171        7
    815  Query_1615341     2GS2_A   29.057             265        171        7
    816  Query_1615341     5FED_A   29.057             265        171        7
    817  Query_1615341     9UBI_A   26.568             271        179        8
    818  Query_1615341     4WKQ_A   29.057             265        171        7
    819  Query_1615341     3D5X_A   30.667             225        139        7
    820  Query_1615341     6JZ0_A   29.057             265        171        7
    821  Query_1615341     6CN9_A   33.716             261        148       11
    822  Query_1615341     7UYV_A   28.470             281        167        9
    823  Query_1615341   7Q6H_AAA   28.470             281        167        9
    824  Query_1615341     2J5E_A   29.057             265        171        7
    825  Query_1615341     8A27_A   29.057             265        171        7
    826  Query_1615341     7AEM_A   29.057             265        171        7
    827  Query_1615341     8PO3_A   29.057             265        171        7
    828  Query_1615341     2RGP_A   29.057             265        171        7
    829  Query_1615341     8PO4_A   29.057             265        171        7
    830  Query_1615341     1XKK_A   29.057             265        171        7
    831  Query_1615341     3ZC6_A   28.470             281        167        9
    832  Query_1615341     4PY1_A   28.986             276        157       10
    833  Query_1615341     5W86_A   28.470             281        167        9
    834  Query_1615341     8H7X_A   29.057             265        171        7
    835  Query_1615341     9N6G_A   29.057             265        171        7
    836  Query_1615341     5NG0_A   27.241             290        192       10
    837  Query_1615341     5XGN_A   29.057             265        171        7
    838  Query_1615341     9BY4_A   29.057             265        171        7
    839  Query_1615341     7XDY_A   26.568             271        179        8
    840  Query_1615341     7XDV_A   26.568             271        179        8
    841  Query_1615341     7XDX_A   26.568             271        179        8
    842  Query_1615341     2GS7_A   29.057             265        171        7
    843  Query_1615341     5TOZ_A   28.315             279        170        8
    844  Query_1615341     6AAM_A   28.986             276        157       10
    845  Query_1615341     5W5J_A   27.241             290        192       10
    846  Query_1615341     3LXK_A   28.315             279        170        8
    847  Query_1615341     9QBG_A   29.434             265        170        6
    848  Query_1615341     9QBF_A   29.434             265        170        6
    849  Query_1615341     6HMX_A   27.437             277        185        9
    850  Query_1615341     4HJO_A   29.057             265        171        7
    851  Query_1615341     4TPT_A   30.986             284        165       11
    852  Query_1615341     8X2O_A   27.241             290        192       10
    853  Query_1615341     9F3V_A   27.241             290        192       10
    854  Query_1615341     6ES0_A   27.891             294        185       11
    855  Query_1615341     6UL8_A   27.241             290        192       10
    856  Query_1615341     5NXD_A   30.986             284        165       11
    857  Query_1615341   9QXN_AAA   29.057             265        171        7
    858  Query_1615341     5AR2_A   27.437             277        185        9
    859  Query_1615341     8AZA_A   27.891             294        185       11
    860  Query_1615341     5NG3_B   27.076             277        186        9
    861  Query_1615341     6HZV_A   28.470             281        167        9
    862  Query_1615341     3LZB_A   29.057             265        171        7
    863  Query_1615341     4C8B_A   28.114             281        178       10
    864  Query_1615341     5ZWJ_A   29.057             265        171        7
    865  Query_1615341     6LUB_A   28.679             265        172        7
    866  Query_1615341     8HV4_A   28.679             265        172        7
    867  Query_1615341     2V5Q_A   32.673             202        119        7
    868  Query_1615341     4E4L_A   25.623             281        178        9
    869  Query_1615341     6N7A_A   25.623             281        178        9
    870  Query_1615341     3IKA_A   28.679             265        172        7
    871  Query_1615341     4GIH_A   28.986             276        157       10
    872  Query_1615341     6GTT_A   26.354             277        176       10
    873  Query_1615341     8WD4_A   28.679             265        172        7
    874  Query_1615341     9D3V_A   28.679             265        172        7
    875  Query_1615341     5TD2_A   24.579             297        188       10
    876  Query_1615341     6ELR_A   25.623             281        178        9
    877  Query_1615341     7TVD_A   28.352             261        175        6
    878  Query_1615341     6GGH_A   25.623             281        178        9
    879  Query_1615341     3PJC_A   28.114             281        168        9
    880  Query_1615341     2YAC_A   32.673             202        119        7
    881  Query_1615341     9S3X_A   28.679             265        172        7
    882  Query_1615341     3PP0_A   28.679             265        172        6
    883  Query_1615341     7PCD_A   28.679             265        172        6
    884  Query_1615341     2QKW_B   27.372             274        177        7
    885  Query_1615341     9QEK_A   30.097             309        170       14
    886  Query_1615341     3ZBF_A   30.097             309        170       14
    887  Query_1615341     7JXH_A   28.679             265        172        6
    888  Query_1615341     3KB7_A   32.673             202        119        7
    889  Query_1615341     2JIU_A   28.679             265        172        7
    890  Query_1615341     4I24_A   28.679             265        172        7
    891  Query_1615341     5KHW_A   25.623             281        178        9
    892  Query_1615341     7SZ0_A   28.571             266        171        7
    893  Query_1615341     5GMP_A   28.679             265        172        7
    894  Query_1615341     5Y9T_A   28.679             265        172        7
    895  Query_1615341     5FEE_A   28.679             265        172        7
    896  Query_1615341     3EYG_A   25.623             281        178        9
    897  Query_1615341     4LQM_A   28.679             265        172        7
    898  Query_1615341     3HGK_A   27.372             274        177        7
    899  Query_1615341     6TPE_A   25.623             281        178        9
    900  Query_1615341     5GNK_A   28.679             265        172        7
    901  Query_1615341     5AJQ_A   28.319             226        137        8
    902  Query_1615341     4G5P_A   28.679             265        172        7
    903  Query_1615341     2J7T_A   28.319             226        137        8
    904  Query_1615341     2JIT_A   28.679             265        172        7
    905  Query_1615341     6EIM_A   28.319             226        137        8
    906  Query_1615341     4BC6_A   28.319             226        137        8
    907  Query_1615341     5XDL_A   28.679             265        172        7
    908  Query_1615341     4YNE_A   27.619             315        173       12
    909  Query_1615341     2OU7_A   32.673             202        119        7
    910  Query_1615341     4NZW_B   29.278             263        169        9
    911  Query_1615341     5J9Z_A   28.679             265        172        7
    912  Query_1615341     5J9Y_A   28.679             265        172        7
    913  Query_1615341     4ZSE_A   28.679             265        172        7
    914  Query_1615341     5NG3_A   27.240             279        183        9
    915  Query_1615341     2EB2_A   28.679             265        172        7
    916  Query_1615341     4QPS_A   28.114             281        168        9
    917  Query_1615341     2ITN_A   28.679             265        172        7
    918  Query_1615341     6JWL_A   28.679             265        172        7
    919  Query_1615341     3THB_A   32.673             202        119        7
    920  Query_1615341     3GGF_A   27.757             263        173        9
    921  Query_1615341     6HXF_A   27.602             221        145        7
    922  Query_1615341     7C3N_A   27.957             279        171        8
    923  Query_1615341     4J52_A   32.353             204        121        7
    924  Query_1615341     2EB3_A   28.679             265        172        7
    925  Query_1615341     6P1D_A   28.679             265        172        7
    926  Query_1615341     2ITT_A   28.679             265        172        7
    927  Query_1615341     4H1J_A   26.022             269        174        8
    928  Query_1615341     9JQ1_A   28.679             265        172        7
    929  Query_1615341     6V5N_A   28.679             265        172        7
    930  Query_1615341     2RKU_A   32.353             204        121        7
    931  Query_1615341     6C9D_A   29.524             210        138        5
    932  Query_1615341     6TFU_A   28.679             265        172        7
    933  Query_1615341     8PEH_A   31.683             202        124        6
    934  Query_1615341     3ET7_A   26.022             269        174        8
    935  Query_1615341     4R3P_A   28.679             265        172        7
    936  Query_1615341     5TO8_A   26.022             269        174        8
    937  Query_1615341     5LWM_A   28.114             281        168        9
    938  Query_1615341     8S9P_C   29.153             295        163       11
    939  Query_1615341     7K1H_A   28.679             265        172        7
    940  Query_1615341     9P9U_A   29.057             265        171        7
    941  Query_1615341     9P9U_A   29.057             265        171        7
    942  Query_1615341     9XU9_A   28.679             265        172        7
    943  Query_1615341     3CC6_A   26.022             269        174        8
    944  Query_1615341     4I20_A   28.679             265        172        7
    945  Query_1615341     4RJ4_A   28.302             265        173        7
    946  Query_1615341     5JFS_A   27.445             317        175       12
    947  Query_1615341     7UYR_A   28.986             276        157       10
    948  Query_1615341     7VKO_A   27.619             315        173       12
    949  Query_1615341     6D22_A   27.445             317        175       12
    950  Query_1615341     8HY7_A   28.302             265        173        7
    951  Query_1615341     4RIY_A   27.407             270        177        7
    952  Query_1615341     7VKM_A   27.619             315        173       12
    953  Query_1615341     5TA6_A   32.673             202        119        7
    954  Query_1615341     3GOP_A   28.679             265        172        7
    955  Query_1615341     9DF3_A   29.104             268        170        8
    956  Query_1615341     9DF2_A   29.104             268        170        8
    957  Query_1615341     4KS7_A   31.658             199        123        7
    958  Query_1615341     6C7Y_A   25.806             279        176        9
    959  Query_1615341     6S9B_A   28.302             265        173        7
    960  Query_1615341     6S9C_A   28.302             265        173        7
    961  Query_1615341     4J7B_A   31.401             207        125        7
    962  Query_1615341     5H3Q_A   27.832             309        168       12
    963  Query_1615341     7MN5_B   28.679             265        172        6
    964  Query_1615341     4AOJ_A   27.832             309        168       12
    965  Query_1615341     7MN6_B   28.679             265        172        6
    966  Query_1615341     8A2B_A   28.679             265        172        7
    967  Query_1615341     5KMI_A   27.832             309        168       12
    968  Query_1615341     3DAK_A   26.950             282        181        9
    969  Query_1615341     5KML_A   27.832             309        168       12
    970  Query_1615341     2HAK_A   29.524             210        138        5
    971  Query_1615341     3NZ0_A   28.623             276        158       10
    972  Query_1615341     2VWI_A   26.950             282        181        9
    973  Query_1615341     2C30_A   31.658             199        123        7
    974  Query_1615341     4ZP5_A   28.571             273        172        9
    975  Query_1615341     8DSW_A   29.104             268        171        8
    976  Query_1615341     4OBO_A   28.571             273        172        9
    977  Query_1615341     8J5W_A   27.832             309        168       12
    978  Query_1615341     5DI1_A   28.571             273        172        9
    979  Query_1615341     5J95_A   28.571             273        172        9
    980  Query_1615341     4F0I_A   27.832             309        168       12
    981  Query_1615341     2QNJ_A   27.600             250        157        8
    982  Query_1615341     4U3Z_A   28.571             273        172        9
    983  Query_1615341     4LRM_A   30.370             270        164        9
    984  Query_1615341     6NSP_A   27.832             309        168       12
    985  Query_1615341     6D1Y_A   27.832             309        168       12
    986  Query_1615341     7XAF_A   27.832             309        168       12
    987  Query_1615341     6IQN_A   27.832             309        168       12
    988  Query_1615341     2A19_B   28.421             285        167        8
    989  Query_1615341     4RVT_A   28.571             273        172        9
    990  Query_1615341     4E1Z_A   27.857             280        163       10
    991  Query_1615341     4E20_A   27.857             280        163       10
    992  Query_1615341     4GT5_A   27.832             309        168       12
    993  Query_1615341     7MN5_A   27.407             270        177        7
    994  Query_1615341     8J5X_A   27.653             311        170       12
    995  Query_1615341     4LL0_A   28.302             265        173        7
    996  Query_1615341     8PO0_A   30.370             270        164        9
    997  Query_1615341     8UOI_A   27.200             250        158        8
    998  Query_1615341     7AAX_A   24.667             300        190       10
    999  Query_1615341     4I5M_A   26.459             257        172        6
    1000 Query_1615341     9U8C_A   30.370             270        164        9
    1001 Query_1615341     7M5Z_A   24.667             300        190       10
    1002 Query_1615341     6NSS_A   27.653             311        170       12
    1003 Query_1615341     4I21_A   28.302             265        173        7
    1004 Query_1615341     8V5I_A   28.571             273        172        9
    1005 Query_1615341     3LMG_A   27.037             270        178        7
    1006 Query_1615341     3UG1_A   28.302             265        173        7
    1007 Query_1615341     7P1L_A   28.455             246        160        9
    1008 Query_1615341     4RIX_A   27.037             270        178        7
    1009 Query_1615341     4RIW_A   27.037             270        178        7
    1010 Query_1615341     7FEH_A   25.697             323        180       12
    1011 Query_1615341     8PO1_A   30.370             270        164        9
    1012 Query_1615341     5F1Z_A   28.261             276        159       10
    1013 Query_1615341     6YAT_A   26.724             232        162        5
    1014 Query_1615341     3NYX_A   28.261             276        159       10
    1015 Query_1615341     9KLW_A   28.302             265        173        7
    1016 Query_1615341     6OP9_A   27.037             270        178        7
    1017 Query_1615341     5SAU_A   25.697             323        180       12
    1018 Query_1615341     7MX3_A   28.519             270        156       10
    1019 Query_1615341   8C12_AAA   30.392             204        123        7
    1020 Query_1615341     8UOH_A   27.200             250        158        8
    1021 Query_1615341     3FE3_A   27.200             250        158        8
    1022 Query_1615341     3COM_A   26.724             232        162        5
    1023 Query_1615341     6CTH_A   28.788             198        132        5
    1024 Query_1615341     7OXB_A   28.302             265        173        7
    1025 Query_1615341     8UOJ_A   27.200             250        158        8
    1026 Query_1615341   9R1W_AAA   32.178             202        120        7
    1027 Query_1615341     3W2O_A   28.302             265        173        7
    1028 Query_1615341     6M0U_A   31.000             200        128        4
    1029 Query_1615341     3KEX_A   27.037             270        178        7
    1030 Query_1615341     8PAV_A   26.724             232        162        5
    1031 Query_1615341     7LGS_A   30.370             270        164        9
    1032 Query_1615341     4I1Z_A   28.302             265        173        7
    1033 Query_1615341     5Y25_A   28.302             265        173        7
    1034 Query_1615341     5WNI_A   28.205             273        182        8
    1035 Query_1615341     8KFQ_A   28.302             265        173        7
    1036 Query_1615341     9FQP_A   30.370             270        164        9
    1037 Query_1615341     7MON_B   28.519             270        156       10
    1038 Query_1615341     2F57_A   30.392             204        123        7
    1039 Query_1615341     8U8X_A   29.520             271        166        8
    1040 Query_1615341     4FZA_B   27.376             263        174        9
    1041 Query_1615341     8VB5_A   29.520             271        166        8
    1042 Query_1615341     9IIC_A   27.004             237        155        6
    1043 Query_1615341     5WR7_A   27.331             311        171       12
    1044 Query_1615341     8D73_A   28.302             265        173        7
    1045 Query_1615341     6NPT_A   27.331             311        171       12
    1046 Query_1615341     4PMM_A   28.155             309        167       13
    1047 Query_1615341     9LFU_A   28.309             272        158       10
    1048 Query_1615341     6D3K_A   27.797             295        163        8
    1049 Query_1615341     2WZJ_A   28.634             227        150        6
    1050 Query_1615341     6PL1_A   27.508             309        169       12
    1051 Query_1615341     4B6L_A   27.907             258        167        8
    1052 Query_1615341     2MSE_D   60.563              71         28        0
    1053 Query_1615341     2R0I_A   29.075             227        149        6
    1054 Query_1615341     4OLI_A   27.737             274        163        9
    1055 Query_1615341     1ZMU_A   29.075             227        149        6
    1056 Query_1615341     5EAK_A   29.439             214        139        6
    1057 Query_1615341     1WXM_A   60.563              71         28        0
    1058 Query_1615341     1RFA_A   54.667              75         31        1
    1059 Query_1615341     6JRK_A   29.057             265        171        7
    1060 Query_1615341     3UIU_A   27.797             295        163        8
    1061 Query_1615341     5DBX_A   26.596             282        182        8
    1062 Query_1615341     8A5J_A   27.830             212        145        5
    1063 Query_1615341     8QEL_B   28.889             270        155        8
    1064 Query_1615341     5BVK_A   25.228             329        180       12
    1065 Query_1615341     3ZOS_A   25.228             329        180       12
    1066 Query_1615341     5D9H_A   26.596             282        182        8
    1067 Query_1615341     1ZMW_A   28.634             227        150        6
    1068 Query_1615341     6S89_A   29.057             265        171        7
    1069 Query_1615341     1ZMV_A   28.634             227        150        6
    1070 Query_1615341     1C1Y_B   55.405              74         30        1
    1071 Query_1615341     4G0N_B   55.405              74         30        1
    1072 Query_1615341     5KZ7_A   29.439             214        139        6
    1073 Query_1615341     6JRJ_A   29.057             265        171        7
    1074 Query_1615341     6VJJ_B   55.405              74         30        1
    1075 Query_1615341     6Y23_A   25.228             329        180       12
    1076 Query_1615341     6BRJ_A   25.228             329        180       12
    1077 Query_1615341     1GUA_B   55.405              74         30        1
    1078 Query_1615341     5E8U_A   30.216             278        148       13
    1079 Query_1615341     6VC0_A   25.758             264        173       10
    1080 Query_1615341     5FDP_A   25.228             329        180       12
    1081 Query_1615341     8JOF_A   53.659              82         33        2
    1082 Query_1615341     4ASZ_A   27.609             297        175       12
    1083 Query_1615341     8TXY_A   28.505             214        141        6
    1084 Query_1615341     1RW8_A   30.686             277        148       13
    1085 Query_1615341     1B6C_B   30.686             277        148       13
    1086 Query_1615341     5USQ_A   30.686             277        148       13
    1087 Query_1615341     9F6X_A   30.686             277        148       13
    1088 Query_1615341     8YHF_A   30.686             277        148       13
    1089 Query_1615341     9J9D_A   30.686             277        148       13
    1090 Query_1615341     1VJY_A   30.686             277        148       13
    1091 Query_1615341     5DH3_A   27.897             233        154        7
    1092 Query_1615341     8A66_B   28.507             221        144        7
    1093 Query_1615341     3TZM_A   30.686             277        148       13
    1094 Query_1615341     1PY5_A   30.686             277        148       13
    1095 Query_1615341     4X0M_A   30.686             277        148       13
    1096 Query_1615341     5FRI_A   30.686             277        148       13
    1097 Query_1615341     8A66_A   28.761             226        137        8
    1098 Query_1615341     2WOT_A   30.686             277        148       13
    1099 Query_1615341     5I8A_A   27.986             293        156       12
    1100 Query_1615341     5KVT_A   27.986             293        156       12
    1101 Query_1615341     5E8T_A   30.686             277        148       13
    1102 Query_1615341     4YNZ_A   28.295             258        164        9
    1103 Query_1615341     4YOM_B   28.295             258        164        9
    1104 Query_1615341     8XFL_A   26.531             245        166        8
    1105 Query_1615341     1RRB_A   54.667              75         31        1
    1106 Query_1615341     5ES1_A   26.531             245        166        8
    1107 Query_1615341     3IEC_A   29.954             217        134        7
    1108 Query_1615341     5LPZ_A   28.239             301        184        9
    1109 Query_1615341     4M68_A   25.177             282        175       11
    1110 Query_1615341     6FER_A   26.174             298        173       11
    1111 Query_1615341     3KUD_B   54.054              74         31        1
    1112 Query_1615341     1Y8G_A   28.634             227        150        6
    1113 Query_1615341     3KUC_B   54.054              74         31        1
    1114 Query_1615341     3V5Q_A   26.689             296        172       10
    1115 Query_1615341     4BFM_A   31.220             205        130        8
    1116 Query_1615341     5LPV_A   28.309             272        166        7
    1117 Query_1615341     5LPB_A   28.309             272        166        7
    1118 Query_1615341     7AYM_A   25.378             331        177       13
    1119 Query_1615341     6N3N_A   30.667             225        115        8
    1120 Query_1615341     4OH4_A   28.309             272        166        7
    1121 Query_1615341     4CQG_A   31.220             205        130        8
    1122 Query_1615341     4LG4_A   28.319             226        138        8
    1123 Query_1615341     5XD6_A   30.556             216        130        6
    1124 Query_1615341     4YUR_A   26.357             258        174        7
    1125 Query_1615341     3MDY_A   30.357             280        147       14
    1126 Query_1615341     5WNO_A   27.027             259        174        6
    1127 Query_1615341     4XUF_A   28.980             245        133        9
    1128 Query_1615341     9Y9N_A   26.357             258        174        7
    1129 Query_1615341     6AO5_A   28.054             221        145        7
    1130 Query_1615341     3COK_A   26.357             258        174        7
    1131 Query_1615341     4JXF_A   26.357             258        174        7
    1132 Query_1615341     4LGD_A   28.054             221        145        7
    1133 Query_1615341     2OO8_X   29.352             293        160       15
    1134 Query_1615341     3TL8_A   31.797             217        132        8
    1135 Query_1615341     6JQR_A   28.980             245        133        9
    1136 Query_1615341     6OYW_A   25.670             261        171       10
    1137 Query_1615341     1FVR_A   29.110             292        162       15
    1138 Query_1615341     4BIB_A   25.692             253        165       10
    1139 Query_1615341     6VRE_A   26.087             253        164       10
    1140 Query_1615341     2BMC_A   24.573             293        182       10
    1141 Query_1615341     3UIM_A   31.797             217        132        8
    1142 Query_1615341     4YMJ_A   26.351             296        173       10
    1143 Query_1615341     1U59_A   25.856             263        175        8
    1144 Query_1615341     3COH_A   24.915             293        181       10
    1145 Query_1615341     3QBN_A   25.338             296        176       11
    1146 Query_1615341     4PRJ_A   24.915             293        181       10
    1147 Query_1615341     8XB1_A   28.980             245        133        9
    1148 Query_1615341     5X02_A   28.980             245        133        9
    1149 Query_1615341     6BFN_A   32.857             210        118        9
    1150 Query_1615341     1RJB_A   28.980             245        133        9
    1151 Query_1615341     3R21_A   24.662             296        184       10
    1152 Query_1615341     4BF2_A   25.670             261        171       10
    1153 Query_1615341     1GZK_A   28.400             250        156       10
    1154 Query_1615341     7T6F_A   25.267             281        179        9
    1155 Query_1615341     3D0E_A   28.400             250        156       10
    1156 Query_1615341     6N3L_A   30.222             225        116        8
    1157 Query_1615341     2XRU_A   24.573             293        182       10
    1158 Query_1615341     1O6K_A   28.400             250        156       10
    1159 Query_1615341     1GZN_A   28.400             250        156       10
    1160 Query_1615341     1MRV_A   28.400             250        156       10
    1161 Query_1615341     6E2M_A   25.670             261        171       10
    1162 Query_1615341     2JED_A   28.155             206        136        5
    1163 Query_1615341     5VIL_A   25.670             261        171       10
    1164 Query_1615341     6IL3_A   28.980             245        133        9
    1165 Query_1615341     2J4Z_A   24.573             293        182       10
    1166 Query_1615341     2JDO_A   28.400             250        156       10
    1167 Query_1615341     5F9E_A   28.155             206        136        5
    1168 Query_1615341     7ZTL_A   24.573             293        182       10
    1169 Query_1615341     3VW6_A   26.087             253        164       10
    1170 Query_1615341     2CLQ_A   25.670             261        171       10
    1171 Query_1615341     4Q9Z_A   28.155             206        136        5
    1172 Query_1615341     6HJJ_A   25.000             296        177       11
    1173 Query_1615341     6SEQ_A   28.175             252        150        9
    1174 Query_1615341     3UNZ_A   24.573             293        182       10
    1175 Query_1615341     3NRM_A   24.662             296        184       10
    1176 Query_1615341     1O6L_A   28.400             250        156       10
    1177 Query_1615341     6XIH_A   26.087             253        164       10
    1178 Query_1615341     5UOR_A   26.087             253        164       10
    1179 Query_1615341     8C1M_A   24.585             301        188       10
    1180 Query_1615341     4PX6_A   28.832             274        161       11
    1181 Query_1615341     5OBJ_A   24.585             301        188       10
    1182 Query_1615341     8Q61_A   28.400             250        156       10
    1183 Query_1615341     8GUW_A   25.000             296        177       11
    1184 Query_1615341     3EMG_A   28.832             274        161       11
    1185 Query_1615341     5DN3_A   24.585             301        188       10
    1186 Query_1615341     4J8M_A   24.573             293        182       10
    1187 Query_1615341     1MQ4_A   24.585             301        188       10
    1188 Query_1615341     2WQB_A   29.010             293        161       15
    1189 Query_1615341     5AAD_A   25.000             296        177       11
    1190 Query_1615341     4RX7_A   28.832             274        161       11
    1191 Query_1615341     2OZO_A   25.475             263        176        8
    1192 Query_1615341     5UOX_A   26.087             253        164       10
    1193 Query_1615341     4RX9_A   28.832             274        161       11
    1194 Query_1615341     6OYT_A   26.087             253        164       10
    1195 Query_1615341     3TUB_A   28.832             274        161       11
    1196 Query_1615341     9BZG_A   25.000             296        177       11
    1197 Query_1615341     1XJD_A   28.155             206        136        5
    1198 Query_1615341     4CEG_A   25.000             296        177       11
    1199 Query_1615341     4Q5J_A   32.178             202        123        5
    1200 Query_1615341     3TUC_A   28.832             274        161       11
    1201 Query_1615341     4YJO_A   28.832             274        161       11
    1202 Query_1615341     6J5T_A   27.018             285        171       11
    1203 Query_1615341     9C1W_A   28.400             250        156       10
    1204 Query_1615341     8OF5_A   25.000             296        177       11
    1205 Query_1615341     5LXM_A   25.000             296        177       11
    1206 Query_1615341     4X3J_A   29.010             293        161       15
    1207 Query_1615341     4K2R_A   25.475             263        176        8
    1208 Query_1615341     6VPJ_A   25.091             275        169       10
    1209 Query_1615341     8X5K_A   28.832             274        161       11
    1210 Query_1615341     2J50_A   24.573             293        182       10
    1211 Query_1615341     4F4P_A   28.832             274        161       11
    1212 Query_1615341     1MUO_A   24.573             293        182       10
    1213 Query_1615341     5ORL_A   24.573             293        182       10
    1214 Query_1615341     3SRV_B   28.832             274        161       11
    1215 Query_1615341     3FDN_A   24.573             293        182       10
    1216 Query_1615341     4BTF_A   25.000             280        178       11
    1217 Query_1615341     3W16_A   24.573             293        182       10
    1218 Query_1615341     4C3P_A   24.573             293        182       10
    1219 Query_1615341     5OS2_A   25.000             296        177       11
    1220 Query_1615341     7QQ6_A   30.804             224        113        9
    1221 Query_1615341     4RSS_A   28.832             274        161       11
    1222 Query_1615341     8SSP_A   24.573             293        182       10
    1223 Query_1615341     8FO7_C   30.435             276        165       11
    1224 Query_1615341     3W10_A   24.573             293        182       10
    1225 Query_1615341     8SMC_C   30.435             276        165       11
    1226 Query_1615341     5OSD_A   24.573             293        182       10
    1227 Query_1615341     7QWK_A   30.804             224        113        9
    1228 Query_1615341     8C1K_A   25.000             284        171       10
    1229 Query_1615341     5TR6_A   28.832             274        161       11
    1230 Query_1615341     6C83_A   24.573             293        182       10
    1231 Query_1615341     1XBA_A   28.832             274        161       11
    1232 Query_1615341     5ZAN_A   24.573             293        182       10
    1233 Query_1615341     6I2U_A   24.573             293        182       10
    1234 Query_1615341     5OS5_A   24.573             293        182       10
    1235 Query_1615341     4PV0_A   28.832             274        161       11
    1236 Query_1615341     2X6D_A   24.573             293        182       10
    1237 Query_1615341     8C15_A   25.000             284        171       10
    1238 Query_1615341     6VNO_A   30.576             278        162       13
    1239 Query_1615341     5C26_A   28.832             274        161       11
    1240 Query_1615341     8TXZ_A   30.576             278        162       13
    1241 Query_1615341     9DMI_A   30.576             278        162       13
    1242 Query_1615341     1OL5_A   24.573             293        182       10
    1243 Query_1615341     2WTW_A   24.573             293        182       10
    1244 Query_1615341     6VOV_A   28.832             274        161       11
    1245 Query_1615341     3LAU_A   24.573             293        182       10
    1246 Query_1615341     2WTV_A   24.573             293        182       10
    1247 Query_1615341     4YJQ_A   28.832             274        161       11
    1248 Query_1615341     2C6E_A   24.573             293        182       10
    1249 Query_1615341     6VP8_A   30.576             278        162       13
    1250 Query_1615341     2C6D_A   24.573             293        182       10
    1251 Query_1615341     2XNG_A   24.573             293        182       10
    1252 Query_1615341     5Y5T_A   28.832             274        161       11
    1253 Query_1615341     4FL3_A   28.832             274        161       11
    1254 Query_1615341     9KDS_A   25.000             296        177       11
    1255 Query_1615341     4DFL_A   28.832             274        161       11
    1256 Query_1615341     8RRQ_A   28.832             274        161       11
    1257 Query_1615341     8PS7_A   27.018             285        146       12
    1258 Query_1615341     7O2V_A   24.573             293        182       10
    1259 Query_1615341     4FL2_A   28.832             274        161       11
    1260 Query_1615341     4BN1_A   24.573             293        182       10
    1261 Query_1615341     8HOD_A   30.396             227        139        6
    1262 Query_1615341     4UZD_A   25.000             296        177       11
    1263 Query_1615341     3E5A_A   24.573             293        182       10
    1264 Query_1615341     7LHT_A   30.576             278        162       13
    1265 Query_1615341     8TZB_A   30.576             278        162       13
    1266 Query_1615341     3HA6_A   24.573             293        182       10
    1267 Query_1615341     4JAI_A   24.252             301        189       10
    1268 Query_1615341     8OKU_A   26.263             297        177       13
    1269 Query_1615341     5DOS_A   24.555             281        176        9
    1270 Query_1615341     3EFW_A   24.573             293        182       10
    1271 Query_1615341     5EW9_A   24.632             272        174        9
    1272 Query_1615341     3H0Y_A   24.573             293        182       10
    1273 Query_1615341     5DT4_A   24.555             281        176        9
    1274 Query_1615341     9CI3_A   30.576             278        162       13
    1275 Query_1615341     5DT3_A   24.555             281        176        9
    1276 Query_1615341     9D8Z_A   30.556             216        120        9
    1277 Query_1615341     8R4O_A   25.728             206        143        5
    1278 Query_1615341     6HM6_A   28.727             275        160       12
    1279 Query_1615341     3SRV_A   28.727             275        160       12
    1280 Query_1615341     6VPG_A   24.632             272        174        9
    1281 Query_1615341     3H9R_A   30.556             216        120        9
    1282 Query_1615341     9RDA_A   30.556             216        120        9
    1283 Query_1615341     6VPH_A   24.632             272        174        9
    1284 Query_1615341     6XKA_A   24.632             272        174        9
    1285 Query_1615341     9KS6_A   24.573             293        182       10
    1286 Query_1615341     9CHO_A   30.576             278        162       13
    1287 Query_1615341     6VPL_A   24.632             272        174        9
    1288 Query_1615341     9D8F_A   30.556             216        120        9
    1289 Query_1615341     6JUX_A   30.556             216        120        9
    1290 Query_1615341     8HO6_A   30.396             227        139        6
    1291 Query_1615341     4DYM_A   30.556             216        120        9
    1292 Query_1615341     5TOS_A   25.566             309        193       10
    1293 Query_1615341   9ESA_AAA   25.912             274        177        9
    1294 Query_1615341     6VPI_A   24.727             275        170       10
    1295 Query_1615341     8JF4_A   25.091             275        169       10
    1296 Query_1615341     9P6A_A   26.996             263        163       11
    1297 Query_1615341     8TZG_A   30.435             276        165       11
    1298 Query_1615341     3W2C_A   25.091             275        169       10
    1299 Query_1615341     2W1C_A   25.091             275        169       10
    1300 Query_1615341     8HOA_A   29.956             227        140        6
    1301 Query_1615341     2H6D_A   29.224             219        142        7
    1302 Query_1615341     2W1D_A   25.091             275        169       10
    1303 Query_1615341     6UNQ_A   30.556             216        120        9
    1304 Query_1615341     6UNR_A   30.556             216        120        9
    1305 Query_1615341     6GVX_A   30.244             205        132        8
    1306 Query_1615341     9F31_A   30.244             205        132        8
    1307 Query_1615341     4BKY_A   30.244             205        132        8
    1308 Query_1615341     5TWU_A   30.244             205        132        8
    1309 Query_1615341     7CTX_A   28.431             204        131        6
    1310 Query_1615341     6XR4_A   30.576             278        162       13
    1311 Query_1615341     6EIX_A   30.093             216        121        9
    1312 Query_1615341     5TWL_A   30.244             205        132        8
    1313 Query_1615341     2WQE_A   25.091             275        169       10
    1314 Query_1615341     3MY0_A   31.019             216        119        9
    1315 Query_1615341     2YZA_A   28.440             218        145        6
    1316 Query_1615341     6VXR_A   30.244             205        132        8
    1317 Query_1615341     3DAJ_A   24.080             299        184       11
    1318 Query_1615341     6KZC_A   26.129             310        170       11
    1319 Query_1615341     4ZHX_A   29.167             216        140        7
    1320 Query_1615341     5TVT_A   30.244             205        132        8
    1321 Query_1615341     5ISO_A   29.167             216        140        7
    1322 Query_1615341     4CFE_A   29.167             216        140        7
    1323 Query_1615341     2Y7J_A   25.357             280        176       13
    1324 Query_1615341     8TZC_A   30.216             278        163       13
    1325 Query_1615341     8BIK_A   29.167             216        140        7
    1326 Query_1615341     7KX8_A   26.471             238        152        9
    1327 Query_1615341     5IH8_A   30.244             205        132        8
    1328 Query_1615341     7F3G_A   26.471             238        152        9
    1329 Query_1615341     7KXW_A   26.471             238        152        9
    1330 Query_1615341     5JZN_A   26.471             238        152        9
    1331 Query_1615341     7LI3_A   30.216             278        163       13
    1332 Query_1615341     4M69_A   27.925             265        160        9
    1333 Query_1615341     4UMQ_A   30.049             203        131        8
    1334 Query_1615341     8TYQ_A   30.216             278        163       13
    1335 Query_1615341     3MTF_A   30.233             215        120        9
    1336 Query_1615341     6VBZ_A   24.275             276        185       11
    1337 Query_1615341     9IC2_A   29.167             216        140        7
    1338 Query_1615341     4UMT_A   30.049             203        131        8
    1339 Query_1615341     7OPO_A   28.270             237        146        9
    1340 Query_1615341     5EZV_A   29.167             216        140        7
    1341 Query_1615341     5D9K_A   28.270             237        146        9
    1342 Query_1615341     3G51_A   28.270             237        146        9
    1343 Query_1615341     5JZJ_A   26.471             238        152        9
    1344 Query_1615341     4QFG_A   29.167             216        140        7
    1345 Query_1615341     4NUS_A   28.270             237        146        9
    1346 Query_1615341     3D14_A   24.080             299        184       11
    1347 Query_1615341     6NPZ_A   28.287             251        156       11
    1348 Query_1615341     8XFY_A   28.270             237        146        9
    1349 Query_1615341     3CQU_A   28.016             257        149       12
    1350 Query_1615341     9S9T_A   29.851             201        129        6
    1351 Query_1615341     6BX6_A   28.372             215        143        6
    1352 Query_1615341     8UWR_A   30.233             215        120        9
    1353 Query_1615341     6KYQ_A   26.267             217        142        7
    1354 Query_1615341     4GV1_A   28.287             251        156       11
    1355 Query_1615341     4RER_A   30.189             212        131        9
    1356 Query_1615341     4REW_A   30.189             212        131        9
    1357 Query_1615341     3OCB_A   28.016             257        149       12
    1358 Query_1615341     6BUU_A   28.016             257        149       12
    1359 Query_1615341     2I0E_A   29.851             201        129        6
    1360 Query_1615341     4EL9_A   28.270             237        146        9
    1361 Query_1615341     6KYR_A   26.267             217        142        7
    1362 Query_1615341     6CCY_A   28.287             251        156       11
    1363 Query_1615341     3UBD_A   28.270             237        146        9
    1364 Query_1615341     4D2P_A   29.557             203        132        8
    1365 Query_1615341     2Z7Q_A   27.397             219        136        7
    1366 Query_1615341     3OHT_A   28.033             239        138        9
    1367 Query_1615341     4CFH_A   29.817             218        136        9
    1368 Query_1615341     4CZU_A   25.118             211        150        6
    1369 Query_1615341     6GR8_A   25.912             274        177        9
    1370 Query_1615341     4RED_A   30.189             212        131        9
    1371 Query_1615341     4CZT_A   25.118             211        150        6
    1372 Query_1615341     7JIJ_A   30.189             212        131        9
    1373 Query_1615341     5DE2_A   26.087             207        138        7
    1374 Query_1615341     8DFP_A   33.696             184        109        5
    1375 Query_1615341     8DFP_A   29.545             132         79        6
    1376 Query_1615341     7JHG_A   30.189             212        131        9
    1377 Query_1615341     3ZDU_A   25.839             298        197       12
    1378 Query_1615341     6C9F_A   30.189             212        131        9
    1379 Query_1615341     8DFQ_A   33.696             184        109        5
    1380 Query_1615341     8DFQ_A   29.545             132         79        6
    1381 Query_1615341     8DFM_A   33.696             184        109        5
    1382 Query_1615341     8DFM_A   29.545             132         79        6
    1383 Query_1615341     6NPY_B   26.087             207        138        7
    1384 Query_1615341     6C9H_A   30.189             212        131        9
    1385 Query_1615341     2DWB_A   24.232             293        183       10
    1386 Query_1615341     1ZYC_A   25.468             267        166        5
    1387 Query_1615341     8SXN_A   25.604             207        139        7
    1388 Query_1615341     1ZY4_A   25.468             267        166        5
    1389 Query_1615341     8E05_A   27.891             294        157       12
    1390 Query_1615341     2FH9_A   29.048             210        135        8
    1391 Query_1615341     3MN3_A   29.048             210        135        8
    1392 Query_1615341     3HYH_A   29.048             210        135        8
    1393 Query_1615341     3O96_A   29.435             248        157       11
    1394 Query_1615341     3MTL_A   28.972             214        140        7
    1395 Query_1615341     4EJN_A   29.435             248        157       11
    1396 Query_1615341     8E04_A   27.891             294        157       12
    1397 Query_1615341     3DAE_A   29.048             210        135        8
    1398 Query_1615341     8FAC_A   27.891             294        157       12
    1399 Query_1615341     6S73_A   26.087             207        138        7
    1400 Query_1615341     4M66_A   27.925             265        160        9
    1401 Query_1615341     9H59_C   25.604             207        139        7
    1402 Query_1615341     9NFQ_C   25.604             207        139        7
    1403 Query_1615341     8WS0_A   25.604             207        139        7
    1404 Query_1615341     4FSN_A   27.854             219        133        9
    1405 Query_1615341     4A4X_A   27.273             220        131        8
    1406 Query_1615341     1IA8_A   27.854             219        133        9
    1407 Query_1615341     4QYE_A   27.854             219        133        9
    1408 Query_1615341     5KCV_A   29.435             248        157       11
    1409 Query_1615341     7AKO_A   28.641             206        129        8
    1410 Query_1615341     2WQM_A   25.604             207        139        7
    1411 Query_1615341     7AKM_A   28.641             206        129        8
    1412 Query_1615341     1ZLT_A   27.854             219        133        9
    1413 Query_1615341     4FSM_A   27.854             219        133        9
    1414 Query_1615341     4FSY_A   27.854             219        133        9
    1415 Query_1615341     2X8E_A   27.854             219        133        9
    1416 Query_1615341     2BR1_A   27.854             219        133        9
    1417 Query_1615341     8QGY_A   25.692             253        165       10
    1418 Query_1615341     2HOG_A   27.982             218        132        9
    1419 Query_1615341     5OQ5_A   28.311             219        132        9
    1420 Query_1615341     2QHM_A   27.982             218        132        9
    1421 Query_1615341     7APJ_A   29.435             248        157       11
    1422 Query_1615341     3ORM_A   28.251             223        135        8
    1423 Query_1615341     4FSW_A   27.854             219        133        9
    1424 Query_1615341     5F4N_A   27.854             219        133        9
    1425 Query_1615341     4D28_A   26.923             208        142        7
    1426 Query_1615341     3JVR_A   27.854             219        133        9
    1427 Query_1615341     6TM5_S   26.941             219        133        7
    1428 Query_1615341     2ACX_A   27.586             232        134        7
    1429 Query_1615341     2E9V_A   27.854             219        133        9
    1430 Query_1615341     4FSZ_A   27.854             219        133        9
    1431 Query_1615341     3OT3_A   27.854             219        133        9
    1432 Query_1615341     8UW7_A   29.435             248        157       11
    1433 Query_1615341     3NYN_A   27.586             232        134        7
    1434 Query_1615341     3IW4_A   28.358             201        132        5
    1435 Query_1615341     2AYP_A   27.854             219        133        9
    1436 Query_1615341     8UVY_A   29.435             248        157       11
    1437 Query_1615341     6TM5_Q   26.941             219        133        7
    1438 Query_1615341     2E9N_A   27.854             219        133        9
    1439 Query_1615341     2YDJ_A   27.854             219        133        9
    1440 Query_1615341     4FT3_A   27.854             219        133        9
    1441 Query_1615341     6TUA_A   26.073             303        189       10
    1442 Query_1615341     8EJ4_K   25.118             211        143        7
    1443 Query_1615341     4FST_A   27.854             219        133        9
    1444 Query_1615341     1MRU_A   28.251             223        135        8
    1445 Query_1615341     2W5A_A   26.941             219        133        7
    1446 Query_1615341     4FR4_A   29.500             200        130        6
    1447 Query_1615341     2JAV_A   26.941             219        133        7
    1448 Query_1615341     8A3T_S   30.286             175        111        7
    1449 Query_1615341     5U94_A   28.251             223        135        8
    1450 Query_1615341     4AF3_A   25.589             297        186       11
    1451 Query_1615341     3F69_A   28.251             223        135        8
    1452 Query_1615341     6B2P_A   28.251             223        135        8
    1453 Query_1615341     4RA4_A   28.358             201        132        5
    1454 Query_1615341     6I2P_A   28.251             223        135        8
    1455 Query_1615341     7PUE_A   26.210             248        165        7
    1456 Query_1615341     3ORI_A   28.251             223        135        8
    1457 Query_1615341     1ZYS_A   27.854             219        133        9
    1458 Query_1615341     3F61_A   28.251             223        135        8
    1459 Query_1615341     4G31_A   29.224             219        121        7
    1460 Query_1615341     1ZXE_A   25.094             267        167        5
    1461 Query_1615341     1O6Y_A   28.251             223        135        8
    1462 Query_1615341     8WS1_A   25.121             207        140        7
    1463 Query_1615341     4FSU_A   27.854             219        133        9
    1464 Query_1615341     2R5T_A   26.210             248        165        7
    1465 Query_1615341     9IWX_A   27.547             265        161        9
    1466 Query_1615341     6G78_A   26.432             227        153        7
    1467 Query_1615341     5I3O_A   24.219             256        172        9
    1468 Query_1615341     4PNI_A   25.616             203        136        5
    1469 Query_1615341     4W9W_A   24.219             256        172        9
    1470 Query_1615341     8VSU_C   23.894             226        162        4
    1471 Query_1615341     4W9X_A   24.219             256        172        9
    1472 Query_1615341     3QC9_A   25.616             203        136        5
    1473 Query_1615341     3C4X_A   25.616             203        136        5
    1474 Query_1615341     3C4W_A   25.616             203        136        5
    1475 Query_1615341     3T8O_A   25.616             203        136        5
    1476 Query_1615341     6G76_A   26.432             227        153        7
    1477 Query_1615341     6G77_A   26.432             227        153        7
    1478 Query_1615341     7MT8_G   25.616             203        136        5
    1479 Query_1615341     4WBO_A   25.616             203        136        5
    1480 Query_1615341     6NTD_B   48.649              74         35        1
    1481 Query_1615341     4E5A_X   27.197             239        140        9
    1482 Query_1615341     3H4J_A   29.717             212        131        8
    1483 Query_1615341     4L9I_A   25.616             203        136        5
    1484 Query_1615341     2PUU_A   27.197             239        140        9
    1485 Query_1615341     3HNG_A   31.176             170        101        5
    1486 Query_1615341     2FSL_X   27.197             239        140        9
    1487 Query_1615341     3P4K_A   27.197             239        140        9
    1488 Query_1615341     8EFJ_A   27.197             239        140        9
    1489 Query_1615341     1LEW_A   27.197             239        140        9
    1490 Query_1615341     2FST_X   27.197             239        140        9
    1491 Query_1615341     4EQM_A   25.391             256        170       11
    1492 Query_1615341     4LOO_A   27.197             239        140        9
    1493 Query_1615341     3VHE_A   30.909             165        102        4
    1494 Query_1615341     3VHK_A   30.909             165        102        4
    1495 Query_1615341     2GTM_A   27.197             239        140        9
    1496 Query_1615341     3NNX_A   27.197             239        140        9
    1497 Query_1615341     1BMK_A   27.197             239        140        9
    1498 Query_1615341     5O90_A   27.197             239        140        9
    1499 Query_1615341     5MRD_A   26.087             207        144        5
    1500 Query_1615341     3TG1_A   27.197             239        140        9
    1501 Query_1615341     3VID_A   30.909             165        102        4
    1502 Query_1615341     2FSO_X   27.197             239        140        9
    1503 Query_1615341     1Y6A_A   30.909             165        102        4
    1504 Query_1615341     2BAQ_A   27.615             239        139        9
    1505 Query_1615341     4TYH_B   27.197             239        140        9
    1506 Query_1615341     8YPE_A   27.197             239        140        9
    1507 Query_1615341     5O8U_A   27.197             239        140        9
    1508 Query_1615341     5O8V_A   27.197             239        140        9
    1509 Query_1615341     3ODZ_X   27.197             239        140        9
    1510 Query_1615341   8ACM_AAA   27.197             239        140        9
    1511 Query_1615341     3PWY_A   27.230             213        144        6
    1512 Query_1615341     3HVC_A   27.197             239        140        9
    1513 Query_1615341     3S3I_A   27.197             239        140        9
    1514 Query_1615341     3OD6_X   27.197             239        140        9
    1515 Query_1615341     3NNU_A   27.197             239        140        9
    1516 Query_1615341     6SO1_A   27.197             239        140        9
    1517 Query_1615341     3ODY_X   27.197             239        140        9
    1518 Query_1615341     5NZZ_E   27.197             239        140        9
    1519 Query_1615341     3FI4_A   27.197             239        140        9
    1520 Query_1615341     3HEC_A   27.197             239        140        9
    1521 Query_1615341     6SOI_A   27.197             239        140        9
    1522 Query_1615341     1A9U_A   27.197             239        140        9
    1523 Query_1615341     9MHB_A   27.197             239        140        9
    1524 Query_1615341     1YW2_A   27.197             239        140        9
    1525 Query_1615341     2YIS_A   27.197             239        140        9
    1526 Query_1615341     1OZ1_A   27.197             239        140        9
    1527 Query_1615341     2GFS_A   27.197             239        140        9
    1528 Query_1615341     3K3I_A   27.197             239        140        9
    1529 Query_1615341     8X3M_A   27.197             239        140        9
    1530 Query_1615341     3GCU_A   27.197             239        140        9
    1531 Query_1615341     3OEF_X   27.197             239        140        9
    1532 Query_1615341     1YWR_A   27.197             239        140        9
    1533 Query_1615341     2GHL_A   27.197             239        140        9
    1534 Query_1615341     3KL8_A   27.317             205        135        7
    1535 Query_1615341     2BAJ_A   27.197             239        140        9
    1536 Query_1615341     3K3J_A   27.197             239        140        9
    1537 Query_1615341     3ZS5_A   27.197             239        140        9
    1538 Query_1615341     1DI9_A   27.197             239        140        9
    1539 Query_1615341     1ZZL_A   27.197             239        140        9
    1540 Query_1615341     2NPQ_A   27.197             239        140        9
    1541 Query_1615341     3KQ7_A   27.197             239        140        9
    1542 Query_1615341     3MPT_A   27.197             239        140        9
    1543 Query_1615341     8VXE_A   27.197             239        140        9
    1544 Query_1615341     5ETC_A   27.197             239        140        9
    1545 Query_1615341     6KA4_A   30.435             184        109        4
    1546 Query_1615341     3PY3_A   27.197             239        140        9
    1547 Query_1615341     2BAL_A   27.197             239        140        9
    1548 Query_1615341     6TCA_B   27.197             239        140        9
    1549 Query_1615341     2OZA_B   27.197             239        140        9
    1550 Query_1615341     3KK9_A   27.317             205        135        7
    1551 Query_1615341     1M7Q_A   27.197             239        140        9
    1552 Query_1615341     3E92_A   27.197             239        140        9
    1553 Query_1615341     3KK8_A   27.317             205        135        7
    1554 Query_1615341     2BDW_A   27.317             205        135        7
    1555 Query_1615341     6ZQS_A   27.197             239        140        9
    1556 Query_1615341     3D7Z_A   27.197             239        140        9
    1557 Query_1615341     9CJ1_A   27.197             239        140        9
    1558 Query_1615341     4F9W_A   27.197             239        140        9
    1559 Query_1615341     4X7H_A   28.378             222        122        7
    1560 Query_1615341     2LGC_A   27.197             239        140        9
    1561 Query_1615341     6Y7W_A   27.197             239        140        9
    1562 Query_1615341     1H1W_A   26.761             213        145        6
    1563 Query_1615341     6SOV_A   27.197             239        140        9
    1564 Query_1615341     4R3C_A   27.197             239        140        9
    1565 Query_1615341     3DT1_A   27.197             239        140        9
    1566 Query_1615341     3D83_A   27.197             239        140        9
    1567 Query_1615341     7PVU_A   27.197             239        140        9
    1568 Query_1615341   7Z6I_AAA   27.197             239        140        9
    1569 Query_1615341     9NYT_A   27.197             239        140        9
    1570 Query_1615341     4F9Y_A   27.197             239        140        9
    1571 Query_1615341     1IAN_A   27.197             239        140        9
    1572 Query_1615341     1UU9_A   26.761             213        145        6
    1573 Query_1615341     1V0O_A   26.190             210        145        6
    1574 Query_1615341     3ION_A   26.761             213        145        6
    1575 Query_1615341     1V0B_A   26.190             210        145        6
    1576 Query_1615341     2Y8O_A   27.197             239        140        9
    1577 Query_1615341     4CT1_A   26.087             207        144        5
    1578 Query_1615341     5ETF_A   26.778             239        141        9
    1579 Query_1615341     6QYX_A   27.848             237        147        8
    1580 Query_1615341     1FOT_A   27.982             218        137        7
    1581 Query_1615341     2R7B_A   26.761             213        145        6
    1582 Query_1615341     5WJJ_A   27.197             239        140        9
    1583 Query_1615341     3H9O_A   26.761             213        145        6
    1584 Query_1615341     4A07_A   26.570             207        143        5
    1585 Query_1615341     3NUS_A   26.761             213        145        6
    1586 Query_1615341     5L4Q_A   24.561             285        188       10
    1587 Query_1615341     3ORX_A   26.761             213        145        6
    1588 Query_1615341     1Z5M_A   26.761             213        145        6
    1589 Query_1615341     5MZ3_A   27.197             239        140        9
    1590 Query_1615341     9YLH_A   27.197             239        140        9
    1591 Query_1615341     5TE0_A   24.561             285        188       10
    1592 Query_1615341     3QC4_A   26.761             213        145        6
    1593 Query_1615341     3HRC_A   26.570             207        143        5
    1594 Query_1615341     5ETI_A   26.778             239        141        9
    1595 Query_1615341     2XCH_A   26.761             213        145        6
    1596 Query_1615341     1OKY_A   26.761             213        145        6
    1597 Query_1615341     1OVE_A   27.197             239        140        9
    1598 Query_1615341     7LVH_A   24.912             285        187       10
    1599 Query_1615341     4XX9_A   26.570             207        143        5
    1600 Query_1615341     4EWQ_A   27.197             239        140        9
    1601 Query_1615341     1OB3_A   26.190             210        145        6
    1602 Query_1615341     3Q4T_A   32.353             204        111        7
    1603 Query_1615341     6M9L_A   27.197             239        140        9
    1604 Query_1615341     3RWQ_A   26.761             213        145        6
    1605 Query_1615341     2BIY_A   26.761             213        145        6
    1606 Query_1615341     6M95_A   27.197             239        140        9
    1607 Query_1615341     3RWP_A   26.761             213        145        6
    1608 Query_1615341     4WSQ_A   24.561             285        188       10
    1609 Query_1615341     7XBR_F   26.733             202        134        7
    1610 Query_1615341     3MH2_A   26.778             239        141        9
    1611 Query_1615341     2WTK_C   23.451             226        163        4
    1612 Query_1615341     4IC7_A   23.881             268        178        8
    1613 Query_1615341     3HKO_A   23.810             294        167       12
    1614 Query_1615341     8A8M_A   27.197             239        140        9
    1615 Query_1615341     9QB5_A   24.561             285        188       10
    1616 Query_1615341     4KIK_A   26.697             221        132        9
    1617 Query_1615341     3NUN_A   26.761             213        145        6
    1618 Query_1615341     4B99_A   23.881             268        178        8
    1619 Query_1615341     2PHK_A   29.412             204        121       10
    1620 Query_1615341     2QLU_A   28.700             223        127        8
    1621 Query_1615341     4GEO_A   26.778             239        141        9
    1622 Query_1615341     9F58_A   28.351             194        120        4
    1623 Query_1615341     3GCP_A   26.778             239        141        9
    1624 Query_1615341     1PHK_A   29.412             204        121       10
    1625 Query_1615341     9CMZ_A   27.119             236        124       12
    1626 Query_1615341     8GMC_A   24.561             285        188       10
    1627 Query_1615341     8OMV_A   26.244             221        133        9
    1628 Query_1615341     4E3C_A   26.244             221        133        9
    1629 Query_1615341     8U2O_A   27.119             236        124       12
    1630 Query_1615341     3NAX_A   26.761             213        145        6
    1631 Query_1615341     2XCK_A   26.291             213        146        6
    1632 Query_1615341     9BF3_A   30.070             143         86        2
    1633 Query_1615341     3NAY_A   26.761             213        145        6
    1634 Query_1615341     2W96_B   26.818             220        133        9
    1635 Query_1615341     3V3V_A   25.676             296        155       12
    1636 Query_1615341     3WE4_A   25.806             248        155       10
    1637 Query_1615341     3A60_A   25.806             248        155       10
    1638 Query_1615341     3GP0_A   26.522             230        140        9
    1639 Query_1615341     3A62_A   25.806             248        155       10
    1640 Query_1615341     4KIK_B   26.244             221        133        9
    1641 Query_1615341     4L45_A   26.210             248        154       10
    1642 Query_1615341     4L46_A   26.210             248        154       10
    1643 Query_1615341     4L43_A   26.210             248        154       10
    1644 Query_1615341     7N91_A   26.210             248        154       10
    1645 Query_1615341     8YGW_A   26.522             230        140        9
    1646 Query_1615341     7XBR_A   26.238             202        135        7
    1647 Query_1615341     1QL6_A   29.412             204        121       10
    1648 Query_1615341     4ZSJ_A   24.793             242        151        7
    1649 Query_1615341     1CM8_A   27.928             222        126        9
    1650 Query_1615341     4L42_A   25.806             248        155       10
    1651 Query_1615341     5BYY_A   24.793             242        151        7
    1652 Query_1615341     6P8E_B   26.818             220        133        9
    1653 Query_1615341     2W99_B   26.818             220        133        9
    1654 Query_1615341     5BYZ_A   24.793             242        151        7
    1655 Query_1615341     6HKM_A   24.793             242        151        7
    1656 Query_1615341     4ZSG_A   24.793             242        151        7
    1657 Query_1615341     8T7T_A   34.043             141         75        4
    1658 Query_1615341     8T7T_A   28.906             128         79        3
    1659 Query_1615341     4L3J_A   27.570             214        136        9
    1660 Query_1615341     5O7I_A   24.793             242        151        7
    1661 Query_1615341     4Y83_A   25.403             248        166       10
    1662 Query_1615341     4AGU_A   25.592             211        148        5
    1663 Query_1615341     6UNA_A   27.727             220        125        9
    1664 Query_1615341     4C57_A   26.016             246        148       11
    1665 Query_1615341     5LOH_A   25.472             212        143        6
    1666 Query_1615341     9CSK_B   26.818             220        133        9
    1667 Query_1615341     7CGA_A   27.928             222        126        9
    1668 Query_1615341     6V6A_A   24.354             271        172        9
    1669 Query_1615341     4RLO_A   27.570             214        136        9
    1670 Query_1615341     3GC8_A   26.522             230        140        9
    1671 Query_1615341     2CN5_A   25.221             226        149        6
    1672 Query_1615341     9LTA_A   25.207             242        150        8
    1673 Query_1615341     2W9F_B   26.818             220        133        9
    1674 Query_1615341     2XK9_A   25.221             226        149        6
    1675 Query_1615341     6LBA_A   29.843             191        109        5
    1676 Query_1615341     5Y8U_A   31.818             154         97        4
    1677 Query_1615341     2YCF_A   25.221             226        149        6
    1678 Query_1615341     2YCR_A   25.221             226        149        6
    1679 Query_1615341     2BFX_A   22.172             221        152        6
    1680 Query_1615341     5Z1D_A   31.818             154         97        4
    1681 Query_1615341     4C2W_A   22.172             221        152        6
    1682 Query_1615341     2W0J_A   25.221             226        149        6
    1683 Query_1615341     4B8M_A   22.172             221        152        6
    1684 Query_1615341     2VRX_A   22.172             221        152        6
    1685 Query_1615341     3TXO_A   28.019             207        137        7
    1686 Query_1615341     4C2V_A   22.172             221        152        6
    1687 Query_1615341     8FP1_A   28.019             207        137        7
    1688 Query_1615341     3COI_A   25.847             236        133       10
    1689 Query_1615341     3GC9_A   26.522             230        140        9
    1690 Query_1615341     4B8L_A   22.172             221        152        6
    1691 Query_1615341     6IB0_A   31.818             154         97        4
    1692 Query_1615341     2DYL_A   31.169             154         98        4
    1693 Query_1615341     5Y90_A   31.818             154         97        4
    1694 Query_1615341     7OVJ_A   31.818             154         97        4
    1695 Query_1615341     2XRW_A   25.161             310        163       13
    1696 Query_1615341     3DLS_A   26.244             221        124        7
    1697 Query_1615341     3NIZ_A   27.273             209        142        6
    1698 Query_1615341     4ZSL_A   24.686             239        149        7
    1699 Query_1615341     6YG4_A   31.818             154         97        4
    1700 Query_1615341     6YFZ_A   31.818             154         97        4
    1701 Query_1615341     6HKN_A   24.561             228        148        6
    1702 Query_1615341     5B2L_A   31.818             154         97        4
    1703 Query_1615341     3ALN_A   26.882             279        176       11
    1704 Query_1615341     3WZU_A   31.818             154         97        4
    1705 Query_1615341     3J4Q_D   27.572             243        152       11
    1706 Query_1615341     3FHI_A   27.572             243        152       11
    1707 Query_1615341     8V5H_A   25.463             216        141        7
    1708 Query_1615341     6YG1_A   31.169             154         98        4
    1709 Query_1615341     4EYJ_A   25.847             236        133       10
    1710 Query_1615341     7DV6_A   27.803             223        126       10
    1711 Query_1615341     6ZR5_A   25.338             296        156       12
    1712 Query_1615341     3G33_A   26.457             223        133        9
    1713 Query_1615341     2QUR_A   27.572             243        152       11
    1714 Query_1615341     2QKR_A   27.273             209        142        6
    1715 Query_1615341     8X23_A   25.847             236        133       10
    1716 Query_1615341     2XS0_A   26.897             290        159       13
    1717 Query_1615341     8YGZ_A   27.803             223        126       10
    1718 Query_1615341     7SJ3_A   26.457             223        133        9
    1719 Query_1615341     5Y7Z_A   25.532             235        145       10
    1720 Query_1615341     5FWK_K   26.457             223        133        9
    1721 Query_1615341     1UKH_A   25.338             296        156       12
    1722 Query_1615341     6YG0_A   31.169             154         98        4
    1723 Query_1615341     8X5M_A   25.338             296        156       12
    1724 Query_1615341     4YR8_A   25.338             296        156       12
    1725 Query_1615341     3QAL_E   27.572             243        152       11
    1726 Query_1615341     1ATP_E   27.572             243        152       11
    1727 Query_1615341     5E8V_A   27.803             223        126       10
    1728 Query_1615341     2ERZ_E   27.572             243        152       11
    1729 Query_1615341     7E0Z_A   27.572             243        152       11
    1730 Query_1615341     1J3H_A   27.572             243        152       11
    1731 Query_1615341     3QA8_A   27.149             221        131       10
    1732 Query_1615341     4YHJ_A   27.273             165        112        4
    1733 Query_1615341     3O17_A   25.753             299        151       13
    1734 Query_1615341     2CPK_E   27.572             243        152       11
    1735 Query_1615341     3RZF_A   27.149             221        131       10
    1736 Query_1615341     4IAC_A   27.572             243        152       11
    1737 Query_1615341     8R99_A   29.557             203        115        9
    1738 Query_1615341     4DG3_E   27.572             243        152       11
    1739 Query_1615341     3X2U_A   27.572             243        152       11
    1740 Query_1615341     4MYG_A   25.847             236        133       10
    1741 Query_1615341     9GLA_A   26.364             220        151        7
    1742 Query_1615341     2G01_A   25.753             299        151       13
    1743 Query_1615341     6CCF_A   26.364             220        149        5
    1744 Query_1615341     8X5L_A   26.939             245        151       10
    1745 Query_1615341     2GNJ_A   27.049             244        152       10
    1746 Query_1615341     3I6U_A   23.667             300        201        9
    1747 Query_1615341     6CMJ_A   25.503             298        177       10
    1748 Query_1615341     1VZO_A   25.481             208        136        7
    1749 Query_1615341     3O7L_B   27.572             243        152       11
    1750 Query_1615341     2GU8_A   28.509             228        139       11
    1751 Query_1615341     1ZRZ_A   25.263             190        129        5
    1752 Query_1615341     3QD2_B   34.862             109         54        3
    1753 Query_1615341     3A8W_A   26.207             145        101        3
    1754 Query_1615341     8R3X_A   26.207             145        101        3
    1755 Query_1615341     3I6W_A   23.667             300        201        9
    1756 Query_1615341     3L9M_A   27.160             243        153       10
    1757 Query_1615341     5LI1_A   26.207             145        101        3
    1758 Query_1615341     2GNF_A   27.049             244        152       10
    1759 Query_1615341     4DC2_A   26.207             145        101        3
    1760 Query_1615341     2PK9_A   25.822             213        145        8
    1761 Query_1615341     8RU8_A   25.551             227        148        9
    1762 Query_1615341     5LI9_A   26.207             145        101        3
    1763 Query_1615341     1CMK_E   27.869             244        150       11
    1764 Query_1615341     1CTP_E   27.869             244        150       11
    1765 Query_1615341     5LIH_A   26.207             145        101        3
    1766 Query_1615341     6FRX_A   26.939             245        151       10
    1767 Query_1615341     1CDK_A   27.869             244        150       11
    1768 Query_1615341     3ZH8_A   26.207             145        101        3
    1769 Query_1615341     3RP9_A   27.459             244        124       11
    1770 Query_1615341     8R9B_A   29.557             203        115        9
    1771 Query_1615341     4NTS_A   27.500             240        150       11
    1772 Query_1615341     9EJK_B   26.207             145        101        3
    1773 Query_1615341     9EJK_B   30.000              50         31        1
    1774 Query_1615341     2V7O_A   26.087             207        132        7
    1775 Query_1615341     4DFX_E   27.500             240        150       11
    1776 Query_1615341     4WNK_A   29.560             159        104        4
    1777 Query_1615341     8R3Y_I   26.207             145        101        3
    1778 Query_1615341     4AE9_A   26.971             241        152       10
    1779 Query_1615341     1SYK_A   27.160             243        153       11
    1780 Query_1615341     10BL_A   26.339             224        147        7
    1781 Query_1615341     4AE6_A   26.971             241        152       10
    1782 Query_1615341     10SL_A   26.339             224        147        7
    1783 Query_1615341     6IB2_A   31.169             154         98        4
    1784 Query_1615341     9HIX_J   29.557             203        115        9
    1785 Query_1615341     5UV4_A   24.413             213        150        6
    1786 Query_1615341     2BFY_A   21.560             218        153        6
    1787 Query_1615341     5OO1_A   25.225             222        155        7
    1788 Query_1615341     4DFY_A   27.160             243        153       11
    1789 Query_1615341     8P7L_J   29.557             203        115        9
    1790 Query_1615341     1JBP_E   28.070             228        140       11
    1791 Query_1615341     8JFK_C   26.961             204        126       10
    1792 Query_1615341     1BKX_A   28.070             228        140       11
    1793 Query_1615341     1L3R_E   28.070             228        140       11
    1794 Query_1615341     5EYK_A   21.560             218        153        6
    1795 Query_1615341     3FI3_A   25.830             271        140       12
    1796 Query_1615341     1APM_E   28.070             228        140       11
    1797 Query_1615341     5X3F_B   28.070             228        140       11
    1798 Query_1615341     8ORM_J   29.557             203        115        9
    1799 Query_1615341     5K3Y_A   21.560             218        153        6
    1800 Query_1615341     2WEL_A   26.087             207        132        7
    1801 Query_1615341     6MM5_E   28.070             228        140       11
    1802 Query_1615341     7UJR_A   25.000             212        132        8
    1803 Query_1615341     8UKP_E   28.070             228        140       11
    1804 Query_1615341     3PVB_A   28.070             228        140       11
    1805 Query_1615341     7E11_A   28.070             228        140       11
    1806 Query_1615341     8UKN_C   28.070             228        140       11
    1807 Query_1615341     1UA2_A   29.557             203        115        9
    1808 Query_1615341     6O9L_8   29.557             203        115        9
    1809 Query_1615341     5VLO_A   25.714             210        129        8
    1810 Query_1615341     1RDQ_E   27.572             243        152       11
    1811 Query_1615341     1XH9_A   26.829             246        150       10
    1812 Query_1615341     3NX8_A   26.749             243        154       10
    1813 Query_1615341     4WB5_A   26.749             243        154       10
    1814 Query_1615341     3TTI_A   25.830             271        140       12
    1815 Query_1615341     3MVJ_A   26.749             243        154       10
    1816 Query_1615341     6C0U_A   26.749             243        154       10
    1817 Query_1615341     8P4Z_A   29.557             203        115        9
    1818 Query_1615341     3AGM_A   26.749             243        154       10
    1819 Query_1615341     3QAM_E   27.160             243        153       11
    1820 Query_1615341     4ERW_A   25.225             222        155        7
    1821 Query_1615341     1Q61_A   26.639             244        153       10
    1822 Query_1615341     4WHZ_A   25.830             271        140       12
    1823 Query_1615341     3KVX_A   25.830             271        140       12
    1824 Query_1615341     7B55_B   25.604             207        133        7
    1825 Query_1615341     4O21_A   28.070             228        140       11
    1826 Query_1615341     9BLH_A   26.087             207        132        7
    1827 Query_1615341     2VN9_A   26.087             207        132        7
    1828 Query_1615341     1SZM_A   26.639             244        153       10
    1829 Query_1615341     1H4L_A   24.091             220        158        5
    1830 Query_1615341     4AU8_A   24.186             215        154        5
    1831 Query_1615341     2GNG_A   26.639             244        153       10
    1832 Query_1615341     2R9S_A   25.830             271        140       12
    1833 Query_1615341     2IW6_A   25.446             224        152        9
    1834 Query_1615341     4Z84_A   26.639             244        153       10
    1835 Query_1615341     3VUL_A   26.667             240        130        9
    1836 Query_1615341     2VZ6_A   25.604             207        133        7
    1837 Query_1615341     7VDP_A   24.186             215        154        5
    1838 Query_1615341     8JPB_G   26.147             218        142        6
    1839 Query_1615341     7B5O_J   29.557             203        115        9
    1840 Query_1615341     3VUM_A   26.667             240        130        9
    1841 Query_1615341     4CFU_A   25.225             222        155        7
    1842 Query_1615341     3AMA_A   26.531             245        152       10
    1843 Query_1615341     3ELJ_A   25.418             299        152       13
    1844 Query_1615341     3PTG_A   25.830             271        140       12
    1845 Query_1615341     6Q4G_A   25.225             222        155        7
    1846 Query_1615341     1GZ8_A   25.225             222        155        7
    1847 Query_1615341     3FV8_A   25.830             271        140       12
    1848 Query_1615341     7ORE_A   25.830             271        140       12
    1849 Query_1615341     1PMN_A   25.830             271        140       12
    1850 Query_1615341     3PZE_A   25.418             299        152       13
    1851 Query_1615341     5N23_A   26.423             246        152       10
    1852 Query_1615341     3OXI_A   25.830             271        140       12
    1853 Query_1615341     1GII_A   25.974             231        142       10
    1854 Query_1615341     1JNK_A   25.830             271        140       12
    1855 Query_1615341     4X21_A   25.830             271        140       12
    1856 Query_1615341     4EOS_A   25.225             222        155        7
    1857 Query_1615341     4UX9_A   25.418             299        152       13
    1858 Query_1615341     1OGU_A   25.225             222        155        7
    1859 Query_1615341     3FI2_A   25.830             271        140       12
    1860 Query_1615341     4W4V_A   25.830             271        140       12
    1861 Query_1615341     2O0U_A   25.830             271        140       12
    1862 Query_1615341     2OK1_A   25.830             271        140       12
    1863 Query_1615341     4QTD_A   25.418             299        152       13
    1864 Query_1615341     9FT9_A   25.418             299        152       13
    1865 Query_1615341     4EOP_A   25.225             222        155        7
    1866 Query_1615341     4BCM_A   25.225             222        155        7
    1867 Query_1615341     4EOQ_A   25.225             222        155        7
    1868 Query_1615341     1VYW_A   25.225             222        155        7
    1869 Query_1615341     2B1P_A   25.830             271        140       12
    1870 Query_1615341     4EOM_A   25.446             224        152        9
    1871 Query_1615341     3PXF_A   25.225             222        155        7
    1872 Query_1615341     1OIT_A   25.225             222        155        7
    1873 Query_1615341     4X3F_A   25.333             225        131        8
    1874 Query_1615341     9CKO_A   28.931             159        105        4
    1875 Query_1615341     5OO0_A   25.225             222        155        7
    1876 Query_1615341     2F7E_E   26.230             244        154       10
    1877 Query_1615341     3VUD_A   26.667             240        130        9
    1878 Query_1615341     2EXC_X   25.830             271        140       12
    1879 Query_1615341     6W4O_A   25.238             210        130        8
    1880 Query_1615341     2IW8_A   25.217             230        145        9
    1881 Query_1615341     4EOO_A   25.225             222        155        7
    1882 Query_1615341     6WJF_A   25.806             248        160       10
    1883 Query_1615341     5LW1_B   25.418             299        152       13
    1884 Query_1615341     4O38_A   26.016             246        148       11
    1885 Query_1615341     4EOJ_A   25.446             224        152        9
    1886 Query_1615341     4X3F_C   25.806             217        124        8
    1887 Query_1615341     8H6P_A   25.225             222        155        7
    1888 Query_1615341     4EON_A   25.446             224        152        9
    1889 Query_1615341     9I9J_K   25.225             222        155        7
    1890 Query_1615341     7E34_A   25.225             222        155        7
    1891 Query_1615341     1H1P_A   25.225             222        155        7
    1892 Query_1615341     1W98_A   25.225             222        155        7
    1893 Query_1615341     4EOK_A   25.446             224        152        9
    1894 Query_1615341     8USO_A   25.714             210        129        8
    1895 Query_1615341     8UV0_A   25.225             222        155        7
    1896 Query_1615341     7NVQ_A   25.225             222        155        7
    1897 Query_1615341     9DC6_A   25.806             248        160       10
    1898 Query_1615341     9NFS_A   25.806             248        160       10
    1899 Query_1615341     3EZR_A   25.225             222        155        7
    1900 Query_1615341     6GUE_A   25.225             222        155        7
    1901 Query_1615341     4OW8_A   25.333             225        131        8
    1902 Query_1615341     8FEC_B   25.806             248        160       10
    1903 Query_1615341     3BHT_A   25.225             222        155        7
    1904 Query_1615341     4TNB_A   28.931             159        105        4
    1905 Query_1615341     6PJX_A   28.931             159        105        4
    1906 Query_1615341     4WB7_A   25.806             248        160       10
    1907 Query_1615341     9HIU_B   25.225             222        155        7
    1908 Query_1615341     4EOI_A   25.217             230        145        9
    1909 Query_1615341     6INL_A   25.225             222        155        7
    1910 Query_1615341     2JGZ_A   25.225             222        155        7
    1911 Query_1615341     5N3N_A   27.160             243        153       11
    1912 Query_1615341     4I3Z_A   25.225             222        155        7
    1913 Query_1615341     1GY3_A   25.225             222        155        7
    1914 Query_1615341     5UQ1_A   25.225             222        155        7
    1915 Query_1615341     4L7F_A   25.418             299        152       13
    1916 Query_1615341     8BZO_A   25.225             222        155        7
    1917 Query_1615341     1E9H_A   25.225             222        155        7
    1918 Query_1615341     9NYQ_A   25.225             222        155        7
    1919 Query_1615341     2JDT_A   26.639             244        153       10
    1920 Query_1615341     6XD3_J   29.064             203        116        9
    1921 Query_1615341     5U6Y_A   25.238             210        130        8
    1922 Query_1615341     5K4J_A   25.225             222        155        7
    1923 Query_1615341     3PJ8_A   25.225             222        155        7
    1924 Query_1615341     1YDR_E   26.230             244        154       10
    1925 Query_1615341     6F14_A   27.160             243        153       11
    1926 Query_1615341     2UVY_A   26.639             244        153       10
    1927 Query_1615341     6Q4I_A   25.225             222        155        7
    1928 Query_1615341     1AQ1_A   25.225             222        155        7
    1929 Query_1615341     1FQ1_B   25.225             222        155        7
    1930 Query_1615341     1B38_A   25.225             222        155        7
    1931 Query_1615341     6OQI_A   25.225             222        155        7
    1932 Query_1615341     1SMH_A   26.423             246        151       10
    1933 Query_1615341     4RT7_A   29.697             165         96        6
    1934 Query_1615341     4RT7_A   27.523             109         68        4
    1935 Query_1615341     6Y05_A   27.160             243        153       11
    1936 Query_1615341     9OB2_A   25.225             222        155        7
    1937 Query_1615341     9PE7_A   24.000             225        152        5
    1938 Query_1615341     8PYR_A   29.064             203        116        9
    1939 Query_1615341     6EM7_A   27.160             243        153       11
    1940 Query_1615341     1SVH_A   26.230             244        154       10
    1941 Query_1615341     5VI9_A   26.230             244        154       10
    1942 Query_1615341     3QHR_A   25.225             222        155        7
    1943 Query_1615341     3G2F_A   27.523             218        117        9
    1944 Query_1615341     4X3F_B   31.250             112         73        2
    1945 Query_1615341     3DND_A   26.230             244        154       10
    1946 Query_1615341     5UUU_A   25.688             218        143        6
    1947 Query_1615341     1Q24_A   26.230             244        154       10
    1948 Query_1615341     1Q8W_A   26.230             244        154       10
    1949 Query_1615341     3AGL_A   26.337             243        155       10
    1950 Query_1615341     2C1A_A   26.230             244        154       10
    1951 Query_1615341     6VZK_A   24.528             212        133        8
    1952 Query_1615341     3SV0_A   26.389             216        140        6
    1953 Query_1615341     3VUG_A   26.667             240        130        9
    1954 Query_1615341     2JDS_A   26.230             244        154       10
    1955 Query_1615341     3VUK_A   26.667             240        130        9
    1956 Query_1615341     6UNP_A   27.523             218        117        9
    1957 Query_1615341     8SF8_A   26.230             244        154       10
    1958 Query_1615341    8GXQ_HI   29.064             203        116        9
    1959 Query_1615341     8YNG_A   23.982             221        159        6
    1960 Query_1615341     6XBZ_J   29.064             203        116        9
    1961 Query_1615341     1BI7_A   23.556             225        153        5
    1962 Query_1615341     3VUI_A   26.667             240        130        9
    1963 Query_1615341     3VUH_A   26.667             240        130        9
    1964 Query_1615341     4C33_A   26.230             244        154       10
    1965 Query_1615341     3KRW_A   25.688             218        143        6
    1966 Query_1615341     4WB8_A   27.193             228        142       10
    1967 Query_1615341     1JOW_B   23.556             225        153        5
    1968 Query_1615341     7PWD_A   25.688             218        143        6
    1969 Query_1615341     6TD3_B   25.991             227        149        8
    1970 Query_1615341     3CIK_A   25.688             218        143        6
    1971 Query_1615341     5NW8_A   28.636             220        133       11
    1972 Query_1615341     6C2Y_A   25.688             218        143        6
    1973 Query_1615341     8BYA_A   25.225             222        155        7
    1974 Query_1615341     5DYK_A   28.571             273        154       14
    1975 Query_1615341     5N1F_A   28.636             220        133       11
    1976 Query_1615341     1OMW_A   25.688             218        143        6
    1977 Query_1615341     3PSC_A   25.688             218        143        6
    1978 Query_1615341     4C0T_A   24.597             248        137       11
    1979 Query_1615341     4WIH_A   28.636             220        133       11
    1980 Query_1615341     6B2Q_A   25.333             225        131        8
    1981 Query_1615341     1UNG_A   23.636             220        159        5
    1982 Query_1615341     8I0M_A   23.556             225        153        5
    1983 Query_1615341     5UKK_A   25.688             218        143        6
    1984 Query_1615341     5HE1_A   25.688             218        143        6
    1985 Query_1615341     3NUP_A   23.556             225        153        5
    1986 Query_1615341     4MK0_A   25.688             218        143        6
    1987 Query_1615341     9I9I_K   26.540             211        142        6
    1988 Query_1615341     5MHQ_A   24.775             222        156        7
    1989 Query_1615341     2UZT_A   26.638             229        142       10
    1990 Query_1615341     5HE3_A   25.688             218        143        6
    1991 Query_1615341     5HE0_A   25.688             218        143        6
    1992 Query_1615341     4IZ5_A   26.471             238        150        9
    1993 Query_1615341     2VO0_A   27.074             229        141       10
    1994 Query_1615341     4WB6_B   26.337             243        155       10
    1995 Query_1615341     4WB6_A   26.337             243        155       10
    1996 Query_1615341     5EFQ_A   26.471             204        134        6
    1997 Query_1615341     2CJM_A   25.225             222        155        7
    1998 Query_1615341     5YV8_A   31.126             151         94        4
    1999 Query_1615341     2ZV2_A   31.126             151         94        4
    2000 Query_1615341     1Q8T_A   25.820             244        155       10
    2001 Query_1615341     6RFP_A   26.471             238        150        9
    2002 Query_1615341     7VDU_A   24.775             222        156        7
    2003 Query_1615341     1XH7_A   25.820             244        155       10
    2004 Query_1615341     3ZO2_A   25.820             244        155       10
    2005 Query_1615341     1H01_A   24.775             222        156        7
    2006 Query_1615341     1OIR_A   24.775             222        156        7
    2007 Query_1615341     3VN9_A   29.747             158        103        5
    2008 Query_1615341     8S79_A   24.380             242        130        7
    2009 Query_1615341     8A8M_B   23.839             323        219       10
    2010 Query_1615341     8XEY_A   24.576             236        139       11
    2011 Query_1615341     5UY6_A   31.126             151         94        4
    2012 Query_1615341     4C34_A   26.638             229        142       10
    2013 Query_1615341     2AC5_A   30.769             143         88        5
    2014 Query_1615341     4XRL_A   26.471             238        150        9
    2015 Query_1615341     2QR8_A   24.576             236        139       11
    2016 Query_1615341     6OQL_A   23.611             216        146        5
    2017 Query_1615341     4JG6_A   24.576             236        139       11
    2018 Query_1615341     3RNY_A   26.033             242        151       12
    2019 Query_1615341     4BCF_A   23.529             204        139        6
    2020 Query_1615341     8WF4_A   26.293             232        143       12
    2021 Query_1615341     4NIF_A   25.941             239        135       13
    2022 Query_1615341     9HW6_B   26.180             233        134        9
    2023 Query_1615341    9IJJ_4Z   26.147             218        142        6
    2024 Query_1615341     8I0L_A   23.529             204        139        6
    2025 Query_1615341     4IZA_A   26.471             238        150        9
    2026 Query_1615341     4S2Z_A   26.471             238        150        9
    2027 Query_1615341     5O1S_A   24.576             236        139       11
    2028 Query_1615341     3KN5_A   25.604             207        125        9
    2029 Query_1615341     6OPG_A   26.471             238        150        9
    2030 Query_1615341     5WP1_A   26.471             238        150        9
    2031 Query_1615341     7UKZ_A   26.415             212        143        6
    2032 Query_1615341     2BHH_A   24.775             222        156        7
    2033 Query_1615341     3NPC_A   27.468             233        137        9
    2034 Query_1615341     4EC8_A   23.529             204        139        6
    2035 Query_1615341     7N8T_A   27.468             233        137        9
    2036 Query_1615341     9HVX_A   26.180             233        134        9
    2037 Query_1615341     2WNT_A   25.941             239        135       13
    2038 Query_1615341     5V62_A   26.471             238        150        9
    2039 Query_1615341     3E7O_A   27.468             233        137        9
    2040 Query_1615341     7OPM_A   26.471             238        150        9
    2041 Query_1615341     3BLH_A   23.529             204        139        6
    2042 Query_1615341     6U2G_A   26.394             269        172       11
    2043 Query_1615341     6OPK_A   26.471             238        150        9
    2044 Query_1615341     7MFD_B   26.394             269        172       11
    2045 Query_1615341     4RZ7_A   37.500             104         54        4
    2046 Query_1615341     2QR7_A   24.891             229        147       10
    2047 Query_1615341     3ZUV_A   26.471             238        150        9
    2048 Query_1615341     5KKR_C   26.394             269        172       11
    2049 Query_1615341     2ERK_A   26.471             238        150        9
    2050 Query_1615341     4QP1_A   26.471             238        150        9
    2051 Query_1615341     5EZR_A   37.500             104         54        4
    2052 Query_1615341     4IZ7_A   26.471             238        150        9
    2053 Query_1615341     7CML_A   27.468             233        137        9
    2054 Query_1615341     7W5O_A   26.471             238        150        9
    2055 Query_1615341     4XJ0_A   26.471             238        150        9
    2056 Query_1615341     9HW6_A   26.180             233        134        9
    2057 Query_1615341     6G9J_A   26.471             238        150        9
    2058 Query_1615341     1PME_A   27.197             239        147       10
    2059 Query_1615341     6GZH_A   23.529             204        139        6
    2060 Query_1615341     6G9K_A   26.471             238        150        9
    2061 Query_1615341     2GPH_A   26.471             238        150        9
    2062 Query_1615341     6W4O_O   24.286             210        132        8
    2063 Query_1615341     5LCK_A   26.471             238        150        9
    2064 Query_1615341     9Z8K_B   25.862             232        134        8
    2065 Query_1615341     5K4I_A   26.471             238        150        9
    2066 Query_1615341     8ZJV_A   26.471             238        150        9
    2067 Query_1615341     4QTA_A   26.471             238        150        9
    2068 Query_1615341     4XOY_A   26.471             238        150        9
    2069 Query_1615341     3FME_A   29.114             158        104        5
    2070 Query_1615341     8RMB_A   26.471             238        150        9
    2071 Query_1615341     8PSR_A   26.471             238        150        9
    2072 Query_1615341     2Y9Q_A   26.471             238        150        9
    2073 Query_1615341     4QP4_A   26.471             238        150        9
    2074 Query_1615341     3O71_A   26.471             238        150        9
    2075 Query_1615341     4FV6_A   26.471             238        150        9
    2076 Query_1615341     5AWM_A   24.370             238        138        7
    2077 Query_1615341     8CHF_E   26.022             269        173       10
    2078 Query_1615341     6W9E_A   23.529             204        139        6
    2079 Query_1615341     7XQK_A   25.328             229        136        9
    2080 Query_1615341     3MI9_A   23.529             204        139        6
    2081 Query_1615341     2OJG_A   26.471             238        150        9
    2082 Query_1615341     3C9W_A   26.471             238        150        9
    2083 Query_1615341     5L1Z_A   23.529             204        139        6
    2084 Query_1615341     4OR5_A   23.529             204        139        6
    2085 Query_1615341     6OTS_A   26.471             238        150        9
    2086 Query_1615341     9AXA_B   26.022             269        173       10
    2087 Query_1615341     1WZY_A   26.471             238        150        9
    2088 Query_1615341     1TVO_A   26.471             238        150        9
    2089 Query_1615341     6RFO_A   25.833             240        143       10
    2090 Query_1615341     10JU_B   25.862             232        134        8
    2091 Query_1615341     2FYS_A   26.471             238        150        9
    2092 Query_1615341     6OT6_A   26.471             238        150        9
    2093 Query_1615341     3ZU7_A   26.471             238        150        9
    2094 Query_1615341     4QYY_A   26.471             238        150        9
    2095 Query_1615341     2Z7L_A   26.471             238        150        9
    2096 Query_1615341     6NBS_A   26.471             238        150        9
    2097 Query_1615341     4IMY_A   23.529             204        139        6
    2098 Query_1615341     4ZZM_A   26.471             238        150        9
    2099 Query_1615341     3R63_A   26.471             238        150        9
    2100 Query_1615341     9SKQ_A   24.268             239        137        8
    2101 Query_1615341     4XOZ_A   26.471             238        150        9
    2102 Query_1615341     9O0U_B   26.022             269        173       10
    2103 Query_1615341     8BW9_C   29.048             210        135        9
    2104 Query_1615341     10JU_A   25.862             232        134        8
    2105 Query_1615341     6PP9_B   26.022             269        173       10
    2106 Query_1615341     4FUX_A   26.471             238        150        9
    2107 Query_1615341     7UGB_A   26.471             238        150        9
    2108 Query_1615341     9TYG_A   26.022             269        173       10
    2109 Query_1615341     8ELC_A   27.468             233        137        9
    2110 Query_1615341     5BUE_A   26.471             238        150        9
    2111 Query_1615341     9QQJ_A   26.471             238        150        9
    2112 Query_1615341     3SOA_A   24.286             210        132        8
    2113 Query_1615341     3QYW_A   26.471             238        150        9
    2114 Query_1615341     4Y72_A   24.268             239        137        8
    2115 Query_1615341     6DMG_A   26.471             238        150        9
    2116 Query_1615341     9Z8K_A   25.862             232        134        8
    2117 Query_1615341     7NJ0_B   24.268             239        137        8
    2118 Query_1615341     6NYB_B   26.022             269        173       10
    2119 Query_1615341     5LCJ_A   26.471             238        150        9
    2120 Query_1615341     6Q0T_C   26.022             269        173       10
    2121 Query_1615341     8RLX_A   26.471             238        150        9
    2122 Query_1615341     9TYG_B   26.471             238        150        9
    2123 Query_1615341     4AN2_A   24.922             321        201       12
    2124 Query_1615341     4YC6_A   24.268             239        137        8
    2125 Query_1615341     6FI6_A   26.471             238        150        9
    2126 Query_1615341     4XNE_A   26.471             238        150        9
    2127 Query_1615341     4XP2_A   26.471             238        150        9
    2128 Query_1615341     8AOC_A   26.471             238        150        9
    2129 Query_1615341     6Z45_A   23.039             204        140        6
    2130 Query_1615341     4NST_A   25.490             204        136        6
    2131 Query_1615341     4O6E_A   26.027             219        144        7
    2132 Query_1615341     4FV7_A   26.471             238        150        9
    2133 Query_1615341     4GSB_A   26.471             238        150        9
    2134 Query_1615341     6FXV_A   26.471             238        150        9
    2135 Query_1615341     3SA0_A   26.471             238        150        9
    2136 Query_1615341     2ZOQ_A   23.853             218        154        5
    2137 Query_1615341     5NHH_A   26.471             238        150        9
    2138 Query_1615341     1GOL_A   26.050             238        151        9
    2139 Query_1615341     6GDM_A   26.471             238        150        9
    2140 Query_1615341     2B9H_A   24.199             281        182       11
    2141 Query_1615341     8RM2_A   26.471             238        150        9
    2142 Query_1615341     5NGU_A   26.471             238        150        9
    2143 Query_1615341     8K5R_A   23.039             204        140        6
    2144 Query_1615341     4N0S_A   26.471             238        150        9
    2145 Query_1615341     4JG8_A   24.153             236        140       11
    2146 Query_1615341     4S30_A   26.471             238        150        9
    2147 Query_1615341     3TEI_A   26.471             238        150        9
    2148 Query_1615341     9TU0_A   23.853             218        154        5
    2149 Query_1615341     4QTB_A   23.853             218        154        5
    2150 Query_1615341     4I5H_A   27.197             239        147       10
    2151 Query_1615341     9JK1_A   25.490             204        136        6
    2152 Query_1615341     2Y4I_C   28.230             209        138        8
    2153 Query_1615341     4CXA_A   25.490             204        136        6
    2154 Query_1615341     4UN0_C   25.490             204        136        6
    2155 Query_1615341     2AC3_A   30.070             143         89        5
    2156 Query_1615341     8XFM_A   30.070             143         89        5
    2157 Query_1615341     6B3E_A   25.490             204        136        6
    2158 Query_1615341     7F2X_A   24.062             320        209       10
    2159 Query_1615341     2B9F_A   24.199             281        182       11
    2160 Query_1615341     7E73_A   26.471             238        150        9
    2161 Query_1615341     6XI8_A   25.676             222        140       10
    2162 Query_1615341     2F9G_A   24.199             281        182       11
    2163 Query_1615341     4H3Q_A   26.471             238        150        9
    2164 Query_1615341     7KUE_A   25.676             222        140       10
    2165 Query_1615341     4CRS_A   27.136             199        130        6
    2166 Query_1615341     3ORN_A   28.230             209        138        8
    2167 Query_1615341     4XHL_A   24.528             265        177        8
    2168 Query_1615341     3N9X_A   26.695             236        125       11
    2169 Query_1615341     7E75_A   26.050             238        151        9
    2170 Query_1615341     2P55_A   28.230             209        138        8
    2171 Query_1615341     5WVD_A   30.128             156         97        7
    2172 Query_1615341     7W5C_A   24.638             207        126        7
    2173 Query_1615341     4U7Z_A   28.230             209        138        8
    2174 Query_1615341     1S9J_A   28.230             209        138        8
    2175 Query_1615341     5EYM_A   28.230             209        138        8
    2176 Query_1615341     2HW6_A   30.128             156         97        7
    2177 Query_1615341     3EQC_A   28.230             209        138        8
    2178 Query_1615341     5YXI_A   40.845              71         39        1
    2179 Query_1615341     7JUQ_C   28.230             209        138        8
    2180 Query_1615341     8YP4_A   28.230             209        138        8
    2181 Query_1615341     3ZLY_A   28.230             209        138        8
    2182 Query_1615341     3ZLS_A   28.230             209        138        8
    2183 Query_1615341     4AW2_A   23.016             252        175       10
    2184 Query_1615341     3DV3_A   28.230             209        138        8
    2185 Query_1615341     5HZE_A   28.230             209        138        8
    2186 Query_1615341     3MBL_A   28.230             209        138        8
    2187 Query_1615341     3W8Q_A   28.230             209        138        8
    2188 Query_1615341     7PQV_A   28.230             209        138        8
    2189 Query_1615341     6PXN_A   26.852             216        133        7
    2190 Query_1615341     5YT3_A   27.751             209        139        7
    2191 Query_1615341     6Z3U_B   25.743             202        133        9
    2192 Query_1615341     6BDL_A   27.184             206        127        8
    2193 Query_1615341     1S9I_A   28.230             209        138        8
    2194 Query_1615341     5V5Y_A   25.701             214        122        9
    2195 Query_1615341     7N3U_A   25.701             214        122        9
    2196 Query_1615341     8WDK_W   25.701             214        122        9
    2197 Query_1615341     2IN6_A   25.701             214        122        9
    2198 Query_1615341     7LV3_A   27.184             206        127        8
    2199 Query_1615341     1X8B_A   25.701             214        122        9
    2200 Query_1615341     3BI6_A   25.701             214        122        9
    2201 Query_1615341     7NAA_A   23.394             218        146        8
    2202 Query_1615341     2Z2W_A   25.701             214        122        9
    2203 Query_1615341     4BGQ_A   23.973             292        199       11
    2204 Query_1615341     8H59_A   28.736             174        103        9
    2205 Query_1615341     7T4T_A   27.184             206        127        8
    2206 Query_1615341     6DTL_A   24.510             204        124        6
    2207 Query_1615341     3ENM_A   28.662             157        104        5
    2208 Query_1615341     5Z33_A   28.736             174        103        9
    2209 Query_1615341     4OTD_A   25.180             278        180       11
    2210 Query_1615341     9AXH_A   28.230             209        138        8
    2211 Query_1615341     5CZO_A   23.755             261        178        7
    2212 Query_1615341     4UAK_A   24.506             253        170       10
    2213 Query_1615341     8ZTC_A   28.324             173        105        8
    2214 Query_1615341     3QFV_A   24.506             253        170       10
    2215 Query_1615341     7JV7_A   25.714             210        143        6
    2216 Query_1615341     5CI6_A   24.510             204        124        6
    2217 Query_1615341     3SLS_A   28.230             209        138        8
    2218 Query_1615341     6P5M_A   23.415             205        146        6
    2219 Query_1615341     3TKU_A   24.506             253        170       10
    2220 Query_1615341     7JNT_A   23.415             205        146        6
    2221 Query_1615341     7B3M_A   28.230             209        138        8
    2222 Query_1615341     5U7Q_A   23.415             205        146        6
    2223 Query_1615341     4WOT_A   23.415             205        146        6
    2224 Query_1615341     4AAA_A   22.936             218        160        4
    2225 Query_1615341     5U7R_A   23.415             205        146        6
    2226 Query_1615341     3NIE_A   26.800             250        130       12
    2227 Query_1615341     8X8X_A   23.415             205        146        6
    2228 Query_1615341     9AXM_A   27.751             209        139        7
    2229 Query_1615341     6ED6_A   23.415             205        146        6
    2230 Query_1615341     3ZLW_A   27.751             209        139        8
    2231 Query_1615341     4L6Q_A   23.415             205        146        6
    2232 Query_1615341     6PXP_A   27.315             216        132        7
    2233 Query_1615341     2F2U_A   22.927             205        147        6
    2234 Query_1615341     4KB8_A   27.074             229        136        8
    2235 Query_1615341     4XH0_A   23.596             267        177        8
    2236 Query_1615341     6RUU_A   25.271             277        180       11
    2237 Query_1615341     2H34_A   29.487             156         97        6
    2238 Query_1615341     6ZIW_I   27.586             203        131        8
    2239 Query_1615341     6E9W_A   22.927             205        147        6
    2240 Query_1615341     9JCU_A   22.927             205        147        6
    2241 Query_1615341     3TV7_A   23.448             145        107        3
    2242 Query_1615341     4W7P_A   23.448             145        107        3
    2243 Query_1615341     4QNY_A   29.139             151         89        8
    2244 Query_1615341     7JOU_A   23.448             145        107        3
    2245 Query_1615341     2ESM_A   23.448             145        107        3
    2246 Query_1615341     8ZH5_A   23.448             145        107        3
    2247 Query_1615341     2V55_A   23.448             145        107        3
    2248 Query_1615341     7S26_A   23.448             145        107        3
    2249 Query_1615341     9B3S_A   26.852             216        133        7
    2250 Query_1615341     7S25_A   23.448             145        107        3
    2251 Query_1615341     4TN6_A   26.852             216        133        7
    2252 Query_1615341     8VXF_A   26.852             216        133        7
    2253 Query_1615341     8D7M_A   26.852             216        133        7
    2254 Query_1615341     5MQV_A   26.852             216        133        7
    2255 Query_1615341     3UYS_A   26.852             216        133        7
    2256 Query_1615341     4JJR_A   26.852             216        133        7
    2257 Query_1615341     8VXD_A   26.852             216        133        7
    2258 Query_1615341     5OKT_A   26.852             216        133        7
    2259 Query_1615341     5IH4_A   26.852             216        133        7
    2260 Query_1615341     5X17_A   26.852             216        133        7
    2261 Query_1615341     8IZC_A   26.852             216        133        7
    2262 Query_1615341     7P7F_A   26.852             216        133        7
    2263 Query_1615341     6RCG_A   26.852             216        133        7
    2264 Query_1615341     1CKI_A   26.852             216        133        7
    2265 Query_1615341     4TW9_A   26.852             216        133        7
    2266 Query_1615341     4HNI_A   26.291             213        138        7
    2267 Query_1615341     5CYZ_A   23.372             261        179        7
    2268 Query_1615341    9FQR_Xr   31.915              94         55        3
    2269 Query_1615341     7UP4_A   25.121             207        106        9
    2270 Query_1615341     2ELI_A   37.037              54         30        1
    2271 Query_1615341     5L2Q_A   25.478             157        110        4
    2272 Query_1615341     2VD5_A   26.180             233        145       10
    2273 Query_1615341     3OZ6_A   23.239             284        172       11
    2274 Query_1615341     7WTT_a   29.936             157         93        6
    2275 Query_1615341     6GZD_A   29.936             157         93        6
    2276 Query_1615341     4FI1_A   24.832             149         95        6
    2277 Query_1615341     8TQ2_A   24.537             216        133        9
    2278 Query_1615341     4O2Z_A   27.523             109         71        3
    2279 Query_1615341     9H8C_A   24.537             216        133        9
    2280 Query_1615341     5X18_A   22.400             250        173        8
    2281 Query_1615341     4CRL_A   24.537             216        133        9
    2282 Query_1615341     3RGF_A   24.537             216        133        9
    2283 Query_1615341     6T41_A   24.537             216        133        9
    2284 Query_1615341     6TPA_A   24.537             216        133        9
    2285 Query_1615341     6QTG_A   24.537             216        133        9
    2286 Query_1615341     9R59_A   23.721             215        148        8
    2287 Query_1615341     8XU4_A   23.721             215        148        8
    2288 Query_1615341     6BXI_A   24.569             232        133        9
    2289 Query_1615341     5IDN_A   24.537             216        133        9
    2290 Query_1615341     5HNB_A   24.537             216        133        9
    2291 Query_1615341     5M4U_A   23.973             146         94        6
    2292 Query_1615341     5XQX_A   24.537             216        133        9
    2293 Query_1615341     5FQD_C   30.189             159         90        6
    2294 Query_1615341     9IHG_A   23.973             146         94        6
    2295 Query_1615341     3OFM_A   23.973             146         94        6
    2296 Query_1615341     5OOI_A   23.973             146         94        6
    2297 Query_1615341     5FGK_A   24.537             216        133        9
    2298 Query_1615341     7KPV_A   32.000             100         57        3
    2299 Query_1615341     6HMD_A   23.973             146         94        6
    2300 Query_1615341     9HK5_A   23.973             146         94        6
    2301 Query_1615341     3KA0_A   24.074             216        146        9
    2302 Query_1615341     6L20_A   23.973             146         94        6
    2303 Query_1615341     6QY8_A   23.973             146         94        6
    2304 Query_1615341     9OTY_C   30.189             159         90        6
    2305 Query_1615341     3E3B_X   23.973             146         94        6
    2306 Query_1615341     6T8X_A   22.120             217        152        8
    2307 Query_1615341     3UFF_A   33.333              51         30        1
    2308 Query_1615341     7PSU_A   24.242             132         86        5
    2309 Query_1615341     3UGD_A   33.333              51         30        1
    2310 Query_1615341     2JBO_A   24.074             216        146        9
    2311 Query_1615341     3GOK_A   24.074             216        146        9
    2312 Query_1615341     6TLL_A   24.242             132         86        5
    2313 Query_1615341     2PZY_A   24.074             216        146        9
    2314 Query_1615341     4MD8_E   24.242             132         86        5
    2315 Query_1615341     6TCA_A   24.074             216        146        9
    2316 Query_1615341     3R2B_A   24.074             216        146        9
    2317 Query_1615341     4MD7_E   24.242             132         86        5
    2318 Query_1615341     5OMY_A   24.242             132         86        5
    2319 Query_1615341     4TYH_A   24.074             216        146        9
    2320 Query_1615341     2OZA_A   24.074             216        146        9
    2321 Query_1615341     3FPM_A   24.074             216        146        9
    2322 Query_1615341     4FKD_A   35.294              51         29        1
    2323 Query_1615341     3R2Y_A   24.074             216        146        9
    2324 Query_1615341     2P3G_X   24.074             216        146        9
    2325 Query_1615341     1KWP_A   24.074             216        146        9
    2326 Query_1615341     3UEJ_A   33.333              51         30        1
    2327 Query_1615341     2ONL_C   24.074             216        146        9
    2328 Query_1615341     1PTQ_A   34.000              50         29        1
    2329 Query_1615341     7KND_A   34.000              50         29        1
    2330 Query_1615341     6O6Q_A   21.888             233        150        9
    2331 Query_1615341     2ENZ_A   36.000              50         28        1
    2332 Query_1615341     3UEY_A   33.333              51         30        1
    2333 Query_1615341     6L22_A   23.288             146         95        6
    2334 Query_1615341     3U87_A   24.242             132         86        5
    2335 Query_1615341     6UWA_A   36.538              52         29        1
    2336 Query_1615341     1NA7_A   22.477             218        140       10
    2337 Query_1615341     6L24_A   23.288             146         95        6
    2338 Query_1615341     5XVU_A   25.455             110         75        4
    2339 Query_1615341     2CSN_A   24.299             214        143        5
    2340 Query_1615341     1CSN_A   24.299             214        143        5
    2341 Query_1615341     2ZJW_A   24.242             132         86        5
    2342 Query_1615341     3JUH_A   24.242             132         86        5
    2343 Query_1615341     1NXK_A   26.027             146         95        6
    2344 Query_1615341     6L23_A   24.242             132         86        5
    2345 Query_1615341     5CSP_A   22.477             218        140       10
    2346 Query_1615341     5OSL_A   22.477             218        140       10
    2347 Query_1615341     6QY7_A   24.242             132         86        5
    2348 Query_1615341     5MOV_A   22.477             218        140       10
    2349 Query_1615341     9HXU_A   24.242             132         86        5
    2350 Query_1615341     4MD9_E   23.358             137         91        5
    2351 Query_1615341     1JWH_A   24.242             132         86        5
    2352 Query_1615341     3Q9W_A   24.242             132         86        5
    2353 Query_1615341     6HME_A   24.242             132         86        5
    2354 Query_1615341     5N1V_A   24.242             132         86        5
    2355 Query_1615341     7A4Q_A   22.477             218        140       10
    2356 Query_1615341     1XA6_A   36.364              66         34        2
    2357 Query_1615341     3W8L_A   24.242             132         86        5
    2358 Query_1615341     5KU8_A   24.242             132         86        5
    2359 Query_1615341     1PJK_A   24.242             132         86        5
    2360 Query_1615341     2R7I_A   24.242             132         86        5
    2361 Query_1615341     3BQC_A   24.242             132         86        5
    2362 Query_1615341   6Z83_AAA   24.242             132         86        5
    2363 Query_1615341     3NSZ_A   24.242             132         86        5
    2364 Query_1615341     3H30_A   24.242             132         86        5
    2365 Query_1615341     3NGA_A   24.242             132         86        5
    2366 Query_1615341     6Z19_B   24.242             132         86        5
    2367 Query_1615341     3MB6_A   24.242             132         86        5
    2368 Query_1615341     5ZN5_A   24.242             132         86        5
    2369 Query_1615341     3Q04_A   24.242             132         86        5
    2370 Query_1615341     6Q38_A   24.242             132         86        5
    2371 Query_1615341     6YPH_A   24.242             132         86        5
    2372 Query_1615341     1TBN_A   35.185              54         31        1
    2373 Query_1615341     5CVG_A   24.242             132         86        5
    2374 Query_1615341     5CLP_A   24.242             132         86        5
    2375 Query_1615341     7QUX_A   24.242             132         86        5
    2376 Query_1615341     6L21_A   24.242             132         86        5
    2377 Query_1615341     6YZH_A   24.242             132         86        5
    2378 Query_1615341     5ZN0_A   23.358             137         91        5
    2379 Query_1615341     3CXL_A   26.000             150         89        5
    2380 Query_1615341     9EDY_A   21.834             229        130        9
    2381 Query_1615341     6U69_A   21.586             227        133        9
    2382 Query_1615341     4DGL_C   24.545             110         76        4
    2383 Query_1615341     6JKK_A   19.910             221        152        5
    2384 Query_1615341     1KBE_A   34.091              44         28        1
    2385 Query_1615341     5M07_A   22.749             211        128        7
    2386 Query_1615341     5M06_A   22.749             211        128        7
    2387 Query_1615341     2PVH_A   22.603             146        103        5
    2388 Query_1615341     4JRN_A   22.086             163        110        3
    2389 Query_1615341     5XKA_A   22.596             208        132        6
    2390 Query_1615341     4ANM_A   22.603             146        103        5
    2391 Query_1615341     1DS5_A   23.636             110         77        4
    2392 Query_1615341     3PVG_A   22.603             146        103        5
    2393 Query_1615341     1M2P_A   23.636             110         77        4
    2394 Query_1615341     1DAW_A   23.636             110         77        4
    2395 Query_1615341     2QC6_A   22.603             146        103        5
    2396 Query_1615341     4DGN_A   22.603             146        103        5
    2397 Query_1615341     3KXG_A   22.603             146        103        5
    2398 Query_1615341     5TS8_A   22.603             146        103        5
    2399 Query_1615341     4DGM_A   22.603             146        103        5
    2400 Query_1615341     1Y8F_A   42.857              42         20        1
    2401 Query_1615341     6RA0_A   39.474              38         19        1
    2402 Query_1615341     3KGA_A   21.860             215        131        8
    2403 Query_1615341     7P5Z_1   26.115             157         88        6
    2404 Query_1615341     2YUU_A   22.388              67         43        2
    2405 Query_1615341     7T7C_A   37.736              53         29        1
    2406 Query_1615341     5UE8_A   37.736              53         29        1
    2407 Query_1615341     6VP8_B   28.358              67         43        2
    2408 Query_1615341     7DG2_A   27.119              59         42        1
    2409 Query_1615341     3UIB_A   23.200             250        132       12
    2410 Query_1615341     3PG1_A   23.200             250        132       12
    2411 Query_1615341     2ENN_A   28.000              50         32        1
    2412 Query_1615341     2E73_A   32.000              50         30        1
    2413 Query_1615341     2DB6_A   29.167              48         30        1
    2414 Query_1615341     9C5F_A   45.833              24         13        0
    2415 Query_1615341     5HNV_A   19.915             236        138        6
    2416 Query_1615341     4IX3_A   23.214             168         98        7
    2417 Query_1615341     4B6D_A   32.000              50         31        1
         q.start q.end s.start s.end    evalue bitscore positives
    1          1   766      28   793  0.00e+00   1591.0     99.87
    2          1   766       1   766  0.00e+00   1590.0     99.74
    3          1   766      28   793  0.00e+00   1590.0     99.87
    4          1   766      28   793  0.00e+00   1589.0     99.74
    5          1   766       2   767  0.00e+00   1589.0     99.74
    6          1   766      63   828  0.00e+00   1588.0     99.87
    7          1   766       3   768  0.00e+00   1587.0     99.74
    8          1   766       2   767  0.00e+00   1587.0     99.74
    9          1   766       2   767  0.00e+00   1586.0     99.61
    10         1   766       2   767  0.00e+00   1586.0     99.61
    11         1   766       2   767  0.00e+00   1583.0     99.48
    12       360   734       1   375  0.00e+00    774.0     99.20
    13       360   734       1   375  0.00e+00    774.0     99.20
    14        97   755      27   669  0.00e+00    711.0     68.65
    15        97   755      27   669  0.00e+00    709.0     68.50
    16        97   755      27   669  0.00e+00    705.0     68.35
    17       157   755      58   647  0.00e+00    703.0     72.27
    18        97   755       6   648  0.00e+00    701.0     68.35
    19        97   755      27   669  0.00e+00    701.0     68.05
    20       157   760      58   652  0.00e+00    700.0     71.84
    21       432   736      25   329  0.00e+00    637.0     99.67
    22       432   726      13   307  0.00e+00    621.0    100.00
    23       432   726      12   306  0.00e+00    621.0    100.00
    24       432   726      13   307  0.00e+00    621.0    100.00
    25       432   726      13   307  0.00e+00    621.0    100.00
    26       433   726       7   300  0.00e+00    620.0    100.00
    27       432   726      13   307  0.00e+00    618.0     99.66
    28       433   727       6   300  0.00e+00    616.0     99.32
    29       432   723      13   304  0.00e+00    614.0    100.00
    30       432   726      13   302  0.00e+00    604.0     98.31
    31       445   726       3   284  0.00e+00    595.0    100.00
    32       446   727       1   282  0.00e+00    594.0    100.00
    33       442   723       2   283  0.00e+00    590.0     99.65
    34       445   723       6   284  0.00e+00    588.0    100.00
    35       445   723       3   281  0.00e+00    588.0    100.00
    36       445   723       2   280  0.00e+00    587.0    100.00
    37       448   722       1   275  0.00e+00    580.0    100.00
    38       448   723       1   276  0.00e+00    580.0     99.64
    39       448   723       1   276  0.00e+00    579.0     99.64
    40       445   723      16   294  0.00e+00    578.0     98.92
    41       449   721       1   273  0.00e+00    576.0    100.00
    42       447   735       2   290  0.00e+00    572.0     95.50
    43       449   720       1   272  0.00e+00    563.0     98.16
    44       442   721      11   290  0.00e+00    554.0     95.36
    45       444   721       5   282  0.00e+00    552.0     95.68
    46       442   721       1   280  0.00e+00    552.0     95.36
    47       444   721       1   278  0.00e+00    552.0     95.68
    48       443   721       1   279  0.00e+00    552.0     95.34
    49       444   721      21   298  0.00e+00    551.0     95.68
    50       442   721      11   290  0.00e+00    551.0     95.00
    51       442   721       2   281  0.00e+00    551.0     95.36
    52       442   721       7   286  0.00e+00    550.0     95.00
    53       444   721       3   280  0.00e+00    550.0     95.32
    54       442   721      11   290  0.00e+00    550.0     95.00
    55       442   721       2   281  0.00e+00    548.0     95.00
    56       442   719      11   288  0.00e+00    546.0     94.96
    57       443   721      17   295  0.00e+00    544.0     94.27
    58       448   721       5   278  0.00e+00    542.0     95.26
    59       443   721      17   295  0.00e+00    541.0     93.91
    60       448   721       3   276  0.00e+00    541.0     95.26
    61       448   719       2   273  0.00e+00    540.0     95.59
    62       442   721       7   281  0.00e+00    536.0     93.57
    63       418   755     224   562  0.00e+00    532.0     86.26
    64       418   755     224   562  0.00e+00    531.0     85.96
    65       432   726      13   307 1.70e-171    496.0     89.15
    66       445   723       2   280 4.39e-169    489.0     89.96
    67       445   723       2   280 6.14e-169    489.0     89.96
    68       442   726       1   285 9.84e-169    489.0     89.12
    69       445   723       2   280 6.04e-163    473.0     90.32
    70       512   717       1   206 1.40e-124    372.0     92.72
    71       157   283       3   132  1.19e-54    186.0     76.92
    72       157   283       4   133  1.23e-54    186.0     76.92
    73       157   283       7   136  1.40e-54    186.0     76.92
    74       449   731      27   310  2.81e-54    191.0     59.04
    75       449   731      27   310  2.95e-54    191.0     59.04
    76       446   731      18   305  3.60e-54    191.0     59.66
    77       149   232      12    95  4.11e-54    182.0    100.00
    78       449   731      50   333  5.21e-54    191.0     59.04
    79       153   237      12    96  5.75e-54    182.0     97.65
    80       157   283       7   136  4.92e-53    181.0     76.15
    81       137   232       6    92  1.58e-52    178.0     89.58
    82       449   731      40   325  3.19e-48    175.0     55.93
    83        21   114       1    88  5.10e-48    166.0     89.36
    84       449   719       6   283  2.53e-46    168.0     56.18
    85       448   727      30   307  1.11e-45    167.0     59.93
    86       448   717       4   267  1.57e-45    166.0     55.20
    87       448   727      30   307  3.58e-45    166.0     59.93
    88        38   110     160   232  2.82e-42    155.0    100.00
    89       446   716      25   277  3.74e-42    157.0     54.95
    90       446   716       3   255  4.05e-42    157.0     54.95
    91       446   715       9   270  4.63e-40    150.0     54.51
    92       446   715       6   267  1.25e-39    149.0     54.15
    93       446   715       2   263  1.30e-39    149.0     54.15
    94       446   715      10   271  1.53e-39    149.0     54.15
    95       446   715      16   277  1.78e-39    149.0     54.15
    96       446   715       9   270  1.85e-39    149.0     54.15
    97       448   715       1   260  2.19e-39    148.0     54.18
    98       448   715       2   261  2.30e-39    148.0     54.18
    99       448   707       6   258  3.07e-39    148.0     54.48
    100      446   715       9   270  3.19e-39    148.0     54.15
    101      446   715       6   267  3.77e-39    148.0     54.15
    102      446   715       9   270  4.40e-39    148.0     53.79
    103      446   715       9   270  4.67e-39    147.0     53.79
    104      446   715       9   270  5.05e-39    147.0     53.79
    105      448   715       1   260  5.07e-39    147.0     53.82
    106      448   715       1   260  5.19e-39    147.0     57.04
    107      446   715       6   267  5.25e-39    147.0     53.99
    108      447   715       1   261  5.36e-39    147.0     53.62
    109      446   715       9   270  5.95e-39    147.0     53.79
    110      447   715      15   275  7.29e-39    147.0     56.83
    111      447   715      10   270  7.39e-39    147.0     56.83
    112      447   715      11   271  7.61e-39    147.0     56.83
    113      450   707       7   257  7.64e-39    147.0     54.51
    114      447   715       5   265  9.43e-39    146.0     56.83
    115      447   715       1   261  1.11e-38    146.0     53.99
    116      446   715       9   270  1.34e-38    146.0     53.43
    117      447   715      13   273  1.76e-38    146.0     56.73
    118      446   715       9   270  1.80e-38    146.0     53.43
    119      446   715       9   270  1.80e-38    146.0     53.43
    120      446   715       9   270  1.91e-38    146.0     53.43
    121      447   715      11   271  2.08e-38    146.0     56.73
    122      447   715       5   265  2.41e-38    145.0     56.73
    123      446   715       9   270  2.46e-38    145.0     53.43
    124      447   715       7   267  2.51e-38    145.0     56.73
    125      447   715      14   274  2.52e-38    145.0     56.73
    126      447   715      14   274  2.59e-38    145.0     56.73
    127      450   715       2   259  2.66e-38    145.0     53.85
    128      446   715       1   262  2.72e-38    145.0     53.43
    129      446   715       9   270  2.91e-38    145.0     53.43
    130      447   715       5   265  2.95e-38    145.0     56.73
    131      447   715       6   266  3.02e-38    145.0     56.73
    132      447   715       1   261  4.02e-38    145.0     53.26
    133      446   715     258   519  4.35e-38    151.0     54.15
    134      446   715     259   520  4.66e-38    151.0     54.15
    135      446   715     174   435  5.44e-38    149.0     54.15
    136      450   707       7   257  5.45e-38    144.0     54.51
    137      446   715     175   436  5.75e-38    149.0     54.15
    138      448   710       8   262  5.90e-38    144.0     53.70
    139      448   710      19   273  6.48e-38    144.0     53.70
    140      446   715     175   436  6.56e-38    149.0     54.15
    141      446   715       9   270  1.12e-37    144.0     53.07
    142      446   715     200   461  1.12e-37    149.0     54.15
    143      446   715       9   270  1.40e-37    143.0     53.45
    144      449   715       3   261  1.40e-37    143.0     56.88
    145      446   715     174   435  1.47e-37    148.0     54.15
    146      446   715     178   439  1.49e-37    148.0     54.15
    147      446   715     173   434  3.83e-37    146.0     53.79
    148      446   715     175   436  4.26e-37    146.0     53.79
    149      447   715       5   265  6.11e-37    141.0     56.00
    150      446   715     176   437  8.24e-37    145.0     53.43
    151      451   732       4   288  9.36e-37    142.0     54.24
    152      451   732      41   325  9.47e-37    144.0     54.24
    153      451   732      12   296  1.09e-36    142.0     53.72
    154      451   732      29   313  1.78e-36    142.0     54.24
    155      448   715     175   424  1.85e-36    144.0     53.68
    156      451   732      41   325  1.85e-36    143.0     54.24
    157      451   732      41   325  2.61e-36    142.0     54.24
    158      451   732      41   325  2.77e-36    142.0     54.24
    159      451   732      41   325  2.77e-36    142.0     54.24
    160      451   732      39   323  2.84e-36    142.0     54.24
    161      451   732      41   325  3.02e-36    142.0     54.24
    162      448   715     181   440  3.06e-36    144.0     53.82
    163      451   732      41   325  3.13e-36    142.0     53.90
    164      448   715     181   440  3.33e-36    144.0     53.82
    165      451   730      15   296  4.41e-36    139.0     54.64
    166      451   732      12   296  4.43e-36    141.0     53.38
    167      451   730      14   295  6.91e-36    139.0     54.64
    168      451   732      41   325  8.88e-36    140.0     53.90
    169      451   732      41   325  1.16e-35    140.0     53.90
    170      447   710       2   258  1.74e-35    137.0     54.10
    171      451   732       4   285  1.99e-35    137.0     54.11
    172      451   732      10   291  2.08e-35    137.0     54.11
    173      447   710       2   258  2.11e-35    136.0     54.10
    174      451   732      25   306  2.18e-35    138.0     54.11
    175      442   729      36   323  2.30e-35    138.0     51.52
    176      451   720      12   283  3.44e-35    137.0     55.16
    177      157   237      14    94  4.52e-35    129.0     83.95
    178      442   729      36   323  4.90e-35    137.0     51.34
    179      451   730      15   296  5.21e-35    136.0     53.95
    180      451   730      10   291  5.39e-35    136.0     53.26
    181      451   720      22   293  5.79e-35    136.0     54.80
    182      451   720      22   293  8.27e-35    136.0     54.80
    183      439   715     176   442  9.64e-35    139.0     52.82
    184      439   715     176   442  1.01e-34    139.0     52.82
    185      451   720       3   274  1.34e-34    135.0     54.80
    186      462   715       9   254  1.38e-34    134.0     53.64
    187      451   720      15   286  1.42e-34    135.0     54.80
    188      451   720      18   289  1.46e-34    135.0     55.67
    189      451   720      18   287  2.66e-34    134.0     53.76
    190      451   732      15   298  4.58e-34    134.0     54.76
    191      451   720      15   286  6.78e-34    133.0     53.38
    192      451   720      29   300  1.79e-33    132.0     53.38
    193      451   710       2   254  1.98e-33    131.0     53.79
    194      451   710       2   254  2.08e-33    131.0     53.79
    195      448   715       3   280  2.33e-33    132.0     50.87
    196      451   710       2   254  2.55e-33    130.0     53.79
    197      448   715      12   289  2.83e-33    132.0     50.87
    198      451   720      13   283  3.68e-33    131.0     53.57
    199      451   720      12   282  3.92e-33    131.0     53.57
    200      448   715      12   289  5.33e-33    131.0     51.19
    201      448   715       9   286  5.86e-33    130.0     51.19
    202      448   715      12   289  6.57e-33    130.0     51.19
    203      448   715       5   282  6.60e-33    130.0     51.19
    204      451   720      40   310  7.23e-33    131.0     53.57
    205      448   715      18   295  8.14e-33    131.0     51.19
    206      448   715      11   288  8.48e-33    130.0     51.19
    207      451   719      30   299  8.65e-33    130.0     52.67
    208      448   715       6   283  9.00e-33    130.0     51.19
    209      448   715      12   289  9.19e-33    130.0     51.19
    210      448   715      11   288  9.48e-33    130.0     51.19
    211      448   715      18   295  9.64e-33    130.0     51.19
    212      448   715      12   289  9.73e-33    130.0     51.19
    213      448   715       5   282  1.09e-32    130.0     51.19
    214      451   719      30   299  1.22e-32    130.0     53.41
    215      451   719      39   308  1.40e-32    130.0     53.41
    216      440   726      97   375  1.42e-32    131.0     50.17
    217      451   720      21   291  2.71e-32    131.0     53.57
    218      440   726      97   375  2.86e-32    130.0     49.83
    219      451   720      21   291  2.89e-32    131.0     53.57
    220      451   720      21   291  3.00e-32    131.0     53.57
    221      448   715      40   317  3.01e-32    129.0     51.19
    222      440   726      99   377  3.06e-32    130.0     49.83
    223      457   720      35   298  6.32e-32    128.0     53.65
    224      461   722      21   271  6.69e-32    128.0     51.84
    225      448   715       8   285  6.82e-32    127.0     50.51
    226      461   722      18   268  7.15e-32    127.0     51.84
    227      448   715      18   295  7.44e-32    128.0     50.51
    228      451   720      29   298  7.67e-32    127.0     53.21
    229      457   720      45   308  7.73e-32    128.0     53.65
    230      451   720      29   298  8.02e-32    127.0     53.21
    231      448   730      32   325  8.45e-32    128.0     50.82
    232      448   730      10   303  1.04e-31    127.0     50.82
    233      440   726      99   377  1.06e-31    129.0     48.81
    234      460   710      27   288  1.40e-31    127.0     51.31
    235      456   714      27   287  1.45e-31    127.0     51.58
    236      460   710      27   288  1.58e-31    126.0     51.31
    237      448   711       7   280  1.91e-31    126.0     50.88
    238      460   710      27   288  1.99e-31    126.0     51.31
    239      450   714       9   265  2.38e-31    125.0     52.38
    240      448   711      10   283  2.42e-31    126.0     50.88
    241      448   715     970  1247  2.76e-31    133.0     51.03
    242      447   714       2   261  2.83e-31    125.0     51.81
    243      450   714      12   268  3.24e-31    125.0     52.01
    244      446   717       2   267  3.31e-31    125.0     50.53
    245      448   730      10   303  3.32e-31    125.0     50.00
    246      421   714       1   285  3.43e-31    125.0     50.83
    247      448   730      10   303  3.48e-31    125.0     51.13
    248      448   717       1   264  3.61e-31    124.0     50.54
    249      448   730      10   303  3.65e-31    125.0     50.49
    250      441   717       3   273  3.72e-31    125.0     50.35
    251      450   714      19   275  3.75e-31    125.0     52.01
    252      450   714      16   272  3.77e-31    125.0     52.01
    253      435   717     105   381  3.94e-31    128.0     50.00
    254      450   714      13   269  4.10e-31    124.0     52.01
    255      450   714      14   270  4.16e-31    124.0     52.01
    256      450   714      19   275  4.54e-31    125.0     52.01
    257      450   714       2   258  4.54e-31    124.0     52.01
    258      450   714      11   267  4.63e-31    124.0     52.01
    259      450   714       3   259  4.65e-31    124.0     52.01
    260      450   714      12   268  4.76e-31    124.0     52.01
    261      450   714      15   271  4.86e-31    124.0     52.01
    262      431   714       2   275  5.09e-31    124.0     50.68
    263      450   714       4   260  5.42e-31    124.0     52.01
    264      450   714       3   259  5.66e-31    124.0     52.01
    265      450   714      10   266  5.71e-31    124.0     52.01
    266      442   717       1   270  5.94e-31    124.0     50.18
    267      448   730      10   303  5.99e-31    125.0     50.81
    268      450   714       1   257  6.17e-31    124.0     52.01
    269      450   714       5   261  6.32e-31    124.0     52.01
    270      450   714       4   260  6.32e-31    124.0     52.01
    271      450   714       4   260  6.42e-31    124.0     52.01
    272      450   714       6   262  6.60e-31    124.0     52.01
    273      450   714       7   263  6.93e-31    124.0     52.38
    274      450   714       7   263  7.13e-31    124.0     52.01
    275      446   717      23   288  7.21e-31    124.0     49.47
    276      448   717      11   274  7.24e-31    124.0     50.54
    277      450   714       7   263  7.37e-31    124.0     52.01
    278      448   717       1   264  7.41e-31    124.0     49.46
    279      450   714       7   263  7.48e-31    124.0     52.01
    280      450   714       9   265  7.80e-31    124.0     52.01
    281      450   714      14   270  8.20e-31    124.0     52.01
    282      448   717       5   268  8.26e-31    124.0     50.54
    283      442   717       1   270  8.33e-31    124.0     50.18
    284      450   714       4   260  8.51e-31    124.0     52.01
    285      448   717      11   274  8.65e-31    124.0     50.54
    286      435   717     158   434  8.88e-31    127.0     50.00
    287      448   717      10   273  8.90e-31    124.0     50.54
    288      460   710      36   297  9.05e-31    124.0     50.94
    289      450   714       4   260  9.27e-31    123.0     52.01
    290      448   717      11   274  9.28e-31    124.0     50.54
    291      450   714       4   260  9.45e-31    123.0     52.01
    292      448   717       6   269  9.98e-31    124.0     50.54
    293      460   710      31   292  1.00e-30    124.0     50.94
    294      448   717       6   269  1.02e-30    124.0     50.54
    295      460   710      30   291  1.04e-30    124.0     50.94
    296      460   710      41   302  1.04e-30    124.0     50.94
    297      450   714       4   260  1.04e-30    124.0     52.01
    298      448   717       5   268  1.05e-30    124.0     50.54
    299      448   717       5   268  1.05e-30    124.0     50.54
    300      450   714       2   258  1.06e-30    123.0     52.75
    301      460   710      36   297  1.06e-30    124.0     50.94
    302      460   710      37   298  1.10e-30    124.0     50.94
    303      448   717       8   271  1.10e-30    124.0     50.54
    304      448   717       5   268  1.10e-30    123.0     50.54
    305      460   710      37   298  1.14e-30    124.0     50.94
    306      450   714      26   282  1.17e-30    124.0     52.01
    307      448   717       6   269  1.17e-30    123.0     50.54
    308      460   710      33   294  1.21e-30    124.0     50.94
    309      460   710      34   295  1.21e-30    124.0     50.94
    310      448   717       6   269  1.25e-30    123.0     50.54
    311      448   717       8   271  1.28e-30    123.0     50.54
    312      450   714       4   260  1.32e-30    123.0     52.01
    313      448   730      10   303  1.32e-30    124.0     51.13
    314      448   717       8   271  1.33e-30    123.0     50.54
    315      448   730      11   304  1.34e-30    124.0     51.47
    316      460   710      34   295  1.36e-30    124.0     50.94
    317      460   705      27   283  1.40e-30    124.0     51.15
    318      451   715       2   258  1.41e-30    123.0     51.82
    319      441   717       3   273  1.44e-30    123.0     50.00
    320      446   717       2   267  1.47e-30    123.0     50.18
    321      460   705      26   282  1.49e-30    124.0     51.15
    322      461   722      20   270  1.62e-30    130.0     51.84
    323      451   715       4   260  1.72e-30    122.0     51.82
    324      460   710      51   312  1.79e-30    124.0     50.94
    325      460   710      32   293  1.80e-30    124.0     50.56
    326      460   710      36   297  1.80e-30    124.0     50.56
    327      460   710      30   291  1.84e-30    123.0     50.56
    328      460   710      52   313  1.89e-30    124.0     50.94
    329      448   717      19   282  1.91e-30    123.0     50.54
    330      460   710      31   292  1.91e-30    123.0     50.56
    331      448   717      33   296  1.95e-30    124.0     50.54
    332      460   710      29   290  2.00e-30    123.0     50.56
    333      460   710      37   298  2.02e-30    124.0     50.56
    334      460   710      33   294  2.06e-30    123.0     50.56
    335      460   710      34   295  2.06e-30    123.0     50.56
    336      460   710      28   289  2.21e-30    123.0     50.56
    337      460   705      26   282  2.26e-30    123.0     50.76
    338      451   715      14   270  2.33e-30    122.0     51.82
    339      460   705      26   282  2.47e-30    123.0     50.76
    340      448   715     990  1267  2.73e-30    130.0     51.19
    341      448   730      12   305  2.74e-30    123.0     50.81
    342      435   717     200   476  2.75e-30    127.0     50.00
    343      448   730      11   304  2.79e-30    123.0     50.81
    344      460   705      26   282  2.83e-30    122.0     50.76
    345      450   714     183   439  2.92e-30    126.0     52.75
    346      450   714      19   275  3.07e-30    122.0     51.65
    347      460   710      52   313  3.25e-30    123.0     50.56
    348      460   710      52   313  3.32e-30    123.0     50.56
    349      448   717      11   274  3.38e-30    122.0     50.18
    350      460   710      46   307  3.43e-30    123.0     50.56
    351      460   705      56   312  3.43e-30    123.0     50.76
    352      448   730    1004  1297  3.74e-30    130.0     50.65
    353      448   717       5   268  3.87e-30    122.0     50.18
    354      435   717     239   515  4.08e-30    127.0     50.00
    355      448   717       5   268  4.35e-30    121.0     50.18
    356      448   715     965  1243  4.41e-30    129.0     51.03
    357      448   730     977  1270  4.47e-30    129.0     50.65
    358      448   717      11   274  4.70e-30    121.0     50.18
    359      435   717     197   473  4.72e-30    126.0     50.00
    360      456   715      43   301  4.84e-30    122.0     52.04
    361      460   710      41   302  4.89e-30    122.0     50.94
    362      456   715      62   320  4.90e-30    123.0     52.42
    363      450   714       2   258  5.11e-30    121.0     51.65
    364      456   715      62   320  5.18e-30    123.0     52.04
    365      460   710      31   292  5.21e-30    122.0     50.94
    366      450   714      15   271  5.37e-30    121.0     51.65
    367      460   710      26   287  5.38e-30    122.0     50.94
    368      448   730     977  1270  5.41e-30    129.0     50.65
    369      450   714     182   438  5.48e-30    125.0     52.01
    370      448   717       6   269  5.54e-30    122.0     50.18
    371      448   717       6   269  6.25e-30    121.0     50.18
    372      460   705     189   445  6.74e-30    125.0     51.53
    373      448   717       1   264  7.06e-30    121.0     50.18
    374      435   717     200   476  1.10e-29    125.0     49.66
    375      448   715     961  1239  1.39e-29    128.0     51.03
    376      451   714       4   259  1.46e-29    120.0     51.84
    377      446   717       2   267  1.71e-29    120.0     49.82
    378      451   714       2   257  1.85e-29    119.0     51.84
    379      448   717       7   270  2.18e-29    120.0     49.82
    380      450   715      18   293  2.56e-29    120.0     48.25
    381      450   715      41   326  3.06e-29    120.0     46.96
    382      450   715      41   326  3.21e-29    120.0     46.96
    383      451   714       2   257  3.23e-29    119.0     52.57
    384      457   714       6   255  4.20e-29    118.0     52.26
    385      450   714     395   651  4.57e-29    125.0     53.11
    386      450   715      18   292  6.20e-29    119.0     48.07
    387      448   717       5   268  8.59e-29    118.0     48.91
    388      448   730     987  1280  1.01e-28    125.0     51.14
    389      448   730    1014  1307  1.02e-28    125.0     51.14
    390      448   730    1002  1295  1.07e-28    125.0     51.14
    391      448   730    1002  1295  1.09e-28    125.0     51.14
    392      450   715      18   307  2.03e-28    118.0     46.33
    393      451   711      12   282  2.04e-28    118.0     50.70
    394      451   714     194   449  2.51e-28    120.0     52.94
    395      450   715      18   307  2.55e-28    117.0     46.33
    396      450   715      41   326  2.74e-28    118.0     46.62
    397      448   730     975  1268  3.05e-28    124.0     51.14
    398      450   715      18   307  5.14e-28    117.0     46.00
    399      446   714       4   269  6.85e-28    115.0     49.29
    400      446   714       2   267  7.00e-28    115.0     49.29
    401      450   715      18   307  7.16e-28    116.0     46.00
    402      446   714       4   269  7.67e-28    115.0     49.29
    403      446   714       3   268  7.72e-28    115.0     49.29
    404      446   714       3   268  7.80e-28    115.0     49.29
    405      450   715      18   307  9.60e-28    116.0     46.00
    406      460   710      27   288  9.83e-28    115.0     49.06
    407      450   715      18   307  9.97e-28    116.0     46.00
    408      446   714       6   271  1.03e-27    115.0     49.29
    409      463   719      46   301  1.06e-27    116.0     49.43
    410      446   714       1   266  1.08e-27    115.0     49.29
    411      446   714       1   266  1.11e-27    115.0     49.29
    412      463   720      23   279  1.46e-27    115.0     49.25
    413      446   714      29   294  1.47e-27    115.0     49.29
    414      450   715      36   319  1.86e-27    115.0     47.46
    415      450   715      35   318  1.86e-27    115.0     47.46
    416      452   707      15   272  1.87e-27    114.0     47.58
    417      451   710       5   280  1.91e-27    114.0     48.44
    418      452   707       5   262  2.10e-27    114.0     47.58
    419      450   715      36   319  2.12e-27    115.0     47.46
    420      452   707      15   272  2.25e-27    114.0     47.58
    421      452   707      15   272  2.43e-27    114.0     47.58
    422      450   720      12   294  2.60e-27    114.0     47.65
    423      452   707      17   274  2.97e-27    114.0     47.58
    424      449   714       1   263  3.03e-27    113.0     49.10
    425      451   710      22   297  3.06e-27    114.0     48.44
    426      452   707       5   262  3.08e-27    113.0     47.58
    427      446   714       1   266  3.12e-27    113.0     48.93
    428      450   715      35   318  3.16e-27    115.0     46.78
    429      452   707      17   274  3.18e-27    114.0     47.58
    430      452   707       5   262  3.24e-27    113.0     47.58
    431      446   714       1   266  3.63e-27    113.0     48.93
    432      452   707       5   262  3.63e-27    113.0     47.58
    433      452   707       9   266  3.64e-27    113.0     47.58
    434      450   715      59   348  3.96e-27    115.0     46.36
    435      450   715      18   307  4.01e-27    114.0     45.67
    436      444   715      12   294  5.99e-27    114.0     48.64
    437      457   714      41   300  7.89e-27    113.0     49.45
    438      457   714      28   287  7.90e-27    112.0     49.45
    439      448   715    1362  1651  8.19e-27    119.0     46.84
    440      452   707       9   266  8.44e-27    114.0     47.58
    441      450   715      38   327  8.64e-27    114.0     46.36
    442      457   714      28   287  9.33e-27    112.0     49.45
    443      450   715      15   304  1.00e-26    113.0     45.67
    444      378   714     309   646  1.02e-26    117.0     46.18
    445      450   720      18   309  1.05e-26    113.0     46.58
    446      450   720      16   307  1.15e-26    112.0     46.58
    447      450   720       7   298  1.20e-26    112.0     46.58
    448      450   720      16   307  1.34e-26    112.0     46.58
    449      450   720      18   309  1.36e-26    112.0     46.58
    450      450   720       8   299  1.40e-26    112.0     46.25
    451      450   720      18   309  1.59e-26    112.0     46.25
    452      450   720      14   305  1.60e-26    112.0     46.25
    453      450   720      30   321  1.61e-26    112.0     46.58
    454      450   720      14   305  1.65e-26    112.0     46.25
    455      450   720      18   309  1.65e-26    112.0     46.58
    456      450   715       7   260  1.66e-26    111.0     48.72
    457      450   720       8   299  1.67e-26    112.0     46.25
    458      450   720      18   309  1.70e-26    112.0     46.25
    459      450   715      18   307  1.74e-26    112.0     45.67
    460      450   720       3   265  1.78e-26    111.0     50.00
    461      450   715      16   269  1.81e-26    111.0     48.72
    462      450   720      16   307  1.83e-26    112.0     46.58
    463      450   720      18   309  1.83e-26    112.0     46.58
    464      450   720      16   307  1.88e-26    112.0     46.58
    465      461   736      14   275  1.92e-26    112.0     47.37
    466      450   715       1   254  1.96e-26    110.0     48.72
    467      457   714      19   278  2.09e-26    111.0     49.08
    468      457   714      28   287  2.23e-26    111.0     49.08
    469      457   714      24   283  2.25e-26    111.0     49.08
    470      450   714       1   262  2.26e-26    110.0     48.91
    471      450   720      30   321  2.27e-26    112.0     46.25
    472      457   714      41   300  2.28e-26    112.0     49.08
    473      457   714      28   287  2.31e-26    111.0     49.08
    474      457   714      41   300  2.36e-26    112.0     49.08
    475      378   714     311   648  2.49e-26    116.0     45.89
    476      450   720      30   321  2.49e-26    112.0     46.58
    477      457   714      24   283  2.52e-26    111.0     49.08
    478      378   714     309   646  2.60e-26    116.0     45.89
    479      450   715      23   307  2.74e-26    112.0     48.30
    480      457   714      31   290  2.86e-26    111.0     49.08
    481      450   720      29   320  2.87e-26    112.0     46.25
    482      457   714      27   286  2.96e-26    111.0     49.08
    483      457   714      50   309  2.98e-26    112.0     49.08
    484      457   714      49   308  3.09e-26    111.0     49.08
    485      457   714      32   291  3.26e-26    111.0     49.08
    486      461   736      15   276  3.37e-26    110.0     47.37
    487      457   714      30   289  3.43e-26    111.0     49.08
    488      461   736      45   306  3.59e-26    112.0     48.40
    489      461   736      22   283  3.68e-26    110.0     48.40
    490      457   714      32   291  3.75e-26    111.0     49.08
    491      457   714      29   288  3.79e-26    111.0     49.08
    492      461   736      21   282  3.79e-26    111.0     48.40
    493      450   720      30   321  3.81e-26    111.0     46.25
    494      453   736       7   276  3.84e-26    110.0     46.76
    495      461   736      23   284  3.87e-26    111.0     48.40
    496      457   714      54   313  3.87e-26    111.0     49.08
    497      461   736      22   283  3.90e-26    111.0     48.40
    498      450   715      35   319  3.93e-26    112.0     47.96
    499      453   736       9   278  4.00e-26    110.0     47.75
    500      457   714      41   300  4.17e-26    111.0     49.08
    501      457   714      51   310  4.18e-26    111.0     49.08
    502      457   714      61   320  4.28e-26    112.0     49.08
    503      450   715      12   296  4.32e-26    110.0     48.30
    504      457   714      31   290  4.41e-26    111.0     49.08
    505      453   736      12   281  4.42e-26    110.0     47.75
    506      450   715      14   298  4.55e-26    111.0     48.30
    507      450   715      13   297  4.58e-26    110.0     48.30
    508      450   715       8   292  4.65e-26    110.0     48.30
    509      453   736       8   277  4.69e-26    110.0     47.75
    510      450   715      14   298  4.80e-26    110.0     48.30
    511      450   715      12   296  4.93e-26    110.0     48.30
    512      457   714      32   291  4.99e-26    110.0     49.08
    513      461   736      22   283  5.01e-26    110.0     48.40
    514      461   736      24   285  5.02e-26    110.0     48.40
    515      450   715      15   299  5.38e-26    110.0     48.30
    516      450   715      17   301  5.54e-26    110.0     48.30
    517      450   715      16   300  5.65e-26    110.0     48.30
    518      457   714      22   281  5.98e-26    110.0     49.08
    519      453   736       8   277  6.01e-26    110.0     46.76
    520      450   715      16   300  6.03e-26    110.0     48.30
    521      450   715      19   303  6.04e-26    110.0     48.30
    522      452   714      16   284  6.04e-26    110.0     50.18
    523      450   715      23   307  6.10e-26    110.0     48.30
    524      463   711      18   266  6.22e-26    109.0     50.95
    525      461   715      19   290  6.32e-26    110.0     51.06
    526      450   715      23   307  6.52e-26    110.0     48.30
    527      447   720       2   267  6.75e-26    109.0     49.47
    528      450   715      16   300  6.89e-26    110.0     47.96
    529      461   736      22   283  7.21e-26    110.0     47.37
    530      450   715      16   300  7.29e-26    110.0     48.30
    531      450   715      14   298  7.35e-26    110.0     47.96
    532      461   736      36   297  7.50e-26    110.0     48.40
    533      450   715      23   307  7.64e-26    110.0     48.30
    534      450   715      17   301  7.64e-26    110.0     47.96
    535      453   736       7   276  7.75e-26    109.0     46.76
    536      450   715      18   299  7.97e-26    110.0     47.44
    537      457   714      33   292  8.02e-26    110.0     49.45
    538      450   715      15   299  8.08e-26    110.0     47.96
    539      461   715      39   310  8.10e-26    110.0     50.53
    540      457   714      31   290  8.15e-26    110.0     48.71
    541      461   736      22   283  8.29e-26    110.0     47.37
    542      450   715      23   307  8.33e-26    110.0     48.30
    543      453   736       6   275  8.67e-26    109.0     46.76
    544      461   736      23   284  9.76e-26    110.0     48.04
    545      461   715      14   285  9.90e-26    109.0     50.53
    546      461   721      10   271  1.00e-25    111.0     48.71
    547      461   736      42   303  1.02e-25    110.0     47.37
    548      461   715      18   289  1.03e-25    109.0     50.53
    549      447   715       2   287  1.07e-25    109.0     49.49
    550      461   715      19   290  1.09e-25    109.0     50.53
    551      461   715      19   290  1.14e-25    109.0     50.71
    552      457   714      28   287  1.14e-25    109.0     48.71
    553      457   714      32   291  1.29e-25    109.0     49.08
    554      457   714      32   291  1.30e-25    109.0     49.08
    555      457   714      91   350  1.36e-25    111.0     49.08
    556      457   714      32   291  1.44e-25    109.0     48.71
    557      450   725      41   349  1.44e-25    110.0     46.08
    558      457   714      37   296  1.45e-25    109.0     49.08
    559      450   714      13   296  1.47e-25    109.0     47.96
    560      450   715      12   284  1.48e-25    109.0     48.58
    561      426   715      17   311  1.48e-25    109.0     48.11
    562      457   714      33   292  1.48e-25    109.0     49.08
    563      450   725      37   345  1.48e-25    110.0     46.08
    564      457   714      30   289  1.55e-25    109.0     49.08
    565      442   715       4   294  1.57e-25    109.0     49.01
    566      461   715      25   296  1.58e-25    109.0     50.71
    567      452   720       4   264  1.58e-25    108.0     50.00
    568      450   715      15   299  1.61e-25    109.0     47.96
    569      461   715      25   296  1.64e-25    109.0     50.71
    570      461   715      37   308  1.64e-25    109.0     50.53
    571      446   653      10   214  1.71e-25    108.0     53.52
    572      461   715      29   300  1.82e-25    109.0     50.53
    573      451   710      15   290  1.84e-25    109.0     48.26
    574      450   715      12   284  1.84e-25    108.0     48.58
    575      442   715       2   292  1.86e-25    108.0     49.01
    576      447   715       2   287  1.94e-25    108.0     49.49
    577      461   715      39   310  1.97e-25    109.0     50.71
    578      450   715     188   441  1.97e-25    112.0     48.72
    579      461   715      19   290  2.08e-25    108.0     50.18
    580      442   715       1   291  2.21e-25    108.0     48.68
    581      450   714      24   305  2.24e-25    109.0     47.95
    582      450   715      64   348  2.25e-25    110.0     48.30
    583      457   714    1078  1337  2.25e-25    114.0     49.08
    584      450   714      23   305  2.35e-25    109.0     47.78
    585      458   653      50   243  2.36e-25    109.0     53.96
    586      446   715      22   308  2.51e-25    108.0     49.33
    587      461   715      15   286  2.58e-25    108.0     50.18
    588      461   715      19   290  2.59e-25    108.0     49.82
    589      461   715      19   290  2.64e-25    108.0     50.18
    590      461   715      37   308  2.66e-25    108.0     50.71
    591      461   715      14   285  2.68e-25    108.0     50.18
    592      461   715      22   293  2.70e-25    108.0     50.18
    593      446   715      23   309  2.73e-25    108.0     49.33
    594      450   714      13   296  2.76e-25    108.0     47.96
    595      450   714      13   296  2.79e-25    108.0     47.96
    596      461   715      19   290  2.79e-25    108.0     49.82
    597      458   653      23   216  2.86e-25    108.0     53.96
    598      446   715       2   288  2.96e-25    108.0     49.66
    599      442   715       4   294  2.97e-25    108.0     48.34
    600      450   714      22   301  2.98e-25    108.0     48.28
    601      450   715      18   303  3.01e-25    108.0     46.80
    602      452   717      16   287  3.06e-25    108.0     49.65
    603      445   715       4   291  3.33e-25    108.0     49.16
    604      461   715      20   291  3.33e-25    108.0     50.18
    605      461   715      47   318  3.35e-25    108.0     50.18
    606      461   715      12   283  3.38e-25    108.0     50.18
    607      447   715      30   310  3.49e-25    108.0     45.70
    608      447   715      38   324  3.51e-25    108.0     45.45
    609      450   714      22   305  3.53e-25    108.0     47.62
    610      450   715      20   305  3.53e-25    108.0     46.80
    611      450   715      17   302  3.53e-25    108.0     46.80
    612      447   715      35   321  3.56e-25    108.0     45.45
    613      461   715      12   283  3.64e-25    107.0     50.71
    614      461   715      16   287  3.76e-25    108.0     50.18
    615      461   715      17   288  3.78e-25    108.0     50.18
    616      452   714      16   284  3.81e-25    107.0     49.82
    617      463   711      34   281  3.82e-25    107.0     50.95
    618      453   711      22   285  3.86e-25    108.0     50.90
    619      447   715      36   322  3.90e-25    108.0     45.45
    620      447   715      41   327  4.03e-25    108.0     45.45
    621      444   715      17   305  4.03e-25    108.0     49.33
    622      457   714      30   289  4.10e-25    108.0     49.08
    623      450   714      22   301  4.14e-25    108.0     48.28
    624      461   736      36   297  4.14e-25    108.0     46.67
    625      450   715      30   314  4.16e-25    108.0     48.64
    626      444   715      17   305  4.38e-25    108.0     49.33
    627      453   711      13   276  4.64e-25    108.0     50.90
    628      449   705      22   293  4.66e-25    108.0     50.35
    629      462   710      25   284  4.66e-25    107.0     46.84
    630      453   711      15   278  4.86e-25    108.0     50.90
    631      453   711      13   276  4.95e-25    107.0     50.90
    632      453   711      12   275  5.09e-25    107.0     50.90
    633      452   720       4   264  5.14e-25    107.0     49.64
    634      462   710      25   284  5.21e-25    107.0     46.84
    635      462   710      25   284  5.21e-25    107.0     46.84
    636      450   714      22   305  5.29e-25    108.0     47.62
    637      463   711      34   281  5.33e-25    107.0     50.95
    638      452   720       2   262  5.36e-25    106.0     49.64
    639      461   715      36   307  5.53e-25    108.0     50.18
    640      462   710      24   283  5.62e-25    107.0     46.84
    641      458   653      23   216  5.70e-25    107.0     53.96
    642      450   715      18   302  5.71e-25    107.0     48.64
    643      463   711      34   281  5.91e-25    107.0     50.95
    644      453   711      11   274  5.92e-25    107.0     50.90
    645      462   710      22   281  5.94e-25    107.0     46.84
    646      463   711      35   282  6.09e-25    107.0     50.95
    647      450   714      13   296  6.10e-25    107.0     47.62
    648      450   714      24   307  6.21e-25    107.0     47.62
    649      463   721      29   289  6.31e-25    107.0     47.79
    650      450   715      18   302  6.33e-25    107.0     50.00
    651      450   715      30   314  6.36e-25    108.0     48.30
    652      463   710      17   275  6.36e-25    107.0     47.01
    653      450   715      30   314  6.48e-25    108.0     48.30
    654      463   710      16   274  6.60e-25    107.0     47.01
    655      463   710      18   276  6.62e-25    107.0     47.01
    656      450   715      30   314  6.86e-25    108.0     48.30
    657      450   715     471   755  6.90e-25    112.0     47.96
    658      450   715      31   315  7.18e-25    107.0     50.00
    659      450   715      34   319  7.36e-25    108.0     46.80
    660      450   715      30   314  7.39e-25    107.0     48.30
    661      463   711      34   281  7.39e-25    107.0     50.95
    662      445   715      22   309  7.42e-25    107.0     49.16
    663      450   714      22   305  7.43e-25    107.0     47.62
    664      447   715      38   318  7.57e-25    107.0     45.70
    665      450   725      40   348  7.78e-25    108.0     45.77
    666      450   725      37   345  7.91e-25    108.0     45.77
    667      463   711      22   269  7.92e-25    106.0     50.95
    668      462   710      31   290  8.08e-25    107.0     46.84
    669      461   712      30   294  8.25e-25    108.0     51.48
    670      450   715      32   317  8.25e-25    107.0     47.81
    671      457   714      30   289  8.81e-25    107.0     48.71
    672      450   715      36   321  8.94e-25    107.0     46.80
    673      450   725      41   349  9.01e-25    108.0     45.77
    674      462   710      25   284  9.74e-25    107.0     46.84
    675      450   714      22   305  9.91e-25    107.0     47.62
    676      445   715       3   290  1.00e-24    107.0     49.16
    677      450   714      14   297  1.01e-24    107.0     47.62
    678      450   715      30   314  1.03e-24    107.0     48.64
    679      462   710      38   297  1.03e-24    107.0     46.84
    680      450   715      30   314  1.04e-24    107.0     48.30
    681      450   715      41   326  1.04e-24    107.0     46.80
    682      462   710      41   300  1.06e-24    107.0     46.84
    683      450   715      36   321  1.08e-24    107.0     46.80
    684      450   725      60   368  1.09e-24    108.0     45.77
    685      450   715      41   326  1.27e-24    107.0     46.80
    686      450   715      30   314  1.32e-24    107.0     48.30
    687      232   283       1    52  1.33e-24     98.6     82.69
    688      450   715      30   314  1.48e-24    107.0     48.30
    689      453   704      20   263  1.54e-24    106.0     50.97
    690      450   715       8   292  1.58e-24    106.0     47.96
    691      450   725      41   349  1.59e-24    107.0     45.45
    692      450   715      30   314  1.65e-24    107.0     48.30
    693      463   711      23   270  1.68e-24    105.0     49.81
    694      448   711      16   271  1.75e-24    106.0     49.63
    695      450   714      59   342  1.81e-24    107.0     47.62
    696      447   715      38   322  1.96e-24    107.0     46.28
    697      448   711      16   271  1.96e-24    106.0     49.63
    698      445   721       2   295  2.00e-24    105.0     48.51
    699      447   715      38   322  2.08e-24    106.0     47.64
    700      450   715      30   314  2.14e-24    106.0     48.30
    701      450   715      36   320  2.76e-24    106.0     48.30
    702      461   711      32   276  2.98e-24    105.0     48.44
    703      448   711      37   292  3.00e-24    106.0     49.63
    704      448   711      37   292  3.05e-24    106.0     49.63
    705      448   711      37   292  3.48e-24    105.0     49.63
    706      450   715      30   314  3.59e-24    105.0     48.30
    707      450   715      30   314  3.59e-24    105.0     48.30
    708      452   720       4   264  3.61e-24    104.0     49.28
    709      450   715      18   302  3.76e-24    105.0     47.28
    710      462   710      31   290  3.85e-24    105.0     46.84
    711      450   715      30   314  4.13e-24    105.0     48.30
    712      450   715      19   303  4.29e-24    105.0     47.28
    713      450   715      10   294  4.55e-24    105.0     47.28
    714      450   715      22   306  4.66e-24    105.0     47.96
    715      450   715      22   306  4.66e-24    105.0     47.96
    716      450   715      19   303  4.89e-24    105.0     47.96
    717      450   715      30   314  5.14e-24    105.0     47.96
    718      450   715      17   301  5.23e-24    105.0     47.96
    719      450   715      22   306  5.67e-24    105.0     47.96
    720      450   715      30   314  5.74e-24    105.0     47.96
    721      450   715      19   303  6.31e-24    104.0     47.28
    722      458   653      22   213  6.39e-24    104.0     55.22
    723      452   720      24   284  6.59e-24    104.0     49.28
    724      399   704     139   430  6.77e-24    107.0     49.52
    725      450   715      14   298  6.83e-24    104.0     47.28
    726      451   717       9   271  6.89e-24    103.0     48.56
    727      450   715      30   314  7.05e-24    105.0     47.96
    728      450   715      30   314  7.60e-24    105.0     47.96
    729      450   715      76   360  7.82e-24    105.0     47.96
    730      450   715      17   301  8.70e-24    104.0     47.96
    731      450   715      30   314  8.92e-24    104.0     47.96
    732      450   715      14   298  8.96e-24    104.0     47.28
    733      450   715      30   314  9.44e-24    104.0     47.96
    734      447   715      42   346  9.51e-24    105.0     45.22
    735      450   715      30   314  9.62e-24    104.0     48.30
    736      399   704     152   443  1.10e-23    107.0     49.52
    737      450   715      30   314  1.32e-23    104.0     47.96
    738      447   715      23   316  1.40e-23    103.0     43.75
    739      450   715      30   314  1.60e-23    103.0     47.62
    740      447   715      31   324  1.85e-23    103.0     43.75
    741      453   704      18   261  1.86e-23    103.0     51.75
    742      447   715      42   335  2.01e-23    103.0     45.21
    743      448   724      12   294  2.08e-23    102.0     48.31
    744      447   715      42   335  2.62e-23    103.0     43.75
    745      453   711       4   253  2.79e-23    102.0     51.89
    746      448   728      29   319  2.81e-23    103.0     47.37
    747      456   722      24   300  3.48e-23    102.0     47.39
    748      463   704      24   264  3.53e-23    102.0     51.16
    749      399   704     220   511  3.84e-23    105.0     49.52
    750      448   724      10   296  3.88e-23    102.0     47.67
    751      450   715      64   348  3.90e-23    104.0     47.62
    752      453   705      13   259  4.32e-23    102.0     49.22
    753      460   721      15   284  4.35e-23    102.0     48.38
    754      453   705      52   298  4.47e-23    103.0     49.22
    755      453   711      27   276  4.61e-23    102.0     52.27
    756      448   724      10   296  4.91e-23    102.0     47.67
    757      448   724      10   296  4.96e-23    102.0     47.67
    758      453   711      20   269  5.14e-23    102.0     52.27
    759      448   724       8   294  5.82e-23    101.0     47.67
    760      453   711      12   261  5.87e-23    101.0     52.27
    761      448   724      28   314  6.03e-23    102.0     47.67
    762      453   711      12   261  6.04e-23    101.0     52.27
    763      453   711      30   279  6.49e-23    101.0     52.27
    764      448   724      27   313  7.56e-23    101.0     47.67
    765      463   704      24   264  7.93e-23    100.0     51.16
    766      452   743      13   316  8.20e-23    102.0     46.89
    767      448   724      29   315  8.60e-23    101.0     47.67
    768      453   711      25   274  8.65e-23    101.0     51.52
    769      463   704      24   264  8.80e-23    100.0     51.16
    770      231   656      99   542  9.33e-23    105.0     41.99
    771      232   285      35    92  1.20e-02     40.8     44.83
    772      231   656      99   542  1.00e-22    105.0     41.99
    773      232   285      35    92  1.20e-02     40.8     44.83
    774      453   711      17   266  1.01e-22    100.0     49.05
    775      453   711       5   254  1.04e-22    101.0     51.52
    776      451   710      43   326  1.05e-22    102.0     44.11
    777      453   711       5   254  1.08e-22    101.0     51.52
    778      453   711      18   267  1.09e-22    100.0     49.05
    779      462   717      12   283  1.10e-22    100.0     46.62
    780      231   656      99   542  1.20e-22    105.0     42.33
    781      232   285      35    92  6.00e-03     42.0     44.83
    782      462   717      20   291  1.36e-22    100.0     46.62
    783      463   719      26   281  1.41e-22    101.0     47.55
    784      463   719      27   282  1.42e-22    101.0     47.55
    785      460   721      19   288  1.46e-22    100.0     48.38
    786      463   711      39   302  1.50e-22    100.0     45.65
    787      463   719      27   282  1.58e-22    101.0     47.55
    788      460   721      19   288  1.61e-22    100.0     48.38
    789      463   711      22   285  1.64e-22    100.0     45.65
    790      463   704      36   276  1.64e-22    100.0     50.97
    791      463   719     718   973  1.69e-22    105.0     46.99
    792      463   704      23   263  1.77e-22    100.0     50.19
    793      462   717      12   283  1.79e-22    100.0     46.62
    794      463   719      25   280  1.82e-22    100.0     47.55
    795      463   719      29   284  1.88e-22    100.0     47.55
    796      463   719      24   279  1.93e-22    100.0     47.55
    797      462   717      14   285  1.95e-22    100.0     46.62
    798      460   676      47   262  1.96e-22    100.0     49.33
    799      462   717      17   288  1.97e-22    100.0     46.59
    800      463   704      33   273  1.99e-22    101.0     51.16
    801      463   719      28   283  2.05e-22    100.0     47.55
    802      462   717      19   290  2.06e-22    100.0     46.62
    803      463   719      25   280  2.06e-22    100.0     47.55
    804      460   676      47   262  2.07e-22    100.0     49.33
    805      463   719      30   285  2.07e-22    100.0     47.55
    806      463   719      20   275  2.14e-22    100.0     47.55
    807      463   719      23   278  2.15e-22    100.0     47.55
    808      463   719      25   280  2.16e-22    100.0     47.55
    809      460   676      47   262  2.19e-22    100.0     49.33
    810      463   719      26   281  2.22e-22    100.0     47.55
    811      463   719      25   280  2.29e-22    100.0     47.55
    812      463   719      26   281  2.29e-22    100.0     47.55
    813      463   719      41   296  2.35e-22    100.0     47.55
    814      463   719      24   279  2.35e-22    100.0     47.55
    815      463   719      26   281  2.37e-22    100.0     47.55
    816      463   719      24   279  2.40e-22    100.0     47.55
    817      463   716      49   316  2.43e-22    100.0     46.86
    818      463   719      26   281  2.43e-22    100.0     47.55
    819      460   676      31   246  2.48e-22    100.0     49.33
    820      463   719      23   278  2.55e-22    100.0     47.55
    821      463   711      34   281  2.56e-22     99.8     48.28
    822      462   717      18   289  2.61e-22     99.8     46.26
    823      462   717      30   301  2.69e-22    100.0     46.62
    824      463   719      23   278  2.74e-22    100.0     47.55
    825      463   719      22   277  2.75e-22    100.0     47.55
    826      463   719      42   297  2.80e-22    100.0     47.92
    827      463   719      26   281  2.85e-22    100.0     47.55
    828      463   719      17   272  2.94e-22    100.0     47.55
    829      463   719      24   279  2.95e-22    100.0     47.55
    830      463   719      48   303  3.12e-22    100.0     47.55
    831      462   717      15   286  3.56e-22     99.4     46.26
    832      463   711      39   302  3.63e-22     99.8     45.65
    833      462   717      14   285  3.63e-22     99.0     46.26
    834      463   719      26   281  3.64e-22    100.0     47.55
    835      463   719      23   278  3.70e-22    100.0     47.92
    836      452   724      17   304  3.70e-22     99.4     50.34
    837      463   719      27   282  3.73e-22    100.0     47.55
    838      463   719      34   289  3.83e-22    100.0     47.55
    839      463   716      39   306  3.98e-22     99.8     46.86
    840      463   716      39   306  4.02e-22     99.8     46.86
    841      463   716      39   306  4.10e-22     99.8     46.86
    842      463   719      26   281  4.22e-22     99.8     47.55
    843      462   717      24   295  4.31e-22     99.8     46.59
    844      463   711      19   282  4.37e-22     99.4     45.65
    845      452   724      12   299  4.48e-22     99.4     50.34
    846      462   717      30   301  4.50e-22     99.8     46.59
    847      463   719     760  1015  4.54e-22    104.0     47.55
    848      463   719     760  1015  4.62e-22    104.0     47.55
    849      462   724      27   301  5.05e-22     99.4     50.18
    850      463   719      33   288  5.08e-22     99.8     47.55
    851      460   721       5   279  5.20e-22     99.0     48.24
    852      452   724      13   300  5.30e-22     99.4     50.34
    853      452   724      17   304  5.41e-22     99.4     50.34
    854      452   724      13   300  5.44e-22     99.4     49.66
    855      452   724       9   296  5.59e-22     99.0     50.00
    856      460   721       7   281  5.60e-22     99.0     48.24
    857      463   719      42   297  5.73e-22     99.8     47.55
    858      462   724      42   316  5.81e-22     99.4     50.18
    859      452   724      13   300  5.98e-22     99.4     49.66
    860      462   724      30   304  6.21e-22     99.0     50.18
    861      462   717      13   284  6.22e-22     98.6     46.26
    862      463   719      23   278  6.48e-22     99.4     47.55
    863      462   724      43   317  6.72e-22     99.4     49.47
    864      463   719      45   300  7.31e-22     99.8     47.55
    865      463   719      25   280  7.35e-22     99.4     47.55
    866      463   719      25   280  7.38e-22     99.4     47.55
    867      463   656      29   221  7.66e-22     99.0     49.50
    868      463   719      29   302  7.66e-22     98.6     47.69
    869      463   719      31   304  7.86e-22     98.6     47.69
    870      463   719      27   282  7.88e-22     99.0     47.17
    871      463   711      22   285  8.03e-22     98.6     45.65
    872      445   711      16   274  8.26e-22     98.6     49.46
    873      463   719      25   280  8.30e-22     99.0     47.55
    874      463   719      26   281  8.35e-22     99.0     47.55
    875      448   721       2   285  8.58e-22     97.8     47.47
    876      463   719      28   301  8.58e-22     98.2     47.69
    877      463   719      34   286  8.71e-22     99.0     47.51
    878      463   719      19   292  8.72e-22     98.2     47.69
    879      462   717      18   289  8.81e-22     98.6     46.26
    880      463   656      25   217  8.96e-22     98.6     49.50
    881      463   719      23   278  9.06e-22     99.0     47.55
    882      463   719      25   280  9.38e-22     99.0     46.79
    883      463   719      24   279  9.51e-22     99.0     46.79
    884      463   719      47   315  9.56e-22     98.6     49.27
    885      452   730       9   301  9.70e-22     98.2     47.25
    886      452   730      30   322  9.72e-22     98.6     47.25
    887      463   719      28   283  9.80e-22     98.6     46.79
    888      463   656      25   217  9.84e-22     98.6     49.50
    889      463   719      24   279  9.94e-22     98.6     47.17
    890      463   719      25   280  1.00e-21     99.0     47.17
    891      463   719      43   316  1.00e-21     98.6     47.69
    892      463   719     718   973  1.03e-21    102.0     46.62
    893      463   719      27   282  1.05e-21     98.6     47.17
    894      463   719      29   284  1.06e-21     99.0     47.17
    895      463   719      24   279  1.08e-21     98.6     47.17
    896      463   719      17   290  1.08e-21     97.8     47.69
    897      463   719      27   282  1.09e-21     98.6     47.55
    898      463   719      47   315  1.09e-21     98.6     49.27
    899      463   719      18   291  1.10e-21     97.8     47.69
    900      463   719      27   282  1.12e-21     97.8     47.17
    901      445   660      13   223  1.14e-21     97.8     52.65
    902      463   719      26   281  1.15e-21     98.6     47.17
    903      445   660      16   226  1.15e-21     98.2     52.65
    904      463   719      23   278  1.17e-21     98.6     47.17
    905      445   660      16   226  1.17e-21     97.8     52.65
    906      445   660       8   218  1.20e-21     97.8     52.65
    907      463   719      27   282  1.24e-21     98.6     47.17
    908      437   722       8   296  1.25e-21     97.8     46.67
    909      463   656      49   241  1.26e-21     98.6     49.50
    910      453   711      18   267  1.28e-21     97.8     48.67
    911      463   719      23   278  1.31e-21     98.2     47.17
    912      463   719      22   277  1.33e-21     98.2     47.17
    913      463   719      27   282  1.37e-21     98.6     47.17
    914      462   724      30   304  1.41e-21     97.8     49.10
    915      463   719      30   285  1.45e-21     98.6     47.17
    916      462   717      17   288  1.49e-21     97.4     46.26
    917      463   719      23   278  1.50e-21     98.2     47.17
    918      463   719      27   282  1.54e-21     98.2     47.17
    919      463   656      47   239  1.55e-21     98.2     49.50
    920      453   711      21   270  1.58e-21     97.4     50.57
    921      445   660      16   226  1.61e-21     97.4     52.94
    922      462   717      18   289  1.62e-21     97.8     45.88
    923      463   658      22   216  1.62e-21     97.4     49.02
    924      463   719      30   285  1.63e-21     98.2     47.17
    925      463   719      23   278  1.65e-21     98.2     47.17
    926      463   719      23   278  1.67e-21     98.2     47.17
    927      457   714      26   280  1.69e-21     97.4     47.58
    928      463   719      27   282  1.69e-21     98.2     47.17
    929      463   719      24   279  1.69e-21     98.2     47.17
    930      463   658      23   217  1.70e-21     97.4     49.02
    931      455   660      21   224  1.75e-21    100.0     51.90
    932      463   719      29   284  1.75e-21     98.2     47.17
    933      463   655      23   219  1.81e-21     97.4     55.45
    934      457   714      10   264  1.84e-21     96.7     47.58
    935      463   719      23   278  1.85e-21     97.8     47.55
    936      457   714      15   269  1.90e-21     97.1     47.58
    937      462   717      18   289  1.90e-21     97.1     46.26
    938      451   710     569   852  1.94e-21    101.0     45.08
    939      463   719      27   282  1.95e-21     97.8     47.17
    940      463   719      23   278  1.95e-21    101.0     47.55
    941      463   719     374   629  1.95e-21    101.0     47.55
    942      463   719      27   282  1.97e-21     97.8     47.17
    943      457   714      14   268  2.03e-21     96.7     47.58
    944      463   719      25   280  2.03e-21     97.8     47.17
    945      463   719      25   280  2.04e-21     97.8     47.17
    946      435   722      11   301  2.12e-21     97.4     46.37
    947      463   711      15   278  2.14e-21     97.1     44.93
    948      437   722      31   319  2.15e-21     97.8     46.67
    949      435   722      10   300  2.22e-21     97.4     46.37
    950      463   719      24   279  2.29e-21     97.8     47.17
    951      463   723      21   280  2.50e-21     97.4     46.30
    952      437   722      31   319  2.54e-21     97.4     46.67
    953      463   656      72   264  2.55e-21     98.2     49.50
    954      463   719      57   312  2.70e-21     98.2     47.17
    955      463   719      26   284  2.71e-21     97.4     47.01
    956      463   719      26   284  2.71e-21     97.4     47.01
    957      462   656      30   219  2.73e-21     96.7     55.28
    958      463   717      14   285  2.74e-21     96.7     47.67
    959      463   719      22   277  2.74e-21     97.4     47.17
    960      463   719      25   280  2.77e-21     97.4     47.17
    961      460   658      27   224  2.79e-21     96.7     50.72
    962      443   722      35   317  2.89e-21     97.4     46.60
    963      463   719     726   981  2.91e-21    101.0     46.79
    964      443   722      40   322  2.98e-21     97.4     46.60
    965      463   719     726   981  3.04e-21    101.0     46.79
    966      463   719      22   277  3.07e-21     97.4     47.17
    967      443   722      34   316  3.27e-21     97.4     46.60
    968      450   711       5   281  3.30e-21     96.3     47.16
    969      443   722      35   317  3.33e-21     97.4     46.60
    970      455   660      15   218  3.35e-21     97.1     51.43
    971      463   711      22   285  3.35e-21     96.7     45.29
    972      450   711      10   286  3.54e-21     96.7     47.16
    973      462   656      52   241  3.75e-21     97.1     55.28
    974      453   711      21   284  3.92e-21     96.7     48.35
    975      463   719      27   286  3.98e-21     97.1     46.64
    976      453   711      22   285  4.02e-21     97.1     48.35
    977      443   722      24   306  4.27e-21     96.7     46.60
    978      453   711      21   284  4.30e-21     96.7     48.35
    979      453   711      20   283  4.36e-21     96.7     48.35
    980      443   722      11   293  4.41e-21     96.3     46.60
    981      455   694      12   247  4.43e-21     96.7     48.40
    982      453   711      22   285  4.45e-21     97.1     48.35
    983      463   719      27   285  4.47e-21     97.1     47.78
    984      443   722       8   290  4.51e-21     96.3     46.60
    985      443   722      31   313  4.53e-21     96.7     46.60
    986      443   722      10   292  4.59e-21     96.3     46.60
    987      443   722       7   289  4.69e-21     96.3     46.60
    988      463   729      19   284  4.92e-21     95.9     46.32
    989      453   711      24   287  5.00e-21     96.7     48.35
    990      463   715      17   284  5.18e-21     95.9     45.00
    991      463   715      16   283  5.23e-21     95.9     45.00
    992      443   722      17   299  5.27e-21     96.3     46.60
    993      463   723     715   974  5.52e-21    100.0     46.30
    994      441   722       8   292  5.97e-21     95.9     46.62
    995      463   719      27   282  6.02e-21     96.7     46.79
    996      463   719      25   283  6.11e-21     96.7     47.78
    997      455   694      12   247  6.19e-21     96.3     48.40
    998      448   724      12   298  6.21e-21     95.9     47.67
    999      460   711      33   277  6.40e-21     95.9     46.69
    1000     463   719      31   289  6.45e-21     96.7     47.41
    1001     448   724      27   313  6.51e-21     95.9     47.67
    1002     441   722      21   305  6.85e-21     95.9     46.62
    1003     463   719      25   280  6.97e-21     96.3     46.79
    1004     453   711      21   284  7.09e-21     95.9     48.35
    1005     463   723      39   298  7.65e-21     96.3     45.56
    1006     463   719      30   285  7.88e-21     96.3     46.79
    1007     455   694       7   242  8.09e-21     95.9     50.00
    1008     463   723      21   280  8.64e-21     95.9     45.56
    1009     463   723      21   280  8.80e-21     95.9     45.56
    1010     451   734      17   318  8.88e-21     95.9     43.34
    1011     463   719      25   283  8.96e-21     95.9     47.78
    1012     463   711      22   285  9.00e-21     95.5     45.29
    1013     439   669      13   237  9.10e-21     95.5     50.00
    1014     463   711      22   285  9.15e-21     95.5     45.29
    1015     463   719      26   281  9.28e-21     95.9     46.79
    1016     463   723      42   301  9.34e-21     95.9     45.56
    1017     451   734      12   313  9.45e-21     95.9     43.34
    1018     463   711      31   284  9.68e-21     95.5     48.89
    1019     462   658      30   221  9.69e-21     95.1     52.45
    1020     455   694      12   247  9.78e-21     95.9     48.40
    1021     455   694      15   250  9.88e-21     95.9     48.40
    1022     439   669      13   237  9.99e-21     95.5     50.00
    1023     461   653      21   214  1.02e-20     95.5     53.54
    1024     463   719      24   279  1.02e-20     95.9     46.79
    1025     455   694      12   247  1.03e-20     95.9     48.40
    1026     463   656      44   236  1.03e-20     95.5     49.01
    1027     463   719      27   282  1.07e-20     95.9     46.79
    1028     463   655      24   220  1.11e-20     95.1     53.00
    1029     463   723      21   280  1.14e-20     95.5     45.56
    1030     439   669      13   237  1.14e-20     95.1     50.00
    1031     463   719      24   282  1.15e-20     95.5     47.78
    1032     463   719      25   280  1.15e-20     95.5     46.79
    1033     463   719      28   283  1.18e-20     95.5     46.79
    1034     455   718      20   287  1.26e-20     95.9     46.89
    1035     463   719      22   277  1.27e-20     95.5     46.79
    1036     463   719      25   283  1.32e-20     95.5     47.41
    1037     463   711      27   280  1.35e-20     95.1     48.89
    1038     462   658      52   243  1.36e-20     95.1     52.45
    1039     463   719      35   294  1.44e-20     95.9     47.60
    1040     453   711       6   255  1.60e-20     94.4     50.19
    1041     463   719      25   284  1.61e-20     95.5     47.60
    1042     439   669       3   227  1.82e-20     94.4     49.79
    1043     441   722      17   301  1.83e-20     94.4     46.30
    1044     463   719      44   299  1.85e-20     95.5     46.79
    1045     441   722      15   299  1.86e-20     94.4     46.30
    1046     441   720       5   287  1.99e-20     94.0     47.25
    1047     461   711      56   311  2.13e-20     95.1     48.90
    1048     463   726      45   320  2.37e-20     94.7     44.41
    1049     442   663       1   220  2.42e-20     94.7     51.10
    1050     443   722      23   305  2.45e-20     94.4     46.28
    1051     460   711      14   258  2.52e-20     93.6     47.67
    1052     157   227       3    73  2.56e-20     87.0     71.83
    1053     442   663       1   220  2.64e-20     94.4     50.66
    1054     463   711     361   624  2.67e-20     97.8     44.89
    1055     442   663       1   220  2.81e-20     94.4     50.66
    1056     455   663      15   221  3.35e-20     94.4     51.40
    1057     157   227      10    80  3.70e-20     87.0     71.83
    1058     157   228       5    79  3.72e-20     87.0     72.00
    1059     463   719      27   282  3.77e-20     94.0     47.92
    1060     463   726      20   295  3.88e-20     93.6     44.41
    1061     450   711       6   282  4.03e-20     94.0     45.74
    1062     459   669       5   209  4.60e-20     92.8     51.42
    1063     463   714      18   268  4.71e-20     92.4     46.67
    1064     451   734      15   322  5.12e-20     93.6     42.55
    1065     451   734       6   313  5.36e-20     93.2     42.55
    1066     450   711      13   289  5.37e-20     94.0     45.74
    1067     442   663       1   220  5.56e-20     93.6     50.66
    1068     463   719      29   284  5.74e-20     93.6     47.92
    1069     442   663       1   220  5.82e-20     93.6     50.66
    1070     157   227       4    77  5.90e-20     86.3     71.62
    1071     157   227       5    78  6.06e-20     86.3     71.62
    1072     455   663      33   239  6.35e-20     93.6     51.40
    1073     463   719      27   282  6.35e-20     93.6     47.92
    1074     157   227       7    80  6.42e-20     86.3     71.62
    1075     451   734      41   348  7.28e-20     93.6     42.55
    1076     451   734      42   349  7.31e-20     93.6     42.55
    1077     157   227       8    81  8.02e-20     85.9     71.62
    1078     457   696       9   278  9.52e-20     92.4     49.28
    1079     466   714      27   282  9.82e-20     92.0     50.76
    1080     451   734      28   335  9.86e-20     93.2     42.55
    1081     151   227       5    86  1.21e-19     85.5     68.29
    1082     452   722      10   292  1.23e-19     92.0     49.49
    1083     455   663       5   211  1.32e-19     92.4     51.40
    1084     457   696       6   275  1.40e-19     91.7     48.74
    1085     457   696      44   313  1.45e-19     92.4     48.74
    1086     457   696       6   275  1.46e-19     91.7     48.74
    1087     457   696       8   277  1.51e-19     91.7     48.74
    1088     457   696      45   314  1.52e-19     92.4     48.74
    1089     457   696       8   277  1.55e-19     91.7     48.74
    1090     457   696       5   274  1.56e-19     91.7     48.74
    1091     441   669      13   235  1.58e-19     92.0     50.64
    1092     453   669       9   219  1.60e-19     91.7     51.58
    1093     457   696      11   280  1.61e-19     91.7     48.74
    1094     457   696      31   300  1.67e-19     92.0     48.74
    1095     457   696       7   276  1.69e-19     91.7     48.74
    1096     457   696       8   277  1.71e-19     91.7     48.74
    1097     453   669       9   219  1.72e-19     91.7     51.33
    1098     457   696       8   277  1.79e-19     91.7     48.74
    1099     441   704       8   274  1.92e-19     91.3     46.76
    1100     441   704       5   271  2.06e-19     91.3     46.42
    1101     457   696       9   278  2.07e-19     91.7     48.74
    1102     455   699      10   259  2.22e-19     92.0     48.84
    1103     455   699      18   267  2.40e-19     92.0     48.84
    1104     455   694      14   249  2.42e-19     91.7     49.80
    1105     157   228      21    95  2.48e-19     85.5     72.00
    1106     455   694      15   250  2.85e-19     91.7     49.80
    1107     455   663       7   213  2.95e-19     91.3     50.69
    1108     463   740      25   316  3.10e-19     91.3     47.18
    1109     462   723      16   281  3.29e-19     90.5     48.23
    1110     451   710       5   293  3.47e-19     90.9     44.97
    1111     157   227       8    81  3.92e-19     84.0     70.27
    1112     442   663       1   220  4.26e-19     90.9     49.78
    1113     157   227       8    81  4.45e-19     84.0     70.27
    1114     457   722      17   297  4.88e-19     90.1     46.62
    1115     461   660      15   213  6.33e-19     90.5     54.15
    1116     463   712      25   289  6.48e-19     89.7     47.79
    1117     463   712      25   289  6.61e-19     89.7     47.79
    1118     451   734      28   335  6.78e-19     90.5     43.20
    1119     461   649      21   240  7.19e-19     90.1     49.78
    1120     463   712      50   314  7.24e-19     90.5     47.79
    1121     461   660      15   213  8.48e-19     90.1     54.15
    1122     453   669       9   219  8.76e-19     89.4     51.33
    1123     463   664      30   239  9.20e-19     89.4     49.07
    1124     459   711      17   263  1.14e-18     88.6     47.67
    1125     456   696      38   308  1.17e-18     89.7     48.93
    1126     456   707      23   273  1.19e-18     89.7     47.10
    1127     450   659       4   242  1.26e-18     89.0     48.16
    1128     459   711      11   257  1.55e-18     87.8     47.67
    1129     453   669       8   218  1.59e-18     89.7     51.58
    1130     459   711      15   261  1.65e-18     88.2     47.67
    1131     459   711      11   257  1.74e-18     87.8     47.67
    1132     453   669      24   234  1.79e-18     89.7     51.58
    1133     449   714      13   285  2.04e-18     88.6     46.08
    1134     463   668      46   257  2.06e-18     89.4     49.77
    1135     450   659      33   271  2.15e-18     89.0     48.16
    1136     463   712      29   277  2.37e-18     88.2     49.04
    1137     449   714      23   295  2.40e-18     88.6     46.58
    1138     463   704      27   267  2.45e-18     88.6     49.80
    1139     463   704      26   266  2.55e-18     87.8     49.80
    1140     450   732      31   294  2.68e-18     88.2     45.39
    1141     463   668      38   249  2.69e-18     88.6     49.31
    1142     457   722      17   297  2.91e-18     88.2     46.28
    1143     454   707       9   260  3.35e-18     87.4     46.01
    1144     450   732       5   268  3.36e-18     87.0     45.05
    1145     450   732       6   269  3.49e-18     87.4     44.93
    1146     450   732       6   269  3.62e-18     87.0     45.05
    1147     450   659      42   280  3.77e-18     88.6     48.16
    1148     450   659      43   281  3.82e-18     88.6     48.16
    1149     462   653      29   233  4.05e-18     88.2     50.95
    1150     450   659      40   278  4.19e-18     88.2     48.16
    1151     447   732       5   271  4.19e-18     86.7     45.27
    1152     463   712      43   291  4.23e-18     88.2     49.04
    1153     463   704      13   247  4.41e-18     87.8     49.20
    1154     463   719     880  1153  4.60e-18     91.3     46.62
    1155     463   704      13   247  5.20e-18     87.8     49.20
    1156     461   649      34   253  5.24e-18     87.8     49.78
    1157     450   732       5   268  5.33e-18     86.7     45.39
    1158     463   704      13   247  5.47e-18     87.8     49.20
    1159     463   704      13   247  5.55e-18     87.8     49.20
    1160     463   704      16   250  5.68e-18     87.8     49.20
    1161     463   712      30   278  5.70e-18     87.0     48.66
    1162     457   656      20   219  5.82e-18     87.8     49.03
    1163     463   712      28   276  5.84e-18     87.0     48.66
    1164     450   659      66   304  5.90e-18     88.2     48.16
    1165     450   732      31   294  5.93e-18     87.0     45.05
    1166     463   704      18   252  6.36e-18     87.8     49.20
    1167     457   656      21   220  6.40e-18     87.8     49.03
    1168     450   732      34   297  6.49e-18     87.0     45.39
    1169     463   704      16   256  6.55e-18     86.3     49.41
    1170     463   712      30   278  6.68e-18     86.7     48.66
    1171     457   656       7   206  6.77e-18     87.4     49.03
    1172     450   732      10   273  7.26e-18     86.3     44.93
    1173     461   688     137   381  7.26e-18     90.5     48.81
    1174     450   732       6   269  7.42e-18     86.3     45.39
    1175     447   732       5   271  7.64e-18     86.3     44.93
    1176     463   704      13   247  7.65e-18     87.4     49.20
    1177     463   704      21   261  7.83e-18     86.7     49.41
    1178     463   704      17   257  8.12e-18     85.9     49.41
    1179     442   732       2   273  8.24e-18     85.9     44.52
    1180     461   720      23   276  8.30e-18     86.3     46.35
    1181     442   732       1   272  8.33e-18     85.9     44.52
    1182     463   704     160   394  8.42e-18     89.0     49.20
    1183     450   732      29   292  8.43e-18     86.7     44.93
    1184     461   720      31   284  8.46e-18     86.3     46.35
    1185     442   732       1   272  8.57e-18     85.9     44.52
    1186     450   732       6   269  8.62e-18     86.3     45.39
    1187     442   732       1   272  8.98e-18     85.9     44.52
    1188     449   714      20   292  9.33e-18     87.0     46.08
    1189     450   732      10   273  9.35e-18     86.3     44.93
    1190     461   720      23   276  9.43e-18     86.3     46.35
    1191     454   707     335   586  9.51e-18     89.4     46.01
    1192     463   704      17   257  9.51e-18     85.9     49.41
    1193     461   720      23   276  9.58e-18     86.3     46.35
    1194     463   704      21   261  9.63e-18     85.9     49.41
    1195     461   720      33   286  9.68e-18     86.3     46.35
    1196     450   732      10   273  9.89e-18     85.9     44.93
    1197     457   656      19   218  1.01e-17     87.0     48.54
    1198     450   732      10   273  1.02e-17     85.9     44.93
    1199     463   655      50   246  1.02e-17     87.0     50.99
    1200     461   720      33   286  1.06e-17     86.3     46.35
    1201     461   720      21   274  1.07e-17     85.9     46.35
    1202     463   717      92   369  1.08e-17     88.2     45.61
    1203     463   704     157   391  1.10e-17     88.2     49.20
    1204     450   732      10   273  1.12e-17     85.9     44.59
    1205     450   732       8   271  1.12e-17     85.9     44.93
    1206     449   714      19   291  1.13e-17     86.7     46.08
    1207     454   707     335   586  1.15e-17     89.4     46.01
    1208     450   711      31   281  1.17e-17     85.9     45.82
    1209     461   720      20   273  1.17e-17     85.9     46.35
    1210     450   732       5   268  1.22e-17     85.5     45.05
    1211     461   720      11   264  1.28e-17     85.5     46.35
    1212     450   732      22   285  1.31e-17     85.9     45.05
    1213     450   732       2   265  1.32e-17     85.1     45.39
    1214     461   720      17   270  1.36e-17     85.5     46.35
    1215     450   732       6   269  1.41e-17     85.5     45.05
    1216     462   723     202   467  1.42e-17     88.2     48.93
    1217     450   732       3   266  1.45e-17     85.5     45.05
    1218     450   732       7   270  1.45e-17     85.5     45.05
    1219     450   732       2   265  1.46e-17     85.1     44.93
    1220     461   649      20   236  1.47e-17     86.3     49.55
    1221     461   720      20   273  1.47e-17     85.5     46.35
    1222     450   732       9   272  1.50e-17     85.5     45.05
    1223     463   723     559   822  1.50e-17     89.7     48.91
    1224     450   732       3   266  1.53e-17     85.5     45.05
    1225     463   723     559   822  1.54e-17     89.4     48.91
    1226     450   732       4   267  1.55e-17     85.1     45.39
    1227     461   649      18   234  1.58e-17     86.3     49.55
    1228     441   711       1   255  1.59e-17     85.1     45.77
    1229     461   720      23   276  1.59e-17     85.5     46.35
    1230     450   732      10   273  1.61e-17     85.5     45.05
    1231     461   720      23   276  1.61e-17     85.5     46.35
    1232     450   732       6   269  1.63e-17     85.5     45.05
    1233     450   732      10   273  1.64e-17     85.5     45.05
    1234     450   732       4   267  1.65e-17     85.1     45.39
    1235     461   720      13   266  1.66e-17     85.1     46.35
    1236     450   732      10   273  1.70e-17     85.5     45.05
    1237     441   711       1   255  1.75e-17     84.7     45.77
    1238     463   723     559   822  1.77e-17     89.4     50.00
    1239     461   720      33   286  1.78e-17     85.5     45.99
    1240     463   723     552   815  1.78e-17     89.4     50.00
    1241     463   723     553   816  1.81e-17     89.4     50.00
    1242     450   732       7   270  1.81e-17     85.1     45.05
    1243     450   732      10   273  1.82e-17     85.1     44.71
    1244     461   720      13   266  1.84e-17     85.1     46.35
    1245     450   732       5   268  1.84e-17     85.1     45.05
    1246     450   732      10   273  1.85e-17     85.1     44.71
    1247     461   720      21   274  1.85e-17     85.1     45.99
    1248     450   732      10   273  1.88e-17     85.1     45.05
    1249     463   723     556   819  1.89e-17     89.4     50.00
    1250     450   732       5   268  1.91e-17     85.1     45.05
    1251     450   732       8   271  1.92e-17     85.1     45.05
    1252     461   720      23   276  1.93e-17     85.5     46.35
    1253     461   720     375   628  1.94e-17     88.6     46.35
    1254     450   732       3   266  1.96e-17     84.7     44.93
    1255     461   720      13   266  1.97e-17     85.1     46.35
    1256     461   720      22   275  1.97e-17     85.1     45.99
    1257     469   711      25   289  1.98e-17     85.5     44.21
    1258     450   732      24   287  1.99e-17     85.5     45.05
    1259     461   720     376   629  1.99e-17     88.6     46.35
    1260     450   732      10   273  2.03e-17     85.1     45.05
    1261     463   675      43   264  2.05e-17     85.9     46.70
    1262     450   732       5   268  2.18e-17     85.1     44.59
    1263     450   732       5   268  2.25e-17     84.7     45.05
    1264     463   723    1885  2148  2.39e-17     89.4     50.00
    1265     463   723    1885  2148  2.43e-17     89.4     50.00
    1266     450   732       5   268  2.50e-17     84.3     45.05
    1267     442   732       8   279  2.51e-17     84.7     44.19
    1268     455   715      18   308  2.51e-17     85.9     47.81
    1269     441   711       1   255  2.55e-17     84.7     45.91
    1270     450   732       4   267  2.55e-17     84.3     45.05
    1271     450   711       9   259  2.58e-17     84.3     45.59
    1272     450   732       5   268  2.64e-17     84.3     45.05
    1273     441   711       1   255  2.72e-17     84.3     45.91
    1274     463   723    1920  2183  2.73e-17     89.0     50.00
    1275     441   711       1   255  2.77e-17     84.3     45.91
    1276     456   649      30   237  2.92e-17     85.5     52.31
    1277     455   656       7   206  2.93e-17     85.5     50.97
    1278     461   720      16   269  3.11e-17     84.3     47.27
    1279     461   720      17   270  3.13e-17     84.3     47.27
    1280     450   711      31   281  3.18e-17     84.7     45.59
    1281     456   649      38   245  3.21e-17     85.5     52.31
    1282     456   649      39   246  3.24e-17     85.5     52.31
    1283     450   711      31   281  3.30e-17     84.7     45.59
    1284     450   711      31   281  3.40e-17     84.7     45.59
    1285     450   732       2   265  3.42e-17     84.0     45.05
    1286     463   723    1343  1606  3.53e-17     88.6     49.64
    1287     450   711      31   281  3.53e-17     84.7     45.59
    1288     456   649      15   222  3.57e-17     85.1     52.31
    1289     456   649       7   214  3.61e-17     84.7     52.31
    1290     463   675      43   264  3.70e-17     85.1     46.26
    1291     456   649       9   216  3.75e-17     84.7     52.31
    1292     463   740      73   375  4.15e-17     85.9     46.93
    1293     443   711      30   282  4.16e-17     84.7     47.08
    1294     450   711      31   281  4.54e-17     84.3     45.82
    1295     450   711       2   252  4.59e-17     83.6     45.09
    1296     459   703      38   289  4.64e-17     84.3     47.91
    1297     463   723     553   816  4.68e-17     87.8     48.55
    1298     450   711       1   251  4.87e-17     83.6     45.09
    1299     450   711       8   258  5.22e-17     83.6     45.09
    1300     463   675      43   264  5.36e-17     84.7     46.70
    1301     446   658       2   213  5.36e-17     83.6     49.32
    1302     450   711       8   258  5.47e-17     83.6     45.09
    1303     456   649      38   245  5.80e-17     84.7     52.31
    1304     456   649      38   245  5.85e-17     84.7     52.31
    1305     461   660      19   217  5.89e-17     84.7     52.20
    1306     461   660      16   214  6.24e-17     84.7     52.20
    1307     461   660      16   214  6.44e-17     84.7     52.20
    1308     461   660      19   217  6.45e-17     84.7     52.20
    1309     463   655      49   248  6.80e-17     84.3     50.98
    1310     463   723    1885  2148  6.83e-17     87.8     49.64
    1311     456   649      38   245  6.94e-17     84.3     52.31
    1312     461   660      16   214  7.22e-17     84.3     52.20
    1313     450   711       2   252  7.22e-17     83.2     45.09
    1314     456   649       9   216  7.70e-17     84.0     50.46
    1315     446   658       2   213  8.41e-17     83.2     47.71
    1316     461   660      16   214  8.55e-17     84.3     52.20
    1317     446   732       5   272  8.57e-17     83.2     44.15
    1318     457   722      11   305  9.15e-17     83.6     45.16
    1319     449   658      21   229  9.46e-17     86.3     49.54
    1320     461   660      16   214  9.85e-17     84.0     52.20
    1321     449   658       8   216  1.00e-16     85.9     49.54
    1322     449   658      27   235  1.05e-16     85.9     49.54
    1323     438   693      74   344  1.09e-16     84.3     47.86
    1324     463   723     553   816  1.11e-16     86.7     49.64
    1325     449   658      15   223  1.12e-16     85.9     49.54
    1326     442   667       3   229  1.14e-16     83.6     50.00
    1327     461   660      14   212  1.18e-16     83.6     52.20
    1328     442   667       5   231  1.25e-16     82.8     50.00
    1329     442   667       3   229  1.28e-16     82.8     50.00
    1330     442   667       3   229  1.37e-16     82.4     50.00
    1331     463   723    1885  2148  1.41e-16     86.7     49.64
    1332     463   707      28   281  1.42e-16     83.2     46.42
    1333     461   658      35   231  1.46e-16     84.0     52.22
    1334     463   723    1885  2148  1.46e-16     86.7     49.64
    1335     457   649      10   216  1.47e-16     82.8     52.09
    1336     462   723      28   293  1.57e-16     82.8     49.28
    1337     449   658       4   212  1.58e-16     82.4     49.54
    1338     461   658      35   231  1.60e-16     83.6     52.22
    1339     456   677      33   260  1.62e-16     83.2     51.05
    1340     449   658      22   230  1.62e-16     85.5     49.54
    1341     456   677      31   258  1.69e-16     83.2     51.05
    1342     456   677      25   252  1.73e-16     83.2     51.05
    1343     442   667      20   246  1.84e-16     82.4     50.00
    1344     449   658      10   218  1.85e-16     84.7     50.00
    1345     456   677      31   258  1.90e-16     82.8     51.05
    1346     446   732       5   272  1.94e-16     82.0     43.81
    1347     463   704      34   269  1.96e-16     83.6     49.00
    1348     456   677      28   255  2.10e-16     82.8     51.05
    1349     463   704      18   253  2.20e-16     83.2     49.42
    1350     463   657      28   222  2.21e-16     83.2     49.25
    1351     449   658       5   213  2.23e-16     82.0     47.91
    1352     457   649      32   238  2.26e-16     82.8     52.09
    1353     459   667      16   222  2.32e-16     82.8     50.69
    1354     463   704      16   251  2.42e-16     82.8     49.00
    1355     455   658       6   208  2.44e-16     84.7     50.94
    1356     455   658       6   208  2.46e-16     84.7     50.94
    1357     463   704      17   252  2.50e-16     82.8     49.42
    1358     463   704      13   248  2.66e-16     82.8     49.42
    1359     463   657      28   222  2.83e-16     82.8     48.76
    1360     456   677      26   253  2.88e-16     82.0     51.05
    1361     459   667      16   222  3.01e-16     82.4     50.69
    1362     463   704      16   251  3.21e-16     82.4     49.00
    1363     456   677      25   252  3.24e-16     82.0     51.05
    1364     461   658      35   231  3.32e-16     82.8     52.22
    1365     463   669      36   243  3.40e-16     82.0     50.68
    1366     450   670      47   269  3.72e-16     82.8     49.79
    1367     449   658      27   235  4.05e-16     83.6     50.46
    1368     455   660      11   218  4.10e-16     83.6     52.13
    1369     443   711       5   257  4.28e-16     80.9     47.08
    1370     455   658      14   216  4.29e-16     82.4     50.94
    1371     455   660      11   218  4.81e-16     83.2     52.13
    1372     455   658       4   206  4.90e-16     83.6     50.94
    1373     459   656      36   236  4.91e-16     81.6     51.21
    1374     556   735     770   944  5.20e-16     84.3     48.91
    1375     450   570     580   708  1.00e-03     44.3     49.24
    1376     455   658       4   206  5.21e-16     83.2     50.94
    1377     458   743       7   292  5.23e-16     81.6     47.99
    1378     455   658      14   216  5.35e-16     83.2     50.94
    1379     556   735     766   940  5.60e-16     84.3     48.91
    1380     450   570     576   704  1.00e-03     43.9     49.24
    1381     556   735     768   942  5.71e-16     84.3     48.91
    1382     450   570     578   706  1.00e-03     43.9     49.24
    1383     459   656      13   213  5.79e-16     80.9     51.69
    1384     455   658      14   216  5.90e-16     83.2     50.94
    1385     450   732      10   273  5.95e-16     80.9     44.03
    1386     463   701      14   275  7.23e-16     80.9     42.32
    1387     459   656      36   236  7.26e-16     80.9     51.21
    1388     463   701      14   275  7.29e-16     80.9     42.32
    1389     460   715    1246  1522  7.76e-16     84.3     47.28
    1390     455   658       8   209  7.79e-16     80.1     50.95
    1391     455   658       4   205  7.86e-16     80.1     50.95
    1392     455   658      13   214  8.06e-16     80.1     50.95
    1393     463   704     159   394  8.49e-16     82.4     49.19
    1394     462   669       9   216  8.65e-16     80.9     49.53
    1395     463   704     156   391  8.88e-16     82.4     49.19
    1396     460   715    1226  1502  8.95e-16     84.0     47.28
    1397     455   658      14   215  8.95e-16     80.1     50.95
    1398     460   715    1278  1554  9.15e-16     84.0     47.28
    1399     459   656      36   236  9.47e-16     80.5     51.21
    1400     463   707      28   281  9.87e-16     80.9     46.42
    1401     459   656      38   238  1.01e-15     80.5     51.21
    1402     459   656      36   236  1.05e-15     80.5     51.21
    1403     459   656      36   236  1.06e-15     80.5     51.21
    1404     448   657       5   207  1.08e-15     79.7     50.23
    1405     463   663      14   223  1.09e-15     79.7     49.55
    1406     448   657       7   209  1.10e-15     80.1     50.23
    1407     448   657       7   209  1.11e-15     80.1     50.23
    1408     463   704     170   405  1.11e-15     82.0     49.19
    1409     461   657      16   212  1.13e-15     80.1     50.97
    1410     459   656      36   236  1.14e-15     80.5     51.21
    1411     461   657      16   212  1.19e-15     80.1     50.97
    1412     448   657       7   209  1.20e-15     80.1     50.23
    1413     448   657       6   208  1.21e-15     79.7     50.23
    1414     448   657       6   208  1.22e-15     79.7     50.23
    1415     448   657       7   209  1.22e-15     79.7     50.23
    1416     448   657       7   209  1.24e-15     80.1     50.23
    1417     463   704     600   840  1.26e-15     83.2     49.41
    1418     448   656       6   207  1.31e-15     80.5     50.46
    1419     448   657       7   209  1.35e-15     80.1     50.68
    1420     448   656       7   208  1.39e-15     80.5     50.46
    1421     463   704     151   386  1.43e-15     81.6     49.19
    1422     447   657      11   220  1.47e-15     80.1     46.19
    1423     448   657       6   208  1.48e-15     79.3     50.23
    1424     448   657       7   209  1.48e-15     79.3     50.23
    1425     455   657       9   211  1.52e-15     81.6     52.40
    1426     448   657       6   208  1.53e-15     79.3     50.23
    1427     463   663      14   223  1.56e-15     81.6     48.86
    1428     500   710     233   451  1.56e-15     82.4     46.55
    1429     448   657       6   208  1.57e-15     79.3     50.23
    1430     448   657       6   208  1.58e-15     79.3     50.23
    1431     448   657       6   208  1.62e-15     79.3     50.23
    1432     463   704     148   383  1.64e-15     81.6     49.19
    1433     500   710     233   451  1.68e-15     82.0     46.55
    1434     463   657      27   221  1.69e-15     80.5     47.26
    1435     448   657       7   209  1.70e-15     79.0     50.23
    1436     463   704     148   383  1.71e-15     81.3     49.19
    1437     463   663      14   223  1.72e-15     81.3     48.86
    1438     448   657       6   208  1.74e-15     79.0     50.23
    1439     448   657       7   209  1.78e-15     79.3     50.23
    1440     448   657       6   208  1.90e-15     79.3     50.23
    1441     434   714      17   306  1.96e-15     79.7     41.91
    1442     455   656      13   217  2.03e-15     79.0     50.24
    1443     448   657       6   208  2.13e-15     79.0     50.23
    1444     447   657      11   220  2.49e-15     79.3     45.74
    1445     463   663      14   223  2.66e-15     78.6     49.32
    1446     463   656      23   217  2.70e-15     80.1     49.00
    1447     463   663      14   223  2.81e-15     78.6     49.32
    1448     501   672     164   330  2.91e-15     82.4     56.00
    1449     447   657      13   222  3.07e-15     78.6     45.74
    1450     423   712       3   271  3.27e-15     78.6     48.15
    1451     447   657      11   220  3.60e-15     79.0     45.74
    1452     447   657       8   217  3.82e-15     78.2     45.74
    1453     463   657      28   222  3.86e-15     79.3     46.77
    1454     447   657       9   218  3.96e-15     78.2     45.74
    1455     463   704      48   283  4.27e-15     79.7     44.76
    1456     447   657      11   220  4.41e-15     78.6     45.74
    1457     448   657       7   209  4.50e-15     77.8     49.77
    1458     447   657      11   220  4.53e-15     78.6     45.74
    1459     461   650      11   224  4.80e-15     78.2     46.58
    1460     463   701      14   275  5.02e-15     78.2     41.20
    1461     447   657      28   237  5.03e-15     78.2     45.74
    1462     459   656      36   236  5.19e-15     78.6     51.21
    1463     448   657       6   208  5.25e-15     77.8     49.77
    1464     463   704      46   281  5.46e-15     79.3     44.76
    1465     463   707      28   281  6.33e-15     78.6     46.04
    1466     455   674      29   248  9.63e-15     77.4     50.66
    1467     456   695       9   258  1.42e-14     77.0     46.88
    1468     463   656     193   389  1.52e-14     79.0     48.28
    1469     456   695      12   261  1.66e-14     77.0     46.88
    1470     455   672      63   286  1.71e-14     78.2     44.69
    1471     456   695      15   264  1.89e-14     76.6     46.88
    1472     463   656     193   389  2.05e-14     78.6     48.28
    1473     463   656     193   389  2.07e-14     78.6     48.28
    1474     463   656     193   389  2.09e-14     78.6     48.28
    1475     463   656     193   389  2.09e-14     78.6     48.28
    1476     455   674      29   248  2.09e-14     76.6     50.22
    1477     455   674      29   248  2.11e-14     76.6     50.22
    1478     463   656     193   389  2.22e-14     78.6     48.28
    1479     463   656     193   389  2.23e-14     78.6     48.28
    1480     157   227       8    81  2.53e-14     70.5     64.86
    1481     450   670      18   240  2.63e-14     77.0     49.79
    1482     455   658       9   210  2.67e-14     76.6     48.11
    1483     463   656     164   360  3.06e-14     77.8     48.28
    1484     450   670      14   236  3.07e-14     76.6     49.79
    1485     550   714     198   356  3.21e-14     76.6     51.18
    1486     450   670      25   247  3.23e-14     76.6     49.79
    1487     450   670      28   250  3.31e-14     76.6     49.79
    1488     450   670      25   247  3.57e-14     76.6     49.79
    1489     450   670      18   240  3.65e-14     76.6     49.79
    1490     450   670      25   247  3.67e-14     76.6     49.79
    1491     459   704      15   259  3.68e-14     75.5     48.83
    1492     450   670      19   241  3.74e-14     76.6     49.79
    1493     553   714     195   350  3.75e-14     76.6     50.30
    1494     553   714     202   357  3.83e-14     76.6     50.30
    1495     450   670      14   236  3.88e-14     76.3     49.79
    1496     450   670      18   240  3.89e-14     76.3     49.79
    1497     450   670      37   259  3.89e-14     76.6     49.79
    1498     450   670      18   240  3.96e-14     76.3     49.79
    1499     460   661      36   238  3.99e-14     75.9     47.34
    1500     450   670      38   260  3.99e-14     76.6     49.79
    1501     553   714     193   348  4.07e-14     76.3     50.30
    1502     450   670      25   247  4.16e-14     76.6     49.79
    1503     553   714     200   355  4.16e-14     76.3     50.30
    1504     450   670      23   245  4.17e-14     76.3     49.79
    1505     450   670      13   235  4.25e-14     76.3     49.79
    1506     450   670      20   242  4.49e-14     76.3     49.79
    1507     450   670      18   240  4.57e-14     76.3     49.79
    1508     450   670      18   240  4.57e-14     76.3     49.79
    1509     450   670      18   240  4.61e-14     76.3     49.79
    1510     450   670      18   240  4.62e-14     76.3     49.79
    1511     460   667      37   243  4.63e-14     75.5     47.89
    1512     450   670      20   242  4.91e-14     76.3     49.79
    1513     450   670      15   237  4.97e-14     75.9     49.79
    1514     450   670      18   240  5.01e-14     76.3     49.79
    1515     450   670      18   240  5.01e-14     75.9     49.79
    1516     450   670      39   261  5.13e-14     76.3     49.79
    1517     450   670      18   240  5.19e-14     76.3     49.79
    1518     450   670      18   240  5.19e-14     76.3     49.37
    1519     450   670      30   252  5.24e-14     76.3     49.79
    1520     450   670      14   236  5.29e-14     75.9     49.79
    1521     450   670      39   261  5.33e-14     76.3     49.79
    1522     450   670      37   259  5.34e-14     76.3     49.79
    1523     450   670      20   242  5.38e-14     75.9     49.79
    1524     450   670      18   240  5.38e-14     75.9     49.79
    1525     450   670      17   239  5.38e-14     75.9     49.79
    1526     450   670      30   252  5.43e-14     76.3     49.79
    1527     450   670      30   252  5.43e-14     76.3     49.79
    1528     450   670      16   238  5.44e-14     75.9     49.79
    1529     450   670      23   245  5.47e-14     76.3     49.79
    1530     450   670      18   240  5.48e-14     75.9     49.79
    1531     450   670      18   240  5.48e-14     75.9     49.79
    1532     450   670      18   240  5.53e-14     75.9     49.79
    1533     450   670      14   236  5.59e-14     75.9     49.79
    1534     459   656      10   207  5.63e-14     74.7     48.29
    1535     450   670      23   245  5.67e-14     75.9     49.79
    1536     450   670      20   242  5.68e-14     75.9     49.79
    1537     450   670      20   242  5.78e-14     75.9     49.79
    1538     450   670      18   240  5.84e-14     75.9     49.79
    1539     450   670      15   237  5.85e-14     75.9     49.79
    1540     450   670      25   247  5.87e-14     75.9     49.79
    1541     450   670      38   260  5.88e-14     76.3     49.79
    1542     450   670      29   251  5.95e-14     75.9     49.79
    1543     450   670      45   267  6.03e-14     76.3     49.79
    1544     450   670      21   243  6.10e-14     75.9     49.79
    1545     554   721     334   514  6.11e-14     77.4     49.46
    1546     450   670      38   260  6.21e-14     75.9     49.37
    1547     450   670      23   245  6.37e-14     75.9     49.79
    1548     450   670      22   244  6.44e-14     75.9     49.37
    1549     450   670      24   246  6.54e-14     75.9     49.79
    1550     459   656       9   206  6.59e-14     74.7     48.29
    1551     450   670      24   246  6.60e-14     75.9     49.79
    1552     450   670      29   251  6.63e-14     75.9     49.79
    1553     459   656      10   207  6.64e-14     74.7     48.29
    1554     459   656      33   230  6.80e-14     75.9     48.29
    1555     450   670      20   242  7.06e-14     75.9     49.37
    1556     450   670      18   240  7.19e-14     75.5     49.79
    1557     450   670      18   240  7.19e-14     75.5     49.37
    1558     450   670      41   263  7.20e-14     75.9     49.79
    1559     461   650      25   241  7.23e-14     75.1     45.95
    1560     450   670      23   245  7.33e-14     75.5     49.79
    1561     450   670      18   240  7.39e-14     75.5     49.79
    1562     460   667      15   221  7.44e-14     74.7     47.42
    1563     450   670      19   241  7.45e-14     75.5     49.79
    1564     450   670      41   263  7.53e-14     75.9     49.79
    1565     450   670      41   263  7.67e-14     75.9     49.79
    1566     450   670      18   240  7.73e-14     75.5     49.79
    1567     450   670      18   240  7.81e-14     75.5     49.79
    1568     450   670      18   240  8.02e-14     75.5     49.79
    1569     450   670      42   264  8.38e-14     75.9     49.79
    1570     450   670      41   263  8.39e-14     75.9     49.79
    1571     450   670      24   246  9.31e-14     75.5     49.79
    1572     460   667      14   220  9.44e-14     74.3     47.42
    1573     461   666       8   211  9.52e-14     74.3     48.57
    1574     460   667      38   244  9.87e-14     74.7     47.42
    1575     461   666       8   211  9.88e-14     74.3     48.57
    1576     450   670      20   242  9.95e-14     75.1     49.79
    1577     460   661      37   239  1.02e-13     74.7     47.34
    1578     450   670      18   240  1.05e-13     75.1     49.79
    1579     450   670      18   246  1.07e-13     75.1     48.52
    1580     450   660       3   207  1.08e-13     74.7     49.54
    1581     460   667      38   244  1.08e-13     74.3     47.42
    1582     450   670      23   245  1.09e-13     75.1     49.79
    1583     460   667      37   243  1.09e-13     74.3     47.42
    1584     460   661      37   239  1.09e-13     74.3     46.86
    1585     460   667      13   219  1.10e-13     73.9     47.42
    1586     456   722      19   294  1.10e-13     75.1     45.26
    1587     460   667      42   248  1.11e-13     74.3     47.42
    1588     460   667      12   218  1.12e-13     73.9     47.42
    1589     450   670      27   249  1.12e-13     75.1     49.79
    1590     450   670      18   240  1.14e-13     75.1     49.37
    1591     456   722      20   295  1.19e-13     74.7     45.26
    1592     460   667      40   246  1.22e-13     74.3     47.42
    1593     460   661      37   239  1.24e-13     74.3     46.86
    1594     450   670      24   246  1.25e-13     75.1     49.79
    1595     460   667      35   241  1.26e-13     74.3     47.42
    1596     460   667      35   241  1.27e-13     74.3     47.42
    1597     450   670      24   246  1.28e-13     75.1     49.79
    1598     456   722      33   308  1.29e-13     74.3     45.26
    1599     460   661      36   238  1.31e-13     74.3     46.86
    1600     450   670      41   263  1.32e-13     75.1     49.37
    1601     461   666       8   211  1.34e-13     73.9     48.57
    1602     466   649      35   231  1.35e-13     74.3     47.55
    1603     450   670      46   268  1.35e-13     75.1     49.79
    1604     460   667      37   243  1.40e-13     73.9     47.42
    1605     460   667      35   241  1.41e-13     73.9     47.42
    1606     450   670      46   268  1.46e-13     75.1     49.79
    1607     460   667      37   243  1.49e-13     73.9     47.42
    1608     456   722      17   292  1.49e-13     73.9     45.26
    1609     462   657      35   228  1.55e-13     73.6     50.50
    1610     450   670      18   240  1.57e-13     74.7     49.79
    1611     455   672       5   228  1.60e-13     73.9     44.25
    1612     420   669      20   279  1.64e-13     75.1     45.90
    1613     463   704      34   322  1.66e-13     74.3     45.24
    1614     450   670      39   261  1.67e-13     74.7     49.37
    1615     456   722      18   293  1.83e-13     73.9     45.26
    1616     455   656      19   228  2.22e-13     75.5     47.51
    1617     460   667      19   225  2.53e-13     73.2     47.42
    1618     420   669      19   278  2.56e-13     74.3     45.90
    1619     479   668      31   225  2.60e-13     72.8     51.96
    1620     452   651      12   225  2.77e-13     73.2     47.98
    1621     450   670      25   247  2.82e-13     73.9     48.95
    1622     522   701     781   969  3.18e-13     75.5     45.88
    1623     450   670      18   240  3.29e-13     73.6     49.37
    1624     479   668      44   238  3.83e-13     72.8     51.96
    1625     463   669      17   233  3.94e-13     72.8     50.42
    1626     456   722      33   308  4.10e-13     72.8     45.26
    1627     455   656      15   224  4.21e-13     74.7     47.51
    1628     455   656      13   222  4.31e-13     74.7     47.51
    1629     463   669      17   233  4.49e-13     72.4     50.42
    1630     460   667      34   240  4.53e-13     72.4     47.42
    1631     460   667      35   241  4.55e-13     72.4     46.95
    1632     520   650     792   932  4.63e-13     75.1     49.65
    1633     460   667      34   240  4.78e-13     72.4     47.42
    1634     463   665      12   220  5.19e-13     72.4     47.73
    1635     451   711      14   279  6.10e-13     73.2     44.59
    1636     463   699      27   256  6.66e-13     72.4     47.18
    1637     463   699      25   254  6.84e-13     72.4     47.18
    1638     450   665      16   230  6.84e-13     72.4     49.57
    1639     463   699      25   254  6.96e-13     72.4     47.18
    1640     455   656      19   228  7.18e-13     73.9     47.06
    1641     463   699      26   255  7.30e-13     72.4     45.97
    1642     463   699      26   255  7.36e-13     72.4     45.97
    1643     463   699      26   255  7.43e-13     72.4     45.97
    1644     463   699      16   245  7.44e-13     72.4     45.97
    1645     450   665      18   232  7.50e-13     72.4     49.57
    1646     462   657      35   228  7.51e-13     71.6     50.00
    1647     479   668      44   238  7.86e-13     71.6     51.47
    1648     445   669       1   228  8.37e-13     72.0     45.87
    1649     461   665      31   235  8.62e-13     72.4     49.10
    1650     463   699      26   255  9.23e-13     72.0     47.18
    1651     445   669       2   229  9.40e-13     72.0     45.87
    1652     463   665      14   222  9.89e-13     71.2     47.27
    1653     463   665      12   220  1.01e-12     71.6     47.27
    1654     445   669       3   230  1.03e-12     72.0     45.87
    1655     445   669       2   229  1.03e-12     72.0     45.87
    1656     445   669       4   231  1.05e-12     72.0     45.87
    1657     525   649     799   937  1.07e-12     73.9     55.32
    1658     493   612     339   462  9.30e-07     54.7     50.78
    1659     463   665      26   231  1.08e-12     71.2     48.13
    1660     445   669       6   233  1.08e-12     72.0     45.87
    1661     463   705      81   314  1.17e-12     71.6     46.77
    1662     462   668      10   215  1.23e-12     71.2     48.82
    1663     463   665      27   229  1.27e-12     72.0     49.55
    1664     444   658      12   254  1.35e-12     71.6     47.97
    1665     456   658      37   242  1.43e-12     71.6     48.58
    1666     463   665      12   220  1.47e-12     71.2     47.27
    1667     461   665      26   230  1.48e-12     71.6     49.55
    1668     459   708      14   272  1.54e-12     71.6     44.65
    1669     463   665      13   218  1.56e-12     70.5     48.13
    1670     450   665      24   238  1.58e-12     71.6     49.57
    1671     456   665      17   238  1.88e-12     70.9     45.13
    1672     445   669       3   230  1.89e-12     71.2     46.69
    1673     463   665      12   220  1.93e-12     70.5     47.27
    1674     456   665      10   231  2.04e-12     70.9     45.13
    1675     548   721     333   515  2.06e-12     72.4     47.12
    1676     561   711     125   273  2.08e-12     70.9     49.35
    1677     456   665      11   232  2.17e-12     70.9     45.13
    1678     456   665      11   232  2.19e-12     70.9     45.13
    1679     443   658       9   214  2.25e-12     70.1     45.70
    1680     561   711     126   274  2.27e-12     70.5     49.35
    1681     443   658       9   214  2.31e-12     70.1     45.70
    1682     456   665      11   232  2.34e-12     70.5     45.13
    1683     443   658      10   215  2.34e-12     70.1     45.70
    1684     443   658      10   215  2.36e-12     70.1     45.70
    1685     463   663      31   231  2.44e-12     70.9     45.89
    1686     443   658      11   216  2.45e-12     70.1     45.70
    1687     463   663      31   231  2.59e-12     70.9     45.89
    1688     450   665      20   233  2.69e-12     70.9     49.15
    1689     450   665      24   238  2.73e-12     70.9     49.57
    1690     443   658      10   215  3.01e-12     69.7     45.70
    1691     561   711     137   285  3.06e-12     70.1     49.35
    1692     561   711     134   282  3.09e-12     70.1     49.35
    1693     561   711     125   273  3.10e-12     70.1     49.35
    1694     561   711     137   285  3.15e-12     70.1     49.35
    1695     441   711       1   280  3.15e-12     70.9     43.87
    1696     463   661      32   235  3.18e-12     70.5     41.63
    1697     461   665      27   229  3.27e-12     70.1     49.76
    1698     448   669       1   225  3.29e-12     70.5     46.03
    1699     561   711     129   277  3.32e-12     70.1     49.35
    1700     561   711     129   277  3.35e-12     70.1     49.35
    1701     459   669       4   224  3.35e-12     70.5     46.49
    1702     561   711     126   274  3.40e-12     70.1     49.35
    1703     450   711      17   284  3.46e-12     70.1     44.44
    1704     561   711     126   274  3.49e-12     70.1     49.35
    1705     427   656       8   239  3.50e-12     70.5     44.86
    1706     427   656       7   238  3.51e-12     70.5     44.86
    1707     456   658      36   244  3.66e-12     70.1     49.07
    1708     561   711     170   318  3.92e-12     70.1     49.35
    1709     450   665      38   251  4.12e-12     70.5     49.15
    1710     457   651      10   225  4.15e-12     69.7     47.09
    1711     451   711      16   281  4.17e-12     70.5     44.26
    1712     463   665      17   228  4.25e-12     69.7     46.19
    1713     427   656       7   238  4.28e-12     70.1     44.86
    1714     461   665      27   229  4.40e-12     69.7     49.76
    1715     450   665      19   232  4.50e-12     70.1     49.15
    1716     451   711      16   281  4.61e-12     70.5     45.52
    1717     457   651       9   224  4.69e-12     69.7     47.09
    1718     463   665      12   223  4.76e-12     69.3     46.19
    1719     451   658      12   243  4.77e-12     69.7     48.51
    1720     463   665      19   230  4.85e-12     69.3     46.19
    1721     451   711      14   279  4.86e-12     70.1     44.26
    1722     561   711     129   277  5.15e-12     69.3     49.35
    1723     451   711      13   278  5.20e-12     70.1     44.26
    1724     451   711      14   279  5.25e-12     70.1     44.26
    1725     427   656       7   238  5.31e-12     69.7     44.44
    1726     427   656       7   238  5.41e-12     69.7     44.44
    1727     457   651      11   226  5.53e-12     69.3     47.09
    1728     427   656       8   239  5.81e-12     69.7     44.44
    1729     427   656      11   242  5.95e-12     69.7     44.44
    1730     427   656       7   238  6.14e-12     69.7     44.44
    1731     455   656      14   223  6.18e-12     70.9     47.96
    1732     500   662     234   392  6.38e-12     70.9     47.88
    1733     451   711      14   279  6.58e-12     69.7     43.48
    1734     427   656       7   238  6.60e-12     69.7     44.44
    1735     455   656      15   224  6.69e-12     70.9     47.96
    1736     427   656      13   244  6.90e-12     69.7     44.44
    1737     463   650      43   232  6.93e-12     69.7     47.78
    1738     427   656      28   259  6.93e-12     69.7     44.86
    1739     427   656      14   245  6.95e-12     69.7     44.44
    1740     450   665      39   252  7.04e-12     69.7     48.73
    1741     463   677      10   223  7.27e-12     68.9     47.73
    1742     451   711      14   279  8.23e-12     69.3     43.81
    1743     495   710      72   282  8.37e-12     68.6     45.45
    1744     427   656       9   240  8.39e-12     69.3     44.90
    1745     427   656       7   238  9.13e-12     68.9     43.85
    1746     388   665      64   357  9.74e-12     69.7     42.67
    1747     456   731      35   309  9.84e-12     68.6     42.28
    1748     463   657      62   263  9.90e-12     68.9     48.56
    1749     427   656       7   238  9.99e-12     68.9     44.44
    1750     442   656       9   225  1.01e-11     68.9     46.05
    1751     519   705      90   269  1.07e-11     68.9     48.95
    1752     554   650     167   270  1.08e-11     68.6     56.88
    1753     519   663      79   217  1.10e-11     68.9     52.41
    1754     519   663      82   220  1.12e-11     68.9     52.41
    1755     388   665      78   371  1.13e-11     69.3     42.67
    1756     427   656       8   239  1.17e-11     68.9     44.44
    1757     519   663      80   218  1.21e-11     68.9     52.41
    1758     427   656       7   238  1.22e-11     68.6     43.85
    1759     519   663     122   260  1.24e-11     68.9     52.41
    1760     461   665      11   218  1.33e-11     68.2     52.11
    1761     461   677      16   231  1.33e-11     68.2     48.46
    1762     519   663      76   214  1.34e-11     68.6     52.41
    1763     427   656       7   238  1.36e-11     68.6     45.90
    1764     427   656       7   238  1.36e-11     68.6     45.90
    1765     519   663      75   213  1.37e-11     68.6     52.41
    1766     427   656       8   239  1.38e-11     68.6     44.49
    1767     427   656       7   238  1.39e-11     68.6     45.90
    1768     519   663      75   213  1.41e-11     68.6     52.41
    1769     449   651      48   279  1.45e-11     69.3     46.31
    1770     463   650      44   233  1.47e-11     68.6     47.78
    1771     430   656      10   238  1.47e-11     68.6     44.58
    1772     519   663     322   460  1.48e-11     69.7     52.41
    1773     235   280     141   190  4.00e+00     32.7     44.00
    1774     461   656      37   233  1.51e-11     68.2     45.41
    1775     430   656      10   238  1.54e-11     68.6     44.58
    1776     500   656     233   385  1.55e-11     69.7     48.43
    1777     519   663      75   213  1.55e-11     68.2     52.41
    1778     429   656       2   231  1.55e-11     68.2     45.23
    1779     427   656       7   238  1.57e-11     68.6     44.44
    1780     501   714      57   272  1.58e-11     68.6     47.77
    1781     429   656       2   231  1.60e-11     68.2     45.23
    1782     501   714      57   272  1.68e-11     68.6     47.77
    1783     561   711     137   285  1.79e-11     67.8     48.70
    1784     463   650      21   210  1.83e-11     68.2     48.28
    1785     463   666      33   243  1.84e-11     67.8     46.95
    1786     451   658       5   214  1.87e-11     67.4     44.95
    1787     461   677      13   228  1.88e-11     67.4     47.75
    1788     427   656      28   259  1.89e-11     68.2     45.27
    1789     463   650      51   240  1.96e-11     68.2     47.78
    1790     442   656      22   238  2.05e-11     68.2     46.05
    1791     479   668      45   239  2.05e-11     68.2     52.94
    1792     442   656      22   238  2.15e-11     67.8     46.05
    1793     442   656      22   238  2.21e-11     67.8     46.05
    1794     451   658       2   211  2.21e-11     67.0     44.95
    1795     451   692      14   252  2.22e-11     68.2     45.39
    1796     442   656      22   238  2.25e-11     67.8     46.05
    1797     442   656      65   281  2.28e-11     68.2     46.49
    1798     463   650      21   210  2.33e-11     67.8     47.78
    1799     451   658       6   215  2.49e-11     67.0     44.95
    1800     461   656      10   206  2.57e-11     67.4     45.89
    1801     442   656      11   227  2.66e-11     67.4     46.05
    1802     461   658      11   209  2.69e-11     66.6     45.75
    1803     442   656      11   227  2.71e-11     67.4     46.05
    1804     442   656      17   233  2.76e-11     67.4     46.05
    1805     442   656      10   226  2.79e-11     67.4     46.05
    1806     442   656      11   227  2.81e-11     67.4     46.05
    1807     463   650      18   207  2.83e-11     67.4     47.78
    1808     463   650      18   207  2.91e-11     67.4     47.78
    1809     461   656      19   215  2.91e-11     67.0     46.19
    1810     427   656       7   238  2.91e-11     67.4     44.44
    1811     427   656       7   238  2.91e-11     67.4     44.31
    1812     427   656       8   239  2.94e-11     67.4     44.03
    1813     427   656       7   238  2.97e-11     67.4     44.44
    1814     451   692      52   290  3.19e-11     68.2     45.39
    1815     427   656      28   259  3.20e-11     67.8     44.86
    1816     427   656      11   242  3.34e-11     67.4     44.44
    1817     463   650      11   200  3.41e-11     66.6     47.78
    1818     427   656       8   239  3.42e-11     67.4     44.44
    1819     427   656       7   238  3.58e-11     67.4     44.03
    1820     461   677      16   231  3.61e-11     66.6     47.75
    1821     427   656       7   238  3.62e-11     67.4     43.85
    1822     451   692      14   252  3.79e-11     67.4     45.39
    1823     451   692      14   252  3.80e-11     67.4     45.39
    1824     461   656      19   215  3.82e-11     67.0     45.41
    1825     442   656      14   230  3.87e-11     67.0     46.05
    1826     461   656       9   205  3.90e-11     66.6     45.89
    1827     461   656      10   206  3.94e-11     66.6     45.89
    1828     427   656       7   238  3.96e-11     67.0     43.85
    1829     461   676       8   222  4.00e-11     66.6     46.82
    1830     461   671      12   221  4.06e-11     66.6     46.98
    1831     427   656       7   238  4.10e-11     67.0     43.44
    1832     451   692       7   245  4.20e-11     67.0     45.39
    1833     461   677      12   227  4.20e-11     66.6     48.66
    1834     427   656       8   239  4.21e-11     67.0     44.26
    1835     451   664      14   233  4.23e-11     67.4     45.83
    1836     461   656      28   224  4.26e-11     66.6     45.41
    1837     461   671       8   217  4.43e-11     66.2     46.98
    1838     457   664     190   398  4.49e-11     68.2     45.41
    1839     463   650      63   252  4.51e-11     67.4     47.78
    1840     451   664      14   233  4.59e-11     67.0     45.83
    1841     461   677      13   228  4.65e-11     66.2     47.75
    1842     427   656       8   239  4.68e-11     67.0     44.08
    1843     451   711      19   284  4.81e-11     67.0     43.48
    1844     451   692      14   252  4.85e-11     67.0     45.39
    1845     461   677      16   231  4.92e-11     66.2     47.75
    1846     461   677       9   224  4.98e-11     66.2     47.75
    1847     451   692      14   252  5.03e-11     67.0     45.39
    1848     451   692      15   253  5.04e-11     67.0     45.39
    1849     451   692      15   253  5.06e-11     67.0     45.39
    1850     451   711       8   273  5.09e-11     67.0     43.48
    1851     427   656       8   240  5.11e-11     66.6     43.90
    1852     451   692      13   251  5.13e-11     67.0     45.39
    1853     461   677       8   223  5.16e-11     66.2     48.05
    1854     451   692      52   290  5.21e-11     67.4     45.39
    1855     451   692      17   255  5.29e-11     67.0     45.39
    1856     461   677      11   226  5.37e-11     66.2     47.75
    1857     451   711      14   279  5.38e-11     67.0     43.48
    1858     461   677      12   227  5.42e-11     66.2     47.75
    1859     451   692      14   252  5.43e-11     66.6     45.39
    1860     451   692      16   254  5.45e-11     67.0     45.39
    1861     451   692      14   252  5.48e-11     66.6     45.39
    1862     451   692      15   253  5.51e-11     66.6     45.39
    1863     451   711      15   280  5.53e-11     66.6     43.48
    1864     451   711      10   275  5.59e-11     66.6     43.48
    1865     461   677      11   226  5.62e-11     66.2     47.75
    1866     461   677      11   226  5.73e-11     66.2     47.75
    1867     461   677      12   227  5.73e-11     66.2     47.75
    1868     461   677      13   228  5.75e-11     66.2     47.75
    1869     451   692       7   245  5.76e-11     66.6     45.39
    1870     461   677      12   227  5.78e-11     66.2     48.21
    1871     461   677      16   231  5.79e-11     66.2     47.75
    1872     461   677       9   224  5.82e-11     65.9     47.75
    1873     455   658      12   220  5.91e-11     66.6     44.89
    1874     500   656     233   385  6.01e-11     67.8     47.80
    1875     461   677      13   228  6.06e-11     65.9     47.75
    1876     427   656       8   239  6.08e-11     66.6     44.26
    1877     451   664      14   233  6.10e-11     66.6     45.42
    1878     451   692       8   246  6.17e-11     66.6     45.39
    1879     461   656      12   208  6.21e-11     67.4     45.71
    1880     461   677      12   227  6.28e-11     65.9     47.83
    1881     461   677      10   225  6.50e-11     65.9     47.75
    1882     422   656      57   293  6.52e-11     67.0     44.76
    1883     451   711      24   289  6.55e-11     66.6     43.48
    1884     444   658      14   256  6.63e-11     66.2     46.34
    1885     461   677      12   227  6.63e-11     65.9     48.21
    1886     463   658      20   220  6.64e-11     66.2     44.70
    1887     461   677       9   224  6.73e-11     65.9     47.75
    1888     461   677      10   225  6.81e-11     65.9     48.21
    1889     461   677      11   226  6.87e-11     65.9     47.75
    1890     461   677       8   223  6.92e-11     65.9     47.75
    1891     461   677      13   228  7.07e-11     65.9     47.75
    1892     461   677       9   224  7.11e-11     65.9     47.75
    1893     461   677      11   226  7.19e-11     65.9     48.21
    1894     461   656      19   215  7.21e-11     67.0     46.19
    1895     461   677       7   222  7.24e-11     65.9     47.75
    1896     461   677      10   225  7.25e-11     65.9     47.75
    1897     422   656      57   293  7.31e-11     66.6     44.35
    1898     422   656      57   293  7.31e-11     66.6     44.35
    1899     461   677      10   225  7.39e-11     65.9     47.75
    1900     461   677      12   227  7.40e-11     65.9     47.30
    1901     455   658      11   219  7.43e-11     65.5     44.89
    1902     422   656      57   293  7.44e-11     66.6     44.35
    1903     461   677      10   225  7.46e-11     65.9     47.75
    1904     500   656     233   385  7.64e-11     67.4     47.80
    1905     500   656     234   386  7.65e-11     67.4     47.80
    1906     422   656      57   293  7.71e-11     66.6     44.35
    1907     461   677       8   223  7.79e-11     65.5     47.75
    1908     461   677       9   224  7.80e-11     65.5     47.83
    1909     461   677       9   224  7.80e-11     65.5     47.75
    1910     461   677       9   224  7.93e-11     65.5     47.75
    1911     427   656       8   239  7.96e-11     66.2     44.86
    1912     461   677       8   223  8.00e-11     65.5     47.75
    1913     461   677       9   224  8.02e-11     65.5     47.75
    1914     461   677      11   226  8.03e-11     65.5     47.75
    1915     451   711      10   275  8.05e-11     66.2     43.48
    1916     461   677       8   223  8.08e-11     65.5     47.75
    1917     461   677       9   224  8.15e-11     65.5     47.75
    1918     461   677       9   224  8.24e-11     65.5     47.75
    1919     427   656       8   239  8.25e-11     66.2     43.85
    1920     463   650      21   210  8.28e-11     66.2     47.29
    1921     461   656       6   202  8.29e-11     67.0     45.71
    1922     461   677       9   224  8.31e-11     65.5     47.75
    1923     461   677       9   224  8.31e-11     65.5     47.75
    1924     427   656       7   238  8.49e-11     66.2     43.44
    1925     427   656       8   239  8.55e-11     66.2     44.86
    1926     427   656       8   239  8.55e-11     66.2     43.85
    1927     461   677      16   231  8.65e-11     65.5     47.75
    1928     461   677       8   223  8.70e-11     65.5     47.75
    1929     461   677       8   223  8.70e-11     65.5     47.75
    1930     461   677       9   224  8.78e-11     65.5     47.75
    1931     461   677       9   224  8.86e-11     65.5     47.75
    1932     427   656       7   238  8.96e-11     65.9     44.31
    1933     553   710     245   396  9.11e-11     66.6     50.30
    1934     450   548      60   167  1.50e-02     40.4     48.62
    1935     427   656      10   241  9.22e-11     65.9     44.86
    1936     461   677      11   226  9.29e-11     65.5     47.75
    1937     454   665      10   228  9.31e-11     65.5     45.33
    1938     463   650      18   207  9.62e-11     65.9     47.29
    1939     427   656      10   241  9.65e-11     65.9     44.86
    1940     427   656       7   238  9.71e-11     65.9     43.85
    1941     427   656       8   239  9.78e-11     65.9     43.85
    1942     461   677      10   225  9.79e-11     65.5     47.75
    1943     463   649      21   228  1.01e-10     65.9     43.58
    1944     547   658     113   220  1.02e-10     65.9     52.68
    1945     427   656       7   238  1.06e-10     65.9     43.85
    1946     457   664     200   408  1.07e-10     66.6     45.87
    1947     427   656       7   238  1.07e-10     65.9     43.85
    1948     427   656       7   238  1.07e-10     65.9     43.85
    1949     427   656       8   239  1.10e-10     65.9     44.03
    1950     427   656       8   239  1.12e-10     65.9     43.85
    1951     461   658      11   209  1.13e-10     64.7     45.75
    1952     456   659       8   216  1.13e-10     66.6     46.30
    1953     451   664      14   233  1.14e-10     65.9     45.42
    1954     427   656       8   239  1.20e-10     65.5     43.85
    1955     451   664      14   233  1.20e-10     65.9     45.42
    1956     463   649      43   250  1.28e-10     65.5     43.58
    1957     427   656      10   241  1.33e-10     65.5     43.85
    1958     463   650      18   207  1.35e-10     65.5     47.29
    1959     461   677      11   226  1.36e-10     65.1     47.96
    1960     463   650      63   252  1.45e-10     65.9     47.29
    1961     454   665      10   228  1.46e-10     65.1     44.89
    1962     451   664      14   233  1.54e-10     65.5     45.42
    1963     451   664      14   233  1.59e-10     65.5     45.42
    1964     427   656       8   239  1.60e-10     65.1     43.85
    1965     457   664     190   398  1.62e-10     66.2     44.95
    1966     442   656       9   225  1.65e-10     65.1     46.05
    1967     454   665      10   228  1.68e-10     64.7     44.89
    1968     457   664     191   399  1.70e-10     66.2     44.95
    1969     442   653       3   225  1.72e-10     65.1     47.58
    1970     457   664     191   399  1.72e-10     66.2     44.95
    1971     450   656      31   239  1.97e-10     65.1     47.27
    1972     457   664     191   399  1.97e-10     66.2     44.95
    1973     461   677       8   223  1.98e-10     64.3     47.30
    1974     441   694     525   775  2.10e-10     66.2     45.05
    1975     450   656      33   241  2.12e-10     65.1     47.27
    1976     457   664     191   399  2.13e-10     65.9     44.95
    1977     457   664     191   399  2.14e-10     65.9     44.95
    1978     438   660     258   480  2.18e-10     66.2     42.34
    1979     450   656      33   241  2.18e-10     64.7     47.27
    1980     455   658      32   240  2.24e-10     64.3     44.89
    1981     461   676       8   222  2.26e-10     64.3     46.82
    1982     454   665      10   228  2.27e-10     64.3     44.89
    1983     457   664     162   370  2.27e-10     65.9     44.95
    1984     457   664     163   371  2.27e-10     65.9     44.95
    1985     454   665      10   228  2.40e-10     64.3     44.89
    1986     457   664     162   370  2.43e-10     65.9     44.95
    1987     462   665      90   294  2.55e-10     65.1     47.87
    1988     461   677      16   231  2.56e-10     64.3     47.30
    1989     442   656       8   224  2.64e-10     64.3     45.85
    1990     457   664     162   370  2.68e-10     65.5     44.95
    1991     457   664     162   370  2.79e-10     65.5     44.95
    1992     451   670       8   238  2.83e-10     64.7     45.38
    1993     442   656      23   239  2.84e-10     64.3     45.41
    1994     427   656       7   238  2.98e-10     64.3     43.62
    1995     427   656       7   238  3.01e-10     64.3     43.62
    1996     463   653      19   219  3.20e-10     64.3     48.53
    1997     461   677       8   223  3.24e-10     63.9     47.30
    1998     562   710     148   290  3.36e-10     63.5     49.67
    1999     562   710     148   290  3.55e-10     63.5     49.67
    2000     427   656       7   238  3.66e-10     64.3     43.44
    2001     451   670      30   260  3.85e-10     64.3     44.96
    2002     461   677      10   225  3.96e-10     63.5     47.30
    2003     427   656       7   238  4.00e-10     63.9     43.44
    2004     427   656       8   239  4.03e-10     63.9     43.44
    2005     461   677       8   223  4.14e-10     63.5     47.30
    2006     461   677       9   224  4.18e-10     63.5     47.30
    2007     556   711     158   309  5.00e-10     63.5     51.90
    2008     501   717      74   287  5.42e-10     62.8     40.91
    2009     397   711      10   313  5.83e-10     63.5     43.34
    2010     459   675      25   240  5.86e-10     63.5     48.31
    2011     562   710     140   282  5.90e-10     62.8     49.67
    2012     442   656      23   239  6.03e-10     63.5     45.41
    2013     524   656      85   226  6.29e-10     63.2     51.75
    2014     451   670       2   232  6.52e-10     63.2     44.96
    2015     459   675      26   241  6.57e-10     63.2     48.31
    2016     463   665       9   218  6.64e-10     62.8     45.83
    2017     459   675      39   254  6.75e-10     63.5     48.31
    2018     444   672      17   243  6.86e-10     63.2     47.93
    2019     462   651      25   225  7.17e-10     63.2     45.10
    2020     454   672       6   222  7.18e-10     63.2     48.71
    2021     454   672       8   224  7.27e-10     63.2     47.70
    2022     459   669      30   246  7.64e-10     63.2     45.49
    2023     457   664     191   399  7.66e-10     64.3     44.95
    2024     462   651      26   226  7.69e-10     63.2     45.10
    2025     451   670       8   238  8.11e-10     63.2     44.54
    2026     451   670      10   240  8.50e-10     63.2     44.96
    2027     459   675      48   263  8.56e-10     63.2     48.31
    2028     463   656      19   209  9.08e-10     62.8     48.79
    2029     451   670       6   236  9.17e-10     62.8     44.54
    2030     451   670       9   239  9.25e-10     62.8     44.54
    2031     461   665      19   224  9.50e-10     62.4     47.64
    2032     461   677       8   223  9.66e-10     62.4     46.85
    2033     451   664      14   233  9.78e-10     62.8     45.92
    2034     462   651      25   225  9.80e-10     63.2     45.10
    2035     451   664       5   224  9.90e-10     62.8     45.92
    2036     459   669       7   223  1.02e-09     62.8     45.49
    2037     454   672      27   243  1.04e-09     62.4     47.70
    2038     451   670       3   233  1.04e-09     62.8     44.54
    2039     451   664      12   231  1.07e-09     62.8     45.92
    2040     451   670      16   246  1.08e-09     62.8     44.54
    2041     462   651      25   225  1.08e-09     62.4     45.10
    2042     391   655      17   263  1.09e-09     63.2     45.72
    2043     451   670       4   234  1.10e-09     62.8     44.54
    2044     391   655      16   262  1.10e-09     62.8     45.72
    2045     563   665     645   738  1.12e-09     63.9     56.73
    2046     459   675      26   241  1.13e-09     62.4     48.91
    2047     451   670      16   246  1.14e-09     62.8     44.54
    2048     391   655      18   264  1.14e-09     62.8     45.72
    2049     451   670      17   247  1.14e-09     62.8     44.54
    2050     451   670      21   251  1.15e-09     62.8     44.54
    2051     563   665     644   737  1.15e-09     63.9     56.73
    2052     451   670       8   238  1.17e-09     62.8     44.54
    2053     451   664      31   250  1.17e-09     62.8     45.92
    2054     451   670      20   250  1.18e-09     62.8     44.54
    2055     451   670       1   231  1.21e-09     62.4     44.54
    2056     459   669      30   246  1.21e-09     62.8     45.06
    2057     451   670      20   250  1.23e-09     62.8     44.54
    2058     451   670      32   262  1.23e-09     62.8     45.19
    2059     462   651      24   224  1.28e-09     62.4     45.10
    2060     451   670      20   250  1.30e-09     62.4     44.54
    2061     451   670      16   246  1.32e-09     62.4     44.54
    2062     461   656      12   208  1.37e-09     63.2     45.24
    2063     451   670      12   242  1.38e-09     62.4     44.54
    2064     501   712      57   270  1.39e-09     62.4     45.26
    2065     451   670       4   234  1.45e-09     62.4     44.54
    2066     451   670       3   233  1.46e-09     62.4     44.54
    2067     451   670      13   243  1.46e-09     62.4     44.54
    2068     451   670       3   233  1.49e-09     62.4     44.54
    2069     556   711     114   265  1.51e-09     61.6     51.90
    2070     451   670      17   247  1.52e-09     62.4     44.96
    2071     451   670      16   246  1.53e-09     62.4     44.54
    2072     451   670      14   244  1.56e-09     62.4     44.54
    2073     451   670      21   251  1.58e-09     62.4     44.54
    2074     451   670      10   240  1.60e-09     62.4     44.54
    2075     451   670      12   242  1.60e-09     62.4     44.54
    2076     451   664      34   253  1.61e-09     62.4     43.70
    2077     391   655      16   262  1.61e-09     62.4     43.87
    2078     462   651      24   224  1.63e-09     62.0     45.10
    2079     461   670     236   448  1.63e-09     62.8     46.72
    2080     462   651      24   224  1.64e-09     62.0     45.10
    2081     451   670      32   262  1.65e-09     62.4     44.54
    2082     451   670       9   239  1.65e-09     62.0     44.54
    2083     462   651      24   224  1.66e-09     62.0     45.10
    2084     462   651      18   218  1.69e-09     62.0     45.10
    2085     451   670      10   240  1.69e-09     62.0     44.54
    2086     391   655      17   263  1.69e-09     62.4     43.87
    2087     451   670      20   250  1.70e-09     62.4     44.54
    2088     451   670      12   242  1.70e-09     62.4     44.54
    2089     451   670      30   254  1.71e-09     62.4     45.83
    2090     501   712      57   270  1.71e-09     62.0     45.26
    2091     451   670      16   246  1.71e-09     62.0     44.54
    2092     451   670      10   240  1.72e-09     62.0     44.54
    2093     451   670      17   247  1.72e-09     62.0     44.54
    2094     451   670      17   247  1.72e-09     62.0     44.54
    2095     451   670      18   248  1.73e-09     62.0     44.54
    2096     451   670      19   249  1.74e-09     62.0     44.54
    2097     462   651      26   226  1.74e-09     62.0     45.10
    2098     451   670       2   232  1.75e-09     62.0     44.54
    2099     451   670      10   240  1.75e-09     62.0     44.54
    2100     461   677      13   229  1.77e-09     61.6     44.77
    2101     451   670       3   233  1.78e-09     62.0     44.54
    2102     391   655      18   264  1.78e-09     62.4     43.87
    2103     451   655      81   281  1.78e-09     62.4     51.43
    2104     501   712      57   270  1.79e-09     62.0     45.26
    2105     391   655      20   266  1.81e-09     62.4     43.87
    2106     451   670      12   242  1.82e-09     62.0     44.54
    2107     451   670      32   262  1.82e-09     62.0     44.54
    2108     391   655      22   268  1.83e-09     62.4     43.87
    2109     451   664      14   233  1.84e-09     62.4     45.06
    2110     451   670      13   243  1.84e-09     62.0     44.54
    2111     451   670      10   240  1.85e-09     62.0     44.54
    2112     461   656      17   213  1.90e-09     62.4     45.24
    2113     451   670      16   246  1.91e-09     62.0     44.54
    2114     461   677      13   229  1.92e-09     61.6     44.77
    2115     451   670       2   232  1.93e-09     62.0     44.54
    2116     501   712      57   270  1.94e-09     62.0     45.26
    2117     461   677       8   224  1.94e-09     61.6     44.77
    2118     391   655      38   284  2.00e-09     62.4     43.87
    2119     451   670      20   250  2.03e-09     62.0     44.54
    2120     391   655      38   284  2.05e-09     62.4     43.87
    2121     451   670      17   247  2.08e-09     62.0     44.54
    2122     451   670      15   245  2.21e-09     62.0     44.54
    2123     451   761       5   295  2.23e-09     61.2     42.68
    2124     461   677       8   224  2.25e-09     61.2     44.77
    2125     451   670      16   246  2.28e-09     61.6     44.54
    2126     451   670       2   232  2.40e-09     61.6     44.54
    2127     451   670       3   233  2.47e-09     61.6     44.54
    2128     451   670      20   250  2.47e-09     61.6     44.54
    2129     462   651      24   224  2.53e-09     61.2     45.59
    2130     463   653      21   221  2.54e-09     61.6     48.53
    2131     463   670      39   250  2.56e-09     61.6     44.75
    2132     451   670      12   242  2.62e-09     61.6     44.54
    2133     451   670      16   246  2.62e-09     61.6     44.54
    2134     451   670      16   246  2.65e-09     61.6     44.54
    2135     451   670      12   242  2.71e-09     61.6     44.54
    2136     461   670      49   262  2.83e-09     61.6     43.58
    2137     451   670      31   261  2.92e-09     61.6     44.54
    2138     451   670      16   246  3.05e-09     61.2     44.54
    2139     451   670      20   250  3.09e-09     61.2     44.54
    2140     463   721      19   290  3.21e-09     61.2     44.48
    2141     451   670      17   247  3.24e-09     61.2     44.54
    2142     451   670      31   261  3.31e-09     61.2     44.54
    2143     462   651      24   224  3.34e-09     60.8     45.10
    2144     451   670      12   242  3.45e-09     61.2     44.54
    2145     459   675      39   254  3.58e-09     61.2     47.88
    2146     451   670      10   240  3.60e-09     61.2     44.54
    2147     451   670      14   244  3.67e-09     61.2     44.12
    2148     461   670      46   259  4.16e-09     61.2     43.58
    2149     461   670      47   260  4.25e-09     60.8     43.58
    2150     451   670      11   241  4.73e-09     60.8     44.77
    2151     463   653      21   221  4.78e-09     60.5     48.53
    2152     451   655      64   264  4.98e-09     60.8     48.80
    2153     463   653      21   221  5.04e-09     60.5     48.53
    2154     463   653      21   221  5.10e-09     60.5     48.53
    2155     524   656      85   226  5.20e-09     60.5     51.05
    2156     524   656      83   224  5.37e-09     60.1     51.05
    2157     463   653      17   217  5.40e-09     60.1     48.53
    2158     451   761      21   315  6.22e-09     60.1     42.81
    2159     463   721      19   290  6.37e-09     60.5     44.48
    2160     451   670      15   245  6.94e-09     60.1     44.12
    2161     456   665       4   212  7.10e-09     59.7     45.50
    2162     463   721      19   290  7.54e-09     60.1     44.48
    2163     451   670      14   244  7.75e-09     60.1     44.12
    2164     456   665       6   214  8.01e-09     59.7     45.50
    2165     466   656      23   214  9.14e-09     59.7     44.22
    2166     451   655       2   202  1.01e-08     59.3     48.80
    2167     456   705       8   264  1.21e-08     59.7     43.02
    2168     452   651      24   247  1.25e-08     59.7     43.64
    2169     451   670      15   245  1.26e-08     59.3     44.12
    2170     451   655       2   202  1.27e-08     59.3     48.80
    2171     512   656      75   229  1.27e-08     58.9     49.36
    2172     463   653      32   224  1.33e-08     59.3     45.89
    2173     451   655       2   202  1.35e-08     59.3     48.80
    2174     451   655       2   202  1.37e-08     59.3     48.80
    2175     451   655      28   228  1.42e-08     59.3     48.80
    2176     512   656      72   226  1.44e-08     58.9     49.36
    2177     451   655      29   229  1.44e-08     59.3     48.80
    2178     157   224       6    76  1.49e-08     53.9     63.38
    2179     451   655      53   253  1.50e-08     59.3     48.80
    2180     451   655      30   230  1.58e-08     59.3     48.80
    2181     451   655      27   227  1.65e-08     58.9     48.80
    2182     451   655      27   227  1.67e-08     58.9     48.80
    2183     463   701      82   327  1.67e-08     59.3     45.24
    2184     451   655       2   202  1.70e-08     58.9     48.80
    2185     451   655      27   227  1.71e-08     58.9     48.80
    2186     451   655       2   202  1.76e-08     58.9     48.80
    2187     451   655      24   224  1.81e-08     58.9     48.80
    2188     451   655      24   224  1.86e-08     58.9     48.80
    2189     459   659      11   216  2.19e-08     58.9     44.91
    2190     451   655      24   224  2.20e-08     58.5     46.41
    2191     460   650      95   290  2.31e-08     58.9     47.52
    2192     463   657      42   235  2.39e-08     58.5     44.17
    2193     451   655      12   212  2.45e-08     58.5     48.80
    2194     461   646      17   221  2.51e-08     57.8     46.73
    2195     461   646      13   217  2.81e-08     57.8     46.73
    2196     461   646      14   218  2.83e-08     57.8     46.73
    2197     461   646      15   219  2.91e-08     57.8     46.73
    2198     463   657     312   505  3.07e-08     58.9     43.69
    2199     461   646      17   221  3.09e-08     57.8     46.73
    2200     461   646      15   219  3.10e-08     57.8     46.73
    2201     457   660      12   222  3.23e-08     57.4     43.12
    2202     461   646      13   217  3.37e-08     57.4     46.73
    2203     463   741      20   301  3.44e-08     57.8     46.23
    2204     559   721     128   291  3.96e-08     58.2     52.87
    2205     463   657     164   357  3.98e-08     58.2     43.69
    2206     463   650      38   227  4.03e-08     57.8     45.10
    2207     556   710     141   291  4.15e-08     57.4     49.68
    2208     559   721     142   305  4.31e-08     58.2     52.87
    2209     463   722      20   287  4.31e-08     57.8     43.53
    2210     451   655      27   227  4.41e-08     57.4     48.80
    2211     459   705      12   265  4.57e-08     57.8     41.76
    2212     463   701      84   329  4.59e-08     58.2     45.85
    2213     559   721     143   306  4.76e-08     57.8     52.60
    2214     463   701      82   327  5.04e-08     57.8     45.85
    2215     462   665      35   237  5.09e-08     57.4     44.76
    2216     463   650      44   233  5.15e-08     57.8     45.10
    2217     451   655      21   221  5.30e-08     57.0     48.80
    2218     459   656      82   282  5.53e-08     57.8     46.83
    2219     463   701      98   343  5.62e-08     57.8     45.85
    2220     459   656      76   276  6.24e-08     57.4     46.83
    2221     451   655      43   243  6.32e-08     57.0     48.80
    2222     459   656      72   272  7.03e-08     57.4     46.83
    2223     459   656      73   273  7.05e-08     57.4     46.83
    2224     463   676      33   246  7.56e-08     57.0     44.50
    2225     459   656      72   272  8.02e-08     57.0     46.83
    2226     451   661      25   260  8.24e-08     57.4     43.60
    2227     459   656      68   268  8.95e-08     57.0     46.83
    2228     451   655      27   227  8.96e-08     56.2     46.41
    2229     459   656      92   292  9.05e-08     57.0     46.83
    2230     451   655      27   227  9.17e-08     56.6     48.33
    2231     459   656      79   279  9.35e-08     57.0     46.83
    2232     459   659      11   216  1.05e-07     56.2     44.91
    2233     459   656      79   279  1.08e-07     56.6     46.83
    2234     452   659      12   230  1.17e-07     56.2     43.67
    2235     456   705       8   264  1.21e-07     56.6     41.95
    2236     463   715      32   305  1.25e-07     55.8     41.52
    2237     552   704     135   280  1.38e-07     55.8     48.08
    2238     463   652      32   231  1.46e-07     55.8     43.84
    2239     459   656     106   306  1.50e-07     56.2     46.34
    2240     459   656      68   268  1.64e-07     56.2     46.83
    2241     513   656     132   273  2.07e-07     55.8     48.97
    2242     513   656     138   279  2.07e-07     55.8     48.97
    2243     559   705     133   269  2.17e-07     55.5     52.32
    2244     513   656     135   276  2.25e-07     55.8     48.97
    2245     513   656     137   278  2.27e-07     55.8     48.97
    2246     513   656     136   277  2.40e-07     55.8     48.97
    2247     513   656     137   278  2.61e-07     55.5     48.97
    2248     513   656     133   274  2.64e-07     55.5     48.97
    2249     459   659      11   216  2.89e-07     55.1     44.91
    2250     513   656     131   272  2.92e-07     55.5     48.97
    2251     459   659      11   216  3.15e-07     54.7     44.91
    2252     459   659      12   217  3.29e-07     54.7     44.91
    2253     459   659      18   223  3.40e-07     54.7     44.91
    2254     459   659      31   236  3.75e-07     54.7     44.91
    2255     459   659      13   218  3.77e-07     54.3     44.91
    2256     459   659      11   216  3.81e-07     54.3     44.91
    2257     459   659      12   217  3.83e-07     54.3     44.91
    2258     459   659      31   236  3.89e-07     54.3     44.91
    2259     459   659      11   216  4.03e-07     54.3     44.91
    2260     459   659      13   218  4.16e-07     54.3     44.91
    2261     459   659       8   213  4.24e-07     54.3     44.91
    2262     459   659      13   218  4.51e-07     54.3     44.91
    2263     459   659      13   218  4.72e-07     53.9     44.91
    2264     459   659      11   216  6.07e-07     53.9     44.91
    2265     459   659      11   216  6.54e-07     53.5     44.44
    2266     459   659      13   218  6.70e-07     53.5     46.01
    2267     459   705      12   265  7.30e-07     53.9     40.61
    2268     559   649     166   253  1.06e-06     54.3     54.26
    2269     463   656      20   190  1.21e-06     52.8     43.96
    2270     231   280      14    67  1.85e-06     48.1     51.85
    2271     540   694     142   293  2.07e-06     52.4     49.04
    2272     437   656      52   270  2.08e-06     52.8     45.92
    2273     461   709      15   287  2.12e-06     52.8     44.01
    2274     452   598      12   161  2.22e-06     52.4     47.13
    2275     452   598      49   198  2.49e-06     52.4     47.13
    2276     526   669     151   287  3.99e-06     51.6     53.69
    2277     460   653      24   231  4.23e-06     52.0     45.37
    2278     559   665     158   260  4.93e-06     51.2     51.38
    2279     460   653      35   242  6.21e-06     51.2     45.37
    2280     459   693      11   254  7.10e-06     50.4     42.00
    2281     460   653      27   234  7.20e-06     50.8     45.37
    2282     460   653      26   233  7.64e-06     50.8     45.37
    2283     460   653      27   234  7.87e-06     50.8     45.37
    2284     460   653      26   233  7.98e-06     50.8     45.37
    2285     460   653      26   233  8.26e-06     50.8     45.37
    2286     452   658      16   222  9.04e-06     50.4     44.65
    2287     452   658      14   220  9.28e-06     50.1     44.65
    2288     463   656      14   241  1.04e-05     50.1     42.67
    2289     460   653      26   233  1.04e-05     50.4     45.37
    2290     460   653      27   234  1.08e-05     50.4     45.37
    2291     526   666     125   258  1.20e-05     50.1     52.05
    2292     460   653      33   240  1.20e-05     50.1     45.37
    2293     452   598      16   165  1.25e-05     50.1     45.91
    2294     526   666     131   264  1.35e-05     50.1     52.05
    2295     526   666     111   244  1.36e-05     49.7     52.05
    2296     526   666     125   258  1.38e-05     50.1     52.05
    2297     460   653      26   233  1.44e-05     49.7     45.37
    2298     556   650     266   359  1.47e-05     50.4     54.00
    2299     526   666     126   259  1.50e-05     49.7     52.05
    2300     526   666     125   258  1.63e-05     49.7     52.05
    2301     452   658      14   220  1.79e-05     49.3     45.37
    2302     526   666     116   249  1.82e-05     49.3     52.05
    2303     526   666     112   245  1.82e-05     49.3     52.05
    2304     452   598       3   152  1.91e-05     48.9     45.91
    2305     526   666     116   249  1.91e-05     49.3     52.05
    2306     448   658      42   247  2.02e-05     49.3     45.16
    2307     234   280       9    59  3.62e-05     43.9     50.98
    2308     537   666     144   263  3.67e-05     48.5     50.76
    2309     234   280       9    59  3.70e-05     43.9     50.98
    2310     452   658      22   228  3.82e-05     48.5     44.91
    2311     452   658      30   236  3.97e-05     48.1     44.91
    2312     537   666     124   243  3.98e-05     48.5     51.52
    2313     452   658      20   226  4.00e-05     48.1     44.91
    2314     537   666     124   243  4.05e-05     48.5     51.52
    2315     452   658      22   228  4.07e-05     48.5     44.91
    2316     452   658      14   220  4.19e-05     48.1     44.91
    2317     537   666     124   243  4.19e-05     48.5     51.52
    2318     537   666     124   243  4.19e-05     48.5     51.52
    2319     452   658      10   216  4.29e-05     48.1     44.91
    2320     452   658      16   222  4.37e-05     48.1     44.91
    2321     452   658      21   227  4.50e-05     48.1     44.91
    2322     234   280       9    59  4.68e-05     43.5     50.98
    2323     452   658      15   221  4.89e-05     48.1     44.91
    2324     452   658      16   222  4.99e-05     48.1     44.91
    2325     452   658      60   266  5.18e-05     48.1     44.91
    2326     234   280       9    59  5.42e-05     43.5     50.98
    2327     452   658      66   272  5.43e-05     48.1     44.91
    2328     235   280       1    50  5.95e-05     42.7     50.00
    2329     235   280       3    52  6.41e-05     42.7     50.00
    2330     453   670      39   254  6.46e-05     47.8     45.49
    2331     235   280      13    62  7.50e-05     43.1     50.00
    2332     234   280       9    59  8.35e-05     42.7     50.98
    2333     526   666     115   248  1.02e-04     47.0     50.68
    2334     537   666     124   243  1.03e-04     47.0     51.52
    2335     233   280       1    52  1.05e-04     45.8     51.92
    2336     461   666      43   243  1.06e-04     47.0     45.87
    2337     526   666     115   248  1.16e-04     47.0     50.68
    2338     559   666     133   237  1.40e-04     46.6     53.64
    2339     459   659      13   220  1.42e-04     46.2     41.59
    2340     459   659      14   221  1.47e-04     46.2     41.59
    2341     537   666     129   248  1.61e-04     46.6     51.52
    2342     537   666     124   243  1.66e-04     46.2     51.52
    2343     519   658     128   266  1.70e-04     46.6     47.26
    2344     537   666     129   248  1.71e-04     46.2     51.52
    2345     461   666      42   242  1.72e-04     46.2     45.41
    2346     461   666      42   242  1.78e-04     46.2     45.41
    2347     537   666     125   244  1.84e-04     46.2     51.52
    2348     461   666      41   241  1.86e-04     46.2     45.41
    2349     537   666     144   263  1.87e-04     46.2     51.52
    2350     532   666     119   243  1.90e-04     46.2     50.36
    2351     537   666     124   243  1.94e-04     46.2     51.52
    2352     537   666     124   243  1.95e-04     46.2     51.52
    2353     537   666     138   257  1.95e-04     46.2     51.52
    2354     537   666     124   243  1.98e-04     46.2     51.52
    2355     461   666      41   241  1.99e-04     46.2     45.41
    2356     234   295     212   273  2.03e-04     46.6     50.00
    2357     537   666     124   243  2.06e-04     46.2     51.52
    2358     537   666     123   242  2.15e-04     45.8     51.52
    2359     537   666     123   242  2.15e-04     45.8     51.52
    2360     537   666     124   243  2.16e-04     45.8     51.52
    2361     537   666     124   243  2.17e-04     45.8     51.52
    2362     537   666     146   265  2.23e-04     46.2     51.52
    2363     537   666     123   242  2.30e-04     45.8     51.52
    2364     537   666     124   243  2.31e-04     45.8     51.52
    2365     537   666     124   243  2.34e-04     45.8     51.52
    2366     537   666     123   242  2.38e-04     45.8     51.52
    2367     537   666     124   243  2.45e-04     45.8     51.52
    2368     537   666     124   243  2.48e-04     45.8     51.52
    2369     537   666     122   241  2.49e-04     45.8     51.52
    2370     537   666     122   241  2.57e-04     45.8     51.52
    2371     537   666     137   256  2.72e-04     45.8     51.52
    2372     231   280       7    60  2.89e-04     42.0     48.15
    2373     537   666     146   265  3.01e-04     45.4     51.52
    2374     537   666     147   266  3.02e-04     45.4     51.52
    2375     537   666     147   266  3.05e-04     45.4     51.52
    2376     537   666     129   248  3.19e-04     45.4     50.76
    2377     537   666     171   290  3.29e-04     45.8     51.52
    2378     532   666     119   243  3.34e-04     45.4     50.36
    2379     235   379     210   342  3.36e-04     45.8     39.33
    2380     459   659      12   219  3.37e-04     45.4     41.48
    2381     459   659      10   217  3.69e-04     45.1     42.73
    2382     559   666     139   243  5.10e-04     44.7     53.64
    2383     451   653      50   263  7.05e-04     44.3     45.25
    2384     234   277       4    46  7.49e-04     39.7     56.82
    2385     458   655      33   221  1.00e-03     43.5     39.81
    2386     458   655      33   221  2.00e-03     42.7     39.34
    2387     559   699     154   294  2.00e-03     43.1     50.00
    2388     496   655     152   300  2.00e-03     43.1     42.33
    2389     458   655      21   209  2.00e-03     42.4     38.94
    2390     559   699     135   275  3.00e-03     42.4     50.00
    2391     559   666     134   238  3.00e-03     42.4     52.73
    2392     559   699     133   273  3.00e-03     42.4     50.00
    2393     559   666     133   237  3.00e-03     42.4     52.73
    2394     559   666     133   237  3.00e-03     42.4     52.73
    2395     559   699     134   274  4.00e-03     42.0     50.00
    2396     559   699     133   273  4.00e-03     42.0     50.00
    2397     559   699     133   273  4.00e-03     42.0     50.00
    2398     559   699     134   274  4.00e-03     42.0     50.00
    2399     559   699     133   273  4.00e-03     42.0     50.00
    2400     235   272      14    55  6.00e-03     37.7     52.38
    2401     240   273     106   143  1.40e-02     38.9     60.53
    2402     452   658      16   201  1.40e-02     40.0     41.40
    2403     462   598      38   186  2.50e-02     39.7     45.86
    2404     218   280       6    67  4.80e-02     35.4     47.76
    2405     233   281      45    97  4.80e-02     39.3     49.06
    2406     233   281      37    89  6.40e-02     38.5     49.06
    2407     463   528     216   278  7.30e-02     37.7     55.22
    2408     216   273     151   209  8.40e-02     37.4     47.46
    2409     458   670      25   251  8.70e-02     37.7     41.60
    2410     458   670      25   251  9.30e-02     37.7     41.60
    2411     235   280      24    73  1.10e-01     34.3     46.00
    2412     235   280       8    57  5.60e-01     32.3     44.00
    2413     235   278      18    65  5.80e-01     32.3     43.75
    2414     257   280       4    27  8.50e-01     33.9     58.33
    2415     456   655      22   242  8.50e-01     34.7     38.56
    2416     463   600      28   194  3.60e+00     32.7     37.50
    2417     235   281       9    58  5.90e+00     29.3     42.00

    $url
                                                                                                                                                                              V61TE0AB016 
    "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=Alignment&ALIGNMENT_VIEW=Tabular&RESULTS_FILE=on&FORMAT_TYPE=CSV&ALIGNMENTS=20000&DESCRIPTIONS=20000&RID=V61TE0AB016" 

    attr(,"class")
    [1] "blast"
