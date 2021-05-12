#!/usr/bin/env Rscript

# 5-methylation_calls.R
#
# Jay Moore, John Innes Centre, Norwich
# jonathan.moore@jic.ac.uk
#
# 29/09/2017
#
# Version 2
#
# Change log
#
# Version 1 made methylation calls at each C individually, then merged calls for CG dinucleotides. This demonstrated a very high level of concordance between adjacent nucleotides in CG dinucleotide pairs.
# Accordingly, this version makes composite methylation calls for pairs of adjacent CG dinucleotides in order to provide more power to the methylation calling algorithm to distinguish between states at a given site.
# In addition, this version develops a generation- and replicate-aware analysis, to assign methylation status to sites in concordance with replicate structure, and to identify sites whose methylation status changes between generations in one or more line.  
# This version also extends the analysis of sites overlapping genes, to generate a set of genic sites whose status remains M, U, or changes in one direction or the other between generations.


# Published data sources (project_id)
# Schmitz et al, 2011 - NCBI SRA accession SRA035939
# Becker et al, 2011 - EMBL ENA accession PRJEB2678

# Summary of functions
#
# Done:
# Read GFF files output from bs_sequel for each sample in a project - each file contains C and T calls at all sites in a meth_context
# Align and merge the C and T calls at sites across samples
# Plot summary data about samples, sites
# Assign methylation status at each site in each sample
#
# This script calls methylation status of a C site based on:
# Fisher's exact test -> "I" (identify sites where coverage is too low to make a clear distinction between methylated and unmethylated status given the conversion efficiency and sequencing and alignment errors for the sample)
# Binomial test with B-H FDH adjusted p-values >= 0.01 -> "U" (identify sites where the number of "C" calls is too few to make a call of methylated for the sample)
# #C >= #T -> differentiate between "M" and "P" (impose an expectation that methylated sites will be represented by at least half of aligned reads being converted, else call Partial methylation)
# Coverage < 5 sigma -> "I" (assume coverage is approxinmately Gaussian distribution, and ignore sites where coverage is too high to be meaningful for the sample)
#
# Once all sites have been called, the script combines adjacent C and G sites identified from the reference and addresses the concordance of calls in adjacent CG dinucleotides.
#
# Bring CG site combined calls into same structure as CHG/CHH calls
#
# Underway:
#
# How many sites with reliable methylation status calls vary between replicates in each sample?
# Which genes are really genes (identify and ignore 'heterochromatic genes'
# What are the clusters of methylated sites across all lines?
#
# To do:
# 
# Address SNPs and remove sites with SNPs in any sample, and adjacent sites where SNP affects a CG site or CH* site
# Estimate forward and reverse rates for M->U and U->M
# Read GFF files for annotations and annotate each meth_context locus
# Carry out annotation-dependent analyses of methylation status
#
#
# Results:
#
# Combining #C and #T from adjacent sites in CG dinucleotides, prior to calling, increases the number of unambiguous (M/U) calls substantially in both experiments
# In the Schmitz data, proportion increases from 83% to 88%, in the Becker data from 51% to 59%
# In both data sets the fraction of sites called as "P" also increases from ~1% to ~1.5%
#
#	Schmitz	Becker				
#C + #T basis	4770529	32929537	I	0.106401687	0.391711504	
#	12914624	17178020	M	0.288047254	0.204340196	
#	810736	1254341	P	0.018082623	0.014920945	
#	26339199	32703892	U	0.587468435	0.389027356	
#	44835088	84065790				
#				0.87551569	0.593367552	Proportion M/U
#						
#Combining calls	7018720	40652504	I	0.156545249	0.483579634	
#	12320603	14225184	M	0.274798234	0.169214897	
#	613851	922635	P	0.013691308	0.010975154	
#	24881914	28265467	U	0.554965209	0.336230314	
#	44835088	84065790				
#				0.829763443	0.505445211	Proportion M/U
#
#
# This analysis, combining #C and #T at adjacent sites in CG dinucleotides before making methylation calls, substantially increases the enrichment of the variable sites for smaller gaps than randomly chosen 'all M' or 'all U' sites. For gaps between sites<50 nt, enrichment is ~ 2x
#
#
# Concordance of calls among reps:
#
# Schmitz data:
# CG
#                   D      I      M    O     P       U
#Line1Gen3     700561  47046 672612 1186 30594 1350194
#Line19Gen3    582266  44225 745651 1747 30437 1397867
#Line29Gen30   278340  40328 792257 2946 35396 1652926
#Line12Gen3     50738  14118 843796 2020 40999 1850522
#Line49Gen30   386491  68604 791920 1659 33344 1520175
#Line59Gen30   500219 153478 750750 2232 27911 1367603
#Line119Gen30 1070180  41982 703868 1829 22375  961959
#
# After converting "P" to "U" with cutoff mCG=0.45:
#                   D      I      M     O P       U
#Line1Gen3     674255  47046 675376 10525 0 1394991
#Line19Gen3    553141  44225 749315 14526 0 1440986
#Line29Gen30   246889  40328 796326 15033 0 1703617
#Line12Gen3     14980  14118 849302 15491 0 1908302
#Line49Gen30   356939  68604 796117 13958 0 1566575
#Line59Gen30   471082 153478 754206 16168 0 1407259
#Line119Gen30 1043287  41982 706586 17784 0  992554
#
# Using p<0.05 in Fisher's Exact test and p<0.005 in Binomial test:
#                  D     I      M    O     P       U
#Line1Gen3    464735 27304 718113 2246 29384 1560411
#Line19Gen3   395339 24927 773073 2765 27442 1578647
#Line29Gen30  178635 24208 815170 4009 32425 1747746
#Line12Gen3    37852 12502 851617 2839 36849 1860534
#Line49Gen30  279957 40318 813520 2606 30468 1635324
#Line59Gen30  360904 80352 783987 3511 24767 1548672
#Line119Gen30 955561 26175 734659 2820 19524 1063454
#
#
# CHG
#                   D       I      M     O     P       U
#Line1Gen3    2397131 1677900 184629 44154   430 1810912
#Line19Gen3   2028733 1697219 409875 30186 26936 1922207
#Line29Gen30  1403778 1449528 429026 41130 32021 2759673
#Line12Gen3    927203  768598 644693 41795 60240 3672627
#Line49Gen30  1304256 1842457 379046 94499   936 2493962
#Line59Gen30  1461453 2325729 359560 41041 16827 1910546
#Line119Gen30 2925630 1482948 382214 26260 11056 1287048
#
#
# Becker untrimmed data:
# CG
#                   D       I      M    O    P       U
#Line4Gen3    1032292  456328 465033 3765 6295  838480
#Line8Gen3     592155  980252 422285 8997 5423  793081
#Line109Gen31  536733  798189 484595 2681 3664  976331
#Line119Gen31  406872  830067 529903 2448 4727 1028176
#Line29Gen31   827466  314559 518567 2573 3826 1135202
#Line39Gen31   756156  599624 395918 2591 1631 1046273
#Line39Gen32   415544 1082574 457008 3010 4528  839529
#Line49Gen31   436353  438682 584406 2702 4057 1335993
#Line49Gen32   285970 1420354 379631 3108 4520  708610
#Line59Gen31   356846  301633 668957 2959 6155 1465643
#Line79Gen31   405627  595338 591730 4002 5568 1199928
#Line89Gen31   831138  608389 480077 2006 5612  874971
#Line99Gen31   890205  696108 426574 2600 2290  784416
#
# Using p<0.05 in Fisher's Exact test and p<0.005 in Binomial test:
#                   D       I      M     O    P       U
#Line4Gen3    1025077  251730 524733  7594 4696  988363
#Line8Gen3     614404  750497 488875 13788 3846  930783
#Line109Gen31  597408  529116 566061  6331 2531 1100746
#Line119Gen31  426723  546761 605729  6113 3225 1213642
#Line29Gen31   677916  165352 609273  6150 2581 1340921
#Line39Gen31   734605  349741 507015  5721  967 1204144
#Line39Gen32   575076  778719 509527  6288 3061  929522
#Line49Gen31   311392  212298 680916  7533 2777 1587277
#Line49Gen32   272700 1277601 447505  6005 3165  795217
#Line59Gen31   258006  157985 730211  7635 4416 1643940
#Line79Gen31   371178  334465 664513  9371 3116 1419550
#Line89Gen31   799590  421540 541453  4603 4077 1030930
#Line99Gen31   961156  421592 493318  5969 1371  918787

# CHG
#                   D       I    M    O  P       U
#Line4Gen3    1125465 1690897 4252 1389  7 1235199
#Line8Gen3     609271 2432588 3156 1352  0 1010842
#Line109Gen31  664280 2354336 2775  410  0 1035408
#Line119Gen31  468511 2389279 5603 1200  4 1192612
#Line29Gen31  1428013 1478126 1613 2509  0 1146948
#Line39Gen31  1251179 1867520 1094  481  0  936935
#Line39Gen32   365539 2621883 3556  659  1 1065571
#Line49Gen31   743526 1672952 2878  988  5 1636860
#Line49Gen32   284362 2883743 3442  943 11  884708
#Line59Gen31   717077 1389237 5650 4335 12 1940898
#Line79Gen31   587909 1998867 6814 2914  7 1460698
#Line89Gen31  1085499 1938008 3252 1747  0 1028703
#Line99Gen31  1004838 2144409 3073  858  0  904031
#
#
# Concordant M or U calls across reps (excluding line 69) range from 59%-96% of sites among the Schmitz data (mean 78%), and 39%-76% among the Becker untrimmed data (mean 53%)
# 
#
# Summarising concordance of consistently replicated M and U calls across all lines (excluding line 69) :
#
# Schmitz data:
#No of CG loci: 2802193 
#No of CG loci with all calls M or U: 1318815 
#No of CG loci with all calls M: 546168 
#No of CG loci with all calls U: 741378 
#No of CG loci with variation in M/U calls: 31269 
#No of CG loci with variation in M/U calls, M at generation 3: 13556 
#No of CG loci with variation in M/U calls, U at generation 3: 11351 
#
# After converting 'P' to 'U' with mCG cutoff=0.45:
#No of CG loci: 2802193 
#No of CG loci with all calls M or U: 1392328 
#No of CG loci with all calls M: 550539 
#No of CG loci with all calls U: 792125 
#No of CG loci with variation in M/U calls: 49664 
#No of CG loci with variation in M/U calls, M at generation 3: 18593 
#No of CG loci with variation in M/U calls, U at generation 3: 16543 
#
# Using p<0.05 in Fisher's Exact test and p<0.005 in Binomial test:
#No of CG loci: 2802193 
#No of CG loci with all calls M or U: 1531054 
#No of CG loci with all calls M: 592062 
#No of CG loci with all calls U: 898177 
#No of CG loci with strict variation in M/U calls (all parentals agree, all lines have agreement between reps): 40815 
#No of CG loci with strict variation in M/U calls, M at generation 3: 17078 
#No of CG loci with strict variation in M/U calls, U at generation 3: 15108 
#No of CG loci with any variation in M/U calls (at least 2 parentals agree, at least one progeny line disagrees and has agreement between reps): 70768 
#No of CG loci with any variation in M/U calls, M at generation 3: 35852 
#No of CG loci with any variation in M/U calls, U at generation 3: 34916 
#
#  We have added the more sensitive version of 'variable' sites, where we use the consensus of parents, rather than unanimity, and where we only ask that one progeny line has enough coverage to make a variant call, rather than requiring all progeny to have clean rep-wise calls.  This increases the number of identified variable sites by 73%.  The other effect of this is to reduce the disparity between counts of M->U and U->M transitions, from 53:47 to 51:49.


#No of CHG loci: 6115156 
#No of CHG loci with all calls M or U: 1014846 
#No of CHG loci with all calls M: 141337 
#No of CHG loci with all calls U: 872486 
#No of CHG loci with variation in M/U calls: 1023 
#No of CHG loci with variation in M/U calls, M at generation 3: 139 
#No of CHG loci with variation in M/U calls, U at generation 3: 197 
#
# Proportion of CG sites with concordant M/U calls per line:
#   Line1Gen3   Line19Gen3  Line29Gen30   Line12Gen3  Line49Gen30  Line59Gen30 Line119Gen30 
#   0.7218653    0.7649430    0.8725962    0.9615034    0.8251020    0.7559626    0.5944726 
# Proportion of CG sites in each line with rep-concordant M/U call within transposons: 
#   Line1Gen3   Line19Gen3  Line29Gen30   Line12Gen3  Line49Gen30  Line59Gen30 Line119Gen30 
#   0.7843881    0.8666747    0.8991857    0.9659413    0.8892960    0.8609134    0.8162608 
# Proportion of CHG sites with concordant M/U calls per line:
#   Line1Gen3   Line19Gen3  Line29Gen30   Line12Gen3  Line49Gen30  Line59Gen30 Line119Gen30 
#   0.3263271    0.3813610    0.5214420    0.7060032    0.4698176    0.3712262    0.2729713 
#
# After converting 'P' to 'U' and using mCG cutoff=0.45:
#   Line1Gen3   Line19Gen3  Line29Gen30   Line12Gen3  Line49Gen30  Line59Gen30 Line119Gen30 
#   0.7388381    0.7816382    0.8921381    0.9840878    0.8431582    0.7713477    0.6063608    
#
# Using p<0.05 for Fisher's Exact test and p<0.005 in Binomial test:
#   Line1Gen3   Line19Gen3  Line29Gen30   Line12Gen3  Line49Gen30  Line59Gen30 Line119Gen30 
#   0.8131217    0.8392427    0.9146108    0.9678673    0.8739027    0.8324405    0.6416806 
#

# Becker data:
#No of CG loci: 2802193 
#No of CG loci with all calls M or U: 765045 
#No of CG loci with all calls M: 211457 
#No of CG loci with all calls U: 530755 
#No of CG loci with variation in M/U calls: 22833 
#No of CG loci with variation in M/U calls, M at generation 3: 11152
#No of CG loci with variation in M/U calls, U at generation 3: 11169 
#No of CG loci with variation in M/U calls, M at generation 3 and 30: 45 
#No of CG loci with variation in M/U calls, U at generation 3 and 30: 46 
#
#No of CHG loci: 6100508 
#No of CHG loci with all calls M or U: 602302 
#No of CHG loci with all calls M: 4333 
#No of CHG loci with all calls U: 597949 
#No of CHG loci with variation in M/U calls: 20 
#No of CHG loci with variation in M/U calls, M at generation 3: 1 
#No of CHG loci with variation in M/U calls, U at generation 3: 19 
#No of CHG loci with variation in M/U calls, M at generation 3 and 30: 1 
#No of CHG loci with variation in M/U calls, U at generation 3 and 30: 0 
#
# Using p<0.05 for Fisher's Exact test and p<0.005 in Binomial test:
#No of CG loci: 2802193 
#No of CG loci with all calls M or U: 933133 
#No of CG loci with all calls M: 283950 
#No of CG loci with all calls U: 614931 
#No of CG loci with strict variation in M/U calls (all parentals agree, all lines have agreement between reps): 34252 
#No of CG loci with strict variation in M/U calls, M at generation 3: 17224 
#No of CG loci with strict variation in M/U calls, U at generation 3: 16138 
#No of CG loci with any variation in M/U calls (at least 2 parentals agree, at least one progeny line disagrees and has agreement between reps): 65864 
#No of CG loci with any variation in M/U calls, M at generation 3: 36166 
#No of CG loci with any variation in M/U calls, U at generation 3: 29698 

# Proportion of CG sites with concordant M/U calls per line:
#   Line4Gen3    Line8Gen3 Line109Gen31 Line119Gen31  Line29Gen31  Line39Gen31  Line39Gen32  Line49Gen31  Line49Gen32  Line59Gen31  Line79Gen31  Line89Gen31  Line99Gen31 
#   0.4651760    0.4337196    0.5213510    0.5560213    0.5901696    0.5146651    0.4626865    0.6853200    0.3883533    0.7617605    0.6393771    0.4835670    0.4321580 
#
# Using p<0.05 for Fisher's Exact test and p<0.005 in Binomial test:
#   Line4Gen3    Line8Gen3 Line109Gen31 Line119Gen31  Line29Gen31  Line39Gen31  Line39Gen32  Line49Gen31  Line49Gen32  Line59Gen31  Line79Gen31  Line89Gen31  Line99Gen31 
#   0.5399685    0.5066239    0.5948223    0.6492668    0.6959528    0.6106499    0.5135439    0.8094350    0.4434819    0.8472475    0.7437257    0.5611259    0.5039285 

# Proportion of CG sites in each line with rep-concordant M/U call within non-heterochromatic genes: 
#   Line4Gen3    Line8Gen3 Line109Gen31 Line119Gen31  Line29Gen31  Line39Gen31  Line39Gen32  Line49Gen31  Line49Gen32  Line59Gen31  Line79Gen31  Line89Gen31  Line99Gen31 
#   0.5408167    0.5081718    0.6071344    0.6454972    0.6885090    0.6146135    0.5376677    0.7834307    0.4536456    0.8516156    0.7368869    0.5548744    0.4965213 
#
# Proportion of CG sites in each line with rep-concordant M/U call within transposons: 
#   Line4Gen3    Line8Gen3 Line109Gen31 Line119Gen31  Line29Gen31  Line39Gen31  Line39Gen32  Line49Gen31  Line49Gen32  Line59Gen31  Line79Gen31  Line89Gen31  Line99Gen31 
#   0.5448761    0.5046653    0.5702040    0.6083041    0.5865765    0.4726979    0.5411500    0.6456594    0.4637271    0.7279808    0.6681512    0.5688265    0.5180915 
#   
# Proportion of CHG sites with concordant M/U calls per line:
#   Line4Gen3    Line8Gen3 Line109Gen31 Line119Gen31  Line29Gen31  Line39Gen31  Line39Gen32  Line49Gen31  Line49Gen32  Line59Gen31  Line79Gen31  Line89Gen31  Line99Gen31 
#   0.2409027    0.1941471    0.1991636    0.2363275    0.2143936    0.1755396    0.2088780    0.3174052    0.1721652    0.3838067    0.2931130    0.2000604    0.1739580 
#
# Proportion of CHG sites in each line with rep-concordant M/U call within non-heterochromatic genes: 
#   Line4Gen3    Line8Gen3 Line109Gen31 Line119Gen31  Line29Gen31  Line39Gen31  Line39Gen32  Line49Gen31  Line49Gen32  Line59Gen31  Line79Gen31  Line89Gen31  Line99Gen31 
#   0.3054935    0.2499250    0.2558860    0.2953299    0.2830914    0.2312006    0.2635129    0.4041542    0.2189066    0.4797751    0.3617048    0.2543510    0.2235783 
#
# Proportion of CHG sites in each line with rep-concordant M/U call within transposons: 
#   Line4Gen3    Line8Gen3 Line109Gen31 Line119Gen31  Line29Gen31  Line39Gen31  Line39Gen32  Line49Gen31  Line49Gen32  Line59Gen31  Line79Gen31  Line89Gen31  Line99Gen31 
#  0.13580409   0.10020552   0.09714683   0.14409826   0.07149915   0.05455186   0.12417623   0.12923320   0.10001737   0.17979043   0.18589045   0.11108174   0.08753462 
#
#
# Gaps between adjacent variable CG sites are very strongly enriched (up to 4fold) for gap sizes below 50nt in both data sets (see plots).  This enrichment increases with decreasing gap size.  The implication of this result is that variability of methylation at CG sites is clustered, to a much greater extent than methylated sites are clustered.
#
#
# Sites which vary between lines (excluding line 69), with consistency between reps, are very strongly enriched in annotated genes:
#
# Gene-wise correlation between average methylation across all samples and count of 'all M' sites as proportion of all sites which are M/U in all samples:
# Schmitz result: r=0.9542264.  After converting P to U/M with 0.45 cutoff:  0.953940455588698. Using p<0.05 Fisher's and p<0.005 Binomial: 0.960349590157647
# Becker result: r=0.873724653313617
#
# Transposon-wise correlation between average methylation across all samples and count of 'all M' sites as proportion of all sites which are M/U in all samples:
# Schmitz result: r=???.  After converting P to U/M with 0.45 cutoff:  0.954046877657042. Using p<0.05 Fisher's and p<0.005 Binomial: 0.963762033310974
# Becker result: 0.970635542056782
#
# Schmitz data for non-heterochromatic genes:
#
#1620783 of 2802193 ( 0.578398061803737 )  CG  sites lie within annotated non-heterochromatic genes.
#91367 of 546168 ( 0.16728735480658 ) 'all M' sites lie within annotated non-heterochromatic genes.
#549934 of 741378 ( 0.74177275290068 ) 'all U' sites lie within annotated non-heterochromatic genes.
#27640 of 31269 ( 0.883942562921744 ) variable sites lie within annotated non-heterochromatic genes.
#12249 of 13556 ( 0.903585128356447 ) M->U transition sites lie within annotated non-heterochromatic genes.
#10152 of 11351 ( 0.894370540040525 ) U->M transition sites lie within annotated non-heterochromatic genes.
#
#516901 of 2802193 ( 0.184463025922911 ) CG sites lie within annotated transposons.
#328323 of 546168 ( 0.601139209913433 ) 'all M' sites lie within annotated transposons.
#20219 of 741378 ( 0.02727218773689 ) 'all U' sites lie within annotated transposons.
#1176 of 31269 ( 0.0376091336467428 ) variable sites lie within annotated transposons.
#434 of 13556 ( 0.032015343759221 ) M->U transition sites lie within annotated transposons.
#360 of 11351 ( 0.0317152673773236 ) U->M transition sites lie within annotated transposons.
#
# How many CG, variable and invariable sites lie within annotated exons?
#1395509 of 2802193 ( 0.498006025994641 )  CG  sites lie within annotated exons.
#72893 of 546168 ( 0.13346259758902 ) 'all M' sites lie within annotated exons.
#457787 of 741378 ( 0.617481230897059 ) 'all U' sites lie within annotated exons.
#21726 of 31269 ( 0.69480955579008 ) variable sites lie within annotated exons.
#9558 of 13556 ( 0.705075243434641 ) M->U transition sites lie within annotated exons.
#8035 of 11351 ( 0.707867148268875 ) U->M transition sites lie within annotated exons.
#
#511281 of 2802193 ( 0.18245745385846 )  CG  sites lie within annotated introns.
#48520 of 546168 ( 0.0888371343615884 ) 'all M' sites lie within annotated introns.
#171042 of 741378 ( 0.230708221716857 ) 'all U' sites lie within annotated introns.
#13675 of 31269 ( 0.437334100866673 ) variable sites lie within annotated introns.
#6527 of 13556 ( 0.48148421363234 ) M->U transition sites lie within annotated introns.
#4619 of 11351 ( 0.406924500044049 ) U->M transition sites lie within annotated introns.
#
# How many CG, variable and invariable sites lie within annotated 5'UTRs?
#177892 of 2802193 ( 0.0634831362436492 )  CG  sites lie within annotated 5'UTRs.
#923 of 546168 ( 0.0016899562039519 ) 'all M' sites lie within annotated 5'UTRs.
#91142 of 741378 ( 0.122935938212356 ) 'all U' sites lie within annotated 5'UTRs.
#333 of 31269 ( 0.010649525088746 ) variable sites lie within annotated 5'UTRs.
#98 of 13556 ( 0.007229271171437 ) M->U transition sites lie within annotated 5'UTRs.
#156 of 11351 ( 0.0137432825301736 ) U->M transition sites lie within annotated 5'UTRs.
#
# How many CG, variable and invariable sites lie within annotated 3'UTRs?
#141368 of 2802193 ( 0.050449059004858 )  CG  sites lie within annotated 3'UTRs.
#4396 of 546168 ( 0.0080488054957449 ) 'all M' sites lie within annotated 3'UTRs.
#50018 of 741378 ( 0.0674662587775737 ) 'all U' sites lie within annotated 3'UTRs.
#1678 of 31269 ( 0.053663372669417 ) variable sites lie within annotated 3'UTRs.
#583 of 13556 ( 0.0430067866627324 ) M->U transition sites lie within annotated 3'UTRs.
#728 of 11351 ( 0.0641353184741432 ) U->M transition sites lie within annotated 3'UTRs.
#
#
# Becker data for non-heterochromatic genes:
#
#
#4057209 of 6100508 ( 0.66506084411331 ) CHG sites lie within annotated genes.
#468 of 4333 ( 0.10800830833141 ) 'all M' sites lie within annotated genes.
#553478 of 597949 ( 0.925627436453611 ) 'all U' sites lie within annotated genes.
#6 of 20 ( 0.3 ) variable sites lie within annotated genes.
#0 of 1 ( 0 ) M->U transition sites lie within annotated genes.
#6 of 19 ( 0.315789473684211 ) U->M transition sites lie within annotated genes.
#
#1036390 of 6100508 ( 0.169885852129036 ) CHG sites lie within annotated transposons.
#3228 of 4333 ( 0.744980383106393 ) 'all M' sites lie within annotated transposons.
#18323 of 597949 ( 0.0306430816006047 ) 'all U' sites lie within annotated transposons.
#9 of 20 ( 0.45 ) variable sites lie within annotated transposons.
#0 of 1 ( 0 ) M->U transition sites lie within annotated transposons.
#9 of 19 ( 0.473684210526316 ) U->M transition sites lie within annotated transposons.
#
# Correlation of 'heterochromatic gene' calls between Schmitz and Becker data:
# There is high correlation between the per-gene 'average methylation' metric between the two studies (r=0.99173)
# 154 genes have either no CG sites or no coverage in either study
# 92 genes have data in one study but not the other (13 have data only in the Becker study, the remainder only in Schmitz)
# 80 genes have opposite calls in each study
# 33015 genes are agreed between the two studies
#
# Heterochromatic genes current threshold is 0.6.  What if we raise this to 0.8?  This maximises agreement between the Schmitz and Becker data sets.  
# In addition, the genes with average methylation between 0.6 and 0.8 have variable sites which look more similar to regular genes than to heterochromatic genes
#
# Schmitz data:
# 	10610 of 33341 genes have one or more sites which vary with concordance among replicates (0.318226807834198).
#	298 of 2287 heterochromatic genes have one or more sites which vary with concordance among replicates (0.130301705290774).
#	192 of 557 heterochromatic genes with average methylation <=0.8 have one or more sites which vary with concordance among replicates (0.344703770197487).
#	106 of 1897 heterochromatic genes with average methylation >0.8 have one or more sites which vary with concordance among replicates (0.0558777016341592).
#
# Becker data:
#   8174 of 33341 genes have one or more sites which vary with concordance among replicates (0.245163612369155).
#   240 of 2276 heterochromatic genes have one or more sites which vary with concordance among replicates (0.105448154657293).
#   149 of 595 heterochromatic genes with average methylation <=0.8 have one or more sites which vary with concordance among replicates (0.250420168067227).
#   91 of 1914 heterochromatic genes with average methylation >0.8 have one or more sites which vary with concordance among replicates (0.0475444096133751).
#
#
# Above, heterochromatic genes were defined as those with average methylation exceeding 0.6.  A further refinement can be considered by comparing the variance among methylation level of CG sites overlapping the gene.  In TEs, typically all CG sites are methylated, so variance in methylation level among sites is close to zero.  In gene-body methylated genes, some sites are highly methylated and some sites not, typically, in a mosaic pattern, leading to higher values of variance of methylation levels among sites (e.g. SRA035939_CG_methylation_level_variance_by_gene_length.PDF).  There is a clear group of genes, relatively short (<2.5kb), with moderate levels of methylation (>0.2, <0.4), but with very low levels of variance in methylation (<0.1) in comparison to other genes of similar length and avaerage methylation (290 genes).  On inspection, these genes predominantly comprise ChrM genes (146 of 152 mitochondrial genes), and genes from the pericentromeric unmethylated region on Chr2 between 3.24 and 3.50 MB (127 of).  Visual inspection of methylation profiles of a number of these genes led us to exclude them from being categorised as 'gene-body methylated', as their patterns of methylation are distinct from the patterns seen in most other 'typical' gene-body-methylated genes. The mitochondrial genes are highlighted in e.g. SRA035939_CG_methylation_level_variance_by_gene_length_CHRM.PDF.  The reason for low variance being detected in methylation of individual sites in mitochondrial genes may be due to large copy numbers leading to an 'averaging out' of methylation levels among copies.  The reasons for the pericentromeric segment of chromosome 2 looking 'mitochondrial' in its methylation status are to be determined. 
#
# Schmitz data on genes and transposons relative to 'heterochromatic' methylation variance threshold:
#33341 annotated genes
#31219 annotated genes with average_methylation<0.6
#30893 annotated genes with average_methylation<0.6 and dispersion index>0.25 where average_methylation>0.2 or average_methylation<0.2
#22234 annotated genes with methylation dispersion index>0.25
#21949 annotated genes with average methylation<0.6 and dispersion index>0.25
#
#31189 annotated transposons
#14324 annotated transposons with average_methylation<0.6
#14064 annotated transposons with average_methylation<0.6 and dispersion index >0.251744745537775 where average_methylation>0.2 or average_methylation<0.2
#14644 annotated transposons with methylation dispersion index>0.251744745537775
#12251 annotated transposons with average methylation<0.6 and dispersion index>0.251744745537775
#
# 2623 genes are assigned as heterochromatic.
#
# After P->U/M with mCGcutoff 0.45:
#33341 annotated genes
#31230 annotated genes with average_methylation<0.6
#31182 annotated genes with average_methylation<0.6 and dispersion index >0.276878180849435 where average_methylation>0.2 or average_methylation<0.2
#28595 annotated genes with methylation dispersion index>0.276878180849435
#28347 annotated genes with average methylation<0.6 and dispersion index>0.276878180849435
#
#31189 annotated transposons
#14380 annotated transposons with average_methylation<0.6
#14188 annotated transposons with average_methylation<0.6 and dispersion index >0.276878180849435 where average_methylation>0.2 or average_methylation<0.2
#15905 annotated transposons with methylation dispersion index>0.276878180849435
#13548 annotated transposons with average methylation<0.6 and dispersion index>0.276878180849435
#
# 2324 genes are assigned as heterochromatic.
#
# Using Fishers p<0.05, Binomial p<0.005:
#
#33341 annotated genes
#31206 annotated genes with average_methylation<0.6
#30890 annotated genes with average_methylation<0.6 and dispersion index >0.249901683365147 where average_methylation>0.2 or average_methylation<0.2
#23582 annotated genes with methylation dispersion index>0.249901683365147
#23287 annotated genes with average methylation<0.6 and dispersion index>0.249901683365147
#
#31189 annotated transposons
#14276 annotated transposons with average_methylation<0.6
#14029 annotated transposons with average_methylation<0.6 and dispersion index >0.249901683365147 where average_methylation>0.2 or average_methylation<0.2
#14880 annotated transposons with methylation dispersion index>0.249901683365147
#12450 annotated transposons with average methylation<0.6 and dispersion index>0.249901683365147
#
# 2627 genes are assigned as heterochromatic.
#
# After splitting P sites between M and U at mCG=0.45:
#12455 of 33341 genes have one or more sites which vary with concordance among replicates (0.37356408026154).
#573 of 2324 heterochromatic genes have one or more sites which vary with concordance among replicates (0.246557659208262).
#275 of 596 heterochromatic genes with average methylation <=0.8 have one or more sites which vary with concordance among replicates (0.461409395973154).
#298 of 1895 heterochromatic genes with average methylation >0.8 have one or more sites which vary with concordance among replicates (0.157255936675462).
#2501 of 31189 transposons have one or more sites which vary with concordance among replicates (0.080188528006669).
#
# Using Fishers p<0.05, Binomial p<0.005:
#12041 of 33341 genes have one or more sites which vary with concordance among replicates (0.361146936204673).
#371 of 2627 heterochromatic genes have one or more sites which vary with concordance among replicates (0.141225732775029).
#228 of 881 heterochromatic genes with average methylation <=0.8 have one or more sites which vary with concordance among replicates (0.258796821793417).
#143 of 1912 heterochromatic genes with average methylation >0.8 have one or more sites which vary with concordance among replicates (0.0747907949790795).
#993 of 31189 transposons have one or more sites which vary with concordance among replicates (0.0318381480650229).
#
# Using more sensitive test for variable sites:
#14455 of 33341 genes have one or more sites which vary with concordance among replicates (0.43355028343481).
#453 of 2627 heterochromatic genes have one or more sites which vary with concordance among replicates (0.172440045679482).
#281 of 881 heterochromatic genes with average methylation <=0.8 have one or more sites which vary with concordance among replicates (0.318955732122588).
#172 of 1912 heterochromatic genes with average methylation >0.8 have one or more sites which vary with concordance among replicates (0.0899581589958159).
#1800 of 31189 transposons have one or more sites which vary with concordance among replicates (0.0577126551027606).


#
# More refined classification:
# Counts of genes by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
#                 347                15095                 1762                 1683                14454 
# Counts of transposons by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
#                8355                 1810                13683                 3742                 3598 

#
#9788 of 15095 genes with gene body methylation (0.648426631334879) have one or more CG sites which vary.
#822 of 18246 genes without gene body methylation (0.045050970075633) have one or more CG sites which vary.
#657 of 14454 unmethylated genes (0.0454545454545455) have one or more CG sites which vary.
#218 of 1811 transposons with gene body methylation (0.120375483158476) have one or more CG sites which vary.
#529 of 29379 transposons without gene body methylation (0.0180060587494469) have one or more CG sites which vary.
#74 of 3599 unmethylated transposons (0.0205612670186163) have one or more CG sites which vary.
#
#
# After conversion of P to U/M with mCG cutoff=0.45:
# Counts of genes by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
#                 347                15110                 1787                 1377                14720 
#Counts of transposons by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
#                8355                 1915                13810                 3251                 3857 
#11192 of 15110 genes with gene body methylation (0.740701522170748) have one or more CG sites which vary.
#1263 of 18231 genes without gene body methylation (0.069277604080961) have one or more CG sites which vary.
#837 of 14720 unmethylated genes (0.0568614130434783) have one or more CG sites which vary.
#449 of 1916 transposons with gene body methylation (0.234342379958246) have one or more CG sites which vary.
#2052 of 29274 transposons without gene body methylation (0.070096331215413) have one or more CG sites which vary.
#123 of 3858 unmethylated transposons (0.0318818040435459) have one or more CG sites which vary.
#
#
# Using Fisher's p<0.05, Binomial p<0.005:
# Counts of genes by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
#                 347                15245                 1783                 1372                14594 
#
#Counts of transposons by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
#                8355                 1930                13896                 3092                 3916 
#
#10946 of 15245 genes with gene body methylation (0.718005903574943) have one or more CG sites which vary.
#1095 of 18096 genes without gene body methylation (0.0605106100795756) have one or more CG sites which vary.
#885 of 14594 unmethylated genes (0.0606413594627929) have one or more CG sites which vary.
#280 of 1930 transposons with gene body methylation (0.145077720207254) have one or more CG sites which vary.
#713 of 29259 transposons without gene body methylation (0.0243685703544209) have one or more CG sites which vary.
#99 of 3916 unmethylated transposons (0.0252808988764045) have one or more CG sites which vary
#
# Using more sensitive detection of variable sites (majority consensus in parental lines, any progeny line showing repwise variation c.f. parental consensus:
# Counts of genes by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
#                 347                15264                 1793                 1333                14604 
# Counts of transposons by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
#               8355                 2001                13986                 2877                 3970 
#
#12971 of 15264 genes with gene body methylation (0.849777253668763) have one or more CG sites which vary.
#1484 of 18077 genes without gene body methylation (0.0820932676882226) have one or more CG sites which vary.
#1213 of 14604 unmethylated genes (0.0830594357710216) have one or more CG sites which vary.
#563 of 2001 transposons with gene body methylation (0.28135932033983) have one or more CG sites which vary.
#1237 of 29188 transposons without gene body methylation (0.0423804303138276) have one or more CG sites which vary.
#206 of 3970 unmethylated transposons (0.0518891687657431) have one or more CG sites which vary.

#

# Analysis of distribution of CG sites in Schmitz data with Fishers p<0.05, Binomial p<0.005:
# 1085545 of 2802193 ( 0.387391232509681 )  CG  sites lie within annotated non-heterochromatic genes.
#113905 of 592062 ( 0.19238694596174 ) 'all M' sites lie within annotated non-heterochromatic genes.
#412585 of 898177 ( 0.45935823339943 ) 'all U' sites lie within annotated non-heterochromatic genes.
#35654 of 40815 ( 0.873551390420189 ) variable sites lie within annotated non-heterochromatic genes.
#15912 of 17078 ( 0.93172502634969 ) M->U transition sites lie within annotated non-heterochromatic genes.
#12647 of 15108 ( 0.83710616891713 ) U->M transition sites lie within annotated non-heterochromatic genes.
# How many CG, variable and invariable sites lie within annotated transposons?
#516901 of 2802193 ( 0.184463025922911 )  CG  sites lie within annotated transposons.
#350754 of 592062 ( 0.592427820059386 ) 'all M' sites lie within annotated transposons.
#24058 of 898177 ( 0.0267853663587467 ) 'all U' sites lie within annotated transposons.
#1644 of 40815 ( 0.040279309077545 ) variable sites lie within annotated transposons.
#571 of 17078 ( 0.0334348284342429 ) M->U transition sites lie within annotated transposons.
#515 of 15108 ( 0.0340879004500927 ) U->M transition sites lie within annotated transposons.
# How many CG, variable and invariable sites lie within annotated exons?
#907817 of 2802193 ( 0.323966621856525 )  CG  sites lie within annotated exons.
#91512 of 592062 ( 0.154564893541555 ) 'all M' sites lie within annotated exons.
#335957 of 898177 ( 0.374043200839033 ) 'all U' sites lie within annotated exons.
#28290 of 40815 ( 0.693127526644616 ) variable sites lie within annotated exons.
#12551 of 17078 ( 0.734922122028341 ) M->U transition sites lie within annotated exons.
#10047 of 15108 ( 0.665011914217633 ) U->M transition sites lie within annotated exons.
# How many CG, variable and invariable sites lie within annotated introns?
#418892 of 2802193 ( 0.149487205199642 )  CG  sites lie within annotated introns.
#59774 of 592062 ( 0.100959021183592 ) 'all M' sites lie within annotated introns.
#150464 of 898177 ( 0.167521546421251 ) 'all U' sites lie within annotated introns.
#17577 of 40815 ( 0.430650496141125 ) variable sites lie within annotated introns.
#8356 of 17078 ( 0.489284459538588 ) M->U transition sites lie within annotated introns.
#5898 of 15108 ( 0.390389197776013 ) U->M transition sites lie within annotated introns.
# How many CG, variable and invariable sites lie within annotated 5'UTRs?
#129019 of 2802193 ( 0.0460421534134159 )  CG  sites lie within annotated 5'UTRs.
#1287 of 592062 ( 0.0021737588293118 ) 'all M' sites lie within annotated 5'UTRs.
#76973 of 898177 ( 0.0856991439326547 ) 'all U' sites lie within annotated 5'UTRs.
#421 of 40815 ( 0.010314835232145 ) variable sites lie within annotated 5'UTRs.
#130 of 17078 ( 0.00761213256821642 ) M->U transition sites lie within annotated 5'UTRs.
#180 of 15108 ( 0.0119142176330421 ) U->M transition sites lie within annotated 5'UTRs.
# How many CG, variable and invariable sites lie within annotated 3'UTRs?
#90969 of 2802193 ( 0.0324635026923556 )  CG  sites lie within annotated 3'UTRs.
#5337 of 592062 ( 0.00901425864183143 ) 'all M' sites lie within annotated 3'UTRs.
#36003 of 898177 ( 0.0400845267692225 ) 'all U' sites lie within annotated 3'UTRs.
#2089 of 40815 ( 0.0511821634203112 ) variable sites lie within annotated 3'UTRs.
#745 of 17078 ( 0.043623375102471 ) M->U transition sites lie within annotated 3'UTRs.
#868 of 15108 ( 0.0574530050304474 ) U->M transition sites lie within annotated 3'UTRs.
#
# Schmitz Fishers p<0.05, Binomial p<0.005, more sensitive method for identifying variable sites:
#1085658 of 2802193 ( 0.387431558068984 )  CG  sites lie within annotated non-heterochromatic genes.
#113905 of 592062 ( 0.19238694596174 ) 'all M' sites lie within annotated non-heterochromatic genes.
#412585 of 898177 ( 0.45935823339943 ) 'all U' sites lie within annotated non-heterochromatic genes.
#60345 of 70768 ( 0.852715916798553 ) variable sites lie within annotated non-heterochromatic genes.
#32607 of 35852 ( 0.90948901037599 ) M->U transition sites lie within annotated non-heterochromatic genes.
#27738 of 34916 ( 0.773680687269887 ) U->M transition sites lie within annotated non-heterochromatic genes.
# How many CG, variable and invariable sites lie within annotated transposons?
#516901 of 2802193 ( 0.184463025922911 )  CG  sites lie within annotated transposons.
#350754 of 592062 ( 0.592427820059386 ) 'all M' sites lie within annotated transposons.
#24058 of 898177 ( 0.0267853663587467 ) 'all U' sites lie within annotated transposons.
#3496 of 70768 ( 0.0494008591453764 ) variable sites lie within annotated transposons.
#1476 of 35852 ( 0.0411692513667299 ) M->U transition sites lie within annotated transposons.
#2020 of 34916 ( 0.0578531332340474 ) U->M transition sites lie within annotated transposons.
# How many CG, variable and invariable sites lie within annotated exons?
#907919 of 2802193 ( 0.324003021918904 )  CG  sites lie within annotated exons.
#91512 of 592062 ( 0.154564893541555 ) 'all M' sites lie within annotated exons.
#335957 of 898177 ( 0.374043200839033 ) 'all U' sites lie within annotated exons.
#49132 of 70768 ( 0.694268595975582 ) variable sites lie within annotated exons.
#26445 of 35852 ( 0.737615753653911 ) M->U transition sites lie within annotated exons.
#22687 of 34916 ( 0.649759422614274 ) U->M transition sites lie within annotated exons.
# How many CG, variable and invariable sites lie within annotated introns?
#418954 of 2802193 ( 0.149509330727755 )  CG  sites lie within annotated introns.
#59774 of 592062 ( 0.100959021183592 ) 'all M' sites lie within annotated introns.
#150464 of 898177 ( 0.167521546421251 ) 'all U' sites lie within annotated introns.
#28286 of 70768 ( 0.399700429572688 ) variable sites lie within annotated introns.
#15977 of 35852 ( 0.445637621332143 ) M->U transition sites lie within annotated introns.
#12309 of 34916 ( 0.352531790583114 ) U->M transition sites lie within annotated introns.
# How many CG, variable and invariable sites lie within annotated 5'UTRs?
#129031 of 2802193 ( 0.0460464357736958 )  CG  sites lie within annotated 5'UTRs.
#1287 of 592062 ( 0.0021737588293118 ) 'all M' sites lie within annotated 5'UTRs.
#76973 of 898177 ( 0.0856991439326547 ) 'all U' sites lie within annotated 5'UTRs.
#721 of 70768 ( 0.0101882206647072 ) variable sites lie within annotated 5'UTRs.
#335 of 35852 ( 0.00934396965301796 ) M->U transition sites lie within annotated 5'UTRs.
#386 of 34916 ( 0.0110551036773972 ) U->M transition sites lie within annotated 5'UTRs.
# How many CG, variable and invariable sites lie within annotated 3'UTRs?
#90984 of 2802193 ( 0.0324688556427056 )  CG  sites lie within annotated 3'UTRs.
#5337 of 592062 ( 0.00901425864183143 ) 'all M' sites lie within annotated 3'UTRs.
#36003 of 898177 ( 0.0400845267692225 ) 'all U' sites lie within annotated 3'UTRs.
#3788 of 70768 ( 0.0535270178611802 ) variable sites lie within annotated 3'UTRs.
#1812 of 35852 ( 0.0505411134664733 ) M->U transition sites lie within annotated 3'UTRs.
#1976 of 34916 ( 0.056592965975484 ) U->M transition sites lie within annotated 3'UTRs.

# Multiple variable sites in transposons are more likely to represent an overall transisiton in one direction, than are those in genes (see *_Mstats_transitions_per_* plots)
#
# Becker data on genes and transposons relative to 'heterochromatic' methylation variance threshold:
#33341 annotated genes
#31296 annotated genes with average_methylation<0.6
#31234 annotated genes with average_methylation<0.6 and dispersion index >0.18968726371385 where average_methylation>0.2 or average_methylation<0.2
#27681 annotated genes with methylation dispersion index>0.18968726371385
#27298 annotated genes with average methylation<0.6 and dispersion index>0.18968726371385
#
#31189 annotated transposons
#15226 annotated transposons with average_methylation<0.6
#15185 annotated transposons with average_methylation<0.6 and dispersion index >0.18968726371385 where average_methylation>0.2 or average_methylation<0.2
#17682 annotated transposons with methylation dispersion index>0.18968726371385
#14867 annotated transposons with average methylation<0.6 and dispersion index>0.18968726371385
#
# 2340 genes are assigned as heterochromatic.
#
# More refined classification:
#
#8174 of 33341 genes have one or more sites which vary with concordance among replicates (0.245163612369155).
#241 of 2340 heterochromatic genes have one or more sites which vary with concordance among replicates (0.102991452991453).
#150 of 659 heterochromatic genes with average methylation <=0.8 have one or more sites which vary with concordance among replicates (0.227617602427921).
#91 of 1914 heterochromatic genes with average methylation >0.8 have one or more sites which vary with concordance among replicates (0.0475444096133751).
#590 of 31189 transposons have one or more sites which vary with concordance among replicates (0.0189169258392382).
#
#Counts of genes by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated
#                 347                15187                 1038                 6542                10226 
#Counts of transposons by classification:
#             CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
#                8355                  407                 5929                15865                  630 
#7499 of 15188 genes with gene body methylation (0.493745061890967) have one or more CG sites which vary.
#675 of 18154 genes without gene body methylation (0.0371818882890823) have one or more CG sites which vary.
#572 of 10227 unmethylated genes (0.0559303803656986) have one or more CG sites which vary.
#117 of 410 transposons with gene body methylation (0.285365853658537) have one or more CG sites which vary.
#473 of 30782 transposons without gene body methylation (0.0153661230589305) have one or more CG sites which vary.
#34 of 633 unmethylated transposons (0.0537124802527646) have one or more CG sites which vary.
#
#

# No. of variable sites in a GBM gene is strongly corellated with gene length (r=0.590166 in Schmitz data), and less strongly with number of exons (r=0.421135 in Schmitz data)
#
#
#
# Analysis w.r.t. Chromatin States (Sequeira-Mendes et al, 2014):
# There is very strong enrichment for genome segments annotated as chromatin state 3,6,7 among sites which vary in one or more lines at generation 30, with rep-wise concordance:
#
# Schmitz data on proportions of CG sites in segments of each chromatin state:
#        CG sites        all_M        all_U   variable  variableMU  variableUM
# [1,] 0.15617353 0.0006463793 2.692862e-01 0.01196073 0.006786663 0.018941062
# [2,] 0.10903776 0.0047022718 1.586336e-01 0.01522274 0.006712895 0.023786451
# [3,] 0.08624345 0.0416960406 8.739324e-02 0.26674982 0.276261434 0.264293895
# [4,] 0.11038159 0.0361258261 1.246767e-01 0.03732131 0.016007672 0.050392036
# [5,] 0.14611876 0.0136325599 1.882913e-01 0.06114682 0.021318973 0.106598538
# [6,] 0.09012404 0.0226104567 1.123065e-01 0.16735425 0.132708764 0.209761255
# [7,] 0.08270483 0.1003169639 5.091104e-02 0.39508779 0.496680437 0.302616510
# [8,] 0.07733529 0.2418264151 8.443355e-03 0.03457098 0.028695780 0.022464981
# [9,] 0.14188075 0.5384430866 5.800675e-05 0.01058556 0.014827383 0.001145274
#
#
#
# There is some corellation between numbers of sites changing in each direction among individual chromatin state segments
# Schmitz data:
#Corellation between numbers of M->U and U->M sites per chromatin state segment (segments containing variable sites only): 0.198037252039741
#Corellation between numbers of M->U and U->M sites per chromatin state segment (all segments, including sites with no variation)): 0.3644361890509
#
#
# DNase-hypersensitivity loci are strongly enriched for Unmethylated sites, and strongly depleted for methylated and variable sites 
# Schmitz data:
#460346 of 2802193 CG sites (0.164280618786786) lie within DHSs.
#330 of 31269 variable CG sites (0.0105535834212799) lie within DHSs.
#20659 of 546168 all M CG sites (0.0378253577653762) lie within DHSs.
#170877 of 741378 all U CG sites (0.230485663183963) lie within DHSs
#
# 89% of CG sites in DHSs with a rep-consistent call are U
# 11% are M
# 0.2% are variable
#
# This compares with these figures for the genome as a whole:
# 56% are U
# 41% are M
# 2% are variable




	
#### This part parses command line options or sets defaults

if(!require(optparse)){
	install.packages("optparse")
	library(optparse)
}

library("optparse")

option_list = list(
	make_option(c("-p", "--project"), action="store", default=NA, type='character',
              help="project code for location of bs_sequel output"),
	make_option(c("-c", "--context"), action="store", default=NA, type='character',
              help="methylation context (CG, CHG or CHH)"),
	make_option(c("-a", "--action"), action="store", default="all", type='character',
              help="actions to perform (load, plot, call, merge, analyse or all)"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Should the program print extra stuff out? [default %default]"),
	make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Make the program not be verbose.")
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$verbose) {
    # show the user which options were chosen
	cat("Project: ")
    cat(opt$project)
    cat("\n")
    
    cat("Context: ")
    cat(opt$context)
    cat("\n")

    cat("Action: ")
    cat(opt$action)
    cat("\n")
}

# project_id and meth_context together are used for finding raw data and naming output files
project_id=NULL
if (is.na(opt$project)) {
	# No project defined - set a default
	#project_id="test"
	project_id="SRA035939"  # Schmitz et al, 2011
	#project_id="PRJEB2678"  # Becker et al, 2011
} else {
	project_id=opt$project
}
	
meth_context=NULL
if (is.na(opt$context)) {
	# No context defined - set a default
	meth_context="CG"
	#meth_context="CHG"
	#meth_context="CHH"
} else {
	meth_context=opt$context
}

action=NULL
if (is.na(opt$action)) {
	# No action defined - set a default
	action="all"
	#action="load"
	#action="plot"
	#action="call"
	#action="merge"
	#action="analyse"
} else {
	action=opt$action
}


# Platform-specific stuff:
#  Base of path to project
#  Where to find bioconductor (unresolved on cluster at the moment)
pathroot = ""
if (.Platform$OS.type=="windows") {
	pathroot="D:"
	source("http://bioconductor.org/biocLite.R")
} else if (.Platform$OS.type=="unix") {
	if(substr(getwd(),1,20)=="/nbi/Research-Groups") {
		# assume we are running on cluster
		pathroot="/nbi/Research-Groups/JIC/Daniel-Zilberman"
	} else {
		# assume we are running in Virtualbox with shared folder /media/sf_D_DRIVE
		pathroot="/media/sf_D_DRIVE"
		source("http://bioconductor.org/biocLite.R")
	}
}

biocLite(c("GenomicRanges", "BSgenome", "BSgenome.Athaliana.TAIR.TAIR9", "MethylSeekR", "karyoploteR"))

#### This part installs packages so will be slow the first time per platform

if(!require(data.table)){
	install.packages("data.table")
	library(data.table)
}
if(!require(stringr)){
	install.packages("stringr")
	library(stringr)
}
if(!require(ggplot2)){
	install.packages("ggplot2")
	library(ggplot2)
}
if(!require(mixtools)){
	install.packages("mixtools")
	library(mixtools)
}
if(!require(rootSolve)){
	install.packages("rootSolve")
	library(rootSolve)
}
install.packages("ggpubr")

	
#### This part sets up libraries and contextual information and metadata about the project

library(data.table)
library(stringr)
library(ggplot2)
library(plyr)
library(ggpubr)


# Source directory containing alignments from bs_sequel
source_dir = paste0(pathroot,"/Projects/Jay-SMP_variation/3-alignments/",project_id,"/")

# Working directory where outputs will be directed to
setwd(dir = paste0(pathroot,"/Projects/Jay-SMP_variation/5-analysis/",project_id,"/"))

reference_fasta = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta")  # TAIR10 genome assembly This never gets used at the moment
reference_CG_sites = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.CG_sites.tsv")   # We prepared this earlier using a perl script to parse the TAIR10 genome assembly
reference_gff = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/Araport11_GFF3_genes_transposons.201606.gff")  # AraPort11 annotation
reference_exons = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-exon.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
reference_introns = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-intron.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
reference_5UTR = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-five_prime_UTR.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
reference_3UTR = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-three_prime_UTR.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
### SNP_loci.txt - we need to import a list of all known SNPs between lines considered for analysis, so SNP loci can inform site exclusion criteria
reference_chromatin_states = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Chromatin_states/Sequeira-Mendes_et_al_2014_tpc124578_SupplementalIDS2.txt")  # genomic ranges assigning each segment of the nuclear genome to one of 9 chromatin states
reference_DHS_loci = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/DHS_scores/TAIR10_DHSs.gff")  # DNase-hypersensitivity sites (indicative of open chromatin)

# Nucleosomes position data from http://plantdhs.org/Download (Jiming Jiang lab) Originally described in Wu Y.F. Zhang W.L. Jiang J.M. Genome-wide nucleosome positioning is orchestrated by genomic regions associated with DNase I hypersensitivity in rice PLoS Genet.  2014 10 e1004378
# For nucleosome positioning data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Nucleosomes/Ath_leaf_NPS.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Nucleosomes/Ath_leaf_NPS.RDS")
#reference_nucleosomes = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Nucleosomes/Ath_leaf_NPS.RDS")  

# Nucleosomes position data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2807196
# Lyons DB, Zilberman D. DDM1 and Lsh remodelers allow methylation of DNA wrapped in nucleosomes. Elife 2017 Nov 15;6. PMID: 29140247
# GSE96994 wt.rep1.mnase_seq results
# GSM2807196_wt.rep1.nucleosomes.bed is all nucleosome calls
# GSM2807196_wt.group_1.bed is well-positioned nucleosome calls

#reference_nucleosomes = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Nucleosomes/GSM2807196_wt.rep1.nucleosomes.bed")
reference_nucleosomes = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Nucleosomes/GSM2807196_wt.group_1.bed")

# H3K9me2 normalised WT data from Dave Lyons
# For H3K9me2 data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/dk9.jacobsen.normalized.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/dk9.jacobsen.normalized.RDS")
reference_H3K9me2 = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/dk9.jacobsen.normalized.RDS")  

# H3K9me1 Plant DNase I hypersensitivity data. Downloaded from http://plantdhs.org/Download (Jiming Jiang lab) on 05/12/2017
# For H3K4me1 data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_buds_H3K4me1.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_buds_H3K4me1.RDS")
reference_H3K4me1 = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_buds_H3K4me1.RDS")

# H3K9me1 Plant DNase I hypersensitivity data. Downloaded from http://plantdhs.org/Download (Jiming Jiang lab) on 05/12/2017
# For H3K27me3 data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_leaf_H3K27me3.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_leaf_H3K27me3.RDS")
reference_H3K27me3 = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_leaf_H3K27me3.RDS")

# H2AW data from Jaemyung
reference_H2AW = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/h2aw.w50.gff")


# get a list of the samples from the subdirectories of the source directory
samples = list.dirs(path = source_dir, full.names = FALSE, recursive = FALSE)
no_samples=length(samples)

# Read the sample metadata in from a text file in the source directory
sample_metadata=read.delim(paste0(source_dir,project_id,"_metadata.tsv"), header=TRUE, sep="\t")
rownames(sample_metadata)=sample_metadata$Identifier

# Read the Schmitz et al, 2011 DMRs data, so this can be used in masking where needed
Schmitz_CG_DMRs = paste0(pathroot,"/Literature/Schmitz_et_al_2011_CG_DMRs.txt")
Schmitz_non_CG_DMRs = paste0(pathroot,"/Literature/Schmitz_et_al_2011_non-CG_DMRs.txt")


#### This part gets the raw alignment/basecall data in, makes some convenient reshapings of the data, and saves it for later before going on to the analyses

if ((action=="load") || (action=="all")) {

	#### This section loads the raw data from GFF files into a data frame
	#### It takes a while, so run it the first time for the project, then subsequently skip to the part after this loop where the data.table is loaded from a cached version

	# read each sample's 'all chromosomes' single-c GFF file
	# append the coverage data to a table with one row per site, and one 'C' column and one 'T' column per sample
	# this table can then be used to access data about sites across samples

	all_samples_ct=NA
	for(sample_no in 1:no_samples) {
		sample_name=strsplit(samples[sample_no],"_")[[1]][[1]]
		if(opt$verbose) {cat(paste0("Loading ",sample_name,"\n"))}
		x=read.table(file=paste0(source_dir,samples[sample_no],"/single-c/",samples[sample_no],".chr.all.single-c.",meth_context,".gff"))

		# parse the awkward "c=2;t=25*" format into two new columns for each sample
		y=do.call(rbind, str_split(do.call(rbind, str_split(x[,9], ';')), '='))
		x[,10]=as.numeric(y[1:nrow(x),2])
		x[,11]=y[(nrow(x)+1):(2*nrow(x)),2]
		x[,12]=as.numeric(ifelse(substr(x[,11],nchar(x[,11]),nchar(x[,11]))=="*",substr(x[,11],1,nchar(x[,11])-1),x[,11]))

		# merge the data for this sample with the accumulated data	
		colnames(x)[10]=paste0(sample_name,"_C")
		colnames(x)[12]=paste0(sample_name,"_T")
	
		if (sample_no==1) {
			all_samples_ct=x[c("V1","V4","V7",paste0(sample_name,"_C"),paste0(sample_name,"_T"))]
		} else {
			all_samples_ct=merge(x[c("V1","V4","V7",paste0(sample_name,"_C"),paste0(sample_name,"_T"))], all_samples_ct, by=c("V1","V4","V7"), all=TRUE)
		}
		#all_samples_ct[,sample_name]=as.character(all_samples_ct[,sample_name])
	
		rm(x)
		rm(y)
	}
	colnames(all_samples_ct)[1]="Chromosome"
	colnames(all_samples_ct)[2]="Locus"
	colnames(all_samples_ct)[3]="Strand"

	all_samples_ct = data.table(all_samples_ct,key=c("Chromosome", "Locus","Strand"))

	#Add a column to identify the number of samples with missing data
	all_samples_ct$missing_count <- apply(all_samples_ct, 1, function(x) sum(is.na(x))/2)

	# Dump the data table to file for a convenience cache
	saveRDS(all_samples_ct, file=paste0(project_id,"_",meth_context,"_all_samples_ct.rds"))


	#### This section develops some summary objects
	#### Coverage data from the data.table are summarised by sample to a data.frame sample_info
	#### Data from each individual sample are appended to a 'long' data.frame coverage_data, for all samples, with one row per site per sample

	sample_info=data.frame(Cov_C=integer(no_samples), Cov_T=integer(no_samples), CT_ratio=numeric(no_samples), ISNA_nonC=integer(no_samples), ISNA_C=integer(no_samples), Cprop_cutoff=numeric(no_samples))
	rownames(sample_info)=substr(samples,1,9)
	# quantile cutoff for C/T ratio in Chloroplast meth_context sites
	Cprop_quantile=0.99

	coverage_data = data.frame("Sample"=as.factor(NULL),"Chromosome"=as.factor(NULL),"Locus"=as.integer(NULL),"Strand"=as.factor(NULL),"Cov_C"=as.integer(NULL),"Cov_T"=as.integer(NULL))

	for(sample_name in rownames(sample_info)) {
		if(opt$verbose) {cat(paste0("Summarising: ",sample_name,"\n"))}
	
		# Summarise the data for the sample
		sample_info[sample_name,"Cov_C"]=sum(all_samples_ct[,paste0(sample_name,"_C"),with=FALSE],na.rm=TRUE)
		sample_info[sample_name,"Cov_T"]=sum(all_samples_ct[,paste0(sample_name,"_T"),with=FALSE],na.rm=TRUE)
		sample_info[sample_name,"CT_ratio"]=sample_info[sample_name,"Cov_C"]/sample_info[sample_name,"Cov_T"]
		sample_info[sample_name,"ISNA_nonC"]=sum(is.na(all_samples_ct[all_samples_ct$Chromosome!="CHRC",paste0(sample_name,"_C"),with=FALSE]))
		sample_info[sample_name,"ISNA_C"]=sum(is.na(all_samples_ct[all_samples_ct$Chromosome=="CHRC",paste0(sample_name,"_C"),with=FALSE]))
		
		
		# In the methylation calling stage later, the C/(C+T) ratio in the chloroplast is used as a proxy estimate of the combined technical errors introduced by sequencing and bisulphite conversion.
		# A threshold C/(C+T) ratio is calculated here for each sample's chloroplast reads, against which to test all the other sites for whether they contain significantly more C reads than an unmethylated site, and/or significantly fewer C reads than a methylated site.
		# Two different methods of calculating the C/(C+T) threshold are proposed.
		# The first method uses the .99 quantile of chloroplast C/(C+T), on the basis that this should be a reasonable upper bound for unmethylated sites. One drawback of this method is that, for low coverage samples, a larger number of sites will fail to be called due to not being significantly more methylated than a U site while neither being significantly less methylkated than an M site.  Such sites ultimately receive an I classification.  A second drawback is that it can potentially increase the number of sites called as U which are in reality methylated at a low level.
		# The second method uses the mean of chloroplast C/(C+T), so makes the test for a U site much more stringent, but at the same time, allows a call to be made for more sites with lower coverage.  This might be justified in cases with low coverage, but where replicates are available to increase the likelihood that sites called incorrectly due to low coverage will be contradicted by the call in their replicate.

		# Method 1
		# Find C/(C+T) cutoff thresholds for the sample from 99% quantile of C/(C+T) distribution on Chloroplast	
		sample_info[sample_name,"Cprop_cutoff"]=quantile(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="CHRC"),paste0(sample_name,"_C"),with=FALSE])[,1]/(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="CHRC"),paste0(sample_name,"_C"),with=FALSE])[,1]+as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="CHRC"),paste0(sample_name,"_T"),with=FALSE])[,1]),Cprop_quantile,na.rm=TRUE)
		
		# Method 2
		# Find C/(C+T) ratio for chloroplast as a whole (estimation of conversion failure rate + sequencing error)
		#sample_info[sample_name,"Cprop_cutoff"]=sum(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="CHRC"),paste0(sample_name,"_C"),with=FALSE])[,1],na.rm=TRUE)/(sum(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="CHRC"),paste0(sample_name,"_C"),with=FALSE])[,1],na.rm=TRUE)+sum(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="CHRC"),paste0(sample_name,"_T"),with=FALSE])[,1],na.rm=TRUE))
		
		# Add the sample coverage data to the long table
		sample_coverage_data=as.data.frame(all_samples_ct[,c("Chromosome","Locus","Strand","Cov_C"=paste0(sample_name,"_C"),"Cov_T"=paste0(sample_name,"_T")),with=FALSE])
		sample_coverage_data=cbind(rep(sample_name,nrow(sample_coverage_data)),sample_coverage_data)
		colnames(sample_coverage_data)=c("Sample","Chromosome","Locus","Strand","Cov_C","Cov_T")	
		coverage_data=rbind(coverage_data,sample_coverage_data)
		rm(sample_coverage_data)
	}

	# Dump the sample info to a file for convenience
	write.table(sample_info,file=paste0(project_id,"_",meth_context,"_sample_info.tsv"),sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)
	saveRDS(sample_info, file=paste0(project_id,"_",meth_context,"_sample_info.rds"))

	# Dump the coverage data table to file for a convenience cache
	saveRDS(coverage_data, file=paste0(project_id,"_",meth_context,"_coverage_data.rds"))

} # end action=="load"

if ((action=="plot") || (action=="all")) {
	
	if (action=="plot") {
		# Reload the data objects developed in earlier actions
		all_samples_ct <- readRDS(paste0(project_id,"_",meth_context,"_all_samples_ct.rds"))
		sample_info <- readRDS(paste0(project_id,"_",meth_context,"_sample_info.rds"))
		coverage_data <- data.table(readRDS(paste0(project_id,"_",meth_context,"_coverage_data.rds")))
	}
	
	#### Plot useful stuff for QC of the BS data
	
	# Plot summary of missing sites data
	pdf(paste0(project_id,"_",meth_context,"_sites_missing_data_summary.pdf"))
	print(ggplot(all_samples_ct[all_samples_ct$missing_count>0,], aes(x=missing_count, colour=Chromosome)) + geom_freqpoly(binwidth = 1) + xlab(paste("Number of samples missing data for a",meth_context,"site")) +ylab("Frequency of sites") + annotate("text", x = 10, y = 1e5, label = paste(nrow(all_samples_ct[all_samples_ct$missing_count==0,]),meth_context,"sites have no missing data")))
	dev.off() 

	# Plot coverage distributions

	# These limits are applied to both C and T coverage at each site, to allow exclusion of lowly or highly covered sites that mess up the plots
	coverage_min=0
	coverage_max=50

	# Plot C/(C+T) distribution by chromosome by sample, and annotate with quantile cutoff value
	# Coverage limits not relevant here
	pdf(paste0(project_id,"_",meth_context,"_site_C_proportion_by_chrom_by_sample.pdf"))
	#print(ggplot(coverage_data[(!is.na(coverage_data$Cov_C)) & (coverage_data$Cov_C>=coverage_min) & (coverage_data$Cov_T>=coverage_min) & (coverage_data$Cov_C<=coverage_max) & (coverage_data$Cov_T<=coverage_max),], aes(x=Cov_C/(Cov_C+Cov_T),colour=Chromosome)) + geom_freqpoly(binwidth = .01) + facet_wrap(~Sample) + xlab(paste(meth_context,"site C/(C+T)")) +ylab("Number of sites") +labs(fill = "Chromosome"))
	# This part is very slow due to the merging and conversion from factor to numeric of the cutoff value. The alternative without the quantile values is commented out above
	print(ggplot(merge(coverage_data[(!is.na(coverage_data$Cov_C)),],cbind("Sample"=rownames(sample_info),Cratio_cutoff=sample_info$Cprop_cutoff)), aes(x=Cov_C/(Cov_C+Cov_T),colour=Chromosome)) + geom_freqpoly(binwidth = .01) + facet_wrap(~Sample) + xlab(paste(meth_context,"site C/(C+T)")) +ylab("Number of sites") +labs(fill = "Chromosome") + geom_vline(aes(xintercept = as.numeric(as.character(Cratio_cutoff)))))
	dev.off()

	# Plot C/(C+T) distribution for chloroplast only, by sample, and annotate with quantile cutoff value
	# Don't apply coverage limits to chloroplast data especially, as coverage will be much higher than other chromosomes
	pdf(paste0(project_id,"_",meth_context,"_site_C_proportion_chloroplast_by_sample.pdf"))
	print(ggplot(merge(coverage_data[(!is.na(coverage_data$Cov_C)) & (coverage_data$Chromosome=="CHRC"),],cbind("Sample"=rownames(sample_info),Cratio_cutoff=sample_info$Cprop_cutoff)), aes(x=Cov_C/(Cov_C+Cov_T),colour=Chromosome)) + geom_freqpoly(binwidth = .01) + facet_wrap(~Sample) + xlab(paste(meth_context,"site C/(C+T)")) +ylab("Number of sites") +labs(fill = "Chromosome") + geom_vline(aes(xintercept = as.numeric(as.character(Cratio_cutoff)))))
	dev.off()

	# plot coverage frequency polys for each sample by strand
	# Coverage limits are applied to make x axis reasonable
	pdf(paste0(project_id,"_",meth_context,"_site_coverage_by_strand_by_sample.pdf"))
	print(ggplot(coverage_data[(!is.na(coverage_data$Cov_C)) & (coverage_data$Cov_C>=coverage_min) & (coverage_data$Cov_T>=coverage_min) & (coverage_data$Cov_C<=coverage_max) & (coverage_data$Cov_T<=coverage_max),], aes(x=Cov_C+Cov_T,colour=Strand)) + geom_freqpoly(binwidth = 1) + facet_wrap(~Sample) + xlab(paste(meth_context,"site coverage")) +ylab("Number of sites") +labs(fill = "Chromosome"))
	dev.off()
	
	# plot coverage frequency polys for each sample by chromosome
	# Coverage limits are applied to make x axis reasonable
	pdf(paste0(project_id,"_",meth_context,"_site_coverage_by_chrom_by_sample.pdf"))
	print(ggplot(coverage_data[(!is.na(coverage_data$Cov_C)) & (coverage_data$Cov_C>=coverage_min) & (coverage_data$Cov_T>=coverage_min) & (coverage_data$Cov_C<=coverage_max) & (coverage_data$Cov_T<=coverage_max),], aes(x=Cov_C+Cov_T,colour=Chromosome)) + geom_freqpoly(binwidth = 1) + facet_wrap(~Sample) + xlab(paste(meth_context,"site coverage")) +ylab("Number of sites") +labs(fill = "Chromosome"))
	dev.off()

	# plot cumulative frequency version
	# Coverage limits are applied to make x axis reasonable
	pdf(paste0(project_id,"_",meth_context,"_site_coverage_by_chrom_by_sample_CumFreq.pdf"))
	print(ggplot(coverage_data[(!is.na(coverage_data$Cov_C)) & (coverage_data$Cov_C>=coverage_min) & (coverage_data$Cov_T>=coverage_min) & (coverage_data$Cov_C<=coverage_max) & (coverage_data$Cov_T<=coverage_max),], aes(x=Cov_C+Cov_T,colour=Chromosome)) + stat_ecdf(geom = "step") + facet_wrap(~Sample) + xlab(paste(meth_context,"site coverage")) +ylab("Cumulative frequency") + labs(fill = "Chromosome"))
	dev.off()
} # end action=="plot"



#### This part classifies each site as one of: Unmethylated; Partially methylated; Methylated; Indeterminate with p-value

if ((action=="call") || (action=="all")) {
	
	if (action=="call") {
		#### Reload the data objects developed in earlier actions
		# all_samples_ct <- readRDS(paste0(project_id,"_",meth_context,"_all_samples_ct.rds")) # Commented out as not needed anymore?
		sample_info <- readRDS(paste0(project_id,"_",meth_context,"_sample_info.rds"))
		coverage_data <- data.table(readRDS(paste0(project_id,"_",meth_context,"_coverage_data.rds")))
	} else {
		coverage_data=data.table(coverage_data)
	}
	
		if (meth_context=="CG") {
		# This part works on pairs of adjacent CG sites (CG dinucleotides in the annotation) to generate a composite methylation call per CG site pair

		# The basis of this process is that in a CG dinucleotide, where adjacent sites both have an unambiguous call (M/U), in the overwhelming majority of cases (>99.9% for both Schmitz and Becker data sets at least) the two calls are found to agree.  We use this information to make inferences about the methylation status of a CG dinucleotide on the basis of combining the #C and #T data from both members of each pair of nucleotides.  Where one nucleotide from a pair has no coverage data, we rely on the coverage at the other nucleotide in the pair (if there is coverage). Where neither member of the pair has coverage, we generate NA values for Cov_C and Cov_T, so the final table will still contain entries for all CG sites in the reference.

		# Record the data on the sample/locus_pair in a new table

		# Read in the list of known CG dinucleotides in the reference, we prepared earlier in perl
		reference_CG_site_pair = data.table(read.table(reference_CG_sites, header=FALSE, sep='\t'))
		colnames(reference_CG_site_pair)=c("Chromosome","Locus","endLocus")
		reference_CG_site_pair$Chromosome = toupper(reference_CG_site_pair$Chromosome)
		#reference_CG_site_pair$Locus=as.integer(levels(reference_CG_site_pair$Locus)[reference_CG_site_pair$Locus])
		#reference_CG_site_pair$endLocus=as.integer(levels(reference_CG_site_pair$endLocus)[reference_CG_site_pair$endLocus])

		sample_no=0
		for(sample_name in rownames(sample_info)) {
			sample_no = sample_no+1
			if(opt$verbose) {cat(paste0("Accumulating #C and #T from CG dinucleotides for sample: ",sample_name,"\n"))}

			# Merge the #C and #T counts for the sample at the first locus of each CG site pair
			sample_combined_CG_coverage_data = merge(reference_CG_site_pair,coverage_data[coverage_data$Sample==sample_name,c("Chromosome","Locus","Strand","Cov_C","Cov_T")], by=c("Chromosome", "Locus"), all.x = TRUE)
			# Rename the merged columns
			colnames(sample_combined_CG_coverage_data) = c("Chromosome", "startLocus", "Locus", "Strand1", "Cov_C1", "Cov_T1")
			# Merge the #C and #T counts for the sample at the second locus of each CG site pair
			sample_combined_CG_coverage_data = merge(sample_combined_CG_coverage_data,coverage_data[coverage_data$Sample==sample_name,c("Chromosome","Locus","Strand","Cov_C","Cov_T")], by=c("Chromosome", "Locus"), all.x = TRUE)
			colnames(sample_combined_CG_coverage_data) = c("Chromosome", "endLocus", "Locus", "Strand1", "Cov_C1", "Cov_T1", "Strand2", "Cov_C2", "Cov_T2")
			sample_combined_CG_coverage_data=data.table(sample_combined_CG_coverage_data)

			# Summarise the #C and #T counts at the two adjacent loci, and make a table of values for the combined locus
			sample_combined_CG_coverage_data$Cov_C=ifelse(is.na(sample_combined_CG_coverage_data$Cov_C1),0,sample_combined_CG_coverage_data$Cov_C1)+ifelse(is.na(sample_combined_CG_coverage_data$Cov_C2),0,sample_combined_CG_coverage_data$Cov_C2)
			sample_combined_CG_coverage_data$Cov_T=ifelse(is.na(sample_combined_CG_coverage_data$Cov_T1),0,sample_combined_CG_coverage_data$Cov_T1)+ifelse(is.na(sample_combined_CG_coverage_data$Cov_T2),0,sample_combined_CG_coverage_data$Cov_T2)

			# Put the sample name and strand in for compatibility with other contexts
			sample_combined_CG_coverage_data$Sample=sample_name
			sample_combined_CG_coverage_data$Strand="+"

 			# Accumulate the results for this sample in the overall table
			if (sample_no == 1) {
				combined_CG_coverage_data=sample_combined_CG_coverage_data[,c("Sample","Chromosome","Locus","Strand","Cov_C","Cov_T")]
			} else {
				combined_CG_coverage_data=rbind(combined_CG_coverage_data,sample_combined_CG_coverage_data[,c("Sample","Chromosome","Locus","Strand","Cov_C","Cov_T")])
			}
			# Clean up sample cache table
			rm(sample_combined_CG_coverage_data)

		} # end for each sample

		# Turn combined_CG_coverage_data from matrix to data frame
		combined_CG_coverage_data=data.frame(combined_CG_coverage_data)

		# Replace the individual site coverage data with the combined data for CG dinucleotides
		individual_CG_coverage_data=coverage_data
		coverage_data=combined_CG_coverage_data
		rm(combined_CG_coverage_data)

		# Dump the pair #C and #T counts to file for a convenience cache
		saveRDS(coverage_data, file=paste0(project_id,"_",meth_context,"_coverage_data.rds"))
		saveRDS(individual_CG_coverage_data, file=paste0(project_id,"_",meth_context,"_individual_CG_coverage_data.rds"))
		coverage_data=data.table(coverage_data)

	} # end meth_context == "CG"


	# Now we have made a special case for CG context by merging adjacent paired loci, we can treat each context the same to do the methylation calling
	
	
	#Values of status are "M"ethylated, "P"artially methylated, "U"nmethylated, "I"ndeterminate
		
	#The algorithm to determine methylation status:
	#
	# Assuming approximately Gaussian distribution of coverage by site (it tends to have a long tail though), if the coverage exceed 5 sigma, then discard the site from further consideration (status is "I"). Otherwise:
	#
	# Use Fisher's exact test to test whether the #C and #T at the site differs significantly from that expected from a completely unmethylated site with a level of sequencing errors and conversion inefficiency similar to that of the Chloroplast (sample_Cprop_cutoff).  Similarly, test whether #C and #T differ significantly from a fully methylated site with a similar level of errors.  If the site cannot be significantly distinguished from either of these two states (coverage too low/intermediate #C/(#C+#T)), then status is "I".  Otherwise:
	#
	# Use binomial test to see whether the site's #C and #T can be called significantly different from that expected of an unmethylated site with errors (sample_Cprop_cutoff). If not, then status is "U". Otherwise:
	#
	# If #C represents more than a given threshold proportion of the reads, then status is "M", else, "P".  The threshold may differ between methylation contexts.

	
	# First define a function that will return a matrix of pairs of p-values using Fisher's exact test, one pair per site - the first p-value is the FDR adjusted p-value that the site is significantly different from a Methylated site with errors, and the second that the site is significanty different from an Unmethylated site with errors
	# We will ultimately use this test to define a sample-specific region of #C,#T space for sites which cannot be reliably distinguished from either Methylated or Unmethylated status
	
	fisherTestMulti <- function(input_matrix, maxN = 100, Mprop_cutoff=0, Uprop_cutoff = 0, bh = TRUE) {
	
		# This function accepts a matrix of pairs of #C and #T values, and a cutoff for errors
		# For each pair it estimates the p-values/FDR of the pair being distinct from those expected from a methylated and from an unmethylated locus
		
		# Confidence level for Fisher test
		fisher_conf_level = 0.95
		
		# Set up two matrices, one for Fisher test results to rule out Unmethylated, and one to rule out Methylated 
        fmatU <- matrix(NA, maxN+1, maxN+1)
        fmatM <- matrix(NA, maxN+1, maxN+1)
		
		# Precompute the test results for C and T values between 0 and 50 and store in matrix
		# Matrix slots will have offset of 1 from relevant value: e.g. values 0-50, slots #1-#51
		
        for (i in 0:maxN) {
            for (j in 0:maxN) {
				Cprop=i/(i+j)
		
				# Minimum number of C calls expected if methylated with errors
				# Truncate to remain on the conservative side (control FP)
				expMC=trunc((i+j)*(1-Mprop_cutoff))
				# Maximum number of T calls expected if methylated with errors
				expMT=i+j-expMC
		
				# Minimum number of T calls expected if unmethylated with errors
				# Truncate to remain on the conservative side (control FP)
				expUT=trunc((i+j)*(1-Uprop_cutoff))
				# Maximum number of C calls expected if unmethylated with errors
				expUC=i+j-expUT

				fmatU[i+1,j+1] = fisher.test(matrix(c(i, j, expUC,expUT), byrow = TRUE, 2, 2), or = 1, alternative = "greater", conf.int = TRUE, conf.level = fisher_conf_level)$p.value
				fmatM[i+1,j+1] = fisher.test(matrix(c(i, j, expMC, expMT), byrow = TRUE, 2, 2), or = 1, alternative = "less", conf.int = TRUE, conf.level = fisher_conf_level)$p.value
			}
		}
        pVfun <- function(x) {
			#is.na test added to cope with NA values in input matrix
			if (is.na(x[2])) {
				# Can't differentiate from either M or U
				return(c(1,1))
			} else if ((x[1] < (maxN + 1)) & (x[2] < (maxN + 1))) { 
				# C and T counts are both small so we can use precomputed p-values
				p_valueU = fmatU[x[1]+1, x[2]+1]
				p_valueM = fmatM[x[1]+1, x[2]+1]
			}  else {
				# We need to calculate the p-values
				Cprop=x[1]/(x[1]+x[2])
		
				# Minimum number of C calls expected if methylated with errors
				# Truncate to remain on the conservative side (control FP)
				expMC=trunc((x[1]+x[2])*(1-Mprop_cutoff))
				# Maximum number of T calls expected if methylated with errors
				expMT=x[1]+x[2]-expMC
		
				# Minimum number of T calls expected if unmethylated with errors
				# Truncate to remain on the conservative side (control FP)
				expUT=trunc((x[1]+x[2])*(1-Uprop_cutoff))
				# Maximum number of C calls expected if unmethylated with errors
				expUC=x[1]+x[2]-expUT

				# p-value of pair being distinct from those expected from a U site
				p_valueU = fisher.test(matrix(c(x[1], x[2], expUC,expUT), byrow = TRUE, 2, 2), or = 1, alternative = "greater", conf.int = TRUE, conf.level = fisher_conf_level)$p.value
				# p-value of pair being distinct from those expected from a M site
				p_valueM = fisher.test(matrix(c(x[1], x[2], expMC, expMT), byrow = TRUE, 2, 2), or = 1, alternative = "less", conf.int = TRUE, conf.level = fisher_conf_level)$p.value
			}
			return(c(p_valueM, p_valueU))
		}
		
		pValues <- apply(input_matrix, 1, pVfun)
        if (bh) {
            pValues[1,] <- p.adjust(pValues[1,], method = "BH")
            pValues[2,] <- p.adjust(pValues[2,], method = "BH")
		}
		#pValues rescaling removed - I think they included this	step to turn p-value into a positive integer for convenience of storage
		#pValues <- -round(10 * log10(pValues))
        return(pValues)
    }
 	
	#sample_Cprop_cutoff is derived from the proportion of C calls in Chloroplast meth_context sites and represents an estimate of the bounds of error
	#Chloroplast should not be methylated, so any C calls (as opposed to T calls) result from one of the various sources of error
	#the sample_Cprop_cutoff is calculated by looking at the Cprop_cutoff quantile of the proportion of C calls among all Chloroplst meth_context sites in the given sample

	meth_data_fisher=matrix(, nrow=nrow(coverage_data), ncol=2)
		
	for(sample_name in rownames(sample_info)) {
		cat(paste0("Fisher's exact tests for sample ",sample_name,"\n"))
		meth_data_fisher[coverage_data$Sample==sample_name,] = t(fisherTestMulti(as.matrix(cbind(coverage_data[coverage_data$Sample==sample_name,]$Cov_C,coverage_data[coverage_data$Sample==sample_name,]$Cov_T)), Mprop_cutoff=sample_info[rownames(sample_info)==sample_name,]$Cprop_cutoff, Uprop_cutoff=sample_info[rownames(sample_info)==sample_name,]$Cprop_cutoff))
	}	

	# We now have a matrix of pairs of FDR-adjusted p-values, one pair per site, indicating whether the site is significantly distinct from the expected #C and #T coverage of a Methylated or from an Unmethylated site with errors
	
	# Save the Fisher's test results to cache file
	saveRDS(meth_data_fisher, file=paste0(project_id,"_",meth_context,"_meth_data_fisher_2.rds"))

	# We have used the Fisher's Exact test to identify sites which can't reliably be distinguished from U or M status
	# Now we use the binomial test to identify sites that can be reliably inferred to be Unmethylated
	
	
	# Function to calculate a matrix of Binomial p-values given an input matrix of C-coverage and C+T coverage values for loci
	# Original function taken from methylPipe$binomTestMulti.  Modifications are commented below.
    binomTestMulti <- function(mat, maxN = 100, p = bc, bh = TRUE) {
        bmat <- matrix(NA, maxN, maxN)
        for (i in 1:maxN) {
            for (j in i:maxN) bmat[i, j] <- binom.test(x = i, 
                n = j, p = p, alternative = "greater")$p.value
        }
        pVfun <- function(x) {
			#is.na test added to cope with NA values in input matrix
			if (is.na(x[2])) {
				return(NA)
			#x[1] == 0 test added after discussion with methylPipe author: https://support.bioconductor.org/p/102762/#102765
			} else if (x[1] == 0) {
				return(1)
			} else if (x[2] < (maxN + 1)) { 
                return(bmat[x[1], x[2]])
			}
            else {
				return(binom.test(x = x[1], n = x[2], p = p, 
                alternative = "greater")$p.value)
			}
		}
        pValues <- apply(mat, 1, pVfun)
        if (bh) {
            pValues <- p.adjust(pValues, method = "BH")
		}
		#pValues rescaling removed - they included this step to turn p-value into a positive integer for convenience of storage
		#pValues <- -round(10 * log10(pValues))
        return(pValues)
    }
 	
	meth_data_binom=matrix(, nrow=nrow(coverage_data), ncol=1)
	
	for(sample_name in rownames(sample_info)) {
		cat(paste0("Binomial tests for sample ",sample_name,"\n"))
		meth_data_binom[coverage_data$Sample==sample_name,] = binomTestMulti(as.matrix(cbind(coverage_data[coverage_data$Sample==sample_name,]$Cov_C,coverage_data[coverage_data$Sample==sample_name,]$Cov_C+coverage_data[coverage_data$Sample==sample_name,]$Cov_T)), p=sample_info[rownames(sample_info)==sample_name,]$Cprop_cutoff)
	}	

	# Save the Binomial test results to cache file
	saveRDS(meth_data_binom, file=paste0(project_id,"_",meth_context,"_meth_data_binom.rds"))

	#View(cbind(coverage_data,meth_data,meth_data_binom)[1:100,])

	
	# p-value/FDR for rejecting hypotheses of 'Methylated' or 'Unmethylated'
	# Making this limit more strict (smaller p-value) has the effect of excluding more sites due to insufficient coverage
	# 0.01 may be a better threshold to use when no replicates are available.  Where replicates are available, errors made in this step will be compensated for by comparison between replicates and discard of sites where replicates disagree
	#fisher_p_value_cutoff = 0.01
	# 0.01 may be too strict in removing low-coverage sites.  We use 0.05 and allow disagreement between replicates to control for errors made in this call.
	fisher_p_value_cutoff = 0.05

	# p-value/FDR for testing whether sites have too much methylation to be considered Unmethylated
	# Making this limit more strict (smaller p-value) has the effect of calling more sites (with a few methylated reads) as Unmethylated, rather than Methylated/Partial
	# Lister et al, 2009 used binomial test with B-H adjusted p-values for FDR<0.01
	# In their work, anything that passed this test was classified as having some methylation, but we will instead use failure to pass this test as a criterion to call sites Unmethylated
	# 0.01 still leaves a mode of mCG ~ 0.1 in the density distribution of 'partially' methylated sites.  0.005 removes this mode, so is likely to return more true positive 'U' sites.
	#binomial_p_value_threshold=0.01	
	binomial_p_value_threshold=0.005	

	# We also impose a strict cutoff on minimum coverage before a site will be considered for calling
	# 3 reads is a very generous minimum, allowing the most sites to be considered.  10 reads would be a conservative minimum, if the experiment contains enough data.  There is literature support for thresholds in this range:
	# Lister et al, 2009 >=3 reads coverage
	# Ziller et al, 2013 >=5 reads coverage
	# Kunde-Ramamoorthy et al, 2013 >=10 reads coverage
	# Cokus et al, 2008 >=5 reads coverage.
	binomial_coverage_cutoff=3

	
	# Using either the Fisher or binomial method alone, a large number of CG dinucleotides are identified as Methylated at one nucleotide and Unmethylated at the other.  This typically happens in cases where one of the sites in the pair has high coverage, and a clear 'U' or 'M' call, and the other has low coverage and a borderline call. For this reason, we use both methods to make a final call - we use Fisher's test to identify sites with ambiguous calls, which either pass or fail both tests (Unmethylated and Methylated).  We use the binomial test to identify sites clearly Unmethylated.  For the remainder (sites not Unmethylated, and with sufficient but not too high coverage), we consider them Methylated if more than a cutoff fraction of reads are methylated, otherwise we consider them partially methylated.  At the moment this fraction is chosen arbitrarily.  Ideally we could optimise this cutoff, however, while reducing it to maximise coverage and minimise sites lost as 'Partially methylated', we would likely also increase the number of adjacent CG dinucleotides showing inconsistent methylation calls.	
	
	# We impose a cutoff for the minimum proportion #C/(#C+#T) we accept to call the locus methylated. This differs from Lister et al, 2009 approach where they allowed the binomial test alone to call site methylation status.
	# According to Gardiner et al, 2015, Cokus et al, 2008 had identified that it is generally possible to impose a far higher proportion as the threshold in the CG context (they used 80%) than in other contexts (they used 25% and 10% for CHG and CHH respectively)
	### We should be able to work out this threshold empirically for CG context using ROC analysis on the basis that methylation status in adjacent CG pairs should be concordant
	### We might also consider whether different thresholds are appropriate dependant on whether each C is considered in isolation, or whether CG dinucleotides are considered as a single site
	
	# New approach to this: we try a range of values from 0.1 to 0.9 for partial_meth_ctuff, and instead of classifying sites below the line as 'P' we classify them as 'U'.  This increases the concordance among reps.  We maximise the mean proportion of sites with concordance among reps in the script 5-methylation_calls_optimise_cutoff_2018-02-23.R.  This comes out with a threshold of mCG=0.45 for the Schmitz data set.  We then use this cutoff, along with a more stringent binomial test threshold to distinguish unmethylated from methylated sites.
	
	partial_meth_cutoff=0
	if (meth_context=="CG") {
		partial_meth_cutoff = 0.45
	} else if (meth_context=="CHG") {
		partial_meth_cutoff = 0.25
	} else if (meth_context=="CHH") {
		partial_meth_cutoff = 0.1
	}
	
	# Create a matrix ready to hold results, for efficiency
	meth_data=matrix(, nrow=nrow(coverage_data), ncol=2)
	
	# Classify each site dependent on which of the two Fisher tests passes the FDR threshold
	fisher_class=ifelse((meth_data_fisher[,1]<fisher_p_value_cutoff),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),"P","U"),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),"M","I"))
	# Classify each site dependent on whether it passes the binomial test
	binom_class=ifelse((coverage_data$Cov_C+coverage_data$Cov_T)<binomial_coverage_cutoff,"I",ifelse(meth_data_binom<binomial_p_value_threshold,"M","U"))

	# Classify each site based on the relevant hybrid assessment
	# Added new test here for binom_class=="I" so that the minimum coverage cutoff is enforced properly
	meth_data[,1] = ifelse(binom_class=="I","I",ifelse(fisher_class=="I","I",ifelse(binom_class=="U","U",ifelse((coverage_data$Cov_C/(coverage_data$Cov_C+coverage_data$Cov_T))>=partial_meth_cutoff,"M","P"))))
	#Alternative approach: consider mC<partial_meth_cutoff calls to be U, rather than P
	#meth_data[,1] = ifelse(fisher_class=="I","I",ifelse(binom_class=="U","U",ifelse((coverage_data$Cov_C/(coverage_data$Cov_C+coverage_data$Cov_T))>=partial_meth_cutoff,"M","U")))

	# Choose the relevant p-value to retain, dependent on which status was chosen
	#### This needs more work
	#### Fisher's method (using a chisquared test) can be used to combine multiple p-values addressing the same general hypothesis.  Implemented in chiCombP function in the methylPipe package.
	meth_data[,2] = ifelse((meth_data_fisher[,1]<fisher_p_value_cutoff),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),pmax(meth_data_fisher[,1],meth_data_fisher[,2]),meth_data_fisher[,1]),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),meth_data_fisher[,2],pmin(meth_data_fisher[,1],meth_data_fisher[,2])))

	# Discard status and replace with "I" if coverage is too high to consider trusting the call (5 sigma)
	# This assumes that coverage is roughly Gaussian distributed (although it typically is long-tailed)
	for(sample_name in rownames(sample_info)) {
		coverage_upper_cutoff = sd(coverage_data[Sample==sample_name,Cov_C+Cov_T],na.rm=TRUE)*5
		meth_data[(coverage_data$Sample==sample_name) & (!is.na(coverage_data$Cov_C)) & ((coverage_data$Cov_C+coverage_data$Cov_T)>coverage_upper_cutoff),1] = "I"
	}
	
	# Organise meth_data from the large matrix to a data frame
	meth_data=data.frame(meth_data)
	colnames(meth_data) = c("meth_status","p_value")
	meth_data$p_value=as.numeric(levels(meth_data$p_value)[meth_data$p_value])
	
	# Plot the distribution of C proportions at sites with #C more or less than half coverage
	pdf(paste0(project_id,"_",meth_context,"_site_C_proportion_ge_0.5_by_sample.pdf"))
	print(ggplot(cbind(coverage_data,meth_data)[(meth_status!="I") & (Cov_C>=Cov_T),], aes(x=Cov_C/(Cov_C+Cov_T),colour=Chromosome)) + geom_freqpoly(binwidth = .01) + facet_wrap(~Sample))
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_site_C_proportion_lt_0.5_by_sample.pdf"))
	print(ggplot(cbind(coverage_data,meth_data)[(meth_status!="I") & (Cov_C<Cov_T),], aes(x=Cov_C/(Cov_C+Cov_T),colour=Chromosome)) + geom_freqpoly(binwidth = .01) + facet_wrap(~Sample))
	dev.off()
	
	# Plot C/T coverage distribution overall, and by methylation status
	plot_coverage_limit=100
	pdf(paste0(project_id,"_",meth_context,"_C_T_coverage_distribution_smooth.pdf"))
	print(smoothScatter(x=coverage_data[(coverage_data$Cov_T<=plot_coverage_limit) & (coverage_data$Cov_C<=plot_coverage_limit),]$Cov_T, y=coverage_data[(coverage_data$Cov_T<=plot_coverage_limit) & (coverage_data$Cov_C<=plot_coverage_limit),]$Cov_C))
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_C_T_coverage_distribution_point.pdf"))
	print(ggplot(cbind(coverage_data,meth_data)[sample(nrow(coverage_data),900000),][!is.na(Cov_C),], aes(x=Cov_T, y=Cov_C)) + geom_point(aes(colour = meth_status), alpha=0.05, size=2) + xlim(0,plot_coverage_limit) + ylim(0,plot_coverage_limit) + theme_minimal())
	dev.off()

	# Plot C/T coverage distribution as classified by Binomial test
	pdf(paste0(project_id,"_",meth_context,"_C_T_coverage_distribution_binomial_sig.pdf"))
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom,binom_class=ifelse(is.na(meth_data_binom),"I",ifelse((coverage_data$Cov_C+coverage_data$Cov_T)<binomial_coverage_cutoff,"I",ifelse(meth_data_binom<binomial_p_value_threshold,"M","U"))))[sample(nrow(coverage_data),300000),][(!is.na(Cov_C)),], aes(x=Cov_T, y=Cov_C)) + geom_point(aes(colour = binom_class.V1), alpha=0.05, size=2) + xlim(0,plot_coverage_limit) + ylim(0,plot_coverage_limit) + theme_minimal())
	dev.off()
	
	# Plot smoothed C/T coverage distribution for sites not classified "U", "M" or "I" by us
	pdf(paste0(project_id,"_",meth_context,"_C_T_coverage_distribution_ambig_sites_smooth.pdf"))
	#print(smoothScatter(cbind(coverage_data,meth_data)[(!is.na(Cov_C)) & (meth_status!="U") & (meth_status!="M") & (meth_status!="I") & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),]$Cov_T,cbind(coverage_data,meth_data)[(!is.na(Cov_C)) & (meth_status!="U") & (meth_status!="M")  & (meth_status!="I") & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),]$Cov_C))
	print(smoothScatter(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="U") & (V1<0.05) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),]$Cov_T,cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="U") & (V1<0.05) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),]$Cov_C))
	dev.off()

	# Plot distribution of C/C+T) for 'Unmethylated' sites (mCG<0.45) for a range of p-value cutoffs of Binomial test
	pdf(paste0(project_id,"_",meth_context,"_C_T_coverage_distribution_ambig_sites.pdf"))
	p_threshold = 1
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="U") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.5
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="U") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.1
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="U") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.05
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="U") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.01
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="U") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.005
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="U") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.001
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="U") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	dev.off()

	# Plot distribution of C/C+T) for 'Partial' sites (mCG<0.45) for a range of p-value cutoffs of Binomial test
	pdf(paste0(project_id,"_",meth_context,"_C_T_coverage_distribution_ambig_sites.pdf"))
	p_threshold = 1
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="P") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.5
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="P") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.1
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="P") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.05
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="P") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.01
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="P") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.005
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="P") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	p_threshold = 0.001
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[(!is.na(Cov_C)) & (meth_status=="P") & (V1<p_threshold) & (Cov_C<=plot_coverage_limit) & (Cov_T<=plot_coverage_limit),], aes(x=Cov_C/(Cov_C+Cov_T))) + geom_density() + annotate("text", x=0, y=-0.1, label= paste0("p<",p_threshold)) +coord_cartesian(xlim=c(0, 0.5)))
	dev.off()

	
	# Plot C/T coverage distribution for sites classified "U" by us, and "M" by Binomial test
	pdf(paste0(project_id,"_",meth_context,"_C_T_coverage_distribution_Usites_binomial_sig.pdf"))
	print(ggplot(cbind(coverage_data,meth_data,meth_data_binom)[sample(nrow(coverage_data),300000),][(!is.na(Cov_C)) & (meth_status=="U") & (V1<binomial_p_value_threshold),], aes(x=Cov_T, y=Cov_C)) + geom_point(aes(colour = meth_status), alpha=0.05, size=2) + xlim(0,plot_coverage_limit) + ylim(0,plot_coverage_limit) + theme_minimal())
	dev.off()

	# Count frequencies and overlaps between methylation classes by method
	count(cbind(coverage_data,meth_data,meth_data_binom,binom_class=ifelse(is.na(meth_data_binom),"I",ifelse(meth_data_binom<binomial_p_value_threshold,"M","U")))[,paste0(meth_status,binom_class.V1)])

	# Dump the methylation calls data frame to file for a convenience cache
	saveRDS(meth_data, file=paste0(project_id,"_",meth_context,"_meth_data.rds"))
	
} # end action=="call"


#### This part brings together methylation calls per site per sample into a per-site matrix for all samples

if ((action=="merge") || (action=="all")) {

	if (action=="merge") {
		#### Reload the data objects developed in earlier actions
		#all_samples_ct <- readRDS(paste0(project_id,"_",meth_context,"_all_samples_ct.rds")) # Commented out as not needed any more?
		sample_info <- readRDS(paste0(project_id,"_",meth_context,"_sample_info.rds"))
		coverage_data <- data.table(readRDS(paste0(project_id,"_",meth_context,"_coverage_data.rds")))
		meth_data <- data.table(readRDS(paste0(project_id,"_",meth_context,"_meth_data.rds"))) # Convert meth_data to a data.table for speed of query
	} else {
		meth_data=data.table(meth_data)
	}

	if (meth_context=="CG") {
		# No need to do anything special for CG context, as this was taken care of in the 'call' section
	}
	
	#### Build per-site tables of methylation calls across all samples
	
	# Create empty tables to accumulate each of the 'all samples' columns data
	all_samples_meth_status=NULL
	all_samples_p_value=NULL

	sample_no=0
	for(sample_name in rownames(sample_info)) {
		sample_no = sample_no+1
		if(opt$verbose) {cat(paste0("Adding methylation calls to all-samples per-site table: ",sample_name,"\n"))}

		sample_meth_status=NULL
		sample_p_value=NULL
		
		# Make separate tables for the sample's meth_call, p_value and meth_call_basis, keyed on locus
		sample_meth_status=cbind(coverage_data,meth_data)[coverage_data$Sample==sample_name,c(2,3,4,7)]
		sample_p_value=cbind(coverage_data,meth_data)[coverage_data$Sample==sample_name,c(2,3,4,8)]
		colnames(sample_meth_status)[4]=sample_name
		colnames(sample_p_value)[4]=sample_name

		# Merge each of the sample's 3 tables into the relevant accumulation table
		if (sample_no==1) {
			all_samples_meth_status=sample_meth_status
			all_samples_p_value=sample_p_value
		} else {
			all_samples_meth_status=merge(all_samples_meth_status, sample_meth_status, by=c("Chromosome","Locus","Strand"), all=TRUE)
			all_samples_p_value=merge(all_samples_p_value, sample_p_value, by=c("Chromosome","Locus","Strand"), all=TRUE)
		}

		# Tidy up the temporary objects
		rm(sample_meth_status, sample_p_value, sample_meth_call_basis)
	} # end for each sample

	#Add a column to identify the number of samples with missing data
	#######CG_pair_meth_status$m_count <- apply(CG_pair_meth_status, 1, function(x) sum(is.na(x))/2)
	
	if(opt$verbose) {
		### Numbers of consistent calls by replicate set:
	}
		
	# Sort the data frames sensibly (locus has received an alpha rather than numerical sort at this point)
	# Dump the site all-samples meth status data table to file for a convenience cache
	all_samples_meth_status=all_samples_meth_status[order(all_samples_meth_status$Chromosome,all_samples_meth_status$Locus),]
	saveRDS(all_samples_meth_status, file=paste0(project_id,"_",meth_context,"_all_samples_meth_status.rds"))
	all_samples_p_value=all_samples_p_value[order(all_samples_p_value$Chromosome,all_samples_p_value$Locus),]
	saveRDS(all_samples_p_value, file=paste0(project_id,"_",meth_context,"_all_samples_p_value.rds"))
} # end action=="merge"


#### Annotation and analysis

if ((action=="analyse") || (action=="all")) {

	if (action=="analyse") {
		#### Reload the data objects developed in earlier actions
		#all_samples_ct <- readRDS(paste0(project_id,"_",meth_context,"_all_samples_ct.rds")) # Commented out as not needed any more?
		sample_info <- readRDS(paste0(project_id,"_",meth_context,"_sample_info.rds"))
		coverage_data <- data.table(readRDS(paste0(project_id,"_",meth_context,"_coverage_data.rds")))
		meth_data <- data.table(readRDS(paste0(project_id,"_",meth_context,"_meth_data.rds")))
		meth_data_binom <- data.table(readRDS(paste0(project_id,"_",meth_context,"_meth_data_binom.rds")))
		meth_data_fisher <- data.table(readRDS(paste0(project_id,"_",meth_context,"_meth_data_fisher_2.rds")))
		# We reload the following as data tables, as the logic tests for all_m, all_u etc fail with data.tables
		all_samples_meth_status <- data.frame(readRDS(paste0(project_id,"_",meth_context,"_all_samples_meth_status.rds")))
		all_samples_p_value <- data.frame(readRDS(paste0(project_id,"_",meth_context,"_all_samples_p_value.rds")))
		if (meth_context=="CG") {
			# nothing special so far?
		}
	} else {
		all_samples_meth_status = data.frame(all_samples_meth_status)
		all_samples_p_value = data.frame(all_samples_p_value)
		if (meth_context=="CG") {
			# nothing special so far?
		}
	}

	
	# What is the highest p-value from the binomial test for a site we identified as "Partial"?
	cat(paste0("Max binomial p-value for any 'P' site:",max(cbind(meth_data,meth_data_binom)[meth_data$meth_status=="P",3]),"\n"))
	# What is the highest p-value from the binomial test for a site we identified as "Methylated"?
	cat(paste0("Max binomial p-value for any 'M' site:",max(cbind(meth_data,meth_data_binom)[meth_data$meth_status=="M",3]),"\n"))
	# What is the lowest p-value from the binomial test for a site we identified as "Indeterminate"?
	cat(paste0("Min binomial p-value for any 'I' site:",min(cbind(meth_data,meth_data_binom)[meth_data$meth_status=="I",3],na.rm=TRUE),"\n"))
	
	# What is the lowest p-value from the binomial test for a site we identified as "Indeterminate" where overall coverage >=10?
	cat(paste0("Min binomial p-value for any 'I' site with coverage>=10:",min(cbind(coverage_data,meth_data,meth_data_binom)[meth_data$meth_status=="I" & ((coverage_data$Cov_C+coverage_data$Cov_T)>=10),9],na.rm=TRUE),"\n"))
	
	# What is the lowest p-value from the binomial test for a site we identified as "Unmethylated"?
	cat(paste0("Min binomial p-value for any 'U' site:",min(cbind(meth_data,meth_data_binom)[meth_data$meth_status=="U",3],na.rm=TRUE),"\n"))
	
	# What is the lowest p-value from the binomial test for a site we identified as "Unmethylated" where overall coverage >=10?
	cat(paste0("Min binomial p-value for any 'U' site with coverage>=10:",min(cbind(coverage_data,meth_data,meth_data_binom)[meth_data$meth_status=="U" & ((coverage_data$Cov_C+coverage_data$Cov_T)>=10),9],na.rm=TRUE),"\n"))


	# Merge the sample metadata with the sample information
	sample_info=cbind(sample_info,sample_metadata)

	# Identify site with either an "M" or "U" call in all samples
	# Identify sites among these without identical calls in all samples

	# Some samples may need to be excluded from this analysis
	# For Becker and Schmitz data sets we want to exclude Line 69 as it is a methylation mutant
	excluded_samples = c()
	valid_samples=c()
	if (project_id=="SRA035939") {
		excluded_samples = c("SRR342380","SRR342390")
	} else if (project_id=="PRJEB2678") {
		excluded_samples = c("ERR046563","ERR046564")
	}
	valid_samples=rownames(sample_info)[!(rownames(sample_info) %in% excluded_samples)]
	valid_samples_meth_status=all_samples_meth_status[,!names(all_samples_meth_status) %in% excluded_samples]

	
	# Find out which sites are generally always methylated, and which unmethylated, across the genome

	# This part sums up methylation calls across lines to give an 'average' methylation level per site irrespective of line
	# Partial calls are assigned an arbitrary methylation level of 0.25
	average_methylation = matrix(0, nrow=nrow(all_samples_meth_status), ncol=2)
	for(sample_name in valid_samples) {
		if(opt$verbose) {cat(paste0("Summing methylation levels in sample ",sample_name,"\n"))}
		average_methylation[,1]=average_methylation[,1]+ifelse(all_samples_meth_status[,sample_name]=="U",0,ifelse(all_samples_meth_status[,sample_name]=="M",1,ifelse(all_samples_meth_status[,sample_name]=="P",0.25,0)))
		average_methylation[,2]=average_methylation[,2]+ifelse(all_samples_meth_status[,sample_name]=="U",1,ifelse(all_samples_meth_status[,sample_name]=="M",1,ifelse(all_samples_meth_status[,sample_name]=="P",1,0)))
	}
	
	
	# Accumulate methylation calls which are concordant among replicates per line_generation

	line_generation_no = 0
	for(this_line in unique(sample_info$Line)) {
		for(this_generation in unique(sample_info[sample_info$Line==this_line,]$Generation)) {
			sample_count=0
			reps_concordant_meth_status=NULL
			for(sample_name in rownames(sample_info[(!rownames(sample_info) %in% excluded_samples) & sample_info$Line==this_line & sample_info$Generation==this_generation,])) {
				cat(paste(this_line, this_generation, sample_name, "\n"))
				sample_count=sample_count+1
				if (sample_count == 1) {
					# This is the first replicate of replicate set
					line_generation_no = line_generation_no + 1
					# Cache the methylation status and move on
					reps_concordant_meth_status=all_samples_meth_status[,sample_name]
				} else {
					# We are doing a second or higher rep - check for concordance between this sample and the consensus of the previous reps_concordant_meth_status
					# If the current sample disagrees with the previous ones, record a status of "D" for disagreement/discordant for the locus, unless status is opposite (M/U) or "O" already found, in which case "O"
					reps_concordant_meth_status=as.factor(ifelse(levels(reps_concordant_meth_status)[reps_concordant_meth_status]==all_samples_meth_status[,sample_name],levels(reps_concordant_meth_status)[reps_concordant_meth_status],ifelse(((levels(reps_concordant_meth_status)[reps_concordant_meth_status]=="M") & (all_samples_meth_status[,sample_name]=="U")) | ((levels(reps_concordant_meth_status)[reps_concordant_meth_status]=="U") & (all_samples_meth_status[,sample_name]=="M")) | ((levels(reps_concordant_meth_status)[reps_concordant_meth_status]=="O") & ((all_samples_meth_status[,sample_name]=="M") | (all_samples_meth_status[,sample_name]=="U"))),"O","D")))
					### Could do something here with the Fisher method to combine p-values from the concordant reps
				}
			}
			# All reps done for this line and generation, report and cache the results
			if (sample_count>0) {
				if (line_generation_no==1) {
					all_reps_meth_status=cbind(all_samples_meth_status[,c("Chromosome","Locus","Strand")],reps_concordant_meth_status)
					#all_reps_p_value=sample_p_value
				} else {
					all_reps_meth_status=cbind(all_reps_meth_status, reps_concordant_meth_status)
					#all_reps_p_value=cbind(all_samples_p_value, sample_p_value, by=c("Chromosome","Locus","Strand"), all=TRUE)
				}
				colnames(all_reps_meth_status)[line_generation_no+3]=paste0("Line",this_line,"Gen",this_generation)
				cat(paste0(sum(all_reps_meth_status[,line_generation_no+3]=="O")," sites have opposite calls (M/U) among reps.\n"))
				# Tidy up the temporary objects
				rm(reps_concordant_meth_status, sample_count)
			}
		}
	}

	saveRDS(all_reps_meth_status, file=paste0(project_id,"_",meth_context,"_all_reps_meth_status.rds"))

	# Summarise the results of the rep concordance analysis
	line_gen_no = 0
	rep_concordance_summary=matrix(, nrow=ncol(all_reps_meth_status)-3, ncol=6)
	rownames(rep_concordance_summary)=colnames(all_reps_meth_status[,4:ncol(all_reps_meth_status)])
	for(line_gen in rownames(rep_concordance_summary)) {
		line_gen_no = line_gen_no + 1
		cat(paste0("Summarising methylation call basis for ",line_gen,"\n"))
		status_counts = count(all_reps_meth_status[,line_gen])
		# Needs a hack for when "P" is zero (nrow should be 6 otherwise):
		if(nrow(status_counts)==5) {
			status_counts$x = as.character(status_counts$x)
			status_counts[6,]=status_counts[5,]
			status_counts[5,]=c("P",0)
			status_counts$x = as.factor(status_counts$x)
			status_counts$freq=as.numeric(status_counts$freq)
		}
		rep_concordance_summary[line_gen_no,]=status_counts$freq
		colnames(rep_concordance_summary)=status_counts$x
	}
	# output the summary of numbers of concordant and discordant sites
	rep_concordance_summary
	# output the proportion of sites in each line with a 'clean' M/U call
	(rep_concordance_summary[,3]+rep_concordance_summary[,6])/sum(rep_concordance_summary[,1:6])*nrow(rep_concordance_summary)
	
	# This part looks for M/U calls across all samples, cognisant of rep structure
	line_gen_no=0
	first_line_gen = TRUE
	for(line_gen in rownames(rep_concordance_summary)) {
		line_gen_no = line_gen_no + 1
		if(opt$verbose) {cat(paste0("Checking variation in methylation calls concordant between reps: ",line_gen,"\n"))}
		if (first_line_gen) {
			all_m_or_u=ifelse(((all_reps_meth_status[,line_gen_no+3]=="U") | (all_reps_meth_status[,line_gen_no+3]=="M")),TRUE,FALSE)
			all_m=ifelse(all_reps_meth_status[,line_gen_no+3]=="M",TRUE,FALSE)
			all_u=ifelse(all_reps_meth_status[,line_gen_no+3]=="U",TRUE,FALSE)
			first_line_gen=FALSE
		} else {
			all_m_or_u=all_m_or_u & ifelse(((all_reps_meth_status[,line_gen_no+3]=="U") | (all_reps_meth_status[,line_gen_no+3]=="M")),TRUE,FALSE)
			all_m=all_m & ifelse(all_reps_meth_status[,line_gen_no+3]=="M",TRUE,FALSE)
			all_u=all_u & ifelse(all_reps_meth_status[,line_gen_no+3]=="U",TRUE,FALSE)
		}
	}
	
	# This version looks for M/U calls across all samples, independent of rep structure
	sample_no=0
	first_sample = TRUE
	for(sample_name in valid_samples) {
		sample_no = sample_no+1
		if(opt$verbose) {cat(paste0("Checking methylation calls in all-samples per-site table: ",sample_name,"\n"))}
		if (first_sample) {
			all_m_or_u_regardless=ifelse(((all_samples_meth_status[,sample_no+3]=="U") | (all_samples_meth_status[,sample_no+3]=="M")),TRUE,FALSE)
			#all_m=ifelse(all_samples_meth_status[,sample_no+3]=="M",TRUE,FALSE)
			#all_u=ifelse(all_samples_meth_status[,sample_no+3]=="U",TRUE,FALSE)
			first_sample=FALSE
		} else {
			all_m_or_u_regardless=all_m_or_u_regardless & ifelse(((all_samples_meth_status[,sample_no+3]=="U") | (all_samples_meth_status[,sample_no+3]=="M")),TRUE,FALSE)
			#all_m=all_m & ifelse(all_samples_meth_status[,sample_no+3]=="M",TRUE,FALSE)
			#all_u=all_u & ifelse(all_samples_meth_status[,sample_no+3]=="U",TRUE,FALSE)
		}
	}
	
	### What proportion of sites which vary between reps also have a clean M/U call in all samples?
	# This is our strict version of identifying variable sites.  It requires that all samples have sufficient coverage at the site, that there is agreement among reps in all lines, and that there are no differences between lines at generation 3.  One effect of this strictness is to filter out all highly labile sites.
		
	# This part looks for variable calls which are invariable at generation 3
	gen3=c()
	if (project_id=="SRA035939") {
		gen3=c("Line1Gen3","Line12Gen3","Line19Gen3")
		gen3_30=c("Line29Gen30","Line49Gen30","Line59Gen30","Line119Gen30")
	} else if(project_id=="PRJEB2678") {
		gen3=c("Line4Gen3","Line8Gen3")
	}
	
	line_gen_no=0
	first_line_gen = TRUE
	for(line_gen in gen3) {
		line_gen_no = line_gen_no + 1
		if(opt$verbose) {cat(paste0("Checking variation in methylation calls concordant within Gen 3: ",line_gen,"\n"))}
		if (first_line_gen) {
			#all_m_or_u_3=ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
			all_m_3=ifelse(all_reps_meth_status[,line_gen]=="M",TRUE,FALSE)
			all_u_3=ifelse(all_reps_meth_status[,line_gen]=="U",TRUE,FALSE)
			first_line_gen=FALSE
		} else {
			#all_m_or_u_3=all_m_or_u_3 & ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
			all_m_3=all_m_3 & ifelse(all_reps_meth_status[,line_gen]=="M",TRUE,FALSE)
			all_u_3=all_u_3 & ifelse(all_reps_meth_status[,line_gen]=="U",TRUE,FALSE)
		}
	}

	# This part looks for variable calls which are invariable at generation 3 and 30/1 
	### We should redo this analysis in a more sensitive way - check each generation 31 line for changes compared to the corresponding generation 30 line, in concert with reps, but ignoring coverage and calls in generation 3 and in all other generation 30/31 lines.  The current analysis does not really add anything useful.

	gen3_30=c()
	if (project_id=="SRA035939") {
		# analysis makes no sense for this project, as there is no generation 31
	} else if(project_id=="PRJEB2678") {
		gen3_30=c("Line4Gen3","Line8Gen3","Line29Gen31","Line39Gen31","Line49Gen31","Line59Gen31","Line79Gen31","Line89Gen31","Line99Gen31","Line109Gen31","Line119Gen31")
	}
	
	line_gen_no=0
	first_line_gen = TRUE
	for(line_gen in gen3_30) {
		line_gen_no = line_gen_no + 1
		if(opt$verbose) {cat(paste0("Checking variation in methylation calls concordant within Gen 3 and 30: ",line_gen,"\n"))}
		if (first_line_gen) {
			#all_m_or_u_3=ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
			all_m_3_30=ifelse(all_reps_meth_status[,line_gen]=="M",TRUE,FALSE)
			all_u_3_30=ifelse(all_reps_meth_status[,line_gen]=="U",TRUE,FALSE)
			first_line_gen=FALSE
		} else {
			#all_m_or_u_3=all_m_or_u_3 & ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
			all_m_3_30=all_m_3_30 & ifelse(all_reps_meth_status[,line_gen]=="M",TRUE,FALSE)
			all_u_3_30=all_u_3_30 & ifelse(all_reps_meth_status[,line_gen]=="U",TRUE,FALSE)
		}
	}

	# New more sensitive approach to calling variable sites:
	# Find a consensus call for parental lines: if at least two parental lines have calls concordant with reps, and they agree with each other, this is the consensus parental call
	
	gen3=c()
	gen3_30=c()
	if (project_id=="SRA035939") {
		gen3=c("Line1Gen3","Line12Gen3","Line19Gen3")
		gen3_30=c("Line29Gen30","Line49Gen30","Line59Gen30","Line119Gen30")
	} else if(project_id=="PRJEB2678") {
		gen3=c("Line4Gen3","Line8Gen3")
		gen3_30=c("Line4Gen3","Line8Gen3","Line29Gen31","Line39Gen31","Line49Gen31","Line59Gen31","Line79Gen31","Line89Gen31","Line99Gen31","Line109Gen31","Line119Gen31")
	}

	# Identify parental consensus calls
	line_gen_no=0
	first_line_gen = TRUE
	for(line_gen in gen3) {
		line_gen_no = line_gen_no + 1
		if(opt$verbose) {cat(paste0("Finding consensus call at Gen 3: ",line_gen,"\n"))}
		if (first_line_gen) {
			#all_m_or_u_3=ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
			no_m_3=ifelse(all_reps_meth_status[,line_gen]=="M",1,0)
			no_u_3=ifelse(all_reps_meth_status[,line_gen]=="U",1,0)
			first_line_gen=FALSE
		} else {
			#all_m_or_u_3=all_m_or_u_3 & ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
			no_m_3=no_m_3 + ifelse(all_reps_meth_status[,line_gen]=="M",1,0)
			no_u_3=no_u_3 + ifelse(all_reps_meth_status[,line_gen]=="U",1,0)
		}
	}
	parental_consensus = ifelse(no_m_3>=2,"M",ifelse(no_u_3>=2,"U","I"))

	# Identify sites where at least one offspring line differs from parental consensus, and also, count how many offspring differ at the locus
	line_gen_no=0
	first_line_gen = TRUE
	site_varies_from_parentals = NA
	num_lines_varying_from_parentals = NA
	num_lines_m_to_u_from_parentals = NA
	num_lines_u_to_m_from_parentals = NA
	num_lines_m_or_u_excl_parentals = NA
	for(line_gen in gen3_30) {
		line_gen_no = line_gen_no + 1
		if(opt$verbose) {cat(paste0("Checking variation in methylation calls from parental consensus at Gen 30: ",line_gen,"\n"))}
		if (first_line_gen) {
			#all_m_or_u_3=ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
			site_varies_from_parentals=ifelse((all_reps_meth_status[,line_gen] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_reps_meth_status[,line_gen]!=parental_consensus),TRUE,FALSE)
			num_lines_varying_from_parentals=ifelse((all_reps_meth_status[,line_gen] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_reps_meth_status[,line_gen]!=parental_consensus),1,0)
			num_lines_m_to_u_from_parentals=ifelse((all_reps_meth_status[,line_gen] == "U") & (parental_consensus == "M"),1,0)
			num_lines_u_to_m_from_parentals=ifelse((all_reps_meth_status[,line_gen] == "M") & (parental_consensus == "U"),1,0)
			num_lines_m_or_u_excl_parentals=ifelse((all_reps_meth_status[,line_gen] == "M") | (all_reps_meth_status[,line_gen] == "U"),1,0)
			first_line_gen=FALSE
		} else {
			#all_m_or_u_3=all_m_or_u_3 & ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
			site_varies_from_parentals=site_varies_from_parentals | ifelse((all_reps_meth_status[,line_gen] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_reps_meth_status[,line_gen]!=parental_consensus),TRUE,FALSE)
			num_lines_varying_from_parentals=num_lines_varying_from_parentals + ifelse((all_reps_meth_status[,line_gen] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_reps_meth_status[,line_gen]!=parental_consensus),1,0)
			num_lines_m_to_u_from_parentals=num_lines_m_to_u_from_parentals + ifelse((all_reps_meth_status[,line_gen] == "U") & (parental_consensus == "M"),1,0)
			num_lines_u_to_m_from_parentals=num_lines_u_to_m_from_parentals + ifelse((all_reps_meth_status[,line_gen] == "M") & (parental_consensus == "U"),1,0)
			num_lines_m_or_u_excl_parentals=num_lines_m_or_u_excl_parentals + ifelse((all_reps_meth_status[,line_gen] == "M") | (all_reps_meth_status[,line_gen] == "U"),1,0)
			#first_line_gen=FALSE
		}
	}

	
	
	if (opt$verbose) {
		cat(paste("No of",meth_context,"loci:",nrow(all_reps_meth_status),"\n"))
		cat(paste("No of",meth_context,"loci with all calls M or U:",nrow(all_reps_meth_status[all_m_or_u==TRUE,]),"\n"))
		cat(paste("No of",meth_context,"loci with all calls M:",nrow(all_reps_meth_status[all_m==TRUE,]),"\n"))
		cat(paste("No of",meth_context,"loci with all calls U:",nrow(all_reps_meth_status[all_u==TRUE,]),"\n"))
		cat(paste("No of",meth_context,"loci with strict variation in M/U calls (all parentals agree, all lines have agreement between reps):",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE)),]),"\n"))
		cat(paste("No of",meth_context,"loci with strict variation in M/U calls, M at generation 3:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3==TRUE)),]),"\n"))
		cat(paste("No of",meth_context,"loci with strict variation in M/U calls, U at generation 3:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3==TRUE)),]),"\n"))
		#cat(paste("No of",meth_context,"loci with strict variation in M/U calls, M at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3_30==TRUE)),]),"\n"))
		#cat(paste("No of",meth_context,"loci with strict variation in M/U calls, U at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3_30==TRUE)),]),"\n"))
		cat(paste("No of",meth_context,"loci with any variation in M/U calls (at least 2 parentals agree, at least one progeny line disagrees and has agreement between reps):",nrow(all_reps_meth_status[site_varies_from_parentals,]),"\n"))
		cat(paste("No of",meth_context,"loci with any variation in M/U calls, M at generation 3:",nrow(all_reps_meth_status[site_varies_from_parentals & (no_m_3>=2),]),"\n"))
		cat(paste("No of",meth_context,"loci with any variation in M/U calls, U at generation 3:",nrow(all_reps_meth_status[site_varies_from_parentals & (no_u_3>=2),]),"\n"))
		cat(paste("No of",meth_context,"loci with a clear call in parental consensus and all offspring:",nrow(cbind.data.frame(parental_consensus, site_varies_from_parentals, num_lines_m_to_u_from_parentals, num_lines_u_to_m_from_parentals, num_lines_m_or_u_excl_parentals)[(num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),]),"\n"))
		# Commented out reporting of generation 31 changes for now.  We should redo this analysis in a more sensitive way - check each generation 31 line for changes compared to the corresponding generation 30 line, in concert with reps, but ignoring coverage and calls in generation 3 and in all other generation 30/31 lines.
		#cat(paste("No of",meth_context,"loci with strict variation in M/U calls, M at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3_30==TRUE)),]),"\n"))
		#cat(paste("No of",meth_context,"loci with strict variation in M/U calls, U at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3_30==TRUE)),]),"\n"))
	}

	# Generate tables of sites with various patterns of calls
	across_samples_meth_status=cbind(all_samples_meth_status[,1:3],"average_methylation"=average_methylation[,1]/average_methylation[,2])[average_methylation[,2]>0,]
	#variant_calls=all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE)),]
	variant_calls=all_samples_meth_status[site_varies_from_parentals,]
	#M_to_U_calls=all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3==TRUE)),]
	#U_to_M_calls=all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3==TRUE)),]
	M_to_U_calls=all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]
	U_to_M_calls=all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]
	all_M_calls=all_samples_meth_status[all_m==TRUE,]
	all_U_calls=all_samples_meth_status[all_u==TRUE,]
	all_clean_calls=all_samples_meth_status[(num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),]
	
	M_parent_calls = all_samples_meth_status[parental_consensus=="M",]
	U_parent_calls = all_samples_meth_status[parental_consensus=="U",]

	# Make loci 2nt long for CG sites
	if (meth_context=="CG") {
		across_samples_meth_status=within(across_samples_meth_status, {Locus2=Locus+1})
		variant_calls=within(variant_calls, {Locus2=Locus+1})
		M_to_U_calls=within(M_to_U_calls, {Locus2=Locus+1})
		U_to_M_calls=within(U_to_M_calls, {Locus2=Locus+1})
		all_M_calls=within(all_M_calls, {Locus2=Locus+1})
		all_U_calls=within(all_U_calls, {Locus2=Locus+1})
		all_clean_calls=within(all_clean_calls, {Locus2=Locus+1})
		
		M_parent_calls=within(M_parent_calls, {Locus2=Locus+1})
		U_parent_calls=within(U_parent_calls, {Locus2=Locus+1})
	} else {
		across_samples_meth_status=within(across_samples_meth_status, {Locus2=Locus})
		variant_calls=within(variant_calls, {Locus2=Locus})
		M_to_U_calls=within(M_to_U_calls, {Locus2=Locus})
		U_to_M_calls=within(M_to_U_calls, {Locus2=Locus})
		all_M_calls=within(all_M_calls, {Locus2=Locus})
		all_U_calls=within(all_U_calls, {Locus2=Locus})
		all_clean_calls=within(all_clean_calls, {Locus2=Locus})

		M_parent_calls=within(M_parent_calls, {Locus2=Locus})
		U_parent_calls=within(U_parent_calls, {Locus2=Locus})
	}

	# Write out the 'average' status of sites across samples to a convenience table for external visualisations
	write.table(across_samples_meth_status[,c("Chromosome","Locus","Locus2","Strand","average_methylation")], file=paste0(project_id,"_",meth_context,"_average_methylation_across_samples.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Write out the variable sites to a convenience table for external visualisations
	write.table(variant_calls[,c("Chromosome","Locus","Locus2","Strand")], file=paste0(project_id,"_",meth_context,"_variable_sites.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Write out the M_to_U sites to a convenience table for external visualisations
	write.table(M_to_U_calls[,c("Chromosome","Locus","Locus2","Strand")], file=paste0(project_id,"_",meth_context,"_M_to_U_sites.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Write out the U_to_M sites to a convenience table for external visualisations
	write.table(U_to_M_calls[,c("Chromosome","Locus","Locus2","Strand")], file=paste0(project_id,"_",meth_context,"_U_to_M_sites.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Write out the all_M sites to a convenience table for external visualisations
	write.table(all_M_calls[,c("Chromosome","Locus","Locus2","Strand")], file=paste0(project_id,"_",meth_context,"_all_M_sites.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Write out the all_U sites to a convenience table for external visualisations
	write.table(all_U_calls[,c("Chromosome","Locus","Locus2","Strand")], file=paste0(project_id,"_",meth_context,"_all_U_sites.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Write out the all_clean sites to a convenience table for external visualisations
	write.table(all_clean_calls[,c("Chromosome","Locus","Locus2","Strand")], file=paste0(project_id,"_",meth_context,"_all_clean_sites.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Plot distribution of separation between variable loci
	# Generate a table of gaps between sites with variable calls
	gaps=matrix(, nrow=nrow(variant_calls), ncol=1)
	no_chroms=0
	prev_chrom="Z"
	for(locus_no in 1:nrow(variant_calls)) {
#	for(locus_no in 1:10) {
		if (prev_chrom==levels(variant_calls[locus_no,]$Chromosome)[variant_calls[locus_no,]$Chromosome]) {
			gaps[locus_no-no_chroms,1]=(variant_calls[locus_no,]$Locus)-(variant_calls[locus_no-1,]$Locus)
		} else {
			no_chroms = no_chroms + 1
		}
		prev_chrom=levels(variant_calls[locus_no,]$Chromosome)[variant_calls[locus_no,]$Chromosome]
	}

	gaps=cbind(data.frame(gaps),"Sites with variable calls among samples")
	colnames(gaps)=c("gap","class")	

	### The gaps analysis only really makes sense for CG sites - too few variable CHG sites to make any sense of gap plots currently, and CHH results unknown so far.
	# Generate a gap table the same size based on random loci which are methylated in all samples
	loci_sample <- sample(1:nrow(all_samples_meth_status[all_m==TRUE,]), nrow(variant_calls),
  	replace=FALSE)
	loci_sample=all_samples_meth_status[all_m==TRUE,][loci_sample,]
	loci_sample=loci_sample[order(loci_sample$Chromosome,loci_sample$Locus),]

	sample_gaps=matrix(, nrow=nrow(loci_sample), ncol=1)
	no_chroms=0
	prev_chrom="Z"
	for(locus_no in 1:nrow(loci_sample)) {
#	for(locus_no in 1:10) {
		if (prev_chrom==levels(loci_sample[locus_no,]$Chromosome)[loci_sample[locus_no,]$Chromosome]) {
			sample_gaps[locus_no-no_chroms,1]=(loci_sample[locus_no,]$Locus)-(loci_sample[locus_no-1,]$Locus)
		} else {
			no_chroms = no_chroms + 1
		}
		prev_chrom=levels(loci_sample[locus_no,]$Chromosome)[loci_sample[locus_no,]$Chromosome]
	}

	sample_gaps_m=cbind(data.frame(sample_gaps),"Sites with M calls in all samples")
	colnames(sample_gaps_m)=c("gap","class")

	# Generate a gap table the same size based on random loci which are unmethylated in all samples
	loci_sample <- sample(1:nrow(all_samples_meth_status[all_u==TRUE,]), nrow(variant_calls),
  	replace=FALSE)
	loci_sample=all_samples_meth_status[all_u==TRUE,][loci_sample,]
	loci_sample=loci_sample[order(loci_sample$Chromosome,loci_sample$Locus),]

	sample_gaps=matrix(, nrow=nrow(loci_sample), ncol=1)
	no_chroms=0
	prev_chrom="Z"
	for(locus_no in 1:nrow(loci_sample)) {
#	for(locus_no in 1:10) {
		if (prev_chrom==levels(loci_sample[locus_no,]$Chromosome)[loci_sample[locus_no,]$Chromosome]) {
			sample_gaps[locus_no-no_chroms,1]=(loci_sample[locus_no,]$Locus)-(loci_sample[locus_no-1,]$Locus)
		} else {
			no_chroms = no_chroms + 1
		}
		prev_chrom=levels(loci_sample[locus_no,]$Chromosome)[loci_sample[locus_no,]$Chromosome]
	}

	sample_gaps_u=cbind(data.frame(sample_gaps),"Sites with U calls in all samples")
	colnames(sample_gaps_u)=c("gap","class")
	
	# Plot distribution of gap sizes for variable sites vs. sites M in all samples or U in all samples
	pdf(paste0(project_id,"_",meth_context,"_variable_sites_gap_frequency.pdf"))
	print(ggplot(rbind(gaps, sample_gaps_m, sample_gaps_u), aes(x=gap, colour=class)) +geom_freqpoly(binwidth=5) + coord_cartesian(xlim = c(0, 5000)))
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_variable_sites_gap_frequency_detailed.pdf"))
	print(ggplot(rbind(gaps, sample_gaps_m, sample_gaps_u), aes(x=gap, colour=class)) +geom_freqpoly(binwidth=2) + coord_cartesian(xlim = c(0, 100)))
	dev.off()

		
	# Test whether the distributions of gap size vary significantly between groups
	ks.test(gaps$gap,sample_gaps_m$gap)
	ks.test(gaps$gap,sample_gaps_u$gap)
	ks.test(sample_gaps_m$gap,sample_gaps_u$gap)


	### This section annotates genomic features with characteristics of observed methylation patterns, and annotates CG sites with observed methylation patterns and details of genomic context

	library(GenomicRanges)
	
	# Load the annotation
	annot_gff = read.delim(reference_gff, header=F, comment.char="#")	
	gff.exons = read.delim(reference_exons, header=F, comment.char="#")
	gff.introns = read.delim(reference_introns, header=F, comment.char="#")
	gff.5UTR = read.delim(reference_5UTR, header=F, comment.char="#")
	gff.3UTR = read.delim(reference_3UTR, header=F, comment.char="#")

	# Grab the portion relating to genes
	gff.genes = annot_gff[annot_gff[,3]=="gene",]
	# Grab the portion relating to transposons
	gff.transposons = annot_gff[annot_gff[,3]=="transposable_element",]

	rm(annot_gff)
	
	# Convert chromosome names to uppercase to match previous objects
	gff.genes$V1=toupper(gff.genes$V1)
	gff.transposons$V1=toupper(gff.transposons$V1)
	gff.exons$V1=toupper(gff.exons$V1)
	gff.introns$V1=toupper(gff.introns$V1)
	gff.5UTR$V1=toupper(gff.5UTR$V1)
	gff.3UTR$V1=toupper(gff.3UTR$V1)
	
	# The following line sometimes fails if stringi has not been installed correctly.  Not sure why.  Fix is to install.packages("stringi") and to allow RStudio to restart R.
	# Grab gene IDs from descriptive field
	gff.genes$gene_ID=str_split_fixed(str_split_fixed(gff.genes$V9, ';',3)[,1],"=",2)[,2]
	gff.transposons$gene_ID=str_split_fixed(str_split_fixed(gff.transposons$V9, ';',3)[,1],"=",2)[,2]
	gff.exons$gene_ID=str_split_fixed(str_split_fixed(gff.exons$V9, ';',3)[,1],"=",2)[,2]
	gff.introns$gene_ID=str_split_fixed(str_split_fixed(gff.introns$V9, ';',3)[,1],"=",2)[,2]
	gff.5UTR$gene_ID=str_split_fixed(str_split_fixed(gff.5UTR$V9, ';',3)[,1],"=",2)[,2]
	gff.3UTR$gene_ID=str_split_fixed(str_split_fixed(gff.3UTR$V9, ';',3)[,1],"=",2)[,2]

	#Gene names are tricky - need to find the element containing "symbol="
	gff.genes=within(gff.genes, {gene_name=str_split_fixed(str_split_fixed(V9, ';',3)[,2],"=",2)[,2]})
	gff.transposons=within(gff.transposons, {gene_name=str_split_fixed(str_split_fixed(V9, ';',3)[,3],"=",2)[,2]})

	# Fix the gene names in exon, intron and UTRs
	gff.exons$exon_ID = gff.exons$gene_ID
	gff.exons$gene_ID = substr(gff.exons$gene_ID,1,9)
	gff.introns$intron_ID = gff.introns$gene_ID
	gff.introns$gene_ID = substr(gff.introns$gene_ID,1,9)
	gff.5UTR$UTR_ID = gff.5UTR$gene_ID
	gff.5UTR$gene_ID = substr(gff.5UTR$gene_ID,1,9)
	gff.3UTR$UTR_ID = gff.3UTR$gene_ID
	gff.3UTR$gene_ID = substr(gff.3UTR$gene_ID,1,9)
	
	# Add the exon and intron numbers
	gff.exons$exon_no = as.numeric(str_split_fixed(gff.exons$exon_ID, ':',3)[,3])
	gff.introns$exon_no = as.numeric(str_split_fixed(gff.introns$intron_ID, ':',3)[,3])
	
	# Find out how many exons each gene has
	install.packages("sqldf")
	library(sqldf)
	gff.genes=merge(gff.genes, sqldf('SELECT genes.gene_id, MAX(exon_no) AS no_exons FROM [gff.exons] exons, [gff.genes] genes WHERE genes.gene_ID=exons.gene_ID GROUP BY genes.gene_id'), by="gene_ID", all=TRUE)
	 
	# Find how many variable sites each gene has
	#gene_info = gff.genes
	#varloc_genes = NULL
	
	# Make a GRanges for gene space
	gene_ranges=makeGRangesFromDataFrame(df = gff.genes, start.field = "V4", end.field = "V5",seqnames.field = "V1")
	# Make a GRanges for transposon space
	transposon_ranges=makeGRangesFromDataFrame(df = gff.transposons, start.field = "V4", end.field = "V5",seqnames.field = "V1")
	# Make a GRanges for exon space
	exon_ranges=makeGRangesFromDataFrame(df = gff.exons, start.field = "V4", end.field = "V5",seqnames.field = "V1")
	# Make a GRanges for intron space
	intron_ranges=makeGRangesFromDataFrame(df = gff.introns, start.field = "V4", end.field = "V5",seqnames.field = "V1")
	# Make a GRanges for 5' UTR space
	UTR5_ranges=makeGRangesFromDataFrame(df = gff.5UTR, start.field = "V4", end.field = "V5",seqnames.field = "V1")
	# Make a GRanges for 3' UTR space
	UTR3_ranges=makeGRangesFromDataFrame(df = gff.3UTR, start.field = "V4", end.field = "V5",seqnames.field = "V1")
	
	#  Make a GRanges for all CG sites
	CG_site_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status, start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	# Make a GRanges for all_M sites
	all_M_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[all_m==TRUE,], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	# Make a GRanges for all_U sites
	all_U_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[all_u==TRUE,], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	# Make a GRanges for variable sites
	variant_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals,], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	# Make a GRanges for M->U sites
	m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	# Make a GRanges for U->M sites
	u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	# Make a GRanges for clean call sites
	clean_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	# Make 'clean' versions of the site ranges, where data is available for all lines
	# No need to make for all_M and all_U as they were already 'clean'
	#clean_all_M_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(all_m==TRUE) & (num_lines_m_or_u_excl_parentals==4),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	#clean_all_U_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(all_u==TRUE) & (num_lines_m_or_u_excl_parentals==4),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	clean_variant_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	clean_m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	clean_u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	# Make single change clean versions of the variable site ranges - sites where data is available for all lines, and where a change occurs but only in one offspring line
	clean_single_variant_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_varying_from_parentals == 1) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	clean_single_m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_varying_from_parentals == 1) & (no_m_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	clean_single_u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_varying_from_parentals == 1) & (no_u_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	
	M_parent_call_ranges = makeGRangesFromDataFrame(df = all_samples_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	U_parent_call_ranges = makeGRangesFromDataFrame(df = all_samples_meth_status[(parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")

	
	# Remove sites from the 'variable' sets, which were identified by Schmitz et al as being part of DMRs
	# We want to do this because we are focusing on spontaneous epimutations at individual sites, and want to remove the effect of coordinated changes at a set of nearby sites (DMRs)
	DMRs.gr = reduce(c(makeGRangesFromDataFrame(df=read.table(file=Schmitz_CG_DMRs, header=FALSE, sep="\t"), start.field = "V2", end.field = "V3",seqnames.field = "V1"),makeGRangesFromDataFrame(df=read.table(file=Schmitz_non_CG_DMRs, header=FALSE, sep="\t"), start.field = "V2", end.field = "V3",seqnames.field = "V1")))
	levels(DMRs.gr@seqnames@values) = toupper(as.character(levels(DMRs.gr@seqnames@values)))
	DMRs.gr@seqinfo@seqnames = levels(DMRs.gr@seqnames@values)

	variant_call_ranges = setdiff(variant_call_ranges, DMRs.gr, ignore.strand=TRUE)
	m_to_u_call_ranges = setdiff(m_to_u_call_ranges, DMRs.gr, ignore.strand=TRUE)
	u_to_m_call_ranges = setdiff(u_to_m_call_ranges, DMRs.gr, ignore.strand=TRUE)
	clean_call_ranges = setdiff(clean_call_ranges, DMRs.gr, ignore.strand = TRUE)
	clean_variant_call_ranges = setdiff(clean_variant_call_ranges, DMRs.gr, ignore.strand = TRUE)
	clean_m_to_u_call_ranges = setdiff(clean_m_to_u_call_ranges, DMRs.gr, ignore.strand = TRUE)
	clean_u_to_m_call_ranges = setdiff(clean_u_to_m_call_ranges, DMRs.gr, ignore.strand = TRUE)
	clean_single_variant_call_ranges = setdiff(clean_variant_call_ranges, DMRs.gr, ignore.strand = TRUE)
	clean_single_m_to_u_call_ranges = setdiff(clean_m_to_u_call_ranges, DMRs.gr, ignore.strand = TRUE)
	clean_single_u_to_m_call_ranges = setdiff(clean_u_to_m_call_ranges, DMRs.gr, ignore.strand = TRUE)
	M_parent_call_ranges = setdiff(M_parent_call_ranges, DMRs.gr, ignore.strand = TRUE)
	U_parent_call_ranges = setdiff(U_parent_call_ranges, DMRs.gr, ignore.strand = TRUE)
	
	# Find the overlaps
	olaps_CG_sites = findOverlaps(CG_site_ranges, gene_ranges)
	olaps_all_M = findOverlaps(all_M_ranges, gene_ranges)
	olaps_all_U = findOverlaps(all_U_ranges, gene_ranges)
	olaps_variable = findOverlaps(variant_call_ranges, gene_ranges)
	olaps_m_to_u = findOverlaps(m_to_u_call_ranges, gene_ranges)
	olaps_u_to_m = findOverlaps(u_to_m_call_ranges, gene_ranges)

	tolaps_CG_sites = findOverlaps(CG_site_ranges, transposon_ranges)
	tolaps_all_M = findOverlaps(all_M_ranges, transposon_ranges)
	tolaps_all_U = findOverlaps(all_U_ranges, transposon_ranges)
	tolaps_variable = findOverlaps(variant_call_ranges, transposon_ranges)
	tolaps_m_to_u = findOverlaps(m_to_u_call_ranges, transposon_ranges)
	tolaps_u_to_m = findOverlaps(u_to_m_call_ranges, transposon_ranges)

	eolaps_CG_sites = findOverlaps(CG_site_ranges, exon_ranges)
	eolaps_all_M = findOverlaps(all_M_ranges, exon_ranges)
	eolaps_all_U = findOverlaps(all_U_ranges, exon_ranges)
	eolaps_variable = findOverlaps(variant_call_ranges, exon_ranges)
	eolaps_m_to_u = findOverlaps(m_to_u_call_ranges, exon_ranges)
	eolaps_u_to_m = findOverlaps(u_to_m_call_ranges, exon_ranges)

	iolaps_CG_sites = findOverlaps(CG_site_ranges, intron_ranges)
	iolaps_all_M = findOverlaps(all_M_ranges, intron_ranges)
	iolaps_all_U = findOverlaps(all_U_ranges, intron_ranges)
	iolaps_variable = findOverlaps(variant_call_ranges, intron_ranges)
	iolaps_m_to_u = findOverlaps(m_to_u_call_ranges, intron_ranges)
	iolaps_u_to_m = findOverlaps(u_to_m_call_ranges, intron_ranges)

	olaps5_CG_sites = findOverlaps(CG_site_ranges, UTR5_ranges)
	olaps5_all_M = findOverlaps(all_M_ranges, UTR5_ranges)
	olaps5_all_U = findOverlaps(all_U_ranges, UTR5_ranges)
	olaps5_variable = findOverlaps(variant_call_ranges, UTR5_ranges)
	olaps5_m_to_u = findOverlaps(m_to_u_call_ranges, UTR5_ranges)
	olaps5_u_to_m = findOverlaps(u_to_m_call_ranges, UTR5_ranges)

	olaps3_CG_sites = findOverlaps(CG_site_ranges, UTR3_ranges)
	olaps3_all_M = findOverlaps(all_M_ranges, UTR3_ranges)
	olaps3_all_U = findOverlaps(all_U_ranges, UTR3_ranges)
	olaps3_variable = findOverlaps(variant_call_ranges, UTR3_ranges)
	olaps3_m_to_u = findOverlaps(m_to_u_call_ranges, UTR3_ranges)
	olaps3_u_to_m = findOverlaps(u_to_m_call_ranges, UTR3_ranges)


	# Generate a count of CG sites per gene
	gff.genes$CG_sites_count=table(gff.genes[subjectHits(olaps_CG_sites),]$gene_ID)[gff.genes$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.genes$CG_sites_count=ifelse(is.na(gff.genes$CG_sites_count),0,gff.genes$CG_sites_count)

	# Generate a count of all_M sites per gene
	gff.genes$all_M_count=table(gff.genes[subjectHits(olaps_all_M),]$gene_ID)[gff.genes$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.genes$all_M_count=ifelse(is.na(gff.genes$all_M_count),0,gff.genes$all_M_count)

	# Generate a count of all_U sites per gene
	gff.genes$all_U_count=table(gff.genes[subjectHits(olaps_all_U),]$gene_ID)[gff.genes$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.genes$all_U_count=ifelse(is.na(gff.genes$all_U_count),0,gff.genes$all_U_count)

	# Generate a count of variable sites per gene
	gff.genes$variable_count=table(gff.genes[subjectHits(olaps_variable),]$gene_ID)[gff.genes$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.genes$variable_count=ifelse(is.na(gff.genes$variable_count),0,gff.genes$variable_count)

	# Generate a count of m->u sites per gene
	gff.genes$m_to_u_count=table(gff.genes[subjectHits(olaps_m_to_u),]$gene_ID)[gff.genes$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.genes$m_to_u_count=ifelse(is.na(gff.genes$m_to_u_count),0,gff.genes$m_to_u_count)

	# Generate a count of u->m sites per gene
	gff.genes$u_to_m_count=table(gff.genes[subjectHits(olaps_u_to_m),]$gene_ID)[gff.genes$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.genes$u_to_m_count=ifelse(is.na(gff.genes$u_to_m_count),0,gff.genes$u_to_m_count)

	
	# Generate a count of CG sites per transposon
	gff.transposons$CG_sites_count=table(gff.transposons[subjectHits(tolaps_CG_sites),]$gene_ID)[gff.transposons$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.transposons$CG_sites_count=ifelse(is.na(gff.transposons$CG_sites_count),0,gff.transposons$CG_sites_count)

	# Generate a count of all_M sites per transposon
	gff.transposons$all_M_count=table(gff.transposons[subjectHits(tolaps_all_M),]$gene_ID)[gff.transposons$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.transposons$all_M_count=ifelse(is.na(gff.transposons$all_M_count),0,gff.transposons$all_M_count)

	# Generate a count of all_U sites per transposon
	gff.transposons$all_U_count=table(gff.transposons[subjectHits(tolaps_all_U),]$gene_ID)[gff.transposons$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.transposons$all_U_count=ifelse(is.na(gff.transposons$all_U_count),0,gff.transposons$all_U_count)

	# Generate a count of variable sites per transposon
	gff.transposons$variable_count=table(gff.transposons[subjectHits(tolaps_variable),]$gene_ID)[gff.transposons$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.transposons$variable_count=ifelse(is.na(gff.transposons$variable_count),0,gff.transposons$variable_count)

	# Generate a count of m->u sites per transposon
	gff.transposons$m_to_u_count=table(gff.transposons[subjectHits(tolaps_m_to_u),]$gene_ID)[gff.transposons$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.transposons$m_to_u_count=ifelse(is.na(gff.transposons$m_to_u_count),0,gff.transposons$m_to_u_count)

	# Generate a count of u->m sites per transposon
	gff.transposons$u_to_m_count=table(gff.transposons[subjectHits(tolaps_u_to_m),]$gene_ID)[gff.transposons$gene_ID]
	# Replace NA values in variant_count with 0s
	gff.transposons$u_to_m_count=ifelse(is.na(gff.transposons$u_to_m_count),0,gff.transposons$u_to_m_count)

	
	# Generate a count of CG sites per exon
	gff.exons$CG_sites_count=table(gff.exons[subjectHits(eolaps_CG_sites),]$exon_ID)[gff.exons$exon_ID]
	# Replace NA values in variant_count with 0s
	gff.exons$CG_sites_count=ifelse(is.na(gff.exons$CG_sites_count),0,gff.exons$CG_sites_count)

	# Generate a count of all_M sites per exon
	gff.exons$all_M_count=table(gff.exons[subjectHits(eolaps_all_M),]$exon_ID)[gff.exons$exon_ID]
	# Replace NA values in variant_count with 0s
	gff.exons$all_M_count=ifelse(is.na(gff.exons$all_M_count),0,gff.exons$all_M_count)

	# Generate a count of all_U sites per exon
	gff.exons$all_U_count=table(gff.exons[subjectHits(eolaps_all_U),]$exon_ID)[gff.exons$exon_ID]
	# Replace NA values in variant_count with 0s
	gff.exons$all_U_count=ifelse(is.na(gff.exons$all_U_count),0,gff.exons$all_U_count)

	# Generate a count of variable sites per exon
	gff.exons$variable_count=table(gff.exons[subjectHits(eolaps_variable),]$exon_ID)[gff.exons$exon_ID]
	# Replace NA values in variant_count with 0s
	gff.exons$variable_count=ifelse(is.na(gff.exons$variable_count),0,gff.exons$variable_count)

	# Generate a count of m->u sites per exon
	gff.exons$m_to_u_count=table(gff.exons[subjectHits(eolaps_m_to_u),]$exon_ID)[gff.exons$exon_ID]
	# Replace NA values in variant_count with 0s
	gff.exons$m_to_u_count=ifelse(is.na(gff.exons$m_to_u_count),0,gff.exons$m_to_u_count)

	# Generate a count of u->m sites per exon
	gff.exons$u_to_m_count=table(gff.exons[subjectHits(eolaps_u_to_m),]$exon_ID)[gff.exons$exon_ID]
	# Replace NA values in variant_count with 0s
	gff.exons$u_to_m_count=ifelse(is.na(gff.exons$u_to_m_count),0,gff.exons$u_to_m_count)

	
	# Generate a count of CG sites per intron
	gff.introns$CG_sites_count=table(gff.introns[subjectHits(iolaps_CG_sites),]$intron_ID)[gff.introns$intron_ID]
	# Replace NA values in variant_count with 0s
	gff.introns$CG_sites_count=ifelse(is.na(gff.introns$CG_sites_count),0,gff.introns$CG_sites_count)

	# Generate a count of all_M sites per intron
	gff.introns$all_M_count=table(gff.introns[subjectHits(iolaps_all_M),]$intron_ID)[gff.introns$intron_ID]
	# Replace NA values in variant_count with 0s
	gff.introns$all_M_count=ifelse(is.na(gff.introns$all_M_count),0,gff.introns$all_M_count)

	# Generate a count of all_U sites per intron
	gff.introns$all_U_count=table(gff.introns[subjectHits(iolaps_all_U),]$intron_ID)[gff.introns$intron_ID]
	# Replace NA values in variant_count with 0s
	gff.introns$all_U_count=ifelse(is.na(gff.introns$all_U_count),0,gff.introns$all_U_count)

	# Generate a count of variable sites per intron
	gff.introns$variable_count=table(gff.introns[subjectHits(iolaps_variable),]$intron_ID)[gff.introns$intron_ID]
	# Replace NA values in variant_count with 0s
	gff.introns$variable_count=ifelse(is.na(gff.introns$variable_count),0,gff.introns$variable_count)

	# Generate a count of m->u sites per intron
	gff.introns$m_to_u_count=table(gff.introns[subjectHits(iolaps_m_to_u),]$intron_ID)[gff.introns$intron_ID]
	# Replace NA values in variant_count with 0s
	gff.introns$m_to_u_count=ifelse(is.na(gff.introns$m_to_u_count),0,gff.introns$m_to_u_count)

	# Generate a count of u->m sites per intron
	gff.introns$u_to_m_count=table(gff.introns[subjectHits(iolaps_u_to_m),]$intron_ID)[gff.introns$intron_ID]
	# Replace NA values in variant_count with 0s
	gff.introns$u_to_m_count=ifelse(is.na(gff.introns$u_to_m_count),0,gff.introns$u_to_m_count)

	
	# Generate a count of CG sites per 5UTR
	gff.5UTR$CG_sites_count=table(gff.5UTR[subjectHits(olaps5_CG_sites),]$UTR_ID)[gff.5UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.5UTR$CG_sites_count=ifelse(is.na(gff.5UTR$CG_sites_count),0,gff.5UTR$CG_sites_count)

	# Generate a count of all_M sites per 5UTR
	gff.5UTR$all_M_count=table(gff.5UTR[subjectHits(olaps5_all_M),]$UTR_ID)[gff.5UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.5UTR$all_M_count=ifelse(is.na(gff.5UTR$all_M_count),0,gff.5UTR$all_M_count)

	# Generate a count of all_U sites per 5UTR
	gff.5UTR$all_U_count=table(gff.5UTR[subjectHits(olaps5_all_U),]$UTR_ID)[gff.5UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.5UTR$all_U_count=ifelse(is.na(gff.5UTR$all_U_count),0,gff.5UTR$all_U_count)

	# Generate a count of variable sites per 5UTR
	gff.5UTR$variable_count=table(gff.5UTR[subjectHits(olaps5_variable),]$UTR_ID)[gff.5UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.5UTR$variable_count=ifelse(is.na(gff.5UTR$variable_count),0,gff.5UTR$variable_count)

	# Generate a count of m->u sites per 5UTR
	gff.5UTR$m_to_u_count=table(gff.5UTR[subjectHits(olaps5_m_to_u),]$UTR_ID)[gff.5UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.5UTR$m_to_u_count=ifelse(is.na(gff.5UTR$m_to_u_count),0,gff.5UTR$m_to_u_count)

	# Generate a count of u->m sites per 5UTR
	gff.5UTR$u_to_m_count=table(gff.5UTR[subjectHits(olaps5_u_to_m),]$UTR_ID)[gff.5UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.5UTR$u_to_m_count=ifelse(is.na(gff.5UTR$u_to_m_count),0,gff.5UTR$u_to_m_count)

	
	# Generate a count of CG sites per 3UTR
	gff.3UTR$CG_sites_count=table(gff.3UTR[subjectHits(olaps_CG_sites),]$UTR_ID)[gff.3UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.3UTR$CG_sites_count=ifelse(is.na(gff.3UTR$CG_sites_count),0,gff.3UTR$CG_sites_count)

	# Generate a count of all_M sites per 3UTR
	gff.3UTR$all_M_count=table(gff.3UTR[subjectHits(olaps_all_M),]$UTR_ID)[gff.3UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.3UTR$all_M_count=ifelse(is.na(gff.3UTR$all_M_count),0,gff.3UTR$all_M_count)

	# Generate a count of all_U sites per 3UTR
	gff.3UTR$all_U_count=table(gff.3UTR[subjectHits(olaps_all_U),]$UTR_ID)[gff.3UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.3UTR$all_U_count=ifelse(is.na(gff.3UTR$all_U_count),0,gff.3UTR$all_U_count)

	# Generate a count of variable sites per 3UTR
	gff.3UTR$variable_count=table(gff.3UTR[subjectHits(olaps_variable),]$UTR_ID)[gff.3UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.3UTR$variable_count=ifelse(is.na(gff.3UTR$variable_count),0,gff.3UTR$variable_count)

	# Generate a count of m->u sites per 3UTR
	gff.3UTR$m_to_u_count=table(gff.3UTR[subjectHits(olaps_m_to_u),]$UTR_ID)[gff.3UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.3UTR$m_to_u_count=ifelse(is.na(gff.3UTR$m_to_u_count),0,gff.3UTR$m_to_u_count)

	# Generate a count of u->m sites per 3UTR
	gff.3UTR$u_to_m_count=table(gff.3UTR[subjectHits(olaps_u_to_m),]$UTR_ID)[gff.3UTR$UTR_ID]
	# Replace NA values in variant_count with 0s
	gff.3UTR$u_to_m_count=ifelse(is.na(gff.3UTR$u_to_m_count),0,gff.3UTR$u_to_m_count)

	
	# Find out which genes are really genes, and which are 'heterochromatic genes' as per Zemach et al, 2013 Figure 6A.
	gff.genes$average_methylation=table(gff.genes[subjectHits(olaps_CG_sites),]$gene_ID)[gff.genes$gene_ID]

	# Find out which transposons are really transposons, and which are 'euchromatic transposons' as per my own method.
	gff.transposons$average_methylation=table(gff.transposons[subjectHits(tolaps_CG_sites),]$gene_ID)[gff.transposons$gene_ID]

	
	# Create a data frame with gene ID and average methylation for each overlap between a site and a gene model
	olaps_average_methylation=data.frame(cbind("gene_ID"=gff.genes[subjectHits(olaps_CG_sites),"gene_ID"],"site_average_methylation"=(average_methylation[,1]/average_methylation[,2])[queryHits(olaps_CG_sites)]))
	tolaps_average_methylation=data.frame(cbind("gene_ID"=gff.transposons[subjectHits(tolaps_CG_sites),"gene_ID"],"site_average_methylation"=(average_methylation[,1]/average_methylation[,2])[queryHits(tolaps_CG_sites)]))
	# Convert average values to numerics
	olaps_average_methylation$site_average_methylation=as.numeric(levels(olaps_average_methylation$site_average_methylation)[olaps_average_methylation$site_average_methylation])
	tolaps_average_methylation$site_average_methylation=as.numeric(levels(tolaps_average_methylation$site_average_methylation)[tolaps_average_methylation$site_average_methylation])

	# Find the mean of sites average methylation levels across each gene, and merge this with the genes table
	#x = merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], mean), by="gene_ID", all=TRUE)
	#gff.genes$average_methylation = x$site_average_methylation[match(gff.genes$gene_ID, x$gene_ID)]
	gff.genes$average_methylation = merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], mean), by="gene_ID", all=TRUE)$site_average_methylation[match(gff.genes$gene_ID, merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], mean), by="gene_ID", all=TRUE)$gene_ID)]
	gff.transposons$average_methylation = merge(gff.transposons, aggregate(. ~ gene_ID, tolaps_average_methylation[], mean), by="gene_ID", all=TRUE)$site_average_methylation[match(gff.transposons$gene_ID, merge(gff.transposons, aggregate(. ~ gene_ID, tolaps_average_methylation[], mean), by="gene_ID", all=TRUE)$gene_ID)]

	# Find the variance of sites average methylation levels across each gene, and merge this with the genes table
	gff.genes$variance_methylation = merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], var), by="gene_ID", all=TRUE)$site_average_methylation[match(gff.genes$gene_ID, merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], var), by="gene_ID", all=TRUE)$gene_ID)]
	gff.transposons$variance_methylation = merge(gff.transposons, aggregate(. ~ gene_ID, tolaps_average_methylation[], var), by="gene_ID", all=TRUE)$site_average_methylation[match(gff.transposons$gene_ID, merge(gff.transposons, aggregate(. ~ gene_ID, tolaps_average_methylation[], var), by="gene_ID", all=TRUE)$gene_ID)]
 

	# Check correlation between our average methylation metric, and the proportion of sites within the gene represented by 'all M' as a proportion of those with a good call in all samples
	cat(paste0("Gene-wise correlation between average methylation across all samples and count of 'all M' sites as proportion of all sites which are M/U in all samples: ",cor(gff.genes$average_methylation,(gff.genes$all_M_count/(gff.genes$all_M_count+gff.genes$all_U_count+gff.genes$variable_count)), use = "pairwise.complete.obs"),"\n"))
	cat(paste0("Transposon-wise correlation between average methylation across all samples and count of 'all M' sites as proportion of all sites which are M/U in all samples: ",cor(gff.transposons$average_methylation,(gff.transposons$all_M_count/(gff.transposons$all_M_count+gff.transposons$all_U_count+gff.transposons$variable_count)), use = "pairwise.complete.obs"),"\n"))
	
	
	cat(paste0("Maximum level of methylation of a ChrM gene: ",max(gff.genes[gff.genes$V1=="CHRM",]$variance_methylation/gff.genes[gff.genes$V1=="CHRM",]$average_methylation,na.rm=TRUE),"\n"))
	cat(paste0("Maximum level of methylation of a ChrM transposon: ",max(gff.transposons[gff.transposons$V1=="CHRM",]$variance_methylation/gff.transposons[gff.transposons$V1=="CHRM",]$average_methylation,na.rm=TRUE),"\n"))

	# Schmitz data:
	#Maximum level of CG methylation of a ChrM gene: 0.172713661981165.  After converting P to U/M with 0.45 mCG cutoff: 0.327320827320827. After setting Fishers p<0.05, Binomial p<0.005: 0.178095863239987
	#Maximum level of CG methylation of a ChrM transposon: -Inf
	# This provides a rationale for heterochromatic_gene_coverage_minimum = 0.2
	
	# Becker data:
	#Maximum level of CG methylation of a ChrM gene: 0.654437081510604
	#Maximum level of CG methylation of a ChrM transposon: -Inf
	
	cat(paste0("Minimum dispersion index of a mitochondrial gene: ",min(gff.genes[gff.genes$V1=="CHRM",]$average_methylation),"\n"))
	cat(paste0("Minimum dispersion index of a mitochondrial transposon: ",min(gff.transposons[gff.transposons$V1=="CHRM",]$average_methylation),"\n"))
	# Schmitz data:
	#Minimum dispersion index of a mitochondrial gene: 0.106599283069871.  After converting P to U/M with 0.45 mCG cutoff:  0.  After setting Fishers p<0.05, Binomial p<0.005: 0.109324499030381
	# This provides supporting evidence for heterochromatic_gene_dispersion_cutoff=0.25
	#Minimum dispersion index of a mitochondrial transposon: Inf
	
	# Becker data:
	#Minimum dispersion index of a mitochondrial gene: 0
	#Minimum dispersion index of a mitochondrial transposon: Inf
	
	
	# Plot distribution of density of dispersion index
	pdf(paste0(project_id,"_",meth_context,"_gene_methylation_dispersion_index_density.pdf"))
	#print(ggplot(gff.genes, aes(x=variance_methylation/average_methylation)) + geom_histogram(aes(y = ..density..), binwidth=0.01) + stat_function(fun = dnorm, colour = "red", args = list(mean = mean(gff.genes[gff.genes$variance_methylation/gff.genes$average_methylation>0.25,]$variance_methylation/gff.genes[gff.genes$variance_methylation/gff.genes$average_methylation>0.25,]$average_methylation, na.rm = TRUE), sd = sd(gff.genes[gff.genes$variance_methylation/gff.genes$average_methylation>0.25,]$variance_methylation/gff.genes[gff.genes$variance_methylation/gff.genes$average_methylation>0.25,]$average_methylation, na.rm = TRUE))) + theme_minimal())
	print(ggplot(gff.genes, aes(x=variance_methylation/average_methylation)) + geom_histogram(aes(y = ..density..), binwidth=0.01) + theme_minimal())
	dev.off()
	# This plot provides a rationale for a heterochromatic_gene_dispersion_cutoff around 0.25 in both Schmitz and Becker data

	pdf(paste0(project_id,"_",meth_context,"_transposon_methylation_dispersion_index_density.pdf"))
	#print(ggplot(gff.genes, aes(x=variance_methylation/average_methylation)) + geom_histogram(aes(y = ..density..), binwidth=0.01) + stat_function(fun = dnorm, colour = "red", args = list(mean = mean(gff.genes[gff.genes$variance_methylation/gff.genes$average_methylation>0.25,]$variance_methylation/gff.genes[gff.genes$variance_methylation/gff.genes$average_methylation>0.25,]$average_methylation, na.rm = TRUE), sd = sd(gff.genes[gff.genes$variance_methylation/gff.genes$average_methylation>0.25,]$variance_methylation/gff.genes[gff.genes$variance_methylation/gff.genes$average_methylation>0.25,]$average_methylation, na.rm = TRUE))) + theme_minimal())
	print(ggplot(gff.transposons, aes(x=variance_methylation/average_methylation)) + geom_histogram(aes(y = ..density..), binwidth=0.01) + theme_minimal())
	dev.off()
	# This plot also provides support for a heterochromatic_gene_dispersion_cutoff around 0.25 in both Schmitz and Becker data

	
	### For CG methylation a mixture of 3 Gaussians looks like a good fit for the gene methylation dispersion index density.
	### For CHG methylation, one of the variances in the models would go to zero, so two Gaussians may be a better fit.  One of the fit curves is nearly flat while the other is a big peak around 1.0.
	
	# For a more hands-off way to set the threshold, fit a mixture of Gaussians to the dispersion index density, and identify the intersection of the fit models
	library(mixtools)
	library(rootSolve)
	# Fit a mixture of 2 Gaussians to the dispersion index density
	#mixmdl = normalmixEM(gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$variance_methylation/gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$average_methylation)

	# Fit a mixture of 3 Gaussians to the dispersion index density
	mixmdl = normalmixEM(gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$variance_methylation/gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$average_methylation, k=ifelse(meth_context=="CG",3,ifelse(meth_context=="CHG",2,2)))
	
	# Mixture of 2 Gaussians fits better once 'P' CG sites are removed from the picture
	#mixmdl = normalmixEM(gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$variance_methylation/gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$average_methylation, k=ifelse(meth_context=="CG",2,ifelse(meth_context=="CHG",2,2)))
	
	# Print the characteristics of the fit curves
	cat(paste0("Fit lambda values (mixture proportions):\n"))
	mixmdl$lambda
	# Schmitz data: 0.2853651 0.1365236 0.5781114.  After converting P to U/M with mCG=0.45 cutoff: 0.2381327 0.7618673. Using Fishers p<0.05, Binomial p<0.005: 0.2413691 0.1474132 0.6112177
	# Becker data:  0.08638043 0.17313572 0.74048385
	
	cat(paste0("Fit mu values (means):\n"))
	mixmdl$mu
	# Schmitz data: 0.01777965 0.10733632 0.64071733.  After converting P to U/M with mCG=0.45 cutoff: 0.08997885 0.65362097. Using Fishers p<0.05, Binomial p<0.005: 0.01632469 0.10906799 0.64090168
	# Becker data:  0.003837026 0.070503558 0.635220579
	
	cat(paste0("Fit sigma values (st.devs):\n"))
	mixmdl$sigma
	# Schmitz data: 0.01006439 0.05852205 0.14541845.  After converting P to U/M with mCG=0.45 cutoff: 0.07725693 0.14292752. Using Fishers p<0.05, Binomial p<0.005: 0.01027012 0.05664957 0.14650136
	# Becker data:  0.004954326 0.045214098 0.165311575
	
	# Find the intersection of the two larger components in the mixture - here we cut off
	mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
	dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }
	low_model=ifelse(meth_context=="CG",2,ifelse(meth_context=="CHG",1,1))
	#low_model=1
	high_model=low_model+1
	heterochromatic_gene_dispersion_cutoff=uniroot.all(mixfn, lower=0, upper=1, m1=mixmdl$mu[low_model], sd1=mixmdl$sigma[low_model], m2=mixmdl$mu[high_model], sd2=mixmdl$sigma[high_model], p1=mixmdl$lambda[low_model], p2=mixmdl$lambda[high_model])
	heterochromatic_gene_dispersion_cutoff
	# Schmitz CG data: 0.2517447, Becker CG data: 0.1896873
	# Schmitz CG data after 'P' sites removed: 0.2768782, Becker:
	# Schmitz data with Fishers p<0.05, Binomial p<0.005: 0.2499017
	# Becker CHG data: 0.9766602
	
	pdf(paste0(project_id,"_",meth_context,"_gene_methylation_dispersion_index_density_fitted.pdf"))
	print(plot(mixmdl,which=2))
	print(lines(density(gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$variance_methylation/gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$average_methylation), lty=2, lwd=2))
	dev.off()
	
	
	
	# Plot gene methylation level vs. gene length, colour by chromosome
	heterochromatic_gene_coverage_cutoff=0.6
	heterochromatic_gene_coverage_minimum=0.2
	#heterochromatic_gene_dispersion_cutoff=0.25  # this is now set above by finding the intersection of the two populations of dispersion index values from the mixture 
	
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_and variance_by_gene.pdf"))
	print(ggplot(gff.genes[,], aes(x=m_class, y=average_methylation, colour=m_class)) + geom_violin(scale="width") + theme_minimal())
	print(ggplot(gff.genes[,], aes(x=m_class, y=variance_methylation, colour=m_class)) + geom_violin(scale="width") + theme_minimal())
	dev.off()
	
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_by_gene_length.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=V1) ,alpha=0.4, size=1) + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Chromosome`")))
	dev.off()

	pdf(paste0(project_id,"_",meth_context,"_methylation_level_by_transposon_length.pdf"))
	print(ggplot(gff.transposons, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=V1) ,alpha=0.4, size=1) + theme_minimal() + xlab("Transposon length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within transposon")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Chromosome`")))
	dev.off()

	# Plot gene methylation level vs. gene length, colour by variance in methylation along length of gene
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_variance_by_gene_length.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variance_methylation) ,alpha=0.4, size=1) + scale_colour_gradient(low="red", high="green") + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff)) + guides(color=guide_legend(title="Variance along gene")))
	dev.off()

	pdf(paste0(project_id,"_",meth_context,"_methylation_level_variance_by_transposon_length.pdf"))
	print(ggplot(gff.transposons, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variance_methylation) ,alpha=0.4, size=1) + scale_colour_gradient(low="red", high="green") + theme_minimal() + xlab("Transposon length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within transposon")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff)) + guides(color=guide_legend(title="Variance along transposon")))
	dev.off()
	
	# Plot gene methylation level vs. gene length, colour by dispersion of methylation along length of gene
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_dispersion_by_gene_length.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variance_methylation/average_methylation) ,alpha=0.4, size=1) + scale_colour_gradient(low="red", high="green") + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Dispersion along gene")))
	dev.off()
	
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_dispersion_by_transposon_length.pdf"))
	print(ggplot(gff.transposons, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variance_methylation/average_methylation) ,alpha=0.4, size=1) + scale_colour_gradient(low="red", high="green") + theme_minimal() + xlab("Transposon length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within transposon")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Dispersion along transposon")))
	dev.off()

	# Plot gene methylation level vs. gene length, highlight mitochondrial genes
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_variance_by_gene_length_CHRM.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=gff.genes$V1=="CHRM") ,alpha=0.4, size=1) + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Mitochondrial genes")))
	dev.off()
	
	# Plot gene methylation level vs. gene length, highlight Chr2 pericentromeric region
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_variance_by_gene_length_CHR2_peri_C.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=(gff.genes$V1=="CHR2" & gff.genes$V4>3240000 & gff.genes$V5<3500000)) ,alpha=0.4, size=1)  + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Chr2 3.24-3.50 MB")))
	dev.off()
	
	# Plot gene methylation level vs. gene length, highlight genes with methylation dispersion index>cutoff (inverse for transposons)
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_variance_by_gene_length_dispersion_cutoff.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variance_methylation/average_methylation>heterochromatic_gene_dispersion_cutoff) ,alpha=0.4, size=1)  + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff)) + guides(color=guide_legend(title="Dispersion index>cutoff")))
	dev.off()
	
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_variance_by_transposon_length_dispersion_cutoff.pdf"))
	print(ggplot(gff.transposons, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variance_methylation/average_methylation<=heterochromatic_gene_dispersion_cutoff) ,alpha=0.4, size=1)  + theme_minimal() + xlab("Transposon length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within transposon")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff)) + guides(color=guide_legend(title="Dispersion index>cutoff")))
	dev.off()
	
	# Plot gene methylation level vs. gene length
	pdf(paste0(project_id,"_",meth_context,"_methylation_proportion_by_gene_length.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=(all_M_count)/(all_M_count+all_U_count+variable_count))) + geom_point(aes(colour=V1) ,alpha=0.4, size=1) + theme_minimal() + xlim(0,10000))
	dev.off()

	pdf(paste0(project_id,"_",meth_context,"_methylation_proportion_by_transposon_length.pdf"))
	print(ggplot(gff.transposons, aes(x=V5-V4, y=(all_M_count)/(all_M_count+all_U_count+variable_count))) + geom_point(aes(colour=V1) ,alpha=0.4, size=1) + theme_minimal() + xlim(0,10000))
	dev.off()

	cat(paste0(nrow(gff.genes), " annotated genes\n"))
	cat(paste0(nrow(gff.genes[gff.genes$average_methylation<heterochromatic_gene_coverage_cutoff,]), " annotated genes with average_methylation<",heterochromatic_gene_coverage_cutoff,"\n"))
	cat(paste0(nrow(gff.genes[gff.genes$average_methylation<heterochromatic_gene_coverage_cutoff & (gff.genes$average_methylation>heterochromatic_gene_coverage_minimum & gff.genes$variance_methylation/gff.genes$average_methylation>heterochromatic_gene_dispersion_cutoff | gff.genes$average_methylation<=heterochromatic_gene_coverage_minimum),]), " annotated genes with average_methylation<",heterochromatic_gene_coverage_cutoff," and dispersion index >",heterochromatic_gene_dispersion_cutoff," where average_methylation>",heterochromatic_gene_coverage_minimum," or average_methylation<",heterochromatic_gene_coverage_minimum,"\n"))
	cat(paste0(	nrow(gff.genes[gff.genes$variance_methylation/gff.genes$average_methylation>heterochromatic_gene_dispersion_cutoff,]), " annotated genes with methylation dispersion index>",heterochromatic_gene_dispersion_cutoff,"\n"))
	cat(paste0(	nrow(gff.genes[gff.genes$average_methylation<heterochromatic_gene_coverage_cutoff & gff.genes$variance_methylation/gff.genes$average_methylation>heterochromatic_gene_dispersion_cutoff,]), " annotated genes with average methylation<",heterochromatic_gene_coverage_cutoff," and dispersion index>",heterochromatic_gene_dispersion_cutoff,"\n"))

	cat(paste0(nrow(gff.transposons), " annotated transposons\n"))
	cat(paste0(nrow(gff.transposons[gff.transposons$average_methylation<heterochromatic_gene_coverage_cutoff,]), " annotated transposons with average_methylation<",heterochromatic_gene_coverage_cutoff,"\n"))
	cat(paste0(nrow(gff.transposons[gff.transposons$average_methylation<heterochromatic_gene_coverage_cutoff & (gff.transposons$average_methylation>heterochromatic_gene_coverage_minimum & gff.transposons$variance_methylation/gff.transposons$average_methylation>heterochromatic_gene_dispersion_cutoff | gff.transposons$average_methylation<=heterochromatic_gene_coverage_minimum),]), " annotated transposons with average_methylation<",heterochromatic_gene_coverage_cutoff," and dispersion index >",heterochromatic_gene_dispersion_cutoff," where average_methylation>",heterochromatic_gene_coverage_minimum," or average_methylation<",heterochromatic_gene_coverage_minimum,"\n"))
	cat(paste0(	nrow(gff.transposons[gff.transposons$variance_methylation/gff.transposons$average_methylation>heterochromatic_gene_dispersion_cutoff,]), " annotated transposons with methylation dispersion index>",heterochromatic_gene_dispersion_cutoff,"\n"))
	cat(paste0(	nrow(gff.transposons[gff.transposons$average_methylation<heterochromatic_gene_coverage_cutoff & gff.transposons$variance_methylation/gff.transposons$average_methylation>heterochromatic_gene_dispersion_cutoff,]), " annotated transposons with average methylation<",heterochromatic_gene_coverage_cutoff," and dispersion index>",heterochromatic_gene_dispersion_cutoff,"\n"))

	# Our new definition of heterochromatic genes is those with average methylation>0.6, or with average_methylation>0.2 and dispersion_index<0.25 ish
	heterochromatic_genes=gff.genes$average_methylation>heterochromatic_gene_coverage_cutoff | (gff.genes$average_methylation>heterochromatic_gene_coverage_minimum & gff.genes$variance_methylation/gff.genes$average_methylation<heterochromatic_gene_dispersion_cutoff)
	cat(paste0(nrow(gff.genes[heterochromatic_genes,])," genes are assigned as heterochromatic.\n"))
	# Schmitz 2627
	
	# Plot the distribution of average methylation by gene length for heterochromatic genes
	pdf(paste0(project_id,"_",meth_context,"_methylation_level_and heterochromatic_status.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=heterochromatic_genes) ,alpha=0.4, size=1) + theme_minimal() + xlab("Gene length") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="'Heterochromatic genes'")))
	dev.off()
		
	# Output the heterochromatic genes to a file for convenience
	write.table(cbind(gff.genes$gene_ID,gff.genes$average_methylation, heterochromatic_genes),file=paste0(project_id,"_",meth_context,"_heterochromatic_genes.tsv"),sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)

	cat(paste0(nrow(gff.genes[gff.genes$variable_count>0,])," of ",nrow(gff.genes)," genes have one or more sites which vary with concordance among replicates (",nrow(gff.genes[gff.genes$variable_count>0,])/nrow(gff.genes),").\n"))
	cat(paste0(nrow(gff.genes[gff.genes$variable_count>0 & heterochromatic_genes,])," of ",nrow(gff.genes[heterochromatic_genes,])," heterochromatic genes have one or more sites which vary with concordance among replicates (",nrow(gff.genes[gff.genes$variable_count>0 & heterochromatic_genes,])/nrow(gff.genes[heterochromatic_genes,]),").\n"))
	cat(paste0(nrow(gff.genes[gff.genes$variable_count>0 & heterochromatic_genes & gff.genes$average_methylation<=0.8,])," of ",nrow(gff.genes[heterochromatic_genes & gff.genes$average_methylation<=0.8,])," heterochromatic genes with average methylation <=0.8 have one or more sites which vary with concordance among replicates (",nrow(gff.genes[gff.genes$variable_count>0 & heterochromatic_genes & gff.genes$average_methylation<=0.8,])/nrow(gff.genes[heterochromatic_genes & gff.genes$average_methylation<=0.8,]),").\n"))
	cat(paste0(nrow(gff.genes[gff.genes$variable_count>0 & heterochromatic_genes & gff.genes$average_methylation>0.8,])," of ",nrow(gff.genes[heterochromatic_genes & gff.genes$average_methylation>0.8,])," heterochromatic genes with average methylation >0.8 have one or more sites which vary with concordance among replicates (",nrow(gff.genes[gff.genes$variable_count>0 & heterochromatic_genes & gff.genes$average_methylation>0.8,])/nrow(gff.genes[heterochromatic_genes & gff.genes$average_methylation>0.8,]),").\n"))
	
	cat(paste0(nrow(gff.transposons[gff.transposons$variable_count>0,])," of ",nrow(gff.transposons)," transposons have one or more sites which vary with concordance among replicates (",nrow(gff.transposons[gff.transposons$variable_count>0,])/nrow(gff.transposons),").\n"))
	
	
	# Check whether number of variable sites per gene is correlated with gene length and number of exons (yes to both)
	
	
	# Plot whether proportion of sites varying changes with gene length (not really)
	cat(paste0("Correlation between proportion of variable sites and length of gene in nt among GBM genes: r=",cor(gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$variable_count/(gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$variable_count+gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$all_M_count+gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$all_U_count), 1+gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$V5-gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$V4),"\n"))
	pdf(paste0(project_id,"_",meth_context,"_proportion_variable_by_GBM_gene_length.pdf"))
	print(ggplot(gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),], aes(x=1+V5-V4, y=variable_count/(variable_count+all_M_count+all_U_count)))  + geom_point(aes(colour=average_methylation) ,alpha=0.4, size=1) + theme_minimal() + scale_colour_gradient(low="green", high="red"))	
	dev.off()

	# Plot whether proportion of sites varying changes with methylation level of gene (yes it does)
	cat(paste0("Correlation between proportion of variable sites and average methylation level gene among GBM genes: r=",cor(gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$average_methylation, gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$variable_count/(gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$variable_count+gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$all_M_count+gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$all_U_count)),"\n"))
	
	# Plot whether proportion of sites varying changes with CG site density among GBM genes (yes it does - slightly negatively corellated)
	cat(paste0("Correlation between proportion of variable sites and CG site density among GBM genes: r=",cor(gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$average_methylation, gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$CG_sites_count/(1+gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$V5-gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$V4)),"\n"))
	pdf(paste0(project_id,"_",meth_context,"_proportion_variable_by_CG_site_density_GBM_genes.pdf"))
	print(ggplot(gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),], aes(x=CG_sites_count/(1+V5-V4), y=variable_count/(variable_count+all_M_count+all_U_count)))  + geom_point(aes(colour=(1+V5-V4)) ,alpha=0.4, size=1) + theme_minimal() + scale_colour_gradient(low="green", high="red", trans="log"))
	dev.off()
	
	
	# Assuming we are in the Schmitz project, read in the Becker heterochromatic genes analysis, and plot the correlation between the two studies' classifications
	becker_heterochromatic_genes=read.table(file=paste0("../PRJEB2678/PRJEB2678_",meth_context,"_heterochromatic_genes.tsv"),sep="\t")
	colnames(becker_heterochromatic_genes)[1:2]=c("gene_ID","average_methylation")
	pdf(paste0(project_id,"_",meth_context,"_average_methylation_by_gene_correlation.pdf"))
	print(ggplot(cbind(becker_heterochromatic_genes,cbind(gff.genes$gene_ID,gff.genes$average_methylation, heterochromatic_genes)), aes(x=becker_heterochromatic_genes$average_methylation, y=gff.genes$average_methylation)) +geom_point(aes(colour=gff.genes$variance_methylation/gff.genes$average_methylation), alpha=0.4, size=1) + scale_colour_gradient(low="red", high="green") + theme_minimal() +geom_hline(aes(yintercept=0.6))+geom_hline(aes(yintercept=0.8)) +geom_vline(aes(xintercept=0.6)) +geom_vline(aes(xintercept=0.8)) + guides(color=guide_legend(title="Dispersion index")))
	dev.off()
	
	pdf(paste0(project_id,"_",meth_context,"_average_methylation_by_gene_correlation_DI0.25.pdf"))
	print(ggplot(cbind(becker_heterochromatic_genes,cbind(gff.genes$gene_ID,gff.genes$average_methylation, heterochromatic_genes)), aes(x=becker_heterochromatic_genes$average_methylation, y=gff.genes$average_methylation)) +geom_point(aes(colour=gff.genes$variance_methylation/gff.genes$average_methylation>0.25), alpha=0.4, size=1)  + theme_minimal() +geom_hline(aes(yintercept=0.6))+geom_hline(aes(yintercept=0.8)) +geom_vline(aes(xintercept=0.6)) +geom_vline(aes(xintercept=0.8)) + guides(color=guide_legend(title="Dispersion index>0.25")))
	dev.off()
	

	
	# Set up a classifier column for genes to summarise each gene's methylation status.
	# Unmethylated means variance is low and average is within that expected from Chloroplast
	# Heterochromatic means variance is low and methylation level is high
	# Gene-body methylated means variance and average methylation level are within expectations
	# CG-poor means there are too few CG sites in the gene to make a determination
	# Unknown means that there is not enough coverage of the CG sites in the gene to make a determination
	
	gff.genes$m_class = ifelse(gff.genes$CG_sites_count<=1,"CG-poor",ifelse(gff.genes$all_M_count+gff.genes$all_U_count+gff.genes$variable_count==0,"Unknown",ifelse(gff.genes$average_methylation==0,"Unmethylated",ifelse((gff.genes$variance_methylation/gff.genes$average_methylation>heterochromatic_gene_dispersion_cutoff),"Gene-body Methylated",ifelse((gff.genes$variance_methylation/gff.genes$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.genes$average_methylation<heterochromatic_gene_coverage_minimum),"Unmethylated",ifelse((gff.genes$variance_methylation/gff.genes$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.genes$average_methylation>=heterochromatic_gene_coverage_minimum),"Heterochromatic","Really Unknown"))))))

	gff.transposons$m_class = ifelse(gff.transposons$CG_sites_count<=1,"CG-poor",ifelse(gff.transposons$all_M_count+gff.transposons$all_U_count+gff.transposons$variable_count==0,"Unknown",ifelse(gff.transposons$average_methylation==0,"Unmethylated",ifelse((gff.transposons$variance_methylation/gff.transposons$average_methylation>heterochromatic_gene_dispersion_cutoff),"Gene-body Methylated",ifelse((gff.transposons$variance_methylation/gff.transposons$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.transposons$average_methylation<heterochromatic_gene_coverage_minimum),"Unmethylated",ifelse((gff.transposons$variance_methylation/gff.transposons$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.transposons$average_methylation>=heterochromatic_gene_coverage_minimum),"Heterochromatic","Really Unknown"))))))

	saveRDS(gff.genes, file=paste0(project_id,"_",meth_context,"_gff.genes.rds"))
	saveRDS(gff.transposons, file=paste0(project_id,"_",meth_context,"_gff.transposons.rds"))
	saveRDS(gff.exons, file=paste0(project_id,"_",meth_context,"_gff.exons.rds"))
	saveRDS(gff.introns, file=paste0(project_id,"_",meth_context,"_gff.introns.rds"))
	saveRDS(gff.5UTR, file=paste0(project_id,"_",meth_context,"_gff.5UTR.rds"))
	saveRDS(gff.3UTR, file=paste0(project_id,"_",meth_context,"_gff.3UTR.rds"))

	
	cat(paste0("Correlation between number of variable sites and number of exons among GBM genes: r=",cor(gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$variable_count, gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$no_exons),"\n"))
	cat(paste0("Correlation between number of variable sites and length of gene in nt among GBM genes: r=",cor(gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$variable_count, 1+gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$V5-gff.genes[(gff.genes$m_class=="Gene-body Methylated") & !is.na(gff.genes$no_exons),]$V4),"\n"))
	# Schmitz:
	#Correlation between number of variable sites and number of exons among GBM genes: r=0.467013088180723
	#Correlation between number of variable sites and length of gene in nt among GBM genes: r=0.659174040597096

	
	# Plot genes coloured by classification
	pdf(paste0(project_id,"_",meth_context,"_gene_methylation_level_and_status.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=m_class) ,alpha=0.4, size=1) + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Gene methylation status\n")))
	dev.off()

	pdf(paste0(project_id,"_",meth_context,"_transposon_methylation_level_and_status.pdf"))
	print(ggplot(gff.transposons, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=m_class) ,alpha=0.4, size=1) + theme_minimal() + xlab("T`ransposon length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within transposon")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Transposon methylation status\n")))
	dev.off()

	# Summarise the breakdown of genes by classification
	cat(paste0("Counts of genes by classification:\n"))
	table(gff.genes$m_class)
	cat(paste0("Counts of transposons by classification:\n"))
	table(gff.transposons$m_class)

	# Count how many gene-body methylated genes have variable sites
	cat(paste0(nrow(gff.genes[(gff.genes$m_class=='Gene-body Methylated') & (gff.genes$variable_count>0),])," of ",nrow(gff.genes[gff.genes$m_class=='Gene-body Methylated',])," genes with gene body methylation (",nrow(gff.genes[(gff.genes$m_class=='Gene-body Methylated') & (gff.genes$variable_count>0),])/nrow(gff.genes[gff.genes$m_class=='Gene-body Methylated',]),") have one or more CG sites which vary.\n"))
	# Schmitz 12806 of 15264 genes with gene body methylation (0.83896750524109) have one or more CG sites which vary.
	
	cat(paste0(nrow(gff.genes[(gff.genes$m_class!='Gene-body Methylated') & (gff.genes$variable_count>0),])," of ",nrow(gff.genes[gff.genes$m_class!='Gene-body Methylated',])," genes without gene body methylation (",nrow(gff.genes[(gff.genes$m_class!='Gene-body Methylated') & (gff.genes$variable_count>0),])/nrow(gff.genes[gff.genes$m_class!='Gene-body Methylated',]),") have one or more CG sites which vary.\n"))
	# Schmitz 1385 of 18077 genes without gene body methylation (0.0766166952481053) have one or more CG sites which vary.
	
	cat(paste0(nrow(gff.genes[(gff.genes$m_class=='Unmethylated') & (gff.genes$variable_count>0),])," of ",nrow(gff.genes[gff.genes$m_class=='Unmethylated',])," unmethylated genes (",nrow(gff.genes[(gff.genes$m_class=='Unmethylated') & (gff.genes$variable_count>0),])/nrow(gff.genes[gff.genes$m_class=='Unmethylated',]),") have one or more CG sites which vary.\n"))
	# Schmitz 1138 of 14604 unmethylated genes (0.0779238564776773) have one or more CG sites which vary.
	
	cat(paste0(nrow(gff.transposons[(gff.transposons$m_class=='Gene-body Methylated') & (gff.transposons$variable_count>0),])," of ",nrow(gff.transposons[gff.transposons$m_class=='Gene-body Methylated',])," transposons with gene body methylation (",nrow(gff.transposons[(gff.transposons$m_class=='Gene-body Methylated') & (gff.transposons$variable_count>0),])/nrow(gff.transposons[gff.transposons$m_class=='Gene-body Methylated',]),") have one or more CG sites which vary.\n"))
	# Schmitz 527 of 2001 transposons with gene body methylation (0.263368315842079) have one or more CG sites which vary.
	
	cat(paste0(nrow(gff.transposons[(gff.transposons$m_class!='Gene-body Methylated') & (gff.transposons$variable_count>0),])," of ",nrow(gff.transposons[gff.transposons$m_class!='Gene-body Methylated',])," transposons without gene body methylation (",nrow(gff.transposons[(gff.transposons$m_class!='Gene-body Methylated') & (gff.transposons$variable_count>0),])/nrow(gff.transposons[gff.transposons$m_class!='Gene-body Methylated',]),") have one or more CG sites which vary.\n"))
	# Schmitz 1137 of 29188 transposons without gene body methylation (0.0389543648074551) have one or more CG sites which vary.
	
	cat(paste0(nrow(gff.transposons[(gff.transposons$m_class=='Unmethylated') & (gff.transposons$variable_count>0),])," of ",nrow(gff.transposons[gff.transposons$m_class=='Unmethylated',])," unmethylated transposons (",nrow(gff.transposons[(gff.transposons$m_class=='Unmethylated') & (gff.transposons$variable_count>0),])/nrow(gff.transposons[gff.transposons$m_class=='Unmethylated',]),") have one or more CG sites which vary.\n"))
	# Schmitz 185 of 3970 unmethylated transposons (0.0465994962216625) have one or more CG sites which vary.
	
	
	# Plot the proportion of variable sites against length for non-heterochromatic genes
	pdf(paste0(project_id,"_",meth_context,"_proportion_of_variable_sites_GBM_genes.pdf"))
	print(ggplot(gff.genes[gff.genes$m_class=="Gene-body Methylated",], aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variable_count/(variable_count+all_M_count+all_U_count)),alpha=0.4, size=1)  + scale_colour_gradient(low="green", high="red", trans="log", breaks=c(0.001,0.01,0.1,1)) + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Proportion of\nvariable sites\n")))
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_proportion_of_variable_sites_all_genes.pdf"))
	print(ggplot(gff.genes[,], aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variable_count/(variable_count+all_M_count+all_U_count)),alpha=0.4, size=1)  + scale_colour_gradient(low="green", high="red", trans="log", breaks=c(0.001,0.01,0.1,1)) + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Proportion of\nvariable sites\n")))
	dev.off()

	pdf(paste0(project_id,"_",meth_context,"_proportion_of_variable_sites_GBM_transposons.pdf"))
	print(ggplot(gff.transposons[gff.transposons$m_class=="Gene-body Methylated",], aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variable_count/(variable_count+all_M_count+all_U_count)),alpha=0.4, size=1)  + scale_colour_gradient(low="green", high="red", trans="log", breaks=c(0.001,0.01,0.1,1)) + theme_minimal() + xlab("Transposon length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within transposon")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Proportion of\nvariable sites\n")))
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_proportion_of_variable_sites_all_transposons.pdf"))
	print(ggplot(gff.transposons[,], aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variable_count/(variable_count+all_M_count+all_U_count)),alpha=0.4, size=1)  + scale_colour_gradient(low="green", high="red", trans="log", breaks=c(0.001,0.01,0.1,1)) + theme_minimal() + xlab("Transposon length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within transposon")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Proportion of\nvariable sites\n")))
	dev.off()

	
	### When CG and CHG are both done, plot the gene-wise correlation
	gff.genes.CG=readRDS(paste0(project_id,"_CG_gff.genes.rds"))
	gff.genes.CHG=readRDS(paste0(project_id,"_CHG_gff.genes.rds"))
	gff.transposons.CG=readRDS(paste0(project_id,"_CG_gff.transposons.rds"))
	gff.transposons.CHG=readRDS(paste0(project_id,"_CHG_gff.transposons.rds"))

	pdf(paste0(project_id,"_",meth_context,"_gene-wise_CHG_vs_CGmethylation.pdf"))
	print(ggplot(merge(gff.genes.CG, gff.genes.CHG, by="gene_ID"), aes(x=average_methylation.x, y=average_methylation.y)) + geom_point(aes(colour=m_class.x) ,alpha=0.4, size=1) + theme_minimal() + xlab("Average CG methylation") + ylab(paste0("Average CHG methylation")) + guides(color=guide_legend(title="Gene methylation status\n")))
	dev.off()
	
	pdf(paste0(project_id,"_",meth_context,"_transposon-wise_CHG_vs_CGmethylation.pdf"))
	print(ggplot(merge(gff.transposons.CG, gff.transposons.CHG, by="gene_ID"), aes(x=average_methylation.x, y=average_methylation.y)) + geom_point(aes(colour=m_class.x) ,alpha=0.4, size=1) + theme_minimal() + xlab("Average CG methylation") + ylab(paste0("Average CHG methylation")) + guides(color=guide_legend(title="Transposon methylation status\n")))
	dev.off()
	
	pdf(paste0(project_id,"_",meth_context,"_gene-wise_CHG_vs_CG_variable_sites.pdf"))
	print(ggplot(merge(gff.genes.CG, gff.genes.CHG, by="gene_ID"), aes(x=variable_count.x/(variable_count.x+all_M_count.x+all_U_count.x),y=average_methylation.y)) + geom_point(aes(colour=variance_methylation.x/average_methylation.x) ,alpha=0.4, size=1) + scale_colour_gradient(low="red", high="green") + theme_minimal() + xlab("Proportion of variable CG sites in gene") + ylab(paste0("Average methylation of CHG sites within gene")) +guides(color=guide_legend(title="Dispersion index\n of CG methylation\n in gene")))
	dev.off()

	pdf(paste0(project_id,"_",meth_context,"_transposon-wise_CHG_vs_CG_variable_sites.pdf"))
	print(ggplot(merge(gff.transposons.CG, gff.transposons.CHG, by="gene_ID"), aes(x=variable_count.x/(variable_count.x+all_M_count.x+all_U_count.x),y=average_methylation.y)) + geom_point(aes(colour=variance_methylation.x/average_methylation.x) ,alpha=0.4, size=1) + scale_colour_gradient(low="red", high="green") + theme_minimal() + xlab("Proportion of variable CG sites in transposon") + ylab(paste0("Average methylation of CHG sites within transposon")) +guides(color=guide_legend(title="Dispersion index\n of CG methylation\n in transposon")))
	dev.off()

	
	# We now want to use CHG methylation to identify a set of additional genes, previously identified as GBM, which we will redefine as Heterochromatic on the basis of CHG methylation.
	### This section is left in for reference, but is superceded by a later analysis.
	### Overall average mCHG over the length of the gene is not necessarily the best way to exclude GBM as the explanation for dispersion of mCG.  Instead, later after we run a segmentation model, we identify arbitrary mCG segments of the genome where mCHG is higher than a threshold defined by comparison with heterochromatic genes and transposons.
	### Neverthelsess, onwards for now:
	
	# First we plot the distributions of CHG methylation level for each class of gene
	pdf(paste0(project_id,"_",meth_context,"_gene_wise_CHG_methylation_density_by_m_class.pdf"))
	print(ggplot(merge(gff.genes.CG, gff.genes.CHG, by="gene_ID")[,], aes(x=average_methylation.y, colour=m_class.x)) + geom_histogram(aes(y=0.1*..density..),binwidth=0.1)+ theme_minimal() + xlab("Average CG methylation") + xlab(paste0("Average CHG methylation")) + guides(color=guide_legend(title="Gene CG methylation status\n")) + facet_grid(m_class.x ~ .))
	dev.off()
	
	# Then we fit a mixture of 2 Gaussians to the combined density of methylation levels of all genes and find where they intersect.  This is where we draw the divide between GBM and Heterochromatic genes
	
	# Fit a mixture of 2 Gaussians to the dispersion index density
	CHG_meth_mixmdl = normalmixEM(merge(gff.genes.CG, gff.genes.CHG[(!is.na(gff.genes.CHG$average_methylation)) & (gff.genes.CHG$average_methylation>0),], by="gene_ID", all.y=TRUE)$average_methylation.y, k=2)
	# Had to exclude sites with zero CHG methylation of the fitting algorithm exploded due to too little variance in one class
	
	
	
	# Print the characteristics of the fit curves
	cat(paste0("Fit lambda values (mixture proportions):\n"))
	CHG_meth_mixmdl$lambda
	# Schmitz data:  0.3730093 0.6269907
	# Becker data: 
	
	cat(paste0("Fit mu values (means):\n"))
	CHG_meth_mixmdl$mu
	# Schmitz data: 0.0005150998 0.3662588779
	# Becker data:  
	
	cat(paste0("Fit sigma values (st.devs):\n"))
	CHG_meth_mixmdl$sigma
	# Schmitz data: 0.0005349217 0.3168317368
	# Becker data:  
	
	# Find the intersection of the two larger components in the mixture - here we cut off

	#packages were installed and mixfn defined earlier
	#install.packages("rootSolve")
	#library(rootSolve)
	#mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
	#dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }

	CHG_meth_low_model=1
	CHG_meth_high_model=low_model+1
	CHG_meth_cutoff=uniroot.all(mixfn, lower=0, upper=1, m1=CHG_meth_mixmdl$mu[CHG_meth_low_model], sd1=CHG_meth_mixmdl$sigma[CHG_meth_low_model], m2=CHG_meth_mixmdl$mu[CHG_meth_high_model], sd2=CHG_meth_mixmdl$sigma[CHG_meth_high_model], p1=CHG_meth_mixmdl$lambda[CHG_meth_low_model], p2=CHG_meth_mixmdl$lambda[CHG_meth_high_model])
	# Schmitz CG data: 0.002423706
	# Becker CG data: 
	
	# In the Schmitz data, cutoff of 0.002423706 removes 1327 GBM genes from the original set of 15095 (8.8%)

	pdf(paste0(project_id,"_",meth_context,"_CHG_methylation_GBM_cutoff_density_fitted.pdf"))
	print(plot(CHG_meth_mixmdl,which=2))
	print(lines(density(merge(gff.genes.CG, gff.genes.CHG[(!is.na(gff.genes.CHG$average_methylation)) & (gff.genes.CHG$average_methylation>0),], by="gene_ID", all.y=TRUE)$average_methylation.y), lty=2, lwd=2))
	dev.off()

	#  This may be a bit much (.2% is a very severe cutoff), so we try instead to look at 1%ile of the CHG methylation among heterochromatic genes
	CHG_meth_quantile_cutoff = 0.01
	CHG_meth_cutoff=quantile(merge(gff.genes.CG[gff.genes.CG$m_class=="Heterochromatic",], gff.genes.CHG[(!is.na(gff.genes.CHG$average_methylation)) & (gff.genes.CHG$average_methylation>0),], by="gene_ID")$average_methylation.y, CHG_meth_quantile_cutoff, na.rm=TRUE)	
	cat(paste0(nrow(merge(gff.genes.CG[gff.genes.CG$m_class=="Gene-body Methylated",], gff.genes.CHG[(!is.na(gff.genes.CHG$average_methylation)) & (gff.genes.CHG$average_methylation>=CHG_meth_cutoff),], by="gene_ID"))," 'GBM' genes have CHG methylation above the ",CHG_meth_quantile_cutoff," quantile of CHG methylation among 'Heterochromatic' genes (",CHG_meth_cutoff,").\n"))
	# Schmitz data:
	# quantile	mCHGcutoff	GBM genes
	# .01	.011	1046
	# .02	.048	681
	# .03	.088	505
	# .04	.118	402
	# .05	.156	315
	

	# Let's look instead at individual CG methylated segments from the segmentation model.  This is done below in MethylSeekR section.
	

	
	# How many CG, variable and invariable sites lie within annotated genes?
	### These summary steps are commented out for now - they double-count for sites that lie within more than one gene model
	#cat(paste(sum(gff.genes$CG_sites_count),"of",nrow(all_samples_meth_status),"(",sum(gff.genes$CG_sites_count)/nrow(all_samples_meth_status),") CG sites lie within annotated genes.\n"))
	#cat(paste(sum(gff.genes$all_M_count),"of",nrow(all_samples_meth_status[all_m==TRUE,]),"(",sum(gff.genes$all_M_count)/nrow(all_samples_meth_status[all_m==TRUE,]),") 'all M' sites lie within annotated genes.\n"))
	#cat(paste(sum(gff.genes$all_U_count),"of",nrow(all_samples_meth_status[all_u==TRUE,]),"(",sum(gff.genes$all_U_count)/nrow(all_samples_meth_status[all_u==TRUE,]),") 'all U' sites lie within annotated genes.\n"))
	#cat(paste(sum(gff.genes$variable_count),"of",nrow(variant_calls),"(",sum(gff.genes$variable_count)/nrow(all_samples_meth_status[(all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE),]),") variable sites lie within annotated genes.\n"))
	#cat(paste(sum(gff.genes$m_to_u_count),"of",nrow(all_samples_meth_status[(all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3==TRUE),]),"(",sum(gff.genes$m_to_u_count)/nrow(all_samples_meth_status[(all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3==TRUE),]),") M->U transition sites lie within annotated genes.\n"))
	#cat(paste(sum(gff.genes$u_to_m_count),"of",nrow(all_samples_meth_status[(all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3==TRUE),]),"(",sum(gff.genes$u_to_m_count)/nrow(all_samples_meth_status[(all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3==TRUE),]),") U->M transition sites lie within annotated genes.\n"))
	
	# Plot distribution of M->U and U->M calls per gene 
	sp=ggscatter(gff.genes, x="u_to_m_count", y="m_to_u_count", colour="all_M_count", size=2, alpha=0.6, palette="jco")
	xplot=ggdensity(gff.genes, "u_to_m_count", fill="all_M_count", palette="jco")
	yplot=ggdensity(gff.genes, "m_to_u_count", fill="all_M_count", palette="jco") + rotate()
	yplot=yplot+clean_theme()
	xplot=xplot+clean_theme()
	pdf(paste0(project_id,"_",meth_context,"_Mstatus_transitions_per_gene.pdf"))
	ggarrange(xplot,NULL,sp,yplot,ncol=2, nrow=2, align="hv", widths=c(2,1), heights=c(1,2))
	dev.off()
	
	# Histogram instead:
	pdf(paste0(project_id,"_",meth_context,"_Mstatus_transitions_per_GBM_gene_var_gt_3.pdf"))
	print(ggplot(gff.genes[gff.genes$m_class=="Gene-body Methylated" & gff.genes$variable_count>3,], aes(x=m_to_u_count/(u_to_m_count+m_to_u_count))) + geom_histogram(aes(y= ..count..), binwidth=0.1) + theme_minimal())
	dev.off()
	
	### In the case of transposons we would expect all the calls to go in one direction in an individual transposon, whereas in the case of genes we might expect a (weak?) correlation between M->U and U->M calls per gene
	# Plot distribution of M->U and U->M calls per transposon
	sp=ggscatter(gff.transposons, x="u_to_m_count", y="m_to_u_count", colour="all_M_count", size=2, alpha=0.6, palette="jco")
	xplot=ggdensity(gff.transposons, "u_to_m_count", fill="all_M_count", palette="jco")
	yplot=ggdensity(gff.transposons, "m_to_u_count", fill="all_M_count", palette="jco") + rotate()
	yplot=yplot+clean_theme()
	xplot=xplot+clean_theme()
	pdf(paste0(project_id,"_",meth_context,"_Mstatus_transitions_per_transposon.pdf"))
	ggarrange(xplot,NULL,sp,yplot,ncol=2, nrow=2, align="hv", widths=c(2,1), heights=c(1,2))
	dev.off()
	
	# Histogram instead:
	pdf(paste0(project_id,"_",meth_context,"_Mstatus_transitions_per_transposon_var_gt_3.pdf"))
	print(ggplot(gff.transposons[gff.transposons$variable_count>3,], aes(x=m_to_u_count/(u_to_m_count+m_to_u_count))) + geom_histogram(aes(y= ..count..), binwidth=0.1) + theme_minimal())
	dev.off()
	
	
	### So far we have considered the characteristics of genes, transposons and other objects in terms of the CG sites overlapping them.  We now go on to consider the characteristics of CG sites in gene space.
	### Many CG sites are double-counted where gene models overlap.  Need to collapse gene space into a compact set of GenomeRanges, then redo the overlap analyses to count numbers of sites in each class within gene space

	# Make collapsed GRanges for gene space (reduce() collapses overlapping granges into a single grange)
	gene_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.genes[(!is.na(heterochromatic_genes)) & !heterochromatic_genes,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	transposon_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.transposons[,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	exon_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.exons[gff.exons$gene_ID %in% gff.genes[(!is.na(heterochromatic_genes)) & !heterochromatic_genes,]$gene_ID,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	intron_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.introns[gff.introns$gene_ID %in% gff.genes[(!is.na(heterochromatic_genes)) & !heterochromatic_genes,]$gene_ID,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	UTR5_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.5UTR[gff.5UTR$gene_ID %in% gff.genes[(!is.na(heterochromatic_genes)) & !heterochromatic_genes,]$gene_ID,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	UTR3_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.3UTR[gff.3UTR$gene_ID %in% gff.genes[(!is.na(heterochromatic_genes)) & !heterochromatic_genes,]$gene_ID,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))

	cat(paste0("The length of non-heterochromatic genes is ",sum(width(gene_ranges_reduced)),"nt\n"))
	cat(paste0("The length of non-heterochromatic exons is ",sum(width(exon_ranges_reduced)),"nt\n"))
	cat(paste0("The length of non-heterochromatic introns is ",sum(width(intron_ranges_reduced)),"nt\n"))
	
	# Alternative version just for gene-body methylated genes
	gene_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.genes[gff.genes$m_class=="Gene-body Methylated",], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	transposon_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.transposons[,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	exon_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.exons[gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	intron_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.introns[gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	UTR5_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.5UTR[gff.5UTR$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))
	UTR3_ranges_reduced=reduce(makeGRangesFromDataFrame(df = gff.3UTR[gff.3UTR$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,], start.field = "V4", end.field = "V5",seqnames.field = "V1"))

	cat(paste0("The length of Gene-body methylated genes is ",sum(width(gene_ranges_reduced)),"nt\n"))
	cat(paste0("The length of exons in gene-body methylated genes is ",sum(width(exon_ranges_reduced)),"nt\n"))
	cat(paste0("The length of introns in gene-body methylated genes is ",sum(width(intron_ranges_reduced)),"nt\n"))
	
	
	# Find the overlaps
	olaps_CG_sites_reduced = findOverlaps(CG_site_ranges, gene_ranges_reduced)
	olaps_all_M_reduced = findOverlaps(all_M_ranges, gene_ranges_reduced)
	olaps_all_U_reduced = findOverlaps(all_U_ranges, gene_ranges_reduced)
	olaps_variable_reduced = findOverlaps(variant_call_ranges, gene_ranges_reduced)
	olaps_m_to_u_reduced = findOverlaps(m_to_u_call_ranges, gene_ranges_reduced)
	olaps_u_to_m_reduced = findOverlaps(u_to_m_call_ranges, gene_ranges_reduced)

	tolaps_CG_sites_reduced = findOverlaps(CG_site_ranges, transposon_ranges_reduced)
	tolaps_all_M_reduced = findOverlaps(all_M_ranges, transposon_ranges_reduced)
	tolaps_all_U_reduced = findOverlaps(all_U_ranges, transposon_ranges_reduced)
	tolaps_variable_reduced = findOverlaps(variant_call_ranges, transposon_ranges_reduced)
	tolaps_m_to_u_reduced = findOverlaps(m_to_u_call_ranges, transposon_ranges_reduced)
	tolaps_u_to_m_reduced = findOverlaps(u_to_m_call_ranges, transposon_ranges_reduced)

	eolaps_CG_sites_reduced = findOverlaps(CG_site_ranges, exon_ranges_reduced)
	eolaps_all_M_reduced = findOverlaps(all_M_ranges, exon_ranges_reduced)
	eolaps_all_U_reduced = findOverlaps(all_U_ranges, exon_ranges_reduced)
	eolaps_variable_reduced = findOverlaps(variant_call_ranges, exon_ranges_reduced)
	eolaps_m_to_u_reduced = findOverlaps(m_to_u_call_ranges, exon_ranges_reduced)
	eolaps_u_to_m_reduced = findOverlaps(u_to_m_call_ranges, exon_ranges_reduced)

	iolaps_CG_sites_reduced = findOverlaps(CG_site_ranges, intron_ranges_reduced)
	iolaps_all_M_reduced = findOverlaps(all_M_ranges, intron_ranges_reduced)
	iolaps_all_U_reduced = findOverlaps(all_U_ranges, intron_ranges_reduced)
	iolaps_variable_reduced = findOverlaps(variant_call_ranges, intron_ranges_reduced)
	iolaps_m_to_u_reduced = findOverlaps(m_to_u_call_ranges, intron_ranges_reduced)
	iolaps_u_to_m_reduced = findOverlaps(u_to_m_call_ranges, intron_ranges_reduced)

	olaps5_CG_sites_reduced = findOverlaps(CG_site_ranges, UTR5_ranges_reduced)
	olaps5_all_M_reduced = findOverlaps(all_M_ranges, UTR5_ranges_reduced)
	olaps5_all_U_reduced = findOverlaps(all_U_ranges, UTR5_ranges_reduced)
	olaps5_variable_reduced = findOverlaps(variant_call_ranges, UTR5_ranges_reduced)
	olaps5_m_to_u_reduced = findOverlaps(m_to_u_call_ranges, UTR5_ranges_reduced)
	olaps5_u_to_m_reduced = findOverlaps(u_to_m_call_ranges, UTR5_ranges_reduced)

	olaps3_CG_sites_reduced = findOverlaps(CG_site_ranges, UTR3_ranges_reduced)
	olaps3_all_M_reduced = findOverlaps(all_M_ranges, UTR3_ranges_reduced)
	olaps3_all_U_reduced = findOverlaps(all_U_ranges, UTR3_ranges_reduced)
	olaps3_variable_reduced = findOverlaps(variant_call_ranges, UTR3_ranges_reduced)
	olaps3_m_to_u_reduced = findOverlaps(m_to_u_call_ranges, UTR3_ranges_reduced)
	olaps3_u_to_m_reduced = findOverlaps(u_to_m_call_ranges, UTR3_ranges_reduced)

	# How many CG, variable and invariable sites lie within annotated genes?
	cat(paste(length(olaps_CG_sites_reduced),"of",nrow(all_samples_meth_status),"(",length(olaps_CG_sites_reduced)/nrow(all_samples_meth_status),") ",meth_context," sites lie within annotated non-heterochromatic genes.\n"))
	cat(paste(length(olaps_all_M_reduced),"of",nrow(all_samples_meth_status[all_m==TRUE,]),"(",length(olaps_all_M_reduced)/nrow(all_samples_meth_status[all_m==TRUE,]),") 'all M' sites lie within annotated non-heterochromatic genes.\n"))
	cat(paste(length(olaps_all_U_reduced),"of",nrow(all_samples_meth_status[all_u==TRUE,]),"(",length(olaps_all_U_reduced)/nrow(all_samples_meth_status[all_u==TRUE,]),") 'all U' sites lie within annotated non-heterochromatic genes.\n"))
	cat(paste(length(olaps_variable_reduced),"of",nrow(variant_calls),"(",length(olaps_variable_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals,]),") variable sites lie within annotated non-heterochromatic genes.\n"))
	cat(paste(length(olaps_m_to_u_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),,]),"(",length(olaps_m_to_u_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),") M->U transition sites lie within annotated non-heterochromatic genes.\n"))
	cat(paste(length(olaps_u_to_m_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),,]),"(",length(olaps_u_to_m_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),") U->M transition sites lie within annotated non-heterochromatic genes.\n"))
	# Schmitz:
	#1085658 of 2802193 ( 0.387431558068984 )  CG  sites lie within annotated non-heterochromatic genes.
	#113905 of 592062 ( 0.19238694596174 ) 'all M' sites lie within annotated non-heterochromatic genes.
	#412585 of 898177 ( 0.45935823339943 ) 'all U' sites lie within annotated non-heterochromatic genes.
	#57641 of 70768 ( 0.814506556635768 ) variable sites lie within annotated non-heterochromatic genes.
	#31706 of 35852 ( 0.884357915876381 ) M->U transition sites lie within annotated non-heterochromatic genes.
	#25935 of 34916 ( 0.723390605823943 ) U->M transition sites lie within annotated non-heterochromatic genes.

	# How many CG, variable and invariable sites lie within annotated transposons?
	cat(paste(length(tolaps_CG_sites_reduced),"of",nrow(all_samples_meth_status),"(",length(tolaps_CG_sites_reduced)/nrow(all_samples_meth_status),") ",meth_context," sites lie within annotated transposons.\n"))
	cat(paste(length(tolaps_all_M_reduced),"of",nrow(all_samples_meth_status[all_m==TRUE,]),"(",length(tolaps_all_M_reduced)/nrow(all_samples_meth_status[all_m==TRUE,]),") 'all M' sites lie within annotated transposons.\n"))
	cat(paste(length(tolaps_all_U_reduced),"of",nrow(all_samples_meth_status[all_u==TRUE,]),"(",length(tolaps_all_U_reduced)/nrow(all_samples_meth_status[all_u==TRUE,]),") 'all U' sites lie within annotated transposons.\n"))
	cat(paste(length(tolaps_variable_reduced),"of",nrow(variant_calls),"(",length(tolaps_variable_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals,]),") variable sites lie within annotated transposons.\n"))
	cat(paste(length(tolaps_m_to_u_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),"(",length(tolaps_m_to_u_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),") M->U transition sites lie within annotated transposons.\n"))
	cat(paste(length(tolaps_u_to_m_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),"(",length(tolaps_u_to_m_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),") U->M transition sites lie within annotated transposons.\n"))
	# Schmitz:
	#516901 of 2802193 ( 0.184463025922911 )  CG  sites lie within annotated transposons.
	#350754 of 592062 ( 0.592427820059386 ) 'all M' sites lie within annotated transposons.
	#24058 of 898177 ( 0.0267853663587467 ) 'all U' sites lie within annotated transposons.
	#2820 of 70768 ( 0.0398485191046801 ) variable sites lie within annotated transposons.
	#1218 of 35852 ( 0.0339730001115698 ) M->U transition sites lie within annotated transposons.
	#1602 of 34916 ( 0.045881544277695 ) U->M transition sites lie within annotated transposons.
	
	# How many CG, variable and invariable sites lie within annotated exons?
	cat(paste(length(eolaps_CG_sites_reduced),"of",nrow(all_samples_meth_status),"(",length(eolaps_CG_sites_reduced)/nrow(all_samples_meth_status),") ",meth_context," sites lie within annotated exons.\n"))
	cat(paste(length(eolaps_all_M_reduced),"of",nrow(all_samples_meth_status[all_m==TRUE,]),"(",length(eolaps_all_M_reduced)/nrow(all_samples_meth_status[all_m==TRUE,]),") 'all M' sites lie within annotated exons.\n"))
	cat(paste(length(eolaps_all_U_reduced),"of",nrow(all_samples_meth_status[all_u==TRUE,]),"(",length(eolaps_all_U_reduced)/nrow(all_samples_meth_status[all_u==TRUE,]),") 'all U' sites lie within annotated exons.\n"))
	cat(paste(length(eolaps_variable_reduced),"of",nrow(variant_calls),"(",length(eolaps_variable_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals,]),") variable sites lie within annotated exons.\n"))
	cat(paste(length(eolaps_m_to_u_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),"(",length(eolaps_m_to_u_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),") M->U transition sites lie within annotated exons.\n"))
	cat(paste(length(eolaps_u_to_m_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),"(",length(eolaps_u_to_m_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),") U->M transition sites lie within annotated exons.\n"))

	# How many CG, variable and invariable sites lie within annotated introns?
	cat(paste(length(iolaps_CG_sites_reduced),"of",nrow(all_samples_meth_status),"(",length(iolaps_CG_sites_reduced)/nrow(all_samples_meth_status),") ",meth_context," sites lie within annotated introns.\n"))
	cat(paste(length(iolaps_all_M_reduced),"of",nrow(all_samples_meth_status[all_m==TRUE,]),"(",length(iolaps_all_M_reduced)/nrow(all_samples_meth_status[all_m==TRUE,]),") 'all M' sites lie within annotated introns.\n"))
	cat(paste(length(iolaps_all_U_reduced),"of",nrow(all_samples_meth_status[all_u==TRUE,]),"(",length(iolaps_all_U_reduced)/nrow(all_samples_meth_status[all_u==TRUE,]),") 'all U' sites lie within annotated introns.\n"))
	cat(paste(length(iolaps_variable_reduced),"of",nrow(variant_calls),"(",length(iolaps_variable_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals,]),") variable sites lie within annotated introns.\n"))
	cat(paste(length(iolaps_m_to_u_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),"(",length(iolaps_m_to_u_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),") M->U transition sites lie within annotated introns.\n"))
	cat(paste(length(iolaps_u_to_m_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),"(",length(iolaps_u_to_m_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),") U->M transition sites lie within annotated introns.\n"))

	# How many CG, variable and invariable sites lie within annotated 5'UTRs?
	cat(paste(length(olaps5_CG_sites_reduced),"of",nrow(all_samples_meth_status),"(",length(olaps5_CG_sites_reduced)/nrow(all_samples_meth_status),") ",meth_context," sites lie within annotated 5'UTRs.\n"))
	cat(paste(length(olaps5_all_M_reduced),"of",nrow(all_samples_meth_status[all_m==TRUE,]),"(",length(olaps5_all_M_reduced)/nrow(all_samples_meth_status[all_m==TRUE,]),") 'all M' sites lie within annotated 5'UTRs.\n"))
	cat(paste(length(olaps5_all_U_reduced),"of",nrow(all_samples_meth_status[all_u==TRUE,]),"(",length(olaps5_all_U_reduced)/nrow(all_samples_meth_status[all_u==TRUE,]),") 'all U' sites lie within annotated 5'UTRs.\n"))
	cat(paste(length(olaps5_variable_reduced),"of",nrow(variant_calls),"(",length(olaps5_variable_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals,]),") variable sites lie within annotated 5'UTRs.\n"))
	cat(paste(length(olaps5_m_to_u_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),"(",length(olaps5_m_to_u_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),") M->U transition sites lie within annotated 5'UTRs.\n"))
	cat(paste(length(olaps5_u_to_m_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),"(",length(olaps5_u_to_m_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),") U->M transition sites lie within annotated 5'UTRs.\n"))

	# How many CG, variable and invariable sites lie within annotated 3'UTRs?
	cat(paste(length(olaps3_CG_sites_reduced),"of",nrow(all_samples_meth_status),"(",length(olaps3_CG_sites_reduced)/nrow(all_samples_meth_status),") ",meth_context," sites lie within annotated 3'UTRs.\n"))
	cat(paste(length(olaps3_all_M_reduced),"of",nrow(all_samples_meth_status[all_m==TRUE,]),"(",length(olaps3_all_M_reduced)/nrow(all_samples_meth_status[all_m==TRUE,]),") 'all M' sites lie within annotated 3'UTRs.\n"))
	cat(paste(length(olaps3_all_U_reduced),"of",nrow(all_samples_meth_status[all_u==TRUE,]),"(",length(olaps3_all_U_reduced)/nrow(all_samples_meth_status[all_u==TRUE,]),") 'all U' sites lie within annotated 3'UTRs.\n"))
	cat(paste(length(olaps3_variable_reduced),"of",nrow(variant_calls),"(",length(olaps3_variable_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals,]),") variable sites lie within annotated 3'UTRs.\n"))
	cat(paste(length(olaps3_m_to_u_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),"(",length(olaps3_m_to_u_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),") M->U transition sites lie within annotated 3'UTRs.\n"))
	cat(paste(length(olaps3_u_to_m_reduced),"of",nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),"(",length(olaps3_u_to_m_reduced)/nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),") U->M transition sites lie within annotated 3'UTRs.\n"))

	
	### These two lines belong later in the code:
	cat(paste0(length(unique(as.data.frame(findOverlaps(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"), variant_call_ranges))$queryHits))," unique GBM genes contain at least one mCG variable site (",100*length(unique(as.data.frame(findOverlaps(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"), variant_call_ranges))$queryHits))/nrow(gff.genes[gff.genes$m_class=="Gene-body Methylated",]),"%)\n"))
	# Schmitz data after filtering variable sites for overlaps with DMRs:
	#12806 unique GBM genes contain at least one mCG variable site (83.896750524109%)
	
	cat(paste0(sum(segment_annotation_summary[segment_annotation_summary$annotation_type=="GBM genes",]$variable_sites)," mCG variable sites lie within GBM genes (",100*sum(segment_annotation_summary[segment_annotation_summary$annotation_type=="GBM genes",]$variable_sites)/length(variant_call_ranges),"%)\n"))
	# Schmitz data after filtering variable sites for overlaps with DMRs:
	#57641 mCG variable sites lie within GBM genes (87.7857480087114%)

	
	
	
	
	
	# Plot net changes per GBM gene vs. total changes
	pdf(paste0(project_id,"_",meth_context,"_net_changes_in_gene_bodies.pdf"))
	print(ggplot(gff.genes[gff.genes$m_class=="Gene-body Methylated",], aes(x=as.factor(m_to_u_count+u_to_m_count), y=u_to_m_count-m_to_u_count)) + geom_boxplot(varwidth=TRUE, alpha = 0.2) + theme_minimal() +xlab("Number of mCG changes in gene body") + ylab("Net change in number of mCG sites in gene body"))
	dev.off()

	# Does number of changes per gene correlate with number of sites per gene?
	cor(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$CG_sites_count, gff.genes[gff.genes$m_class=="Gene-body Methylated",]$variable_count)
	# Schmitz: 0.6431103

	# Plot total changes per gene and net changes histograms/freqpolys
	pdf(paste0(project_id,"_",meth_context,"_net_and_total_changes_in_gene_bodies.pdf"))
	print(ggplot(rbind.data.frame(cbind.data.frame("count"=gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count+gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count,"metric"=rep("Number of sites where mCG changes",nrow(gff.genes[gff.genes$m_class=="Gene-body Methylated",]))), cbind("count"=abs(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count-gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count),"metric"=rep("Net change in number of mCG sites",nrow(gff.genes[gff.genes$m_class=="Gene-body Methylated",])))), aes(x=as.numeric(count), fill=metric)) + geom_histogram(position="dodge", binwidth=1) + theme_minimal() +xlab("mCG changes in gene body") + ylab("Number of gene models") +xlim(0,32))
	print(ggplot(rbind.data.frame(cbind.data.frame("count"=gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count+gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count,"metric"=rep("Number of sites where mCG changes",nrow(gff.genes[gff.genes$m_class=="Gene-body Methylated",]))), cbind("count"=abs(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count-gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count),"metric"=rep("Net change in number of mCG sites",nrow(gff.genes[gff.genes$m_class=="Gene-body Methylated",])))), aes(x=as.numeric(count), fill=metric)) + geom_density(bw=1, alpha=0.2) + theme_minimal() +xlab("mCG changes in gene body") + ylab("Density") +xlim(0,32))
	dev.off()
	
	### Exon/intron analysis  
	
	pdf(paste0(project_id,"_",meth_context,"_GBM_introns_exons_variation.pdf"))
	# Within exons, is there a correlation between exon number and number of variable sites? (spoiler: NO)
	print(ggplot(gff.exons, aes(x=as.factor(exon_no), y=variable_count))  + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Variable sites per exon, all genes"))
	cor(gff.exons$exon_no, gff.exons$variable_count)
	# Schmitz r=0.04615776
	
	# Is the same true for introns? (spoiler: YES)
	print(ggplot(gff.introns, aes(x=as.factor(exon_no), y=variable_count))  + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Variable sites per intron, all genes"))
	cor(gff.introns$exon_no, gff.introns$variable_count)	
	# Schmitz r=0.0792817
	
	# Does the number of CG sites increase with exon Number? (decreases actually, but in exons, not introns):
	exons_CG_sites_model=lm(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$CG_sites_count ~ gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	print(ggplot(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=CG_sites_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_abline(intercept = coef(exons_CG_sites_model)[1], slope = coef(exons_CG_sites_model)[2], color="black", size=.1) + ggtitle("CG sites per exon, GBM genes"))
	cat(paste0("Relationship between number of CG sites and exon no in exons:\n"))
	coef(exons_CG_sites_model)
	cat(paste0("r=",cor(gff.exons$exon_no, gff.exons$CG_sites_count),"\n"))
	# Schmitz: Intercept 11.9973933, Slope -0.3651917, -0.208990998917695
	
	introns_CG_sites_model=lm(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$CG_sites_count ~ gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	print(ggplot(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=CG_sites_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_abline(intercept = coef(introns_CG_sites_model)[1], slope = coef(introns_CG_sites_model)[2], color="black", size=.1) + ggtitle("CG sites per intron, GBM genes"))
	cat(paste0("Relationship between number of CG sites and exon no in introns:\n"))
	coef(introns_CG_sites_model)
	cat(paste0("r=",cor(gff.introns$exon_no, gff.introns$CG_sites_count),"\n"))
	# Schmitz: Intercept 5.49547183, Slope 0.02744224, r=0.0212447963946072

	
	# Does the density, as well as the number, of CG sites decrease with increasing exon number? (YES in exons but not introns)
	exons_CG_density_model=lm(gff.exons[(gff.exons$V5-gff.exons$V4>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$CG_sites_count/(gff.exons[(gff.exons$V5-gff.exons$V4>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$V5-gff.exons[(gff.exons$V5-gff.exons$V4>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$V4+1) ~ gff.exons[(gff.exons$V5-gff.exons$V4>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	print(ggplot(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=CG_sites_count/(V5-V4+1))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_abline(intercept = coef(exons_CG_density_model)[1], slope = coef(exons_CG_density_model)[2], color="black", size=.1) + ggtitle("CG site density per exon, GBM genes"))
	cat(paste0("Relationship between CG site density and exon no in exons:\n"))
	coef(exons_CG_density_model)
	cat(paste0("r=",cor(gff.exons$exon_no, gff.exons$CG_sites_count/(gff.exons$V5-gff.exons$V4+1)),"\n"))
	# Schmitz: Intercept 0.0317658269, Slope -0.0004554816, r=-0.186781022106755

	introns_CG_density_model=lm(gff.introns[(gff.introns$V5-gff.introns$V4>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$CG_sites_count/(gff.introns[(gff.introns$V5-gff.introns$V4>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$V5-gff.introns[(gff.introns$V5-gff.introns$V4>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$V4+1) ~ gff.introns[(gff.introns$V5-gff.introns$V4>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	print(ggplot(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=CG_sites_count/(V5-V4+1))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_abline(intercept = coef(introns_CG_density_model)[1], slope = coef(introns_CG_density_model)[2], color="black", size=.1) + ggtitle("CG site density per intron, GBM genes"))
	cat(paste0("Relationship between CG site density and exon no in introns:\n"))
	coef(introns_CG_density_model)
	cat(paste0("r=",cor(gff.introns$exon_no, gff.introns$CG_sites_count/(gff.introns$V5-gff.introns$V4+1)),"\n"))
	# Schmitz: Intercept 1.360824e-02, Slope -7.389975e-05, r=-0.0586177065072982

	
	# Length of exons and introns in nt do also show these correlations as strongly, but still do:
	cor(gff.exons$V5-gff.exons$V4+1, gff.exons$exon_no)
	print(ggplot(gff.exons, aes(x=exon_no, y=V5-V4+1)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Exon length, all genes"))
	# Schmitz r=-0.1609663
	
	cor(gff.introns$V5-gff.introns$V4+1, gff.introns$exon_no)
	print(ggplot(gff.introns, aes(x=exon_no, y=V5-V4+1)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Intron length, all genes"))
	# Schmitz r=0.03356494
	
	# Are CG sites more likely to vary as exon number increases? (spoiler: YES but weakly)
	print(ggplot(gff.exons, aes(x=as.factor(exon_no), y=variable_count/CG_sites_count))  + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Proportion of CG sites varying by exon, all genes"))
	cor(gff.exons[gff.exons$CG_sites_count>0,]$exon_no, gff.exons[gff.exons$CG_sites_count>0,]$variable_count/gff.exons[gff.exons$CG_sites_count>0,]$CG_sites_count)
	# Schmitz r=0.1495166
	
	print(ggplot(gff.introns, aes(x=as.factor(exon_no), y=variable_count/CG_sites_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Proportion of CG sites varying by intron, all genes"))
	cor(gff.introns[gff.introns$CG_sites_count>0,]$exon_no, gff.introns[gff.introns$CG_sites_count>0,]$variable_count/gff.introns[gff.introns$CG_sites_count>0,]$CG_sites_count)
	# Schmitz r=0.1147976
	
	# Do we still find these correlations if we only consider GBM genes? (even weaker, but looks like clearly related from the boxplot)
	print(ggplot(gff.exons[gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,], aes(x=as.factor(exon_no), y=variable_count/CG_sites_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Proportion of CG sites varying by exon, GBM genes"))
	cor(gff.exons[(gff.exons$CG_sites_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no, gff.exons[(gff.exons$CG_sites_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$variable_count/gff.exons[(gff.exons$CG_sites_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$CG_sites_count)
	# Schmitz r=0.1026963
	
	print(ggplot(gff.introns[gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,], aes(x=as.factor(exon_no), y=variable_count/CG_sites_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Proportion of CG sites varying by intron, GBM genes"))
	cor(gff.introns[(gff.introns$CG_sites_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no, gff.introns[(gff.introns$CG_sites_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$variable_count/gff.introns[(gff.introns$CG_sites_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$CG_sites_count)
	# Schmitz r=0.07989526
	
	# We also see a positive correlation between proportion of methylated sites and exon number
	print(ggplot(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=all_M_count/CG_sites_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Proportion of CG sites methylated by exon, GBM genes"))
	cor(gff.exons[(gff.exons$CG_sites_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no, gff.exons[(gff.exons$CG_sites_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count/gff.exons[(gff.exons$CG_sites_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$CG_sites_count)
	# Schmitz r=0.1726383
	
	print(ggplot(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=all_M_count/CG_sites_count)) +  theme_minimal() + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + ggtitle("Proportion of CG sites methylated by intron, GBM genes"))
	cor(gff.introns[(gff.introns$CG_sites_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no, gff.introns[(gff.introns$CG_sites_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count/gff.introns[(gff.introns$CG_sites_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$CG_sites_count)
	# Schmitz r=0.1481438
	
	#Does the relationship between no variable sites and no methylated sites change with exon no? (declines but not overall correlated)
	print(ggplot(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=as.factor(exon_no), y=variable_count/all_M_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Ratio of varying to methylated CG sites by exon, GBM genes"))
	cor(gff.exons[(gff.exons$all_M_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no, gff.exons[(gff.exons$all_M_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$variable_count/gff.exons[(gff.exons$all_M_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count)
	# Schmitz r=-0.1067306
	
	print(ggplot(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=as.factor(exon_no), y=variable_count/all_M_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + ggtitle("Ratio of varying to methylated CG sites by intron, GBM genes"))
	cor(gff.introns[(gff.introns$all_M_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no, gff.introns[(gff.introns$all_M_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$variable_count/gff.introns[(gff.introns$all_M_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count)
	# Schmitz r=-0.07679508
	
	# Does the direction of methylation change vary with exon number? (YES - M->U > U->M increasingly as exon No. increases, but not in introns)
	print(ggplot(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=m_to_u_count-u_to_m_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_hline(aes(yintercept = 0)) + ggtitle("Imbalance in direction of methylation change by exon, GBM genes"))
		
	print(ggplot(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=m_to_u_count-u_to_m_count)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_hline(aes(yintercept = 0)) + ggtitle("Imbalance in direction of methylation change by intron, GBM genes"))
	
	# Are the numbers of transitions in each direction in proportion with the overall state of methylation of the exon/intron? (YES for M->U, but with a slight effect of exon no for U->M)
	exons_m_trans_model=lm(gff.exons[(gff.exons$all_M_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$m_to_u_count/(gff.exons[(gff.exons$all_M_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count + gff.exons[(gff.exons$all_M_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$m_to_u_count) ~ gff.exons[(gff.exons$all_M_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	print(ggplot(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=m_to_u_count/(all_M_count+m_to_u_count))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_abline(intercept = coef(exons_m_trans_model)[1], slope = coef(exons_m_trans_model)[2], color="black", size=.1) + ggtitle("Proportion of methylated sites changing M->U by exon, GBM genes"))

	exons_u_trans_model=lm(gff.exons[(gff.exons$all_U_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$u_to_m_count/(gff.exons[(gff.exons$all_U_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_U_count + gff.exons[(gff.exons$all_U_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$u_to_m_count) ~ gff.exons[(gff.exons$all_U_count>0) & (gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	print(ggplot(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=u_to_m_count/(all_U_count+u_to_m_count))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_abline(intercept = coef(exons_u_trans_model)[1], slope = coef(exons_u_trans_model)[2], color="black", size=.1) + ggtitle("Proportion of unmethylated sites changing U->M by exon, GBM genes"))

	introns_m_trans_model=lm(gff.introns[(gff.introns$all_M_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$m_to_u_count/(gff.introns[(gff.introns$all_M_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count + gff.introns[(gff.introns$all_M_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$m_to_u_count) ~ gff.introns[(gff.introns$all_M_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	print(ggplot(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=m_to_u_count/(all_M_count+m_to_u_count))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_abline(intercept = coef(introns_m_trans_model)[1], slope = coef(introns_m_trans_model)[2], color="black", size=.1) + ggtitle("Proportion of methylated sites changing M->U by intron, GBM genes"))

	introns_u_trans_model=lm(gff.introns[(gff.introns$all_U_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$u_to_m_count/(gff.introns[(gff.introns$all_U_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_U_count + gff.introns[(gff.introns$all_U_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$u_to_m_count) ~ gff.introns[(gff.introns$all_U_count>0) & (gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	print(ggplot(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=u_to_m_count/(all_U_count+u_to_m_count))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal() + geom_abline(intercept = coef(introns_u_trans_model)[1], slope = coef(introns_u_trans_model)[2], color="black", size=.1) + ggtitle("Proportion of unmethylated sites changing U->M by intron, GBM genes"))

	dev.off()
	
	
	# What is the best model for numbers of transitions in each direction per exon/intron?
	# Overdispersed would indicate negative binomial, neutral Poisson, underdispersed Binomial
	mean(gff.exons[,]$m_to_u_count)
	# Schmitz: 0.0673076. Fishers p<0.05, Binomial p<0.005: 0.08548833. More sensitive change detection: 0.1813037
	var(gff.exons[,]$m_to_u_count)
	# Schmitz: 0.1201823. Fishers p<0.05, Binomial p<0.005: 0.1629113. More sensitive change detection: 0.4579262
	mean(gff.exons[,]$u_to_m_count)
	# Schmitz: 0.05603315. Fishers p<0.05, Binomial p<0.005: 0.0753508. More sensitive change detection: 0.1678402
	var(gff.exons[,]$u_to_m_count)
	# Schmitz: 0.08842758. Fishers p<0.05, Binomial p<0.005: 0.1285242. More sensitive change detection: 0.3814964
	mean(gff.introns[,]$m_to_u_count)
	# Schmitz: 0.07781043. Fishers p<0.05, Binomial p<0.005: 0.0952231. More sensitive change detection: 0.179759
	var(gff.introns[,]$m_to_u_count)
	# Schmitz: 0.1892647. Fishers p<0.05, Binomial p<0.005: 0.2557544. More sensitive change detection: 0.6730327
	mean(gff.introns[,]$u_to_m_count)
	# Schmitz: 0.05254412. Fishers p<0.05, Binomial p<0.005: 0.06765668. More sensitive change detection: 0.1400885
	var(gff.introns[,]$u_to_m_count)
	# Schmitz: 0.08994424. Fishers p<0.05, Binomial p<0.005: 0.1260843. More sensitive change detection: 0.3697017
	
	# Overdispersed on all fronts - Negative Binomial applies (but is it a good fit? fit the models and check the plots for good ness of fit.  Looks decent for Schmitz data)
	transition_model_class="Negative Binomial"
	goodness_of_fit_class="nbinomial"
	
	# Accumulate results from the model fitting to each group of data (no's of m->u/u->m transitions per exon/intron)
	transition_no_counts=data.frame(feature_type=character(),transition_type=character(),feature_no=integer(),mu=numeric(),size=numeric(), stringsAsFactors=FALSE)

	# We are going to use the Ricci, 2005 document to lead us in fitting a negative binomial model to the methylation status change count data, and testing for goodness of fit using Chi Squared test (https://cran.r-project.org/doc/contrib/Ricci-distributions-en.pdf)
	
	#install.packages("MASS")  # for fitting negative binomial models to data
	library(MASS)
	#install.packages("vcd")  # for goodness of fit tests
	library(vcd)
		
		
	pdf(paste0(project_id,"_",meth_context,"_Mstatus_transitions_per_exon_intron_by_exon_no.pdf"))
	cat("M->U transitions per exon, by exon No:\n")
	for (exon_no in 1:max(gff.exons$exon_no)) {
		transition_model=fitdistr(gff.exons[gff.exons$exon_no==exon_no,]$m_to_u_count, transition_model_class)
		transition_no_counts=rbind(transition_no_counts,data.frame(feature_type="exon",transition_type="M_U",feature_no=exon_no,mu=transition_model$estimate[2],size=transition_model$estimate[1]))
		cat(paste0("Exon No: ",exon_no," mu: ",coef(transition_model)[2],"\n"))
		summary(goodfit(as.integer(gff.exons[gff.exons$exon_no==exon_no,]$m_to_u_count), type=goodness_of_fit_class, method="MinChisq"))
		barplot(table(gff.exons[gff.exons$exon_no==exon_no,]$m_to_u_count)/nrow(gff.exons[gff.exons$exon_no==exon_no,]), xlab="No. of M->U transitions", ylab=paste0("Proportion of exon No. ",exon_no,"s"), main=paste0("Numbers of transitions per exon\nwith Negative Binomial model fitted to count data"))
		transition_model.density=dnbinom(0:21,mu=transition_model$estimate[2],size=transition_model$estimate[1])
		lines(transition_model.density)
		rm(transition_model, transition_model.density)
	}
	cat("U->M transitions per exon, by exon No:\n")
	for (exon_no in 1:max(gff.exons$exon_no)) {
		transition_model=fitdistr(gff.exons[gff.exons$exon_no==exon_no,]$u_to_m_count, transition_model_class)
		transition_no_counts=rbind(transition_no_counts,data.frame(feature_type="exon",transition_type="U_M",feature_no=exon_no,mu=transition_model$estimate[2],size=transition_model$estimate[1]))
		cat(paste0("Exon No: ",exon_no," mu: ",coef(transition_model)[2],"\n"))
		summary(goodfit(as.integer(gff.exons[gff.exons$exon_no==exon_no,]$u_to_m_count), type=goodness_of_fit_class, method="MinChisq"))
		barplot(table(gff.exons[gff.exons$exon_no==exon_no,]$u_to_m_count)/nrow(gff.exons[gff.exons$exon_no==exon_no,]), xlab="No. of U->M transitions", ylab=paste0("Proportion of exon No. ",exon_no,"s"), main=paste0("Numbers of transitions per exon\nwith Negative Binomial model fitted to count data"))
		transition_model.density=dnbinom(0:21,mu=transition_model$estimate[2],size=transition_model$estimate[1])
		lines(transition_model.density)
		rm(transition_model, transition_model.density)
	}
	cat("M->U transitions per intron, by exon No:\n")
	for (exon_no in 1:max(gff.exons$exon_no)) {
		transition_model=fitdistr(gff.introns[gff.introns$exon_no==exon_no,]$m_to_u_count, transition_model_class)
		transition_no_counts=rbind(transition_no_counts,data.frame(feature_type="intron",transition_type="M_U",feature_no=exon_no,mu=transition_model$estimate[2],size=transition_model$estimate[1]))
		cat(paste0("Exon No: ",exon_no," mu: ",coef(transition_model)[2],"\n"))
		summary(goodfit(as.integer(gff.introns[gff.introns$exon_no==exon_no,]$m_to_u_count), type=goodness_of_fit_class, method="MinChisq"))
		barplot(table(gff.introns[gff.introns$exon_no==exon_no,]$m_to_u_count)/nrow(gff.introns[gff.introns$exon_no==exon_no,]), xlab="No. of M->U transitions", ylab=paste0("Proportion of intron No. ",exon_no,"s"), main=paste0("Numbers of transitions per intron\nwith Negative Binomial model fitted to count data"))
		transition_model.density=dnbinom(0:21,mu=transition_model$estimate[2],size=transition_model$estimate[1])
		lines(transition_model.density)
		rm(transition_model, transition_model.density)
	}
	cat("U->M transitions per intron, by exon No:\n")
	for (exon_no in 1:max(gff.exons$exon_no)) {
		transition_model=fitdistr(gff.introns[gff.introns$exon_no==exon_no,]$u_to_m_count, transition_model_class)
		transition_no_counts=rbind(transition_no_counts,data.frame(feature_type="intron",transition_type="U_M",feature_no=exon_no,mu=transition_model$estimate[2],size=transition_model$estimate[1]))
		cat(paste0("Exon No: ",exon_no," mu: ",coef(transition_model)[2],"\n"))
		summary(goodfit(as.integer(gff.introns[gff.introns$exon_no==exon_no,]$u_to_m_count), type=goodness_of_fit_class, method="MinChisq"))
		barplot(table(gff.introns[gff.introns$exon_no==exon_no,]$u_to_m_count)/nrow(gff.introns[gff.introns$exon_no==exon_no,]), xlab="No. of U->M transitions", ylab=paste0("Proportion of intron No. ",exon_no,"s"), main=paste0("Numbers of transitions per intron\nwith Negative Binomial model fitted to count data"))
		transition_model.density=dnbinom(0:21,mu=transition_model$estimate[2],size=transition_model$estimate[1])
		lines(transition_model.density)
		rm(transition_model, transition_model.density)
	}

	#  In the majority of cases above (M-U/U-M per exon/intron number) we can't reject the hypothesis that the count data are distributed according to a negative binomial distribution model (p>0.05).  However there are a few cases where this does not appear to hold.

	# When we take the counts data as a whole, ignoring the exon number, only the U->M counts in introns cannot be distinguished from a negative binomial model.  It might be the case that there is underlying structure, related to exon number, which distorts the data from a negative binomial distribution, or that aggregating the counts together gives large enough numbers to distinguish easily from the idealised distribution:
	summary(goodfit(as.integer(gff.exons[,]$m_to_u_count), type=goodness_of_fit_class, method="MinChisq"))
	summary(goodfit(as.integer(gff.exons[,]$u_to_m_count), type=goodness_of_fit_class, method="MinChisq"))
	summary(goodfit(as.integer(gff.introns[,]$m_to_u_count), type=goodness_of_fit_class, method="MinChisq"))
	summary(goodfit(as.integer(gff.introns[,]$u_to_m_count), type=goodness_of_fit_class, method="MinChisq"))
	
	# add a factor to the accumulated model fit data for plotting
	transition_no_counts=cbind(transition_no_counts,"plotgroup"=paste0(transition_no_counts$feature_type," ",transition_no_counts$transition_type))
	print(ggplot(transition_no_counts, aes(x=feature_no, y=mu)) +geom_point(aes(colour=plotgroup))  +geom_line(aes(colour=plotgroup)) +theme_minimal() +xlab("Exon/Intron number") +ylab("mu transitions per exon/intron\nfrom fitted Negative Binomial model"))
	# pointrange version:
	print(ggplot() +geom_pointrange(data=transition_no_counts[transition_no_counts$feature_no<25,], aes(x=feature_no, y=mu, ymin=mu-size/2, ymax=mu+size/2, colour=plotgroup))   +theme_minimal() +xlab("Exon/Intron number") +ylab("mu transitions per exon/intron\nfrom fitted Negative Binomial model") + facet_wrap(~plotgroup))
	
	dev.off()
	
	
	
	# Fit linear model to relationship between proportion of sites methylated and exon number:
	### This might be better done with average methylation for the exon if we go back and calculate that later
	exons_m_prop_model=lm(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count/(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count+gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_U_count+gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$variable_count) ~ gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	
	cat(paste0("Relationship between all_M_count/(all_M_count+all_U_count+variable_count) and exon no in exons:\n"))
	coef(exons_m_prop_model)
	
	ggplot(gff.exons[(gff.exons$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=all_M_count/(all_M_count+all_U_count+variable_count))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal()  + geom_abline(intercept = coef(exons_m_prop_model)[1], slope = coef(exons_m_prop_model)[2], color="black", size=.1)
	
	introns_m_prop_model=lm(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count/(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_M_count+gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$all_U_count+gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$variable_count) ~ gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),]$exon_no)
	
	cat(paste0("Relationship between all_M_count/(all_M_count+all_U_count+variable_count) and exon no in introns:\n"))
	coef(introns_m_prop_model)
	
	ggplot(gff.introns[(gff.introns$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID),], aes(x=exon_no, y=all_M_count/(all_M_count+all_U_count+variable_count))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") + theme_minimal()  + geom_abline(intercept = coef(introns_m_prop_model)[1], slope = coef(introns_m_prop_model)[2], color="black", size=.1)
	

	
	# Summarise the results of the rep concordance analysis
	library(plyr)
	line_gen_no = 0
	rep_concordance_summary_genes=matrix(, nrow=ncol(all_reps_meth_status)-3, ncol=6)
	rownames(rep_concordance_summary_genes)=colnames(all_reps_meth_status[,4:ncol(all_reps_meth_status)])
	for(line_gen in rownames(rep_concordance_summary_genes)) {
		line_gen_no = line_gen_no + 1
		cat(paste0("Summarising methylation call basis for ",line_gen,"\n"))
		status_counts = count(all_reps_meth_status[queryHits(olaps_CG_sites_reduced) ,line_gen])
		# Needs a hack for when "P" is zero (nrow should be 6 otherwise):
		if(nrow(status_counts)==5) {
			status_counts$x = as.character(status_counts$x)
			status_counts[6,]=status_counts[5,]
			status_counts[5,]=c("P",0)
			status_counts$x = as.factor(status_counts$x)
			status_counts$freq=as.numeric(status_counts$freq)
		}
		rep_concordance_summary_genes[line_gen_no,]=status_counts$freq
		colnames(rep_concordance_summary_genes)=status_counts$x
	}
	# output the summary of numbers of concordant and discordant sites within gene space
	rep_concordance_summary_genes
	# output the proportion of sites in each line with a 'clean' M/U call within gene space
	cat(paste0("Proportion of ",meth_context," sites in each line with rep-concordant M/U call within non-heterochromatic genes: \n"))
	(rep_concordance_summary_genes[,3]+rep_concordance_summary_genes[,6])/sum(rep_concordance_summary_genes[,1:6])*nrow(rep_concordance_summary_genes)

	line_gen_no = 0
	rep_concordance_summary_transposons=matrix(, nrow=ncol(all_reps_meth_status)-3, ncol=6)
	rownames(rep_concordance_summary_transposons)=colnames(all_reps_meth_status[,4:ncol(all_reps_meth_status)])
	for(line_gen in rownames(rep_concordance_summary_transposons)) {
		line_gen_no = line_gen_no + 1
		cat(paste0("Summarising methylation call basis for ",line_gen,"\n"))
		status_counts = count(all_reps_meth_status[queryHits(tolaps_CG_sites_reduced) ,line_gen])
		# Needs a hack for when "P" is zero (nrow should be 6 otherwise):
		if(nrow(status_counts)==5) {
			status_counts$x = as.character(status_counts$x)
			status_counts[6,]=status_counts[5,]
			status_counts[5,]=c("P",0)
			status_counts$x = as.factor(status_counts$x)
			status_counts$freq=as.numeric(status_counts$freq)
		}
		rep_concordance_summary_transposons[line_gen_no,]=status_counts$freq
		colnames(rep_concordance_summary_transposons)=status_counts$x
	}
	# output the summary of numbers of concordant and discordant sites within transposon space
	rep_concordance_summary_transposons
	# output the proportion of sites in each line with a 'clean' M/U call within transposon space
	cat(paste0("Proportion of ",meth_context," sites in each line with rep-concordant M/U call within transposons: \n"))
	(rep_concordance_summary_transposons[,3]+rep_concordance_summary_transposons[,6])/sum(rep_concordance_summary_transposons[,1:6])*nrow(rep_concordance_summary_transposons)

### Control point 1
	save.image(file=paste0(project_id,"_",meth_context,"_control_point_1.RData"))
	# load (paste0(project_id,"_",meth_context,"_control_point_1.RData"))

	
	# Load the BSGenome version of the reference 
	library(BSgenome)

	# List the available genome assemblies
	available.genomes()

	# TAIR9 and TAIR10 correspond to the same genome assembly so there is no need for a BSgenome pkg for TAIR10 :-)
	# https://stat.ethz.ch/pipermail/bioconductor/2010-December/036938.html
	library(BSgenome.Athaliana.TAIR.TAIR9)

	# Find the sequence lengths for each chromosome
	sLengths=seqlengths(Athaliana)


	# Use the MethylSeekR package to segment the methylome
	library(MethylSeekR)

	# MethylSeekR expects an input file with 4 columns:  Chrom, Locus, T and M.  T is total reads, M is methylated reads.
	
	# Create an input file for MethylSeekR based on the merged methylomes of all valid samples
	
	# Sum the M (Cov_C) and T (Cov_C+Cov_T) counts for all valid samples 

	# Create empty tables to accumulate each of the 'across samples' columns data
	across_samples_TM=matrix(0, nrow=nrow(all_samples_meth_status), ncol=2)

	for(sample_name in valid_samples) {
		if(opt$verbose) {cat(paste0("Adding methylation read counts to across-samples per-site table: ",sample_name,"\n"))}
		# Merge each of the sample counts into the relevant accumulation table, accounting for possible NA values

		across_samples_TM[,1]=ifelse(is.na(coverage_data[coverage_data$Sample==sample_name,]$Cov_C+coverage_data[coverage_data$Sample==sample_name,]$Cov_T),across_samples_TM[,1],across_samples_TM[,1]+coverage_data[coverage_data$Sample==sample_name,]$Cov_C+coverage_data[coverage_data$Sample==sample_name,]$Cov_T)

		across_samples_TM[,2]=ifelse(is.na(coverage_data[coverage_data$Sample==sample_name,]$Cov_C),across_samples_TM[,2],across_samples_TM[,2]+coverage_data[coverage_data$Sample==sample_name,]$Cov_C)

	} # end for each sample

	# Write out a copy of the methylome in the format that MethylSeekR likes to read in
	write.table(cbind(paste0(substr(as.character(all_samples_meth_status$Chromosome),1,1),tolower(substr(as.character(all_samples_meth_status$Chromosome),2,3)),substr(as.character(all_samples_meth_status$Chromosome),4,4)),all_reps_meth_status$Locus,across_samples_TM),file=paste0(project_id,"_",meth_context,"_TM_read_counts_across_samples.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	
	# Write out a copy of the methylome in Bismark .cov format for visualising in SeqMonk
	write.table(cbind(as.character(all_reps_meth_status$Chromosome),all_reps_meth_status$Locus,ifelse(meth_context=="CG",1,0)+all_reps_meth_status$Locus,across_samples_TM[,2]/across_samples_TM[,1],across_samples_TM[,2],across_samples_TM[,1]-across_samples_TM[,2]),file=paste0(project_id,"_",meth_context,"_CT_read_counts_across_samples.cov"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	
	
	# Read the methylome data in as a Granges object
	meth.gr <- readMethylome(FileName=paste0(project_id,"_",meth_context,"_TM_read_counts_across_samples.tsv"), seqLengths=sLengths)
	
	### Read in the SNPs data
	#snps.gr <- readSNPTable(FileName=snpFname, seqLengths=sLengths)

	### Remove the CpG sites from the SNP loci
	#meth.gr <- removeSNPs(meth.gr, snps.gr)
	
	plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr1", num.cores=1)
	plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr2", num.cores=1)
	plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr3", num.cores=1)
	plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr4", num.cores=1)
	plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr5", num.cores=1)
	plotAlphaDistributionOneChr(m=meth.gr, chr.sel="ChrM", num.cores=1)
	plotAlphaDistributionOneChr(m=meth.gr, chr.sel="ChrC", num.cores=1)

	# Posterior means of alpha<1 indicate a bias towards very low or very high methylation.  Onlt ChrM has a bimodal distribution, with significant density above 1.0, indicative of likely presence of partially methylated domains (PMDs) (although in our case it is more likely due to multiple copies of differently methylated mitochondrial genomes being analysed together).  Accoringly we skip the next part.

#	library(parallel)
#	detectCores()
	
#	PMDsegments.gr1 <- segmentPMDs(m=meth.gr, chr.sel="Chr1", seqLengths=sLengths, num.cores=1)
#	PMDsegments.gr2 <- segmentPMDs(m=meth.gr, chr.sel="Chr2", seqLengths=sLengths, num.cores=1)
#	PMDsegments.gr3 <- segmentPMDs(m=meth.gr, chr.sel="Chr3", seqLengths=sLengths, num.cores=1)
#	PMDsegments.gr4 <- segmentPMDs(m=meth.gr, chr.sel="Chr4", seqLengths=sLengths, num.cores=1)
#	PMDsegments.gr5 <- segmentPMDs(m=meth.gr, chr.sel="Chr5", seqLengths=sLengths, num.cores=1)
#	PMDsegments.grM <- segmentPMDs(m=meth.gr, chr.sel="ChrM", seqLengths=sLengths, num.cores=1)
#	PMDsegments.grC <- segmentPMDs(m=meth.gr, chr.sel="ChrC", seqLengths=sLengths, num.cores=1)
	
	# Confirm that allowing for PMDs has removed bimodality from distribution of posterior of alpha:
#	plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.gr, PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), chr.sel="Chr1", num.cores=1)
#	plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.gr, PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), chr.sel="Chr2", num.cores=1)
#	plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.gr, PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), chr.sel="Chr3", num.cores=1)
#	plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.gr, PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), chr.sel="Chr4", num.cores=1)
#	plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.gr, PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), chr.sel="Chr5", num.cores=1)
	#plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.gr, PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), chr.sel="ChrM", num.cores=1)
	#plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.gr, PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), chr.sel="ChrC", num.cores=1)
	
#	plotPMDSegmentation(m=meth.gr, segs=PMDsegments.gr1, pdfFilename="PMDs1.pdf")
#	plotPMDSegmentation(m=meth.gr, segs=PMDsegments.gr2, pdfFilename="PMDs2.pdf")
#	plotPMDSegmentation(m=meth.gr, segs=PMDsegments.gr3, pdfFilename="PMDs3.pdf")
#	plotPMDSegmentation(m=meth.gr, segs=PMDsegments.gr4, pdfFilename="PMDs4.pdf")
#	plotPMDSegmentation(m=meth.gr, segs=PMDsegments.gr5, pdfFilename="PMDs5.pdf")
#	plotPMDSegmentation(m=meth.gr, segs=PMDsegments.grM, pdfFilename="PMDsM.pdf")
	#plotPMDSegmentation(m=meth.gr, segs=PMDsegments.grC, pdfFilename="PMDsC.pdf")
	
#	savePMDSegments(PMDs=PMDsegments.gr1, GRangesFilename="PMDs.gr1.rds", TableFilename="PMDs1.tab")
#	savePMDSegments(PMDs=PMDsegments.gr2, GRangesFilename="PMDs.gr2.rds", TableFilename="PMDs2.tab")
#	savePMDSegments(PMDs=PMDsegments.gr3, GRangesFilename="PMDs.gr3.rds", TableFilename="PMDs3.tab")
#	savePMDSegments(PMDs=PMDsegments.gr4, GRangesFilename="PMDs.gr4.rds", TableFilename="PMDs4.tab")
#	savePMDSegments(PMDs=PMDsegments.gr5, GRangesFilename="PMDs.gr5.rds", TableFilename="PMDs5.tab")
#	savePMDSegments(PMDs=PMDsegments.grM, GRangesFilename="PMDs.grM.rds", TableFilename="PMDsM.tab")
#	#savePMDSegments(PMDs=PMDsegments.grC, GRangesFilename="PMDs.grC.rds", TableFilename="PMDsC.tab")
	
	# The PMDs identified by the model fitting do not look at all convincing in the plots of sample genomic regions.
	# Putative PMDs are masked out before the clustering of the genome into low and highly methylated regions
	# The next steps with this package would be to estimate FDR by masking out CpG islands to identify background levels of methylation outside CpG islands.  This may or may not be appropriate in plants


	# The next section is supposed to retrieve CPG islands from UCSC genome browser, but no point as they don't have Ath, and Ath doen's have CPG islands
	#library(rtracklayer)
	#session <- browserSession()
	#genome(session) <- "hg18"
	#query <- ucscTableQuery(session, "cpgIslandExt")
	#CpGislands.gr <- track(query)
	#genome(CpGislands.gr) <- NA

	
	# Plants don't have CPG islands, but we will use the mask_loci data frame in place of CpG_islands to construct a set of genomic regions to be masked (ignored).
	# This will be the chloroplast and mitochondrial molecules, and genes identified as 'Unknown' (including the mitochondrial genes from chromosome 2)
	mask_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	mask_loci=rbind(mask_loci,list(Chromosome="ChrM",start_site=1,end_site=sLengths[6]),list(Chromosome="ChrC",start_site=1,end_site=sLengths[7]))
	#colnames(mask_loci)=c("Chromosome","start_site","end_site")
	
	# Leave unknown genes out of the mask now - they seem to increase numbers of segments unnecessarily
	#unknown_gene_loci=gff.genes[gff.genes$m_class=="Unknown",c("V1","V4","V5")]
	#colnames(unknown_gene_loci)=colnames(mask_loci)
	#unknown_gene_loci$Chromosome=paste0(substr(unknown_gene_loci$Chromosome,1,1),tolower(substr(unknown_gene_loci$Chromosome,2,3)),substr(unknown_gene_loci$Chromosome,4,4))
	#mask_loci=rbind.data.frame(mask_loci,unknown_gene_loci)
	#rm(unknown_gene_loci)

	mask_loci.gr=makeGRangesFromDataFrame(df=mask_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome")

	# Define genomicranges representing the whole genome
	genome_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	levels(genome_loci$Chromosome)=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")
	for (chrom_no in 1:length(sLengths)) {
		genome_loci[chrom_no,]=list(Chromosome=toupper(names(sLengths[chrom_no])),start_site=1,end_site=sLengths[chrom_no])
	}
	genome_loci.gr=makeGRangesFromDataFrame(df=genome_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome")


	# Calculate FDRs
	# This code is taken from the package MethylSeekR
	# It doesn't seem to work properly here, because all FDRs end up being several 100%
	calculateFDRs <- function (m, CGIs, PMDs = NA, pdfFilename = NULL, num.cores = 1, 
    nCpG.smoothing = 3, meth.cutoffs = seq(0.3, 0.7, by = 0.1), 
    nCpG.cutoffs = seq(1, 6, by = 1), minCover = 5) 
{
    m = m[values(m)[, 1] >= minCover]
    message("calculating false discovery rate")
    n <- length(meth.cutoffs) * length(nCpG.cutoffs)
    parameters <- vector("list", n)
    count = 1
    for (meth.cutoff in meth.cutoffs) {
        for (k in nCpG.cutoffs) {
            parameters[[count]] <- c(meth.cutoff, k)
            count <- count + 1
        }
    }
    meth.notrand <- m
    if (class(PMDs) == "GRanges") {
        message("removing PMDs for randomization")
        meth.notrand <- subsetByOverlaps(m, PMDs[values(PMDs)$type == 
            "notPMD"])
    }
    meth.rand <- meth.notrand
    ov <- findOverlaps(meth.rand, CGIs)
    meth.rand <- meth.rand[-unique(queryHits(ov))]
    values(meth.rand) <- values(meth.rand)[sample(length(meth.rand)), 
        ]
    res <- mclapply(parameters, function(params) {
        meth.cutoff <- params[1]
        k <- params[2]
        mean.meth <- runmean(Rle(values(meth.notrand)[, 2]/values(meth.notrand)[, 
            1]), k = nCpG.smoothing, endrule = "constant")
        indx <- mean.meth < meth.cutoff
        nSeg <- sum(runValue(indx) == TRUE & runLength(indx) >= 
            k)
        mean.meth = runmean(Rle(values(meth.rand)[, 2]/values(meth.rand)[, 
            1]), k = nCpG.smoothing, endrule = "constant")
        indx <- mean.meth < meth.cutoff
        nSeg.rand <- sum(runValue(indx) == TRUE & runLength(indx) >= 
            k)
        c(nSeg, nSeg.rand)
    }, mc.cores = num.cores)
    FDRs = 100 * sapply(res, function(x) {
        x[2]/x[1]
    })
    tmp = matrix(NA, nrow = length(meth.cutoffs), ncol = length(nCpG.cutoffs))
    rownames(tmp) = meth.cutoffs
    colnames(tmp) = nCpG.cutoffs
    count = 1
    for (meth.cutoff in meth.cutoffs) {
        for (k in nCpG.cutoffs) {
            tmp[as.character(meth.cutoff), as.character(k)] = FDRs[count]
            count = count + 1
        }
    }
    FDRs = tmp
    numSegments = sapply(res, function(x) {
        x[1]
    })
    tmp = matrix(NA, nrow = length(meth.cutoffs), ncol = length(nCpG.cutoffs))
    rownames(tmp) = meth.cutoffs
    colnames(tmp) = nCpG.cutoffs
    count = 1
    for (meth.cutoff in meth.cutoffs) {
        for (k in nCpG.cutoffs) {
            tmp[as.character(meth.cutoff), as.character(k)] = numSegments[count]
            count = count + 1
        }
    }
    numSegments = tmp
    rownames(FDRs) = as.character(100 * as.numeric(rownames(FDRs)))
    rownames(numSegments) = as.character(100 * as.numeric(rownames(numSegments)))
    if (!is.null(pdfFilename)) {
        pdf(pdfFilename, width = 9, height = 4.5)
    }
    par(mfrow = c(1, 2))
    barplot(pmin((t(FDRs)), 80), beside = TRUE, ylab = "FDR (%)", 
        ylim = c(0, 80), xlab = "methylation cut-off (%)")
    barplot(t(numSegments), beside = TRUE, ylab = "number of segments", 
        xlab = "methylation cut-off (%)")
    legend("topleft", legend = paste(colnames(FDRs), c("CpG", 
        rep("CpGs", ncol(FDRs) - 1)), sep = " "), fill = grey.colors(ncol(FDRs)), 
        bty = "n")
    if (!is.null(pdfFilename)) 
        dev.off()
    rownames(FDRs) = as.character(as.numeric(rownames(FDRs))/100)
    rownames(numSegments) = as.character(as.numeric(rownames(numSegments))/100)
    list(FDRs = FDRs, numSegments = numSegments)
}
	stats <- calculateFDRs(m=meth.gr, CGIs=mask_loci.gr, num.cores=1, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_FDRs.pdf"))
	stats

	#$FDRs
	#		1        2        3        4        5        6
	#0.3 711.3040 627.3805 551.1079 467.4624 385.3993 314.4064
	#0.4 570.7304 514.6336 566.9219 568.0249 539.4143 510.4416
	#0.5 544.9517 502.1575 554.8679 566.1457 547.7668 525.5805
	#0.6 396.7911 386.9440 420.8281 457.6987 468.9901 471.3245
	#0.7 139.0559 143.5250 147.0059 178.0157 194.3320 205.6945

	#$numSegments
	#		1     2     3     4     5     6
	#0.3 44179 37797 33034 29741 27526 25801
	#0.4 46304 43113 36698 31781 28888 26787
	#0.5 45738 42085 35970 31033 28076 25969
	#0.6 46839 42356 37507 31843 28507 26228
	#0.7 44188 42718 41682 33351 29852 27518	
		

	# For human, from the vignette https://bioconductor.org/packages/release/bioc/vignettes/MethylSeekR/inst/doc/MethylSeekR.pdf, , m.sel=0.5 and n.sel=4 are chosen from the FDR calculation.
	# For Arabidopsis, FDR figures are all >100 (supposed to be percentages), so it is not clear how to estimate m.sel and n.sel to minimise FDR.
	# For ranges of values tried (n=1,2,3,4,5,6, m=0.5), numbers of segements are:

	#m=0.5, n=1  45739
	#m=0.5, n=2  39508
	#m=0.5, n=3  28878
	#m=0.5, n=4  22703
	#m=0.5, n=5  19583
	#m=0.5, n=6  17677
	#m=0.5, n=7  16142

	# There is not such a clear separation between unmethylated and low-methylated regions as in vertebrates. In vertebrates, there are shorter, low-methylated segments, and longer unmethylated segments.  In Arabidopsis, the length range of unmethylated segments overlaps that of low-methylated segments, whichever m and n values are chosen for the segmentation.
		
	# m.sel seems to be the maximum level of methylation which is entertained to be either unmethylated or low-methylated, n.sel seems to represent possibly the minimum length of segments in CG sites. Assuming this, I could choose 0.75 as the max level (this is the highest percentage of methylation found in a non-heterochromatic gene model), or alternatively choose 0.6 as per the Zilberman grou paper, where they used 0.6 methylation as upper limit of genes to be considered not heterochromatic.

	# Evaluate parameter space by looping over values of m and n.  For each model, calculate overlap of positive controls (gene loci classed as unmethylated or GBM) with segments classed as unmethylated or lowly methylated.  Compare this with overlap of negative controls (transposons and genes classified as heterochromatic). Plot ROC curves, and estimate AUC.
	# Create data frames for positive and negative 'controls' to evaluate performance of segmentation

	methylated_loci.gr=reduce(makeGRangesFromDataFrame(df=rbind(gff.transposons[,c("V1","V4","V5")],gff.genes[gff.genes$m_class=="Heterochromatic",c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1"))
	gene_body_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))
	unmethylated_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))

	# First we make the comparison between Unmethylated + GBM genes and Transposons + Heterochromatic genes
	
	# Create data frame to store the results of the optimisation:
	# For each value of m and n parameters, how many variable sites are overlapping with UMR and LMR segments? What proportion of supposed Transposons and Heterochromatic gene space overlaps with UMR or LMR segemtns (fp)?  What proportion of supposed Unmenthylated or GBM methylated gene space overlaps with UMR or LMR segments (tp)? 
	m_n_optimisation=data.frame(m=numeric(), n=numeric(), UMR=integer(), LMR=integer(), fp=numeric(), tp=numeric())

	# Loop through values for m and n parameters 
	m.seq=seq(0.05,0.95, by=0.05)
	n.seq=seq(1,6, by=1)
	for (m.sel in m.seq) {
		for (n.sel in n.seq) {
			#n.sel=4
			UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)

			# Capitalise chromosome names in segmentation object
			levels(UMRLMRsegments.gr@seqnames)=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")
			UMRLMRsegments.gr@seqinfo@seqnames=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")

			# Find overlaps between segments and positive and negative 'control' loci (annotated genes and transposons)
			fp_hits=findOverlaps(methylated_loci.gr, UMRLMRsegments.gr)
			fp_overlaps <- pintersect(methylated_loci.gr[queryHits(fp_hits)], UMRLMRsegments.gr[subjectHits(fp_hits)])
			fp_overlap_prop <- sum(width(fp_overlaps)) / sum(width(methylated_loci.gr))
			tp_hits=findOverlaps(c(gene_body_loci.gr,unmethylated_loci.gr), UMRLMRsegments.gr)
			tp_overlaps <- pintersect(c(gene_body_loci.gr,unmethylated_loci.gr)[queryHits(tp_hits)], UMRLMRsegments.gr[subjectHits(tp_hits)])
			tp_overlap_prop <- sum(width(tp_overlaps)) / sum(width(c(gene_body_loci.gr,unmethylated_loci.gr)))
			
			# Find overlaps between segments and variable loci
			variable_sites_UMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="UMR"])
			variable_sites_LMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="LMR"])
			m_n_optimisation=rbind(m_n_optimisation, list(m=m.sel, n=n.sel, UMR=length(variable_sites_UMR_olaps@from), LMR=length(variable_sites_LMR_olaps@from), fp=fp_overlap_prop, tp=tp_overlap_prop))
			# How many variable sites were captured?
			cat(paste0("m.sel=",m.sel," n.sel=",n.sel," Variable sites in UMRs: ",length(variable_sites_UMR_olaps@from),"  LMRs: ",length(variable_sites_LMR_olaps@from),"\n"))
		}
	}
	m_n_optimisation$variable_sites_captured=m_n_optimisation$UMR+m_n_optimisation$LMR

	# Plot ROC curves for each value of n, and estimate AUC
	pdf(paste0(project_id,"_",meth_context,"_MethylSeekR_m_n_optimisation_ROC.pdf"))
	print(ggplot(m_n_optimisation, aes(x=fp, y=tp)) + geom_point() +geom_text(aes(label=paste0(m)),hjust=0, vjust=0) +geom_line() +facet_wrap(~ n))
	dev.off()
	
	# Estimate AUC from trapezium approximation
	cat(paste0("Estimated AUC values from ROC curves generated by varying m parameter for each value of n in MethylSeekR:\n"))
	for (n.sel in n.seq) {
		m_n_optimisation_auc=0
		prev_tp=0
		prev_fp=0
		for (m.sel in m.seq) {
			m_n_optimisation_auc = m_n_optimisation_auc + (1-m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$fp)*(m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$tp-prev_tp) + ((m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$tp-prev_tp)*(m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$fp-prev_fp))/2
			prev_tp=m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$tp
			prev_fp=m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$fp
		}
		cat(paste0(" n=",n.sel," AUC=",m_n_optimisation_auc,"\n"))
	}
	
	# Schmitz data:
	#n=1 AUC=0.864135028898874
	#n=2 AUC=0.865661778246867
	#n=3 AUC=0.867474307742376
	#n=4 AUC=0.862142542459426
	#n=5 AUC=0.858323878986377
	#n=6 AUC=0.854695249151202

	# n=3 gives largest AUC for Schmitz data.  m=0.65:0.75 is the inflexion point of the ROC curve.

	# The following function is from the package MethylSeekR
	# The following changes have been made:
	# nCG.classification threshold on no CG sites to differentiate between LMR and UMR segments replaced by an array of 'median methylation' cutoff levels (may want to change this to use pmeth instead)
	# This was done because in Arabidopsis, there is not a clear separation in segment size between unmethylated and low-methylated segments, but there is a clear separation in methylation level, in contrast to mammals.
	# In Arabidopsis, CG sites tend to be fully methylated, or unmethylated, so using the median values introduces sampling artefacts
	# median.meth replaced by pmeth in scatterplot
	# nCG replaced by nSites for generality across contexts
	segmentUMRsLMRs <- function (m, meth.cutoff = 0.5, nCpG.cutoff = 3, PMDs = NA, pdfFilename = NULL, 
    num.cores = 1, myGenomeSeq, seqLengths, nCpG.smoothing = 3, 
    minCover = 5) 
	{
		#nCG.classification <- 30
		#mMeth.classification=c(0.2,0.4,0.75)
		mMeth.classification=c(0.2)
		message("identifying UMRs and LMRs")
		m = m[values(m)[, 1] >= minCover]
		nCGsPerChr = table(as.character(seqnames(m)))
		chrs = names(nCGsPerChr)[nCGsPerChr >= nCpG.smoothing]
		res <- mclapply(chrs, function(chr) {
			sel <- which(as.character(seqnames(m)) == chr)
			mean.meth <- runmean(Rle(values(m)[sel, 2]/values(m)[sel,1]), k = nCpG.smoothing, endrule = "constant")
			indx <- mean.meth < meth.cutoff
			runValue(indx)[runLength(indx) < nCpG.cutoff & runValue(indx) == TRUE] = FALSE
			runValue(indx)[runLength(indx) < nCpG.cutoff & runValue(indx) == FALSE] = TRUE
			tmp.ids <- rep(1:length(runLength(indx)), runLength(indx))
			tmp <- split(1:length(sel), tmp.ids)
			tmp <- tmp[runValue(indx) == TRUE]
			if (length(tmp) > 0) {
				coords <- cbind(sapply(tmp, min), sapply(tmp, max))
				starts <- round((start(m)[sel[pmax(1, coords[, 1] - 1)]] + start(m)[sel[coords[, 1]]])/2)
				ends <- round((start(m)[sel[coords[, 2]]] + start(m)[sel[pmin(length(sel), coords[, 2] + 1)]])/2)
				hmr.gr = GRanges(seqnames = unique(seqnames(m[sel])), 
					strand = "*", ranges = IRanges(starts, ends), 
					seqlengths = seqLengths)
			}
			else {
				hmr.gr = GRanges(, seqlengths = seqLengths)
			}
			hmr.gr
		}, mc.cores = num.cores)
		segments.gr = do.call(c, unname(res))
		if (class(PMDs) == "GRanges") {
			segments.gr = subsetByOverlaps(segments.gr, PMDs[values(PMDs)$type == "notPMD"])
		}
		nSites = vcountPattern("CG", getSeq(myGenomeSeq, resize(segments.gr, width(segments.gr), fix = "start"), as.character = FALSE))
		ov <- findOverlaps(m, segments.gr)
		T = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), sum)
		M = tapply(values(m[queryHits(ov)])[, 2], subjectHits(ov), sum)
		nSites.segmentation = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), length)
		median.meth = tapply(as.vector(runmean(Rle(values(m[queryHits(ov)])[,2]/values(m[queryHits(ov)])[, 1]), nCpG.smoothing, endrule = "constant")), subjectHits(ov), median)
		median.meth = pmax(0, median.meth)
		if (!all.equal(as.numeric(names(T)), 1:length(segments.gr))) {
			message("error in calculating methylation levels for PMDs")
		}
		type=ifelse(median.meth<mMeth.classification[1],"UMR", "LMR")
		#type=ifelse(median.meth<mMeth.classification[1],0,ifelse(median.meth<mMeth.classification[2],1,ifelse(median.meth<mMeth.classification[3],2,3)))
		#type = c("UMR", "LMR")[as.integer(nCG < nCG.classification) + 1]
		values(segments.gr) = DataFrame(nSites.segmentation, nSites, T, M, pmeth = M/T, median.meth = median.meth, type)
		jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
		if (!is.null(pdfFilename)) {
			pdf(pdfFilename, height = 5, width = 5)
		}
		smoothScatter(log2(values(segments.gr)$nSites), 100 * values(segments.gr)$pmeth, 
			colramp = jet.colors, xlab = "log2 number of sites in segment", 
			ylab = "mean methylation (%)")
		#abline(v = log2(nCG.classification), lty = 5)
		abline(h = 100 * mMeth.classification[1], lty = 5)
		abline(h = 100 * mMeth.classification[2], lty = 5)
		abline(h = 100 * mMeth.classification[3], lty = 5)
		if (!is.null(pdfFilename)) 
			dev.off()
		segments.gr
	}
	
	
	
	# Execute the preferred segmentation model, visualise and save the results
	#UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=0.7, nCpG.cutoff=3, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape.pdf"))

	# Although 0.7,3 works well to differentiate between methylated and GBM annotations, we try a different approach.
	# Using a cutoff of 0.4 we try to differentiate between unmethylated and weakly methylated GBM regions on the one hand, and methylated regions on the other.
	# 0.2 cutoff generated unmethylated regions which were always one CG site short at each end, whereas 0.4 threshold does not.
	# n=1 allows for single CG loci which are unmethylated or methylated among a larger contrasting region.
	
	m.sel=0.4
	n.sel=1
	UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape_",m.sel,"_",n.sel,".pdf"))
	head(UMRLMRsegments.gr)
	plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename="UMRsLMRs.gr.rds", TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_",m.sel,"_",n.sel,".tsv"))

	# Now we have defined a segmentation of the genome into UMR, LMR and FMR, let's overlap each of these sets of segments with each of the principal classifications of gene/transposon annotations and estimate the proportions of each class of annotation represented by each class of overlapped mCG segment.
	
	# Capitalise chromosome names in segmentation object
	levels(UMRLMRsegments.gr@seqnames)=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")
	UMRLMRsegments.gr@seqinfo@seqnames=levels(UMRLMRsegments.gr@seqnames)

	# Find the fully methylated regions by subtracting the UMR/LMR from the whole genome
	FMRsegments.gr=setdiff(genome_loci.gr,UMRLMRsegments.gr)
	
	# Find the overlaps of each set of genomic segments with each set of annotated loci
	M_L_hits=findOverlaps(methylated_loci.gr, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="LMR"])
	M_L_overlaps <- pintersect(methylated_loci.gr[queryHits(M_L_hits)], UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="LMR"][subjectHits(M_L_hits)])
	M_L_overlap_prop <- sum(width(M_L_overlaps)) / sum(width(methylated_loci.gr))
	M_M_hits=findOverlaps(methylated_loci.gr, FMRsegments.gr)
	M_M_overlaps <- pintersect(methylated_loci.gr[queryHits(M_M_hits)], FMRsegments.gr[subjectHits(M_M_hits)])
	M_M_overlap_prop <- sum(width(M_M_overlaps)) / sum(width(methylated_loci.gr))
	M_U_hits=findOverlaps(methylated_loci.gr, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="UMR"])
	M_U_overlaps <- pintersect(methylated_loci.gr[queryHits(M_U_hits)], UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="UMR"][subjectHits(M_U_hits)])
	M_U_overlap_prop <- sum(width(M_U_overlaps)) / sum(width(methylated_loci.gr))
	
	L_L_hits=findOverlaps(gene_body_loci.gr, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="LMR"])
	L_L_overlaps <- pintersect(gene_body_loci.gr[queryHits(L_L_hits)], UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="LMR"][subjectHits(L_L_hits)])
	L_L_overlap_prop <- sum(width(L_L_overlaps)) / sum(width(gene_body_loci.gr))
	L_M_hits=findOverlaps(gene_body_loci.gr, FMRsegments.gr)
	L_M_overlaps <- pintersect(gene_body_loci.gr[queryHits(L_M_hits)], FMRsegments.gr[subjectHits(L_M_hits)])
	L_M_overlap_prop <- sum(width(L_M_overlaps)) / sum(width(gene_body_loci.gr))
	L_U_hits=findOverlaps(gene_body_loci.gr, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="UMR"])
	L_U_overlaps <- pintersect(gene_body_loci.gr[queryHits(L_U_hits)], UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="UMR"][subjectHits(L_U_hits)])
	L_U_overlap_prop <- sum(width(L_U_overlaps)) / sum(width(gene_body_loci.gr))
	
	U_L_hits=findOverlaps(unmethylated_loci.gr, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="LMR"])
	U_L_overlaps <- pintersect(unmethylated_loci.gr[queryHits(U_L_hits)], UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="LMR"][subjectHits(U_L_hits)])
	U_L_overlap_prop <- sum(width(U_L_overlaps)) / sum(width(unmethylated_loci.gr))
	U_M_hits=findOverlaps(unmethylated_loci.gr, FMRsegments.gr)
	U_M_overlaps <- pintersect(unmethylated_loci.gr[queryHits(U_M_hits)], FMRsegments.gr[subjectHits(U_M_hits)])
	U_M_overlap_prop <- sum(width(U_M_overlaps)) / sum(width(unmethylated_loci.gr))
	U_U_hits=findOverlaps(unmethylated_loci.gr, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="UMR"])
	U_U_overlaps <- pintersect(unmethylated_loci.gr[queryHits(U_U_hits)], UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="UMR"][subjectHits(U_U_hits)])
	U_U_overlap_prop <- sum(width(U_U_overlaps)) / sum(width(unmethylated_loci.gr))
	
	# Schmitz data:
	M_L_overlap_prop
	#[1] 0.002267616
	M_M_overlap_prop
	#[1] 0.8008056
	M_U_overlap_prop
	#[1] 0.1969268

	# 80% of transposon/heterochromatic gene space is covered by  fully methylated, 20% by unmethylated and 0.2% by lowly methylated segments

	L_L_overlap_prop
	#[1] 0.02057531
	L_M_overlap_prop
	#[1] 0.248998
	L_U_overlap_prop
	#[1] 0.7304267

	# 73% of GBM gene space is covered by unmethylated segments, 24% by fully methylated and 2% by lowly methylated segments
	
	U_L_overlap_prop
	#[1] 7.258314e-05
	U_M_overlap_prop
	#[1] 0.0004609684
	U_U_overlap_prop
	#[1] 0.9994664
	
	# 99.9% of unmethylated gene space is covered by unmethylated segments

		
	# We are interested to exclude potential GBM genes as classified by mCG, where they contain a segment which is mCHG:
	# L_M_overlaps is the intersection of Lowly methylated genes (i.e. genes annotated as GBM) and methylated segments (mCG>0.4)
	# We should also check out L_L_overlaps for completeness sake (intersection of genes annotated as GBM and segments with 0.2<mCG<0.4)
	
	# Define a function to annotate an arbitrary genomicRanges object with methylation values from an arbitrary methylome
	# Function adapted from MethylSeekR package
	meth_by_segment <- function (m, segment_model, meth_context, num.cores = 1, myGenomeSeq, seqLengths, nSite.smoothing = 3, mMeth.classification=c(0.7), mMeth.classes=c("LMR","FMR")) 
	{
		# Replaced nCG with nSites for generality
		
		# m is a genomic ranges object representing a methylome
		# segment_model is an arbitrary genomic ranges object representing a set of query segments of interest

		# Count CG sites in each segment
		# Added fixed=FALSE so that ambiguity codes will work properly
		nSites = vcountPattern(meth_context, getSeq(Athaliana, resize(segment_model, width(segment_model), fix = "start"), as.character = FALSE), fixed=FALSE)
	
		# Overlap the methylome with the segments
		ov <- findOverlaps(m, segment_model)

		# Count Ts and Ms per segment
		T = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), sum)
		M = tapply(values(m[queryHits(ov)])[, 2], subjectHits(ov), sum)
		nSites.segmentation = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), length)
		# The next step carries out a robust smoothing of the methylation values across the segment.  It takes k adjacent sites (default=3) and calculates M/T for a sliding window of that size, then takes the median value of this proportion across the segment. k=3 give 'minimal' robust smoothing eliminating isolated outliers.
		# This might be problematic for us, as in cases where each site is either fully methylated or unmethylated, it will tend to generate sequences like this for a window size of 3:
		# 0,0,0,1/3,1/3,1/3,1/3,2/3,2/3,2/3,1,1,1,1,1  - this will lead to clusters of median methylation at 1/3 or 2/3 except for mostly unmethylated or methylated segments.
		median.meth = tapply(as.vector(runmean(Rle(values(m[queryHits(ov)])[,2]/values(m[queryHits(ov)])[, 1]), nSite.smoothing, endrule = "constant")), subjectHits(ov), median)
		median.meth = pmax(0, median.meth)
		if (!all.equal(as.numeric(names(T)), 1:length(segment_model))) {
			message("error in calculating methylation levels for PMDs")
		}	
		# Original code used median.meth to make a 'type' determination for the segment.  We will use pmeth instead
		type=ifelse(M/T<mMeth.classification[1],mMeth.classes[1], mMeth.classes[2])
		#type=ifelse(median.meth<mMeth.classification[1],mMeth.classes[1], mMeth.classes[2])
		#type=ifelse(median.meth<mMeth.classification[1],0,ifelse(median.meth<mMeth.classification[2],1,ifelse(median.meth<mMeth.classification[3],2,3)))
		#type = c("UMR", "LMR")[as.integer(nCG < nCG.classification) + 1]

		DataFrame(nSites.segmentation, nSites, T, M, pmeth = M/T, median.meth = median.meth, type)
	}
	
	# Define a corresponding convenience function to save an annotated genomic ranges object
	# Function adapted from MethylSeekR package
	saveUMRLMRSegments <- function (segs, GRangesFilename = NULL, TableFilename = NULL) 
	{
		# Replaced nCG with nSites for generality
		if (!is.null(GRangesFilename)) 
			saveRDS(segs, GRangesFilename)
		if (!is.null(TableFilename)) {
			D = data.frame(chr = as.character(seqnames(segs)), start = start(segs), 
				end = end(segs), type = values(segs)$type, nSites.segmentation = values(segs)$nSites.segmentation, 
				nSites.seq = values(segs)$nSites, mean.meth = 100 * values(segs)$pmeth, 
				median.meth = 100 * values(segs)$median.meth)
			write.table(D, file = TableFilename, quote = FALSE, sep = "\t", 
			row.names = FALSE)
		}
	}
	
	# We need to relate the overlap segments with the reference genome, so first, uncapitalise the sequence names
	levels(M_M_overlaps@seqnames)=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM")
	M_M_overlaps@seqinfo@seqnames=levels(M_M_overlaps@seqnames)
	levels(L_M_overlaps@seqnames)=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
	L_M_overlaps@seqinfo@seqnames=levels(L_M_overlaps@seqnames)
	levels(L_L_overlaps@seqnames)=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
	L_L_overlaps@seqinfo@seqnames=levels(L_L_overlaps@seqnames)

	# For each mCG overlap between an anotation from one of the principal classes, and a mCG classified segment, estimate its average CHG methylation
	# Load in the all_samples CHG methylome
	meth_CHG.gr <- readMethylome(FileName=paste0(project_id,"_CHG_TM_read_counts_across_samples.tsv"), seqLengths=sLengths)

	# There is a small chance that some short segments from the end of a gene will contain 0 CG/CHG sites (e.g. segments from HMM segmentation may start and finish midway between CG sites, whereas gene boundaries may be arbitrary)

	# adjust the segment model to remove any segments which do not overlap with CG sites
	M_M_overlaps = M_M_overlaps[unique(subjectHits(findOverlaps(meth.gr, M_M_overlaps)))]
	L_M_overlaps = L_M_overlaps[unique(subjectHits(findOverlaps(meth.gr, L_M_overlaps)))]
	L_L_overlaps = L_L_overlaps[unique(subjectHits(findOverlaps(meth.gr, L_L_overlaps)))]

	# adjust the segment model to remove any segments which do not overlap with CHG sites
	M_M_overlaps = M_M_overlaps[unique(subjectHits(findOverlaps(meth_CHG.gr, M_M_overlaps)))]
	L_M_overlaps = L_M_overlaps[unique(subjectHits(findOverlaps(meth_CHG.gr, L_M_overlaps)))]
	L_L_overlaps = L_L_overlaps[unique(subjectHits(findOverlaps(meth_CHG.gr, L_L_overlaps)))]
	
	# Add no. CG and CHG sites, and methylation levels, as annotations to *_*_overlaps
	values(M_M_overlaps) = cbind(values(M_M_overlaps),meth_by_segment(meth.gr, segment_model=M_M_overlaps, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths),meth_by_segment(meth_CHG.gr, segment_model=M_M_overlaps, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths))
	values(L_M_overlaps) = cbind(values(L_M_overlaps),meth_by_segment(meth.gr, segment_model=L_M_overlaps, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths),meth_by_segment(meth_CHG.gr, segment_model=L_M_overlaps, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths))
	values(L_L_overlaps) = cbind(values(L_L_overlaps),meth_by_segment(meth.gr, segment_model=L_L_overlaps, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths),meth_by_segment(meth_CHG.gr, segment_model=L_L_overlaps, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths))
	
	ggplot(data.frame(cbind(M_M_overlaps@elementMetadata@listData$median.meth, M_M_overlaps@elementMetadata@listData$median.meth.1))) + geom_point(aes(x=X1, y=X2))
	ggplot(data.frame(cbind(L_M_overlaps@elementMetadata@listData$median.meth, L_M_overlaps@elementMetadata@listData$median.meth.1))) + geom_point(aes(x=X1, y=X2))
	ggplot(data.frame(cbind(L_L_overlaps@elementMetadata@listData$median.meth, L_L_overlaps@elementMetadata@listData$median.meth.1))) + geom_point(aes(x=X1, y=X2))
	
	# This leads us back to thresholds.  Can we separate mCG GBM segments from mCG heterochromatic segments on the basis of average mCHG?

	# look at density distribution of mCHG proportions by genic segment
	
	# First we plot the distributions of CHG methylation level for each class of segment
	pdf(paste0(project_id,"_",meth_context,"_segment_wise_CHG_methylation_by_m_class.pdf"))
	print(ggplot(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))) + geom_point(aes(x=mCG, y=mCHG, colour=m_class)))
	dev.off()
	
	# Then we fit a mixture of 2 Gaussians to the combined density of mCHG levels of all segments from GBM genes and Heterochromatic genes/Transposons and find where they intersect.  This is where we draw the divide between GBM and Heterochromatic segments
	
	# Fit a mixture of 2 Gaussians to the mCHG density
	mCHG_mixmdl = normalmixEM(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG[!is.na(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG)], k=2)
	# Had to exclude sites with NA CHG methylation or the fitting algorithm exploded due to too little variance in one class
	
	
	# Print the characteristics of the fit curves
	cat(paste0("Fit lambda values (mixture proportions):\n"))
	mCHG_mixmdl$lambda
	# Schmitz data:  0.5899235 0.410.  Fishers p<0.05, Binomial p<0.005: 0.6828615 0.3171385
	# Becker data: 
	
	cat(paste0("Fit mu values (means):\n"))
	mCHG_mixmdl$mu
	# Schmitz data: 0.006048313 0.394139568.  Fishers p<0.05, Binomial p<0.005: 0.006111586 0.375212735
	# Becker data:  
	
	cat(paste0("Fit sigma values (st.devs):\n"))
	mCHG_mixmdl$sigma
	# Schmitz data: 0.002239083 0.213967850.  Fishers p<0.05, Binomial p<0.005: 0.002688758 0.221452772
	# Becker data:  
	
	# Find the intersection of the two larger components in the mixture - here we cut off

	#packages were installed and mixfn defined earlier
	#install.packages("rootSolve")
	#library(rootSolve)
	#mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
	#dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }

	mCHG_low_model=1
	mCHG_high_model=mCHG_low_model+1
	mCHG_cutoff=uniroot.all(mixfn, lower=0, upper=1, m1=mCHG_mixmdl$mu[mCHG_low_model], sd1=mCHG_mixmdl$sigma[mCHG_low_model], m2=mCHG_mixmdl$mu[mCHG_high_model], sd2=mCHG_mixmdl$sigma[mCHG_high_model], p1=mCHG_mixmdl$lambda[mCHG_low_model], p2=mCHG_mixmdl$lambda[mCHG_high_model])
	mCHG_cutoff
	# Schmitz CG data: 0.01414456.  Fishers p<0.05, Binomial p<0.005: 0.01577781
	# Becker CG data: 

	cat(paste0(length(L_M_overlaps@elementMetadata@listData$median.meth.1[L_M_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(L_M_overlaps@elementMetadata@listData$median.meth.1)," mCG segments among GBM genes have mCHG above ",mCHG_cutoff,".\n"))
	cat(paste0(length(L_L_overlaps@elementMetadata@listData$median.meth.1[L_L_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(L_L_overlaps@elementMetadata@listData$median.meth.1)," low-mCG segments among GBM genes have mCHG above ",mCHG_cutoff,".\n"))
	# Schmitz:
	#2053 of 22926 mCG segments among GBM genes have mCHG above 0.014144555663187.
	#63 of 2514 low-mCG segments among GBM genes have mCHG above 0.0141445571752365.
	# Fishers p<0.05, Binomial p<0.005: 
	#2768 of 36687 mCG segments among GBM genes have mCHG above 0.0157778100719126.
	#279 of 5825 low-mCG segments among GBM genes have mCHG above 0.0157778100719126.

	cat(paste0(length(M_M_overlaps@elementMetadata@listData$median.meth.1[M_M_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(M_M_overlaps@elementMetadata@listData$median.meth.1)," mCG segments among heterochromatic genes and transposons have mCHG above ",mCHG_cutoff,".\n"))
	#Schmitz:
	#12566 of 12869 mCG segments among heterochromatic genes and transposons have mCHG above 0.014144555663187.
	# Fishers p<0.05, Binomial p<0.005: 
	#13023 of 13436 mCG segments among heterochromatic genes and transposons have mCHG above 0.0157778100719126.
	
	pdf(paste0(project_id,"_",meth_context,"_mCHG_GBM_cutoff_density_fitted.pdf"))
	print(plot(mCHG_mixmdl,which=2, breaks=seq(0,1,0.01)))
	print(lines(density(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG[!is.na(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG)], lty=2, lwd=2)))
	dev.off()

	# Gaussians do not look like a great fit, and intersection of the two fitted curves is low, and excludes possibly too many segments.
	# Tried gamma distribution instead, and ignoring points with zero mCHG values (which cause gamma fitting to fail apparently), but couldn't get a nice result
	
	# How many GBM gene models are represented by the segments excluded through the fit Gaussians discriminator?

	# Define a set of genomic ranges combining low-mCHG and mCHG ranges from GBM annotated genes which have mCHG above the threshold
	Lo_CG_Hi_CHG_segments = c(L_M_overlaps[(!is.na(L_M_overlaps@elementMetadata@listData$median.meth.1)) & (L_M_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff),],L_L_overlaps[(!is.na(L_L_overlaps@elementMetadata@listData$median.meth.1)) & (L_L_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff),])
	
	# Capitalise chromosome names in overlap object to bring back in line with gene_body_loci.gr
	levels(Lo_CG_Hi_CHG_segments@seqnames)=c("CHR1", "CHR2", "CHR3", "CHR4", "CHR5")
	Lo_CG_Hi_CHG_segments@seqinfo@seqnames=levels(Lo_CG_Hi_CHG_segments@seqnames)

	
	# Find the overlap between these segments and GBM annotated genes
	GBM_Lo_hi_hits = findOverlaps(gene_body_loci.gr, Lo_CG_Hi_CHG_segments)
	gene_body_loci.gr[unique(GBM_Lo_hi_hits@from),]
	
	GBM_Lo_hi_hits = findOverlaps(gene_ranges, Lo_CG_Hi_CHG_segments)
	#gff.genes[unique(GBM_Lo_hi_hits@from),gene_ID]
	
	# These genes should predominantly be GBM genes (as we started with them) but there may be a few which happen to overlap a relevant GBM gene.  Check the distribution of m_class as a sense check:
	table(gff.genes[unique(GBM_Lo_hi_hits@from),]$m_class)

	# Schmitz data:
	#CG-poor Gene-body Methylated      Heterochromatic              Unknown         Unmethylated 
    #  2                 1960                   43                    4                 5 
	#Fishers p<0.05, Binomial p<0.005:
    #  3                 2662                   44                    4                	4 
	# Looks good
	
	# Reannotate GBM genes which overlap with Lo_CG_Hi_CHG_segments as Lo_CG_Hi_CHG
	gff.genes[unique(GBM_Lo_hi_hits@from),]$m_class = ifelse(gff.genes[unique(GBM_Lo_hi_hits@from),]$m_class=="Gene-body Methylated","Lo_CG_Hi_CHG",gff.genes[unique(GBM_Lo_hi_hits@from),]$m_class)
	
	# Check the CG/CHG profile now of the various gene classes
	pdf(paste0(project_id,"_",meth_context,"_gene-wise_CHG_vs_CGmethylation_corrected.pdf"))
	print(ggplot(merge(gff.genes, gff.genes.CHG, by="gene_ID"), aes(x=average_methylation.x, y=average_methylation.y)) + geom_point(aes(colour=m_class.x) ,alpha=0.4, size=1) + theme_minimal() + xlab("Average CG methylation") + ylab(paste0("Average CHG methylation")) + guides(color=guide_legend(title="Gene methylation status\n")))
	dev.off()

	#pdf(paste0(project_id,"_",meth_context,"_gene-wise_CHG_vs_CG_variable_sites_corrected.pdf"))
	#print(ggplot(merge(gff.genes, gff.genes.CHG, by="gene_ID"), aes(x=variable_count.x/(variable_count.x+all_M_count.x+all_U_count.x),y=average_methylation.y)) + geom_point(aes(colour=variance_methylation.x/average_methylation.x) ,alpha=0.4, size=1) + scale_colour_gradient(low="red", high="green") + theme_minimal() + xlab("Proportion of variable CG sites in gene") + ylab(paste0("Average methylation of CHG sites within gene")) +guides(color=guide_legend(title="Dispersion index\n of CG methylation\n in gene")))
	#dev.off()

	# Redo plot of the distributions of CHG methylation level for each class of gene
	pdf(paste0(project_id,"_",meth_context,"_gene_wise_CHG_methylation_density_by_m_class_corrected.pdf"))
	print(ggplot(merge(gff.genes, gff.genes.CHG, by="gene_ID")[,], aes(x=average_methylation.y, colour=m_class.x)) + geom_histogram(aes(y=0.1*..density..),binwidth=0.1)+ theme_minimal() + xlab("Average CHG methylation") + guides(color=guide_legend(title="Gene CG methylation status\n")) + facet_grid(m_class.x ~ .))
	dev.off()

	pdf(paste0(project_id,"_",meth_context,"_gene_length_by_m_class_corrected.pdf"))
	print(ggplot(merge(gff.genes, gff.genes.CHG, by="gene_ID")[,], aes(x=V5.x-V4.x+1, colour=m_class.x)) + geom_histogram(aes(y=100*..density..),binwidth=100)+ theme_minimal() + xlab("Gene-model length (nt)") + ylab("Density") + guides(color=guide_legend(title="Gene CG methylation status\n")) + facet_grid(m_class.x ~ .))
	dev.off()
	
	# Redo plot genes coloured by classification
	pdf(paste0(project_id,"_",meth_context,"_gene_methylation_level_and_status_corrected.pdf"))
	print(ggplot(gff.genes, aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=m_class) ,alpha=0.4, size=1) + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Gene methylation status\n")))
	dev.off()

	# Plot proportion of CG sites with changes by mCG status
	pdf(paste0(project_id,"_",meth_context,"_proportion_of_CG_sites_with_mCG_changes_by_m_class.pdf"))
	print(ggplot(gff.genes[,], aes(x=variable_count/(variable_count+all_M_count+all_U_count), colour=m_class)) + geom_histogram(aes(y=0.1*..density..),binwidth=.01)+ theme_minimal() + xlab("Proportion of CG sites with mCG changes") + ylab("Proportion of genes") + guides(color=guide_legend(title="Gene CG methylation status\n")) + facet_grid(m_class ~ .))
	dev.off()
	
	# Summarise the breakdown of genes by classification
	cat(paste0("Counts of genes by classification:\n"))
	table(gff.genes$m_class)

	# Count how many gene-body methylated genes have variable sites
	cat(paste0(nrow(gff.genes[(gff.genes$m_class=='Gene-body Methylated') & (gff.genes$variable_count>0),])," of ",nrow(gff.genes[gff.genes$m_class=='Gene-body Methylated',])," genes with gene body methylation (",nrow(gff.genes[(gff.genes$m_class=='Gene-body Methylated') & (gff.genes$variable_count>0),])/nrow(gff.genes[gff.genes$m_class=='Gene-body Methylated',]),") have one or more CG sites which vary.\n"))
	# Schmitz: 8840 of 13135 genes with gene body methylation (0.673011039208222) have one or more CG sites which vary.
	### This proportion increased slightly by removing Lo_CG_Hi_CHG genes
	# Fishers p<0.05, Binomial p<0.005: 
	# 9271 of 12583 genes with gene body methylation (0.736787729476277) have one or more CG sites which vary.
	# More sensitive method of detecting changing sites:
	# 12971 of 15264 genes with gene body methylation (0.849777253668763) have one or more CG sites which vary.
	
	cat(paste0(nrow(gff.genes[(gff.genes$m_class!='Gene-body Methylated') & (gff.genes$variable_count>0),])," of ",nrow(gff.genes[gff.genes$m_class!='Gene-body Methylated',])," genes without gene body methylation (",nrow(gff.genes[(gff.genes$m_class!='Gene-body Methylated') & (gff.genes$variable_count>0),])/nrow(gff.genes[gff.genes$m_class!='Gene-body Methylated',]),") have one or more CG sites which vary.\n"))
	#Schmitz: 1770 of 20206 genes without gene body methylation (0.0875977432445808) have one or more CG sites which vary.
	### This proportion also increased (~doubled) by removing Lo_CG_Hi_CHG genes, indicating that the filter is likely conservative (which we suspected)
	# Fishers p<0.05, Binomial p<0.005: 
	# 2770 of 20758 genes without gene body methylation (0.133442528181906) have one or more CG sites which vary.
	# More sensitive method of detecting changing sites:
	# 1484 of 18077 genes without gene body methylation (0.0820932676882226) have one or more CG sites which vary.
	
	cat(paste0(nrow(gff.genes[(gff.genes$m_class=='Unmethylated') & (gff.genes$variable_count>0),])," of ",nrow(gff.genes[gff.genes$m_class=='Unmethylated',])," unmethylated genes (",nrow(gff.genes[(gff.genes$m_class=='Unmethylated') & (gff.genes$variable_count>0),])/nrow(gff.genes[gff.genes$m_class=='Unmethylated',]),") have one or more CG sites which vary.\n"))
	# Schmitz: 657 of 14454 unmethylated genes (0.0454545454545455) have one or more CG sites which vary.
	# Fishers p<0.05, Binomial p<0.005: 
	# 885 of 14594 unmethylated genes (0.0606413594627929) have one or more CG sites which vary.
	# More sensitive method of detecting changing sites:
	# 1213 of 14604 unmethylated genes (0.0830594357710216) have one or more CG sites which vary.
	
	# Plot the proportion of variable sites against length for non-heterochromatic genes
	pdf(paste0(project_id,"_",meth_context,"_proportion_of_variable_sites_GBM_genes_corrected.pdf"))
	print(ggplot(gff.genes[gff.genes$m_class=="Gene-body Methylated",], aes(x=V5-V4, y=average_methylation)) + geom_point(aes(colour=variable_count/(variable_count+all_M_count+all_U_count)),alpha=0.4, size=1)  + scale_colour_gradient(low="green", high="red", trans="log", breaks=c(0.001,0.01,0.1,1)) + theme_minimal() + xlab("Gene length (nt)") + ylab(paste0("Average methylation of ",meth_context," sites within gene")) + xlim(0,10000) + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff))+ guides(color=guide_legend(title="Proportion of\nvariable sites\n")))
	dev.off()


	# Reannotate genome into segments UMR, GBM, FMR, and calculate mean and variance of mCG and mCHG of each segment

	# UMRLMRsegments.gr contains unmethylated (UMR) and low-methylated (LMR) regions.  LMR are mostly part of GBM genes (I think!)
	# FMRsegments.gr contains fully methylated (FMR) and moderately methylated (GBM) regions. Moderately methylated genes are mostly part of GBM genes (I think!) 
	
	#Characteristics of UMRs:
	#low mean levels of mCG (<0.2ish)
	#low variance/mean of mCG (<0.25ish)
	#low mean levels of mCHG (<?)
	#low variance/mean of mCHG (<?)

	#Characteristics of FMRs:
	#high levels of mCG (>0.75ish)
	#low variance/mean of mCG (<0.25ish)
	#high mean levels of mCHG (>?)
	#low variance/mean of mCHG (<?)

	#Characteristics of GBMs:
	#moderate levels of mCG (0.2-0.75ish)
	#high variance/mean of mCG (>0.25ish)
	#low mean levels of mCHG (<.01ish?)
	#low variance/mean of mCHG (<?)
	
	# Remake UMRLMRsegments.gr
	m.sel=0.4
	n.sel=3
	UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape_",m.sel,"_",n.sel,".pdf"))
	levels(UMRLMRsegments.gr@seqnames@values) = paste0(substr(as.character(levels(UMRLMRsegments.gr@seqnames@values)),1,1),toupper(substr(as.character(levels(UMRLMRsegments.gr@seqnames@values)),2,3)),substr(as.character(levels(UMRLMRsegments.gr@seqnames@values)),4,4))
	UMRLMRsegments.gr@seqinfo@seqnames = levels(UMRLMRsegments.gr@seqnames@values)
	cat(paste0("UMRLMRsegments.gr contains ",length(UMRLMRsegments.gr)," segments.\n"))
	# 30355
	
	# Remake FMRsegments.gr
	FMRsegments.gr=setdiff(genome_loci.gr,UMRLMRsegments.gr)
	cat(paste0("FMRsegments.gr contains ",length(FMRsegments.gr)," segments.\n"))
	# 30362
	
	# Reset chromosome names to mixed case
	levels(FMRsegments.gr@seqnames@values) = paste0(substr(as.character(levels(FMRsegments.gr@seqnames@values)),1,1),tolower(substr(as.character(levels(FMRsegments.gr@seqnames@values)),2,3)),substr(as.character(levels(FMRsegments.gr@seqnames@values)),4,4))
	FMRsegments.gr@seqinfo@seqnames = levels(FMRsegments.gr@seqnames@values)
	levels(UMRLMRsegments.gr@seqnames@values) = paste0(substr(as.character(levels(UMRLMRsegments.gr@seqnames@values)),1,1),tolower(substr(as.character(levels(UMRLMRsegments.gr@seqnames@values)),2,3)),substr(as.character(levels(UMRLMRsegments.gr@seqnames@values)),4,4))
	UMRLMRsegments.gr@seqinfo@seqnames = levels(UMRLMRsegments.gr@seqnames@values)

	# Adjust segments to remove mask_loci.gr
	FMRsegments.gr=setdiff(FMRsegments.gr, mask_loci.gr)
	cat(paste0("FMRsegments.gr contains ",length(FMRsegments.gr)," segments after removing ChrM, ChrC and 'Unknown' genes.\n"))
	# 30371
	UMRLMRsegments.gr=setdiff(UMRLMRsegments.gr, mask_loci.gr)
	cat(paste0("UMRLMRsegments.gr contains ",length(UMRLMRsegments.gr)," segments after removing ChrM, ChrC and 'Unknown' genes.\n"))
	# 31390

	# Adjust segments to remove any segments which don't overlap with any CG sites
	FMRsegments.gr = FMRsegments.gr[unique(subjectHits(findOverlaps(meth.gr, FMRsegments.gr)))]
	cat(paste0("FMRsegments.gr contains ",length(FMRsegments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 30344
	UMRLMRsegments.gr = UMRLMRsegments.gr[unique(subjectHits(findOverlaps(meth.gr, UMRLMRsegments.gr)))]
	cat(paste0("UMRLMRsegments.gr contains ",length(UMRLMRsegments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 31333
	
	# Adjust segments to remove any segments which don't overlap with any CHG sites
	FMRsegments.gr = FMRsegments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, FMRsegments.gr)))]
	cat(paste0("FMRsegments.gr contains ",length(FMRsegments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 30151
	UMRLMRsegments.gr = UMRLMRsegments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, UMRLMRsegments.gr)))]
	cat(paste0("UMRLMRsegments.gr contains ",length(UMRLMRsegments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 31333

	# Annotate segments with CG methylation values
	values(FMRsegments.gr) = cbind(values(FMRsegments.gr),meth_by_segment(meth.gr, segment_model=FMRsegments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.7), mMeth.classes=c("LMR","FMR")))
	values(UMRLMRsegments.gr) = cbind(values(UMRLMRsegments.gr),meth_by_segment(meth.gr, segment_model=UMRLMRsegments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","LMR")))
	
	# Annotate segments with CHG methylation values
	values(FMRsegments.gr) = cbind(values(FMRsegments.gr),meth_by_segment(meth_CHG.gr, segment_model=FMRsegments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.5), mMeth.classes=c("LMR","HMR")))
	values(UMRLMRsegments.gr) = cbind(values(UMRLMRsegments.gr),meth_by_segment(meth_CHG.gr, segment_model=UMRLMRsegments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.5), mMeth.classes=c("LMR","HMR")))
	
	# Merge the two sets of segments
	MRsegments.gr = c(FMRsegments.gr, UMRLMRsegments.gr)
	
	# Plot mCHG vs mCG for each mCG class of segment
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_vs_mCG_median_by_mCG_class.pdf"))
	print(ggplot(data.frame(values(MRsegments.gr)@listData), aes(x=median.meth, y=median.meth.1)) + geom_point(aes(colour=type)) + theme_minimal())
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_vs_mCG_mean_by_mCG_class.pdf"))
	print(ggplot(data.frame(values(MRsegments.gr)@listData), aes(x=pmeth, y=pmeth.1)) + geom_point(aes(colour=type)) + theme_minimal())
	dev.off()

	# Plot density of segment mCHG 
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_density.pdf"))
	print(ggplot(data.frame(values(MRsegments.gr)@listData), aes(x=median.meth.1)) + geom_density(aes(colour=type)) + theme_minimal())
	dev.off()
	
	# Fit a mixture of 2 Gaussians to the mCHG density
	 mCHG_mixmdl = normalmixEM(data.frame(values(MRsegments.gr)@listData)[!is.na(data.frame(values(MRsegments.gr)@listData)$median.meth.1),]$median.meth.1, k=2)
	# Had to exclude sites with NA CHG methylation or the fitting algorithm exploded	
	
	# Print the characteristics of the fit curves
	cat(paste0("Fit lambda values (mixture proportions):\n"))
	mCHG_mixmdl$lambda
	# Schmitz data:  0.8235399 0.1764601
	# Fishers p<0.05, Binomial p<0.005: 0.8232132 0.1767868
	# Becker data: 
	
	cat(paste0("Fit mu values (means):\n"))
	mCHG_mixmdl$mu
	# Schmitz data: 0.005611094 0.315166202
	# Fishers p<0.05, Binomial p<0.005: 0.005613074 0.315730014
	# Becker data:  
	
	cat(paste0("Fit sigma values (st.devs):\n"))
	mCHG_mixmdl$sigma
	# Schmitz data: 0.001962522 0.244004872
	# Fishers p<0.05, Binomial p<0.005: 0.001963114 0.243960731
	# Becker data:  
	
	# Find the intersection of the two larger components in the mixture - here we cut off

	#packages were installed and mixfn defined earlier
	#install.packages("rootSolve")
	#library(rootSolve)
	#mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
	#dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }

	mCHG_low_model=1
	mCHG_high_model=mCHG_low_model+1
	mCHG_cutoff=uniroot.all(mixfn, lower=0, upper=1, m1=mCHG_mixmdl$mu[mCHG_low_model], sd1=mCHG_mixmdl$sigma[mCHG_low_model], m2=mCHG_mixmdl$mu[mCHG_high_model], sd2=mCHG_mixmdl$sigma[mCHG_high_model], p1=mCHG_mixmdl$lambda[mCHG_low_model], p2=mCHG_mixmdl$lambda[mCHG_high_model])
	mCHG_cutoff
	# Schmitz CG data: 0.01299478
	# Fishers p<0.05, Binomial p<0.005: 0.01299959
	# Becker CG data: 

	cat(paste0(length(L_M_overlaps@elementMetadata@listData$median.meth.1[L_M_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(L_M_overlaps@elementMetadata@listData$median.meth.1)," mCG segments among GBM genes have mCHG above ",mCHG_cutoff,".\n"))
	cat(paste0(length(L_L_overlaps@elementMetadata@listData$median.meth.1[L_L_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(L_L_overlaps@elementMetadata@listData$median.meth.1)," low-mCG segments among GBM genes have mCHG above ",mCHG_cutoff,".\n"))
	# Schmitz:
	#2053 of 22926 mCG segments among GBM genes have mCHG above 0.014144555663187.
	#63 of 2514 low-mCG segments among GBM genes have mCHG above 0.0141445571752365.
	# Fishers p<0.05, Binomial p<0.005:
	#3509 of 36687 mCG segments among GBM genes have mCHG above 0.0129995872248819.
	#419 of 5825 low-mCG segments among GBM genes have mCHG above 0.0129995872248819.

	
	cat(paste0(length(M_M_overlaps@elementMetadata@listData$median.meth.1[M_M_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(M_M_overlaps@elementMetadata@listData$median.meth.1)," mCG segments among heterochromatic genes and transposons have mCHG above ",mCHG_cutoff,".\n"))
	#Schmitz:
	#12566 of 12869 mCG segments among heterochromatic genes and transposons have mCHG above 0.014144555663187.
	# Fishers p<0.05, Binomial p<0.005:
	#13047 of 13436 mCG segments among heterochromatic genes and transposons have mCHG above 0.0129995872248819.
	
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_cutoff_density_fitted.pdf"))
	print(plot(mCHG_mixmdl,which=2, breaks=seq(0,1,0.001)))
	print(lines(density(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG[!is.na(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG)], lty=2, lwd=2)))
	dev.off()

	# Plot mCG vs. segment size, with segments coloured by mCHG 
	
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCG_vs_length.pdf"))
	print(ggplot(cbind(data.frame(values(MRsegments.gr)@listData),mCHG=ifelse(values(MRsegments.gr)@listData$median.meth.1>mCHG_cutoff,"High","Low"), seg_width=MRsegments.gr@ranges@width), aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=mCHG)) + coord_trans(x = "log10") + geom_hline(aes(yintercept = heterochromatic_gene_coverage_cutoff)) + theme_minimal())
	
	print(ggplot(cbind(data.frame(values(MRsegments.gr)@listData),mCHG=values(MRsegments.gr)@listData$median.meth.1, seg_width=MRsegments.gr@ranges@width), aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=mCHG), alpha=0.1) + coord_trans(x = "log10") + theme_minimal() + geom_hline(aes(yintercept = 0.2)) + facet_wrap(~ ifelse(mCHG>mCHG_cutoff*2,1,2), ncol=1) +geom_density2d() +scale_colour_gradient(low="green", high="red", trans="log"))
	dev.off()
	

	### New direction for this analysis:
	# First build a combined methylome for the parental lines and use this to identify methylated segments within genes.  
	
	across_parents_TM=matrix(0, nrow=nrow(all_samples_meth_status), ncol=2)

	for(sample_name in sample_metadata[sample_metadata$Generation==3,]$Identifier) {
		if(opt$verbose) {cat(paste0("Adding methylation read counts to across-parents per-site table: ",sample_name,"\n"))}
		# Merge each of the sample counts into the relevant accumulation table, accounting for possible NA values

		across_parents_TM[,1]=ifelse(is.na(coverage_data[coverage_data$Sample==sample_name,]$Cov_C+coverage_data[coverage_data$Sample==sample_name,]$Cov_T),across_parents_TM[,1],across_parents_TM[,1]+coverage_data[coverage_data$Sample==sample_name,]$Cov_C+coverage_data[coverage_data$Sample==sample_name,]$Cov_T)

		across_parents_TM[,2]=ifelse(is.na(coverage_data[coverage_data$Sample==sample_name,]$Cov_C),across_parents_TM[,2],across_parents_TM[,2]+coverage_data[coverage_data$Sample==sample_name,]$Cov_C)

	} # end for each sample

	# Write out a copy of the methylome in the format that MethylSeekR likes to read in
	write.table(cbind(paste0(substr(as.character(all_reps_meth_status$Chromosome),1,1),tolower(substr(as.character(all_reps_meth_status$Chromosome),2,3)),substr(as.character(all_reps_meth_status$Chromosome),4,4)),all_reps_meth_status$Locus,across_parents_TM),file=paste0(project_id,"_",meth_context,"_TM_read_counts_across_parents.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

	# coverage_data version: Write out a copy of the methylome in the format that MethylSeekR likes to read in
	### NB. this version has a problem I think - coverage_data has one row per site per sample so the output ends up being repeated N times.  
	#write.table(cbind(paste0(substr(as.character(coverage_data$Chromosome),1,1),tolower(substr(as.character(coverage_data$Chromosome),2,3)),substr(as.character(coverage_data$Chromosome),4,4)),coverage_data$Locus,across_parents_TM),file=paste0(project_id,"_",meth_context,"_TM_read_counts_across_parents.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	### Here is a bodged version with the first Schmitz sample name hard coded (so that we can get meth_parents_CHH.gr without needing to make all_reps_meth_status):
	#write.table(cbind(paste0(substr(as.character(coverage_data[coverage_data$Sample=="SRR342347",]$Chromosome),1,1),tolower(substr(as.character(coverage_data[coverage_data$Sample=="SRR342347",]$Chromosome),2,3)),substr(as.character(coverage_data[coverage_data$Sample=="SRR342347",]$Chromosome),4,4)),coverage_data[coverage_data$Sample=="SRR342347",]$Locus,across_parents_TM),file=paste0(project_id,"_",meth_context,"_TM_read_counts_across_parents.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	
	
	
	# Write out a copy of the methylome in Bismark .cov format for visualising in SeqMonk
	#write.table(cbind(as.character(all_reps_meth_status$Chromosome),all_reps_meth_status$Locus,ifelse(meth_context=="CG",1,0)+all_reps_meth_status$Locus,across_parents_TM[,2]/across_parents_TM[,1],across_parents_TM[,2],across_parents_TM[,1]-across_parents_TM[,2]),file=paste0(project_id,"_",meth_context,"_CT_read_counts_across_parents.cov"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

	# coverage_data version: Write out a copy of the methylome in Bismark .cov format for visualising in SeqMonk
	write.table(cbind(as.character(coverage_data$Chromosome),coverage_data$Locus,ifelse(meth_context=="CG",1,0)+coverage_data$Locus,across_parents_TM[,2]/across_parents_TM[,1],across_parents_TM[,2],across_parents_TM[,1]-across_parents_TM[,2]),file=paste0(project_id,"_",meth_context,"_CT_read_counts_across_parents.cov"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	
	
	#### At this point, once we have written out parental methylomes in all 3 contexts, we depart from this script and do the actual segmentation in the script 5-methylation_calls_optimise_non_CG_segmentation_2018-03-14.R
	#### The following code is left in here but may not function quite the same.  However, it does illustrate how mixtures of Gaussians were fitted to segment methylation density distributions to pick cutoffs for mCG and mCHG to differentiate between TEM, GBM and UMR segments.  The values chosen below were then used in the above script.
	#### We pick up our analysis in this script a long way below.
	
	# Define a convenience function to convert chromosome names to mixed-case (e.g. Chr9, ChrC)
	mixedCaseChr <- function(s, strict = FALSE) {
		paste0(toupper(substr(s,1,1)),tolower(substr(s,2,3)),toupper(substr(s,4,4)))
	}

	# Set chromosome names to mixed case in genome segments object
	levels(genome_loci.gr@seqnames@values) = mixedCaseChr(as.character(levels(genome_loci.gr@seqnames@values)))
	genome_loci.gr@seqinfo@seqnames = levels(genome_loci.gr@seqnames@values)
	levels(mask_loci.gr@seqnames@values) = mixedCaseChr(as.character(levels(mask_loci.gr@seqnames@values)))
	mask_loci.gr@seqinfo@seqnames = levels(mask_loci.gr@seqnames@values)

	# Set up control loci based on the annotation that we can use for optimising segmentation parameters (TP/FP)
	
	methylated_loci.gr=reduce(makeGRangesFromDataFrame(df=rbind(gff.transposons[,c("V1","V4","V5")],gff.genes[gff.genes$m_class=="Heterochromatic",c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1"))
	gene_body_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))
	unmethylated_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))

	
	### Control point 2
	# This control point is so that mCHG and mCHH can be analysed, and parental methylomes saved, before the segmentation process proceeds
	save.image(file=paste0(project_id,"_",meth_context,"_control_point_2.RData"))
	# load (paste0(project_id,"_",meth_context,"_control_point_2.RData"))
		
		
	
	# We first want to segment purely on the basis of CHG methylation
	### N.B. This section is now obsolete, we went on, instead, to segment jointly in CHG and CHH segment to identify first-pass TEM segments

	# Read the CHG and CHH methylomes we prepared earlier using meth_context="CHG" and "CHH" as Granges objects, and sort the concatenated object so that it can be used to build valid segments
	meth_parents_CHG.gr <- readMethylome(FileName=paste0(project_id,"_CHG_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)
	#meth_parents_CHH.gr <- readMethylome(FileName=paste0(project_id,"_CHH_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)	
	#meth_parents_CHG_CHH.gr <- sort(sortSeqlevels(c(readMethylome(FileName=paste0(project_id,"_CHG_TM_read_counts_across_parents.tsv"), seqLengths=sLengths), readMethylome(FileName=paste0(project_id,"_CHH_TM_read_counts_across_parents.tsv"), seqLengths=sLengths))))
	
	# First we make the comparison between Unmethylated + GBM genes and Transposons + Heterochromatic genes
	
	# Create data frame to store the results of the optimisation:
	# For each value of m and n parameters, how many variable sites are overlapping with UMR and LMR segments? What proportion of supposed Transposons and Heterochromatic gene space overlaps with UMR or LMR segemtns (fp)?  What proportion of supposed Unmenthylated or GBM methylated gene space overlaps with UMR or LMR segments (tp)? 
	m_n_CHG_optimisation=data.frame(m=numeric(), n=numeric(), UMR=integer(), LMR=integer(), fp=numeric(), tp=numeric())

	# Loop through values for m and n parameters 
	#m_CHG.seq=seq(0.002,0.1, by=0.002)
	m_CHG.seq=seq(0.01,0.1, by=0.01)
	#n_CHG.seq=seq(1,20, by=1)
	n_CHG.seq=seq(20,200, by=20)
	for (m_CHG.sel in m_CHG.seq) {
		for (n_CHG.sel in n_CHG.seq) {
			#n.sel=4
			#UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=m_CHG.sel, nCpG.cutoff=n_CHG.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)
			#UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHH.gr, meth.cutoff=m_CHG.sel, nCpG.cutoff=n_CHG.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)
			UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=sort(sortSeqlevels(c(meth_parents_CHG.gr, meth_parents_CHH.gr))), meth.cutoff=m_CHG.sel, nCpG.cutoff=n_CHG.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)

			# Capitalise chromosome names in segmentation object
			levels(UMRLMRsegments_CHG.gr@seqnames)=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")
			UMRLMRsegments_CHG.gr@seqinfo@seqnames=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")

			# Find overlaps between segments and positive and negative 'control' loci (annotated genes and transposons)
			fp_hits=findOverlaps(methylated_loci.gr, UMRLMRsegments_CHG.gr)
			fp_overlaps <- pintersect(methylated_loci.gr[queryHits(fp_hits)], UMRLMRsegments_CHG.gr[subjectHits(fp_hits)])
			fp_overlap_prop <- sum(width(fp_overlaps)) / sum(width(methylated_loci.gr))
			tp_hits=findOverlaps(unmethylated_loci.gr, UMRLMRsegments_CHG.gr)
			tp_overlaps <- pintersect(unmethylated_loci.gr[queryHits(tp_hits)], UMRLMRsegments_CHG.gr[subjectHits(tp_hits)])
			tp_overlap_prop <- sum(width(tp_overlaps)) / sum(width(unmethylated_loci.gr))
			
			# Find overlaps between segments and variable loci
#			variable_sites_UMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments_CHG.gr[UMRLMRsegments_CHG.gr@elementMetadata$type=="UMR"])
#			variable_sites_LMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments_CHG.gr[UMRLMRsegments_CHG.gr@elementMetadata$type=="LMR"])
#			m_n_CHG_optimisation=rbind(m_n_CHG_optimisation, list(m=m_CHG.sel, n=n_CHG.sel, UMR=length(variable_sites_UMR_olaps@from), LMR=length(variable_sites_LMR_olaps@from), fp=fp_overlap_prop, tp=tp_overlap_prop))
			# How many variable sites were captured?
#			cat(paste0("m.sel=",m_CHG.sel," n.sel=",n_CHG.sel," Variable sites in UMRs: ",length(variable_sites_UMR_olaps@from),"  LMRs: ",length(variable_sites_LMR_olaps@from),"\n"))
		}
	}
	m_n_CHG_optimisation$variable_sites_captured=m_n_CHG_optimisation$UMR+m_n_CHG_optimisation$LMR

	# Plot ROC curves for each value of n, and estimate AUC
	pdf(paste0(project_id,"_",meth_context,"_MethylSeekR_m_n_CHG_optimisation_parents_ROC.pdf"))
	print(ggplot(m_n_CHG_optimisation, aes(x=fp, y=tp)) + geom_point() +geom_text(aes(label=paste0(m)),hjust=0, vjust=0) +geom_line() +facet_wrap(~ n))
	dev.off()
	
	# Estimate AUC from trapezium approximation
	cat(paste0("Estimated AUC values from ROC curves generated by varying m parameter for each value of n in MethylSeekR:\n"))
	best_n_CHG = 0
	prev_best_auc_CHG = 0
	for (n_CHG.sel in n_CHG.seq) {
		m_n_CHG_optimisation_auc=0
		prev_tp=0
		prev_fp=0
		for (m_CHG.sel in m_CHG.seq) {
			m_n_CHG_optimisation_auc = m_n_CHG_optimisation_auc + (1-m_n_CHG_optimisation[(m_n_CHG_optimisation$m==m_CHG.sel) & (m_n_CHG_optimisation$n==n_CHG.sel),]$fp)*(m_n_CHG_optimisation[(m_n_CHG_optimisation$m==m_CHG.sel) & (m_n_CHG_optimisation$n==n_CHG.sel),]$tp-prev_tp) + ((m_n_CHG_optimisation[(m_n_CHG_optimisation$m==m_CHG.sel) & (m_n_CHG_optimisation$n==n_CHG.sel),]$tp-prev_tp)*(m_n_CHG_optimisation[(m_n_CHG_optimisation$m==m_CHG.sel) & (m_n_CHG_optimisation$n==n_CHG.sel),]$fp-prev_fp))/2
			prev_tp=m_n_CHG_optimisation[(m_n_CHG_optimisation$m==m_CHG.sel) & (m_n_CHG_optimisation$n==n_CHG.sel),]$tp
			prev_fp=m_n_CHG_optimisation[(m_n_CHG_optimisation$m==m_CHG.sel) & (m_n_CHG_optimisation$n==n_CHG.sel),]$fp
		}
		if (m_n_CHG_optimisation_auc > prev_best_auc_CHG) {
			best_n_CHG = n_CHG.sel
			prev_best_auc_CHG = m_n_CHG_optimisation_auc
		}
		cat(paste0(" n=",n_CHG.sel," AUC=",m_n_CHG_optimisation_auc,"\n"))
	}
	cat(paste0("n=",best_n_CHG," maximises AUC (",prev_best_auc_CHG,").\n"))
	
	# Schmitz data:
	#n=1 AUC=0.915774648823515
	#n=2 AUC=0.917097379474582
	#n=3 AUC=0.9188342431869
	#n=4 AUC=0.921845282039692
	#n=5 AUC=0.92306137528743
	#n=6 AUC=0.923972510806274
	#n=7 AUC=0.924920580611287
	#n=8 AUC=0.925526574512847
	#n=9 AUC=0.926096549181427
	#n=10 AUC=0.926450397107005
	#n=11 AUC=0.926820666617607
	#n=12 AUC=0.927119227810704
	#n=13 AUC=0.927757065154047
	#n=14 AUC=0.928303328166869
	#n=15 AUC=0.928354627629721
	#n=16 AUC=0.92872082408907
	#n=17 AUC=0.928816451405229
	#n=18 AUC=0.928850563061825
	#n=19 AUC=0.929141101049031 ***
	#n=20 AUC=0.929047131061828
 
	# n=19 generates the largest AUC. 
	
	### What is the best threshold for mCHG cutoff between methylated and unmethylated?
	
	# Jaemyung Choi indicated that 10% is a good 'rule of thumb' cutoff for mCHG
	# We aim to identify a good cutoff visually, then programatically
		
	# Segment first using 0.9 threshold with various minimum segment lengths to visualise overall landscape.
	UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=0.9, nCpG.cutoff=60, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_60_0_9.pdf"))
	UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=0.9, nCpG.cutoff=best_n_CHG, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_",best_n_CHG,"_0_9.pdf"))
	UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=0.9, nCpG.cutoff=6, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_6_0_9.pdf"))
	UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=0.9, nCpG.cutoff=2, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_2_0_9.pdf"))
	
	# Plot density distribution of segment mCHG
	pdf(paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_",best_n_CHG,"_0_9_density.pdf"))
	print(ggplot(as.data.frame(UMRLMRsegments_CHG.gr), aes(x=pmeth)) + geom_histogram(binwidth=0.01))
	dev.off()

	# Visual inspection shows a clear minimum in the density at mCHG=0.5.  Rerun segmentation with 0.5 threshold
	UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=0.5, nCpG.cutoff=best_n_CHG, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_",best_n_CHG,"_0_5.pdf"))
	pdf(paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_",best_n_CHG,"_0_5_density.pdf"))
	print(ggplot(as.data.frame(UMRLMRsegments_CHG.gr), aes(x=pmeth)) + geom_histogram(binwidth=0.01))
	dev.off()

	# Fit mixture of two Gaussians to mCHG density to identify cutoff
	mCHG_mixmdl = normalmixEM(as.data.frame(UMRLMRsegments_CHG.gr)[(!is.na(as.data.frame(UMRLMRsegments_CHG.gr)$pmeth)),]$pmeth, k=3)
	# Had to exclude sites with NA CHG methylation or the fitting algorithm exploded	
	
	# Print the characteristics of the fit curves
	cat(paste0("Fit lambda values (mixture proportions):\n"))
	mCHG_mixmdl$lambda
	# Schmitz data:   
	#[1] 0.3646905 0.6353095
	# Becker data: 
	
	cat(paste0("Fit mu values (means):\n"))
	mCHG_mixmdl$mu
	# Schmitz data:
	#[1] 0.02390335 0.24271846
	# Becker data:  
	
	cat(paste0("Fit sigma values (st.devs):\n"))
	mCHG_mixmdl$sigma
	# Schmitz data: 
	#[1] 0.01584327 0.08767327
	# Becker data:  
	
	# Find the intersection of the two larger components in the mixture - here we cut off

	#packages were installed and mixfn defined earlier
	#install.packages("rootSolve")
	#library(rootSolve)
	#mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
	#dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }

	mCHG_low_model=2
	mCHG_high_model=mCHG_low_model+1
	segment_mCHG_cutoff=uniroot.all(mixfn, lower=0, upper=1, m1=mCHG_mixmdl$mu[mCHG_low_model], sd1=mCHG_mixmdl$sigma[mCHG_low_model], m2=mCHG_mixmdl$mu[mCHG_high_model], sd2=mCHG_mixmdl$sigma[mCHG_high_model], p1=mCHG_mixmdl$lambda[mCHG_low_model], p2=mCHG_mixmdl$lambda[mCHG_high_model])
	segment_mCHG_cutoff
	# Schmitz data: 0.1262584
	# Becker data: 

#	cat(paste0(length(L_M_overlaps@elementMetadata@listData$median.meth.1[L_M_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(L_M_overlaps@elementMetadata@listData$median.meth.1)," mCG segments among GBM genes have mCHG above ",mCHG_cutoff,".\n"))
#	cat(paste0(length(L_L_overlaps@elementMetadata@listData$median.meth.1[L_L_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(L_L_overlaps@elementMetadata@listData$median.meth.1)," low-mCG segments among GBM genes have mCHG above ",mCHG_cutoff,".\n"))
	# Schmitz:
	#2053 of 22926 mCG segments among GBM genes have mCHG above 0.014144555663187.
	#63 of 2514 low-mCG segments among GBM genes have mCHG above 0.0141445571752365.
	
#	cat(paste0(length(M_M_overlaps@elementMetadata@listData$median.meth.1[M_M_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(M_M_overlaps@elementMetadata@listData$median.meth.1)," mCG segments among heterochromatic genes and transposons have mCHG above ",mCHG_cutoff,".\n"))
	#Schmitz:
	#12566 of 12869 mCG segments among heterochromatic genes and transposons have mCHG above 0.014144555663187.
	
	pdf(paste0(project_id,"_",meth_context,"_mCHG_segmentation_mCHG_cutoff_density_fitted.pdf"))
	print(plot(mCHG_mixmdl,which=2, breaks=seq(0,1,0.001)))
	#print(lines(density(data.frame(rbind(cbind.data.frame(mCHG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG[!is.na(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG)], lty=2, lwd=2)))
	dev.off()
	
	# Segment parental mCHG using selected optimal values
	m.sel=segment_mCHG_cutoff
	n.sel=best_n_CHG
	UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_",best_n_CHG,"_",segment_mCHG_cutoff,".pdf"))
	saveUMRLMRSegments(segs=UMRLMRsegments_CHG.gr, GRangesFilename=paste0("UMRsLMRs_parents_CHG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_CHG_",m.sel,"_",n.sel,".tsv"))

	# Segment parental mCHG again, using 1/3 optimum minimum segment length to catch shorter segments
	m.sel=segment_mCHG_cutoff
	n.sel=6
	UMRLMRsegments_medium_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_",best_n_CHG,"_",segment_mCHG_cutoff,".pdf"))
	saveUMRLMRSegments(segs=UMRLMRsegments_CHG.gr, GRangesFilename=paste0("UMRsLMRs_parents_CHG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_CHG_",m.sel,"_",n.sel,".tsv"))

	# Segment parental mCHG again, using 1/9 optimum minimum segment length to catch shorter segments
	m.sel=segment_mCHG_cutoff
	n.sel=2
	UMRLMRsegments_fine_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_",best_n_CHG,"_",segment_mCHG_cutoff,".pdf"))
	saveUMRLMRSegments(segs=UMRLMRsegments_CHG.gr, GRangesFilename=paste0("UMRsLMRs_parents_CHG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_CHG_",m.sel,"_",n.sel,".tsv"))

	### We now combine the long, medium and short mCHG segments to form a composite mCHGome
	
	MG_segments.gr = setdiff(genome_loci.gr, mask_loci.gr)
	#LMR_segments.gr = setdiff(MG_segments.gr, UMRLMRsegments.gr)
	#MMR_segments.gr = setdiff(setdiff(MG_segments.gr, UMRLMRsegments_medium.gr), LMR_segments.gr)
	#SMR_segments.gr = setdiff(setdiff(MG_segments.gr, UMRLMRsegments_fine.gr), c(LMR_segments.gr, MMR_segments.gr))
	#UMR_segments.gr = setdiff(UMRLMRsegments.gr, c(SMR_segments.gr, MMR_segments.gr))
	#MR_segments.gr = setdiff(MG_segments.gr, UMR_segments.gr)

	LMR_CHG_segments.gr = setdiff(setdiff(genome_loci.gr, UMRLMRsegments_CHG.gr), mask_loci.gr)
	MMR_CHG_segments.gr = setdiff(setdiff(setdiff(genome_loci.gr, UMRLMRsegments_medium_CHG.gr), mask_loci.gr), LMR_CHG_segments.gr)
	SMR_CHG_segments.gr = setdiff(setdiff(setdiff(setdiff(genome_loci.gr, UMRLMRsegments_fine.gr), mask_loci.gr), LMR_CHG_segments.gr), MMR_CHG_segments.gr)
	UMR_CHG_segments.gr = setdiff(setdiff(setdiff(setdiff(genome_loci.gr, SMR_CHG_segments.gr), MMR_CHG_segments.gr), LMR_CHG_segments.gr), mask_loci.gr)
	mCHG_segments.gr = setdiff(MG_segments.gr, UMR_CHG_segments.gr)

	# Adjust segments to remove any segments which don't overlap with any CG sites
	# Counts here are first, using Fishers p<0.01, Binomial p<0.01; second, using  Fishers p<0.05, Binomial p<0.005; third, masking mCHG segments prior to segmenting mCG
	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 1065
	LMR_CHG_segments.gr = LMR_CHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, LMR_CHG_segments.gr)))]
	cat(paste0("LMR_CHG_segments.gr contains ",length(LMR_CHG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 4453
	MMR_CHG_segments.gr = MMR_CHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MMR_CHG_segments.gr)))]
	cat(paste0("MMR_CHG_segments.gr contains ",length(MMR_CHG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 2244
	SMR_CHG_segments.gr = SMR_CHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, SMR_CHG_segments.gr)))]
	cat(paste0("SMR_CHG_segments.gr contains ",length(SMR_CHG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 81215
	UMR_CHG_segments.gr = UMR_CHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, UMR_CHG_segments.gr)))]
	cat(paste0("UMR_CHG_segments.gr contains ",length(UMR_CHG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 85240
	mCHG_segments.gr = mCHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, mCHG_segments.gr)))]
	cat(paste0("mCHG_segments.gr contains ",length(mCHG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 83040
	
	# Annotate segments with CG methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=LMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))
	#values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MMR")))
	#values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=SMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","SMR")))
	#values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=UMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#replace values previously added: 
	#values(mCHG_segments.gr) = cbind(meth_by_segment(meth_parents_CG.gr, segment_model=mCHG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(segment_mCHG_cutoff), mMeth.classes=c("UMR","TEM")))
	
	# Adjust segments to remove any segments which don't overlap with any CHG sites
	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 1047
	LMR_CHG_segments.gr = LMR_CHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, LMR_CHG_segments.gr)))]
	cat(paste0("LMR_CHG_segments.gr contains ",length(LMR_CHG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 4453
	MMR_CHG_segments.gr = MMR_CHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MMR_CHG_segments.gr)))]
	cat(paste0("MMR_CHG_segments.gr contains ",length(MMR_CHG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 2244
	SMR_CHG_segments.gr = SMR_CHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, SMR_CHG_segments.gr)))]
	cat(paste0("SMR_CHG_segments.gr contains ",length(SMR_CHG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 81215
	UMR_CHG_segments.gr = UMR_CHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, UMR_CHG_segments.gr)))]
	cat(paste0("UMR_CHG_segments.gr contains ",length(UMR_CHG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 85240
	mCHG_segments.gr = mCHG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, mCHG_segments.gr)))]
	cat(paste0("mCHG_segments.gr contains ",length(mCHG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 83040
	
	# Annotate segments with CHG methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MG_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=LMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=SMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=UMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(mCHG_segments.gr) = cbind(values(mCHG_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=mCHG_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	
	# Merge the various distinct classes of segment
	segmentation_model_CHG.gr = c(LMR_CHG_segments.gr, MMR_CHG_segments.gr, SMR_CHG_segments.gr, UMR_CHG_segments.gr)
	
	# Define a convenience function to save a segmentation model. This is adapted from MethylSeekR package
	saveSegmentationModel = function (segs, GRangesFilename = NULL, TableFilename = NULL) 
	{
		if (!is.null(GRangesFilename)) 
			saveRDS(segs, GRangesFilename)
		if (!is.null(TableFilename)) {
			D = data.frame(chr = as.character(seqnames(segs)), start = start(segs), 
				end = end(segs), type = values(segs)$type)
        write.table(D, file = TableFilename, quote = FALSE, sep = "\t", 
            row.names = FALSE)
		}
	}
	
	#saveSegmentationModel(segs=LMR_segments.gr, GRangesFilename=paste0("LMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_LMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	#saveSegmentationModel(segs=MMR_segments.gr, GRangesFilename=paste0("MMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_MMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	#saveSegmentationModel(segs=SMR_segments.gr, GRangesFilename=paste0("SMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_SMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	#saveSegmentationModel(segs=UMR_segments.gr, GRangesFilename=paste0("UMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_UMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	#saveSegmentationModel(segs=MR_segments.gr, GRangesFilename=paste0("MRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_MR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))

	
	
	### Segmenting in CHG context was not bad, but we should get better resolution of segment boundaries by segmenting in CHH and CHG contexts jointly.  This script does that: 5-methylation_calls_optimise_non_CG_segmentationf_2018-03-14.R
	# We reload the segmentation model constructed in the joint CHG/CHH contexts:
	m_non_CG.sel=0.15
	best_n_non_CG=500
	
	UMR_non_CG_segments.gr = readRDS(paste0("UMRsLMRs_parents_non_CG_",m_non_CG.sel,"_",best_n_non_CG,".gr.rds"))	
		
	# Now segment the genomic regions in CG context, for segments not annotated as m_non_CG segments.
	# First segment with n>=1 site.  We can use the overlap between this and the m_non_CG segmentation to extend the m_non_CG segments to the nearest mCG site(s), where necessary.
	
	# Read the CG methylome data in as a Granges object
	meth_parents_CG.gr <- readMethylome(FileName=paste0(project_id,"_",meth_context,"_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)
	
	# Identify the parts of the parental mCGome which overlap with unmethylated parental mCHGome
	#potential_CG_segments.gr = meth_parents_CG.gr[unique(subjectHits(findOverlaps(UMRLMRsegments_CHG.gr, meth_parents_CG.gr)))]
	potential_CG_segments.gr = meth_parents_CG.gr[unique(subjectHits(findOverlaps(UMR_non_CG_segments.gr, meth_parents_CG.gr)))]
	
	# Segment CG methylome in non-m_non_CG loci
	
	# What should m.sel be?
	# Set m.sel to 99ile of C/(C+T) for "U" sites
	#m.sel=quantile(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T),.99)
	#cat(paste0("mCG proportion which captures 99% of sites 'U' across all samples is ",m.sel,".\n"))
	# Schmitz: mCG proportion which captures 99% of sites 'U' across all samples is 0.0769230769230769.
	# Fishers p<0.05, Binomial p<0.005: mCG proportion which captures 99% of sites 'U' across all samples is 0.0833333333333333.

	# Using pmeth, rather than median.meth with .99 threshold for unmethylated segments produces some bleed into UMR from SMR particularly.  .95 is a better cutoff:
	#m.sel=quantile(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T),.95)
	cat(paste0("mCG proportion which captures 95% of sites 'U' across all samples is ",m.sel,".\n"))
	# Schmitz: mCG proportion which captures 95% of sites 'U' across all samples is 0.0377358490566038.
	# Fishers p<0.05, Binomial p<0.005: mCG proportion which captures 95% of sites 'U' across all samples is 0.0384615384615385.

	#m.sel=0.08748551  # This value was derived previously (further on in the script here) as a sensible mCG cutoff for unmethylated short segments, by fitting mixture of Gaussians to distribution of mCG of short parental segments
	
	# Try m.sel=0.95 to check out whole methylome landscape
	m.sel=0.95
	n.sel=1  # We segment using minimum segment size of 1 so that we can extend our mCHG segments by a single CG site only, if necessary (minimise likelihood of extending a m_non_CG segment into an adjacent mCG but non-m_non_CG segment).
	UMRLMRsegments_CG.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_",n.sel,"_",m.sel,".pdf"))
	
	# Save a copy for visualisation
	saveUMRLMRSegments(segs=UMRLMRsegments_CG.gr, GRangesFilename=paste0("UMRsLMRs_parents_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_",m.sel,"_",n.sel,".tsv"))

	# Visual inspection of landscape shows clear cutoff at 0.4 between methylated and unmethylated segments among non-mCHG loci
	m.sel=0.4
	n.sel=1  # We segment using minimum segment size of 1 so that we can extend our mCHG segments by a single CG site only, if necessary (minimise likelihood of extending a mCHG segment into an adjacent mCG but non-mCHG segment).
	UMRLMRsegments_CG.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=UMRLMRsegments_CG.gr, GRangesFilename=paste0("UMRsLMRs_parents_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_",m.sel,"_",n.sel,".tsv"))
	
	# Find which mCG segments overlap with mCHG segments, and merge these segments with mCHG segments for a final set of mCHG segments congnicent of CG site loci/boundary of TEM segment in CG site-space:
	
	# Find the mCG segments, rather than the unmethylated, and filter by masked genome:
	mCG_segments_1.gr = setdiff(setdiff(genome_loci.gr, mask_loci.gr),UMRLMRsegments_CG.gr)
	length(mCG_segments_1.gr)
	# 41863 mCG segments

	# Adjust segments to remove any segments which don't overlap with any CG sites
	mCG_segments_1.gr = mCG_segments_1.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, mCG_segments_1.gr)))]
	length(mCG_segments_1.gr)
	# 41844 mCG segments
	
	# Add back the mCG values for the mCG segments, so it will save properly
	values(mCG_segments_1.gr) = cbind(values(mCG_segments_1.gr),meth_by_segment(meth_parents_CG.gr, segment_model=mCG_segments_1.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=mCG_segments_1.gr, GRangesFilename=paste0("UMRsLMRs_parents_mCG_masked_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_mCG_masked_",m.sel,"_",n.sel,".tsv"))
	
	# Find the mCG segments which overlap with m_non_CG segments:
	#mCHG_mCG_segments.gr=mCG_segments_1.gr[unique(subjectHits(findOverlaps(setdiff(setdiff(genome_loci.gr, mask_loci.gr),UMRLMRsegments_CHG.gr), mCG_segments_1.gr)))]
	m_non_CG_mCG_segments.gr=mCG_segments_1.gr[unique(subjectHits(findOverlaps(m_non_CG_segments.gr, mCG_segments_1.gr)))]
	length(m_non_CG_mCG_segments.gr)
	# 2683 mCG segments have an overlap with m_non_CG segments

	# Add back the mCG values for the mCG overlap segments, so it will save properly
	values(m_non_CG_mCG_segments.gr) = cbind(values(m_non_CG_mCG_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=m_non_CG_mCG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=m_non_CG_mCG_segments.gr, GRangesFilename=paste0("UMRsLMRs_parents_non_CG_CG_overlap_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_non_CG_CG_overlap_",m.sel,"_",n.sel,".tsv"))

	# Take the Union of the m_non_CG best fit segments (subtracted from masked genome) and those mCG 1 0.4 segments which overlap them, to extend the m_non_CG best fit segments to the next closest CG site where appropriate:
	m_non_CG_segments.gr = reduce(union(setdiff(setdiff(genome_loci.gr, UMR_non_CG_segments.gr), mask_loci.gr), m_non_CG_mCG_segments.gr))
	length(m_non_CG_segments.gr)
	# 8492 segments
	
	# Adjust segments to remove any segments which don't overlap with any CG sites
	m_non_CG_segments.gr = m_non_CG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, m_non_CG_segments.gr)))]
	length(m_non_CG_segments.gr)
	# 8286 m_non_CG_segments

	# Add back the mCG values for the mCHG overlap segments, so it will save properly
	values(m_non_CG_segments.gr) = cbind(values(m_non_CG_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=m_non_CG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=m_non_CG_segments.gr, GRangesFilename=paste0("UMRsLMRs_parents_non_CG_segments.gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_non_CG_.tsv"))
	
	# Rebuild the non-m_non_CG fraction of the masked genome for subsequent segmentation in the CG context
	potential_CG_segments.gr = meth_parents_CG.gr[unique(subjectHits(findOverlaps(setdiff(setdiff(genome_loci.gr, m_non_CG_segments.gr), mask_loci.gr), meth_parents_CG.gr)))]
	
	
	### This is where we may have optimised the minimum segment length for segmenting the mCGome.  This code moved to separate script
	
	# Schmitz numbers:
	best_n = 9
	prev_best_auc = 0.936370415538498
			
			
	# Execute the preferred segmentation model, visualise and save the results
	#UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=0.012, nCpG.cutoff=9, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape.pdf"))

	
	# Here are a couple of other ideas for setting cutoffs between U and M, especially for short segments:
	#m.sel=0.025
	# Set m.sel to 99ile of C/(C+T) for "U" sites
	m.sel=quantile(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T),.99)
	cat(paste0("mCG proportion which captures 99% of sites 'U' across all samples is ",m.sel,".\n"))
	# Schmitz: mCG proportion which captures 99% of sites 'U' across all samples is 0.0769230769230769.
	# Fishers p<0.05, Binomial p<0.005: mCG proportion which captures 99% of sites 'U' across all samples is 0.0833333333333333.

	# Using pmeth, rather than median.meth with .99 threshold for unmethylated segments produces some bleed into UMR from SMR particularly.  .95 is a better cutoff:
	m.sel=quantile(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T),.95)
	cat(paste0("mCG proportion which captures 95% of sites 'U' across all samples is ",m.sel,".\n"))
	# Schmitz: mCG proportion which captures 95% of sites 'U' across all samples is 0.0377358490566038.
	# Fishers p<0.05, Binomial p<0.005: mCG proportion which captures 95% of sites 'U' across all samples is 0.0384615384615385.

	# Finally, we choose m.sel iteratively.  By building a finished model of long, medium and short segments in CG context, and assessing density distribution of mCG per segment, we identify a lot of short, lowly methylated segments, that are not credible.  Further on in this script, we fit a mixture of Gaussians to this distribution and use their crossing point to set the mCG threshold here:
	
	m.sel=0.09309094
	
	# First segment using the minimum segment which maximises AUC for ROC curve
	n.sel=best_n
	
	UMRLMRsegments.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, nCpG.smoothing = 1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape_parents_CG_",m.sel,"_",n.sel,".pdf"))
	#head(UMRLMRsegments.gr)
	#plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename=paste0("UMRsLMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))

	# We observe that n.sel=9 generates reasonable segments, not breaking up the majority of GBM loci.  However, it does miss short isolated regions of mCG which look like real GBM.
	# As a strategy to address this we carry out a new segmentation with the same m.sel, but with n.sel=3 then 1
	# We will use these finer segmentations to split Unmethylated segments, where necessary, but will not use it to further divide mCG segments.

	# Segment using minimum segment length of 3
	n.sel=3
	UMRLMRsegments_medium.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, nCpG.smoothing = 1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape_parents_CG_",m.sel,"_",n.sel,".pdf"))
	#head(UMRLMRsegments.gr)
	#plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments_medium.gr, GRangesFilename=paste0("UMRsLMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))

	# Segment using minimum segment length of 1
	n.sel=1
	UMRLMRsegments_fine.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, nCpG.smoothing = 1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape_parents_CG_",m.sel,"_",n.sel,".pdf"))
	#head(UMRLMRsegments.gr)
	#plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments_fine.gr, GRangesFilename=paste0("UMRsLMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))

	# We have generated Unmethylated Regions with minimum lengths of 9, 3 and 1 CG sites respectively (UMR9, UMR3 and UMR1)
	# Masked Genome (MG) = Genome - Mask
	# Long Methylated regions (LMR) = MG - (UMR9 + mCHG)
	# Medium Methylated regions (MMR) = (MG - (UMR3 + mCHG)) - LMR
	# Short methylated regions (SMR) = (MG - (UMR1 + mCHG)) - (LMR + MMR)
	# Unmethylated regions (UMR) = UMR9 - (SMR + MMR +mCHG)
	# Methylated regions (MR) = MG - UMR

	
	# Capitalise chromosome names in segmentation
	#levels(UMRLMRsegments.gr@seqnames)=toupper(levels(UMRLMRsegments.gr@seqnames))
	#UMRLMRsegments.gr@seqinfo@seqnames=levels(UMRLMRsegments.gr@seqnames)
	#levels(UMRLMRsegments_medium.gr@seqnames)=toupper(levels(UMRLMRsegments_medium.gr@seqnames))
	#UMRLMRsegments_medium.gr@seqinfo@seqnames=levels(UMRLMRsegments_medium.gr@seqnames)
	#levels(UMRLMRsegments_fine.gr@seqnames)=toupper(levels(UMRLMRsegments_fine.gr@seqnames))
	#UMRLMRsegments_fine.gr@seqinfo@seqnames=levels(UMRLMRsegments_fine.gr@seqnames)

	# Capitalise chromosome names in mask_loci.gr
	#levels(mask_loci.gr@seqnames)=toupper(levels(mask_loci.gr@seqnames))
	#mask_loci.gr@seqinfo@seqnames=levels(mask_loci.gr@seqnames)
	
	MG_segments.gr = setdiff(genome_loci.gr, mask_loci.gr)
	#LMR_segments.gr = setdiff(MG_segments.gr, UMRLMRsegments.gr)
	#MMR_segments.gr = setdiff(setdiff(MG_segments.gr, UMRLMRsegments_medium.gr), LMR_segments.gr)
	#SMR_segments.gr = setdiff(setdiff(MG_segments.gr, UMRLMRsegments_fine.gr), c(LMR_segments.gr, MMR_segments.gr))
	#UMR_segments.gr = setdiff(UMRLMRsegments.gr, c(SMR_segments.gr, MMR_segments.gr))
	#MR_segments.gr = setdiff(MG_segments.gr, UMR_segments.gr)
	MG_non_CG_segments.gr = setdiff(setdiff(genome_loci.gr, m_non_CG_segments.gr), mask_loci.gr)
	#LMR_segments.gr = setdiff(MG_non_CG_segments.gr, UMRLMRsegments.gr)
	LMR_segments.gr = setdiff(setdiff(setdiff(genome_loci.gr, UMRLMRsegments.gr), m_non_CG_segments.gr), mask_loci.gr)
	#MMR_segments.gr = setdiff(setdiff(MG_non_CG_segments.gr, UMRLMRsegments_medium.gr), LMR_segments.gr)
	MMR_segments.gr = setdiff(setdiff(setdiff(setdiff(genome_loci.gr, UMRLMRsegments_medium.gr), m_non_CG_segments.gr), mask_loci.gr), LMR_segments.gr)
	#SMR_segments.gr = setdiff(setdiff(MG_non_CG_segments.gr, UMRLMRsegments_fine.gr), c(LMR_segments.gr, MMR_segments.gr))
	SMR_segments.gr = setdiff(setdiff(setdiff(setdiff(setdiff(genome_loci.gr, UMRLMRsegments_fine.gr), m_non_CG_segments.gr), mask_loci.gr), LMR_segments.gr), MMR_segments.gr)
	#UMR_segments.gr = setdiff(UMRLMRsegments.gr, c(SMR_segments.gr, MMR_segments.gr))
	#UMR_segments.gr = setdiff(setdiff(setdiff(setdiff(UMRLMRsegments.gr, SMR_segments.gr), MMR_segments.gr), mCHG_segments.gr), mask_loci.gr)
	MR_segments.gr = setdiff(MG_non_CG_segments.gr, UMR_segments.gr)

	# Adjust segments to remove any segments which don't overlap with any CG sites
	# Counts here are first, using Fishers p<0.01, Binomial p<0.01; second, using  Fishers p<0.05, Binomial p<0.005; third, masking mCHG segments prior to segmenting mCG; fourth, segmenting CHG and CHH together, then masking
	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 1362, 1065, 1065, 1029
	LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, LMR_segments.gr)))]
	cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 16500, 16909, 13659, 12401
	MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MMR_segments.gr)))]
	cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 5294, 5315, 5898, 5053
	SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, SMR_segments.gr)))]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 5635, 8385, 9949, 7657
	MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MR_segments.gr)))]
	cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 27429, 30609, 29506, 25111

	# Annotate segments with CG methylation values
	values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=LMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))
	values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MMR")))
	values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=SMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","SMR")))
	values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#replace values previously added, and annotate these segments as TEM: 
	values(m_non_CG_segments.gr) = cbind(meth_by_segment(meth_parents_CG.gr, segment_model=m_non_CG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m_non_CG.sel), mMeth.classes=c("UMR","TEM")))

	# Adjust SMRs to remove any with too low mCG, after visual inspection of distributions found that SMRs were over-enriched for very short low-mCG segments
	# Visual observation of the distribution of mCG density, and good sense (we are averaging over six samples) suggests 0.25
	SMR_segments.gr = SMR_segments.gr[ as.vector((!is.na(SMR_segments.gr$pmeth)) & (SMR_segments.gr$pmeth>0.25)) ]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after removing segments with mCG<0.25.\n"))

	# Build set of UMRs corresponding to new set of LMRs, MMRs, SMRs, TEMs:
	UMR_segments.gr = setdiff(setdiff(setdiff(setdiff(setdiff(genome_loci.gr, SMR_segments.gr), MMR_segments.gr), LMR_segments.gr), m_non_CG_segments.gr), mask_loci.gr)
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments.\n"))
	# 26495
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 27188, 30383, 29026, 26389
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=UMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	meth_parents_CHG.gr <- readMethylome(FileName=paste0(project_id,"_CHG_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)
	
	# Adjust segments to remove any segments which don't overlap with any CHG sites
	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 1344, 1047, 1047, 1012
	LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, LMR_segments.gr)))]
	cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 16483, 16892, 13621, 12307
	MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MMR_segments.gr)))]
	cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 5205, 5221, 5750, 4874
	SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, SMR_segments.gr)))]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 4677, 6672, 7964, 5992
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 27185, 30374, 29000, 26268
	MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MR_segments.gr)))]
	cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 26365, 28785, 27335, 23173
	m_non_CG_segments.gr = m_non_CG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, m_non_CG_segments.gr)))]
	cat(paste0("m_non_CG_segments.gr contains ",length(m_non_CG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	#  ,  , 4437, 8193
	
	# Annotate segments with CHG methylation values
	values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MG_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=LMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=SMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=UMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(m_non_CG_segments.gr) = cbind(values(m_non_CG_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=m_non_CG_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Now the SLOW bit
	meth_parents_CHH.gr <- readMethylome(FileName=paste0(project_id,"_CHH_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)

	# Adjust segments to remove any segments which don't overlap with any CHH sites
	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 1012
	LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, LMR_segments.gr)))]
	cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 12303
	MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, MMR_segments.gr)))]
	cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 4870
	SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, SMR_segments.gr)))]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 3726
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 26268
	MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, MR_segments.gr)))]
	cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 23116
	m_non_CG_segments.gr = m_non_CG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, m_non_CG_segments.gr)))]
	cat(paste0("m_non_CG_segments.gr contains ",length(m_non_CG_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 8193
	
	# Annotate segments with CHH methylation values
	values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=MG_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=LMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=MMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=SMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=UMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=MR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(m_non_CG_segments.gr) = cbind(values(m_non_CG_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=m_non_CG_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	
	# Merge the various distinct classes of segment
	segmentation_model.gr = c(LMR_segments.gr, MMR_segments.gr, SMR_segments.gr, UMR_segments.gr, m_non_CG_segments.gr)
	
	# Define a convenience function to save a segmentation model. This is adapted from MethylSeekR package
	saveSegmentationModel = function (segs, GRangesFilename = NULL, TableFilename = NULL) 
	{
		if (!is.null(GRangesFilename)) 
			saveRDS(segs, GRangesFilename)
		if (!is.null(TableFilename)) {
			D = data.frame(chr = as.character(seqnames(segs)), start = start(segs), 
				end = end(segs), type = values(segs)$type)
        write.table(D, file = TableFilename, quote = FALSE, sep = "\t", 
            row.names = FALSE)
		}
	}
	
	saveSegmentationModel(segs=LMR_segments.gr, GRangesFilename=paste0("LMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_LMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	saveSegmentationModel(segs=MMR_segments.gr, GRangesFilename=paste0("MMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_MMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	saveSegmentationModel(segs=SMR_segments.gr, GRangesFilename=paste0("SMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_SMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	saveSegmentationModel(segs=UMR_segments.gr, GRangesFilename=paste0("UMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_UMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	saveSegmentationModel(segs=MR_segments.gr, GRangesFilename=paste0("MRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_MR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))


	# Plot mCHG vs mCG for each mCG class of segment
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_vs_mCG_median_by_mCG_class.pdf"))
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=median.meth, y=median.meth.1)) + geom_point(aes(colour=type)) + theme_minimal())
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_vs_mCG_mean_by_mCG_class.pdf"))
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=pmeth, y=pmeth.1)) + geom_point(aes(colour=type)) + theme_minimal())
	dev.off()

	# Plot density of segment mCHG 
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_median_mCHG_density.pdf"))
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=median.meth.1)) + geom_density(aes(colour=type)) + theme_minimal())
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mean_mCHG_density.pdf"))
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=pmeth.1)) + geom_density(aes(colour=type)) + theme_minimal())
	ggplot(data.frame(values(segmentation_model.gr)@listData)[data.frame(values(segmentation_model.gr)@listData)$type=="LMR",], aes(x=pmeth.1)) + geom_density(aes(colour=type)) + theme_minimal()
	dev.off()
	
	# Fit a mixture of 2 Gaussians to the mCHG density
	mCHG_mixmdl_median = normalmixEM(data.frame(values(segmentation_model.gr)@listData)[!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1),]$median.meth.1, k=2)
	mCHG_mixmdl_mean = normalmixEM(data.frame(values(segmentation_model.gr)@listData)[!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1),]$pmeth.1, k=2)
	mCHG_mixmdl_mean = normalmixEM(data.frame(values(segmentation_model.gr)@listData)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$type)) & (data.frame(values(segmentation_model.gr)@listData)$type=="LMR"),]$pmeth.1, k=2)

	# Had to exclude sites with NA CHG methylation or the fitting algorithm exploded	
	
	# Print the characteristics of the fit curves
	cat(paste0("Fit lambda values (mixture proportions) median then mean models:\n"))
	mCHG_mixmdl_median$lambda
	mCHG_mixmdl_mean$lambda
	# Schmitz data:   
	#[1] 0.8245796 0.1754204
	#[1] 0.8229448 0.1770552

	# Becker data: 
	
	cat(paste0("Fit mu values (means) median then mean models:\n"))
	mCHG_mixmdl_median$mu
	mCHG_mixmdl_mean$mu
	# Schmitz data:
	#[1] 0.005231171 0.272684494
	#[1] 0.006069077 0.288362640
	# Becker data:  
	
	cat(paste0("Fit sigma values (st.devs) median then mean models:\n"))
	mCHG_mixmdl_median$sigma
	mCHG_mixmdl_mean$sigma
	# Schmitz data: 
	#[1] 0.001853741 0.243537961
	#[1] 0.001818791 0.242231344

	# Becker data:  
	
	# Find the intersection of the two larger components in the mixture - here we cut off

	#packages were installed and mixfn defined earlier
	#install.packages("rootSolve")
	#library(rootSolve)
	#mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
	#dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }

	mCHG_low_model=1
	mCHG_high_model=mCHG_low_model+1
	mCHG_cutoff_median=uniroot.all(mixfn, lower=0, upper=1, m1=mCHG_mixmdl_median$mu[mCHG_low_model], sd1=mCHG_mixmdl_median$sigma[mCHG_low_model], m2=mCHG_mixmdl_median$mu[mCHG_high_model], sd2=mCHG_mixmdl_median$sigma[mCHG_high_model], p1=mCHG_mixmdl_median$lambda[mCHG_low_model], p2=mCHG_mixmdl_median$lambda[mCHG_high_model])
	mCHG_cutoff_mean=uniroot.all(mixfn, lower=0, upper=1, m1=mCHG_mixmdl_mean$mu[mCHG_low_model], sd1=mCHG_mixmdl_mean$sigma[mCHG_low_model], m2=mCHG_mixmdl_mean$mu[mCHG_high_model], sd2=mCHG_mixmdl_mean$sigma[mCHG_high_model], p1=mCHG_mixmdl_mean$lambda[mCHG_low_model], p2=mCHG_mixmdl_mean$lambda[mCHG_high_model])
	mCHG_cutoff_median
	mCHG_cutoff_mean
	# Schmitz CG data:  0.010539, 0.01293227
	# Becker CG data: 

#	cat(paste0(length(L_M_overlaps@elementMetadata@listData$median.meth.1[L_M_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(L_M_overlaps@elementMetadata@listData$median.meth.1)," mCG segments among GBM genes have mCHG above ",mCHG_cutoff,".\n"))
#	cat(paste0(length(L_L_overlaps@elementMetadata@listData$median.meth.1[L_L_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(L_L_overlaps@elementMetadata@listData$median.meth.1)," low-mCG segments among GBM genes have mCHG above ",mCHG_cutoff,".\n"))
	# Schmitz:
	#2053 of 22926 mCG segments among GBM genes have mCHG above 0.014144555663187.
	#63 of 2514 low-mCG segments among GBM genes have mCHG above 0.0141445571752365.
	
#	cat(paste0(length(M_M_overlaps@elementMetadata@listData$median.meth.1[M_M_overlaps@elementMetadata@listData$median.meth.1>mCHG_cutoff])," of ",length(M_M_overlaps@elementMetadata@listData$median.meth.1)," mCG segments among heterochromatic genes and transposons have mCHG above ",mCHG_cutoff,".\n"))
	#Schmitz:
	#12566 of 12869 mCG segments among heterochromatic genes and transposons have mCHG above 0.014144555663187.
	
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_cutoff_density_fitted_median.pdf"))
	print(plot(mCHG_mixmdl_median,which=2, breaks=seq(0,1,0.001)))
	print(lines(density(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG[!is.na(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$median.meth, mCHG=M_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$median.meth, mCHG=L_M_overlaps@elementMetadata@listData$median.meth.1, m_class="Gene-body Methylated")))$mCHG)], lty=2, lwd=2)))
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_cutoff_density_fitted_mean.pdf"))
	print(plot(mCHG_mixmdl_mean,which=2, breaks=seq(0,1,0.001)))
	print(lines(density(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$pmeth, mCHG=M_M_overlaps@elementMetadata@listData$pmeth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$pmeth, mCHG=L_M_overlaps@elementMetadata@listData$pmeth.1, m_class="Gene-body Methylated")))$mCHG[!is.na(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$pmeth, mCHG=M_M_overlaps@elementMetadata@listData$pmeth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$pmeth, mCHG=L_M_overlaps@elementMetadata@listData$pmeth.1, m_class="Gene-body Methylated")))$mCHG)], lty=2, lwd=2)))
	dev.off()

	# Show a table of counts of segments by type

	#table(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),]$type)
	table(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$type)
	# Schmitz data:
	#LMR   MMR   SMR   UMR 
	#16546  5201  6769 30811 
	#16508  5205     4  6668 30643

	# Using TEM segments separately derived from segmenting mCHG:
	#LMR   MMR    MR   SMR   TEM   UMR 
	#13146  5726     1  2940  4343 25645 
	# Using TEM segments separately derived from CHG/CHH combined:
	#11813  4850    23  3724  7823 26994 

	# 4,1 segments ended up with type="MR" (short segments with relatively high levels of methylation) so we set these to be "SMR":
	segmentation_model.gr@elementMetadata@listData$type[(!is.na(segmentation_model.gr@elementMetadata@listData$type)) & (segmentation_model.gr@elementMetadata@listData$type=="MR")]="SMR"

	### This section makes refinements to the segmentation model based on observed 'bugs'
	
	### This part was pre-empted by going back to the start and setting the m.sel threshold for CG segmentations as 0.09309094.  This setting was made on the basis of having previously gone through the code segment below, to find an empirical cutoff for excluding very low mCG short segments
	# 1. SMRs are strongly enriched for very short regions methylated at low levels (<~.25).  At this level, this means that the sites in question are fully methylated in less than two replicates among the parental set of 6.  When viewed in context these do not look convincing enough to allow to break a longer unmethylated region.  From the segment Mean mCG density plots, a similar but far smaller enrichment can be seen in MMRs, but on inspection of a couple of these, due to the minimum length requirement of 3 sites, they seem to be more convincing and are found in more of the lines.
	
	# Fit mixture of 3 Gaussians to mCG density on non-UMR segments (all methylated segments). The intended three components are:  TEMs, GBMs, noise SMRs
	mCG_mixmdl_mean = normalmixEM(as.data.frame(segmentation_model.gr)[(as.data.frame(segmentation_model.gr)$type!="UMR") & (!is.na(as.data.frame(segmentation_model.gr)$type)),]$pmeth, k=2)
	
	# Print the characteristics of the fit curves
	cat(paste0("Fit lambda values (mixture proportions):\n"))
	mCG_mixmdl_mean$lambda
	# Schmitz data:   
	#[1] 0.1043745 0.7456716 0.1499538

	# Becker data: 
	
	cat(paste0("Fit mu values (means):\n"))
	mCG_mixmdl_mean$mu
	# Schmitz data:
	#[1]  0.05220705 0.51706026 0.86424723
	# Becker data:  
	
	cat(paste0("Fit sigma values (st.devs):\n"))
	mCG_mixmdl_mean$sigma
	# Schmitz data: 
	#[1]  0.01350256 0.21413227 0.04437924

	# Becker data:  
	
	# Find the intersection of the two smaller mCG components in the mixture (GBMs and noise SMRs) - here we cut off

	#packages were installed and mixfn defined earlier
	#install.packages("rootSolve")
	#library(rootSolve)
	#mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
	#dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }

	mCG_low_model=1
	mCG_high_model=mCG_low_model+1
	mCG_cutoff_mean=uniroot.all(mixfn, lower=0, upper=1, m1=mCG_mixmdl_mean$mu[mCG_low_model], sd1=mCG_mixmdl_mean$sigma[mCG_low_model], m2=mCG_mixmdl_mean$mu[mCG_high_model], sd2=mCG_mixmdl_mean$sigma[mCG_high_model], p1=mCG_mixmdl_mean$lambda[mCG_low_model], p2=mCG_mixmdl_mean$lambda[mCG_high_model])
	mCG_cutoff_mean
	# Schmitz CG data:  0.01631153 0.08439491; 0.01556350 0.08748551; 0.01380286 0.09309094
	# Becker CG data: 

	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCG_cutoff_density_fitted_mean.pdf"))
	print(plot(mCG_mixmdl_mean,which=2, breaks=seq(0,1,0.01)))
	#print(lines(density(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$pmeth, mCHG=M_M_overlaps@elementMetadata@listData$pmeth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$pmeth, mCHG=L_M_overlaps@elementMetadata@listData$pmeth.1, m_class="Gene-body Methylated")))$mCHG[!is.na(data.frame(rbind(cbind.data.frame(mCG=M_M_overlaps@elementMetadata@listData$pmeth, mCHG=M_M_overlaps@elementMetadata@listData$pmeth.1, m_class="Heterochromatic"),cbind.data.frame(mCG=L_M_overlaps@elementMetadata@listData$pmeth, mCHG=L_M_overlaps@elementMetadata@listData$pmeth.1, m_class="Gene-body Methylated")))$mCHG)], lty=2, lwd=2)))
	dev.off()
	# Fit is not fantastic, but looks good enough for a relatively conservative cutoff to ignore low-mCG SMRs
	
	# N.B. The curves cross twice and we want the upper cross, so mCG_cutoff_mean[2]
	
	# Rebuild the segmentation model, removing SMRs with mCG<threshold
	#MG_segments.gr = setdiff(genome_loci.gr, mask_loci.gr)
	#LMR_segments.gr = setdiff(MG_segments.gr, UMRLMRsegments.gr)
	#MMR_segments.gr = setdiff(setdiff(MG_segments.gr, UMRLMRsegments_medium.gr), LMR_segments.gr)
	#SMR_segments.gr = setdiff(setdiff(MG_segments.gr, UMRLMRsegments_fine.gr), c(LMR_segments.gr, MMR_segments.gr))
	# Adjust segments to remove any segments which don't overlap with any CG sites
	#SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, SMR_segments.gr)))]
	#cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 8535
	#MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MMR_segments.gr)))]
	#cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 5312

	# Add the mCG values to SMRs so they can be used as a filter
	#values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=SMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","SMR")))

	
	# We set mCG_cutoff_mean here to the value previously found, and reiterate mCHG cutoffs to be used for classifying finished segments
	mCG_cutoff_mean = 0.09309094
	segment_mCHG_cutoff = 0.15
	mCHG_cutoff_mean = 0.15
	mCHG_cutoff_median = 0.15
	
	# Remove SMR segments with mCG below threshold
	#  Visual observation of the distribution of mCG density suggests 0.25
	SMR_segments.gr = SMR_segments.gr[SMR_segments.gr$pmeth>mCG_cutoff_mean]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after removing segments with mCG<",mCG_cutoff_mean[2],".\n"))
	# 4619, 5255

	# Remove the segment annotations again from SMR_segments.gr, so they don't interfere with the setdiff process
	#SMR_segments.gr@elementMetadata@listData = list()

	UMR_segments.gr = setdiff(UMRLMRsegments.gr, c(SMR_segments.gr, MMR_segments.gr))
	MR_segments.gr = setdiff(MG_segments.gr, UMR_segments.gr)

	# Adjust segments to remove any segments which don't overlap with any CG sites
	#MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MG_segments.gr)))]
	#cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 1362
	#LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, LMR_segments.gr)))]
	#cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 16952
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 26602, 24178
	MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 26883, 24397

	# Annotate segments with CG methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=LMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))
	#values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MMR")))
	#values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=SMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","SMR")))
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=UMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Adjust segments to remove any segments which don't overlap with any CHG sites
	#MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MG_segments.gr)))]
	#cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 1344, 1047
	#LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, LMR_segments.gr)))]
	#cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 16935, 16892
	#MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MMR_segments.gr)))]
	#cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 5218, 5221
	#SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, SMR_segments.gr)))]
	#cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 3944, 3896
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 26601, 26527, 24177
	MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MR_segments.gr)))]
	cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 26097, 26009, 24395

	# Annotate segments with CHG methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MG_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=LMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=SMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=UMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	
	# Merge the various distinct classes of segment
	segmentation_model.gr = c(LMR_segments.gr, MMR_segments.gr, SMR_segments.gr, UMR_segments.gr, mCHG_segments.gr)

	# Show a table of counts of segments by type

	table(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$type)
	Schmitz data:
	#LMR   MMR   SMR   UMR 
	#16546  5201  3942 26881    SMRs have reduced by almost 3000 and UMRs have also reduced correspondingly 
	#16508  5205  3894 26801	Changed test p-value thresholds do not change balance overmuch
	#LMR   MMR    SMR   TEM   UMR 
	#11813  4850  3747  7823 26994    Separated out TEM segments by segmenting mCHG and mCHH combined
	
		
	# Plot mCG vs. segment size, with segments coloured by type 
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCG_vs_length_by_type_revised.pdf"))

	# Plot distribution of segment mCG and mCHG by type
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Median mCG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Mean mCG"))

	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=type), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Median mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=type), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth.1)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth.1)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Mean mCHG"))

	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth, y=median.meth.1)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCG") + ylab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth, y=pmeth.1)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCG") + ylab("Segment Mean mCHG"))
	
	### Cut-offs to separate UMR, GBM, TEM:
	# What threshold of mCG would capture 80% of TMR LMRs?
	#mCHG_quantile = 0.2
	#TMR_mCG_cutoff = quantile(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)) & (values(segmentation_model.gr)@listData$median.meth.1 > mCHG_cutoff) & (values(segmentation_model.gr)@listData$type=="LMR"),]$median.meth,mCHG_quantile)
	#cat(paste0(100*(1-mCHG_quantile),"% of 'transposon-like' Long Methylated Regions have mCG>",TMR_mCG_cutoff,".\n"))
	# Schmitz: 80% of 'transposon-like' Long Methylated Regions have mCG>0.754069511898385.

	# What threshold of mCG would capture 80% of GBM LMRs?
	#mCG_quantile = 0.8
	#GBM_mCG_cutoff = quantile(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)) & (values(segmentation_model.gr)@listData$median.meth.1 <= mCHG_cutoff) & (values(segmentation_model.gr)@listData$median.meth > m.sel) & (values(segmentation_model.gr)@listData$type=="LMR"),]$median.meth,mCG_quantile)
	#cat(paste0(100*mCG_quantile,"% of 'gene-body methylated-like' Long Methylated Regions have mCG<",GBM_mCG_cutoff,".\n"))
	# Schmitz: 80% of 'gene-body methylated-like' Long Methylated Regions have mCG<0.673452805081518.

	# Try to fit 3 Gaussians model to pick mCHG cutoff to separate TEMs from GBMs:
	#plot(normalmixEM(data.frame(values(segmentation_model.gr)@listData)[(data.frame(values(segmentation_model.gr)@listData)$pmeth>m.sel) & (data.frame(values(segmentation_model.gr)@listData)$pmeth.1>0.01) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)),]$pmeth.1,k=3), which=2, breaks=seq(0,1,0.01))
	
	#m.sel=0.2
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=mCHG), alpha=0.1) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>mCHG_cutoff*2,1,ifelse(median.meth>m.sel,2,3)), ncol=1) +geom_density2d() +scale_colour_gradient(low="green", high="red"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=type), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>mCHG_cutoff_median,paste0("1. 'Transposon-like' segments:\n median mCHG > ",round(mCHG_cutoff_median, digits=4)),ifelse(median.meth>m.sel,paste0("2. 'Gene-body-like' segments:\n median mCG > ",round(m.sel, digits=4),"\n median mCHG <= ",round(mCHG_cutoff_median, digits=4)),paste0("3. Unmethylated segments:\n median mCG <= ",round(m.sel, digits=4),"\n median mCHG <= ",round(mCHG_cutoff_median, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Median mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=type), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>segment_mCHG_cutoff,paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
		
	dev.off()
	
	# Crystalise segment methylation status as a column in the segmentation object
	# mCHG cutoff is segment_mCHG_cutoff
	# mCG cutoff is m.sel
	segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$mCHG>segment_mCHG_cutoff,"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$pmeth>m.sel,"GBM","UMR"))
	# For Schmitz data, two of these segments have NA pmeth values
	#segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$mCHG>mCHG_cutoff_mean,"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$pmeth>m.sel,"GBM","UMR"))

	# Remove any unclassified segments
	length(segmentation_model.gr)
	# 55360
	
	cat(paste0("There are ",nrow(as.data.frame(segmentation_model.gr)[is.na(as.data.frame(segmentation_model.gr)$segment.mStatus),]), " segments with no methylation status classification.\n"))
	# Schmitz data: 5
	segmentation_model.gr = segmentation_model.gr[!is.na(as.data.frame(segmentation_model.gr)$segment.mStatus)]
	length(segmentation_model.gr)
	# 55355
	
	# Use the reduce function on each segment type, separately, to concatenate any adjoining segments of the same type
	segmentation_model.gr = c(reduce(segmentation_model.gr[segmentation_model.gr@elementMetadata@listData$segment.mStatus == "TEM"]), reduce(segmentation_model.gr[segmentation_model.gr@elementMetadata@listData$segment.mStatus == "GBM"]), reduce(segmentation_model.gr[segmentation_model.gr@elementMetadata@listData$segment.mStatus == "UMR"]))
	length(segmentation_model.gr)
	# 53507
	
	# Reannotate the segmentation model (again!)
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_CG.gr, segment_model=segmentation_model.gr.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	
	# Output a text version of the segmentation model for visualisation
	write.table(cbind("Chromosome"=as.character(as.data.frame(segmentation_model.gr)$seqnames),"Start"=as.data.frame(segmentation_model.gr)$start,"End"=as.data.frame(segmentation_model.gr)$end,"Type"=as.data.frame(segmentation_model.gr)$segment.mStatus),file=paste0(project_id,"_",meth_context,"_segmentation_model_revised.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	
	
	#### After carrying out the segmentation in all three contexts, we pick up here to analyse the results of the segmentation.

	# Load in the segmentation model we prepared earlier:
	
	# This version loads in the model prepared using the annotation-independent segmentation of mC data only:
	segmentation_model.gr = readRDS(paste0(project_id,"_",meth_context,"_segmentation_model_draft.rds"))

	# This version loads in the model prepared using the annotation to convert GBM segments to TEM if they don't overlap a gene model:
	#segmentation_model.gr = readRDS(paste0(project_id,"_",meth_context,"_segmentation_model_draft.rds"))
	
	# Reiterate cutoffs just in case
	# We set mCG_cutoff_mean here to the value previously found, and reiterate mCHG cutoffs to be used later for classifying finished segments
	mCG_cutoff_mean = 0.09309094    # derived by fitting mixture of Gaussians to density distribution of segment mCG of whole parental methylome segmented in CG context
	m.sel = 0.09309094
	segment_mCHG_cutoff = 0.1262584 # derived by fitting mixture of Gaussians to density distribution of segment mCHG of whole parental methylome segmented in CHG context
	mCHG_cutoff_mean = 0.1262584
	mCHG_cutoff_median = 0.1262584
	segment_mCHH_cutoff = 0.1      # this is somewhat arbitrary, and is only used to classify the segments as UMR/MR in CHH.  May come back to this, and give it a better cutoff later, if this is thought to be useful
	mCHH_cutoff_mean = 0.1
	
	# Projection of M->U and U->M sites into parental methylation segment model
	# Find overlaps between segments and CG sites with various patterns of variation
	
	# Capitalise chromosome names in segmentation_model.gr
	levels(segmentation_model.gr@seqnames)=toupper(levels(segmentation_model.gr@seqnames))
	segmentation_model.gr@seqinfo@seqnames=levels(segmentation_model.gr@seqnames)

	levels(variant_call_ranges@seqnames)=toupper(levels(variant_call_ranges@seqnames))
	variant_call_ranges@seqinfo@seqnames=levels(variant_call_ranges@seqnames)

	segmentation_model.gr$ID=row.names(as.data.frame(segmentation_model.gr))
	CG_segment_olaps = findOverlaps(CG_site_ranges, segmentation_model.gr)
	segmentation_model.gr$CG_site_count=table(segmentation_model.gr[subjectHits(CG_segment_olaps),]$ID)[segmentation_model.gr$ID]
	segmentation_model.gr$CG_site_count[is.na(segmentation_model.gr$CG_site_count)]=0
	variant_segment_olaps = findOverlaps(variant_call_ranges, segmentation_model.gr)
	segmentation_model.gr$variant_count=table(segmentation_model.gr[subjectHits(variant_segment_olaps),]$ID)[segmentation_model.gr$ID]
	segmentation_model.gr$variant_count[is.na(segmentation_model.gr$variant_count)]=0
	all_M_segment_olaps = findOverlaps(all_M_ranges, segmentation_model.gr)
	segmentation_model.gr$all_M_count=table(segmentation_model.gr[subjectHits(all_M_segment_olaps),]$ID)[segmentation_model.gr$ID]
	segmentation_model.gr$all_M_count[is.na(segmentation_model.gr$all_M_count)]=0
	all_U_segment_olaps = findOverlaps(all_U_ranges, segmentation_model.gr)
	segmentation_model.gr$all_U_count=table(segmentation_model.gr[subjectHits(all_U_segment_olaps),]$ID)[segmentation_model.gr$ID]
	segmentation_model.gr$all_U_count[is.na(segmentation_model.gr$all_U_count)]=0
	m_to_u_segment_olaps = findOverlaps(m_to_u_call_ranges, segmentation_model.gr)
	segmentation_model.gr$m_to_u_count=table(segmentation_model.gr[subjectHits(m_to_u_segment_olaps),]$ID)[segmentation_model.gr$ID]
	segmentation_model.gr$m_to_u_count[is.na(segmentation_model.gr$m_to_u_count)]=0
	u_to_m_segment_olaps = findOverlaps(u_to_m_call_ranges, segmentation_model.gr)
	segmentation_model.gr$u_to_m_count=table(segmentation_model.gr[subjectHits(u_to_m_segment_olaps),]$ID)[segmentation_model.gr$ID]
	segmentation_model.gr$u_to_m_count[is.na(segmentation_model.gr$u_to_m_count)]=0
	M_parent_segment_olaps = findOverlaps(M_parent_call_ranges, segmentation_model.gr)
	segmentation_model.gr$M_parent_count=table(segmentation_model.gr[subjectHits(M_parent_segment_olaps),]$ID)[segmentation_model.gr$ID]
	segmentation_model.gr$M_parent_count[is.na(segmentation_model.gr$M_parent_count)]=0
	U_parent_segment_olaps = findOverlaps(U_parent_call_ranges, segmentation_model.gr)
	segmentation_model.gr$U_parent_count=table(segmentation_model.gr[subjectHits(U_parent_segment_olaps),]$ID)[segmentation_model.gr$ID]
	segmentation_model.gr$U_parent_count[is.na(segmentation_model.gr$U_parent_count)]=0
	
	# Summarise the numbers of variant, all_M and all_U sites by segment methylation status:
	
	# Quick barplot version:
	barplot(prop.table(cbind("CG sites"=aggregate(CG_site_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$CG_site_count.Freq,"all_M"=aggregate(all_M_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$all_M_count.Freq,"all_U"=aggregate(all_U_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$all_U_count.Freq, "variable"=aggregate(variant_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$variant_count.Freq,"variableMU"=aggregate(m_to_u_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$m_to_u_count.Freq,"variableUM"=aggregate(u_to_m_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$u_to_m_count.Freq), 2))
	
	# Summarise in a table for nicer plotting:
	segmentation_model_site_summary=data.frame(rbind(cbind("Methylation"="All CG sites","mStatus"=aggregate(CG_site_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"Sites"=aggregate(CG_site_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$CG_site_count.Freq),cbind("Methylation"="Variable","mStatus"=aggregate(variant_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"Sites"=aggregate(variant_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$variant_count.Freq),cbind("Methylation"="VariableMU","mStatus"=aggregate(m_to_u_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"Sites"=aggregate(m_to_u_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$m_to_u_count.Freq),cbind("Methylation"="VariableUM","mStatus"=aggregate(u_to_m_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"Sites"=aggregate(u_to_m_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$u_to_m_count.Freq),cbind("Methylation"="Methylated","mStatus"=aggregate(all_M_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"Sites"=aggregate(all_M_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$all_M_count.Freq),cbind("Methylation"="Unmethylated","mStatus"=aggregate(all_U_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"Sites"=aggregate(all_U_count.Freq ~ segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$all_U_count.Freq)))

	segmentation_model_site_summary$Sites=as.numeric(as.character(segmentation_model_site_summary$Sites))

	# Plot the distribution of CG sites with different methylation status among the segment methylation states
	# This version does percentages
	pdf(paste0(project_id,"_",meth_context,"_methylation_status_by_segment_mStatus.pdf"))
	print(ggplot() + geom_bar(data=ddply(segmentation_model_site_summary, "Methylation", mutate, percent_Sites=Sites/sum(Sites) * 100), aes(x=Methylation, y=percent_Sites, fill=mStatus), stat="identity") + labs(x=paste0(meth_context," Sites (",project_id,")"), y="Percentage of Sites", fill="Segment Methylation Status\n") + theme_minimal())

	# This version does actual counts of sites
	print(ggplot() + geom_bar(data=segmentation_model_site_summary, aes(x=Methylation, y=Sites, fill=mStatus), stat="identity") + labs(x=paste0(meth_context," Sites (Schmitz et al, 2011)"), y="Number of Sites", fill="Segment Methylation Status\n") + facet_wrap(~Methylation, scales="free", ncol=6) + theme_minimal() + theme(strip.text.x = element_blank())) 

	# Rebuild the summary table with detail by segment type (as originally segmented) and by methylation status
	segmentation_model_site_summary=data.frame(rbind(cbind("Methylation"="All CG sites","mStatus"=aggregate(CG_site_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"segType"=aggregate(CG_site_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$type,"Sites"=aggregate(CG_site_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$CG_site_count.Freq),cbind("Methylation"="Variable","mStatus"=aggregate(variant_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"segType"=aggregate(variant_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$type,"Sites"=aggregate(variant_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$variant_count.Freq),cbind("Methylation"="VariableMU","mStatus"=aggregate(m_to_u_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"segType"=aggregate(m_to_u_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$type,"Sites"=aggregate(m_to_u_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$m_to_u_count.Freq),cbind("Methylation"="VariableUM","mStatus"=aggregate(u_to_m_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"segType"=aggregate(u_to_m_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$type,"Sites"=aggregate(u_to_m_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$u_to_m_count.Freq),cbind("Methylation"="Methylated","mStatus"=aggregate(all_M_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"segType"=aggregate(all_M_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$type,"Sites"=aggregate(all_M_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$all_M_count.Freq),cbind("Methylation"="Unmethylated","mStatus"=aggregate(all_U_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$segment.mStatus,"segType"=aggregate(all_U_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$type,"Sites"=aggregate(all_U_count.Freq ~type+segment.mStatus, data=as.data.frame(segmentation_model.gr), FUN=sum)$all_U_count.Freq)))

	segmentation_model_site_summary$Sites=as.numeric(as.character(segmentation_model_site_summary$Sites))

	#print(ggplot(segmentation_model_site_summary, aes(x=Methylation, y=Sites)) + geom_bar(aes(fill=mStatus), stat="identity") + labs(x=paste0(meth_context," Sites (Schmitz et al, 2011)"), y="Number of Sites", fill="Segment\nMethylation\nStatus\n") + facet_wrap(~segType+Methylation, scales="free", ncol=6) + theme_minimal())
	print(ggplot() + geom_bar(data=segmentation_model_site_summary, aes(x=segType, y=Sites, fill=mStatus), stat="identity") + labs(x=paste0(meth_context," Sites (Schmitz et al, 2011)"), y="Number of Sites", fill="Segment\nMethylation\nStatus\n") + facet_grid(Methylation~segType, scales="free") + theme_minimal())
	
	# Replot the segment lanscape coloured by numbers of variant sites
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width, as.data.frame(segmentation_model.gr)$variant_count.Freq)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=variant_count.Freq), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>segment_mCHG_cutoff,paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)) + scale_colour_gradient(low="green", high="red", trans="log"))
	dev.off()

	pdf(paste0(project_id,"_segmentation_model_summary_2019_02_18.pdf"))
	# Replot the segmentation landscape in terms of CG site length
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=nSites, y=pmeth)) + geom_point(aes(), size=0.1, alpha=0.25) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1) + xlab("Segment length (CG sites)") +ylab("Segment Mean mCG") + scale_x_continuous(limits=c(0,500)))
	# Replot a quantised version showing proportion of 'M' sites in segment
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=nSites, y=M_parent_count.Freq/(M_parent_count.Freq+U_parent_count.Freq))) + geom_point(aes(), size=0.1, alpha=0.25) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1) + xlab("Segment length (CG sites)") +ylab("Proportion of mCG sites in segment") + scale_x_continuous(limits=c(0,500)))
	
	# Plot the segment length and segment methylation spectra
	# Plot length density
	print(ggplot(data.frame(values(segmentation_model.gr)@listData)) +geom_density(aes(x=nSites)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1) + xlab("Segment length (CG sites)") + scale_x_continuous(limits=c(0,100)))	
	# Plot length histogram
	print(ggplot(data.frame(values(segmentation_model.gr)@listData)) +geom_histogram(aes(x=nSites), binwidth=1) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1) + xlab("Segment length (CG sites)") + scale_x_continuous(limits=c(0,20)))	
	# Plot pmeth density
	print(ggplot(data.frame(values(segmentation_model.gr)@listData)) +geom_density(aes(x=pmeth)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment mCG proportion") + scale_x_continuous(limits=c(0,1)))	
	# Plot prop meth sites density
	print(ggplot(data.frame(values(segmentation_model.gr)@listData)) +geom_density(aes(x=M_parent_count.Freq/(M_parent_count.Freq+U_parent_count.Freq))) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Proportion of mCG sites in segment") + scale_x_continuous(limits=c(0,1)))	
	# Plot prop meth sites histogram
	print(ggplot(data.frame(values(segmentation_model.gr)@listData)) +geom_histogram(aes(x=M_parent_count.Freq/(M_parent_count.Freq+U_parent_count.Freq))) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Proportion of mCG sites in segment") + scale_x_continuous(limits=c(-0.1,1.1)))	
	dev.off()
	 


	
	#ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=nSites, y=pmeth)) + geom_point(aes(), size=0.1, alpha=0.25) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1) + xlab("Segment length (CG sites)") +ylab("Segment Mean mCG") + scale_x_continuous(limits=c(0,500)) + scale_colour_gradient(low="green", high="red", trans="log")
	
	
	


	### Show mCHG distribution among variable sites in TEM segments in comparison with 'random sites'

	# These are the ranges corresponding to variable sites in TEM segments:
	View(as.data.frame(variant_call_ranges)[as.data.frame(variant_segment_olaps)[as.data.frame(variant_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM",]),]$queryHits,])
	
	CHG_window_size = 50 # window length in nt
	CHG_step_size=25 # No. nt to slide window per step
	CHG_window_range = 100 # No steps to take each side of central site
	
	window_summary=data.frame()

	central_sites=as.data.frame(variant_call_ranges)[as.data.frame(variant_segment_olaps)[as.data.frame(variant_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM",]),]$queryHits,]

	first_window = TRUE
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		if (first_window) {
			window_summary=rbind.data.frame(window_summary, c("Variable TEM CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)))
			names(window_summary) = c("category","position","mean_mCHG","var_mCHG")
			window_summary$category = as.character(window_summary$category)
			window_summary$position = as.integer(levels(window_summary$position))[window_summary$position]
			window_summary$mean_mCHG = as.numeric(levels(window_summary$mean_mCHG))[window_summary$mean_mCHG]
			window_summary$var_mCHG = as.numeric(levels(window_summary$var_mCHG))[window_summary$var_mCHG]
			first_window = FALSE
		} else {
			window_summary=rbind.data.frame(window_summary, setNames(c("Variable TEM CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		}
		#colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	# Sample an equal number of CG sites randomly, from among TEM segments, and repeat the analysis
	
	sample_size=nrow(central_sites)
	
	# Random 'all_M' sites in TEM segments:
	#central_sites=as.data.frame(all_M_ranges)[as.data.frame(all_M_segment_olaps)[as.data.frame(all_M_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM",]),]$queryHits,][sample(nrow(as.data.frame(all_M_ranges)[as.data.frame(all_M_segment_olaps)[as.data.frame(all_M_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM",]),]$queryHits,]), sample_size),]

	# Random CG sites in TEM segments:
	central_sites=as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM",]),]$queryHits,][sample(nrow(as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM",]),]$queryHits,]), sample_size),]
	
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		window_summary=rbind.data.frame(window_summary, setNames(c("Random TEM CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	# Random CG sites:
	central_sites=as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[,]),]$queryHits,][sample(nrow(as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[,]),]$queryHits,]), sample_size),]
	
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		window_summary=rbind.data.frame(window_summary, setNames(c("Random CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	#window_summary$category = as.character(window_summary$category)
	window_summary$position = as.integer(window_summary$position)
	window_summary$mean_mCHG = as.numeric(window_summary$mean_mCHG)
	window_summary$var_mCHG = as.numeric(window_summary$var_mCHG)

	# Plot mean and variance of mCHG in windows centred at sites in TEM segments where mCG varies
	pdf(paste0(project_id,"_",meth_context,"_TEM_variable_CG_sites_mCHG_windows_",CHG_window_size,".pdf"))
	#ggplot(window_summary, aes(x=position, y=mean_mCHG)) + geom_point() 
	print(ggplot() +geom_pointrange(data=window_summary, aes(x=position, y=mean_mCHG, ymin=mean_mCHG-var_mCHG/2, ymax=mean_mCHG+var_mCHG/2, colour=category)))
	dev.off()
	
	
	### Show mCHG distribution among variable sites in GBM segments in comparison with 'random sites'

	# These are the ranges corresponding to variable sites in GBM segments:
	View(as.data.frame(variant_call_ranges)[as.data.frame(variant_segment_olaps)[as.data.frame(variant_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM",]),]$queryHits,])
	
	CHG_window_size = 50
	CHG_step_size=25
	
	window_summary=data.frame()

	central_sites=as.data.frame(variant_call_ranges)[as.data.frame(variant_segment_olaps)[as.data.frame(variant_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM",]),]$queryHits,]

	first_window = TRUE
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		if (first_window) {
			window_summary=rbind.data.frame(window_summary, c("Variable GBM CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)))
			names(window_summary) = c("category","position","mean_mCHG","var_mCHG")
			window_summary$category = as.character(window_summary$category)
			window_summary$position = as.integer(levels(window_summary$position))[window_summary$position]
			window_summary$mean_mCHG = as.numeric(levels(window_summary$mean_mCHG))[window_summary$mean_mCHG]
			window_summary$var_mCHG = as.numeric(levels(window_summary$var_mCHG))[window_summary$var_mCHG]
			first_window = FALSE
		} else {
			window_summary=rbind.data.frame(window_summary, setNames(c("Variable GBM CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		}
		#colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	# Sample an equal number of CG sites randomly, from among GBM segments, and repeat the analysis
	
	sample_size=nrow(central_sites)
	
	# Random CG sites in GBM segments:
	central_sites=as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM",]),]$queryHits,][sample(nrow(as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM",]),]$queryHits,]), sample_size),]
	
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		window_summary=rbind.data.frame(window_summary, setNames(c("Random GBM CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	# Random CG sites:
	central_sites=as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[,]),]$queryHits,][sample(nrow(as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[,]),]$queryHits,]), sample_size),]
	
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		window_summary=rbind.data.frame(window_summary, setNames(c("Random CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	#window_summary$category = as.character(window_summary$category)
	window_summary$position = as.integer(window_summary$position)
	window_summary$mean_mCHG = as.numeric(window_summary$mean_mCHG)
	window_summary$var_mCHG = as.numeric(window_summary$var_mCHG)

	# Plot mean and variance of mCHG in windows centred at sites in TEM segments where mCG varies
	pdf(paste0(project_id,"_",meth_context,"_GBM_variable_CG_sites_mCHG_windows_",CHG_window_size,".pdf"))
	#ggplot(window_summary, aes(x=position, y=mean_mCHG)) + geom_point() 
	print(ggplot() +geom_pointrange(data=window_summary, aes(x=position, y=mean_mCHG, ymin=mean_mCHG-var_mCHG/2, ymax=mean_mCHG+var_mCHG/2, colour=category)))
	dev.off()
	
	
	### Show mCHG distribution among variable sites in GBM-like segments in comparison with 'random sites'

	# These are the ranges corresponding to variable sites in GBM-like segments:
	View(as.data.frame(variant_call_ranges)[as.data.frame(variant_segment_olaps)[as.data.frame(variant_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like",]),]$queryHits,])
	
	CHG_window_size = 50
	CHG_step_size=25
	
	window_summary=data.frame()

	central_sites=as.data.frame(variant_call_ranges)[as.data.frame(variant_segment_olaps)[as.data.frame(variant_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like",]),]$queryHits,]

	first_window = TRUE
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		if (first_window) {
			window_summary=rbind.data.frame(window_summary, c("Variable GBM-like CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)))
			names(window_summary) = c("category","position","mean_mCHG","var_mCHG")
			window_summary$category = as.character(window_summary$category)
			window_summary$position = as.integer(levels(window_summary$position))[window_summary$position]
			window_summary$mean_mCHG = as.numeric(levels(window_summary$mean_mCHG))[window_summary$mean_mCHG]
			window_summary$var_mCHG = as.numeric(levels(window_summary$var_mCHG))[window_summary$var_mCHG]
			first_window = FALSE
		} else {
			window_summary=rbind.data.frame(window_summary, setNames(c("Variable GBM-like CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		}
		#colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	# Sample an equal number of CG sites randomly, from among GBM-like segments, and repeat the analysis
	
	sample_size=nrow(central_sites)
	
	# Random CG sites in GBM-like segments:
	central_sites=as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like",]),]$queryHits,][sample(nrow(as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like",]),]$queryHits,]), sample_size),]
	
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		window_summary=rbind.data.frame(window_summary, setNames(c("Random GBM-like CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	# Random CG sites:
	central_sites=as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[,]),]$queryHits,][sample(nrow(as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[,]),]$queryHits,]), sample_size),]
	
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		window_summary=rbind.data.frame(window_summary, setNames(c("Random CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	#window_summary$category = as.character(window_summary$category)
	window_summary$position = as.integer(window_summary$position)
	window_summary$mean_mCHG = as.numeric(window_summary$mean_mCHG)
	window_summary$var_mCHG = as.numeric(window_summary$var_mCHG)

	# Plot mean and variance of mCHG in windows centred at sites in TEM segments where mCG varies
	pdf(paste0(project_id,"_",meth_context,"_GBM_like_variable_CG_sites_mCHG_windows_",CHG_window_size,".pdf"))
	#ggplot(window_summary, aes(x=position, y=mean_mCHG)) + geom_point() 
	print(ggplot() +geom_pointrange(data=window_summary, aes(x=position, y=mean_mCHG, ymin=mean_mCHG-var_mCHG/2, ymax=mean_mCHG+var_mCHG/2, colour=category)))
	dev.off()
	
	
	### Show mCHG distribution among variable sites in UMR segments in comparison with 'random sites'

	# These are the ranges corresponding to variable sites in UMR segments:
	View(as.data.frame(variant_call_ranges)[as.data.frame(variant_segment_olaps)[as.data.frame(variant_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR",]),]$queryHits,])
	
	CHG_window_size = 50
	CHG_step_size=25
	
	window_summary=data.frame()

	central_sites=as.data.frame(variant_call_ranges)[as.data.frame(variant_segment_olaps)[as.data.frame(variant_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR",]),]$queryHits,]

	first_window = TRUE
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		if (first_window) {
			window_summary=rbind.data.frame(window_summary, c("Variable UMR CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)))
			names(window_summary) = c("category","position","mean_mCHG","var_mCHG")
			window_summary$category = as.character(window_summary$category)
			window_summary$position = as.integer(levels(window_summary$position))[window_summary$position]
			window_summary$mean_mCHG = as.numeric(levels(window_summary$mean_mCHG))[window_summary$mean_mCHG]
			window_summary$var_mCHG = as.numeric(levels(window_summary$var_mCHG))[window_summary$var_mCHG]
			first_window = FALSE
		} else {
			window_summary=rbind.data.frame(window_summary, setNames(c("Variable UMR CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		}
		#colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	# Sample an equal number of CG sites randomly, from among UMR segments, and repeat the analysis
	
	sample_size=nrow(central_sites)
	
	# Random CG sites in UMR segments:
	central_sites=as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR",]),]$queryHits,][sample(nrow(as.data.frame(CG_site_ranges)[as.data.frame(CG_segment_olaps)[as.data.frame(CG_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR",]),]$queryHits,]), sample_size),]
	
	for (CHG_window_no in -CHG_window_range:CHG_window_range) {
		# for upstream and downstream:
		# make a genomic ranges for the window
		# calculate the mean mCHG for the window
		# store the resulting mCHG values 
		
		window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
		window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
		window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
		# Remove any windows which overlap ends of chromosomes
		if(CHG_window_no<0) {
			window_sites = window_sites[window_sites$end>=1,]
			window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
		} else {
			window_sites = window_sites[window_sites$start>=1,]
			window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
		}
		
		window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")
		
		# ensure all window sites overlap with one or more CHG sites
		window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

		window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
		cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), "variance: ",var(window_meths$pmeth, na.rm=TRUE),"\n"))
		window_summary=rbind.data.frame(window_summary, setNames(c("Random UMR CG sites", CHG_window_no*CHG_step_size, mean(window_meths$pmeth, na.rm=TRUE), var(window_meths$pmeth, na.rm=TRUE)), names(window_summary)))
		colnames(window_summary) = c("category","position","mean_mCHG","var_mCHG")
	}
	
	#window_summary$category = as.character(window_summary$category)
	window_summary$position = as.integer(window_summary$position)
	window_summary$mean_mCHG = as.numeric(window_summary$mean_mCHG)
	window_summary$var_mCHG = as.numeric(window_summary$var_mCHG)

	# Plot mean and variance of mCHG in windows centred at sites in TEM segments where mCG varies
	pdf(paste0(project_id,"_",meth_context,"_UMR_variable_CG_sites_mCHG_windows_",CHG_window_size,".pdf"))
	#ggplot(window_summary, aes(x=position, y=mean_mCHG)) + geom_point() 
	print(ggplot() +geom_pointrange(data=window_summary, aes(x=position, y=mean_mCHG, ymin=mean_mCHG-var_mCHG/2, ymax=mean_mCHG+var_mCHG/2, colour=category)))
	dev.off()
	
	
	
	### Making metaplots by sliding windows across aligned segment-class boundaries
	
	window_size = 100 # window length in nt
	step_size=10 # No. nt to slide window per step
	window_range = 500 # No steps to take each side of central site
	
	window_site_summary=data.frame()
	first_window = TRUE

	for (focus_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
		for (adjacent_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
	
			if (focus_mStatus != adjacent_mStatus) {
				#focus_mStatus = "UMR"
				#adjacent_mStatus = "TEM"
	
				# Create set of central sites to be the segments on the right hand side of the requested boundaries
				central_sites = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[segmentation_model.gr$segment.mStatus==adjacent_mStatus], maxgap=1L), end(first) == start(second) - 1L)@second)
				
				# Create corresponding set of left hand segments for each boundary
				left_hand_segments = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[segmentation_model.gr$segment.mStatus==adjacent_mStatus], maxgap=1L), end(first) == start(second) - 1L)@first)

				# We want to avoid accumulating data from outside of the segments under consideration.  Accordingly, create genomicranges corresponding to the adjoining halfs of each of the two adjoining segments.  We only use adjoining halfs of the segments to avoid edge effects from the other end of the segments
				relevant_sites_end = central_sites$end-round((central_sites$end-central_sites$start)/2)
				relevant_sites_start = left_hand_segments$end-round((left_hand_segments$end-left_hand_segments$start)/2)

				cat(paste0(nrow(central_sites)," boundaries between ",focus_mStatus," and ",adjacent_mStatus," segments.\n"))
				#	levels(central_sites@seqnames)=toupper(levels(central_sites@seqnames))
				#	central_sites@seqinfo@seqnames=levels(central_sites@seqnames)


				for (window_no in -window_range:window_range) {
					# for upstream and downstream:
					# make a genomic ranges for the window
					# calculate the No. CG sites and No. variable sites for the window
					# store the resulting values 
		
					window_sites=data.frame(cbind("seqnames"=as.character(central_sites$seqnames), "start"=central_sites$start + (window_no*step_size) - round((window_size/2) - 1.25), "end"=central_sites$start + (window_no*step_size) +round((window_size/2) + 0.25), "width"=window_size, "strand"=as.character(central_sites$strand)))

					window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
					window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
					window_sites$start = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$start<relevant_sites_start, relevant_sites_start, window_sites$start)))
					window_sites$end = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$end>relevant_sites_end, relevant_sites_end, window_sites$end)))

					# Remove any windows which did not overlap relvant sites at all
					window_sites = window_sites[!is.na(window_sites$start),]
					
					# Remove any windows which overlap ends of chromosomes
					if(window_no<0) {
						window_sites = window_sites[window_sites$end>=1,]
						window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
					} else {
						window_sites = window_sites[window_sites$start>=1,]
						window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
					}
		
					window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end", seqnames.field = "seqnames")

					window_sites.gr$ID=row.names(as.data.frame(window_sites.gr))
					#CG_segment_olaps = findOverlaps(CG_site_ranges, window_sites.gr)
					window_sites.gr$CG_site_count=table(window_sites.gr[subjectHits(findOverlaps(CG_site_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
					window_sites.gr$CG_site_count[is.na(window_sites.gr$CG_site_count)]=0
					#variant_segment_olaps = findOverlaps(variant_call_ranges, window_sites.gr)
					window_sites.gr$variant_count=table(window_sites.gr[subjectHits(findOverlaps(variant_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
					window_sites.gr$variant_count[is.na(window_sites.gr$variant_count)]=0
					#all_M_segment_olaps = findOverlaps(all_M_ranges, window_sites.gr)
					window_sites.gr$all_M_count=table(window_sites.gr[subjectHits(findOverlaps(all_M_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
					window_sites.gr$all_M_count[is.na(window_sites.gr$all_M_count)]=0
					#all_U_segment_olaps = findOverlaps(all_U_ranges, window_sites.gr)
					window_sites.gr$all_U_count=table(window_sites.gr[subjectHits(findOverlaps(all_U_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
					window_sites.gr$all_U_count[is.na(window_sites.gr$all_U_count)]=0
					#m_to_u_segment_olaps = findOverlaps(m_to_u_call_ranges, window_sites.gr)
					window_sites.gr$m_to_u_count=table(window_sites.gr[subjectHits(findOverlaps(m_to_u_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
					window_sites.gr$m_to_u_count[is.na(window_sites.gr$m_to_u_count)]=0
					#u_to_m_segment_olaps = findOverlaps(u_to_m_call_ranges, window_sites.gr)
					window_sites.gr$u_to_m_count=table(window_sites.gr[subjectHits(findOverlaps(u_to_m_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
					window_sites.gr$u_to_m_count[is.na(window_sites.gr$u_to_m_count)]=0
		
					# ensure all window sites overlap with one or more CHG sites
					#window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

					#window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
					cat(paste0("Window: ",window_no," mean variants: ",mean(window_sites.gr$variant_count, na.rm=TRUE), "variance: ",var(window_sites.gr$variant_count, na.rm=TRUE),"\n"))
			
					if (first_window) {
						window_site_summary=rbind.data.frame(window_site_summary, c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE), mean(window_sites.gr$m_to_u_count, na.rm=TRUE), var(window_sites.gr$m_to_u_count, na.rm=TRUE), mean(window_sites.gr$u_to_m_count, na.rm=TRUE), var(window_sites.gr$u_to_m_count, na.rm=TRUE), mean(window_sites.gr$CG_site_count, na.rm=TRUE), var(window_sites.gr$CG_site_count, na.rm=TRUE), nrow(window_sites)))
						names(window_site_summary) = c("category","position","mean_vmCG","var_vmCG","mean_vMU","var_vMU","mean_vUM","var_vUM","mean_nCG","var_nCG","seg_count")
						window_site_summary$category = as.character(window_site_summary$category)
						window_site_summary$position = as.integer(levels(window_site_summary$position))[window_site_summary$position]
						window_site_summary$mean_vmCG = as.numeric(levels(window_site_summary$mean_vmCG))[window_site_summary$mean_vmCG]
						window_site_summary$var_vmCG = as.numeric(levels(window_site_summary$var_vmCG))[window_site_summary$var_vmCG]
						window_site_summary$mean_vMU = as.numeric(levels(window_site_summary$mean_vMU))[window_site_summary$mean_vMU]
						window_site_summary$var_vMU = as.numeric(levels(window_site_summary$var_vMU))[window_site_summary$var_vMU]
						window_site_summary$mean_vUM = as.numeric(levels(window_site_summary$mean_vUM))[window_site_summary$mean_vUM]
						window_site_summary$var_vUM = as.numeric(levels(window_site_summary$var_vUM))[window_site_summary$var_vUM]
						window_site_summary$mean_nCG = as.numeric(levels(window_site_summary$mean_nCG))[window_site_summary$mean_nCG]
						window_site_summary$var_nCG = as.numeric(levels(window_site_summary$var_nCG))[window_site_summary$var_nCG]
						window_site_summary$seg_count = as.integer(levels(window_site_summary$seg_count))[window_site_summary$seg_count]
						first_window = FALSE
					} else {
						window_site_summary=rbind.data.frame(window_site_summary, setNames(c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE), mean(window_sites.gr$m_to_u_count, na.rm=TRUE), var(window_sites.gr$m_to_u_count, na.rm=TRUE), mean(window_sites.gr$u_to_m_count, na.rm=TRUE), var(window_sites.gr$u_to_m_count, na.rm=TRUE), mean(window_sites.gr$CG_site_count, na.rm=TRUE), var(window_sites.gr$CG_site_count, na.rm=TRUE), nrow(window_sites)), names(window_site_summary)))
					}
			#		window_summary=rbind.data.frame(window_summary, setNames(c("Random UMR CG sites", window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE)), names(window_summary)))
			#		colnames(window_summary) = c("category","position","mean_vmCG","var_vmCG")

				}
	
				#window_summary$category = as.character(window_summary$category)
				window_site_summary$position = as.integer(window_site_summary$position)
				window_site_summary$mean_vmCG = as.numeric(window_site_summary$mean_vmCG)
				window_site_summary$var_vmCG = as.numeric(window_site_summary$var_vmCG)
				window_site_summary$mean_vMU = as.numeric(window_site_summary$mean_vMU)
				window_site_summary$var_vMU = as.numeric(window_site_summary$var_vMU)
				window_site_summary$mean_vUM = as.numeric(window_site_summary$mean_vUM)
				window_site_summary$var_vUM = as.numeric(window_site_summary$var_vUM)
				window_site_summary$mean_nCG = as.numeric(window_site_summary$mean_nCG)
				window_site_summary$var_nCG = as.numeric(window_site_summary$var_nCG)
				window_site_summary$seg_count = as.integer(window_site_summary$seg_count)

				# Plot mean and variance of vmCG in windows centred at segment boundaries
				#pdf(paste0(project_id,"_",meth_context,"_",focus_mStatus,"_",adjacent_mStatus,"_variable_CG_sites_windows_",CHG_window_size,".pdf"))
				#print(ggplot() +geom_pointrange(data=window_site_summary, aes(x=position, y=mean_vmCG, ymin=mean_vmCG-var_vmCG/2, ymax=mean_vmCG+var_vmCG/2, colour=category)))
				#dev.off()
			}
		}
	}	
	
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
		
		

	# Reshape data into one row per quantity plotted per window (count data)
	paneltext=array()
	paneltext[1]=paste0("Mean increase in number of mCG sites per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Change"=rep("U->M sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$mean_vUM, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary))), cbind("Change"=rep("M->U sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary))), cbind("Change"=rep("Net changes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"= window_site_summary$mean_vUM-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary))), cbind("Change"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary))))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$change_count=as.numeric(levels(window_site_detail$change_count)[window_site_detail$change_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
	
	# Plot counts of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_count_by_segment boundary_windows_",window_size,".pdf"))
	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable UMR−GBM CG sites"),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab("UMR Segments | GBM Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		
	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable UMR−TEM CG sites"),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab("UMR Segments | TEM Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))

	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable UMR−GBM-like CG sites"),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab("UMR Segments | GBM-like Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()

	# Plot counts of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_count_by_segment boundary_windows_mirror_",window_size,".pdf"))
	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable GBM-UMR CG sites"),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab("GBM Segments | UMR Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		
	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable TEM-UMR CG sites"),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab("TEM Segments | UMR Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))

	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable GBM-like-UMR CG sites"),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab("GBM-like Segments | UMR Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()

	
	
	# Reshape data into one row per quantity plotted per window (rate data)
	paneltext=array()
	paneltext[1]=paste0("Mean increase in proportion of mCG sites per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Change"=rep("U->M sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$mean_vUM, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary))), cbind("Change"=rep("M->U sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary))), cbind("Change"=rep("Net changes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"= window_site_summary$mean_vUM-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary))), cbind("Change"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary))))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$change_count=as.numeric(levels(window_site_detail$change_count)[window_site_detail$change_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
	
	# Plot rates of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_prop_by_segment_boundary_windows_",window_size,".pdf"))
	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable UMR−GBM CG sites"),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab("UMR Segments | GBM Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	
	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable UMR−TEM CG sites"),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab("UMR Segments | TEM Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))

	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable UMR−GBM-like CG sites"),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab("UMR Segments | GBM-like Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()

	# Plot counts of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_prop_by_segment_boundary_windows_mirror_",window_size,".pdf"))
	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable GBM-UMR CG sites"),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab("GBM Segments | UMR Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	
	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable TEM-UMR CG sites"),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab("TEM Segments | UMR Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))

	print(ggplot() +geom_line(data=window_site_detail[window_site_detail$category %in% c("Variable GBM-like-UMR CG sites"),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab("GBM-like Segments | UMR Segments") + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()


	# Remake the boundary plots by dividing GBM, TEM and GBM-like segments into 3 quantiles each (short, medium, long), and finding boundaries 

	window_size = 100 # window length in nt
	step_size=10 # No. nt to slide window per step
	window_range = 500 # No steps to take each side of central site
	
	window_site_summary=data.frame()
	first_window = TRUE

	no_width_quantiles = 3
	for (width_quantile in 1:no_width_quantiles) {
		for (focus_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
			min_width = quantile(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_mStatus,]$width,(width_quantile-1)/no_width_quantiles)
			max_width = quantile(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_mStatus,]$width,(width_quantile)/no_width_quantiles)
			
			for (adjacent_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
	
				if (focus_mStatus != adjacent_mStatus) {
					#focus_mStatus = "UMR"
					#adjacent_mStatus = "TEM"
	
					# Create set of central sites to be the segments on the right hand side of the requested boundaries
					central_sites = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[(segmentation_model.gr$segment.mStatus==adjacent_mStatus) & (segmentation_model.gr@ranges@width>min_width) & (segmentation_model.gr@ranges@width<=max_width)], maxgap=1L), end(first) == start(second) - 1L)@second)
				
					# Create corresponding set of left hand segments for each boundary
					left_hand_segments = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[(segmentation_model.gr$segment.mStatus==adjacent_mStatus) & (segmentation_model.gr@ranges@width>min_width) & (segmentation_model.gr@ranges@width<=max_width)], maxgap=1L), end(first) == start(second) - 1L)@first)

					# We want to avoid accumulating data from outside of the segments under consideration.  Accordingly, create genomicranges corresponding to the adjoining halfs of each of the two adjoining segments.  We only use adjoining halfs of the segments to avoid edge effects from the other end of the segments
					relevant_sites_end = central_sites$end-round((central_sites$end-central_sites$start)/2)
					relevant_sites_start = left_hand_segments$end-round((left_hand_segments$end-left_hand_segments$start)/2)

					cat(paste0(nrow(central_sites)," boundaries between ",focus_mStatus," and ",adjacent_mStatus," segments.\n"))
					#	levels(central_sites@seqnames)=toupper(levels(central_sites@seqnames))
					#	central_sites@seqinfo@seqnames=levels(central_sites@seqnames)


					for (window_no in -window_range:window_range) {
						# for upstream and downstream:
						# make a genomic ranges for the window
						# calculate the No. CG sites and No. variable sites for the window
						# store the resulting values 
		
						window_sites=data.frame(cbind("seqnames"=as.character(central_sites$seqnames), "start"=central_sites$start + (window_no*step_size) - round((window_size/2) - 1.25), "end"=central_sites$start + (window_no*step_size) +round((window_size/2) + 0.25), "width"=window_size, "strand"=as.character(central_sites$strand)))

						window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
						window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
						window_sites$start = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$start<relevant_sites_start, relevant_sites_start, window_sites$start)))
						window_sites$end = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$end>relevant_sites_end, relevant_sites_end, window_sites$end)))

						# Remove any windows which did not overlap relvant sites at all
						window_sites = window_sites[!is.na(window_sites$start),]
					
						# Remove any windows which overlap ends of chromosomes
						if(window_no<0) {
							window_sites = window_sites[window_sites$end>=1,]
							window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
						} else {
							window_sites = window_sites[window_sites$start>=1,]
							window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
						}
		
						window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end", seqnames.field = "seqnames")

						window_sites.gr$ID=row.names(as.data.frame(window_sites.gr))
						#CG_segment_olaps = findOverlaps(CG_site_ranges, window_sites.gr)
						window_sites.gr$CG_site_count=table(window_sites.gr[subjectHits(findOverlaps(CG_site_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$CG_site_count[is.na(window_sites.gr$CG_site_count)]=0
						#variant_segment_olaps = findOverlaps(variant_call_ranges, window_sites.gr)
						window_sites.gr$variant_count=table(window_sites.gr[subjectHits(findOverlaps(variant_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$variant_count[is.na(window_sites.gr$variant_count)]=0
						#all_M_segment_olaps = findOverlaps(all_M_ranges, window_sites.gr)
						window_sites.gr$all_M_count=table(window_sites.gr[subjectHits(findOverlaps(all_M_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$all_M_count[is.na(window_sites.gr$all_M_count)]=0
						#all_U_segment_olaps = findOverlaps(all_U_ranges, window_sites.gr)
						window_sites.gr$all_U_count=table(window_sites.gr[subjectHits(findOverlaps(all_U_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$all_U_count[is.na(window_sites.gr$all_U_count)]=0
						#m_to_u_segment_olaps = findOverlaps(m_to_u_call_ranges, window_sites.gr)
						window_sites.gr$m_to_u_count=table(window_sites.gr[subjectHits(findOverlaps(m_to_u_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$m_to_u_count[is.na(window_sites.gr$m_to_u_count)]=0
						#u_to_m_segment_olaps = findOverlaps(u_to_m_call_ranges, window_sites.gr)
						window_sites.gr$u_to_m_count=table(window_sites.gr[subjectHits(findOverlaps(u_to_m_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$u_to_m_count[is.na(window_sites.gr$u_to_m_count)]=0
		
						# ensure all window sites overlap with one or more CHG sites
						#window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

						#window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
						cat(paste0("Window: ",window_no," mean variants: ",mean(window_sites.gr$variant_count, na.rm=TRUE), "variance: ",var(window_sites.gr$variant_count, na.rm=TRUE),"\n"))
			
						if (first_window) {
							window_site_summary=rbind.data.frame(window_site_summary, c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE), mean(window_sites.gr$m_to_u_count, na.rm=TRUE), var(window_sites.gr$m_to_u_count, na.rm=TRUE), mean(window_sites.gr$u_to_m_count, na.rm=TRUE), var(window_sites.gr$u_to_m_count, na.rm=TRUE), mean(window_sites.gr$CG_site_count, na.rm=TRUE), var(window_sites.gr$CG_site_count, na.rm=TRUE), nrow(window_sites), width_quantile))
							names(window_site_summary) = c("category","position","mean_vmCG","var_vmCG","mean_vMU","var_vMU","mean_vUM","var_vUM","mean_nCG","var_nCG","seg_count","width_quantile")
							window_site_summary$category = as.character(window_site_summary$category)
							window_site_summary$position = as.integer(levels(window_site_summary$position))[window_site_summary$position]
							window_site_summary$mean_vmCG = as.numeric(levels(window_site_summary$mean_vmCG))[window_site_summary$mean_vmCG]
							window_site_summary$var_vmCG = as.numeric(levels(window_site_summary$var_vmCG))[window_site_summary$var_vmCG]
							window_site_summary$mean_vMU = as.numeric(levels(window_site_summary$mean_vMU))[window_site_summary$mean_vMU]
							window_site_summary$var_vMU = as.numeric(levels(window_site_summary$var_vMU))[window_site_summary$var_vMU]
							window_site_summary$mean_vUM = as.numeric(levels(window_site_summary$mean_vUM))[window_site_summary$mean_vUM]
							window_site_summary$var_vUM = as.numeric(levels(window_site_summary$var_vUM))[window_site_summary$var_vUM]
							window_site_summary$mean_nCG = as.numeric(levels(window_site_summary$mean_nCG))[window_site_summary$mean_nCG]
							window_site_summary$var_nCG = as.numeric(levels(window_site_summary$var_nCG))[window_site_summary$var_nCG]
							window_site_summary$seg_count = as.integer(levels(window_site_summary$seg_count))[window_site_summary$seg_count]
							window_site_summary$width_quantile = as.integer(levels(window_site_summary$width_quantile))[window_site_summary$width_quantile]
							first_window = FALSE
						} else {
							window_site_summary=rbind.data.frame(window_site_summary, setNames(c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE), mean(window_sites.gr$m_to_u_count, na.rm=TRUE), var(window_sites.gr$m_to_u_count, na.rm=TRUE), mean(window_sites.gr$u_to_m_count, na.rm=TRUE), var(window_sites.gr$u_to_m_count, na.rm=TRUE), mean(window_sites.gr$CG_site_count, na.rm=TRUE), var(window_sites.gr$CG_site_count, na.rm=TRUE), nrow(window_sites), width_quantile), names(window_site_summary)))
						}
				#		window_summary=rbind.data.frame(window_summary, setNames(c("Random UMR CG sites", window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE)), names(window_summary)))
				#		colnames(window_summary) = c("category","position","mean_vmCG","var_vmCG")

					}
	
					#window_summary$category = as.character(window_summary$category)
					window_site_summary$position = as.integer(window_site_summary$position)
					window_site_summary$mean_vmCG = as.numeric(window_site_summary$mean_vmCG)
					window_site_summary$var_vmCG = as.numeric(window_site_summary$var_vmCG)
					window_site_summary$mean_vMU = as.numeric(window_site_summary$mean_vMU)
					window_site_summary$var_vMU = as.numeric(window_site_summary$var_vMU)
					window_site_summary$mean_vUM = as.numeric(window_site_summary$mean_vUM)
					window_site_summary$var_vUM = as.numeric(window_site_summary$var_vUM)
					window_site_summary$mean_nCG = as.numeric(window_site_summary$mean_nCG)
					window_site_summary$var_nCG = as.numeric(window_site_summary$var_nCG)
					window_site_summary$seg_count = as.integer(window_site_summary$seg_count)
					window_site_summary$width_quantile = as.integer(window_site_summary$width_quantile)

					# Plot mean and variance of vmCG in windows centred at segment boundaries
					#pdf(paste0(project_id,"_",meth_context,"_",focus_mStatus,"_",adjacent_mStatus,"_variable_CG_sites_windows_",CHG_window_size,".pdf"))
					#print(ggplot() +geom_pointrange(data=window_site_summary, aes(x=position, y=mean_vmCG, ymin=mean_vmCG-var_vmCG/2, ymax=mean_vmCG+var_vmCG/2, colour=category)))
					#dev.off()
				}
			}
		}	
	} # for each quantile
	
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
		
		

	# Reshape data into one row per quantity plotted per window (count data)
	paneltext=array()
	paneltext[1]=paste0("Mean increase in number of mCG sites per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Change"=rep("U->M sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$mean_vUM, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("M->U sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("Net changes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"= window_site_summary$mean_vUM-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$change_count=as.numeric(levels(window_site_detail$change_count)[window_site_detail$change_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
		
	# Plot counts of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_count_by_segment_boundary_by_width_windows_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}	
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	# Plot counts of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_count_by_segment_boundary_by_width_windows_mirror_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	
	
	# Reshape data into one row per quantity plotted per window (rate data)
	paneltext=array()
	paneltext[1]=paste0("Mean increase in proportion of mCG sites per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Change"=rep("U->M sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$mean_vUM, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("M->U sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("Net changes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"= window_site_summary$mean_vUM-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$change_count=as.numeric(levels(window_site_detail$change_count)[window_site_detail$change_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
	
	# Plot rates of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_prop_by_segment_boundary_by_width_windows_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	# Plot rates of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_prop_by_segment boundary_by_width_windows_mirror_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	
	# Remake the boundary plots by dividing GBM, TEM and GBM-like segments into 3 quantiles each (short, medium, long), and finding boundaries 
	# This time, when we calculate rates we will use changes per M/U site rather than changes per site as a whole

	window_size = 100 # window length in nt
	step_size=10 # No. nt to slide window per step
	window_range = 250 # No steps to take each side of central site
	
	window_site_summary=data.frame()
	first_window = TRUE

	no_width_quantiles = 3
	for (width_quantile in 1:no_width_quantiles) {
		for (focus_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
			min_width = quantile(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_mStatus,]$width,(width_quantile-1)/no_width_quantiles)
			max_width = quantile(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_mStatus,]$width,(width_quantile)/no_width_quantiles)
			
			for (adjacent_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
	
				if (focus_mStatus != adjacent_mStatus) {
					#focus_mStatus = "UMR"
					#adjacent_mStatus = "TEM"
	
					# Create set of central sites to be the segments on the right hand side of the requested boundaries
					central_sites = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[(segmentation_model.gr$segment.mStatus==adjacent_mStatus) & (segmentation_model.gr@ranges@width>min_width) & (segmentation_model.gr@ranges@width<=max_width)], maxgap=1L), end(first) == start(second) - 1L)@second)
				
					# Create corresponding set of left hand segments for each boundary
					left_hand_segments = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[(segmentation_model.gr$segment.mStatus==adjacent_mStatus) & (segmentation_model.gr@ranges@width>min_width) & (segmentation_model.gr@ranges@width<=max_width)], maxgap=1L), end(first) == start(second) - 1L)@first)

					# We want to avoid accumulating data from outside of the segments under consideration.  Accordingly, create genomicranges corresponding to the adjoining halfs of each of the two adjoining segments.  We only use adjoining halfs of the segments to avoid edge effects from the other end of the segments
					relevant_sites_end = central_sites$end-round((central_sites$end-central_sites$start)/2)
					relevant_sites_start = left_hand_segments$end-round((left_hand_segments$end-left_hand_segments$start)/2)

					cat(paste0(nrow(central_sites)," boundaries between ",focus_mStatus," and ",adjacent_mStatus," segments.\n"))
					#	levels(central_sites@seqnames)=toupper(levels(central_sites@seqnames))
					#	central_sites@seqinfo@seqnames=levels(central_sites@seqnames)


					for (window_no in -window_range:window_range) {
						# for upstream and downstream:
						# make a genomic ranges for the window
						# calculate the No. CG sites and No. variable sites for the window
						# store the resulting values 
		
						window_sites=data.frame(cbind("seqnames"=as.character(central_sites$seqnames), "start"=central_sites$start + (window_no*step_size) - round((window_size/2) - 1.25), "end"=central_sites$start + (window_no*step_size) +round((window_size/2) + 0.25), "width"=window_size, "strand"=as.character(central_sites$strand)))

						window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
						window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
						window_sites$start = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$start<relevant_sites_start, relevant_sites_start, window_sites$start)))
						window_sites$end = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$end>relevant_sites_end, relevant_sites_end, window_sites$end)))

						# Remove any windows which did not overlap relvant sites at all
						window_sites = window_sites[!is.na(window_sites$start),]
					
						# Remove any windows which overlap ends of chromosomes
						if(window_no<0) {
							window_sites = window_sites[window_sites$end>=1,]
							window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
						} else {
							window_sites = window_sites[window_sites$start>=1,]
							window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
						}
		
						window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end", seqnames.field = "seqnames")

						window_sites.gr$ID=row.names(as.data.frame(window_sites.gr))
						#CG_segment_olaps = findOverlaps(CG_site_ranges, window_sites.gr)
						window_sites.gr$CG_site_count=table(window_sites.gr[subjectHits(findOverlaps(CG_site_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$CG_site_count[is.na(window_sites.gr$CG_site_count)]=0
						#variant_segment_olaps = findOverlaps(variant_call_ranges, window_sites.gr)
						window_sites.gr$variant_count=table(window_sites.gr[subjectHits(findOverlaps(variant_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$variant_count[is.na(window_sites.gr$variant_count)]=0
						#all_M_segment_olaps = findOverlaps(all_M_ranges, window_sites.gr)
						window_sites.gr$all_M_count=table(window_sites.gr[subjectHits(findOverlaps(all_M_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$all_M_count[is.na(window_sites.gr$all_M_count)]=0
						#all_U_segment_olaps = findOverlaps(all_U_ranges, window_sites.gr)
						window_sites.gr$all_U_count=table(window_sites.gr[subjectHits(findOverlaps(all_U_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$all_U_count[is.na(window_sites.gr$all_U_count)]=0
						#m_to_u_segment_olaps = findOverlaps(m_to_u_call_ranges, window_sites.gr)
						window_sites.gr$m_to_u_count=table(window_sites.gr[subjectHits(findOverlaps(m_to_u_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$m_to_u_count[is.na(window_sites.gr$m_to_u_count)]=0
						#u_to_m_segment_olaps = findOverlaps(u_to_m_call_ranges, window_sites.gr)
						window_sites.gr$u_to_m_count=table(window_sites.gr[subjectHits(findOverlaps(u_to_m_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$u_to_m_count[is.na(window_sites.gr$u_to_m_count)]=0
		
						# ensure all window sites overlap with one or more CHG sites
						#window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

						#window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
						cat(paste0("Window: ",window_no," mean variants: ",mean(window_sites.gr$variant_count, na.rm=TRUE), "variance: ",var(window_sites.gr$variant_count, na.rm=TRUE),"\n"))
			
						if (first_window) {
							window_site_summary=rbind.data.frame(window_site_summary, c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE), mean(window_sites.gr$m_to_u_count, na.rm=TRUE), var(window_sites.gr$m_to_u_count, na.rm=TRUE), mean(window_sites.gr$u_to_m_count, na.rm=TRUE), var(window_sites.gr$u_to_m_count, na.rm=TRUE), mean(window_sites.gr$CG_site_count, na.rm=TRUE), var(window_sites.gr$CG_site_count, na.rm=TRUE), mean(window_sites.gr$all_M_count, na.rm=TRUE), var(window_sites.gr$all_M_count, na.rm=TRUE), mean(window_sites.gr$all_U_count, na.rm=TRUE), var(window_sites.gr$all_U_count, na.rm=TRUE), nrow(window_sites), width_quantile))
							names(window_site_summary) = c("category","position","mean_vmCG","var_vmCG","mean_vMU","var_vMU","mean_vUM","var_vUM","mean_nCG","var_nCG","mean_all_M","var_all_M","mean_all_U","var_all_U","seg_count","width_quantile")
							window_site_summary$category = as.character(window_site_summary$category)
							window_site_summary$position = as.integer(levels(window_site_summary$position))[window_site_summary$position]
							window_site_summary$mean_vmCG = as.numeric(levels(window_site_summary$mean_vmCG))[window_site_summary$mean_vmCG]
							window_site_summary$var_vmCG = as.numeric(levels(window_site_summary$var_vmCG))[window_site_summary$var_vmCG]
							window_site_summary$mean_vMU = as.numeric(levels(window_site_summary$mean_vMU))[window_site_summary$mean_vMU]
							window_site_summary$var_vMU = as.numeric(levels(window_site_summary$var_vMU))[window_site_summary$var_vMU]
							window_site_summary$mean_vUM = as.numeric(levels(window_site_summary$mean_vUM))[window_site_summary$mean_vUM]
							window_site_summary$var_vUM = as.numeric(levels(window_site_summary$var_vUM))[window_site_summary$var_vUM]
							window_site_summary$mean_nCG = as.numeric(levels(window_site_summary$mean_nCG))[window_site_summary$mean_nCG]
							window_site_summary$var_nCG = as.numeric(levels(window_site_summary$var_nCG))[window_site_summary$var_nCG]
							window_site_summary$mean_all_M = as.numeric(levels(window_site_summary$mean_all_M))[window_site_summary$mean_all_M]
							window_site_summary$var_all_M = as.numeric(levels(window_site_summary$var_all_M))[window_site_summary$var_all_M]
							window_site_summary$mean_all_U = as.numeric(levels(window_site_summary$mean_all_U))[window_site_summary$mean_all_U]
							window_site_summary$var_all_U = as.numeric(levels(window_site_summary$var_all_U))[window_site_summary$var_all_U]
							window_site_summary$seg_count = as.integer(levels(window_site_summary$seg_count))[window_site_summary$seg_count]
							window_site_summary$width_quantile = as.integer(levels(window_site_summary$width_quantile))[window_site_summary$width_quantile]
							first_window = FALSE
						} else {
							window_site_summary=rbind.data.frame(window_site_summary, setNames(c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE), mean(window_sites.gr$m_to_u_count, na.rm=TRUE), var(window_sites.gr$m_to_u_count, na.rm=TRUE), mean(window_sites.gr$u_to_m_count, na.rm=TRUE), var(window_sites.gr$u_to_m_count, na.rm=TRUE), mean(window_sites.gr$CG_site_count, na.rm=TRUE), var(window_sites.gr$CG_site_count, na.rm=TRUE), mean(window_sites.gr$all_M_count, na.rm=TRUE), var(window_sites.gr$all_M_count, na.rm=TRUE), mean(window_sites.gr$all_U_count, na.rm=TRUE), var(window_sites.gr$all_U_count, na.rm=TRUE), nrow(window_sites), width_quantile), names(window_site_summary)))
						}
				#		window_summary=rbind.data.frame(window_summary, setNames(c("Random UMR CG sites", window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE)), names(window_summary)))
				#		colnames(window_summary) = c("category","position","mean_vmCG","var_vmCG")

					}
	
					#window_summary$category = as.character(window_summary$category)
					window_site_summary$position = as.integer(window_site_summary$position)
					window_site_summary$mean_vmCG = as.numeric(window_site_summary$mean_vmCG)
					window_site_summary$var_vmCG = as.numeric(window_site_summary$var_vmCG)
					window_site_summary$mean_vMU = as.numeric(window_site_summary$mean_vMU)
					window_site_summary$var_vMU = as.numeric(window_site_summary$var_vMU)
					window_site_summary$mean_vUM = as.numeric(window_site_summary$mean_vUM)
					window_site_summary$var_vUM = as.numeric(window_site_summary$var_vUM)
					window_site_summary$mean_nCG = as.numeric(window_site_summary$mean_nCG)
					window_site_summary$var_nCG = as.numeric(window_site_summary$var_nCG)
					window_site_summary$mean_all_M = as.numeric(window_site_summary$mean_all_M)
					window_site_summary$var_all_M = as.numeric(window_site_summary$var_all_M)
					window_site_summary$mean_all_U = as.numeric(window_site_summary$mean_all_U)
					window_site_summary$var_all_U = as.numeric(window_site_summary$var_all_U)
					window_site_summary$seg_count = as.integer(window_site_summary$seg_count)
					window_site_summary$width_quantile = as.integer(window_site_summary$width_quantile)

					# Plot mean and variance of vmCG in windows centred at segment boundaries
					#pdf(paste0(project_id,"_",meth_context,"_",focus_mStatus,"_",adjacent_mStatus,"_variable_CG_sites_windows_",CHG_window_size,".pdf"))
					#print(ggplot() +geom_pointrange(data=window_site_summary, aes(x=position, y=mean_vmCG, ymin=mean_vmCG-var_vmCG/2, ymax=mean_vmCG+var_vmCG/2, colour=category)))
					#dev.off()
				}
			}
		}	
	} # for each quantile
	
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
		
		

	# Reshape data into one row per quantity plotted per window (count data)
	paneltext=array()
	paneltext[1]=paste0("Mean increase in number of mCG sites per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Change"=rep("U->M sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$mean_vUM, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("M->U sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("Net changes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"= window_site_summary$mean_vUM-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_nCG, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$change_count=as.numeric(levels(window_site_detail$change_count)[window_site_detail$change_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
		
	# Plot counts of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_count2_by_segment_boundary_by_width_windows_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}	
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	# Plot counts of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_count2_by_segment_boundary_by_width_windows_mirror_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=change_count, colour=Change)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()
	
	
	
	# Reshape data into one row per quantity plotted per window (rate data - U->M per U site and M->U per M site)
	paneltext=array()
	paneltext[1]=paste0("Mean proportion of CG sites becoming methylated per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Change"=rep("U->M sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$mean_vUM, "site_count"=window_site_summary$mean_all_U+window_site_summary$mean_vUM, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("M->U sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_all_M+window_site_summary$mean_vMU, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("Net changes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"= window_site_summary$mean_vUM-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_all_M+window_site_summary$mean_all_U+window_site_summary$mean_vUM+window_site_summary$mean_vMU, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$change_count=as.numeric(levels(window_site_detail$change_count)[window_site_detail$change_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
	
	# Plot rates of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_prop2_by_segment_boundary_by_width_windows_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	# Plot rates of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_prop2_by_segment boundary_by_width_windows_mirror_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	

	# Reshape data into one row per quantity plotted per window (rate data - U->M per U site and M->U per M site), but without seg_count data
	paneltext=array()
	paneltext[1]=paste0("Mean proportion of CG sites becoming methylated per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Change"=rep("U->M sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=window_site_summary$mean_vUM, "site_count"=window_site_summary$mean_all_U+window_site_summary$mean_vUM, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("M->U sites", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"=-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_all_M+window_site_summary$mean_vMU, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Change"=rep("Net changes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "change_count"= window_site_summary$mean_vUM-window_site_summary$mean_vMU, "site_count"=window_site_summary$mean_all_M+window_site_summary$mean_all_U+window_site_summary$mean_vUM+window_site_summary$mean_vMU, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$change_count=as.numeric(levels(window_site_detail$change_count)[window_site_detail$change_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
	
	# Plot rates of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_prop4_by_segment_boundary_by_width_windows_",window_size,".pdf"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nby GBM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nby TEM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nby GBM-like segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()

	# Plot rates of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_mCG_change_prop4_by_segment boundary_by_width_windows_mirror_",window_size,".pdf"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nby GBM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nby TEM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")),], aes(x=position, y=ifelse(site_count==0,0,change_count/site_count), colour=Change)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nby GBM-like segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()

	
	
	
	#### Comparison with nucleosome positions
	### Read in the nucleosome positions as a dump of genomicRanges which was previously prepared in a VM (Windows cannot import BigWig files)
	nucleosome_positions.gr = makeGRangesFromDataFrame(df = read.table(reference_nucleosomes), start.field = "V2", end.field = "V3", seqnames.field = "V1")

	
	# Plot nuclosome coverage +- 1kb relative to variable sites
	
	# Remake the boundary plots by calculating the number of times a relevant window overlaps with a nucleosome
	# This time, when we calculate rates we will use No. nucleosome overlaps / No. hemi-segments in boundary

	levels(nucleosome_positions.gr@seqnames)=toupper(levels(nucleosome_positions.gr@seqnames))
	nucleosome_positions.gr@seqinfo@seqnames=levels(nucleosome_positions.gr@seqnames)

	
	window_size = 1 # window length in nt
	step_size=window_size # No. nt to slide window per step
	window_range = 250 # No steps to take each side of central site
	
	window_site_summary=data.frame()
	first_window = TRUE

	no_width_quantiles = 3
	for (width_quantile in 1:no_width_quantiles) {
		for (focus_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
			min_width = quantile(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_mStatus,]$width,(width_quantile-1)/no_width_quantiles)
			max_width = quantile(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_mStatus,]$width,(width_quantile)/no_width_quantiles)
			
			for (adjacent_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
	
				if (focus_mStatus != adjacent_mStatus) {
					#focus_mStatus = "UMR"
					#adjacent_mStatus = "TEM"
	
					# Create set of central sites to be the segments on the right hand side of the requested boundaries
					central_sites = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[(segmentation_model.gr$segment.mStatus==adjacent_mStatus) & (segmentation_model.gr@ranges@width>min_width) & (segmentation_model.gr@ranges@width<=max_width)], maxgap=1L), end(first) == start(second) - 1L)@second)
				
					# Create corresponding set of left hand segments for each boundary
					left_hand_segments = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[(segmentation_model.gr$segment.mStatus==adjacent_mStatus) & (segmentation_model.gr@ranges@width>min_width) & (segmentation_model.gr@ranges@width<=max_width)], maxgap=1L), end(first) == start(second) - 1L)@first)

					# We want to avoid accumulating data from outside of the segments under consideration.  Accordingly, create genomicranges corresponding to the adjoining halfs of each of the two adjoining segments.  We only use adjoining halfs of the segments to avoid edge effects from the other end of the segments
					relevant_sites_end = central_sites$end-round((central_sites$end-central_sites$start)/2)
					relevant_sites_start = left_hand_segments$end-round((left_hand_segments$end-left_hand_segments$start)/2)

					cat(paste0(nrow(central_sites)," boundaries between ",focus_mStatus," and ",adjacent_mStatus," segments.\n"))
					#	levels(central_sites@seqnames)=toupper(levels(central_sites@seqnames))
					#	central_sites@seqinfo@seqnames=levels(central_sites@seqnames)


					for (window_no in -window_range:window_range) {
						# for upstream and downstream:
						# make a genomic ranges for the window
						# calculate the No. times the window is present in a relevant segment boundary, and No. times the window overlaps a nucleosome
						# store the resulting values 
		
						window_sites=data.frame(cbind("seqnames"=as.character(central_sites$seqnames), "start"=central_sites$start + (window_no*step_size) - round((window_size/2) - 1.25), "end"=central_sites$start + (window_no*step_size) +round((window_size/2) + 0.25), "width"=window_size, "strand"=as.character(central_sites$strand)))

						window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
						window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
						window_sites$start = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$start<relevant_sites_start, relevant_sites_start, window_sites$start)))
						window_sites$end = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$end>relevant_sites_end, relevant_sites_end, window_sites$end)))

						# Remove any windows which did not overlap relvant sites at all
						window_sites = window_sites[!is.na(window_sites$start),]
					
						# Remove any windows which overlap ends of chromosomes
						if(window_no<0) {
							window_sites = window_sites[window_sites$end>=1,]
							window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
						} else {
							window_sites = window_sites[window_sites$start>=1,]
							window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
						}
		
						window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end", seqnames.field = "seqnames")

						window_sites.gr$ID=row.names(as.data.frame(window_sites.gr))
						#window_sites.gr$CG_site_count=table(window_sites.gr[subjectHits(findOverlaps(CG_site_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$CG_site_count[is.na(window_sites.gr$CG_site_count)]=0
						#window_sites.gr$variant_count=table(window_sites.gr[subjectHits(findOverlaps(variant_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$variant_count[is.na(window_sites.gr$variant_count)]=0
						window_sites.gr$nucleosome_count=table(window_sites.gr[subjectHits(findOverlaps(nucleosome_positions.gr, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$nucleosome_count[is.na(window_sites.gr$nucleosome_count)]=0
						#window_sites.gr$all_M_count=table(window_sites.gr[subjectHits(findOverlaps(all_M_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$all_M_count[is.na(window_sites.gr$all_M_count)]=0
						#window_sites.gr$all_U_count=table(window_sites.gr[subjectHits(findOverlaps(all_U_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$all_U_count[is.na(window_sites.gr$all_U_count)]=0
						#window_sites.gr$m_to_u_count=table(window_sites.gr[subjectHits(findOverlaps(m_to_u_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$m_to_u_count[is.na(window_sites.gr$m_to_u_count)]=0
						#window_sites.gr$u_to_m_count=table(window_sites.gr[subjectHits(findOverlaps(u_to_m_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$u_to_m_count[is.na(window_sites.gr$u_to_m_count)]=0
		
						# ensure all window sites overlap with one or more CHG sites
						#window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

						#window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
						cat(paste0("Window: ",window_no," nucleosome proportion: ",mean(window_sites.gr$nucleosome_count, na.rm=TRUE), "variance: ",var(window_sites.gr$nucleosome_count, na.rm=TRUE),"\n"))
			
						if (first_window) {
							#window_site_summary=rbind.data.frame(window_site_summary, c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$nucleosome_count, na.rm=TRUE), var(window_sites.gr$nucleosome_count, na.rm=TRUE), mean(window_sites.gr$m_to_u_count, na.rm=TRUE), var(window_sites.gr$m_to_u_count, na.rm=TRUE), mean(window_sites.gr$u_to_m_count, na.rm=TRUE), var(window_sites.gr$u_to_m_count, na.rm=TRUE), mean(window_sites.gr$CG_site_count, na.rm=TRUE), var(window_sites.gr$CG_site_count, na.rm=TRUE), mean(window_sites.gr$all_M_count, na.rm=TRUE), var(window_sites.gr$all_M_count, na.rm=TRUE), mean(window_sites.gr$all_U_count, na.rm=TRUE), var(window_sites.gr$all_U_count, na.rm=TRUE), nrow(window_sites), width_quantile))
							#names(window_site_summary) = c("category","position","mean_vmCG","var_vmCG","mean_vMU","var_vMU","mean_vUM","var_vUM","mean_nCG","var_nCG","mean_all_M","var_all_M","mean_all_U","var_all_U","seg_count","width_quantile")
							window_site_summary=rbind.data.frame(window_site_summary, c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$nucleosome_count, na.rm=TRUE), var(window_sites.gr$nucleosome_count, na.rm=TRUE), sum(window_sites.gr$nucleosome_count, na.rm=TRUE), nrow(window_sites), width_quantile))
							names(window_site_summary) = c("category","position","mean_nucl","var_nucl","nucl_count","seg_count","width_quantile")
							window_site_summary$category = as.character(window_site_summary$category)
							window_site_summary$position = as.integer(levels(window_site_summary$position))[window_site_summary$position]
							window_site_summary$mean_nucl = as.numeric(levels(window_site_summary$mean_nucl))[window_site_summary$mean_nucl]
							window_site_summary$var_nucl = as.numeric(levels(window_site_summary$var_nucl))[window_site_summary$var_nucl]
							#window_site_summary$mean_vmCG = as.numeric(levels(window_site_summary$mean_vmCG))[window_site_summary$mean_vmCG]
							#window_site_summary$var_vmCG = as.numeric(levels(window_site_summary$var_vmCG))[window_site_summary$var_vmCG]
							#window_site_summary$mean_vMU = as.numeric(levels(window_site_summary$mean_vMU))[window_site_summary$mean_vMU]
							#window_site_summary$var_vMU = as.numeric(levels(window_site_summary$var_vMU))[window_site_summary$var_vMU]
							#window_site_summary$mean_vUM = as.numeric(levels(window_site_summary$mean_vUM))[window_site_summary$mean_vUM]
							#window_site_summary$var_vUM = as.numeric(levels(window_site_summary$var_vUM))[window_site_summary$var_vUM]
							#window_site_summary$mean_nCG = as.numeric(levels(window_site_summary$mean_nCG))[window_site_summary$mean_nCG]
							#window_site_summary$var_nCG = as.numeric(levels(window_site_summary$var_nCG))[window_site_summary$var_nCG]
							#window_site_summary$mean_all_M = as.numeric(levels(window_site_summary$mean_all_M))[window_site_summary$mean_all_M]
							#window_site_summary$var_all_M = as.numeric(levels(window_site_summary$var_all_M))[window_site_summary$var_all_M]
							#window_site_summary$mean_all_U = as.numeric(levels(window_site_summary$mean_all_U))[window_site_summary$mean_all_U]
							#window_site_summary$var_all_U = as.numeric(levels(window_site_summary$var_all_U))[window_site_summary$var_all_U]
							window_site_summary$nucl_count = as.integer(levels(window_site_summary$nucl_count))[window_site_summary$nucl_count]
							window_site_summary$seg_count = as.integer(levels(window_site_summary$seg_count))[window_site_summary$seg_count]
							window_site_summary$width_quantile = as.integer(levels(window_site_summary$width_quantile))[window_site_summary$width_quantile]
							first_window = FALSE
						} else {
							window_site_summary=rbind.data.frame(window_site_summary, setNames(c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$nucleosome_count, na.rm=TRUE), var(window_sites.gr$nucleosome_count, na.rm=TRUE), sum(window_sites.gr$nucleosome_count, na.rm=TRUE), nrow(window_sites), width_quantile), names(window_site_summary)))
						}
				#		window_summary=rbind.data.frame(window_summary, setNames(c("Random UMR CG sites", window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE)), names(window_summary)))
				#		colnames(window_summary) = c("category","position","mean_vmCG","var_vmCG")

					}
	
					#window_summary$category = as.character(window_summary$category)
					window_site_summary$position = as.integer(window_site_summary$position)
					window_site_summary$mean_nucl = as.numeric(window_site_summary$mean_nucl)
					window_site_summary$var_nucl = as.numeric(window_site_summary$var_nucl)
					#window_site_summary$mean_vmCG = as.numeric(window_site_summary$mean_vmCG)
					#window_site_summary$var_vmCG = as.numeric(window_site_summary$var_vmCG)
					#window_site_summary$mean_vMU = as.numeric(window_site_summary$mean_vMU)
					#window_site_summary$var_vMU = as.numeric(window_site_summary$var_vMU)
					#window_site_summary$mean_vUM = as.numeric(window_site_summary$mean_vUM)
					#window_site_summary$var_vUM = as.numeric(window_site_summary$var_vUM)
					#window_site_summary$mean_nCG = as.numeric(window_site_summary$mean_nCG)
					#window_site_summary$var_nCG = as.numeric(window_site_summary$var_nCG)
					#window_site_summary$mean_all_M = as.numeric(window_site_summary$mean_all_M)
					#window_site_summary$var_all_M = as.numeric(window_site_summary$var_all_M)
					#window_site_summary$mean_all_U = as.numeric(window_site_summary$mean_all_U)
					#window_site_summary$var_all_U = as.numeric(window_site_summary$var_all_U)
					window_site_summary$nucl_count = as.integer(window_site_summary$nucl_count)
					window_site_summary$seg_count = as.integer(window_site_summary$seg_count)
					window_site_summary$width_quantile = as.integer(window_site_summary$width_quantile)

					# Plot mean and variance of vmCG in windows centred at segment boundaries
					#pdf(paste0(project_id,"_",meth_context,"_",focus_mStatus,"_",adjacent_mStatus,"_variable_CG_sites_windows_",CHG_window_size,".pdf"))
					#print(ggplot() +geom_pointrange(data=window_site_summary, aes(x=position, y=mean_vmCG, ymin=mean_vmCG-var_vmCG/2, ymax=mean_vmCG+var_vmCG/2, colour=category)))
					#dev.off()
				}
			}
		}	
	} # for each quantile
	
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
		
		

	# Reshape data into one row per quantity plotted per window (count data)
	paneltext=array()
	paneltext[1]=paste0("Proportion of segments overlapping nucleosomes per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Feature"=rep("Nuclosomes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$nucl_count, "site_count"=window_site_summary$seg_count, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Feature"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$feature_count=as.numeric(levels(window_site_detail$feature_count)[window_site_detail$feature_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
		
	# Plot counts of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_count_by_segment_boundary_by_width_windows_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}	
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	# Plot counts of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_count_by_segment_boundary_by_width_windows_mirror_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()
	
	
	
	# Reshape data into one row per quantity plotted per window (rate data - nucleosomes per segment)
	window_site_detail = rbind.data.frame(cbind("Feature"=rep("Nuclosomes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$nucl_count/window_site_summary$seg_count, "site_count"=window_site_summary$seg_count, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Feature"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$feature_count=as.numeric(levels(window_site_detail$feature_count)[window_site_detail$feature_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
		
	# Plot counts of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_prop_by_segment_boundary_by_width_windows_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}	
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	# Plot counts of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_prop_by_segment_boundary_by_width_windows_mirror_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()
	
	

	# Reshape data into one row per quantity plotted per window (rate data - nucleosomes per segment), but without seg_count data
	paneltext=array()
	paneltext[1]=paste0("Proportion of segments overlapping nucleosomes per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Feature"=rep("Nuclosomes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$nucl_count, "site_count"=window_site_summary$seg_count, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$feature_count=as.numeric(levels(window_site_detail$feature_count)[window_site_detail$feature_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])

	# Plot rates of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_prop4_by_segment_boundary_by_width_windows_",window_size,".pdf"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nby GBM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nby TEM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nby GBM-like segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()

	# Plot rates of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_prop4_by_segment boundary_by_width_windows_mirror_",window_size,".pdf"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nby GBM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nby TEM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nby GBM-like segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()

	
	### Correlate mCG-variable sites with histone modification data
	
	### Read in the H3K* data as a dump of genomicRanges which was previously prepared in a VM (Windows cannot import BigWig files)
	H3K9me2.gr = readRDS(reference_H3K9me2)	# This data is normalised and smoothed over contiguous 100bp windows
	H3K4me1.gr = readRDS(reference_H3K4me1) # This data is not normalised and is over 50nt windows
	H3K27me3.gr = readRDS(reference_H3K27me3) # This data is not normalised and is over 50nt windows
	H2AW = read.table(reference_H2AW, header=FALSE, sep="\t") # This data is normalised over contiguous 50bp windows
	H2AW = H2AW[substr(H2AW$V1,1,3) == "chr",] # drop chloroplast and mitochondria
	H2AW$V1 = mixedCaseChr(H2AW$V1)
	H2AW.gr = makeGRangesFromDataFrame(df=H2AW, start.field = "V4", end.field = "V5",seqnames.field = "V1")
	H2AW.gr@elementMetadata@listData$score = H2AW$V6


	# Compare histone methylation at the 50bp/100bp windows overlapping each M->U site and U->M site, divided by segment type
	
	#for (focus_segment_type in c("GBM", "GBM-like", "TEM", "UMR")) {
	for (focus_segment_type in c("GBM", "GBM-like", "TEM")) {  # leave out UMR as not enough methylation
	#focus_segment_type = "TEM"
	window_summary=data.frame()

	first_window = TRUE

	direction_strings = c("M->U", "U->M")
	
	for (direction in direction_strings) {
	
	for (num_lines_varying in 0:4) {
		
		if (num_lines_varying == 0) {
			#sample_size=nrow(central_sites)
			sample_size = 1000
			
			# Random all_M or all_U sites in TEM segments:
			if (direction == "M->U") {
				central_sites.gr=makeGRangesFromDataFrame(df=as.data.frame(all_M_ranges)[as.data.frame(all_M_segment_olaps)[as.data.frame(all_M_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_segment_type,]),]$queryHits,][sample(nrow(as.data.frame(all_M_ranges)[as.data.frame(all_M_segment_olaps)[as.data.frame(all_M_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_segment_type,]),]$queryHits,]), sample_size, replace=TRUE),])
			} else {
				central_sites.gr=makeGRangesFromDataFrame(df=as.data.frame(all_U_ranges)[as.data.frame(all_U_segment_olaps)[as.data.frame(all_U_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_segment_type,]),]$queryHits,][sample(nrow(as.data.frame(all_U_ranges)[as.data.frame(all_U_segment_olaps)[as.data.frame(all_U_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_segment_type,]),]$queryHits,]), sample_size, replace=TRUE),])
			}
		} else {

			# Set central_sites to be the CG sites where mCG varies in num_lines_varying offspring, and which occur in TEM segments
			if (direction == "M->U") {
				central_sites.gr=makeGRangesFromDataFrame(df=cbind.data.frame(as.data.frame(m_to_u_call_ranges),"num_lines"=num_lines_m_to_u_from_parentals[num_lines_m_to_u_from_parentals>0])[(as.data.frame(m_to_u_segment_olaps)[as.data.frame(m_to_u_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_segment_type,]),]$queryHits),][cbind.data.frame(as.data.frame(m_to_u_call_ranges),"num_lines"=num_lines_m_to_u_from_parentals[num_lines_m_to_u_from_parentals>0])[(as.data.frame(m_to_u_segment_olaps)[as.data.frame(m_to_u_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_segment_type,]),]$queryHits),]$num_lines==num_lines_varying,])
			} else {
				central_sites.gr=makeGRangesFromDataFrame(df=cbind.data.frame(as.data.frame(u_to_m_call_ranges),"num_lines"=num_lines_u_to_m_from_parentals[num_lines_u_to_m_from_parentals>0])[(as.data.frame(u_to_m_segment_olaps)[as.data.frame(u_to_m_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_segment_type,]),]$queryHits),][cbind.data.frame(as.data.frame(u_to_m_call_ranges),"num_lines"=num_lines_u_to_m_from_parentals[num_lines_u_to_m_from_parentals>0])[(as.data.frame(u_to_m_segment_olaps)[as.data.frame(u_to_m_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_segment_type,]),]$queryHits),]$num_lines==num_lines_varying,])
			}
		#central_sites.gr=makeGRangesFromDataFrame(df=as.data.frame(variant_call_ranges)[as.data.frame(variant_segment_olaps)[as.data.frame(variant_segment_olaps)$subjectHits %in% rownames(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM",]),]$queryHits,])
		}

		central_sites = as.data.frame(setdiff(central_sites.gr, DMRs.gr, ignore.strand=TRUE))
	
		CHG_window_size = 1 # window length in nt
		CHG_step_size=50 # No. nt to slide window per step
		CHG_window_range = 0 # No steps to take each side of central site
	
		for (CHG_window_no in -CHG_window_range:CHG_window_range) {
			# for upstream and downstream:
			# make a genomic ranges for the window
			# calculate the mean mCHG for the window
			# store the resulting mCHG values 
		
			window_sites=data.frame(cbind("seqnames"=mixedCaseChr(as.character(central_sites$seqnames)), "start"=central_sites$start + (CHG_window_no*CHG_step_size) - round((CHG_window_size/2) - 1.25), "end"=central_sites$start + (CHG_window_no*CHG_step_size) +round((CHG_window_size/2) + 0.25), "width"=CHG_window_size, "strand"=as.character(central_sites$strand)))
		
			window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
			window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
			# Remove any windows which overlap ends of chromosomes
			if(CHG_window_no<0) {
				window_sites = window_sites[window_sites$end>=1,]
				window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
			} else {
				window_sites = window_sites[window_sites$start>=1,]
				window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
			}

			# First do H3K9me2
			#browser()
			window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")

			# ensure all window sites overlap with one or more H3K9me2 windows
			window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(H3K9me2.gr, window_sites.gr)))]

			#window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
			# Read off the H3K9me2 values from the 100nt windows (used in the H3K9me2 file) overlapping the sites of interest
			window_meths = as.data.frame(H3K9me2.gr[unique(queryHits(findOverlaps(H3K9me2.gr, window_sites.gr)))])
		
			#cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), " variance: ",var(window_meths$pmeth, na.rm=TRUE)," prop_meth_0.2: ",sum(window_meths$pmeth>0.2)/nrow(window_meths)," prop_meth_0.4: ",sum(window_meths$pmeth>0.4)/nrow(window_meths)," prop_meth_0.6: ",sum(window_meths$pmeth>0.6)/nrow(window_meths)," prop_meth_0.8: ",sum(window_meths$pmeth>0.8)/nrow(window_meths),"\n"))
			cat(paste0("Window: ",CHG_window_no," mean H3K9me2: ",mean(window_meths$score, na.rm=TRUE), " variance: ",var(window_meths$score, na.rm=TRUE)," prop_meth_0.2: ",sum(window_meths$score>0.2)/nrow(window_meths)," prop_meth_0.4: ",sum(window_meths$score>0.4)/nrow(window_meths)," prop_meth_0.6: ",sum(window_meths$score>0.6)/nrow(window_meths)," prop_meth_0.8: ",sum(window_meths$score>0.8)/nrow(window_meths),"\n"))
			#print(ggplot(as.data.frame(window_meths)) + geom_histogram(aes(x=pmeth, y=(..count..)/sum(..count..)),bins=10) + xlim(0:1) + scale_y_continuous(labels = scales::percent, limits = c(0,0.3)))

			# Store up the results for plotting later
			if (first_window) {
				window_summary=rbind.data.frame(window_summary, cbind.data.frame(rep(direction, nrow(window_meths)), rep("H3K9me2", nrow(window_meths)), rep(num_lines_varying, nrow(window_meths)), window_meths$score))
				names(window_summary) = c("direction", "context", "num_lines", "m_nonCG")
				first_window = FALSE
			} else {
				window_summary=rbind.data.frame(window_summary, setNames(cbind.data.frame(rep(direction, nrow(window_meths)), rep("H3K9me2", nrow(window_meths)), rep(num_lines_varying, nrow(window_meths)), window_meths$score), names(window_summary)))
			}

			# Next do H3K4me1
			#browser()
			window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")

			# ensure all window sites overlap with one or more H3K4me1 windows
			window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(H3K4me1.gr, window_sites.gr)))]

			#window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
			# Read off the H3K4me1 values from the 100nt windows (used in the H3K4me1 file) overlapping the sites of interest
			window_meths = as.data.frame(H3K4me1.gr[unique(queryHits(findOverlaps(H3K4me1.gr, window_sites.gr)))])
			window_meths$score=log(window_meths$score)  # log for scaling
			
			#cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), " variance: ",var(window_meths$pmeth, na.rm=TRUE)," prop_meth_0.2: ",sum(window_meths$pmeth>0.2)/nrow(window_meths)," prop_meth_0.4: ",sum(window_meths$pmeth>0.4)/nrow(window_meths)," prop_meth_0.6: ",sum(window_meths$pmeth>0.6)/nrow(window_meths)," prop_meth_0.8: ",sum(window_meths$pmeth>0.8)/nrow(window_meths),"\n"))
			cat(paste0("Window: ",CHG_window_no," mean H3K4me1: ",mean(window_meths$score, na.rm=TRUE), " variance: ",var(window_meths$score, na.rm=TRUE)," prop_meth_0.2: ",sum(window_meths$score>0.2)/nrow(window_meths)," prop_meth_0.4: ",sum(window_meths$score>0.4)/nrow(window_meths)," prop_meth_0.6: ",sum(window_meths$score>0.6)/nrow(window_meths)," prop_meth_0.8: ",sum(window_meths$score>0.8)/nrow(window_meths),"\n"))
			#print(ggplot(as.data.frame(window_meths)) + geom_histogram(aes(x=pmeth, y=(..count..)/sum(..count..)),bins=10) + xlim(0:1) + scale_y_continuous(labels = scales::percent, limits = c(0,0.3)))

			# Store up the results for plotting later
			if (first_window) {
				window_summary=rbind.data.frame(window_summary, cbind.data.frame(rep(direction, nrow(window_meths)), rep("H3K4me1", nrow(window_meths)), rep(num_lines_varying, nrow(window_meths)), window_meths$score))
				names(window_summary) = c("direction", "context", "num_lines", "m_nonCG")
				first_window = FALSE
			} else {
				window_summary=rbind.data.frame(window_summary, setNames(cbind.data.frame(rep(direction, nrow(window_meths)), rep("H3K4me1", nrow(window_meths)), rep(num_lines_varying, nrow(window_meths)), window_meths$score), names(window_summary)))
			}

			# Next do H3K27me3
			#browser()
			window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")

			# ensure all window sites overlap with one or more H3K27me3 windows
			window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(H3K27me3.gr, window_sites.gr)))]

			#window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
			# Read off the H3K27me3 values from the 100nt windows (used in the H3K27me3 file) overlapping the sites of interest
			window_meths = as.data.frame(H3K27me3.gr[unique(queryHits(findOverlaps(H3K27me3.gr, window_sites.gr)))])
			window_meths$score=log(window_meths$score)  # log for scaling
			
			#cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), " variance: ",var(window_meths$pmeth, na.rm=TRUE)," prop_meth_0.2: ",sum(window_meths$pmeth>0.2)/nrow(window_meths)," prop_meth_0.4: ",sum(window_meths$pmeth>0.4)/nrow(window_meths)," prop_meth_0.6: ",sum(window_meths$pmeth>0.6)/nrow(window_meths)," prop_meth_0.8: ",sum(window_meths$pmeth>0.8)/nrow(window_meths),"\n"))
			cat(paste0("Window: ",CHG_window_no," mean H3K27me3: ",mean(window_meths$score, na.rm=TRUE), " variance: ",var(window_meths$score, na.rm=TRUE)," prop_meth_0.2: ",sum(window_meths$score>0.2)/nrow(window_meths)," prop_meth_0.4: ",sum(window_meths$score>0.4)/nrow(window_meths)," prop_meth_0.6: ",sum(window_meths$score>0.6)/nrow(window_meths)," prop_meth_0.8: ",sum(window_meths$score>0.8)/nrow(window_meths),"\n"))
			#print(ggplot(as.data.frame(window_meths)) + geom_histogram(aes(x=pmeth, y=(..count..)/sum(..count..)),bins=10) + xlim(0:1) + scale_y_continuous(labels = scales::percent, limits = c(0,0.3)))

			# Store up the results for plotting later
			if (first_window) {
				window_summary=rbind.data.frame(window_summary, cbind.data.frame(rep(direction, nrow(window_meths)), rep("H3K27me3", nrow(window_meths)), rep(num_lines_varying, nrow(window_meths)), window_meths$score))
				names(window_summary) = c("direction", "context", "num_lines", "m_nonCG")
				first_window = FALSE
			} else {
				window_summary=rbind.data.frame(window_summary, setNames(cbind.data.frame(rep(direction, nrow(window_meths)), rep("H3K27me3", nrow(window_meths)), rep(num_lines_varying, nrow(window_meths)), window_meths$score), names(window_summary)))
			}

			# Next do H2AW
			#browser()
			window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end",seqnames.field = "seqnames")

			# ensure all window sites overlap with one or more H2AW windows
			window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(H2AW.gr, window_sites.gr)))]

			#window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
			# Read off the H2AW values from the 50nt windows (used in the H2AW file) overlapping the sites of interest
			window_meths = as.data.frame(H2AW.gr[unique(queryHits(findOverlaps(H2AW.gr, window_sites.gr)))])
			
			#cat(paste0("Window: ",CHG_window_no," mean mCHG: ",mean(window_meths$pmeth, na.rm=TRUE), " variance: ",var(window_meths$pmeth, na.rm=TRUE)," prop_meth_0.2: ",sum(window_meths$pmeth>0.2)/nrow(window_meths)," prop_meth_0.4: ",sum(window_meths$pmeth>0.4)/nrow(window_meths)," prop_meth_0.6: ",sum(window_meths$pmeth>0.6)/nrow(window_meths)," prop_meth_0.8: ",sum(window_meths$pmeth>0.8)/nrow(window_meths),"\n"))
			cat(paste0("Window: ",CHG_window_no," mean H2AW: ",mean(window_meths$score, na.rm=TRUE), " variance: ",var(window_meths$score, na.rm=TRUE)," prop_meth_0.2: ",sum(window_meths$score>0.2)/nrow(window_meths)," prop_meth_0.4: ",sum(window_meths$score>0.4)/nrow(window_meths)," prop_meth_0.6: ",sum(window_meths$score>0.6)/nrow(window_meths)," prop_meth_0.8: ",sum(window_meths$score>0.8)/nrow(window_meths),"\n"))
			#print(ggplot(as.data.frame(window_meths)) + geom_histogram(aes(x=pmeth, y=(..count..)/sum(..count..)),bins=10) + xlim(0:1) + scale_y_continuous(labels = scales::percent, limits = c(0,0.3)))

			# Store up the results for plotting later
			if (first_window) {
				window_summary=rbind.data.frame(window_summary, cbind.data.frame(rep(direction, nrow(window_meths)), rep("H2AW", nrow(window_meths)), rep(num_lines_varying, nrow(window_meths)), window_meths$score))
				names(window_summary) = c("direction", "context", "num_lines", "m_nonCG")
				first_window = FALSE
			} else {
				window_summary=rbind.data.frame(window_summary, setNames(cbind.data.frame(rep(direction, nrow(window_meths)), rep("H2AW", nrow(window_meths)), rep(num_lines_varying, nrow(window_meths)), window_meths$score), names(window_summary)))
			}
			} # for CHG_window_no
	} # for num_lines_varying
	} # for direction
	
	
	pdf(paste0(project_id,"_",meth_context,"_",focus_segment_type,"_variable_CG_sites_noDMRs_H3K9me2_boxplots_",CHG_window_size,".pdf"))
	print(ggplot(window_summary, aes(y=m_nonCG, x = as.factor(num_lines))) + geom_boxplot(varwidth=TRUE, alpha = 0.2) + theme_minimal() + xlab("Number of independant offspring lines with site mCG change") + ylab("Histone methylation in 50/100nt window within which site mCG changes")+ coord_flip() +facet_grid(context ~ direction, scales="fixed"))
	dev.off()
	
	for (nonCG_context in c("H3K9me2","H3K4me1","H3K27me3","H2AW")) {
		for (direction in direction_strings) {
			for (num_lines_varying in 0:1) {
				cat(paste0(focus_segmnent_type," ",nonCG_context,", ",direction,", ",num_lines_varying," lines  "))
				cat("less: ",wilcox.test(window_summary[(window_summary$direction==direction) & (window_summary$context==nonCG_context) & (window_summary$num_lines<=num_lines_varying),]$m_nonCG, window_summary[(window_summary$direction==direction) & (window_summary$context==nonCG_context) & (window_summary$num_lines>num_lines_varying),]$m_nonCG, alternative="less")$p.value," more: ",wilcox.test(window_summary[(window_summary$direction==direction) & (window_summary$context==nonCG_context) & (window_summary$num_lines<=num_lines_varying),]$m_nonCG, window_summary[(window_summary$direction==direction) & (window_summary$context==nonCG_context) & (window_summary$num_lines>num_lines_varying),]$m_nonCG, alternative="greater")$p.value,"\n")
			}
		}
	}

	} # for focus_segment_type

	# Below results for Wilcoxon rank sum tests with continuity correction.
	# 'less' reports unadjusted p-value testing whether number of mCG changes is associated with higher histone modification level
	# 'more' reports unadjusted p-value testing whether number of mCG changes is associated with lower histone modification level
# GBM H3K9me2, M->U, 0 lines  less:  0.5683068  more:  0.4316946 
# GBM H3K9me2, M->U, 1 lines  less:  0.9038439  more:  0.09615635 
# GBM H3K9me2, U->M, 0 lines  less:  0.2372345  more:  0.7627674 
# GBM H3K9me2, U->M, 1 lines  less:  0.4408107  more:  0.5591906 
# GBM H3K4me1, M->U, 0 lines  less:  0.9999999  more:  7.522657e-08 
# GBM H3K4me1, M->U, 1 lines  less:  1  more:  3.771793e-39 
# GBM H3K4me1, U->M, 0 lines  less:  0.9688957  more:  0.0311047 
# GBM H3K4me1, U->M, 1 lines  less:  0.9993433  more:  0.0006566611 
# GBM H3K27me3, M->U, 0 lines  less:  0.01612021  more:  0.98388 
# GBM H3K27me3, M->U, 1 lines  less:  1.025732e-13  more:  1 
# GBM H3K27me3, U->M, 0 lines  less:  0.4025407  more:  0.5974626 
# GBM H3K27me3, U->M, 1 lines  less:  0.01202818  more:  0.987972 
# GBM H2AW, M->U, 0 lines  less:  0.6726855  more:  0.3273158 
# GBM H2AW, M->U, 1 lines  less:  0.2224516  more:  0.7775489 
# GBM H2AW, U->M, 0 lines  less:  0.1716095  more:  0.8283922 
# GBM H2AW, U->M, 1 lines  less:  5.625005e-07  more:  0.9999994 
# GBM-like H3K9me2, M->U, 0 lines  less:  0.9822253  more:  0.01777929 
# GBM-like H3K9me2, M->U, 1 lines  less:  0.3709837  more:  0.6290696 
# GBM-like H3K9me2, U->M, 0 lines  less:  0.178639  more:  0.8214038 
# GBM-like H3K9me2, U->M, 1 lines  less:  0.1405542  more:  0.8594892 
# GBM-like H3K4me1, M->U, 0 lines  less:  3.120205e-15  more:  1 
# GBM-like H3K4me1, M->U, 1 lines  less:  0.08640084  more:  0.9136409 
# GBM-like H3K4me1, U->M, 0 lines  less:  0.9999146  more:  8.544241e-05 
# GBM-like H3K4me1, U->M, 1 lines  less:  0.9999992  more:  7.616209e-07 
# GBM-like H3K27me3, M->U, 0 lines  less:  0.3335062  more:  0.666538 
# GBM-like H3K27me3, M->U, 1 lines  less:  0.01796117  more:  0.9820459 
# GBM-like H3K27me3, U->M, 0 lines  less:  6.393251e-05  more:  0.9999361 
# GBM-like H3K27me3, U->M, 1 lines  less:  0.01105749  more:  0.9889489 
# GBM-like H2AW, M->U, 0 lines  less:  1  more:  1.857123e-27 
# GBM-like H2AW, M->U, 1 lines  less:  1  more:  3.707524e-13 
# GBM-like H2AW, U->M, 0 lines  less:  7.586186e-08  more:  0.9999999 
# GBM-like H2AW, U->M, 1 lines  less:  3.242278e-09  more:  1 
# TEM H3K9me2, M->U, 0 lines  less:  1  more:  6.054511e-86 
# TEM H3K9me2, M->U, 1 lines  less:  1  more:  1.09205e-33 
# TEM H3K9me2, U->M, 0 lines  less:  6.183755e-17  more:  1 
# TEM H3K9me2, U->M, 1 lines  less:  1.134785e-14  more:  1 
# TEM H3K4me1, M->U, 0 lines  less:  1.620178e-06  more:  0.9999984 
# TEM H3K4me1, M->U, 1 lines  less:  1.512131e-10  more:  1 
# TEM H3K4me1, U->M, 0 lines  less:  0.9924358  more:  0.007570642 
# TEM H3K4me1, U->M, 1 lines  less:  0.9998549  more:  0.0001452943 
# TEM H3K27me3, M->U, 0 lines  less:  9.145327e-20  more:  1 
# TEM H3K27me3, M->U, 1 lines  less:  1.407202e-06  more:  0.9999986 
# TEM H3K27me3, U->M, 0 lines  less:  0.05702767  more:  0.9429872 
# TEM H3K27me3, U->M, 1 lines  less:  0.9788822  more:  0.02112446 
# TEM H2AW, M->U, 0 lines  less:  1  more:  1.948109e-91 
# TEM H2AW, M->U, 1 lines  less:  1  more:  2.015542e-29 
# TEM H2AW, U->M, 0 lines  less:  4.466346e-58  more:  1 
# TEM H2AW, U->M, 1 lines  less:  2.291465e-36  more:  1 
	
	
	
	
	
	
	
	
	
	
	### Plotting karyotype-level data
	library(karyoploteR)
	
	# Reload segmentation model to be sure
	#segmentation_model.gr = readRDS(paste0(project_id,"_",meth_context,"_segmentation_model_draft.rds"))

	levels(segmentation_model.gr@seqnames@values) = mixedCaseChr(as.character(levels(segmentation_model.gr@seqnames@values)))
	segmentation_model.gr@seqinfo@seqnames = levels(segmentation_model.gr@seqnames@values)

	levels(variant_call_ranges@seqnames@values) = mixedCaseChr(as.character(levels(variant_call_ranges@seqnames@values)))
	variant_call_ranges@seqinfo@seqnames = levels(variant_call_ranges@seqnames@values)

	karyoplot_params <- getDefaultPlotParams(plot.type=2)
	karyoplot_params$ideogramheight <- 80
	karyoplot_params$data1height <- 50
	karyoplot_params$data2height <- karyoplot_params$ideogramheight
	karyoplot_params$data2inmargin <- -1 * karyoplot_params$ideogramheight



	# Create version of segmentation model with 'Unknown' genes masked out	
	#segment_mask_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	#segment_mask_loci=rbind(segment_mask_loci,list(Chromosome="ChrM",start_site=1,end_site=sLengths[6]),list(Chromosome="ChrC",start_site=1,end_site=sLengths[7]))
	#unknown_gene_loci=gff.genes[gff.genes$m_class=="Unknown",c("V1","V4","V5")]
	#colnames(unknown_gene_loci)=colnames(segment_mask_loci)
	#unknown_gene_loci$Chromosome=paste0(substr(unknown_gene_loci$Chromosome,1,1),tolower(substr(unknown_gene_loci$Chromosome,2,3)),substr(unknown_gene_loci$Chromosome,4,4))
	#segment_mask_loci=rbind.data.frame(segment_mask_loci,unknown_gene_loci)
	#rm(unknown_gene_loci)
	#segment_mask_loci.gr=makeGRangesFromDataFrame(df=segment_mask_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
	#masked_segmentation_model.gr = segmentation_model.gr[unique(queryHits(findOverlaps(segmentation_model.gr, setdiff(segmentation_model.gr,segment_mask_loci.gr))))]
	
	# Create version of segmentation model with Chromosome 2 mitochondrial section masked out. The limits of the masked area are defined here as the start and end of the first and last gene model in the section which have the status 'Unknown', largely as a result of abnormally high coverage (>5sigma).
	segment_mask_loci.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"="Chr2", "start_site"=3239693, "end_site"=3505260))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
	
		
	masked_segmentation_model.gr = segmentation_model.gr[unique(queryHits(findOverlaps(segmentation_model.gr, setdiff(segmentation_model.gr,segment_mask_loci.gr))))]

	# Track of GBM genes	
	#segment_mask_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	#segment_mask_loci=rbind(segment_mask_loci,list(Chromosome="ChrM",start_site=1,end_site=sLengths[6]),list(Chromosome="ChrC",start_site=1,end_site=sLengths[7]))
	GBM_gene_loci=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")]
	colnames(GBM_gene_loci)=c("Chromosome", "start_site", "end_site")
	GBM_gene_loci$Chromosome=mixedCaseChr(GBM_gene_loci$Chromosome)
	GBM_gene_loci.gr=reduce(makeGRangesFromDataFrame(df=GBM_gene_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome"))
	
	# Whole genome:
	karyoplot_title="A. Arabidopsis mCG-variable sites predominantly occur in Gene-Body Methylated genes\n"

	pdf_plot_width=40
	pdf_plot_height=10
	pdf(paste0(project_id,"_genome_segmentation_and_mCG_variable_sites.pdf"), width=pdf_plot_width, height=pdf_plot_height)
	kp <- plotKaryotype(genome="BSgenome.Athaliana.TAIR.TAIR9", plot.type=2, chromosomes=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), zoom=NULL, main=karyoplot_title, plot.params = karyoplot_params)
	kpAddBaseNumbers(kp, tick.dist = 1e6, tick.len = 7, tick.col="black", cex = 1, minor.tick.dist = 1e5, minor.tick.len=4, minor.tick.col = "gray")

	kpDataBackground(kp, data.panel = 1, col="white")
	kpDataBackground(kp, data.panel = 2, col="white")

	# Plot individual segment types separately:
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"], data.panel = 1, col="#0000AA")
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"], data.panel = 1, col="#AA0000")
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"], data.panel = 1, col="#00AA00")

	# Plot segmentation model, coloured by segment type:
	GBM_col = "#00FF00"
	GBM_like_col = "#AAAA00"
	TEM_col = "#FF0000"
	UMR_col = "#0000FF"
	UMR_txp = "AA"  # Render UMR with low opacity so the various types of methylation show through strongly
	other_col="#FFFFFF"
	seg_txp = "CC"

	kpPlotRegions(kp, masked_segmentation_model.gr, data.panel = 2, col = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM-like",paste0(GBM_like_col,seg_txp),paste0(other_col,seg_txp))))), border = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,UMR_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM-like",paste0(GBM_like_col,seg_txp),paste0(other_col,seg_txp))))))
	
	#kpText(kp, chr="Chr1", x=28e6, y=1.25, labels="Genomic segments", data.panel = 2)	
	kpRect(kp, chr="Chr1", x0=26e6, x1=30e6, y0=0.1, y1=0.9, col="white", border="black", data.panel = 2)
	kpText(kp, chr="Chr1", x=28e6, y=0.25, labels="Gene-body-like methylation", data.panel = 2, col="green")	
	kpText(kp, chr="Chr1", x=28e6, y=0.5, labels="TE-like methylation", data.panel = 2, col="red")	
	kpText(kp, chr="Chr1", x=28e6, y=0.75, labels="Unmethylated", data.panel = 2, col="blue")	
	#kpAddLabels(kp, "Methylation segments", r0=NULL, r1=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0)
	#kpAddLabels(kp, "Methylation segments", r0=NULL, r1=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0)

	# Plot GBM genes:
	GBM_gene_col = "#CCCC00"
	GBM_gene_txp = "44"
	kpPlotRegions(kp, GBM_gene_loci.gr, data.panel=1, r0=0, r1=0.5, col=paste0(GBM_gene_col,GBM_gene_txp), border=paste0(GBM_gene_col,GBM_gene_txp))
	
	# Plot variable sites:
	#kpPlotDensity(kp, variant_call_ranges, window_size=1e+04, ymin=NULL, ymax=NULL, data.panel=1, r0=0.5, r1=1, col="#00CCCC")
	var_col = "#00CCCC"
	var_txp = "22"
	kpPlotRegions(kp, variant_call_ranges, data.panel=1, r0=0.5, r1=1, col=paste0(var_col,var_txp), border=paste0(var_col,var_txp))

	kpRect(kp, chr="Chr1", x0=26e6, x1=30e6, y0=0.1, y1=0.9, col="white", border="black", data.panel = 1)
	kpText(kp, chr="Chr1", x=28e6, y=0.3, labels="GBM-annotated genes", col=GBM_gene_col, data.panel = 1)
	kpText(kp, chr="Chr1", x=28e6, y=0.7, labels="mCG-variable sites", col=var_col, data.panel = 1)
	
	#dev.off()


	# Zoomed euchromatic fragment:
	
	karyoplot_title="B. mCG-variable sites occur in the majority of euchromatic Gene-Body Methylated genes\n"

	#pdf(paste0(project_id,"_genome_segmentation_and_mCG_variable_sites_zoom1.pdf"), width=pdf_plot_width, height=pdf_plot_height)
	zoom.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"="Chr1", "start_site"=1.8e6, "end_site"=1.9e6))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
	kp <- plotKaryotype(genome="BSgenome.Athaliana.TAIR.TAIR9", plot.type=2, chromosomes=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), zoom=zoom.gr, main=karyoplot_title, plot.params = karyoplot_params)

	kpAddBaseNumbers(kp, tick.dist = 1e5, tick.len = 5, tick.col="black", cex = 1, minor.tick.dist = 1e4, minor.tick.len=2, minor.tick.col = "gray")

	kpDataBackground(kp, data.panel = 1, col="white")
	kpDataBackground(kp, data.panel = 2, col="white")

	# Plot individual segment types separately:
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"], data.panel = 1, col="#0000AA")
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"], data.panel = 1, col="#AA0000")
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"], data.panel = 1, col="#00AA00")

	# Plot segmentation model, coloured by segment type:
	GBM_col = "#00FF00"
	TEM_col = "#FF0000"
	UMR_col = "#0000FF"
	other_col="#FFFFFF"
	seg_txp = "FF"

	kpPlotRegions(kp, masked_segmentation_model.gr, data.panel = 2, col = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),paste0(other_col,seg_txp)))), border = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),paste0(other_col,seg_txp)))))
	
	#kpText(kp, chr="Chr1", x=28e6, y=1.25, labels="Genomic segments", data.panel = 2)	
	kpRect(kp, chr="Chr1", x0=26e6, x1=30e6, y0=0.1, y1=0.9, col="white", border="black", data.panel = 2)
	kpText(kp, chr="Chr1", x=28e6, y=0.25, labels="Gene-body-like methylation", data.panel = 2, col="green")	
	kpText(kp, chr="Chr1", x=28e6, y=0.5, labels="TE-like methylation", data.panel = 2, col="red")	
	kpText(kp, chr="Chr1", x=28e6, y=0.75, labels="Unmethylated", data.panel = 2, col="blue")	
	#kpAddLabels(kp, "Methylation segments", r0=NULL, r1=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0)
	#kpAddLabels(kp, "Methylation segments", r0=NULL, r1=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0)

	# Plot GBM genes:
	GBM_gene_col = "#CCCC00"
	GBM_gene_txp = "88"
	kpPlotRegions(kp, GBM_gene_loci.gr, data.panel=1, r0=0, r1=0.5, col=paste0(GBM_gene_col,GBM_gene_txp), border=paste0(GBM_gene_col,GBM_gene_txp))
	
	# Plot variable sites:
	#kpPlotDensity(kp, variant_call_ranges, window_size=1e+04, ymin=NULL, ymax=NULL, data.panel=1, r0=0.5, r1=1, col="#00CCCC")
	var_col = "#00CCCC"
	var_txp = "88"
	kpPlotRegions(kp, variant_call_ranges, data.panel=1, r0=0.5, r1=1, col=paste0(var_col,var_txp), border=paste0(var_col,var_txp))

	kpRect(kp, chr="Chr1", x0=26e6, x1=30e6, y0=0.1, y1=0.9, col="white", border="black", data.panel = 1)
	kpText(kp, chr="Chr1", x=28e6, y=0.3, labels="GBM-annotated genes", col=GBM_gene_col, data.panel = 1)
	kpText(kp, chr="Chr1", x=28e6, y=0.7, labels="mCG-variable sites", col=var_col, data.panel = 1)
	
	#dev.off()
	
	# Zoomed heterochromatic fragment:
	
	karyoplot_title="C. mCG-variable sites are rare in heterochromatin\n"

	#pdf(paste0(project_id,"_genome_segmentation_and_mCG_variable_sites_zoom2.pdf"), width=pdf_plot_width, height=pdf_plot_height)
	zoom.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"="Chr1", "start_site"=15.8e6, "end_site"=15.9e6))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
	kp <- plotKaryotype(genome="BSgenome.Athaliana.TAIR.TAIR9", plot.type=2, chromosomes=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), zoom=zoom.gr, main=karyoplot_title, plot.params = karyoplot_params)

	kpAddBaseNumbers(kp, tick.dist = 1e5, tick.len = 5, tick.col="black", cex = 1, minor.tick.dist = 1e4, minor.tick.len=2, minor.tick.col = "gray")

	kpDataBackground(kp, data.panel = 1, col="white")
	kpDataBackground(kp, data.panel = 2, col="white")

	# Plot individual segment types separately:
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"], data.panel = 1, col="#0000AA")
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"], data.panel = 1, col="#AA0000")
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"], data.panel = 1, col="#00AA00")

	# Plot segmentation model, coloured by segment type:
	GBM_col = "#00FF00"
	TEM_col = "#FF0000"
	UMR_col = "#0000FF"
	other_col="#FFFFFF"
	seg_txp = "FF"

	kpPlotRegions(kp, masked_segmentation_model.gr, data.panel = 2, col = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),paste0(other_col,seg_txp)))), border = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),paste0(other_col,seg_txp)))))
	
	#kpText(kp, chr="Chr1", x=28e6, y=1.25, labels="Genomic segments", data.panel = 2)	
	kpRect(kp, chr="Chr1", x0=26e6, x1=30e6, y0=0.1, y1=0.9, col="white", border="black", data.panel = 2)
	kpText(kp, chr="Chr1", x=28e6, y=0.25, labels="Gene-body-like methylation", data.panel = 2, col="green")	
	kpText(kp, chr="Chr1", x=28e6, y=0.5, labels="TE-like methylation", data.panel = 2, col="red")	
	kpText(kp, chr="Chr1", x=28e6, y=0.75, labels="Unmethylated", data.panel = 2, col="blue")	
	#kpAddLabels(kp, "Methylation segments", r0=NULL, r1=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0)
	#kpAddLabels(kp, "Methylation segments", r0=NULL, r1=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0)

	# Plot GBM genes:
	GBM_gene_col = "#CCCC00"
	GBM_gene_txp = "88"
	kpPlotRegions(kp, GBM_gene_loci.gr, data.panel=1, r0=0, r1=0.5, col=paste0(GBM_gene_col,GBM_gene_txp), border=paste0(GBM_gene_col,GBM_gene_txp))
	
	# Plot variable sites:
	#kpPlotDensity(kp, variant_call_ranges, window_size=1e+04, ymin=NULL, ymax=NULL, data.panel=1, r0=0.5, r1=1, col="#00CCCC")
	var_col = "#00CCCC"
	var_txp = "88"
	kpPlotRegions(kp, variant_call_ranges, data.panel=1, r0=0.5, r1=1, col=paste0(var_col,var_txp), border=paste0(var_col,var_txp))

	kpRect(kp, chr="Chr1", x0=26e6, x1=30e6, y0=0.1, y1=0.9, col="white", border="black", data.panel = 1)
	kpText(kp, chr="Chr1", x=28e6, y=0.3, labels="GBM-annotated genes", col=GBM_gene_col, data.panel = 1)
	kpText(kp, chr="Chr1", x=28e6, y=0.7, labels="mCG-variable sites", col=var_col, data.panel = 1)
	
	dev.off()
	

	
	# Plot karyotype view of variable site density
	
	karyoplot_title = "Density of mCG-variable sites (100kb window)"
	
	karyoplot_params <- getDefaultPlotParams(plot.type=1)
 	karyoplot_params$ideogramheight <- 0
 	karyoplot_params$data1height <- 50
 
	karyoplot_params$topmargin <- 10
  	karyoplot_params$bottommargin <- 60
 	
 	pdf_plot_width=10
 	pdf_plot_height=10
 
 	pdf(paste0(project_id,"_mCG_variable_sites_density_100kb.pdf"), width=pdf_plot_width, height=pdf_plot_height)

 	kp <- plotKaryotype(genome="BSgenome.Athaliana.TAIR.TAIR9", plot.type=1, chromosomes=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), zoom=NULL, main=karyoplot_title, plot.params = karyoplot_params)
 	kpAddBaseNumbers(kp, tick.dist = 5e6, tick.len = 7, tick.col="black", cex = 1, minor.tick.dist = 1e6, minor.tick.len=4, minor.tick.col = "gray")
 
 	kpDataBackground(kp, data.panel = 1, col="white")
 
	kpPlotDensity(kp, variant_call_ranges, window.size=1e+05, data.panel=1, r0=-0.4, r1=1, col="#FFFFFF")
	
	dev.off()
	
	

	#library(ggbio)
	
	
	
	
	
	
	# How large can C/C+T get for an individual site that is classified as "U"?
	cat(paste0("CG site with largest proportion of C/(C+T) which is classified as 'U': ",max(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T)),"\n"))
	# Schmitz data: 0.1666667
	
	# What does the distribution of C/(C+T) look like for "U" sites?
	ggplot(data.frame("prop_c"=cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T)), aes(x=prop_c)) +geom_density()
	
	# 99ile of C/(C+T) for "U" sites:
	quantile(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T),.99)
	# Schmitz: 0.07692308
	
	
	
	
	# Now we have defined a segmentation of the genome into UMR, GBM and TEM, let's overlap each of these sets of segments with each of the principal classifications of gene/transposon annotations and estimate the proportions of each class of annotation represented by each class of overlapped mCG segment.  Then we can analyse the distribution of sites whose mCG changes among these segment classes.


	segmentation_model.gr = masked_segmentation_model.gr
	
	mask_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	mask_loci=rbind(mask_loci,list(Chromosome="ChrM",start_site=1,end_site=sLengths[6]),list(Chromosome="ChrC",start_site=1,end_site=sLengths[7]))
	mask_loci.gr=makeGRangesFromDataFrame(df=mask_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome")


	#colnames(mask_loci)=c("Chromosome","start_site","end_site")

	# Commented out the part which masks out genes annotated as mCG 'Unknown'
	#unknown_gene_loci=gff.genes[gff.genes$m_class=="Unknown",c("V1","V4","V5")]
	#colnames(unknown_gene_loci)=colnames(mask_loci)
	#unknown_gene_loci$Chromosome=paste0(substr(unknown_gene_loci$Chromosome,1,1),tolower(substr(unknown_gene_loci$Chromosome,2,3)),substr(unknown_gene_loci$Chromosome,4,4))
	#mask_loci=rbind.data.frame(mask_loci,unknown_gene_loci)
	#rm(unknown_gene_loci)


	
	# Define genomicranges representing the whole genome
	genome_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	levels(genome_loci$Chromosome)=c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrM","ChrC")
	for (chrom_no in 1:length(sLengths)) {
		genome_loci[chrom_no,]=list(Chromosome=names(sLengths[chrom_no]),start_site=1,end_site=sLengths[chrom_no])
	}
	genome_loci.gr=makeGRangesFromDataFrame(df=genome_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome")

	gff.genes$Chromosome = mixedCaseChr(gff.genes$V1)
	gff.transposons$Chromosome = mixedCaseChr(gff.transposons$V1)
	
	
	GBM_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	Het_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Heterochromatic",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	UM_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Unmethylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	other_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.genes[!(gff.genes$m_class %in% c("Gene-body Methylated","Heterochromatic","Unmethylated")),c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	GBM_TE_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.transposons[gff.transposons$m_class=="Gene-body Methylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	Het_TE_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.transposons[gff.transposons$m_class=="Heterochromatic",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	UM_TE_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.transposons[gff.transposons$m_class=="Unmethylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	other_TE_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.transposons[!(gff.transposons$m_class %in% c("Gene-body Methylated","Heterochromatic","Unmethylated")),c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	
	Unannotated_loci.gr = setdiff(genome_loci.gr, c(mask_loci.gr, segment_mask_loci.gr, GBM_gene_loci.gr, Het_gene_loci.gr, UM_gene_loci.gr, other_gene_loci.gr, GBM_TE_loci.gr, Het_TE_loci.gr, UM_TE_loci.gr, other_TE_loci.gr))
	
	# What proportion of GBM segments is encompassed by GBM-annotated genes?
	# What proportion of TEM segments is encompassed by Het-annotated TEs?
	# What proportion of UMR segments is encompassed by UM-annotated loci?
	
	# Make a summary table to hold proportions and variable site counts
	segment_annotation_summary = data.frame()
	
	# Make an overlaps table to capture each class of segment overlap containing variable sites
	segment_annotation_detail = data.frame()
	#colnames(segment_annotation_detail) = c("seqnames","start","end","width","strand","hit")

	# Set levels of single clean call ranges to mixed case
	levels(clean_single_m_to_u_call_ranges@seqnames@values) = mixedCaseChr(as.character(levels(clean_single_m_to_u_call_ranges@seqnames@values)))
	clean_single_m_to_u_call_ranges@seqinfo@seqnames = levels(clean_single_m_to_u_call_ranges@seqnames@values)
	levels(clean_single_u_to_m_call_ranges@seqnames@values) = mixedCaseChr(as.character(levels(clean_single_u_to_m_call_ranges@seqnames@values)))
	clean_single_u_to_m_call_ranges@seqinfo@seqnames = levels(clean_single_u_to_m_call_ranges@seqnames@values)
	levels(all_M_ranges@seqnames@values) = mixedCaseChr(as.character(levels(all_M_ranges@seqnames@values)))
	all_M_ranges@seqinfo@seqnames = levels(all_M_ranges@seqnames@values)
	levels(all_U_ranges@seqnames@values) = mixedCaseChr(as.character(levels(all_U_ranges@seqnames@values)))
	all_U_ranges@seqinfo@seqnames = levels(all_U_ranges@seqnames@values)
	levels(m_to_u_call_ranges@seqnames@values) = mixedCaseChr(as.character(levels(m_to_u_call_ranges@seqnames@values)))
	m_to_u_call_ranges@seqinfo@seqnames = levels(m_to_u_call_ranges@seqnames@values)
	levels(u_to_m_call_ranges@seqnames@values) = mixedCaseChr(as.character(levels(u_to_m_call_ranges@seqnames@values)))
	u_to_m_call_ranges@seqinfo@seqnames = levels(u_to_m_call_ranges@seqnames@values)
	levels(segmentation_model.gr@seqnames@values) = mixedCaseChr(as.character(levels(segmentation_model.gr@seqnames@values)))
	segmentation_model.gr@seqinfo@seqnames = levels(segmentation_model.gr@seqnames@values)

	# Find the overlaps of each set of genomic segments with each set of annotated loci.  Find how many variable sites are encompassed by each class of overlap.
	G_G_G_hits=findOverlaps(GBM_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_G_G_overlaps <- pintersect(GBM_gene_loci.gr[queryHits(G_G_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_G_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM", "annotation_type"="GBM gene", as.data.frame(G_G_G_overlaps)))
	G_G_G_overlap_prop <- sum(width(G_G_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	G_G_G_v_hits = length(unique(subjectHits(findOverlaps(G_G_G_overlaps, variant_call_ranges))))
	G_G_G_mu_hits = length(unique(subjectHits(findOverlaps(G_G_G_overlaps, m_to_u_call_ranges))))
	G_G_G_um_hits = length(unique(subjectHits(findOverlaps(G_G_G_overlaps, u_to_m_call_ranges))))
	G_G_G_csm_hits = length(unique(subjectHits(findOverlaps(G_G_G_overlaps, clean_single_m_to_u_call_ranges))))
	G_G_G_csu_hits = length(unique(subjectHits(findOverlaps(G_G_G_overlaps, clean_single_u_to_m_call_ranges))))
	G_G_G_all_M_hits = length(unique(subjectHits(findOverlaps(G_G_G_overlaps, all_M_ranges))))
	G_G_G_all_U_hits = length(unique(subjectHits(findOverlaps(G_G_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM segments", "annotation_type" = "GBM genes", "overlap" = G_G_G_overlap_prop, "variable_sites" = G_G_G_v_hits, "m_to_u_sites" = G_G_G_mu_hits, "u_to_m_sites" = G_G_G_um_hits, "csm_sites" = G_G_G_csm_hits, "csu_sites" = G_G_G_csu_hits, "all_M_sites" = G_G_G_all_M_hits, "all_U_sites" = G_G_G_all_U_hits))
	names(segment_annotation_summary) = c("segment_type", "annotation_type", "overlap", "variable_sites", "m_to_u_sites", "u_to_m_sites", "csm_sites", "csu_sites", "all_M_sites", "all_U_sites")
	segment_annotation_summary$segment_type = as.character(segment_annotation_summary$segment_type)
	segment_annotation_summary$annotation_type = as.character(segment_annotation_summary$annotation_type)
	segment_annotation_summary$overlap = as.numeric(levels(segment_annotation_summary$overlap))[segment_annotation_summary$overlap]
	segment_annotation_summary$variable_sites = as.integer(levels(segment_annotation_summary$variable_sites))[segment_annotation_summary$variable_sites]
	segment_annotation_summary$m_to_u_sites = as.integer(levels(segment_annotation_summary$m_to_u_sites))[segment_annotation_summary$m_to_u_sites]
	segment_annotation_summary$u_to_m_sites = as.integer(levels(segment_annotation_summary$u_to_m_sites))[segment_annotation_summary$u_to_m_sites]
	segment_annotation_summary$csm_sites = as.integer(levels(segment_annotation_summary$csm_sites))[segment_annotation_summary$csm_sites]
	segment_annotation_summary$csu_sites = as.integer(levels(segment_annotation_summary$csu_sites))[segment_annotation_summary$csu_sites]
	segment_annotation_summary$all_M_sites = as.integer(levels(segment_annotation_summary$all_M_sites))[segment_annotation_summary$all_M_sites]
	segment_annotation_summary$all_U_sites = as.integer(levels(segment_annotation_summary$all_U_sites))[segment_annotation_summary$all_U_sites]

	G_H_G_hits=findOverlaps(Het_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_H_G_overlaps <- pintersect(Het_gene_loci.gr[queryHits(G_H_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_H_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM", "annotation_type"="Het. gene", as.data.frame(G_H_G_overlaps)))
	G_H_G_overlap_prop <- sum(width(G_H_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	G_H_G_v_hits = length(unique(subjectHits(findOverlaps(G_H_G_overlaps, variant_call_ranges))))
	G_H_G_mu_hits = length(unique(subjectHits(findOverlaps(G_H_G_overlaps, m_to_u_call_ranges))))
	G_H_G_um_hits = length(unique(subjectHits(findOverlaps(G_H_G_overlaps, u_to_m_call_ranges))))
	G_H_G_csm_hits = length(unique(subjectHits(findOverlaps(G_H_G_overlaps, clean_single_m_to_u_call_ranges))))
	G_H_G_csu_hits = length(unique(subjectHits(findOverlaps(G_H_G_overlaps, clean_single_u_to_m_call_ranges))))
	G_H_G_all_M_hits = length(unique(subjectHits(findOverlaps(G_H_G_overlaps, all_M_ranges))))
	G_H_G_all_U_hits = length(unique(subjectHits(findOverlaps(G_H_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM segments", "annotation_type" = "Het. genes", "overlap" = G_H_G_overlap_prop, "variable_sites" = G_H_G_v_hits, "m_to_u_sites" = G_H_G_mu_hits, "u_to_m_sites" = G_H_G_um_hits, "csm_sites" = G_H_G_csm_hits, "csu_sites" = G_H_G_csu_hits, "all_M_sites" = G_H_G_all_M_hits, "all_U_sites" = G_H_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	G_U_G_hits=findOverlaps(UM_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_U_G_overlaps <- pintersect(UM_gene_loci.gr[queryHits(G_U_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_U_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM", "annotation_type"="Unmeth. gene", as.data.frame(G_U_G_overlaps)))
	G_U_G_overlap_prop <- sum(width(G_U_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	G_U_G_v_hits = length(unique(subjectHits(findOverlaps(G_U_G_overlaps, variant_call_ranges))))
	G_U_G_mu_hits = length(unique(subjectHits(findOverlaps(G_U_G_overlaps, m_to_u_call_ranges))))
	G_U_G_um_hits = length(unique(subjectHits(findOverlaps(G_U_G_overlaps, u_to_m_call_ranges))))
	G_U_G_csm_hits = length(unique(subjectHits(findOverlaps(G_U_G_overlaps, clean_single_m_to_u_call_ranges))))
	G_U_G_csu_hits = length(unique(subjectHits(findOverlaps(G_U_G_overlaps, clean_single_u_to_m_call_ranges))))
	G_U_G_all_M_hits = length(unique(subjectHits(findOverlaps(G_U_G_overlaps, all_M_ranges))))
	G_U_G_all_U_hits = length(unique(subjectHits(findOverlaps(G_U_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM segments", "annotation_type" = "Unmeth. genes", "overlap" = G_U_G_overlap_prop, "variable_sites" = G_U_G_v_hits, "m_to_u_sites" = G_U_G_mu_hits, "u_to_m_sites" = G_U_G_um_hits, "csm_sites" = G_U_G_csm_hits, "csu_sites" = G_U_G_csu_hits, "all_M_sites" = G_U_G_all_M_hits, "all_U_sites" = G_U_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	G_O_G_hits=findOverlaps(other_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_O_G_overlaps <- pintersect(other_gene_loci.gr[queryHits(G_O_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_O_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM", "annotation_type"="Other gene", as.data.frame(G_O_G_overlaps)))
	G_O_G_overlap_prop <- sum(width(G_O_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	G_O_G_v_hits = length(unique(subjectHits(findOverlaps(G_O_G_overlaps, variant_call_ranges))))
	G_O_G_mu_hits = length(unique(subjectHits(findOverlaps(G_O_G_overlaps, m_to_u_call_ranges))))
	G_O_G_um_hits = length(unique(subjectHits(findOverlaps(G_O_G_overlaps, u_to_m_call_ranges))))
	G_O_G_csm_hits = length(unique(subjectHits(findOverlaps(G_O_G_overlaps, clean_single_m_to_u_call_ranges))))
	G_O_G_csu_hits = length(unique(subjectHits(findOverlaps(G_O_G_overlaps, clean_single_u_to_m_call_ranges))))
	G_O_G_all_M_hits = length(unique(subjectHits(findOverlaps(G_O_G_overlaps, all_M_ranges))))
	G_O_G_all_U_hits = length(unique(subjectHits(findOverlaps(G_O_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM segments", "annotation_type" = "Other genes", "overlap" = G_O_G_overlap_prop, "variable_sites" = G_O_G_v_hits, "m_to_u_sites" = G_O_G_mu_hits, "u_to_m_sites" = G_O_G_um_hits, "csm_sites" = G_O_G_csm_hits, "csu_sites" = G_O_G_csu_hits, "all_M_sites" = G_O_G_all_M_hits, "all_U_sites" = G_O_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	G_G_T_hits=findOverlaps(GBM_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_G_T_overlaps <- pintersect(GBM_TE_loci.gr[queryHits(G_G_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_G_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM", "annotation_type"="GBM TE", as.data.frame(G_G_T_overlaps)))
	G_G_T_overlap_prop <- sum(width(G_G_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	G_G_T_v_hits = length(unique(subjectHits(findOverlaps(G_G_T_overlaps, variant_call_ranges))))
	G_G_T_mu_hits = length(unique(subjectHits(findOverlaps(G_G_T_overlaps, m_to_u_call_ranges))))
	G_G_T_um_hits = length(unique(subjectHits(findOverlaps(G_G_T_overlaps, u_to_m_call_ranges))))
	G_G_T_csm_hits = length(unique(subjectHits(findOverlaps(G_G_T_overlaps, clean_single_m_to_u_call_ranges))))
	G_G_T_csu_hits = length(unique(subjectHits(findOverlaps(G_G_T_overlaps, clean_single_u_to_m_call_ranges))))
	G_G_T_all_M_hits = length(unique(subjectHits(findOverlaps(G_G_T_overlaps, all_M_ranges))))
	G_G_T_all_U_hits = length(unique(subjectHits(findOverlaps(G_G_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM segments", "annotation_type" = "GBM TEs", "overlap" = G_G_T_overlap_prop, "variable_sites" = G_G_T_v_hits, "m_to_u_sites" = G_G_T_mu_hits, "u_to_m_sites" = G_G_T_um_hits, "csm_sites" = G_G_T_csm_hits, "csu_sites" = G_G_T_csu_hits, "all_M_sites" = G_G_T_all_M_hits, "all_U_sites" = G_G_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)
	
	G_H_T_hits=findOverlaps(Het_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_H_T_overlaps <- pintersect(Het_TE_loci.gr[queryHits(G_H_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_H_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM", "annotation_type"="Het. TE", as.data.frame(G_H_T_overlaps)))
	G_H_T_overlap_prop <- sum(width(G_H_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	G_H_T_v_hits = length(unique(subjectHits(findOverlaps(G_H_T_overlaps, variant_call_ranges))))
	G_H_T_mu_hits = length(unique(subjectHits(findOverlaps(G_H_T_overlaps, m_to_u_call_ranges))))
	G_H_T_um_hits = length(unique(subjectHits(findOverlaps(G_H_T_overlaps, u_to_m_call_ranges))))
	G_H_T_csm_hits = length(unique(subjectHits(findOverlaps(G_H_T_overlaps, clean_single_m_to_u_call_ranges))))
	G_H_T_csu_hits = length(unique(subjectHits(findOverlaps(G_H_T_overlaps, clean_single_u_to_m_call_ranges))))
	G_H_T_all_M_hits = length(unique(subjectHits(findOverlaps(G_H_T_overlaps, all_M_ranges))))
	G_H_T_all_U_hits = length(unique(subjectHits(findOverlaps(G_H_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM segments", "annotation_type" = "Het. TEs", "overlap" = G_H_T_overlap_prop, "variable_sites" = G_H_T_v_hits, "m_to_u_sites" = G_H_T_mu_hits, "u_to_m_sites" = G_H_T_um_hits, "csm_sites" = G_H_T_csm_hits, "csu_sites" = G_H_T_csu_hits, "all_M_sites" = G_H_T_all_M_hits, "all_U_sites" = G_H_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	G_U_T_hits=findOverlaps(UM_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_U_T_overlaps <- pintersect(UM_TE_loci.gr[queryHits(G_U_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_U_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM", "annotation_type"="Unmeth. TE", as.data.frame(G_U_T_overlaps)))
	G_U_T_overlap_prop <- sum(width(G_U_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	G_U_T_v_hits = length(unique(subjectHits(findOverlaps(G_U_T_overlaps, variant_call_ranges))))
	G_U_T_mu_hits = length(unique(subjectHits(findOverlaps(G_U_T_overlaps, m_to_u_call_ranges))))
	G_U_T_um_hits = length(unique(subjectHits(findOverlaps(G_U_T_overlaps, u_to_m_call_ranges))))
	G_U_T_csm_hits = length(unique(subjectHits(findOverlaps(G_U_T_overlaps, clean_single_m_to_u_call_ranges))))
	G_U_T_csu_hits = length(unique(subjectHits(findOverlaps(G_U_T_overlaps, clean_single_u_to_m_call_ranges))))
	G_U_T_all_M_hits = length(unique(subjectHits(findOverlaps(G_U_T_overlaps, all_M_ranges))))
	G_U_T_all_U_hits = length(unique(subjectHits(findOverlaps(G_U_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM segments", "annotation_type" = "Unmeth. TEs", "overlap" = G_U_T_overlap_prop, "variable_sites" = G_U_T_v_hits, "m_to_u_sites" = G_U_T_mu_hits, "u_to_m_sites" = G_U_T_um_hits, "csm_sites" = G_U_T_csm_hits, "csu_sites" = G_U_T_csu_hits, "all_M_sites" = G_U_T_all_M_hits, "all_U_sites" = G_U_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	G_O_T_hits=findOverlaps(other_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_O_T_overlaps <- pintersect(other_TE_loci.gr[queryHits(G_O_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_O_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM", "annotation_type"="Other TE", as.data.frame(G_O_T_overlaps)))
	G_O_T_overlap_prop <- sum(width(G_O_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	G_O_T_v_hits = length(unique(subjectHits(findOverlaps(G_O_T_overlaps, variant_call_ranges))))
	G_O_T_mu_hits = length(unique(subjectHits(findOverlaps(G_O_T_overlaps, m_to_u_call_ranges))))
	G_O_T_um_hits = length(unique(subjectHits(findOverlaps(G_O_T_overlaps, u_to_m_call_ranges))))
	G_O_T_csm_hits = length(unique(subjectHits(findOverlaps(G_O_T_overlaps, clean_single_m_to_u_call_ranges))))
	G_O_T_csu_hits = length(unique(subjectHits(findOverlaps(G_O_T_overlaps, clean_single_u_to_m_call_ranges))))
	G_O_T_all_M_hits = length(unique(subjectHits(findOverlaps(G_O_T_overlaps, all_M_ranges))))
	G_O_T_all_U_hits = length(unique(subjectHits(findOverlaps(G_O_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM segments", "annotation_type" = "Other TEs", "overlap" = G_O_T_overlap_prop, "variable_sites" = G_O_T_v_hits, "m_to_u_sites" = G_O_T_mu_hits, "u_to_m_sites" = G_O_T_um_hits, "csm_sites" = G_O_T_csm_hits, "csu_sites" = G_O_T_csu_hits, "all_M_sites" = G_O_T_all_M_hits, "all_U_sites" = G_O_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	G_Unann_hits=findOverlaps(Unannotated_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_Unann_overlaps <- pintersect(Unannotated_loci.gr[queryHits(G_Unann_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_Unann_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM", "annotation_type"="Unannotated", as.data.frame(G_Unann_overlaps)))
	G_Unann_overlap_prop <- sum(width(G_Unann_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	G_Unann_v_hits = length(unique(subjectHits(findOverlaps(G_Unann_overlaps, variant_call_ranges))))
	G_Unann_mu_hits = length(unique(subjectHits(findOverlaps(G_Unann_overlaps, m_to_u_call_ranges))))
	G_Unann_um_hits = length(unique(subjectHits(findOverlaps(G_Unann_overlaps, u_to_m_call_ranges))))
	G_Unann_csm_hits = length(unique(subjectHits(findOverlaps(G_Unann_overlaps, clean_single_m_to_u_call_ranges))))
	G_Unann_csu_hits = length(unique(subjectHits(findOverlaps(G_Unann_overlaps, clean_single_u_to_m_call_ranges))))
	G_Unann_all_M_hits = length(unique(subjectHits(findOverlaps(G_Unann_overlaps, all_M_ranges))))
	G_Unann_all_U_hits = length(unique(subjectHits(findOverlaps(G_Unann_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM segments", "annotation_type" = "Unannotated", "overlap" = G_Unann_overlap_prop, "variable_sites" = G_Unann_v_hits, "m_to_u_sites" = G_Unann_mu_hits, "u_to_m_sites" = G_Unann_um_hits, "csm_sites" = G_Unann_csm_hits, "csu_sites" = G_Unann_csu_hits, "all_M_sites" = G_Unann_all_M_hits, "all_U_sites" = G_Unann_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	L_G_G_hits=findOverlaps(GBM_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"])
	L_G_G_overlaps <- pintersect(GBM_gene_loci.gr[queryHits(L_G_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"][subjectHits(L_G_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM-like", "annotation_type"="GBM gene", as.data.frame(L_G_G_overlaps)))
	L_G_G_overlap_prop <- sum(width(L_G_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]))
	L_G_G_v_hits = length(unique(subjectHits(findOverlaps(L_G_G_overlaps, variant_call_ranges))))
	L_G_G_mu_hits = length(unique(subjectHits(findOverlaps(L_G_G_overlaps, m_to_u_call_ranges))))
	L_G_G_um_hits = length(unique(subjectHits(findOverlaps(L_G_G_overlaps, u_to_m_call_ranges))))
	L_G_G_csm_hits = length(unique(subjectHits(findOverlaps(L_G_G_overlaps, clean_single_m_to_u_call_ranges))))
	L_G_G_csu_hits = length(unique(subjectHits(findOverlaps(L_G_G_overlaps, clean_single_u_to_m_call_ranges))))
	L_G_G_all_M_hits = length(unique(subjectHits(findOverlaps(L_G_G_overlaps, all_M_ranges))))
	L_G_G_all_U_hits = length(unique(subjectHits(findOverlaps(L_G_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM-like segments", "annotation_type" = "GBM genes", "overlap" = L_G_G_overlap_prop, "variable_sites" = L_G_G_v_hits, "m_to_u_sites" = L_G_G_mu_hits, "u_to_m_sites" = L_G_G_um_hits, "csm_sites" = L_G_G_csm_hits, "csu_sites" = L_G_G_csu_hits, "all_M_sites" = L_G_G_all_M_hits, "all_U_sites" = L_G_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)
	
	L_H_G_hits=findOverlaps(Het_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"])
	L_H_G_overlaps <- pintersect(Het_gene_loci.gr[queryHits(L_H_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"][subjectHits(L_H_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM-like", "annotation_type"="Het. gene", as.data.frame(L_H_G_overlaps)))
	L_H_G_overlap_prop <- sum(width(L_H_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]))
	L_H_G_v_hits = length(unique(subjectHits(findOverlaps(L_H_G_overlaps, variant_call_ranges))))
	L_H_G_mu_hits = length(unique(subjectHits(findOverlaps(L_H_G_overlaps, m_to_u_call_ranges))))
	L_H_G_um_hits = length(unique(subjectHits(findOverlaps(L_H_G_overlaps, u_to_m_call_ranges))))
	L_H_G_csm_hits = length(unique(subjectHits(findOverlaps(L_H_G_overlaps, clean_single_m_to_u_call_ranges))))
	L_H_G_csu_hits = length(unique(subjectHits(findOverlaps(L_H_G_overlaps, clean_single_u_to_m_call_ranges))))
	L_H_G_all_M_hits = length(unique(subjectHits(findOverlaps(L_H_G_overlaps, all_M_ranges))))
	L_H_G_all_U_hits = length(unique(subjectHits(findOverlaps(L_H_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM-like segments", "annotation_type" = "Het. genes", "overlap" = L_H_G_overlap_prop, "variable_sites" = L_H_G_v_hits, "m_to_u_sites" = L_H_G_mu_hits, "u_to_m_sites" = L_H_G_um_hits, "csm_sites" = L_H_G_csm_hits, "csu_sites" = L_H_G_csu_hits, "all_M_sites" = L_H_G_all_M_hits, "all_U_sites" = L_H_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	L_U_G_hits=findOverlaps(UM_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"])
	L_U_G_overlaps <- pintersect(UM_gene_loci.gr[queryHits(L_U_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"][subjectHits(L_U_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM-like", "annotation_type"="Unmeth. gene", as.data.frame(L_U_G_overlaps)))
	L_U_G_overlap_prop <- sum(width(L_U_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]))
	L_U_G_v_hits = length(unique(subjectHits(findOverlaps(L_U_G_overlaps, variant_call_ranges))))
	L_U_G_mu_hits = length(unique(subjectHits(findOverlaps(L_U_G_overlaps, m_to_u_call_ranges))))
	L_U_G_um_hits = length(unique(subjectHits(findOverlaps(L_U_G_overlaps, u_to_m_call_ranges))))
	L_U_G_csm_hits = length(unique(subjectHits(findOverlaps(L_U_G_overlaps, clean_single_m_to_u_call_ranges))))
	L_U_G_csu_hits = length(unique(subjectHits(findOverlaps(L_U_G_overlaps, clean_single_u_to_m_call_ranges))))
	L_U_G_all_M_hits = length(unique(subjectHits(findOverlaps(L_U_G_overlaps, all_M_ranges))))
	L_U_G_all_U_hits = length(unique(subjectHits(findOverlaps(L_U_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM-like segments", "annotation_type" = "Unmeth. genes", "overlap" = L_U_G_overlap_prop, "variable_sites" = L_U_G_v_hits, "m_to_u_sites" = L_U_G_mu_hits, "u_to_m_sites" = L_U_G_um_hits, "csm_sites" = L_U_G_csm_hits, "csu_sites" = L_U_G_csu_hits, "all_M_sites" = L_U_G_all_M_hits, "all_U_sites" = L_U_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	L_O_G_hits=findOverlaps(other_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"])
	L_O_G_overlaps <- pintersect(other_gene_loci.gr[queryHits(L_O_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"][subjectHits(L_O_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM-like", "annotation_type"="Other gene", as.data.frame(L_O_G_overlaps)))
	L_O_G_overlap_prop <- sum(width(L_O_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]))
	L_O_G_v_hits = length(unique(subjectHits(findOverlaps(L_O_G_overlaps, variant_call_ranges))))
	L_O_G_mu_hits = length(unique(subjectHits(findOverlaps(L_O_G_overlaps, m_to_u_call_ranges))))
	L_O_G_um_hits = length(unique(subjectHits(findOverlaps(L_O_G_overlaps, u_to_m_call_ranges))))
	L_O_G_csm_hits = length(unique(subjectHits(findOverlaps(L_O_G_overlaps, clean_single_m_to_u_call_ranges))))
	L_O_G_csu_hits = length(unique(subjectHits(findOverlaps(L_O_G_overlaps, clean_single_u_to_m_call_ranges))))
	L_O_G_all_M_hits = length(unique(subjectHits(findOverlaps(L_O_G_overlaps, all_M_ranges))))
	L_O_G_all_U_hits = length(unique(subjectHits(findOverlaps(L_O_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM-like segments", "annotation_type" = "Other genes", "overlap" = L_O_G_overlap_prop, "variable_sites" = L_O_G_v_hits, "m_to_u_sites" = L_O_G_mu_hits, "u_to_m_sites" = L_O_G_um_hits, "csm_sites" = L_O_G_csm_hits, "csu_sites" = L_O_G_csu_hits, "all_M_sites" = L_O_G_all_M_hits, "all_U_sites" = L_O_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	L_G_T_hits=findOverlaps(GBM_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"])
	L_G_T_overlaps <- pintersect(GBM_TE_loci.gr[queryHits(L_G_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"][subjectHits(L_G_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM-like", "annotation_type"="GBM TE", as.data.frame(L_G_T_overlaps)))
	L_G_T_overlap_prop <- sum(width(L_G_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]))
	L_G_T_v_hits = length(unique(subjectHits(findOverlaps(L_G_T_overlaps, variant_call_ranges))))
	L_G_T_mu_hits = length(unique(subjectHits(findOverlaps(L_G_T_overlaps, m_to_u_call_ranges))))
	L_G_T_um_hits = length(unique(subjectHits(findOverlaps(L_G_T_overlaps, u_to_m_call_ranges))))
	L_G_T_csm_hits = length(unique(subjectHits(findOverlaps(L_G_T_overlaps, clean_single_m_to_u_call_ranges))))
	L_G_T_csu_hits = length(unique(subjectHits(findOverlaps(L_G_T_overlaps, clean_single_u_to_m_call_ranges))))
	L_G_T_all_M_hits = length(unique(subjectHits(findOverlaps(L_G_T_overlaps, all_M_ranges))))
	L_G_T_all_U_hits = length(unique(subjectHits(findOverlaps(L_G_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM-like segments", "annotation_type" = "GBM TEs", "overlap" = L_G_T_overlap_prop, "variable_sites" = L_G_T_v_hits, "m_to_u_sites" = L_G_T_mu_hits, "u_to_m_sites" = L_G_T_um_hits, "csm_sites" = L_G_T_csm_hits, "csu_sites" = L_G_T_csu_hits, "all_M_sites" = L_G_T_all_M_hits, "all_U_sites" = L_G_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)
	
	L_H_T_hits=findOverlaps(Het_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"])
	L_H_T_overlaps <- pintersect(Het_TE_loci.gr[queryHits(L_H_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"][subjectHits(L_H_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM-like", "annotation_type"="Het. TE", as.data.frame(L_H_T_overlaps)))
	L_H_T_overlap_prop <- sum(width(L_H_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]))
	L_H_T_v_hits = length(unique(subjectHits(findOverlaps(L_H_T_overlaps, variant_call_ranges))))
	L_H_T_mu_hits = length(unique(subjectHits(findOverlaps(L_H_T_overlaps, m_to_u_call_ranges))))
	L_H_T_um_hits = length(unique(subjectHits(findOverlaps(L_H_T_overlaps, u_to_m_call_ranges))))
	L_H_T_csm_hits = length(unique(subjectHits(findOverlaps(L_H_T_overlaps, clean_single_m_to_u_call_ranges))))
	L_H_T_csu_hits = length(unique(subjectHits(findOverlaps(L_H_T_overlaps, clean_single_u_to_m_call_ranges))))
	L_H_T_all_M_hits = length(unique(subjectHits(findOverlaps(L_H_T_overlaps, all_M_ranges))))
	L_H_T_all_U_hits = length(unique(subjectHits(findOverlaps(L_H_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM-like segments", "annotation_type" = "Het. TEs", "overlap" = L_H_T_overlap_prop, "variable_sites" = L_H_T_v_hits, "m_to_u_sites" = L_H_T_mu_hits, "u_to_m_sites" = L_H_T_um_hits, "csm_sites" = L_H_T_csm_hits, "csu_sites" = L_H_T_csu_hits, "all_M_sites" = L_H_T_all_M_hits, "all_U_sites" = L_H_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	L_U_T_hits=findOverlaps(UM_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"])
	L_U_T_overlaps <- pintersect(UM_TE_loci.gr[queryHits(L_U_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"][subjectHits(L_U_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM-like", "annotation_type"="Unmeth. TE", as.data.frame(L_U_T_overlaps)))
	L_U_T_overlap_prop <- sum(width(L_U_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]))
	L_U_T_v_hits = length(unique(subjectHits(findOverlaps(L_U_T_overlaps, variant_call_ranges))))
	L_U_T_mu_hits = length(unique(subjectHits(findOverlaps(L_U_T_overlaps, m_to_u_call_ranges))))
	L_U_T_um_hits = length(unique(subjectHits(findOverlaps(L_U_T_overlaps, u_to_m_call_ranges))))
	L_U_T_csm_hits = length(unique(subjectHits(findOverlaps(L_U_T_overlaps, clean_single_m_to_u_call_ranges))))
	L_U_T_csu_hits = length(unique(subjectHits(findOverlaps(L_U_T_overlaps, clean_single_u_to_m_call_ranges))))
	L_U_T_all_M_hits = length(unique(subjectHits(findOverlaps(L_U_T_overlaps, all_M_ranges))))
	L_U_T_all_U_hits = length(unique(subjectHits(findOverlaps(L_U_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM-like segments", "annotation_type" = "Unmeth. TEs", "overlap" = L_U_T_overlap_prop, "variable_sites" = L_U_T_v_hits, "m_to_u_sites" = L_U_T_mu_hits, "u_to_m_sites" = L_U_T_um_hits, "csm_sites" = L_U_T_csm_hits, "csu_sites" = L_U_T_csu_hits, "all_M_sites" = L_U_T_all_M_hits, "all_U_sites" = L_U_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	L_O_T_hits=findOverlaps(other_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"])
	L_O_T_overlaps <- pintersect(other_TE_loci.gr[queryHits(L_O_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"][subjectHits(L_O_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM-like", "annotation_type"="Other TE", as.data.frame(L_O_T_overlaps)))
	L_O_T_overlap_prop <- sum(width(L_O_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]))
	L_O_T_v_hits = length(unique(subjectHits(findOverlaps(L_O_T_overlaps, variant_call_ranges))))
	L_O_T_mu_hits = length(unique(subjectHits(findOverlaps(L_O_T_overlaps, m_to_u_call_ranges))))
	L_O_T_um_hits = length(unique(subjectHits(findOverlaps(L_O_T_overlaps, u_to_m_call_ranges))))
	L_O_T_csm_hits = length(unique(subjectHits(findOverlaps(L_O_T_overlaps, clean_single_m_to_u_call_ranges))))
	L_O_T_csu_hits = length(unique(subjectHits(findOverlaps(L_O_T_overlaps, clean_single_u_to_m_call_ranges))))
	L_O_T_all_M_hits = length(unique(subjectHits(findOverlaps(L_O_T_overlaps, all_M_ranges))))
	L_O_T_all_U_hits = length(unique(subjectHits(findOverlaps(L_O_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM-like segments", "annotation_type" = "Other TEs", "overlap" = L_O_T_overlap_prop, "variable_sites" = L_O_T_v_hits, "m_to_u_sites" = L_O_T_mu_hits, "u_to_m_sites" = L_O_T_um_hits, "csm_sites" = L_O_T_csm_hits, "csu_sites" = L_O_T_csu_hits, "all_M_sites" = L_O_T_all_M_hits, "all_U_sites" = L_O_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	L_Unann_hits=findOverlaps(Unannotated_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"])
	L_Unann_overlaps <- pintersect(Unannotated_loci.gr[queryHits(L_Unann_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"][subjectHits(L_Unann_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="GBM-like", "annotation_type"="Unannotated", as.data.frame(L_Unann_overlaps)))
	L_Unann_overlap_prop <- sum(width(L_Unann_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]))
	L_Unann_v_hits = length(unique(subjectHits(findOverlaps(L_Unann_overlaps, variant_call_ranges))))
	L_Unann_mu_hits = length(unique(subjectHits(findOverlaps(L_Unann_overlaps, m_to_u_call_ranges))))
	L_Unann_um_hits = length(unique(subjectHits(findOverlaps(L_Unann_overlaps, u_to_m_call_ranges))))
	L_Unann_csm_hits = length(unique(subjectHits(findOverlaps(L_Unann_overlaps, clean_single_m_to_u_call_ranges))))
	L_Unann_csu_hits = length(unique(subjectHits(findOverlaps(L_Unann_overlaps, clean_single_u_to_m_call_ranges))))
	L_Unann_all_M_hits = length(unique(subjectHits(findOverlaps(L_Unann_overlaps, all_M_ranges))))
	L_Unann_all_U_hits = length(unique(subjectHits(findOverlaps(L_Unann_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "GBM-like segments", "annotation_type" = "Unannotated", "overlap" = L_Unann_overlap_prop, "variable_sites" = L_Unann_v_hits, "m_to_u_sites" = L_Unann_mu_hits, "u_to_m_sites" = L_Unann_um_hits, "csm_sites" = L_Unann_csm_hits, "csu_sites" = L_Unann_csu_hits, "all_M_sites" = L_Unann_all_M_hits, "all_U_sites" = L_Unann_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	T_G_G_hits=findOverlaps(GBM_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"])
	T_G_G_overlaps <- pintersect(GBM_gene_loci.gr[queryHits(T_G_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"][subjectHits(T_G_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="TEM", "annotation_type"="GBM gene", as.data.frame(T_G_G_overlaps)))
	T_G_G_overlap_prop <- sum(width(T_G_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]))
	T_G_G_v_hits = length(unique(subjectHits(findOverlaps(T_G_G_overlaps, variant_call_ranges))))
	T_G_G_mu_hits = length(unique(subjectHits(findOverlaps(T_G_G_overlaps, m_to_u_call_ranges))))
	T_G_G_um_hits = length(unique(subjectHits(findOverlaps(T_G_G_overlaps, u_to_m_call_ranges))))
	T_G_G_csm_hits = length(unique(subjectHits(findOverlaps(T_G_G_overlaps, clean_single_m_to_u_call_ranges))))
	T_G_G_csu_hits = length(unique(subjectHits(findOverlaps(T_G_G_overlaps, clean_single_u_to_m_call_ranges))))
	T_G_G_all_M_hits = length(unique(subjectHits(findOverlaps(T_G_G_overlaps, all_M_ranges))))
	T_G_G_all_U_hits = length(unique(subjectHits(findOverlaps(T_G_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "TEM segments", "annotation_type" = "GBM genes", "overlap" = T_G_G_overlap_prop, "variable_sites" = T_G_G_v_hits, "m_to_u_sites" = T_G_G_mu_hits, "u_to_m_sites" = T_G_G_um_hits, "csm_sites" = T_G_G_csm_hits, "csu_sites" = T_G_G_csu_hits, "all_M_sites" = T_G_G_all_M_hits, "all_U_sites" = T_G_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	T_H_G_hits=findOverlaps(Het_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"])
	T_H_G_overlaps <- pintersect(Het_gene_loci.gr[queryHits(T_H_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"][subjectHits(T_H_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="TEM", "annotation_type"="Het. gene", as.data.frame(T_H_G_overlaps)))
	T_H_G_overlap_prop <- sum(width(T_H_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]))
	T_H_G_v_hits = length(unique(subjectHits(findOverlaps(T_H_G_overlaps, variant_call_ranges))))
	T_H_G_mu_hits = length(unique(subjectHits(findOverlaps(T_H_G_overlaps, m_to_u_call_ranges))))
	T_H_G_um_hits = length(unique(subjectHits(findOverlaps(T_H_G_overlaps, u_to_m_call_ranges))))
	T_H_G_csm_hits = length(unique(subjectHits(findOverlaps(T_H_G_overlaps, clean_single_m_to_u_call_ranges))))
	T_H_G_csu_hits = length(unique(subjectHits(findOverlaps(T_H_G_overlaps, clean_single_u_to_m_call_ranges))))
	T_H_G_all_M_hits = length(unique(subjectHits(findOverlaps(T_H_G_overlaps, all_M_ranges))))
	T_H_G_all_U_hits = length(unique(subjectHits(findOverlaps(T_H_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "TEM segments", "annotation_type" = "Het. genes", "overlap" = T_H_G_overlap_prop, "variable_sites" = T_H_G_v_hits, "m_to_u_sites" = T_H_G_mu_hits, "u_to_m_sites" = T_H_G_um_hits, "csm_sites" = T_H_G_csm_hits, "csu_sites" = T_H_G_csu_hits, "all_M_sites" = T_H_G_all_M_hits, "all_U_sites" = T_H_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	T_U_G_hits=findOverlaps(UM_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"])
	T_U_G_overlaps <- pintersect(UM_gene_loci.gr[queryHits(T_U_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"][subjectHits(T_U_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="TEM", "annotation_type"="Unmeth. gene", as.data.frame(T_U_G_overlaps)))
	T_U_G_overlap_prop <- sum(width(T_U_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]))
	T_U_G_v_hits = length(unique(subjectHits(findOverlaps(T_U_G_overlaps, variant_call_ranges))))
	T_U_G_mu_hits = length(unique(subjectHits(findOverlaps(T_U_G_overlaps, m_to_u_call_ranges))))
	T_U_G_um_hits = length(unique(subjectHits(findOverlaps(T_U_G_overlaps, u_to_m_call_ranges))))
	T_U_G_csm_hits = length(unique(subjectHits(findOverlaps(T_U_G_overlaps, clean_single_m_to_u_call_ranges))))
	T_U_G_csu_hits = length(unique(subjectHits(findOverlaps(T_U_G_overlaps, clean_single_u_to_m_call_ranges))))
	T_U_G_all_M_hits = length(unique(subjectHits(findOverlaps(T_U_G_overlaps, all_M_ranges))))
	T_U_G_all_U_hits = length(unique(subjectHits(findOverlaps(T_U_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "TEM segments", "annotation_type" = "Unmeth. genes", "overlap" = T_U_G_overlap_prop, "variable_sites" = T_U_G_v_hits, "m_to_u_sites" = T_U_G_mu_hits, "u_to_m_sites" = T_U_G_um_hits, "csm_sites" = T_U_G_csm_hits, "csu_sites" = T_U_G_csu_hits, "all_M_sites" = T_U_G_all_M_hits, "all_U_sites" = T_U_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	T_O_G_hits=findOverlaps(other_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"])
	T_O_G_overlaps <- pintersect(other_gene_loci.gr[queryHits(T_O_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"][subjectHits(T_O_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="TEM", "annotation_type"="Other gene", as.data.frame(T_O_G_overlaps)))
	T_O_G_overlap_prop <- sum(width(T_O_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]))
	T_O_G_v_hits = length(unique(subjectHits(findOverlaps(T_O_G_overlaps, variant_call_ranges))))
	T_O_G_mu_hits = length(unique(subjectHits(findOverlaps(T_O_G_overlaps, m_to_u_call_ranges))))
	T_O_G_um_hits = length(unique(subjectHits(findOverlaps(T_O_G_overlaps, u_to_m_call_ranges))))
	T_O_G_csm_hits = length(unique(subjectHits(findOverlaps(T_O_G_overlaps, clean_single_m_to_u_call_ranges))))
	T_O_G_csu_hits = length(unique(subjectHits(findOverlaps(T_O_G_overlaps, clean_single_u_to_m_call_ranges))))
	T_O_G_all_M_hits = length(unique(subjectHits(findOverlaps(T_O_G_overlaps, all_M_ranges))))
	T_O_G_all_U_hits = length(unique(subjectHits(findOverlaps(T_O_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "TEM segments", "annotation_type" = "Other genes", "overlap" = T_O_G_overlap_prop, "variable_sites" = T_O_G_v_hits, "m_to_u_sites" = T_O_G_mu_hits, "u_to_m_sites" = T_O_G_um_hits, "csm_sites" = T_O_G_csm_hits, "csu_sites" = T_O_G_csu_hits, "all_M_sites" = T_O_G_all_M_hits, "all_U_sites" = T_O_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	T_G_T_hits=findOverlaps(GBM_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"])
	T_G_T_overlaps <- pintersect(GBM_TE_loci.gr[queryHits(T_G_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"][subjectHits(T_G_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="TEM", "annotation_type"="GBM TE", as.data.frame(T_G_T_overlaps)))
	T_G_T_overlap_prop <- sum(width(T_G_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]))
	T_G_T_v_hits = length(unique(subjectHits(findOverlaps(T_G_T_overlaps, variant_call_ranges))))
	T_G_T_mu_hits = length(unique(subjectHits(findOverlaps(T_G_T_overlaps, m_to_u_call_ranges))))
	T_G_T_um_hits = length(unique(subjectHits(findOverlaps(T_G_T_overlaps, u_to_m_call_ranges))))
	T_G_T_csm_hits = length(unique(subjectHits(findOverlaps(T_G_T_overlaps, clean_single_m_to_u_call_ranges))))
	T_G_T_csu_hits = length(unique(subjectHits(findOverlaps(T_G_T_overlaps, clean_single_u_to_m_call_ranges))))
	T_G_T_all_M_hits = length(unique(subjectHits(findOverlaps(T_G_T_overlaps, all_M_ranges))))
	T_G_T_all_U_hits = length(unique(subjectHits(findOverlaps(T_G_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "TEM segments", "annotation_type" = "GBM TEs", "overlap" = T_G_T_overlap_prop, "variable_sites" = T_G_T_v_hits, "m_to_u_sites" = T_G_T_mu_hits, "u_to_m_sites" = T_G_T_um_hits, "csm_sites" = T_G_T_csm_hits, "csu_sites" = T_G_T_csu_hits, "all_M_sites" = T_G_T_all_M_hits, "all_U_sites" = T_G_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)
	
	T_H_T_hits=findOverlaps(Het_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"])
	T_H_T_overlaps <- pintersect(Het_TE_loci.gr[queryHits(T_H_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"][subjectHits(T_H_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="TEM", "annotation_type"="Het. TE", as.data.frame(T_H_T_overlaps)))
	T_H_T_overlap_prop <- sum(width(T_H_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]))
	T_H_T_v_hits = length(unique(subjectHits(findOverlaps(T_H_T_overlaps, variant_call_ranges))))
	T_H_T_mu_hits = length(unique(subjectHits(findOverlaps(T_H_T_overlaps, m_to_u_call_ranges))))
	T_H_T_um_hits = length(unique(subjectHits(findOverlaps(T_H_T_overlaps, u_to_m_call_ranges))))
	T_H_T_csm_hits = length(unique(subjectHits(findOverlaps(T_H_T_overlaps, clean_single_m_to_u_call_ranges))))
	T_H_T_csu_hits = length(unique(subjectHits(findOverlaps(T_H_T_overlaps, clean_single_u_to_m_call_ranges))))
	T_H_T_all_M_hits = length(unique(subjectHits(findOverlaps(T_H_T_overlaps, all_M_ranges))))
	T_H_T_all_U_hits = length(unique(subjectHits(findOverlaps(T_H_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "TEM segments", "annotation_type" = "Het. TEs", "overlap" = T_H_T_overlap_prop, "variable_sites" = T_H_T_v_hits, "m_to_u_sites" = T_H_T_mu_hits, "u_to_m_sites" = T_H_T_um_hits, "csm_sites" = T_H_T_csm_hits, "csu_sites" = T_H_T_csu_hits, "all_M_sites" = T_H_T_all_M_hits, "all_U_sites" = T_H_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	T_U_T_hits=findOverlaps(UM_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"])
	T_U_T_overlaps <- pintersect(UM_TE_loci.gr[queryHits(T_U_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"][subjectHits(T_U_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="TEM", "annotation_type"="Unmeth. TE", as.data.frame(T_U_T_overlaps)))
	T_U_T_overlap_prop <- sum(width(T_U_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]))
	T_U_T_v_hits = length(unique(subjectHits(findOverlaps(T_U_T_overlaps, variant_call_ranges))))
	T_U_T_mu_hits = length(unique(subjectHits(findOverlaps(T_U_T_overlaps, m_to_u_call_ranges))))
	T_U_T_um_hits = length(unique(subjectHits(findOverlaps(T_U_T_overlaps, u_to_m_call_ranges))))
	T_U_T_csm_hits = length(unique(subjectHits(findOverlaps(T_U_T_overlaps, clean_single_m_to_u_call_ranges))))
	T_U_T_csu_hits = length(unique(subjectHits(findOverlaps(T_U_T_overlaps, clean_single_u_to_m_call_ranges))))
	T_U_T_all_M_hits = length(unique(subjectHits(findOverlaps(T_U_T_overlaps, all_M_ranges))))
	T_U_T_all_U_hits = length(unique(subjectHits(findOverlaps(T_U_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "TEM segments", "annotation_type" = "Unmeth. TEs", "overlap" = T_U_T_overlap_prop, "variable_sites" = T_U_T_v_hits, "m_to_u_sites" = T_U_T_mu_hits, "u_to_m_sites" = T_U_T_um_hits, "csm_sites" = T_U_T_csm_hits, "csu_sites" = T_U_T_csu_hits, "all_M_sites" = T_U_T_all_M_hits, "all_U_sites" = T_U_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	T_O_T_hits=findOverlaps(other_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"])
	T_O_T_overlaps <- pintersect(other_TE_loci.gr[queryHits(T_O_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"][subjectHits(T_O_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="TEM", "annotation_type"="Other TE", as.data.frame(T_O_T_overlaps)))
	T_O_T_overlap_prop <- sum(width(T_O_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]))
	T_O_T_v_hits = length(unique(subjectHits(findOverlaps(T_O_T_overlaps, variant_call_ranges))))
	T_O_T_mu_hits = length(unique(subjectHits(findOverlaps(T_O_T_overlaps, m_to_u_call_ranges))))
	T_O_T_um_hits = length(unique(subjectHits(findOverlaps(T_O_T_overlaps, u_to_m_call_ranges))))
	T_O_T_csm_hits = length(unique(subjectHits(findOverlaps(T_O_T_overlaps, clean_single_m_to_u_call_ranges))))
	T_O_T_csu_hits = length(unique(subjectHits(findOverlaps(T_O_T_overlaps, clean_single_u_to_m_call_ranges))))
	T_O_T_all_M_hits = length(unique(subjectHits(findOverlaps(T_O_T_overlaps, all_M_ranges))))
	T_O_T_all_U_hits = length(unique(subjectHits(findOverlaps(T_O_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "TEM segments", "annotation_type" = "Other TEs", "overlap" = T_O_T_overlap_prop, "variable_sites" = T_O_T_v_hits, "m_to_u_sites" = T_O_T_mu_hits, "u_to_m_sites" = T_O_T_um_hits, "csm_sites" = T_O_T_csm_hits, "csu_sites" = T_O_T_csu_hits, "all_M_sites" = T_O_T_all_M_hits, "all_U_sites" = T_O_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	T_Unann_hits=findOverlaps(Unannotated_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"])
	T_Unann_overlaps <- pintersect(Unannotated_loci.gr[queryHits(T_Unann_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"][subjectHits(T_Unann_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="TEM", "annotation_type"="Unannotated", as.data.frame(T_Unann_overlaps)))
	T_Unann_overlap_prop <- sum(width(T_Unann_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]))
	T_Unann_v_hits = length(unique(subjectHits(findOverlaps(T_Unann_overlaps, variant_call_ranges))))
	T_Unann_mu_hits = length(unique(subjectHits(findOverlaps(T_Unann_overlaps, m_to_u_call_ranges))))
	T_Unann_um_hits = length(unique(subjectHits(findOverlaps(T_Unann_overlaps, u_to_m_call_ranges))))
	T_Unann_csm_hits = length(unique(subjectHits(findOverlaps(T_Unann_overlaps, clean_single_m_to_u_call_ranges))))
	T_Unann_csu_hits = length(unique(subjectHits(findOverlaps(T_Unann_overlaps, clean_single_u_to_m_call_ranges))))
	T_Unann_all_M_hits = length(unique(subjectHits(findOverlaps(T_Unann_overlaps, all_M_ranges))))
	T_Unann_all_U_hits = length(unique(subjectHits(findOverlaps(T_Unann_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "TEM segments", "annotation_type" = "Unannotated", "overlap" = T_Unann_overlap_prop, "variable_sites" = T_Unann_v_hits, "m_to_u_sites" = T_Unann_mu_hits, "u_to_m_sites" = T_Unann_um_hits, "csm_sites" = T_Unann_csm_hits, "csu_sites" = T_Unann_csu_hits, "all_M_sites" = T_Unann_all_M_hits, "all_U_sites" = T_Unann_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	U_G_G_hits=findOverlaps(GBM_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"])
	U_G_G_overlaps <- pintersect(GBM_gene_loci.gr[queryHits(U_G_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"][subjectHits(U_G_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="UMR", "annotation_type"="GBM gene", as.data.frame(U_G_G_overlaps)))
	U_G_G_overlap_prop <- sum(width(U_G_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]))
	U_G_G_v_hits = length(unique(subjectHits(findOverlaps(U_G_G_overlaps, variant_call_ranges))))
	U_G_G_mu_hits = length(unique(subjectHits(findOverlaps(U_G_G_overlaps, m_to_u_call_ranges))))
	U_G_G_um_hits = length(unique(subjectHits(findOverlaps(U_G_G_overlaps, u_to_m_call_ranges))))
	U_G_G_csm_hits = length(unique(subjectHits(findOverlaps(U_G_G_overlaps, clean_single_m_to_u_call_ranges))))
	U_G_G_csu_hits = length(unique(subjectHits(findOverlaps(U_G_G_overlaps, clean_single_u_to_m_call_ranges))))
	U_G_G_all_M_hits = length(unique(subjectHits(findOverlaps(U_G_G_overlaps, all_M_ranges))))
	U_G_G_all_U_hits = length(unique(subjectHits(findOverlaps(U_G_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "UMR segments", "annotation_type" = "GBM genes", "overlap" = U_G_G_overlap_prop, "variable_sites" = U_G_G_v_hits, "m_to_u_sites" = U_G_G_mu_hits, "u_to_m_sites" = U_G_G_um_hits, "csm_sites" = U_G_G_csm_hits, "csu_sites" = U_G_G_csu_hits, "all_M_sites" = U_G_G_all_M_hits, "all_U_sites" = U_G_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	U_H_G_hits=findOverlaps(Het_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"])
	U_H_G_overlaps <- pintersect(Het_gene_loci.gr[queryHits(U_H_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"][subjectHits(U_H_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="UMR", "annotation_type"="Het. gene", as.data.frame(U_H_G_overlaps)))
	U_H_G_overlap_prop <- sum(width(U_H_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]))
	U_H_G_v_hits = length(unique(subjectHits(findOverlaps(U_H_G_overlaps, variant_call_ranges))))
	U_H_G_mu_hits = length(unique(subjectHits(findOverlaps(U_H_G_overlaps, m_to_u_call_ranges))))
	U_H_G_um_hits = length(unique(subjectHits(findOverlaps(U_H_G_overlaps, u_to_m_call_ranges))))
	U_H_G_csm_hits = length(unique(subjectHits(findOverlaps(U_H_G_overlaps, clean_single_m_to_u_call_ranges))))
	U_H_G_csu_hits = length(unique(subjectHits(findOverlaps(U_H_G_overlaps, clean_single_u_to_m_call_ranges))))
	U_H_G_all_M_hits = length(unique(subjectHits(findOverlaps(U_H_G_overlaps, all_M_ranges))))
	U_H_G_all_U_hits = length(unique(subjectHits(findOverlaps(U_H_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "UMR segments", "annotation_type" = "Het. genes", "overlap" = U_H_G_overlap_prop, "variable_sites" = U_H_G_v_hits, "m_to_u_sites" = U_H_G_mu_hits, "u_to_m_sites" = U_H_G_um_hits, "csm_sites" = U_H_G_csm_hits, "csu_sites" = U_H_G_csu_hits, "all_M_sites" = U_H_G_all_M_hits, "all_U_sites" = U_H_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	U_U_G_hits=findOverlaps(UM_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"])
	U_U_G_overlaps <- pintersect(UM_gene_loci.gr[queryHits(U_U_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"][subjectHits(U_U_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="UMR", "annotation_type"="Unmeth. gene", as.data.frame(U_U_G_overlaps)))
	U_U_G_overlap_prop <- sum(width(U_U_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]))
	U_U_G_v_hits = length(unique(subjectHits(findOverlaps(U_U_G_overlaps, variant_call_ranges))))
	U_U_G_mu_hits = length(unique(subjectHits(findOverlaps(U_U_G_overlaps, m_to_u_call_ranges))))
	U_U_G_um_hits = length(unique(subjectHits(findOverlaps(U_U_G_overlaps, u_to_m_call_ranges))))
	U_U_G_csm_hits = length(unique(subjectHits(findOverlaps(U_U_G_overlaps, clean_single_m_to_u_call_ranges))))
	U_U_G_csu_hits = length(unique(subjectHits(findOverlaps(U_U_G_overlaps, clean_single_u_to_m_call_ranges))))
	U_U_G_all_M_hits = length(unique(subjectHits(findOverlaps(U_U_G_overlaps, all_M_ranges))))
	U_U_G_all_U_hits = length(unique(subjectHits(findOverlaps(U_U_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "UMR segments", "annotation_type" = "Unmeth. genes", "overlap" = U_U_G_overlap_prop, "variable_sites" = U_U_G_v_hits, "m_to_u_sites" = U_U_G_mu_hits, "u_to_m_sites" = U_U_G_um_hits, "csm_sites" = U_U_G_csm_hits, "csu_sites" = U_U_G_csu_hits, "all_M_sites" = U_U_G_all_M_hits, "all_U_sites" = U_U_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	U_O_G_hits=findOverlaps(other_gene_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"])
	U_O_G_overlaps <- pintersect(other_gene_loci.gr[queryHits(U_O_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"][subjectHits(U_O_G_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="UMR", "annotation_type"="Other gene", as.data.frame(U_O_G_overlaps)))
	U_O_G_overlap_prop <- sum(width(U_O_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]))
	U_O_G_v_hits = length(unique(subjectHits(findOverlaps(U_O_G_overlaps, variant_call_ranges))))
	U_O_G_mu_hits = length(unique(subjectHits(findOverlaps(U_O_G_overlaps, m_to_u_call_ranges))))
	U_O_G_um_hits = length(unique(subjectHits(findOverlaps(U_O_G_overlaps, u_to_m_call_ranges))))
	U_O_G_csm_hits = length(unique(subjectHits(findOverlaps(U_O_G_overlaps, clean_single_m_to_u_call_ranges))))
	U_O_G_csu_hits = length(unique(subjectHits(findOverlaps(U_O_G_overlaps, clean_single_u_to_m_call_ranges))))
	U_O_G_all_M_hits = length(unique(subjectHits(findOverlaps(U_O_G_overlaps, all_M_ranges))))
	U_O_G_all_U_hits = length(unique(subjectHits(findOverlaps(U_O_G_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "UMR segments", "annotation_type" = "Other genes", "overlap" = U_O_G_overlap_prop, "variable_sites" = U_O_G_v_hits, "m_to_u_sites" = U_O_G_mu_hits, "u_to_m_sites" = U_O_G_um_hits, "csm_sites" = U_O_G_csm_hits, "csu_sites" = U_O_G_csu_hits, "all_M_sites" = U_O_G_all_M_hits, "all_U_sites" = U_O_G_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	U_G_T_hits=findOverlaps(GBM_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"])
	U_G_T_overlaps <- pintersect(GBM_TE_loci.gr[queryHits(U_G_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"][subjectHits(U_G_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="UMR", "annotation_type"="GBM TE", as.data.frame(U_G_T_overlaps)))
	U_G_T_overlap_prop <- sum(width(U_G_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]))
	U_G_T_v_hits = length(unique(subjectHits(findOverlaps(U_G_T_overlaps, variant_call_ranges))))
	U_G_T_mu_hits = length(unique(subjectHits(findOverlaps(U_G_T_overlaps, m_to_u_call_ranges))))
	U_G_T_um_hits = length(unique(subjectHits(findOverlaps(U_G_T_overlaps, u_to_m_call_ranges))))
	U_G_T_csm_hits = length(unique(subjectHits(findOverlaps(U_G_T_overlaps, clean_single_m_to_u_call_ranges))))
	U_G_T_csu_hits = length(unique(subjectHits(findOverlaps(U_G_T_overlaps, clean_single_u_to_m_call_ranges))))
	U_G_T_all_M_hits = length(unique(subjectHits(findOverlaps(U_G_T_overlaps, all_M_ranges))))
	U_G_T_all_U_hits = length(unique(subjectHits(findOverlaps(U_G_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "UMR segments", "annotation_type" = "GBM TEs", "overlap" = U_G_T_overlap_prop, "variable_sites" = U_G_T_v_hits, "m_to_u_sites" = U_G_T_mu_hits, "u_to_m_sites" = U_G_T_um_hits, "csm_sites" = U_G_T_csm_hits, "csu_sites" = U_G_T_csu_hits, "all_M_sites" = U_G_T_all_M_hits, "all_U_sites" = U_G_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)
	
	U_H_T_hits=findOverlaps(Het_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"])
	U_H_T_overlaps <- pintersect(Het_TE_loci.gr[queryHits(U_H_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"][subjectHits(U_H_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="UMR", "annotation_type"="Het. TE", as.data.frame(U_H_T_overlaps)))
	U_H_T_overlap_prop <- sum(width(U_H_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]))
	U_H_T_v_hits = length(unique(subjectHits(findOverlaps(U_H_T_overlaps, variant_call_ranges))))
	U_H_T_mu_hits = length(unique(subjectHits(findOverlaps(U_H_T_overlaps, m_to_u_call_ranges))))
	U_H_T_um_hits = length(unique(subjectHits(findOverlaps(U_H_T_overlaps, u_to_m_call_ranges))))
	U_H_T_csm_hits = length(unique(subjectHits(findOverlaps(U_H_T_overlaps, clean_single_m_to_u_call_ranges))))
	U_H_T_csu_hits = length(unique(subjectHits(findOverlaps(U_H_T_overlaps, clean_single_u_to_m_call_ranges))))
	U_H_T_all_M_hits = length(unique(subjectHits(findOverlaps(U_H_T_overlaps, all_M_ranges))))
	U_H_T_all_U_hits = length(unique(subjectHits(findOverlaps(U_H_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "UMR segments", "annotation_type" = "Het. TEs", "overlap" = U_H_T_overlap_prop, "variable_sites" = U_H_T_v_hits, "m_to_u_sites" = U_H_T_mu_hits, "u_to_m_sites" = U_H_T_um_hits, "csm_sites" = U_H_T_csm_hits, "csu_sites" = U_H_T_csu_hits, "all_M_sites" = U_H_T_all_M_hits, "all_U_sites" = U_H_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	U_U_T_hits=findOverlaps(UM_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"])
	U_U_T_overlaps <- pintersect(UM_TE_loci.gr[queryHits(U_U_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"][subjectHits(U_U_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="UMR", "annotation_type"="Unmeth. TE", as.data.frame(U_U_T_overlaps)))
	U_U_T_overlap_prop <- sum(width(U_U_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]))
	U_U_T_v_hits = length(unique(subjectHits(findOverlaps(U_U_T_overlaps, variant_call_ranges))))
	U_U_T_mu_hits = length(unique(subjectHits(findOverlaps(U_U_T_overlaps, m_to_u_call_ranges))))
	U_U_T_um_hits = length(unique(subjectHits(findOverlaps(U_U_T_overlaps, u_to_m_call_ranges))))
	U_U_T_csm_hits = length(unique(subjectHits(findOverlaps(U_U_T_overlaps, clean_single_m_to_u_call_ranges))))
	U_U_T_csu_hits = length(unique(subjectHits(findOverlaps(U_U_T_overlaps, clean_single_u_to_m_call_ranges))))
	U_U_T_all_M_hits = length(unique(subjectHits(findOverlaps(U_U_T_overlaps, all_M_ranges))))
	U_U_T_all_U_hits = length(unique(subjectHits(findOverlaps(U_U_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "UMR segments", "annotation_type" = "Unmeth. TEs", "overlap" = U_U_T_overlap_prop, "variable_sites" = U_U_T_v_hits, "m_to_u_sites" = U_U_T_mu_hits, "u_to_m_sites" = U_U_T_um_hits, "csm_sites" = U_U_T_csm_hits, "csu_sites" = U_U_T_csu_hits, "all_M_sites" = U_U_T_all_M_hits, "all_U_sites" = U_U_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	U_O_T_hits=findOverlaps(other_TE_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"])
	U_O_T_overlaps <- pintersect(other_TE_loci.gr[queryHits(U_O_T_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"][subjectHits(U_O_T_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="UMR", "annotation_type"="Other TE", as.data.frame(U_O_T_overlaps)))
	U_O_T_overlap_prop <- sum(width(U_O_T_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]))
	U_O_T_v_hits = length(unique(subjectHits(findOverlaps(U_O_T_overlaps, variant_call_ranges))))
	U_O_T_mu_hits = length(unique(subjectHits(findOverlaps(U_O_T_overlaps, m_to_u_call_ranges))))
	U_O_T_um_hits = length(unique(subjectHits(findOverlaps(U_O_T_overlaps, u_to_m_call_ranges))))
	U_O_T_csm_hits = length(unique(subjectHits(findOverlaps(U_O_T_overlaps, clean_single_m_to_u_call_ranges))))
	U_O_T_csu_hits = length(unique(subjectHits(findOverlaps(U_O_T_overlaps, clean_single_u_to_m_call_ranges))))
	U_O_T_all_M_hits = length(unique(subjectHits(findOverlaps(U_O_T_overlaps, all_M_ranges))))
	U_O_T_all_U_hits = length(unique(subjectHits(findOverlaps(U_O_T_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "UMR segments", "annotation_type" = "Other TEs", "overlap" = U_O_T_overlap_prop, "variable_sites" = U_O_T_v_hits, "m_to_u_sites" = U_O_T_mu_hits, "u_to_m_sites" = U_O_T_um_hits, "csm_sites" = U_O_T_csm_hits, "csu_sites" = U_O_T_csu_hits, "all_M_sites" = U_O_T_all_M_hits, "all_U_sites" = U_O_T_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	U_Unann_hits=findOverlaps(Unannotated_loci.gr, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"])
	U_Unann_overlaps <- pintersect(Unannotated_loci.gr[queryHits(U_Unann_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"][subjectHits(U_Unann_hits)])
	segment_annotation_detail = rbind.data.frame(segment_annotation_detail, cbind("segment_type"="UMR", "annotation_type"="Unannotated", as.data.frame(U_Unann_overlaps)))
	U_Unann_overlap_prop <- sum(width(U_Unann_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]))
	U_Unann_v_hits = length(unique(subjectHits(findOverlaps(U_Unann_overlaps, variant_call_ranges))))
	U_Unann_mu_hits = length(unique(subjectHits(findOverlaps(U_Unann_overlaps, m_to_u_call_ranges))))
	U_Unann_um_hits = length(unique(subjectHits(findOverlaps(U_Unann_overlaps, u_to_m_call_ranges))))
	U_Unann_csm_hits = length(unique(subjectHits(findOverlaps(U_Unann_overlaps, clean_single_m_to_u_call_ranges))))
	U_Unann_csu_hits = length(unique(subjectHits(findOverlaps(U_Unann_overlaps, clean_single_u_to_m_call_ranges))))
	U_Unann_all_M_hits = length(unique(subjectHits(findOverlaps(U_Unann_overlaps, all_M_ranges))))
	U_Unann_all_U_hits = length(unique(subjectHits(findOverlaps(U_Unann_overlaps, all_U_ranges))))
	segment_annotation_summary=rbind.data.frame(segment_annotation_summary, c("segment_type" = "UMR segments", "annotation_type" = "Unannotated", "overlap" = U_Unann_overlap_prop, "variable_sites" = U_Unann_v_hits, "m_to_u_sites" = U_Unann_mu_hits, "u_to_m_sites" = U_Unann_um_hits, "csm_sites" = U_Unann_csm_hits, "csu_sites" = U_Unann_csu_hits, "all_M_sites" = U_Unann_all_M_hits, "all_U_sites" = U_Unann_all_U_hits))
	segment_annotation_summary$overlap = as.numeric(segment_annotation_summary$overlap)
	segment_annotation_summary$variable_sites = as.integer(segment_annotation_summary$variable_sites)
	segment_annotation_summary$m_to_u_sites = as.integer(segment_annotation_summary$m_to_u_sites)
	segment_annotation_summary$u_to_m_sites = as.integer(segment_annotation_summary$u_to_m_sites)
	segment_annotation_summary$csm_sites = as.integer(segment_annotation_summary$csm_sites)
	segment_annotation_summary$csu_sites = as.integer(segment_annotation_summary$csu_sites)
	segment_annotation_summary$all_M_sites = as.integer(segment_annotation_summary$all_M_sites)
	segment_annotation_summary$all_U_sites = as.integer(segment_annotation_summary$all_U_sites)

	# Write out segment/annotation overlaps for visual checking
	write.table(segment_annotation_detail, file=paste0(project_id,"_segment_annotation_overlaps.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
	
	pdf(paste0(project_id,"_",meth_context,"_segment_annotation_overlaps.pdf"))
	print(ggplot(segment_annotation_summary, aes(x=segment_type, y=overlap)) + geom_bar(aes(fill=annotation_type), stat="identity") + labs(x=paste0("Segment classification"), y="Segment proportion overlapping annotated loci", fill="Annotation\nclass") + facet_grid(~segment_type, scales="free") + theme_minimal())
	print(ggplot(segment_annotation_summary, aes(x=segment_type, y=variable_sites)) + geom_bar(aes(fill=annotation_type), stat="identity") + labs(x=paste0("Segment classification"), y="Number of mCG-variable sites overlapping annotated loci", fill="Annotation\nclass") + facet_grid(~segment_type, scales="free") + theme_minimal())
	#print(ggplot(segment_annotation_summary, aes(x=annotation_type, y=overlap)) + geom_bar(aes(fill=segment_type), stat="identity") + labs(x=paste0("Annotation classification"), y="Segment proportion overlapping annotated loci", fill="Segment\nclass") + facet_grid(~annotation_type, scales="free") + theme_minimal())
	print(ggplot(segment_annotation_summary, aes(x=annotation_type, y=variable_sites)) + geom_bar(aes(fill=segment_type), stat="identity") + labs(x=paste0("Annotation classification"), y="Number of mCG-variable sites overlapping annotated loci", fill="Segment\nclass") + facet_grid(~annotation_type, scales="free") + theme_minimal())
	dev.off()

	# Summarise proportions of sites changing in each direction by segment type
	library(dplyr)
	
	z = rbind.data.frame(cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(csm_sites)/sum(csm_sites+all_M_sites))),"metric"=rep("Fraction of mCG sites losing methylation",4)),cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(csu_sites)/sum(csu_sites+all_U_sites))),"metric"=rep("Fraction of uCG sites gaining methylation",4)))

	z[z$segment_type == "GBM segments",]$segment_type = "gbM segments"
	z[z$segment_type == "GBM-like segments",]$segment_type = "gbM-like segments"
	z[z$segment_type == "UMR segments",]$segment_type = "UM segments"
	
		
	# NB. plotting this as a PDF reintroduces the legend and messes up the facet y axis limits, so need to plot it and save the PDF separately
	pdf(paste0(project_id,"_",meth_context,"_delta_mCG_fractions_by_segment_type.pdf"))
	print(ggplot(z) + geom_bar(aes(x=segment_type, y=score, fill=segment_type), stat="identity") + facet_grid(metric ~ ., scales="fixed") + labs(x="", y="") + theme_minimal())
	
	
	plot1 = ggplot(z[1:4,]) + geom_bar(aes(x=factor(segment_type, levels=c("UM segments", "gbM segments", "gbM-like segments", "TEM segments")), y=score, fill=segment_type), stat="identity") + facet_grid(metric ~ .) + labs(x=NULL, y=NULL) + ylim(0,0.3) + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme(legend.position="none")
	plot2 = ggplot(z[5:8,]) + geom_bar(aes(x=factor(segment_type, levels=c("UM segments", "gbM segments", "gbM-like segments", "TEM segments")), y=score, fill=segment_type), stat="identity") + facet_grid(metric ~ .) + labs(x=NULL, y=NULL) + ylim(-0.3,0) + theme_minimal() + theme(legend.position="none")
	gridExtra::grid.arrange(plot1,plot2,nrow=2)
	dev.off()

	
	
#segment_type      mean_csm mean_csu
#  <chr>                <dbl>    <dbl>
#1 GBM-like segments  0.0771   0.186  
#2 GBM segments       0.168    0.159  
#3 TEM segments       0.00214  0.205  
#4 UMR segments       0.259    0.00909


	segment_annotation_summary %>% group_by(annotation_type) %>% summarise(mean_csm = sum(csm_sites)/sum(csm_sites+all_M_sites), mean_csu = sum(csu_sites)/sum(csu_sites+all_U_sites))

#annotation_type mean_csm mean_csu
#  <chr>              <dbl>    <dbl>
#1 GBM genes        0.157    0.0371 
#2 GBM TEs          0.0355   0.0373 
#3 Het. genes       0.00579  0.194  
#4 Het. TEs         0.00129  0.235  
#5 Other genes      0        0      
#6 Other TEs        0.0115   0.0303 
#7 Unannotated      0.00610  0.00479
#8 Unmeth. genes    1        0.00403
#9 Unmeth. TEs      1        0.00805


	
	### Here is where we print out the figures for the paper
	
	# Figure Xa - karyoplot of density of variable sites in chrom 1 and zoom in to a few genes
	
	pdf_plot_height = 10
	pdf_plot_width = 10
	pdf(paste0(project_id,"_",meth_context,"_fig_Xa.pdf"), width=pdf_plot_width, height=pdf_plot_height)
	karyoplot_title = "Fig. Xa mCG-variable sites are associated with euchromatin(density of mCG-variable sites per 100kb window)"
	
	karyoplot_params <- getDefaultPlotParams(plot.type=1)
 	karyoplot_params$ideogramheight <- 0
 	karyoplot_params$data1height <- 50
 
	karyoplot_params$topmargin <- 10
  	karyoplot_params$bottommargin <- 60
 	
 	pdf_plot_width=10
 	pdf_plot_height=10
 
 	kp <- plotKaryotype(genome="BSgenome.Athaliana.TAIR.TAIR9", plot.type=1, chromosomes=c("Chr1","Chr2","Chr3","Chr4","Chr5"), zoom=NULL, main=karyoplot_title, plot.params = karyoplot_params)
 	kpAddBaseNumbers(kp, tick.dist = 5e6, tick.len = 7, tick.col="black", cex = 1, minor.tick.dist = 1e6, minor.tick.len=4, minor.tick.col = "gray")
 
 	kpDataBackground(kp, data.panel = 1, col="white")
	
	# Plot density of variable sites as a line and plot density of genes as a band
	kpPlotDensity(kp, GenomicRanges::setdiff(gene_ranges, segment_mask_loci.gr, ignore.strand=TRUE), window.size=1e+05, data.panel=1, r0=-0.4, r1=1, col="#00FF0044", border="#00FF0066")
	kpPlotDensity(kp, GenomicRanges::setdiff(variant_call_ranges, segment_mask_loci.gr, ignore.strand=TRUE), window.size=1e+05, data.panel=1, r0=-0.4, r1=1, col="#FFFFFF00", border="#000000FF")
	kpAxis(kp, data.panel=1, r0=-0.4, r1=1)
	kpText(kp, chr="Chr1", x=15e6, y=1, col="black", labels="Density of mCG changes", cex=0.8, clipping=FALSE)
	kpText(kp, chr="Chr1", x=15e6, y=0.8, col="#00FF00FF", labels="Density of annotated genes", cex=0.8, clipping=FALSE)
	
	dev.off()

	# Figure Xb - zoomed segment of chromosome 1 showing gene bodies, segmentation and variable sites
	pdf(paste0(project_id,"_",meth_context,"_fig_Xb1.pdf"))

	levels(m_to_u_call_ranges@seqnames@values) = mixedCaseChr(as.character(levels(m_to_u_call_ranges@seqnames@values)))
	m_to_u_call_ranges@seqinfo@seqnames = levels(m_to_u_call_ranges@seqnames@values)
	levels(u_to_m_call_ranges@seqnames@values) = mixedCaseChr(as.character(levels(u_to_m_call_ranges@seqnames@values)))
	u_to_m_call_ranges@seqinfo@seqnames = levels(u_to_m_call_ranges@seqnames@values)
	levels(gene_ranges@seqnames@values) = mixedCaseChr(as.character(levels(gene_ranges@seqnames@values)))
	gene_ranges@seqinfo@seqnames = levels(gene_ranges@seqnames@values)

	karyoplot_title="Fig. Xb. mCG-variable sites are associated with gbM\n"

	karyoplot_params <- getDefaultPlotParams(plot.type=2)
	karyoplot_params$ideogramheight <- 10
	karyoplot_params$data1height <- 100
	karyoplot_params$data2height <- 10
	karyoplot_params$data2inmargin <- -1 * karyoplot_params$ideogramheight

	#pdf(paste0(project_id,"_genome_segmentation_and_mCG_variable_sites_zoom1.pdf"), width=pdf_plot_width, height=pdf_plot_height)
	zoom.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"="Chr1", "start_site"=19.8705e6, "end_site"=19.890e6))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
	kp <- plotKaryotype(genome="BSgenome.Athaliana.TAIR.TAIR9", plot.type=2, chromosomes=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), zoom=zoom.gr, main=karyoplot_title, plot.params = karyoplot_params)

	kpAddBaseNumbers(kp, tick.dist = 1e5, tick.len = 5, tick.col="black", cex = 1, minor.tick.dist = 1e4, minor.tick.len=2, minor.tick.col = "gray")

	kpDataBackground(kp, data.panel = 1, col="white")
	kpDataBackground(kp, data.panel = 2, col="white")

	# Plot individual segment types separately:
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"], data.panel = 1, col="#0000AA")
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"], data.panel = 1, col="#AA0000")
	#kpPlotRegions(kp, segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"], data.panel = 1, col="#00AA00")

	# Plot segmentation model, coloured by segment type:
	GBM_col = "#00FF00"
	TEM_col = "#FF0000"
	UMR_col = "#0000FF"
	other_col="#FFFFFF"
	seg_txp = "FF"

	kpPlotRegions(kp, masked_segmentation_model.gr, data.panel = 2, col = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),paste0(other_col,seg_txp)))), border = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),paste0(other_col,seg_txp)))))

	#kpText(kp, chr="Chr1", x=28e6, y=1.25, labels="Genomic segments", data.panel = 2)	
	kpRect(kp, chr="Chr1", x0=26e6, x1=30e6, y0=0.1, y1=0.9, col="white", border="black", data.panel = 2)
	kpText(kp, chr="Chr1", x=28e6, y=0.25, labels="Gene-body-like methylation", data.panel = 2, col="green")	
	kpText(kp, chr="Chr1", x=28e6, y=0.5, labels="TE-like methylation", data.panel = 2, col="red")	
	kpText(kp, chr="Chr1", x=28e6, y=0.75, labels="Unmethylated", data.panel = 2, col="blue")	
	#kpAddLabels(kp, "Methylation segments", r0=NULL, r1=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0)
	#kpAddLabels(kp, "Methylation segments", r0=NULL, r1=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0)

	# Plot GBM genes:
	GBM_gene_col = "#CCCC00"
	GBM_gene_txp = "88"
	#kpPlotRegions(kp, GBM_gene_loci.gr, data.panel=1, r0=0, r1=0.1, col=paste0(GBM_gene_col,GBM_gene_txp), border=paste0(GBM_gene_col,GBM_gene_txp))
	kpPlotRegions(kp, gene_ranges, data.panel=1, r0=0.1, r1=0.5, col=paste0(GBM_gene_col,GBM_gene_txp), border=paste0(GBM_gene_col,GBM_gene_txp))
	kpPlotMarkers(kp, data=gene_ranges, labels=gff.genes$gene_ID, data.panel=1, r0=0, r1=0.3, text.orientation="vertical", line.color="#FFFFFF00", clipping=FALSE)

	# Plot variable sites:
	#kpPlotDensity(kp, variant_call_ranges, window_size=1e+04, ymin=NULL, ymax=NULL, data.panel=1, r0=0.5, r1=1, col="#00CCCC")
	#kpPlotDensity(kp, variant_call_ranges, window_size=1e+04, ymin=NULL, ymax=NULL, data.panel=1, r0=0.5, r1=1, col="#00CCCC")
	#var_col = "#00CCCC"
	var_col = UMR_col
	var_txp = "88"
	#kpPlotRegions(kp, variant_call_ranges, data.panel=1, r0=0, r1=0.1, col=paste0(var_col,var_txp), border=paste0(var_col,var_txp))
	kpPlotRegions(kp, m_to_u_call_ranges, data.panel=1, r0=0, r1=0.1, col=paste0(var_col,var_txp), border=paste0(var_col,var_txp))
	#var_col = "#00CCCC"
	var_col = TEM_col
	var_txp = "88"
	#kpPlotRegions(kp, variant_call_ranges, data.panel=1, r0=0, r1=0.1, col=paste0(var_col,var_txp), border=paste0(var_col,var_txp))
	kpPlotRegions(kp, u_to_m_call_ranges, data.panel=1, r0=0, r1=0.1, col=paste0(var_col,var_txp), border=paste0(var_col,var_txp))

	kpRect(kp, chr="Chr1", x0=26e6, x1=30e6, y0=0.1, y1=0.9, col="white", border="black", data.panel = 1)
	kpText(kp, chr="Chr1", x=28e6, y=0.3, labels="GBM-annotated genes", col=GBM_gene_col, data.panel = 1)
	kpText(kp, chr="Chr1", x=28e6, y=0.7, labels="mCG-variable sites", col=var_col, data.panel = 1)
	dev.off()
	
	
	# Revised Fig. Xb:
	# Showing the pattern of parental and offspring methylation at sample loci
	pdf(paste0(project_id,"_",meth_context,"_fig_Xb.pdf"), width=pdf_plot_width, height=pdf_plot_height)

	# Query the examples based on changes and coverage
	#genes_of_interest = gff.genes[(gff.genes$variable_count>6) & (gff.genes$m_to_u_count>=2) & (gff.genes$u_to_m_count>=2) & ((gff.genes$all_M_count+gff.genes$all_U_count+gff.genes$variable_count)/gff.genes$CG_sites_count > 0.8),]

	# Query the examples based on a gene list
	genes_of_interest = gff.genes[(gff.genes$gene_ID %in% c("AT5G37190","AT5G40450","AT5G52280","AT1G69070","AT2G35030","AT3G22790","AT4G14760","AT5G39680","AT5G65560","AT1G77580")) ,]
	#genes_of_interest = gff.genes[(gff.genes$gene_ID %in% c("AT5G52280")) ,]

	
	locus_chr = genes_of_interest$Chromosome
	locus_start = genes_of_interest$V4
	locus_end = genes_of_interest$V5
	locus_name = genes_of_interest$gene_ID
	
	UM_once_col = "#FF0000"
	UM_many_col = "#AA0000"
	MU_once_col = "#0000FF"
	MU_many_col = "#0000AA"
	M_no_change_col = "#00FF00"
	U_no_change_col = "#00FF00"

	meth_col = "#CC0000"
	unmeth_col = "#0000CC"
	na_col = "#CCCCCC"
	
	var_txp = "FF"
		
	# Do we want to scale the plot per bp or per CG site?
	#plot_scale = "bp"  # this version plots the locus scaled per bp (chromosome coordinates)
	plot_scale = "cg"  # this version plots the locus scaled per CG site (rowname from all_reps_meth_status table stands as proxy for coordinates - is there a problem with higher chromosomes that they do not have as many bp as there are CG sites in the whole genome?)
	
	for (example_locus in 1:length(locus_name)) {
		locus_data=cbind.data.frame(all_reps_meth_status, parental_consensus, num_lines_m_to_u_from_parentals, num_lines_u_to_m_from_parentals)[(all_reps_meth_status$Chromosome==toupper(locus_chr[example_locus])) & (all_reps_meth_status$Locus>=locus_start[example_locus]) & (all_reps_meth_status$Locus<=locus_end[example_locus]),]
	
		if (plot_scale == "bp") {
			# This version turns the locus_data into a genomic range on bp scale
			locus_data.gr = makeGRangesFromDataFrame(df=cbind.data.frame("Chromosome" = mixedCaseChr(locus_data$Chromosome), "start_site" = locus_data$Locus, "end_site" = locus_data$Locus), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
		} else {
			# This version turns the locus_data into a genomic range scaled per CG site
			locus_data.gr = makeGRangesFromDataFrame(df=cbind.data.frame("Chromosome" = mixedCaseChr(locus_data$Chromosome), "start_site" = rownames(locus_data), "end_site" = rownames(locus_data)), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
		}
		
		cat(paste0(nrow(locus_data),"\n"))
		
		#levels(m_to_u_call_ranges@seqnames@values) = mixedCaseChr(as.character(levels(m_to_u_call_ranges@seqnames@values)))
		#m_to_u_call_ranges@seqinfo@seqnames = levels(m_to_u_call_ranges@seqnames@values)
		#levels(u_to_m_call_ranges@seqnames@values) = mixedCaseChr(as.character(levels(u_to_m_call_ranges@seqnames@values)))
		#u_to_m_call_ranges@seqinfo@seqnames = levels(u_to_m_call_ranges@seqnames@values)
		#levels(gene_ranges@seqnames@values) = mixedCaseChr(as.character(levels(gene_ranges@seqnames@values)))
		#gene_ranges@seqinfo@seqnames = levels(gene_ranges@seqnames@values)

		karyoplot_title=paste0("Fig. Xb. Spontaneous losses and non-templated gains of mCG at individual sites\n are characteristic of gbM loci (locus example ",locus_name[example_locus],")\n")

		karyoplot_params <- getDefaultPlotParams(plot.type=2)
		karyoplot_params$ideogramheight <- 0
		karyoplot_params$data1height <- 100
		karyoplot_params$data2height <- 100
		karyoplot_params$data2inmargin <- -1 * karyoplot_params$ideogramheight

		if (plot_scale == "bp") {
			# This version makes a zoom scaled per bp
			zoom.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"=locus_chr[example_locus], "start_site"=locus_start[example_locus], "end_site"=locus_end[example_locus]))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
		} else {
			# This version expresses the genomic range as adjacent 1bp segments, each representing an adjacent CG site, numbered by their position in the original CG sites table
			zoom.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"=locus_chr[example_locus], "start_site"=min(rownames(locus_data)), "end_site"=max(rownames(locus_data))))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
		}

		kp <- plotKaryotype(genome="BSgenome.Athaliana.TAIR.TAIR9", plot.type=2, chromosomes=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), zoom=zoom.gr, main=karyoplot_title, plot.params = karyoplot_params, ideogram.plotter=NULL, labels.plotter = NULL)

		#kpAddBaseNumbers(kp, tick.dist = 100, tick.len = 5, tick.col="black", cex = 1, minor.tick.dist = 1e4, minor.tick.len=2, minor.tick.col = "gray")

		kpDataBackground(kp, data.panel = 1, col="white")
		kpDataBackground(kp, data.panel = 2, col="white")

		# Plot site status
		var_col = ifelse(locus_data$num_lines_m_to_u_from_parentals==1,MU_once_col,ifelse(locus_data$num_lines_m_to_u_from_parentals>1,MU_many_col,ifelse(locus_data$num_lines_u_to_m_from_parentals==1,UM_once_col,ifelse(locus_data$num_lines_u_to_m_from_parentals>1,UM_many_col,ifelse(locus_data$parental_consensus=="M", M_no_change_col, ifelse(locus_data$parental_consensus=="U", U_no_change_col, na_col))))))
		var_txp = "FF"
		kpPlotRegions(kp, locus_data.gr, data.panel=1, r0=0, r1=0.1, col=paste0(var_col,var_txp), border=paste0(var_col,var_txp))

		# Plot segmentation model, coloured by segment type:
		#GBM_col = "#00FF00"
		#TEM_col = "#FF0000"
		#UMR_col = "#0000FF"
		#other_col="#FFFFFF"
		#seg_txp = "FF"

		#kpPlotRegions(kp, masked_segmentation_model.gr, data.panel = 2, col = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),paste0(other_col,seg_txp)))), border = ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="UMR",paste0(UMR_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="TEM",paste0(TEM_col,seg_txp),ifelse(as.data.frame(masked_segmentation_model.gr)$segment.mStatus=="GBM",paste0(GBM_col,seg_txp),paste0(other_col,seg_txp)))))

		#kpRect(kp, chr="Chr1", x0=26e6, x1=30e6, y0=0.1, y1=0.9, col="white", border="black", data.panel = 2)
		#kpText(kp, chr="Chr1", x=28e6, y=0.25, labels="Gene-body-like methylation", data.panel = 2, col="green")	
		#kpText(kp, chr="Chr1", x=28e6, y=0.5, labels="TE-like methylation", data.panel = 2, col="red")	
		#kpText(kp, chr="Chr1", x=28e6, y=0.75, labels="Unmethylated", data.panel = 2, col="blue")	
	
		# Plot GBM genes:
		#GBM_gene_col = "#CCCC00"
		#GBM_gene_txp = "88"
		#kpPlotRegions(kp, gene_ranges, data.panel=1, r0=0.1, r1=0.5, col=paste0(GBM_gene_col,GBM_gene_txp), border=paste0(GBM_gene_col,GBM_gene_txp))
		#kpPlotMarkers(kp, data=gene_ranges, labels=gff.genes$gene_ID, data.panel=1, r0=0, r1=0.3, text.orientation="vertical", line.color="#FFFFFF00", clipping=FALSE)

		
		# Plot site details
		line_data = locus_data$parental_consensus
		kpPlotRegions(kp, locus_data.gr, data.panel=2, r0=0, r1=0.15, col=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp), border=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp))
		line_data = locus_data$Line29Gen30
		kpPlotRegions(kp, locus_data.gr, data.panel=2, r0=0.25, r1=0.4, col=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp), border=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp))
		line_data = locus_data$Line49Gen30
		kpPlotRegions(kp, locus_data.gr, data.panel=2, r0=0.45, r1=0.6, col=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp), border=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp))
		line_data = locus_data$Line59Gen30
		kpPlotRegions(kp, locus_data.gr, data.panel=2, r0=0.65, r1=0.8, col=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp), border=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp))
		line_data = locus_data$Line119Gen30
		kpPlotRegions(kp, locus_data.gr, data.panel=2, r0=0.85, r1=1, col=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp), border=paste0(ifelse(line_data == "M", meth_col, ifelse(line_data == "U", unmeth_col, na_col)),var_txp))

		# Add axis labels to site details
		# Axes on data panels seem to fail for zoomed regions
		#kpAxis(kp, r0=0, r1=0.15, tick.pos = c(0.5), labels = c("Parental"), col="#000000", cex=0.5, data.panel=2)


		# Offset in numbers of bp/CG sites from first site at which to put labels
		line_label_offset = -3
		
		# Offset in numbers of bp/CG sites from last site at which to put legend labels
		line_label_offset2 = 0.6

		# Add status label
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(min(rownames(locus_data)))+line_label_offset-1, y=0.05, labels="Change in mCG", cex=0.5, data.panel = 1, col="black", clipping=FALSE)
		# Add status legend
		kpLines(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=rep(as.numeric(max(rownames(locus_data)))-1,2), y=c(0.15,0.25), col=na_col, clipping=FALSE, data.panel = 1)
		kpLines(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=rep(as.numeric(max(rownames(locus_data)))-1,2), y=c(0.27,0.37), col=UM_once_col, clipping=FALSE, data.panel = 1)
		kpLines(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=rep(as.numeric(max(rownames(locus_data)))-1,2), y=c(0.39,0.49), col=MU_once_col, clipping=FALSE, data.panel = 1)
		kpLines(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=rep(as.numeric(max(rownames(locus_data)))-1,2), y=c(0.51,0.61), col=M_no_change_col, clipping=FALSE, data.panel = 1)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_label_offset2, y=0.2, labels="No data", cex=0.5, data.panel = 1, col=na_col, clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_label_offset2, y=0.32, labels="U->M", cex=0.5, data.panel = 1, col=UM_once_col, clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_label_offset2, y=0.44, labels="M->U", cex=0.5, data.panel = 1, col=MU_once_col, clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_label_offset2, y=0.56, labels="Unchanged", cex=0.5, data.panel = 1, col=M_no_change_col, clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_label_offset2, y=0.7, labels="Change in mCG:", cex=0.5, data.panel = 1, col="black", clipping=FALSE)
		

		# Add line labels
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(min(rownames(locus_data)))+line_label_offset-1, y=0, labels="mCG in Line:", cex=0.5, data.panel = 2, col="black", clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(min(rownames(locus_data)))+line_label_offset, y=0.1, labels="Parental", cex=0.5, data.panel = 2, col="black", clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(min(rownames(locus_data)))+line_label_offset, y=.35, labels="Line29Gen30", cex=0.5, data.panel = 2, col="black", clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(min(rownames(locus_data)))+line_label_offset, y=.55, labels="Line49Gen30", cex=0.5, data.panel = 2, col="black", clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(min(rownames(locus_data)))+line_label_offset, y=.75, labels="Line59Gen30", cex=0.5, data.panel = 2, col="black", clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(min(rownames(locus_data)))+line_label_offset, y=.95, labels="Line119Gen30", cex=0.5, data.panel = 2, col="black", clipping=FALSE)

		# Add line status legend
		line_legend_offset = -7
		#kpLines(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=rep(as.numeric(max(rownames(locus_data)))+line_legend_offset-1,2), y=c(0.15,0.25), col=na_col, clipping=FALSE, data.panel = 1)
		kpLines(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=rep(as.numeric(max(rownames(locus_data)))+line_legend_offset-1.5,2), y=c(0.27,0.37), col=na_col, clipping=FALSE, data.panel = 1)
		kpLines(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=rep(as.numeric(max(rownames(locus_data)))+line_legend_offset-1.5,2), y=c(0.39,0.49), col=UM_once_col, clipping=FALSE, data.panel = 1)
		kpLines(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=rep(as.numeric(max(rownames(locus_data)))+line_legend_offset-1.5,2), y=c(0.51,0.61), col=MU_once_col, clipping=FALSE, data.panel = 1)
		#kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_legend_offset+line_label_offset2, y=0.2, labels="No data", cex=0.5, data.panel = 1, col=na_col, clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_legend_offset+line_label_offset2, y=0.32, labels="No data", cex=0.5, data.panel = 1, col=na_col, clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_legend_offset+line_label_offset2, y=0.44, labels="Methylated", cex=0.5, data.panel = 1, col=UM_once_col, clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_legend_offset+line_label_offset2, y=0.56, labels="Unmethylated", cex=0.5, data.panel = 1, col=MU_once_col, clipping=FALSE)
		kpText(kp, chr=mixedCaseChr(locus_data$Chromosome[1]), x=as.numeric(max(rownames(locus_data)))+line_legend_offset+line_label_offset2, y=0.7, labels="mCG in Line:", cex=0.5, data.panel = 1, col="black", clipping=FALSE)
		
	}
	dev.off()
	
	# Revised Fig. Xb:
	# Showing the pattern of parental and offspring methylation at sample loci plotted using ggplot rather than karyoploteR
	pdf(paste0(project_id,"_",meth_context,"_fig_Xb.pdf"), width=9, height=1.2)

	# Query the examples based on changes and coverage
	#genes_of_interest = gff.genes[(gff.genes$variable_count>6) & (gff.genes$m_to_u_count>=2) & (gff.genes$u_to_m_count>=2) & ((gff.genes$all_M_count+gff.genes$all_U_count+gff.genes$variable_count)/gff.genes$CG_sites_count > 0.8),]

	# Query the examples based on a gene list
	genes_of_interest = gff.genes[(gff.genes$gene_ID %in% c("AT5G37190","AT5G40450","AT5G52280","AT1G69070","AT2G35030","AT3G22790","AT4G14760","AT5G39680","AT5G65560","AT1G77580")) ,]
	plot_starts = c()
	plot_ends = c()
	
	#genes_of_interest = gff.genes[(gff.genes$gene_ID %in% c("AT5G52280")) ,]

	
	locus_chr = genes_of_interest$Chromosome
	locus_start = genes_of_interest$V4
	locus_end = genes_of_interest$V5
	locus_name = genes_of_interest$gene_ID

	MU_once_col = 'blue'
	MU_many_col = 'blue'
	UM_once_col = 'red'
	UM_many_col = 'red'
	M_no_change_col = 'green'
	U_no_change_col = 'green'
	na_col = 'grey'
	na_alpha = 0.3
		
	colour_vector = structure(c("blue","green","grey","red"))
		
	for (example_locus in 1:length(locus_name)) {
		locus_data=cbind.data.frame(all_reps_meth_status, parental_consensus, num_lines_m_to_u_from_parentals, num_lines_u_to_m_from_parentals)[(all_reps_meth_status$Chromosome==toupper(locus_chr[example_locus])) & (all_reps_meth_status$Locus>=locus_start[example_locus]) & (all_reps_meth_status$Locus<=locus_end[example_locus]),]
	
		locus_data$Line1Gen3 = as.character(locus_data$Line1Gen3)
		locus_data$Line19Gen3 = as.character(locus_data$Line19Gen3)
		locus_data$Line29Gen30 = as.character(locus_data$Line29Gen30)
		locus_data$Line12Gen3 = as.character(locus_data$Line12Gen3)
		locus_data$Line49Gen30 = as.character(locus_data$Line49Gen30)
		locus_data$Line59Gen30 = as.character(locus_data$Line59Gen30)
		locus_data$Line119Gen30 = as.character(locus_data$Line119Gen30)
		locus_data$parental_consensus = as.character(locus_data$parental_consensus)

		# D, I, O, P are all valid status codes for sites not classified as M or U.  Replace all these with Z so they can be cleanly represented in the plot
		locus_data[locus_data$Line1Gen3 %in% c("D","I","O","P"),"Line1Gen3"] = "Z"
		locus_data[locus_data$Line19Gen3 %in% c("D","I","O","P"),"Line19Gen3"] = "Z"
		locus_data[locus_data$Line29Gen30 %in% c("D","I","O","P"),"Line29Gen30"] = "Z"
		locus_data[locus_data$Line12Gen3 %in% c("D","I","O","P"),"Line12Gen3"] = "Z"
		locus_data[locus_data$Line49Gen30 %in% c("D","I","O","P"),"Line49Gen30"] = "Z"
		locus_data[locus_data$Line59Gen30 %in% c("D","I","O","P"),"Line59Gen30"] = "Z"
		locus_data[locus_data$Line119Gen30 %in% c("D","I","O","P"),"Line119Gen30"] = "Z"
		locus_data[locus_data$parental_consensus %in% c("D","I","O","P"),"parental_consensus"] = "Z"
	
		cat(paste0(nrow(locus_data),"\n"))
		var_col = ifelse(locus_data$num_lines_m_to_u_from_parentals==1,MU_once_col,ifelse(locus_data$num_lines_m_to_u_from_parentals>1,MU_many_col,ifelse(locus_data$num_lines_u_to_m_from_parentals==1,UM_once_col,ifelse(locus_data$num_lines_u_to_m_from_parentals>1,UM_many_col,ifelse(locus_data$parental_consensus=="M", M_no_change_col, ifelse(locus_data$parental_consensus=="U", U_no_change_col, na_col))))))

		locus_plot = ggplot() + geom_point(data = cbind.data.frame("xpos" = seq(from=0, to=nrow(locus_data)-1), "ypos" = rep(7,nrow(locus_data)), "status" = rep("A",nrow(locus_data))), aes(x=xpos, y=ypos, shape = status, colour=var_col)) + scale_color_manual("Change in mCG",values=colour_vector, labels=c("mCG->uCG","Unchanged","No data","uCG->mCG"), guide="none") + scale_shape_manual("Site mCG by line",values = c(15, 16, 21, 21), guide = "none") + scale_alpha_continuous(NULL, guide="none")

		#cat("1\n")
		#print(locus_plot)
		line_data = locus_data$parental_consensus
#		locus_plot = locus_plot + geom_point(data = cbind.data.frame("xpos" = seq(from=0, to=length(line_data)-1), "ypos" = rep(5,length(line_data)), "status" = line_data), aes(x=xpos, y=ypos, shape = status, alpha=ifelse(status %in% c("M","U"),1,na_alpha)), fill='gray')
		locus_plot = locus_plot + geom_point(data = cbind.data.frame("xpos" = seq(from=0, to=length(line_data)-1), "ypos" = rep(5,length(line_data)), "status" = line_data), aes(x=xpos, y=ypos, shape = status), alpha=ifelse(line_data %in% c("M","U"),1,na_alpha), fill=ifelse(line_data=="A",'gray','white'))
		line_data = locus_data$Line29Gen30
		locus_plot = locus_plot + geom_point(data = cbind.data.frame("xpos" = seq(from=0, to=length(line_data)-1), "ypos" = rep(3,length(line_data)), "status" = line_data), aes(x=xpos, y=ypos, shape = status), alpha=ifelse(line_data %in% c("M","U"),1,na_alpha), fill=ifelse(line_data=="A",'gray','white'))
		line_data = locus_data$Line49Gen30
		locus_plot = locus_plot + geom_point(data = cbind.data.frame("xpos" = seq(from=0, to=length(line_data)-1), "ypos" = rep(2,length(line_data)), "status" = line_data), aes(x=xpos, y=ypos, shape = status), alpha=ifelse(line_data %in% c("M","U"),1,na_alpha), fill=ifelse(line_data=="A",'gray','white'))
		line_data = locus_data$Line59Gen30
		locus_plot = locus_plot + geom_point(data = cbind.data.frame("xpos" = seq(from=0, to=length(line_data)-1), "ypos" = rep(1,length(line_data)), "status" = line_data), aes(x=xpos, y=ypos, shape = status), alpha=ifelse(line_data %in% c("M","U"),1,na_alpha), fill=ifelse(line_data=="A",'gray','white'))
		line_data = locus_data$Line119Gen30
		locus_plot = locus_plot + geom_point(data = cbind.data.frame("xpos" = seq(from=0, to=length(line_data)-1), "ypos" = rep(0,length(line_data)), "status" = line_data), aes(x=xpos, y=ypos, shape = status), alpha=ifelse(line_data %in% c("M","U"),1,na_alpha), fill=ifelse(line_data=="A",'gray','white')) 
		#cat("2\n")
		#print(locus_plot)
		# Add in the final formatting elements
		locus_plot = locus_plot + theme_minimal() + labs(x=NULL, y=NULL) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_y_continuous(minor_breaks = c(0,1,2,3,5,7), breaks = c(0,1,2,3,5,7), labels=c("Line119Gen30", "Line59Gen30", "Line49Gen30", "Line29Gen30", "Parental", "Change in mCG")) +scale_x_continuous(breaks=seq(from=0, to=length(line_data)-1), minor_breaks = seq(from=0, to=length(line_data)-1), limits = c(0,100)) +ggtitle(locus_name[example_locus])
		#+ scale_colour_manual(name="Error Bars",values=cols) + scale_fill_manual(name="Bar",values=cols)
		#cat("3\n")
		
		print(locus_plot)
	}
	dev.off()

	# We use ggpubr to arrange a number of ggplots on the page
	#library(ggpubr)

	# Fig Xc - partitioning of variable sites between methylation segments and loci annotations
	pdf(paste0(project_id,"_",meth_context,"_fig_Xc.pdf"), width=pdf_plot_width, height=pdf_plot_height)
	print(ggplot(segment_annotation_summary, aes(x=segment_type, y=variable_sites)) + geom_bar(aes(fill=annotation_type), stat="identity") + labs(x=paste0("Segments classified by methylation pattern"), y="Number of mCG-variable sites", fill="Annotated loci") + facet_grid(~segment_type, scales="free") + theme_minimal() +ggtitle("Fig. Xc Sites of spontaneous mCG change occur predominantly in gbM genes"))

	dev.off()

	# New version of Xc
	# Plot distribution of variable CG sites among different methylation-type segments
	
	z = segmentation_model_site_summary[segmentation_model_site_summary$mStatus!="GBM-like",]
	levels(z$mStatus) = c("gbM segments","gbM-like segments","TEM segments","UM segments")
	levels(z$Methylation) = c("All CG sites","mCG sites","uCG sites","sites where\nmCG varies","mCG->uCG","uCG->mCG")
	#z$Methylation = as.character(z$Methylation)

	pdf(paste0(project_id,"_",meth_context,"_fig_Xc.pdf"), width=pdf_plot_width, height=pdf_plot_height)
	print(ggplot() + geom_bar(data=ddply(z, "Methylation", mutate, percent_Sites=Sites/sum(Sites) * 100), aes(x=Methylation, y=percent_Sites, fill=mStatus), stat="identity") + labs(x="Distribution of invariable and variable CG sites by methylation segment", y="Percentage of Sites", fill="Segment Methylation Status\n") + theme_minimal())
	dev.off()
	
	# New version of Figure Xd
	# Identify sites which have a clear call in parental consensus and all offspring lines
	
	
	# Report on sites which have a transition in each direction, as a proportion of those which could potentially have made the same transition
	
	#locus_data=cbind.data.frame(all_reps_meth_status, parental_consensus, num_lines_m_to_u_from_parentals, num_lines_u_to_m_from_parentals)[(all_reps_meth_status$Chromosome==toupper(locus_chr[example_locus])) & (all_reps_meth_status$Locus>=locus_start[example_locus]) & (all_reps_meth_status$Locus<=locus_end[example_locus]),]
	
	#if (plot_scale == "bp") {
		# This version turns the locus_data into a genomic range on bp scale
	#	locus_data.gr = makeGRangesFromDataFrame(df=cbind.data.frame("Chromosome" = mixedCaseChr(locus_data$Chromosome), "start_site" = locus_data$Locus, "end_site" = locus_data$Locus), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
	#} else {
		# This version turns the locus_data into a genomic range scaled per CG site
	#	locus_data.gr = makeGRangesFromDataFrame(df=cbind.data.frame("Chromosome" = mixedCaseChr(locus_data$Chromosome), "start_site" = rownames(locus_data), "end_site" = rownames(locus_data)), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
	#}


	# This version is restricted to sites which have a clean call in all samples and only change in one gen30 line
	#z = rbind.data.frame(cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(csm_sites)/sum(csm_sites+all_M_sites))),"metric"=rep("-Fraction of mCG sites losing methylation",4)),cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(csu_sites)/sum(csu_sites+all_U_sites))),"metric"=rep("+Fraction of uCG sites gaining methylation",4)))
	
	# This version reflects all sites, regardless of missing data and number of lines showing change
	z = rbind.data.frame(cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(m_to_u_sites)/sum(m_to_u_sites+all_M_sites))),"metric"=rep("-Fraction of mCG sites losing methylation",4)),cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(u_to_m_sites)/sum(u_to_m_sites+all_U_sites))),"metric"=rep("+Fraction of uCG sites gaining methylation",4)))

	z = z[z$segment_type != "GBM-like segments",]
	z[z$segment_type == "GBM segments",]$segment_type = "gbM segments"
	#z[z$segment_type == "GBM-like segments",]$segment_type = "gbM-like segments"
	z[z$segment_type == "UMR segments",]$segment_type = "UM segments"
		
	# NB. plotting this as a PDF reintroduces the legend and messes up the facet y axis limits, so need to plot it and save the PDF separately - 4"x5.2" landscape
	pdf(paste0(project_id,"_",meth_context,"_fig_Xd.pdf"), width=pdf_plot_width, height=pdf_plot_height)
	#print(ggplot(z) + geom_bar(aes(x=segment_type, y=score, fill=segment_type), stat="identity") + facet_grid(metric ~ ., scales="fixed") + labs(x="", y="") + theme_minimal())
	
	#plot1 = ggplot(z[4:6,]) + geom_bar(aes(x=factor(segment_type, levels=c("UM segments", "gbM segments", "TEM segments")), y=score, fill=segment_type), stat="identity") + facet_grid(metric ~ .) + labs(x=NULL, y=NULL) + ylim(0,0.3) + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme(legend.position="none")
	#plot2 = ggplot(z[1:3,]) + geom_bar(aes(x=factor(segment_type, levels=c("UM segments", "gbM segments", "TEM segments")), y=-score, fill=segment_type), stat="identity") + facet_grid(metric ~ .) + labs(x=NULL, y=NULL) + ylim(-0.3,0) + theme_minimal() + theme(legend.position="none")
	#gridExtra::grid.arrange(plot1,plot2,nrow=2)

	
	# This version is restricted to sites which have a clean call in all samples and only change in one gen30 line
	#z = rbind.data.frame(cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(csm_sites)/sum(csm_sites+all_M_sites))),"metric"=rep("Sites losing methylation",4)),cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(csu_sites)/sum(csu_sites+all_U_sites))),"metric"=rep("Sites gaining methylation",4)))
	
	# This version reflects all sites, regardless of missing data and number of lines showing change
	z = rbind.data.frame(cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(m_to_u_sites)/sum(m_to_u_sites+all_M_sites))),"metric"=rep("Sites losing methylation",4)),cbind.data.frame(as.data.frame(segment_annotation_summary %>% group_by(segment_type) %>% summarise(score = sum(u_to_m_sites)/sum(u_to_m_sites+all_U_sites))),"metric"=rep("Sites gaining methylation",4)))

	z = z[z$segment_type != "GBM-like segments",]
	z[z$segment_type == "GBM segments",]$segment_type = "gbM segments"
	#z[z$segment_type == "GBM-like segments",]$segment_type = "gbM-like segments"
	z[z$segment_type == "UMR segments",]$segment_type = "UM segments"

	plot1 = ggplot(z[,]) + geom_bar(aes(x=factor(segment_type, levels=c("UM segments", "gbM segments", "TEM segments")), y=score, fill=segment_type), stat="identity") + facet_grid(. ~ metric) + labs(x=NULL, y="Fraction of sites with given change") + theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme(legend.position="none")
	#plot2 = ggplot(z[1:3,]) + geom_bar(aes(x=factor(segment_type, levels=c("UM segments", "gbM segments", "TEM segments")), y=score, fill=segment_type), stat="identity") + facet_grid(metric ~ .) + labs(x=NULL, y=NULL) + ylim(0,0.3) + theme_minimal() + theme(legend.position="none")
	print(plot1)

	dev.off()

	# What are all the methylated sites in UM segments?
	cat(sum(segment_annotation_summary[segment_annotation_summary$segment_type == "UMR segments","all_M_sites"]) + sum(segment_annotation_summary[segment_annotation_summary$segment_type == "UMR segments","m_to_u_sites"]))
	# 65
	
	# Checking why there are m->u sites in UM segments 
	# e.g. Show the parental read data for the M->U sites which overlap with the Unannotated regions
	meth_parents_CG.gr[queryHits(findOverlaps(meth_parents_CG.gr, m_to_u_call_ranges[subjectHits(findOverlaps(U_Unann_overlaps, m_to_u_call_ranges))]))]
#	GRanges object with 13 ranges and 2 metadata columns:
#       seqnames    ranges strand |         T         M
#          <Rle> <IRanges>  <Rle> | <numeric> <numeric>
#   [1]     Chr1  12971502      * |       181       130
#   [2]     Chr1  13403672      * |       434       284
#   [3]     Chr1  24913335      * |       259       196
#   [4]     Chr2   8952993      * |       508       427
#   [5]     Chr2   8952997      * |       512       437
#   ...      ...       ...    ... .       ...       ...
#   [9]     Chr3   8585079      * |       441       228
#  [10]     Chr3   8585095      * |       453       288
#  [11]     Chr3  15171113      * |       205       112
#  [12]     Chr5   3958299      * |       389       245
#  [13]     Chr5   7728519      * |       160       117
#	
#	All these loci are methylated in well below 50% of reads, but, in all likelihood 2 out of 3 parental lines were above the threshold to be called as methylated, whereas the consensus at the site was too low to break the valid unmethylated region surrounding the locus
	
	# Fig Xe - tendency of individual gbM loci to lose and gain methylation at independent sites
	pdf(paste0(project_id,"_",meth_context,"_fig_Xe.pdf"), width=pdf_plot_width, height=pdf_plot_height)
	print(ggplot(rbind.data.frame(cbind.data.frame("count"=gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count+gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count,"metric"=rep("Number of sites in gene\nwhere mCG changes",nrow(gff.genes[gff.genes$m_class=="Gene-body Methylated",]))), cbind.data.frame("count"=abs(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count-gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count),"metric"=rep("Overall net change in mCG in gene",nrow(gff.genes[gff.genes$m_class=="Gene-body Methylated",])))), aes(x=as.numeric(count), fill=metric)) + geom_histogram(position="dodge", binwidth=1) + theme_minimal() +xlab("mCG changes per gene") + ylab("Number of annotated gbM genes") +xlim(0,32) +ggtitle("Fig Xe Gains and losses of mCG in individual genes tend towards balance"))
	
	
	#ggarrange(figXa, figXb, figXc, labels = c("A","B","C"), ncol=2, nrow=2)
	dev.off()
	
	
	# Fig Ya
	# Script 5-methylation_calls_metaplot_mCWG_around_mCG_var_loci_2018-05-10.R executes this plot
	
	# Fig Yb
	
	
	
	
	
	# First pass results correlating annotation with segmentation class:
	#GBM segments	GBM genes	0.935591851	49971
	#GBM segments	Het. genes	0.005780837	276
	#GBM segments	Unmeth. genes	0.003273897	72
	#GBM segments	Other genes	0.001123589	1
	#GBM segments	GBM TEs	0.018209808	696
	#GBM segments	Het. TEs	0.011099579	262
	#GBM segments	Unmeth. TEs	0.001311192	21
	#GBM segments	Other TEs	0.005625418	11
	#GBM segments	Unannotated	0.026934275	1086
	#TEM segments	GBM genes	0.02483327	500
	#TEM segments	Het. genes	0.056177657	396
	#TEM segments	Unmeth. genes	0.001030929	7
	#TEM segments	Other genes	0.002353038	0
	#TEM segments	GBM TEs	0.026433546	210
	#TEM segments	Het. TEs	0.652372367	1065
	#TEM segments	Unmeth. TEs	0.000946794	11
	#TEM segments	Other TEs	0.000946794	30
	#TEM segments	Unannotated	0.233613784	1527
	#UMR segments	GBM genes	0.375435916	9869
	#UMR segments	Het. genes	0.00105929	131
	#UMR segments	Unmeth. genes	0.255903738	1837
	#UMR segments	Other genes	0.002442081	1
	#UMR segments	GBM TEs	0.015618335	501
	#UMR segments	Het. TEs	0.004185256	344
	#UMR segments	Unmeth. TEs	0.025360416	313
	#UMR segments	Other TEs	0.01564349	43
	#UMR segments	Unannotated	0.322267066	2344

	# sample of each type of overlap was reviewed visually to diagnose the likely cause of anomalous annotations/mis-segmentations:
	
	# Unmethylated segments are generally correct.
	
	# TEM segments are generally annotated correctly.  16 of the longer overlaps between unmethylated genes and TEM segments are mis-segmentations.
	# TEs annotated as GBM, but overlapping TEM segments should be annotated as het. TEs
	# TEs annotated as 'other', but overlapping TEM segments should be annotated as het. TEs
	# Some genes annotated as GBM but overlapping TEM segment should be annotated as het. or part het.
	
	#<5% of genome which is annotated as GBM segments is misannotated or mis-segmented (Het. TEs, Unmeth. TEs, Other TEs, Unannotated) - we will reannotate these segments if they don't overlap with a gene body, and if they do, we will split the segments and reannotate the part extending past the limit of the gene body.
	
	# GBM segments overlapping TEs, but not genes, should be reannotated as TEM, but we will call GBM-like for now
	# GBM segments overlapping no genes or TEs are mostly TEM, but we can call 'GBM-like' for now
	# Remaining GBM segments overlapping an unannotated region are likely segmentation anomalies
	# Het. genes, if overlapping a GBP segment, are likely GBM genes
	# 'other' genes overlapping a GBM segment should be annotated as GBM genes
	# Unmeth genes with GBM segment overlap should be annotated as GBM genes
	
	
	# After reannotating segments based on overlap with gene model annotations, we have this:
	#GBM segments	GBM genes	0.991971754	49971
	#GBM segments	Het. genes	0.006015511	276
	#GBM segments	Unmeth. genes	0.003471186	72
	#GBM segments	Other genes	0.000975244	1
	#GBM segments	GBM TEs	0.003729965	214
	#GBM segments	Het. TEs	0.001578532	67
	#GBM segments	Unmeth. TEs	0.000269474	2
	#GBM segments	Other TEs	0.001181848	2
	#GBM segments	Unannotated	0	0
	#GBM-like segments	GBM genes	0	0
	#GBM-like segments	Het. genes	0	0
	#GBM-like segments	Unmeth. genes	0	0
	#GBM-like segments	Other genes	0	0
	#GBM-like segments	GBM TEs	0.272720886	480
	#GBM-like segments	Het. TEs	0.175798179	189
	#GBM-like segments	Unmeth. TEs	0.019015788	18
	#GBM-like segments	Other TEs	0.08032981	9
	#GBM-like segments	Unannotated	0.456986139	1054
	#TEM segments	GBM genes	0.024832672	500
	#TEM segments	Het. genes	0.056019107	396
	#TEM segments	Unmeth. genes	0.001030904	7
	#TEM segments	Other genes	0.002194251	0
	#TEM segments	GBM TEs	0.026525792	212
	#TEM segments	Het. TEs	0.652491568	1066
	#TEM segments	Unmeth. TEs	0.000967936	11
	#TEM segments	Other TEs	0.000967936	30
	#TEM segments	Unannotated	0.233686608	1531
	#UMR segments	GBM genes	0.375191562	9869
	#UMR segments	Het. genes	0.0010586	131
	#UMR segments	Unmeth. genes	0.255737182	1837
	#UMR segments	Other genes	0.002440492	1
	#UMR segments	GBM TEs	0.015622706	501
	#UMR segments	Het. TEs	0.004200362	349
	#UMR segments	Unmeth. TEs	0.025347898	314
	#UMR segments	Other TEs	0.015668219	43
	#UMR segments	Unannotated	0.322637765	2372

	
	
	# At this point, segmentation_model.gr == masked_segmentation_model.gr
	
	# We will make amendments to segmentation_model.gr which we could potentially rollback by setting segmentation_model.gr=masked_segmentation_model.gr
	
	# Capture initial state of segments from model, so that x_x_x_hits will continue to work properly
	GBM_segments.gr = segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]
	TEM_segments.gr = segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]
	UMR_segments.gr = segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]
	length(GBM_segments.gr)
	length(TEM_segments.gr)
	length(UMR_segments.gr)
	# 18730
	# 7447
	# 25648
	
	# GBM segments overlapping an annotated TE:
	G_x_T_segments = unique(c(subjectHits(G_G_T_hits),subjectHits(G_H_T_hits),subjectHits(G_U_T_hits),subjectHits(G_O_T_hits)))
	length(G_x_T_segments)
	# 1525 segments from the ones currently annotated as GBM

	# GBM segments overlapping an annotated TE, but no annotated genes:
	# (find which of the G_x_T segments overlaps with any of the gene loci ranges, and take the set difference to exclude these from the G_x_T segments)
	G_x_T_no_gene_segments.gr = setdiff(GBM_segments.gr[G_x_T_segments], GBM_segments.gr[G_x_T_segments][unique(subjectHits(findOverlaps(c(GBM_gene_loci.gr, Het_gene_loci.gr, UM_gene_loci.gr, other_gene_loci.gr),GBM_segments.gr[G_x_T_segments])))])
	length(G_x_T_no_gene_segments.gr)
	# 1014 segments
	
	# Set such segments to be GBM-like rather than GBM
	segmentation_model.gr[unique(subjectHits(findOverlaps(G_x_T_no_gene_segments.gr, segmentation_model.gr)))]$segment.mStatus = rep("GBM-like", length(unique(subjectHits(findOverlaps(G_x_T_no_gene_segments.gr, segmentation_model.gr)))))
	
	
	# GBM segments overlapping an annotated gene or TE:
	G_x_x_segments = unique(c(G_x_T_segments, subjectHits(G_G_G_hits),subjectHits(G_H_G_hits),subjectHits(G_U_G_hits),subjectHits(G_O_G_hits)))
	length(G_x_x_segments)
	# 18015 segments from the ones currently annotated as GBM

	# GBM segments not overlapping an annotated TE, nor an annotated gene:
	# (find which non-G_x_G segments overlaps with any of the gene loci ranges, and take the set difference to exclude these from the G_x_T segments)
	G_x_x_no_annotation_segments.gr = setdiff(GBM_segments.gr, GBM_segments.gr[G_x_x_segments])
	length(G_x_x_no_annotation_segments.gr)
	# 714 segments
	
	# Set such segments to be GBM-like rather than GBM
	segmentation_model.gr[unique(subjectHits(findOverlaps(G_x_x_no_annotation_segments.gr, segmentation_model.gr)))]$segment.mStatus = rep("GBM-like", length(unique(subjectHits(findOverlaps(G_x_x_no_annotation_segments.gr, segmentation_model.gr)))))
	
	# Recapture modified state of segments from model, so that segmentation model can be rebuilt after reduce operation
	GBM_segments.gr = segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]
	TEM_segments.gr = segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="TEM"]
	UMR_segments.gr = segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="UMR"]
	GBM_like_segments.gr = segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM-like"]
	length(GBM_segments.gr)
	length(TEM_segments.gr)
	length(UMR_segments.gr)
	length(GBM_like_segments.gr)
	# 17001
	# 7447
	# 25648
	# 1729

	# We still have some GBM segments which overflow past the bounds of the annotated genes.  We will trim these to end at the end of the gene body, and reannotate accordingly.  We will also annotate the remaining segments which lie outside gene bodies.
	
	G_x_G_hits = findOverlaps(reduce(c(GBM_gene_loci.gr, Het_gene_loci.gr, UM_gene_loci.gr, other_gene_loci.gr)), reduce(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))
	# Define set of segments representing the overlap between GBM segments and gene models
	true_GBM_segments.gr = setdiff(reduce(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])[subjectHits(G_x_G_hits)],setdiff(reduce(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])[subjectHits(G_x_G_hits)],reduce(c(GBM_gene_loci.gr, Het_gene_loci.gr, UM_gene_loci.gr, other_gene_loci.gr))[queryHits(G_x_G_hits)]))
	length(true_GBM_segments.gr)
	# 17285
	# Define set of segments representing the parts of GBM segments which are external of gene models
	ext_GBM_segments.gr = setdiff(reduce(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]), true_GBM_segments.gr) 
	length(ext_GBM_segments.gr)
	# 892
	
	
	# Read the parental methylome data in as Granges objects, to use to annotate new segments
	meth_parents_CG.gr <- readMethylome(FileName=paste0(project_id,"_",meth_context,"_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)
	meth_parents_CHG.gr <- readMethylome(FileName=paste0(project_id,"_CHG_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)
	meth_parents_CHH.gr <- readMethylome(FileName=paste0(project_id,"_CHH_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)

	# First the 'real_GBM' segments
	# Set aside any which overlap no CG sites
	set_aside_no_CG.gr = setdiff(true_GBM_segments.gr, true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, true_GBM_segments.gr)))])
	cat(paste0("true_GBM_segments.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 70
	true_GBM_segments.gr = true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, true_GBM_segments.gr)))]
	cat(paste0("true_GBM_segments.gr contains ",length(true_GBM_segments.gr)," segments after trimming segments with no CG sites.\n"))
	# 17215
	# Annotate the remainder with their mean mCG values
	values(true_GBM_segments.gr) = cbind(values(true_GBM_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=true_GBM_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set aside any of the remainder which overlap no CHG sites 
	set_aside_no_CHG.gr = setdiff(true_GBM_segments.gr, true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, true_GBM_segments.gr)))])
	cat(paste0("true_GBM_segments.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 389
	true_GBM_segments.gr = true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, true_GBM_segments.gr)))]
	cat(paste0("true_GBM_segments.gr contains ",length(true_GBM_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 16826
	# Annotate the remainder with their mean mCHG values
	values(true_GBM_segments.gr) = cbind(values(true_GBM_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=true_GBM_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set aside any of the remainder which overlap no CHH sites 
	set_aside_no_CHH.gr = setdiff(true_GBM_segments.gr, true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, true_GBM_segments.gr)))])
	cat(paste0("true_GBM_segments.gr contains ",length(set_aside_no_CHH.gr)," segments which contain no CHH sites.\n"))
	# 10
	true_GBM_segments.gr = true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, true_GBM_segments.gr)))]
	cat(paste0("true_GBM_segments.gr contains ",length(true_GBM_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 16816
	# Annotate the remainder with their mean mCHH values
	values(true_GBM_segments.gr) = cbind(values(true_GBM_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=true_GBM_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	true_GBM_segments.gr$segment.mStatus = rep("GBM", length(true_GBM_segments.gr))
			
	# Function to return NAs for cases where segmentation model doen's overlap any sites in the given context
	NA_by_segment <- function (segment_model)
	{
		segment_count = length(segment_model)
		nSites.segmentation = rep(as.integer(NA),segment_count)
		nSites = rep(as.integer(NA),segment_count)
		T = array(NA,segment_count)
		M = array(NA,segment_count)
		pmeth = array(NA,segment_count)
		median.meth = rep(as.numeric(NA),segment_count)
		type = array(NA,segment_count)

		DataFrame(nSites.segmentation, nSites, T, M, pmeth = M/T, median.meth = median.meth, type)
	}
	
	# Define a corresponding convenience function to save an annotated genomic ranges object
	# Function adapted from MethylSeekR package
	saveUMRLMRSegments <- function (segs, GRangesFilename = NULL, TableFilename = NULL) 
	{
		# Replaced nCG with nSites for generality
		if (!is.null(GRangesFilename)) 
			saveRDS(segs, GRangesFilename)
		if (!is.null(TableFilename)) {
			D = data.frame(chr = as.character(seqnames(segs)), start = start(segs), 
				end = end(segs), type = values(segs)$type, nSites.segmentation = values(segs)$nSites.segmentation, 
				nSites.seq = values(segs)$nSites, mean.meth = 100 * values(segs)$pmeth, 
				median.meth = 100 * values(segs)$median.meth)
			write.table(D, file = TableFilename, quote = FALSE, sep = "\t", 
			row.names = FALSE)
		}
	}

	# For now, we will annotate mean methylation values of the set-aside segments as NA, and annotate them all as GBM
	values(set_aside_no_CG.gr) = cbind(NA_by_segment(set_aside_no_CG.gr), NA_by_segment(set_aside_no_CG.gr), NA_by_segment(set_aside_no_CG.gr))
	values(set_aside_no_CHG.gr) = cbind(NA_by_segment(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr))
	values(set_aside_no_CHH.gr) = cbind(NA_by_segment(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr))
	set_aside_segments.gr = c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr)
	set_aside_segments.gr$segment.mStatus = rep("GBM",length(set_aside_segments.gr))
	true_GBM_segments.gr = c(true_GBM_segments.gr, set_aside_segments.gr)

	
	# Now tackle the GBM_leftovers.gr segments, which are not true GBM, but may not be GBM-like either
	# Set aside any which overlap no CG sites
	set_aside_no_CG.gr = setdiff(ext_GBM_segments.gr, ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, ext_GBM_segments.gr)))])
	cat(paste0("ext_GBM_segments.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 242
	ext_GBM_segments.gr = ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, ext_GBM_segments.gr)))]
	cat(paste0("ext_GBM_segments.gr contains ",length(ext_GBM_segments.gr)," segments after trimming segments with no CG sites.\n"))
	# 650
	# Annotate the remainder with their mean mCG values
	values(ext_GBM_segments.gr) = cbind(values(ext_GBM_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=ext_GBM_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set aside any which overlap no CHG sites
	set_aside_no_CHG.gr = setdiff(ext_GBM_segments.gr, ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, ext_GBM_segments.gr)))])
	cat(paste0("ext_GBM_segments.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 47
	ext_GBM_segments.gr = ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, ext_GBM_segments.gr)))]
	cat(paste0("ext_GBM_segments.gr contains ",length(ext_GBM_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 603
	# Annotate the remainder with their mean mCHG values
	values(ext_GBM_segments.gr) = cbind(values(ext_GBM_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=ext_GBM_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set aside any which overlap no CHH sites
	set_aside_no_CHH.gr = setdiff(ext_GBM_segments.gr, ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, ext_GBM_segments.gr)))])
	cat(paste0("ext_GBM_segments.gr contains ",length(set_aside_no_CHH.gr)," segments which contain no CHH sites.\n"))
	# 0
	ext_GBM_segments.gr = ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, ext_GBM_segments.gr)))]
	cat(paste0("ext_GBM_segments.gr contains ",length(ext_GBM_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 603
	# Annotate the remainder with their mean mCHH values
	values(ext_GBM_segments.gr) = cbind(values(ext_GBM_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=ext_GBM_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set the segment class of the ext_GBM_segments.gr based on their methylation means
	ext_GBM_segments.gr@elementMetadata@listData$segment.mStatus = ifelse((cbind(data.frame(values(ext_GBM_segments.gr)@listData),mCHG=values(ext_GBM_segments.gr)@listData$pmeth.1, seg_width=ext_GBM_segments.gr@ranges@width)[,]$mCHG>segment_mCHG_cutoff) | (cbind(data.frame(values(ext_GBM_segments.gr)@listData),mCHH=values(ext_GBM_segments.gr)@listData$pmeth.2, seg_width=ext_GBM_segments.gr@ranges@width)[,]$mCHH>segment_mCHH_cutoff),"TEM",ifelse(cbind(data.frame(values(ext_GBM_segments.gr)@listData),mCHG=values(ext_GBM_segments.gr)@listData$pmeth.1, seg_width=ext_GBM_segments.gr@ranges@width)[,]$pmeth>m.sel,"GBM-like","UMR"))
	### We should then reduce each class of the model, since we have introduced a number of new segments with annotations that already exist in the model, and may be adjacent

	# For now, we will annotate mean methylation values of the set-aside segments as NA, and annotate them all as GBM-ext
	### A better strategy would be to allow them to take on the annotation of their neighbour, then merge them, then reannotate the neighbours
	values(set_aside_no_CG.gr) = cbind(NA_by_segment(set_aside_no_CG.gr), NA_by_segment(set_aside_no_CG.gr), NA_by_segment(set_aside_no_CG.gr))
	values(set_aside_no_CHG.gr) = cbind(NA_by_segment(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr))
	values(set_aside_no_CHH.gr) = cbind(NA_by_segment(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr))
	set_aside_segments.gr = c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr)
	
	#set_aside_segments.gr$segment.mStatus = rep("GBM-ext",length(set_aside_segments.gr))
	### Setting segment.mStatus to "UMR" is almost always the right thing to do here since a genuine GBM segment within a gene body almost invariably results in the gene having a UMR region surrounding it. In fact these set-aside segments are often an indicator that there is a larger unmethylated segment that should expand from the seed of the set-aside segment
	set_aside_segments.gr$segment.mStatus = rep("UMR",length(set_aside_segments.gr))
	
	### We would then ideally need to reduce the UMR segments in the model, and reannotate them, but the numbers are so small as to be inconsequential
	ext_GBM_segments.gr = c(ext_GBM_segments.gr, set_aside_segments.gr)
	
	# We may have now left some GBM segments without any CG sites, or having had their mCG sites shipped out to an ext_GBM_segments.gr segment.
	# Reannotate all GBM sites 
	
	GBM_like_segments.gr$ID=NULL
	GBM_like_segments.gr$CG_site_count=NULL
	GBM_like_segments.gr$variant_count=NULL
	GBM_like_segments.gr$all_M_count=NULL
	GBM_like_segments.gr$all_U_count=NULL
	GBM_like_segments.gr$m_to_u_count=NULL
	GBM_like_segments.gr$u_to_m_count=NULL
	TEM_segments.gr$ID=NULL
	TEM_segments.gr$CG_site_count=NULL
	TEM_segments.gr$variant_count=NULL
	TEM_segments.gr$all_M_count=NULL
	TEM_segments.gr$all_U_count=NULL
	TEM_segments.gr$m_to_u_count=NULL
	TEM_segments.gr$u_to_m_count=NULL
	UMR_segments.gr$ID=NULL
	UMR_segments.gr$CG_site_count=NULL
	UMR_segments.gr$variant_count=NULL
	UMR_segments.gr$all_M_count=NULL
	UMR_segments.gr$all_U_count=NULL
	UMR_segments.gr$m_to_u_count=NULL
	UMR_segments.gr$u_to_m_count=NULL
	
	
	segmentation_model.gr = sort(sortSeqlevels(c(true_GBM_segments.gr, ext_GBM_segments.gr, GBM_like_segments.gr, TEM_segments.gr, UMR_segments.gr)))
	length(segmentation_model.gr)
	# 53001

	write.table(cbind("Chromosome"=as.character(as.data.frame(segmentation_model.gr)$seqnames),"Start"=as.data.frame(segmentation_model.gr)$start,"End"=as.data.frame(segmentation_model.gr)$end,"Type"=as.data.frame(segmentation_model.gr)$segment.mStatus),file=paste0(project_id,"_",meth_context,"_segmentation_model_draft2.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Save a copy in case we want to come back to it later
	saveUMRLMRSegments(segs=segmentation_model.gr, GRangesFilename=paste0(project_id,"_",meth_context,"_segmentation_model_draft2.rds"), TableFilename=paste0(project_id,"_",meth_context,"_segmentation_model_detailed_draft2.tsv"))

	# Check whether there are any gaps in the segmentation
	View(as.data.frame(setdiff(setdiff(genome_loci.gr, c(mask_loci.gr, segment_mask_loci.gr)), segmentation_model.gr)))
	
	# 50 gaps identified - on inspection, all are very short (<600bp) UMR segments, typically in heterochromatic regions, which do not contain any CHG sites.  Often they are in between a TEM and GBM-like segment, and seem enriched for variable sites.
	#Instantiate these as new UMR segments in the model:
	missing_segments.gr = setdiff(setdiff(genome_loci.gr, c(mask_loci.gr, segment_mask_loci.gr)), segmentation_model.gr)
	values(missing_segments.gr) = cbind(NA_by_segment(missing_segments.gr),NA_by_segment(missing_segments.gr),NA_by_segment(missing_segments.gr))
	missing_segments.gr$segment.mStatus = rep("UMR",length(missing_segments.gr))
	segmentation_model.gr = c(segmentation_model.gr, missing_segments.gr)
	
	write.table(cbind("Chromosome"=as.character(as.data.frame(segmentation_model.gr)$seqnames),"Start"=as.data.frame(segmentation_model.gr)$start,"End"=as.data.frame(segmentation_model.gr)$end,"Type"=as.data.frame(segmentation_model.gr)$segment.mStatus),file=paste0(project_id,"_",meth_context,"_segmentation_model_draft3.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Save a copy in case we want to come back to it later
	saveUMRLMRSegments(segs=segmentation_model.gr, GRangesFilename=paste0(project_id,"_",meth_context,"_segmentation_model_draft3.rds"), TableFilename=paste0(project_id,"_",meth_context,"_segmentation_model_detailed_draft3.tsv"))
	
	### THE SEGMENTATION IS FINISHED!  NOW NEED TO GO BACK ABOVE TO REGENERATE REPORTS ON OVERLAPS BETWEEN SEGMENTS AND VARIABLE SITES, SEGMENTS AND ANNOTATED MODELS, AND WHOLE GENOME PLOTS ###
	
	# What proportion of remaining GBM segments lays outside of annotated genes?

	G_x_G_hits=findOverlaps(reduce(c(GBM_gene_loci.gr, Het_gene_loci.gr, UM_gene_loci.gr, other_gene_loci.gr)), segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"])
	G_x_G_overlaps <- pintersect(reduce(c(GBM_gene_loci.gr, Het_gene_loci.gr, UM_gene_loci.gr, other_gene_loci.gr))[queryHits(G_x_G_hits)], segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"][subjectHits(G_x_G_hits)])
	G_x_G_overlap_prop <- sum(width(G_x_G_overlaps)) / sum(width(segmentation_model.gr[as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"]))

	
	
	
	# What is the distribution of CG site coverage across all samples?
	# How does this stack up against genome coverage?
	pdf(paste0(project_id,"_",meth_context,"_site_coverage_density_by_sample.pdf"))
	CG_sample_cov_means <- ddply(coverage_data[coverage_data$Chromosome %in% c("CHR1","CHR2","CHR3","CHR4","CHR5"),][sample(nrow(coverage_data[coverage_data$Chromosome %in% c("CHR1","CHR2","CHR3","CHR4","CHR5"),]), 10000),], "Sample", summarise, grp.mean=mean(Cov_C+Cov_T))
	print(ggplot(coverage_data[coverage_data$Chromosome %in% c("CHR1","CHR2","CHR3","CHR4","CHR5"),][sample(nrow(coverage_data[coverage_data$Chromosome %in% c("CHR1","CHR2","CHR3","CHR4","CHR5"),]), 10000),], aes(x=Cov_C+Cov_T, colour=Sample)) +geom_density() + geom_vline(data=CG_sample_cov_means, aes(xintercept=grp.mean, colour=Sample)) + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) + theme_bw() + annotation_logticks() +ggtitle(paste0(project_id," ",meth_context," site coverage density by sample")))
	dev.off()
	
	
	
	
	
	# Sequeira-Mendes et al, 2014 http://www.plantcell.org/content/26/6/2351/ described 9 distinct chromatin states in the Arabidopsis genome and their topographic relationships.
	# Load the Sequeira-Mendes classifications as a set of genomic ranges for useful stuff later
	# These ranges appear to be contiguous, covering the entire nuclear genome
	
	chromatin_state_ranges = read.table(file=reference_chromatin_states, sep="\t", header=TRUE)
	chromatin_state_ranges$Chrom=paste0("CHR",chromatin_state_ranges$Chrom)
	# Plot distribution of chromatin state range sizes
	pdf(paste0(project_id,"_",meth_context,"_chromatin_state_segment_lengths.pdf"))
	print(ggplot(chromatin_state_ranges[chromatin_state_ranges$To-chromatin_state_ranges$From<10000,], aes(x=To-From)) + geom_histogram(aes(y= ..count..), binwidth=150) + theme_minimal())
	dev.off()
	
	# looks binomial with maximum at ~700nt
	### Fit a model to this histogram
	
	# Turn the chromatin states into GenomicRanges:
	chromatin_state_granges=makeGRangesFromDataFrame(df = chromatin_state_ranges, start.field = "From", end.field = "To",seqnames.field = "Chrom")

	# Find overlaps between chromatin states and CG sites with various patterns of variation
	chromatin_state_ranges$ID=row.names(chromatin_state_ranges)
	CG_chromatin_olaps = findOverlaps(CG_site_ranges, chromatin_state_granges)
	chromatin_state_ranges$CG_site_count=table(chromatin_state_ranges[subjectHits(CG_chromatin_olaps),]$ID)[chromatin_state_ranges$ID]
	chromatin_state_ranges$CG_site_count[is.na(chromatin_state_ranges$CG_site_count)]=0
	variant_chromatin_olaps = findOverlaps(variant_call_ranges, chromatin_state_granges)
	chromatin_state_ranges$variant_count=table(chromatin_state_ranges[subjectHits(variant_chromatin_olaps),]$ID)[chromatin_state_ranges$ID]
	chromatin_state_ranges$variant_count[is.na(chromatin_state_ranges$variant_count)]=0
	all_M_chromatin_olaps = findOverlaps(all_M_ranges, chromatin_state_granges)
	chromatin_state_ranges$all_M_count=table(chromatin_state_ranges[subjectHits(all_M_chromatin_olaps),]$ID)[chromatin_state_ranges$ID]
	chromatin_state_ranges$all_M_count[is.na(chromatin_state_ranges$all_M_count)]=0
	all_U_chromatin_olaps = findOverlaps(all_U_ranges, chromatin_state_granges)
	chromatin_state_ranges$all_U_count=table(chromatin_state_ranges[subjectHits(all_U_chromatin_olaps),]$ID)[chromatin_state_ranges$ID]
	chromatin_state_ranges$all_U_count[is.na(chromatin_state_ranges$all_U_count)]=0
	m_to_u_chromatin_olaps = findOverlaps(m_to_u_call_ranges, chromatin_state_granges)
	chromatin_state_ranges$m_to_u_count=table(chromatin_state_ranges[subjectHits(m_to_u_chromatin_olaps),]$ID)[chromatin_state_ranges$ID]
	chromatin_state_ranges$m_to_u_count[is.na(chromatin_state_ranges$m_to_u_count)]=0
	u_to_m_chromatin_olaps = findOverlaps(u_to_m_call_ranges, chromatin_state_granges)
	chromatin_state_ranges$u_to_m_count=table(chromatin_state_ranges[subjectHits(u_to_m_chromatin_olaps),]$ID)[chromatin_state_ranges$ID]
	chromatin_state_ranges$u_to_m_count[is.na(chromatin_state_ranges$u_to_m_count)]=0

	# Summarise the numbers of variant, all_M and all_U sites by chromatin state:
	
	# Quick barplot version:
	barplot(prop.table(cbind("CG sites"=aggregate(CG_site_count ~ State, data=chromatin_state_ranges, FUN=sum)$CG_site_count,"all_M"=aggregate(all_M_count ~ State, data=chromatin_state_ranges, FUN=sum)$all_M_count,"all_U"=aggregate(all_U_count ~ State, data=chromatin_state_ranges, FUN=sum)$all_U_count, "variable"=aggregate(variant_count ~ State, data=chromatin_state_ranges, FUN=sum)$variant_count,"variableMU"=aggregate(m_to_u_count ~ State, data=chromatin_state_ranges, FUN=sum)$m_to_u_count,"variableUM"=aggregate(u_to_m_count ~ State, data=chromatin_state_ranges, FUN=sum)$u_to_m_count), 2))
	
	# Summarise in a table for nicer plotting:
	chromatin_state_site_summary=data.frame(rbind(cbind("Methylation"="All CG sites","State"=aggregate(CG_site_count ~ State, data=chromatin_state_ranges, FUN=sum)$State,"Sites"=aggregate(CG_site_count ~ State, data=chromatin_state_ranges, FUN=sum)$CG_site_count),cbind("Methylation"="Variable","State"=aggregate(variant_count ~ State, data=chromatin_state_ranges, FUN=sum)$State,"Sites"=aggregate(variant_count ~ State, data=chromatin_state_ranges, FUN=sum)$variant_count),cbind("Methylation"="VariableMU","State"=aggregate(m_to_u_count ~ State, data=chromatin_state_ranges, FUN=sum)$State,"Sites"=aggregate(m_to_u_count ~ State, data=chromatin_state_ranges, FUN=sum)$m_to_u_count),cbind("Methylation"="VariableUM","State"=aggregate(u_to_m_count ~ State, data=chromatin_state_ranges, FUN=sum)$State,"Sites"=aggregate(u_to_m_count ~ State, data=chromatin_state_ranges, FUN=sum)$u_to_m_count),cbind("Methylation"="Methylated","State"=aggregate(all_M_count ~ State, data=chromatin_state_ranges, FUN=sum)$State,"Sites"=aggregate(all_M_count ~ State, data=chromatin_state_ranges, FUN=sum)$all_M_count),cbind("Methylation"="Unmethylated","State"=aggregate(all_U_count ~ State, data=chromatin_state_ranges, FUN=sum)$State,"Sites"=aggregate(all_U_count ~ State, data=chromatin_state_ranges, FUN=sum)$all_U_count)))
	chromatin_state_site_summary$Sites=as.numeric(as.character(chromatin_state_site_summary$Sites))
	 
	# Plot the distribution of CG sites with different methylation status among the chromatin states
	# This version does percentages
	pdf(paste0(project_id,"_",meth_context,"_methylation_status_by_chromatin_state.pdf"))
	print(ggplot(ddply(chromatin_state_site_summary, "Methylation", mutate, percent_Sites=Sites/sum(Sites) * 100), aes(x=Methylation, y=percent_Sites)) + geom_bar(aes(fill=State), stat="identity") + labs(x=paste0(meth_context," Sites (",project_id,")"), y="Percentage of Sites", fill="Chromatin State\n(Sequeira-Mendes et al, 2014)") + theme_minimal())
	dev.off()
	# This version does actual counts of sites
	pdf(paste0(project_id,"_",meth_context,"_methylation_status_by_chromatin_state.pdf"))
	print(ggplot(chromatin_state_site_summary, aes(x=Methylation, y=Sites)) + geom_bar(aes(fill=State), stat="identity") + labs(x=paste0(meth_context," Sites (Schmitz et al, 2011)"), y="Number of Sites", fill="Chromatin State\n(Sequeira-Mendes\n et al, 2014)") + facet_wrap(~Methylation, scales="free", ncol=6) + theme_minimal() + theme(strip.text.x = element_blank())) 
	dev.off()
	
	# In the Schmitz data there is some correlation among chromatin state segments between numbers of M->U sites and numbers of U->M sites:
	
	cat(paste0("Correlation between numbers of M->U and U->M sites per chromatin state segment (segments containing variable sites only): ",cor(chromatin_state_ranges[subjectHits(variant_chromatin_olaps),]$m_to_u_count,chromatin_state_ranges[subjectHits(variant_chromatin_olaps),]$u_to_m_count),"\n"))
	#[1] 0.1980373
	cat(paste0("Correlation between numbers of M->U and U->M sites per chromatin state segment (all segments, including sites with no variation)): ",cor(chromatin_state_ranges$m_to_u_count,chromatin_state_ranges$u_to_m_count)))
	#[1] 0.3644362
	
	
	
	
	# Zhang et al, 2015 https://academic.oup.com/nar/article/44/D1/D1148/2503132 present a database of DNase-hypersensitivity sites in the Arabidopsis genome. These sites should represent open chromatin, to allow to test the hypothesis that heritable changes in methylation status are confined or enriched in open chromatin
	
	DHS_ranges = read.table(file=reference_DHS_loci, sep="\t", header=FALSE)
	DHS_ranges=DHS_ranges[,c(1,4,5)]
	colnames(DHS_ranges) = c("Chrom","From","To")
	DHS_ranges$Chrom=paste0("CHR",substr(DHS_ranges$Chrom,4,4))

	# Plot distribution of DHS range sizes
	pdf(paste0(project_id,"_",meth_context,"_DHS_segment_lengths.pdf"))
	print(ggplot(DHS_ranges[DHS_ranges$To-DHS_ranges$From<10000,], aes(x=To-From)) + geom_histogram(aes(y= ..count..), binwidth=150) + theme_minimal())
	dev.off()
	
	# looks binomial with maximum at ~300nt
	### Fit a model to this histogram
	
	# Turn the DHS loci into GenomicRanges:
	DHS_granges=makeGRangesFromDataFrame(df = DHS_ranges, start.field = "From", end.field = "To",seqnames.field = "Chrom")

	# Find overlaps between DHS and CG sites with various patterns of variation
	DHS_ranges$ID=row.names(DHS_ranges)
	CG_DHS_olaps = findOverlaps(CG_site_ranges, DHS_granges)
	DHS_ranges$CG_site_count=table(DHS_ranges[subjectHits(CG_DHS_olaps),]$ID)[DHS_ranges$ID]
	DHS_ranges$CG_site_count[is.na(DHS_ranges$CG_site_count)]=0
	variant_DHS_olaps = findOverlaps(variant_call_ranges, DHS_granges)
	DHS_ranges$variant_count=table(DHS_ranges[subjectHits(variant_DHS_olaps),]$ID)[DHS_ranges$ID]
	DHS_ranges$variant_count[is.na(DHS_ranges$variant_count)]=0
	all_M_DHS_olaps = findOverlaps(all_M_ranges, DHS_granges)
	DHS_ranges$all_M_count=table(DHS_ranges[subjectHits(all_M_DHS_olaps),]$ID)[DHS_ranges$ID]
	DHS_ranges$all_M_count[is.na(DHS_ranges$all_M_count)]=0
	all_U_DHS_olaps = findOverlaps(all_U_ranges, DHS_granges)
	DHS_ranges$all_U_count=table(DHS_ranges[subjectHits(all_U_DHS_olaps),]$ID)[DHS_ranges$ID]
	DHS_ranges$all_U_count[is.na(DHS_ranges$all_U_count)]=0
	m_to_u_DHS_olaps = findOverlaps(m_to_u_call_ranges, DHS_granges)
	DHS_ranges$m_to_u_count=table(DHS_ranges[subjectHits(m_to_u_DHS_olaps),]$ID)[DHS_ranges$ID]
	DHS_ranges$m_to_u_count[is.na(DHS_ranges$m_to_u_count)]=0
	u_to_m_DHS_olaps = findOverlaps(u_to_m_call_ranges, DHS_granges)
	DHS_ranges$u_to_m_count=table(DHS_ranges[subjectHits(u_to_m_DHS_olaps),]$ID)[DHS_ranges$ID]
	DHS_ranges$u_to_m_count[is.na(DHS_ranges$u_to_m_count)]=0

	# Summarise the numbers of variant, all_M and all_U sites by DHS state:
	cat(paste0(sum(DHS_ranges$CG_site_count)," of ",length(CG_site_ranges)," ",meth_context," sites (",sum(DHS_ranges$CG_site_count)/length(CG_site_ranges),") lie within DHSs.\n"))
	cat(paste0(sum(DHS_ranges$variant_count)," of ",length(variant_call_ranges)," variable ",meth_context," sites (",sum(DHS_ranges$variant_count)/length(variant_call_ranges),") lie within DHSs.\n"))
	cat(paste0(sum(DHS_ranges$all_M_count)," of ",length(all_M_ranges)," all M ",meth_context," sites (",sum(DHS_ranges$all_M_count)/length(all_M_ranges),") lie within DHSs.\n"))
	cat(paste0(sum(DHS_ranges$all_U_count)," of ",length(all_U_ranges)," all U ",meth_context," sites (",sum(DHS_ranges$all_U_count)/length(all_U_ranges),") lie within DHSs.\n"))
	
	
	 
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	# Use changepoint detection to segment the genome based on levels of methylation at adjacent sites:
	# http://www.sciencedirect.com/science/article/pii/S0168165617315936
	# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1809-5
	
	install.packages("changepoint")
	library(changepoint)
	
	# Estimate change points in methylation proportion
	# Higher penalties result in larger segments
	# Try a range of penalty values and calculate the difference in the fragment length spectra between the changepoint analysis and the Sequeira-Mendes chromatin states
	
	for (cp_penalty in c(0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5)) {
		#cp_penalty=1.0
		cp_method="PELT"
		changepoints.amp = data.frame()
		changepoints.var = data.frame()
		for (this_chromosome in unique(all_reps_meth_status$Chromosome)) {
			changepoints.amp=rbind(changepoints.amp,cbind("Chromosome"=this_chromosome, "changepoint"=cpt.mean(cbind(all_reps_meth_status$Locus,across_samples_TM[,2]/across_samples_TM[,1])[all_reps_meth_status$Chromosome==this_chromosome & !is.na(across_samples_TM[,2]/across_samples_TM[,1]),2], method=cp_method, penalty="Manual", pen.value=cp_penalty, class=FALSE)))
			#cat(paste0(cp_method, " amplitude changepoint detection with penalty ", cp_penalty, " detected ", nrow(changepoints.amp[changepoints.amp$Chromosome==this_chromosome,]), " changepoints in chromosome ", this_chromosome, ".\n"))
			changepoints.var=rbind(changepoints.var,cbind("Chromosome"=this_chromosome, "changepoint"=cpt.var(cbind(all_reps_meth_status$Locus,across_samples_TM[,2]/across_samples_TM[,1])[all_reps_meth_status$Chromosome==this_chromosome & !is.na(across_samples_TM[,2]/across_samples_TM[,1]),2], method=cp_method, penalty="Manual", pen.value=cp_penalty, class=FALSE)))
			#cat(paste0(cp_method, " variance changepoint detection with penalty ", cp_penalty, " detected ", nrow(changepoints.var[changepoints.var$Chromosome==this_chromosome,]), " changepoints in chromosome ", this_chromosome, ".\n"))
		}
		changepoints.amp$changepoint=as.numeric(levels(changepoints.amp$changepoint)[changepoints.amp$changepoint])
		changepoints.var$changepoint=as.numeric(levels(changepoints.var$changepoint)[changepoints.var$changepoint])
		cat(paste0("Changepoint penalty ",cp_penalty," found ",nrow(changepoints.amp), " amplitude changepoints in the genome.\n"))
		cat(paste0("Changepoint penalty ",cp_penalty," found ",nrow(changepoints.var), " variance changepoints in the genome.\n"))

		# changepoints.amp and .var now contain pointers to the k'th entries of across_samples_TM, at which point changes occur, for each chromosome (which in turn correspond to rows of 	all_reps_meth_status)
		# Need to populate new locus columns with the actual loci of the changepoints, by resolving the pointers through a lookup to all_reps_meth_status
		changepoints.amp$Locus=NA
		changepoints.var$Locus=NA
		changepoint_segments.amp=data.frame()
		changepoint_segments.var=data.frame()
		for (this_chromosome in unique(all_reps_meth_status$Chromosome)) {
			changepoints.amp[changepoints.amp$Chromosome==this_chromosome,]$Locus=all_reps_meth_status[all_reps_meth_status$Chromosome==this_chromosome,][changepoints.amp[changepoints.amp$Chromosome==this_chromosome,]$changepoint,]$Locus
			changepoints.var[changepoints.var$Chromosome==this_chromosome,]$Locus=all_reps_meth_status[all_reps_meth_status$Chromosome==this_chromosome,][changepoints.var[changepoints.var$Chromosome==this_chromosome,]$changepoint,]$Locus
			cp_segment_start=c(1,changepoints.amp[changepoints.amp$Chromosome==this_chromosome,]$Locus+1+ifelse(meth_context=="CG",1,0))
			cp_segment_end=c(changepoints.amp[changepoints.amp$Chromosome==this_chromosome,]$Locus+ifelse(meth_context=="CG",1,0),sLengths[[paste0(substr(this_chromosome,1,1),tolower(substr(this_chromosome,2,3)),substr(this_chromosome,4,4))]])
		
			# Work out start and end points for the segment
			# Changepoints occur at the last site of each segment, so:
			# Segments start at 1, and at next site after end of previous segment. 
			# Segments end at changepoints (+1 for CG sites as they are 2nt long), and at end of chromosomes
			# This means segments are not contiguous - they have gaps between adjacent sites
			cp_segment_start=c(1,all_reps_meth_status[all_reps_meth_status$Chromosome==this_chromosome,][changepoints.amp[changepoints.amp$Chromosome==this_chromosome,]$changepoint+1,]$Locus)
			cp_segment_end=c(changepoints.amp[changepoints.amp$Chromosome==this_chromosome,]$Locus+ifelse(meth_context=="CG",1,0),sLengths[[paste0(substr(this_chromosome,1,1),tolower(substr(this_chromosome,2,3)),substr(this_chromosome,4,4))]])
		
			changepoint_segments.amp=rbind(changepoint_segments.amp,cbind(this_chromosome,cp_segment_start,cp_segment_end))
			rm(cp_segment_start, cp_segment_end)
		}	
		# Convert factors to numeric
		changepoint_segments.amp$cp_segment_start=as.numeric(levels(changepoint_segments.amp$cp_segment_start)[changepoint_segments.amp$cp_segment_start])
		changepoint_segments.amp$cp_segment_end=as.numeric(levels(changepoint_segments.amp$cp_segment_end)[changepoint_segments.amp$cp_segment_end])

	
		# Calculate the sum of absolute difference between the two densities, and find the penalty which minimises the difference in segment legnth spectra
		cat(paste0("For changepoint penalty of ",cp_penalty," difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is ", sum(abs(density(chromatin_state_ranges[,]$To-chromatin_state_ranges[,]$From)$y-density(changepoint_segments.amp[changepoint_segments.amp$this_chromosome %in% c("CHR1","CHR2","CHR3","CHR4","CHR5"),]$cp_segment_end-changepoint_segments.amp[changepoint_segments.amp$this_chromosome %in% c("CHR1","CHR2","CHR3","CHR4","CHR5"),]$cp_segment_start)$y)),"\n"))
	#	cat(paste0("For changepoint penalty of ",cp_penalty," difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is ", sum(abs(density(chromatin_state_ranges[chromatin_state_ranges$To-chromatin_state_ranges$From<10000,]$To-chromatin_state_ranges[chromatin_state_ranges$To-chromatin_state_ranges$From<10000,]$From)$y-density(changepoint_segments.amp[changepoint_segments.amp$cp_segment_end-changepoint_segments.amp$cp_segment_start<10000,]$cp_segment_end-changepoint_segments.amp[changepoint_segments.amp$cp_segment_end-changepoint_segments.amp$cp_segment_start<10000,]$cp_segment_start)$y)),"\n"))
	}
	
	
	# For Schmitz data, cp_penalty=1.0 minimises difference in segment length spectra between amplitude changepoint analysis and Sequeira-Mendes chromatin states
	#For changepoint penalty of 0.5 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00482043506838301
	#For changepoint penalty of 0.6 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00399490933083712
	#For changepoint penalty of 0.7 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.0040534050780842
	#For changepoint penalty of 0.8 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00370984112985995
	#For changepoint penalty of 0.9 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00358097534604698
	#For changepoint penalty of 1 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00357496896463758
	#For changepoint penalty of 1.1 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00358983816692342
	#For changepoint penalty of 1.2 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00369665748113977
	#For changepoint penalty of 1.3 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00395712689866591
	#For changepoint penalty of 1.4 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00396698779372347
	#For changepoint penalty of 1.5 difference in segment length spectra between changepoint analysis and Sequeira-Mendes chromatin state segments is 0.00399866524563563

	# Plot segment length spectra
	pdf(paste0(project_id,"_",meth_context,"_segment_length_spectra.pdf"))
	print(ggplot(data.frame(rbind.data.frame(cbind.data.frame("segment_length"=changepoint_segments.amp$cp_segment_end-changepoint_segments.amp$cp_segment_start,"analysis"="Changepoints"),cbind.data.frame("segment_length"=chromatin_state_ranges[,]$To-chromatin_state_ranges[,]$From,"analysis"="ChromatinStates"))), aes(x=segment_length, color=analysis)) + geom_density() + theme_minimal() +xlim(0,10000))
	dev.off()
	#ggplot(changepoint_segments.amp[changepoint_segments.amp$cp_segment_end-changepoint_segments.amp$cp_segment_start<10000,], aes(x=cp_segment_end-cp_segment_start)) + geom_histogram(aes(y= ..density..), binwidth=150) + theme_minimal()
	
	# The peak in density of segment lengths is slightly lower from changepoint analysis (~400nt) than for the chromatin states analysis (~700nt).  However, there are considerably more considerably longer segments from the changepoints analysis.
	
	# Work out the mean methylation level for each changepoint-defined segment
	
	# First the row-wise version - takes ages (half hour?)
	#changepoint_segments.amp$mean_methylation=NA
	#for (i in 1:nrow(changepoint_segments.amp)) {
	#	changepoint_segments.amp[i,"mean_methylation"]=mean(across_samples_meth_status[across_samples_meth_status$Chromosome==changepoint_segments.amp[i,]$this_chromosome & across_samples_meth_status$Locus>=changepoint_segments.amp[i,]$cp_segment_start & across_samples_meth_status$Locus<=changepoint_segments.amp[i,]$cp_segment_end,]$average_methylation)
	#}
	
	# Now the sqldf version.  Cleaner code but relies on external package and still takes ages
	#library(sqldf)
	changepoint_segments.amp = sqldf("SELECT cp.this_chromosome, cp.cp_segment_start, cp.cp_segment_end, AVG(ms.average_methylation) mean_methylation, POWER(STDEV(ms.average_methylation),2) variance_methylation  FROM [changepoint_segments.amp] cp JOIN across_samples_meth_status ms ON cp.this_chromosome = ms.Chromosome AND ms.Locus BETWEEN cp.cp_segment_start AND cp.cp_segment_end GROUP BY cp.this_chromosome, cp.cp_segment_start")
	
	# Add convenience columns to show integer counts from 0-10 for mean and variance of segment methylation level
	changepoint_segments.amp$methylation_count=round(10*changepoint_segments.amp$mean_methylation)
	changepoint_segments.amp$meth_var_count=round(20*changepoint_segments.amp$variance_methylation)
	 
	# Save the changepoint segment methylation averages so we don't need to recalculate this ever
	saveRDS(changepoint_segments.amp, file=paste0(project_id,"_",meth_context,"_changepoint_segments.amp.rds"))

	# Output changepoint segments with mean methylation values to a convenience file for visualisation
	write.table(changepoint_segments.amp,file=paste0(project_id,"_",meth_context,"_changepoint_segments_",cp_penalty,".tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
	
	# Plot segment mean and variance methylation in genomic context
	pdf(paste0(project_id,"_",meth_context,"_segments_in_genome.pdf"))
	print(ggplot(changepoint_segments.amp, aes(x=cp_segment_start, y=mean_methylation)) + geom_point(aes(colour=variance_methylation), alpha=0.4, size=0.1) + scale_colour_gradient(low="green", high="red") + facet_grid(this_chromosome ~ .) + theme_minimal())
	dev.off()
	
	# Plot methylation variation landscape by segment
	pdf(paste0(project_id,"_",meth_context,"_segments_methylation_variation.pdf"))
	print(ggplot(changepoint_segments.amp, aes(x=mean_methylation, y=variance_methylation/mean_methylation)) + geom_point(aes(colour=variance_methylation/mean_methylation), alpha=0.4, size=0.1) + scale_colour_gradient(low="green", high="red") + theme_minimal() + geom_density_2d())
	dev.off()
	
	# Plot methylation variation by segment length
	pdf(paste0(project_id,"_",meth_context,"_segments_methylation_variation_by_length.pdf"))
	print(ggplot(changepoint_segments.amp, aes(x=log(cp_segment_end+1-cp_segment_start), y=variance_methylation/mean_methylation)) + geom_point(aes(colour=variance_methylation/mean_methylation), alpha=0.4, size=0.1) + scale_colour_gradient(low="green", high="red") + theme_minimal() + geom_density_2d())
	dev.off()
	
	# Plot methylation level and variation by segment length
	pdf(paste0(project_id,"_",meth_context,"_segments_methylation_and_variation_by_length.pdf"))
	print(ggplot(changepoint_segments.amp, aes(x=log(cp_segment_end+1-cp_segment_start), y=mean_methylation)) + geom_point(aes(colour=variance_methylation/mean_methylation), alpha=0.4, size=0.1) + scale_colour_gradient(low="green", high="red") + theme_minimal() + geom_density_2d())
	dev.off()

	# Corresponding plot for genes
	pdf(paste0(project_id,"_",meth_context,"_genes_methylation_and_variation_by_length.pdf"))
	print(ggplot(gff.genes, aes(x=log(V5+1-V4), y=average_methylation)) + geom_point(aes(colour=variance_methylation/average_methylation), alpha=0.4, size=0.1)+ scale_colour_gradient(low="green", high="red") + theme_minimal() + geom_density_2d())
	dev.off()
	
	# Corresponding plot for transposons
	pdf(paste0(project_id,"_",meth_context,"_transposons_methylation_and_variation_by_length.pdf"))
	print(ggplot(gff.transposons, aes(x=log(V5+1-V4), y=average_methylation)) + geom_point(aes(colour=variance_methylation/average_methylation), alpha=0.4, size=0.1)+ scale_colour_gradient(low="green", high="red") + theme_minimal() + geom_density_2d())
	dev.off()
	
	# Turn changepoint segments into GenomicRanges
	changepoint_granges=makeGRangesFromDataFrame(df = changepoint_segments.amp, start.field = "cp_segment_start", end.field = "cp_segment_end",seqnames.field = "this_chromosome")

	### Find overlaps between changepoint segments and variable sites
	
	### Find overlaps between changepoint segments and chromatin states
	
	

	# Non-parametric changepont detection:
	
	install.packages("changepoint.np")
	library(changepoint.np)

	cp_penalty=1.5
	cp_method="PELT"
	changepoints.amp = data.frame()
	#changepoints.var = data.frame()
	for (this_chromosome in unique(all_reps_meth_status$Chromosome)) {
		changepoints.amp=rbind(changepoints.amp,cbind("Chromosome"=this_chromosome, "changepoint"=cpt.np(cbind(all_reps_meth_status$Locus,across_samples_TM[,2]/across_samples_TM[,1])[all_reps_meth_status$Chromosome==this_chromosome & !is.na(across_samples_TM[,2]/across_samples_TM[,1]),2], method=cp_method, penalty="Manual", pen.value=cp_penalty, class=FALSE)))
		cat(paste0(cp_method, " amplitude changepoint detection with penalty ", cp_penalty, " detected ", nrow(changepoints.amp[changepoints.amp$Chromosome==this_chromosome,]), " changepoints in chromosome ", this_chromosome, ".\n"))
		#changepoints.var=rbind(changepoints.var,cbind("Chromosome"=this_chromosome, "changepoint"=cpt.var(cbind(all_reps_meth_status$Locus,across_samples_TM[,2]/across_samples_TM[,1])[all_reps_meth_status$Chromosome==this_chromosome & !is.na(across_samples_TM[,2]/across_samples_TM[,1]),2], method=cp_method, penalty="Manual", pen.value=cp_penalty, class=FALSE)))
		#cat(paste0(cp_method, " variance changepoint detection with penalty ", cp_penalty, " detected ", nrow(changepoints.var[changepoints.var$Chromosome==this_chromosome,]), " changepoints in chromosome ", this_chromosome, ".\n"))
	}
	changepoints.amp$changepoint=as.numeric(levels(changepoints.amp$changepoint)[changepoints.amp$changepoint])
	#changepoints.var$changepoint=as.numeric(levels(changepoints.var$changepoint)[changepoints.var$changepoint])
	cat(paste0("Found ",nrow(changepoints.amp), " amplitude changepoints in the genome.\n"))
	#cat(paste0("Found ",nrow(changepoints.var), " variance changepoints in the genome.\n"))




	
	
	### Check what is going on at the centromeres
	
	# Centromeres (according to https://www.biostars.org/p/18782/)
	#Chr1 TAIR BAC_cloned_genomic_insert 15086046 15087045 . + . ID=AssemblyUnit_CEN1;Name=CEN1
	#Chr2 TAIR BAC_cloned_genomic_insert 3607930 3608929 . + . ID=AssemblyUnit_CEN2;Name=CEN2
	#Chr3 TAIR BAC_cloned_genomic_insert 13799418 13800417 . + . "ID=AssemblyUnit_""CEN3, PT.2OF3"";Name=""CEN3, PT.2OF3"""
	#Chr3 TAIR BAC_cloned_genomic_insert 14208953 14209952 . + . "ID=AssemblyUnit_""CEN3, PT.3OF3"";Name=""CEN3, PT.3OF3"""
	#Chr4 TAIR BAC_cloned_genomic_insert 3956022 3957021 . + . ID=AssemblyUnit_CEN4;Name=CEN4
	#Chr5 TAIR BAC_cloned_genomic_insert 11725025 11726024 . + . ID=AssemblyUnit_CEN5;Name=CEN5

	# Zhang et al, 2006 map methylation in Arabidopsis genome
	
	# Kawabe et al, 2006 identified centromeric regions: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1526690/
	#Chr1: 14.16:15.93
	#Chr2: 3.05:4.32
	#Chr3: 13.49:14.24
	#Chr4: 3.67:4.07
	#Chr5: 11.19:12.57

	
	
	
	
	
	
	# Find how many coding bases the primary transcript has


	# Generate a set of variable calls for the "all M" and "all U" sites randomly selected previously
	loci_sample <- sample(1:nrow(all_samples_meth_status[all_m==TRUE,]), nrow(variant_calls),
  	replace=FALSE)
	loci_sample=all_samples_meth_status[all_m==TRUE,][loci_sample,]
	loci_sample=loci_sample[order(loci_sample$Chromosome,loci_sample$Locus),]
	all_M_ranges=makeGRangesFromDataFrame(df = loci_sample, start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	all_M_olaps = findOverlaps(all_M_ranges, gene_ranges)
	gff.genes$all_M_variant_count=table(gff.genes[subjectHits(all_M_olaps),]$gene_ID)[gff.genes$gene_ID]
	gff.genes$all_M_variant_count=ifelse(is.na(gff.genes$all_M_variant_count),0,gff.genes$all_M_variant_count)
	
	loci_sample <- sample(1:nrow(all_samples_meth_status[all_u==TRUE,]), nrow(variant_calls),
  	replace=FALSE)
	loci_sample=all_samples_meth_status[all_u==TRUE,][loci_sample,]
	loci_sample=loci_sample[order(loci_sample$Chromosome,loci_sample$Locus),]
	all_U_ranges=makeGRangesFromDataFrame(df = loci_sample, start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	all_U_olaps = findOverlaps(all_U_ranges, gene_ranges)
	gff.genes$all_U_variant_count=table(gff.genes[subjectHits(all_U_olaps),]$gene_ID)[gff.genes$gene_ID]
	gff.genes$all_U_variant_count=ifelse(is.na(gff.genes$all_U_variant_count),0,gff.genes$all_U_variant_count)
	
	# How many invariable sites lie within annotated genes?
	cat(paste(sum(gff.genes$all_M_variant_count),"of",nrow(variant_calls),"randomly sampled 'all M' sites lie within annotated genes.\n"))
	cat(paste(sum(gff.genes$all_U_variant_count),"of",nrow(variant_calls),"randomly sampled 'all U' sites lie within annotated genes.\n"))

	# Compare the distributions of site counts per gene between variant sites, random all M and random all U sites.
	#### This needs redoing: Resample randomly chosen 'all M' and 'all U' sites to find the same number that lie within genes as the variable sites
	sites_per_gene=data.frame(rbind(cbind(Count=gff.genes$variant_count,type="Variable sites"),cbind(Count=gff.genes$all_M_variant_count,type="Methylated sites"),cbind(Count=gff.genes$all_U_variant_count,type="Unmethylated sites")))
	rownames(sites_per_gene) = NULL
	sites_per_gene$Count=as.integer(levels(sites_per_gene$Count)[sites_per_gene$Count])

	pdf(paste0(project_id,"_",meth_context,"_variable_sites_per gene_frequency.pdf"))
	print(ggplot(sites_per_gene, aes(x=Count, colour=type)) +geom_freqpoly(binwidth=1))
	dev.off()

	#ggplot(rbind(gaps, sample_gaps_m, sample_gaps_u), aes(x=gap, colour=class)) + geom_density() + xlim(0, 20000) + coord_trans(x = "log10")
	
	# Test whether the distributions of gap size vary significantly between groups
	ks.test(sites_per_gene[sites_per_gene$type=="Variable sites",]$Count,sites_per_gene[sites_per_gene$type=="Methylated sites",]$Count)
	ks.test(sites_per_gene[sites_per_gene$type=="Variable sites",]$Count,sites_per_gene[sites_per_gene$type=="Unmethylated sites",]$Count)
	ks.test(sites_per_gene[sites_per_gene$type=="Methylated sites",]$Count,sites_per_gene[sites_per_gene$type=="Unmethylated sites",]$Count)
	

	# Add column for sample names to sample_info
	sample_info$name=paste0(sample_info$Identifier,"_CG_Gen",sample_info$Generation,"_Line",sample_info$Line,"_Rep",sample_info$Rep)

	### Load each sample's single-c files to seqmonk
	### Rename the samples in seqmonk according to sample_info$name
	### Load the variant_calls file to seqmonk as an annotation track
	
	
	# From the sites which vary between lines/generations with consistency among reps, which align with genes, identify which are changing from M->U and vice versa
	
	
	

	# Automate annotation and selection of sites
	#library(sqldf)
	
	
	# Plot things in genomic context
	biocLite("ggbio")
	library(ggbio)
	
	autoplot(variant_call_ranges, layout="linear")
	autoplot(variant_call_ranges, layout="circle")
	autoplot(variant_call_ranges, layout="karyogram")
	# Assign some metadata to ranges:
	mcols(variant_call_ranges)$Up <- 0.5
	plotGrandLinear(variant_call_ranges, aes(y=score, col=up))
	
	# To do:

	# Remove SNPs
	# How to 'normalise' for undetected data/sites?
	
	# How does stability vary among sites in different genomic and methylation contexts?
	# What is going on with sites with Mismatched status?  Is it real?  Is it inherited?
	# Is a transposon wholly methylated?
	# Does a whole transposon lose its methylation?
	# On-rate and off-rate of methylation changes
	
	# Annotate each CG locus with lots of binary and numeric variables:
	
	# Methylation status of the site:
	# M in all samples (precise)
	# U in all samples (precise)
	# M/U in all samples (fuzzy - allow I and P calls in some samples to be ignored in the assessment)
	# Variable among samples precise (all samples have M or U)
	# Variable fuzzy (some samples have I or P calls)

	# Methylation context:
	# In methylated gene (how is this defined? 60% cutoff?)
	# In unmethylated gene (ditto)
	# In methylated transposon
	# In unmethylated transposon
	# In methylated intergenic locus
	# In unmethylated intergenic locus
	# Neighbourhood methylation status (1 site either side, 1 site 5', 1 site 3')
	# Neighbourhood methylation status (2 sites either side, 2 sites 5', 2 sites 3'))
	# Neighbourhood methylation status (3 sites either side, 2 sites 5', 3 sites 3'))

	# Incorporate Sequeira-Mendes 9 chromatin states model
	
	# Genomic status:
	# In CpG island (where can I find these annotated?)
	# In gene
	# In mRNA
	# In CDS
	# In intron
	# In 5' UTR
	# In 3' UTR
	# In exon
	# In transposon
	# In pseudogene
	# In lncRNA
	
	
	# Generate a diversity index score for each locus, check the mean diversity index of each gene
#	sum(all_samples_meth_status[,4:19]=="U")
#	for(locus_no in 1:nrow(loci_sample)) {
#	for(locus_no in 1:10) {
#		if (prev_chrom==levels(loci_sample[locus_no,]$Chromosome)[loci_sample[locus_no,]$Chromosome]) {
#			sample_gaps[locus_no-no_chroms,1]=(loci_sample[locus_no,]$Locus)-(loci_sample[locus_no-1,]$Locus)
#		} else {
#			no_chroms = no_chroms + 1
#		}
#		prev_chrom=levels(loci_sample[locus_no,]$Chromosome)[loci_sample[locus_no,]$Chromosome]
#	}
	
	
	
	#source('http://bioconductor.org/biocLite.R')
	#biocLite('phyloseq')

	#library(CluMix)
	#mix.heatmap(mixdata, rowmar=7, legend.mat=TRUE)

#
# To do:

# Flexible analysis of M/U calls across groups of samples
# Assign project metadata to samples, and compare methylation status among groups
# Read GFF files for annotations and annotate each meth_context locus
# Carry out annotation-dependent analyses of methylation status
# Visualisations in genomic context
# Shiny browser?



# Installing Gviz or methylKit failed on Windows and Ubuntu
# On Ubuntu:
#sudo apt-get install aptitude
#sudo apt-get install libcurl4-openssl-dev
#sudo apt-get install libxml2-dev


# Plot stuff in genomic context
#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
#biocLite("GenomeGraphs")

#biocLite("ensembldb")
#library(AnnotationHub)
#ah <- AnnotationHub()

#biocLite("TxDb.Athaliana.BioMart.plantsmart28")

#library(ensembldb)
 
#	gff <-rtracklayer::import(reference_gff)
#useDevel(TRUE)
#biocLite("Gviz")
#useDevel(FALSE)


 
} # end action=="analyse"
