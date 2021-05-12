d:
cd Projects\Ancestral_Methylation_Analysis
# check if alignment can be read by RAxML
..\..\Software\RAxML\raxmlHPC.exe -f c -m BINCAT -s 4-summary\tree_samples.fas -n draft_tree 
..\..\Software\RAxML\raxmlHPC.exe -f c -m BINCAT -s 4-summary\tree_samples.fas.reduced -n draft_tree -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# check that all settings are correct
..\..\Software\RAxML\raxmlHPC.exe -f a -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples.fas.reduced -n draft_tree -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis -y --flag-check

# get rapid parsimony tree plus rapid Bootstrap analysis and search for best-scoring ML tree
..\..\Software\RAxML\raxmlHPC.exe -f a -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples.fas.reduced -n draft_tree -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis -y -p 12345

# final thorough optimisation of ML tree from rapid bootstrap search
#..\..\Software\RAxML\raxmlHPC.exe -f T -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples.fas.reduced -t 5-analysis\draft_tree.nwk -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# simple rooting based on most balanced tree
..\..\Software\RAxML\raxmlHPC.exe -f I -t D:\Projects\Ancestral_Methylation_Analysis\5-analysis\RAxML_bipartitions.draft_tree -m BINCAT -n draft_tree_rooted -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# compute marginal ancestral states on a ROOTED reference tree provided with -t
..\..\Software\RAxML\raxmlHPC.exe -f A -t D:\Projects\Ancestral_Methylation_Analysis\5-analysis\RAxML_rootedTree.draft_tree_rooted -m BINCAT -n draft_tree_rooted_anc -s 4-summary\tree_samples.fas.reduced -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis



## Schmitz data only, using gamma model
# check if alignment can be read by RAxML
..\..\Software\RAxML\raxmlHPC.exe -f c -m BINGAMMA -s 4-summary\tree_samples_schmitz.fas -n draft_tree_schmitz 
..\..\Software\RAxML\raxmlHPC.exe -f c -m BINGAMMA -s 4-summary\tree_samples_schmitz.fas.reduced -n draft_tree_schmitz -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# check that all settings are correct
..\..\Software\RAxML\raxmlHPC.exe -f a -m BINGAMMA -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples_schmitz.fas.reduced -n draft_tree_schmitz -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis -y --flag-check

# get rapid parsimony tree plus rapid Bootstrap analysis and search for best-scoring ML tree
..\..\Software\RAxML\raxmlHPC.exe -f a -m BINGAMMA -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples_schmitz.fas.reduced -n draft_tree_schmitz -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis -y -p 12345

# final thorough optimisation of ML tree from rapid bootstrap search
#..\..\Software\RAxML\raxmlHPC.exe -f T -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples.fas.reduced -t 5-analysis\draft_tree.nwk -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# simple rooting based on most balanced tree
..\..\Software\RAxML\raxmlHPC.exe -f I -t D:\Projects\Ancestral_Methylation_Analysis\5-analysis\RAxML_bipartitions.draft_tree_schmitz -m BINGAMMA -n draft_tree_schmitz_rooted -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# compute marginal ancestral states on a ROOTED reference tree provided with -t
..\..\Software\RAxML\raxmlHPC.exe -f A -K MK -t D:\Projects\Ancestral_Methylation_Analysis\5-analysis\RAxML_rootedTree.draft_tree_schmitz_rooted -m BINGAMMA -n draft_tree_schmitz_rooted_anc -s 4-summary\tree_samples_schmitz.fas.reduced -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis


## Schmitz data, BINCAT model, multifurcating constraint tree (this doesn't really work probably - RAxML like bifurcating trees)
# check if alignment can be read by RAxML
..\..\Software\RAxML\raxmlHPC.exe -f c -m BINCAT -s 4-summary\tree_samples_schmitz.fas -n draft_tree_schmitz 
..\..\Software\RAxML\raxmlHPC.exe -f c -m BINCAT -s 4-summary\tree_samples_schmitz.fas.reduced -n draft_tree_schmitz -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# check that all settings are correct
..\..\Software\RAxML\raxmlHPC.exe -f a -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples_schmitz.fas.reduced -g 4-summary\multifurcating_constraint_tree_schmitz.nwk -n draft_tree_schmitz -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis --flag-check

# get rapid parsimony tree plus rapid Bootstrap analysis and search for best-scoring ML tree
..\..\Software\RAxML\raxmlHPC.exe -f a -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples_schmitz.fas.reduced -g 4-summary\multifurcating_constraint_tree_schmitz.nwk -n draft_tree_schmitz -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis -p 12345
..\..\Software\RAxML\raxmlHPC.exe -f d -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples_schmitz.fas.reduced -g 4-summary\multifurcating_constraint_tree_schmitz.nwk -n draft_tree_schmitz -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis -p 12345


# final thorough optimisation of ML tree from rapid bootstrap search
#..\..\Software\RAxML\raxmlHPC.exe -f T -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples.fas.reduced -t 5-analysis\draft_tree.nwk -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# simple rooting based on most balanced tree
..\..\Software\RAxML\raxmlHPC.exe -f I -t D:\Projects\Ancestral_Methylation_Analysis\5-analysis\RAxML_bipartitions.draft_tree_schmitz -m BINCAT -g 4-summary\multifurcating_constraint_tree_schmitz.nwk -n draft_tree_schmitz_rooted -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# compute marginal ancestral states on a ROOTED reference tree provided with -t
..\..\Software\RAxML\raxmlHPC.exe -f A -t D:\Projects\Ancestral_Methylation_Analysis\5-analysis\RAxML_rootedTree.draft_tree_schmitz_rooted -m BINCAT -n draft_tree_schmitz_rooted_anc -s 4-summary\tree_samples_schmitz.fas.reduced -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis


## All data, BINCAT model
# check if alignment can be read by RAxML
..\..\Software\RAxML\raxmlHPC.exe -f c -m BINCAT -s 4-summary\tree_samples.fas -n draft_tree 
..\..\Software\RAxML\raxmlHPC.exe -f c -m BINCAT -s 4-summary\tree_samples.fas.reduced -n draft_tree -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# check that all settings are correct
..\..\Software\RAxML\raxmlHPC.exe -f a -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples.fas.reduced -n draft_tree -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis -y --flag-check

# get rapid parsimony tree plus rapid Bootstrap analysis and search for best-scoring ML tree
..\..\Software\RAxML\raxmlHPC.exe -f a -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples.fas.reduced -n draft_tree -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis -y -p 12345

# final thorough optimisation of ML tree from rapid bootstrap search
#..\..\Software\RAxML\raxmlHPC.exe -f T -m BINCAT -N 100 -x 12345 -B 0.03 -K MK -O -s 4-summary\tree_samples.fas.reduced -t 5-analysis\draft_tree.nwk -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# simple rooting based on most balanced tree
..\..\Software\RAxML\raxmlHPC.exe -f I -t D:\Projects\Ancestral_Methylation_Analysis\5-analysis\RAxML_bipartitions.draft_tree -m BINCAT -n draft_tree_rooted -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis

# compute marginal ancestral states on a ROOTED reference tree provided with -t
..\..\Software\RAxML\raxmlHPC.exe -f A -K MK -t D:\Projects\Ancestral_Methylation_Analysis\5-analysis\RAxML_rootedTree.draft_tree_rooted -m BINCAT -n draft_tree_rooted_anc -s 4-summary\tree_samples.fas.reduced -w D:\Projects\Ancestral_Methylation_Analysis\5-analysis
