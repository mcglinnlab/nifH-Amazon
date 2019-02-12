## create phylo tree at given cutoff of phylo distance and/or base pair distance

library('ape')
library('seqinr')

## read in genetic information
fasta = read.alignment('./nifh_DNA_prot_data/NifH_phylogentic_trees/508_Sub1_FPS_ref_comb_maligned.fasta',
                           format = "fasta")
tree = read.tree('./nifh_DNA_prot_data/NifH_phylogentic_trees/508_Sub1_FPS_ref_comb_tree.tree')


## read in treatment and abundance data
trt = read.csv('./nifh_DNA_prot_data/Dist_Decay_analysis/treatment_forest_pasture,secForest_dd_coords.csv')
abu = read.csv('./nifh_DNA_prot_data/NifH_phylogentic_trees/90%_FSP_S1_prot_named_abu_data.csv')

## NOTE:
## tree files have 507 tips but the abundance file has only 430 datapoints
## need to include those names in the tree that are not in the abundance file as zeros
ref_tips = tree$tip.label[!(tree$tip.label %in% abu[ , 'Seq_Names'])]
#ref_dat = data.frame(Seq_Names = ref_tips)

## create new abundance file
#abu = merge(abu, ref_dat, all=T)

all(names(abu)[-(1:4)] == trt$sample)

habitat = unique(trt$system)
hab_sum = matrix(NA, ncol= length(habitat), nrow= nrow(abu))
colnames(hab_sum) = habitat
rownames(hab_sum) = abu$Seq_Names
for(i in seq_along(habitat)) {
  hab_sum[ , i] = apply(abu[ , which(names(abu) %in% trt$sample[trt$system == habitat[i]])], 1, sum, na.rm=T)
}

hab_prop = hab_sum / apply(hab_sum, 1, sum)

dist_align = as.matrix(dist.alignment(fasta, matrix = "similarity"))

dist_tree = cophenetic.phylo(tree)

## if less distant then the cutoff tips are lumped to common node
cutoff = .1 ## 10% distant 

merge_align = ifelse(dist_align < cutoff, T, F)
merge_tree = ifelse(dist_tree < cutoff, T, F)

diag(merge_align) = NA
diag(merge_tree) = NA

## find any sequence with a TRUE this is a tip to be dropped

tips_to_drop_align = names(which(apply(merge_align, 2, sum, na.rm=T) > 0))
tips_to_drop_tree = names(which(apply(merge_tree, 2, sum, na.rm=T) > 0))


ptree_align = drop.tip(tree, tip=tips_to_drop_align) 
ptree_tree = drop.tip(tree, tip=tips_to_drop_tree) 

par(mfrow=c(1,3))
plot(tree, type='fan', show.tip.label=F)
plot(ptree_align, type='fan', show.tip.label=F)
plot(ptree_tree,type='fan', show.tip.label=F)

par(mfrow=c(1,1))
plot(ptree_align, type='unrooted', show.tip.label=F)

par(mfrow=c(1,1))
plot(ptree_align, type='fan', show.tip.label=F, use.edge.length=T)

## phylogram with pie charts
## only use pies for tips still on reduced tree and with given amt of abundance
abu_cutoff = 40
pie_tips = rownames(hab_sum)
pie_tips = pie_tips[which(apply(hab_sum, 1, sum) >= abu_cutoff)]
pie_tips = pie_tips[pie_tips %in% ptree_align$tip.label]

hab_sum_subset = hab_sum[match(pie_tips, row.names(hab_sum)), ]

pie_size = log10(apply(hab_sum_subset, 1, sum))
pie_size = pie_size / max(pie_size) * .5

plot(ptree_align, type='fan', show.tip.label=F, use.edge.length=T)
tiplabels(pie=hab_prop, tip= match(pie_tips, ptree_align$tip.label),
          cex=pie_size, piecol=c('blue', 'purple', 'red'))

plot(ptree_align, type='p', show.tip.label=F, use.edge.length=T)
tiplabels(pie=hab_prop, cex=pie_size, piecol=c('blue', 'purple', 'red'))

## attempt to make 3 tree figure: one for each habitat
for(i in 1:3) {
 tips = c(ref_tips, row.names(hab_sum)[hab_sum[ , i] > 20])
 subtree = drop.tip(tree, which(!(tree$tip.label %in% tips)))
 hab_sum_sub = hab_sum[row.names(hab_sum) %in% subtree$tip.label, ]
 hab_prop_sub =  hab_prop[row.names(hab_prop) %in% subtree$tip.label, ]
 pie_size = log10(apply(hab_sum_sub, 1, sum) + 1)
 pie_size = pie_size / max(pie_size) * 1
 pie_tips = match(names(pie_size), subtree$tip.label)
 if(i == 1) {
   par(mfrow=c(1,3))
 }
 plot(subtree, type='f', show.tip.label=F, use.edge.length=T)
 tiplabels(pie=hab_prop_sub, tip=pie_tips,
           cex=pie_size, piecol=c('blue', 'purple', 'red'))
}

## attempt to make 3 tree figure: one for each habitat
## in this case each tree is identical only colors change
cols = c('blue', 'purple', 'red')
for(i in 1:3) {
  tips = c(ref_tips, row.names(hab_sum)[hab_sum[ , i] > 20])
  hab_sum_sub = hab_sum[row.names(hab_sum) %in% tree$tip.label, ]
  hab_prop_sub =  hab_prop[row.names(hab_prop) %in% tree$tip.label, ]
  pie_size = log10(apply(hab_sum_sub, 1, sum) + 1)
  #pie_size = pie_size / max(pie_size) * 1
  pie_tips = match(names(pie_size), tree$tip.label)
  if(i == 1) {
    par(mfrow=c(1,3))
  }
  plot(tree, type='f', show.tip.label=F, use.edge.length=T)
  tiplabels(pie=rep(1, length(pie_tips)), tip=pie_tips,
            cex=pie_size, piecol=cols[i])
}




## convert to ultrametric tree
chronotree = chronos(tree) 
pchronotree_align = drop.tip(chronotree, tip=tips_to_drop_align) 

plot(pchronotree_align, show.tip.label=F,type='f')
tiplabels(pie=hab_prop, cex=.25, piecol=c('blue', 'purple', 'red'))
