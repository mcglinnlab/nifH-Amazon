library(phyloseq)

compute_unifrac = function(abu_file, tree_file) {
  tree = read.tree(file.path('./phylogenetic_trees',tree_file))
  root_tip = root_tip = 'FIB10_D58U3'
  tree = root(tree, outgroup=root_tip, resolve.root=TRUE)
  abu = read.csv(file.path('./abundance_data',abu_file))
  abu = as.data.frame(t(abu))                     ##transpose pgylogenetic data
  names(abu) = as.matrix(abu[1,  ])               ## change the data as matrix
  rownames_abu = as.matrix(rownames(abu))         ##change the row names as matrix
  abu = abu[-1, ]                                 ## ----except first row
  abu = apply(abu, 2, as.numeric)                 ## set all columns to numeric class
  rownames(abu) = rownames_abu[-1,]               ## Adjust the row names
  abu_data = otu_table(abu, taxa_are_rows=F)
  tree_abu = phyloseq(abu_data, tree)
  result = UniFrac(tree_abu, weighted=TRUE)
  return(result)
}

abu_files = dir('./abundance_data/')

for(i in seq_along(abu_files)) {
  split_name = strsplit(abu_files[i], '_')       ## split up the file name into its components
  percent = sub('%', '', split_name[[1]][1])     ## identify the % similarity
  marker = sub('.csv', '', split_name[[1]][4])   ## identify the marker
  tree_file = paste(percent, '%_FSP_S1_', marker, '_tree.tree', sep='')
  start = proc.time()
  result = compute_unifrac(abu_files[i], tree_file)
  end = proc.time()
  print(end - start)
  export_file = paste('./unifrac_dist/unifrac_dist', percent, marker, 'abu.txt',
                      sep='_')
  write.table(as.matrix(result), file=export_file, row.names=T, col.names=T)
}


