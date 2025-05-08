library(tidyverse)
install.packages('pheatmap')
library(pheatmap)


#Set File names and types so the cleaning function can create a new cleaned object for
#each data set

names = c('gene_reads_brain_amygdala',
          'GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS')
types = c('.gct',
          '.txt')
data_names = c()

#cleaning function, takes file name and type as input, dynamically names object
#based on file name
cleaning = function(file_name, file_type){
  
  #read file and clean any comment lines, then store in global variable
  file = readLines(paste(file_name, file_type, sep =''))
  file[2] = paste('#',file[2])
  
  #clean the file of any comments
  clean_file = file[!grepl('^[#]', file)]
  
  #temporarily store the cleaned data for assignment
  temp_file = tempfile()
  writeLines(clean_file, temp_file)
  
  #assign data to global matrix
  name_match = sub(".*v[0-9]+(?:\\.[0-9]+)*_", "", file_name)
  assign(name_match, read_delim(temp_file, delim = '\t'), envir = .GlobalEnv)
  
  
  
}



#loop to iterate through all provided files for cleaning
for (i in 1:length(names))
  cleaning(names[i], types[i])

gene_reads_brain_amygdala = as.data.frame(t(gene_reads_brain_amygdala))
#Annotations_SubjectPhenotypesDS = as.data.frame(t(Annotations_SubjectPhenotypesDS))

GTEXids = row.names(gene_reads_brain_amygdala)[c(-1,-2)]
gene_reads_brain_amygdala = add_column(gene_reads_brain_amygdala, 
                                       SUBJID = c('na', 'na', str_extract(GTEXids, '^GTEX-[a-zA-Z0-9]+')), 
                                       .before = 1)

joined = left_join(gene_reads_brain_amygdala, Annotations_SubjectPhenotypesDS, by = 'SUBJID')
joined = relocate(joined, SEX, AGE, DTHHRDY,.before = 2)
joined = joined[-2,]

new_colnames = as.character(joined[1, 5:ncol(joined)])
colnames(joined)[5:ncol(joined)] = new_colnames
joined = joined[-1, ]

joined = na.omit(joined, target.colnames = AGE)


age_groups = joined$AGE
gene_counts = joined[, 5:ncol(joined)]
gene_by_age = cbind(AGE = age_groups, gene_counts)
gene_by_age$AGE = factor(gene_by_age$AGE)
avg_expression = aggregate(. ~ AGE, data = gene_by_age, FUN = function(x) mean(as.numeric(x), na.rm = TRUE))


rownames(avg_expression) = avg_expression$AGE
avg_expression = avg_expression[, -1]

#change the indexing here to change which genes you are looking at
#current set up is for genes 1 through 100
short = avg_expression[, 1:100]
scaled_expr = t(scale(t(as.matrix(short))))
pheatmap(scaled_expr,
         main = "Average Gene Expression by Age Range",
         cluster_cols = FALSE,
         cluster_rows = FALSE)


