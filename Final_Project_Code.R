library(tidyverse)


#Set File names and types so the cleaning function can create a new cleaned object for
#each data set

names = c('GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_median_tpm',
          'gene_reads_brain_amygdala')
types = c('.gct',
          '.gct')

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