library(shiny)
library(stringr) #needed to split columns in dataframe
library(dplyr) 
library(ggplot2)
library(ggpubr)

################################################################################
############################# FUNCTIONS ########################################
################################################################################

# Creates sequence based on start and stop values
create_sequence <- function (start, end){
  return (seq(start, end))
}

# Returns the number of the strand, 
# e.g. virus_2020 returns 2020
get_ordered <- function (name){
  name <- strsplit(name, split='_')
  return(name[[1]][2])
}
# Trims leading white space 
trim_leading <- function (x)  sub("^\\s+", "", x)

# Add leading white space
add_leading_spaces <- function (x, y){
  length_of_asterisk <- nchar(x)
  length_of_sequence <- y
  leading_spaces <- paste(replicate(length_of_sequence - length_of_asterisk, " "), collapse = "")
  new_string <- paste(leading_spaces, x)
  new_string <- substr(new_string, 2, 61)
  return (new_string)
}
# Intiailizes dataframe with given number of strands (rows) for NxN matrix
# and strand names (to change the names of rows)
initialize_df <- function (number_of_strands, strand_names){
  df <- data.frame(matrix(ncol = number_of_strands, 
                          nrow = number_of_strands), 
                   row.names = strand_names)
  colnames(df) <- strand_names
  return (df)
}
# Gets indices of columns of metrics for the entire sequence, CDS and nonCDS
# e.g. Similarity returns 1, 10, 19
get_selected_columns <- function(x){
  if (x < 10){
    y <- x + 9
    z <- y + 9
  }
  else if (x < 19){
    x <- x - 9
    y <- x + 9
    z <- y + 9
  }
  else {
    x <- x - 9 - 9
    y <- x + 9
    z <- y + 9
  }
  return (c(x,y,z))
}
# Initializes x dataframes (empty dataframes) 
# with dimension x dimension dimensions, with row_names as row names 
initialize_dataframes <- function(dimension, row_names, x) {
  dataframes <- c()

  for (i in seq(1,x)){
    dataframe <- list(initialize_df(dimension, row_names))
    dataframes <- c(dataframes, dataframe)
  }
  return (dataframes)
}
##### FILLS THE METRICS

fill_metrics <- function(metrics, df_clustal, indices_with_mutations, sequence){
  number_of_strands <- nrow(metrics[[1]])
  for (i in seq(1, number_of_strands)){
    sequence1 <- df_clustal[i, 2]
    length <- nchar(sequence1)
    for (j in seq(1, number_of_strands)){
      sequence2 <- df_clustal[j, 2]
      mutations <- 0
      transitions <- 0
      transversions <- 0
      gaps <- 0
      insertions <- 0
      deletions <- 0
      Ns <- 0
      # iterate over all characters in a strand, all strands have same length
      for (k in sequence){
        if (!(k %in% indices_with_mutations)){
          next
        }
        if (substr(sequence1, k, k) == 'N' || substr(sequence2, k, k) == 'N'){
          Ns <- Ns + 1
        }
        if(substr(sequence1, k, k) != substr(sequence2, k, k)){
          mutations <- mutations + 1
          if (substr(sequence1, k, k) == '-'){
            gaps <- gaps + 1
            insertions <- insertions + 1
          }
          else if(substr(sequence2, k, k) == '-'){
            gaps <- gaps + 1
            deletions <- deletions + 1
          }
          else if ((substr(sequence1, k, k) == 'A' && substr(sequence2, k, k) == 'G') || 
                   (substr(sequence1, k, k) == 'G' && substr(sequence2, k, k) == 'A')){
            transitions <- transitions + 1
          }
          else if ((substr(sequence1, k, k) == 'T' && substr(sequence2, k, k) == 'C') || 
                   (substr(sequence1, k, k) == 'C' && substr(sequence2, k, k) == 'T')){
            transitions <- transitions + 1
          }
          else{
            transversions <- transversions + 1
          }
        }
      }
      metrics[[1]][i,j] <- 1 - (mutations/length)
      metrics[[2]][i,j] <- mutations
      metrics[[3]][i,j] <- metrics[[2]][i,j]/length
      metrics[[4]][i,j] <- transitions
      metrics[[5]][i,j] <- transversions
      if(transversions != 0){
        metrics[[6]][i,j] <- transitions/transversions
      }
      else{
        metrics[[6]][i,j] <- 0
      }
      metrics[[7]][i,j] <- gaps
      metrics[[8]][i,j] <- insertions
      metrics[[9]][i,j] <- deletions
      # MSA_similarity[i, j] <- 1 - (mutations/length)
      # MSA_mutation_num[i, j] <- mutations
      # MSA_transition_num[i, j] <- transitions
      # MSA_transversion_num [i, j] <- transversions
      # if(transversions != 0){
      #   MSA_tt_ratio[i, j] <- transitions/transversions
      # }
      # else{
      #   MSA_tt_ratio[i, j] <- 0
      # }
      # MSA_gaps [i, j] <- gaps
      # MSA_insertions[i,j] <- insertions
      # MSA_deletions[i, j] <- deletions
    }
  }
  return (metrics)
}

# Creates and fills summary pairwise dataframe
fill_summary_pairwise <- function(MSA_metrics, CDS_metrics, nonCDS_metrics, number_of_strands, strand_names){
  summary_pairwise <- data.frame(matrix(ncol = 27, nrow = number_of_strands-1), #except first strand 
                                 row.names = strand_names[2:length(strand_names)])
  for (i in seq(1,length(strand_names)-1)){
    summary_pairwise[i, 1] = MSA_metrics[[1]][i+1, i]
    summary_pairwise[i, 2] = MSA_metrics[[2]][i+1, i]
    summary_pairwise[i, 3] = MSA_metrics[[3]][i+1, i]
    summary_pairwise[i, 4] = MSA_metrics[[4]][i+1, i]
    summary_pairwise[i, 5] = MSA_metrics[[5]][i+1, i]
    summary_pairwise[i, 6] = MSA_metrics[[6]][i+1, i]
    summary_pairwise[i, 7] = MSA_metrics[[7]][i+1, i]
    summary_pairwise[i, 8] = MSA_metrics[[8]][i+1, i]
    summary_pairwise[i, 9] = MSA_metrics[[9]][i+1, i]
    summary_pairwise[i, 10] = CDS_metrics[[1]][i+1, i]
    summary_pairwise[i, 11] = CDS_metrics[[2]][i+1, i]
    summary_pairwise[i, 12] = CDS_metrics[[3]][i+1, i]
    summary_pairwise[i, 13] = CDS_metrics[[4]][i+1, i]
    summary_pairwise[i, 14] = CDS_metrics[[5]][i+1, i]
    summary_pairwise[i, 15] = CDS_metrics[[6]][i+1, i]
    summary_pairwise[i, 16] = CDS_metrics[[7]][i+1, i]
    summary_pairwise[i, 17] = CDS_metrics[[8]][i+1, i]
    summary_pairwise[i, 18] = CDS_metrics[[9]][i+1, i]
    summary_pairwise[i, 19] = nonCDS_metrics[[1]][i+1, i]
    summary_pairwise[i, 20] = nonCDS_metrics[[2]][i+1, i]
    summary_pairwise[i, 21] = nonCDS_metrics[[3]][i+1, i]
    summary_pairwise[i, 22] = nonCDS_metrics[[4]][i+1, i]
    summary_pairwise[i, 23] = nonCDS_metrics[[5]][i+1, i]
    summary_pairwise[i, 24] = nonCDS_metrics[[6]][i+1, i]
    summary_pairwise[i, 25] = nonCDS_metrics[[7]][i+1, i]
    summary_pairwise[i, 26] = nonCDS_metrics[[8]][i+1, i]
    summary_pairwise[i, 27] = nonCDS_metrics[[9]][i+1, i]
  }
  return(summary_pairwise)
}
##### VARIABLES

metrics <- c('Similarity', 'Mutations', 'Mutation_Frequency', 'Transitions', 
             'Transversions', 'TT_ratio', 'Gaps', 'Insertions', 'Deletions',
             'CDS_Similarity', 'CDS_Mutations', 'CDS_Mutation_Frequency', 
             'CDS_Transitions', 'CDS_Transversions', 'CDS_TT_ratio',
             'CDS_Gaps', 'CDS_Insertions', 'CDS_Deletions',
             'nonCDS_Similarity', 'nonCDS_Mutations', 
             'nonCDS_Mutation_Frequency', 'nonCDS_Transitions', 
             'nonCDS_Transversions', 'nonCDS_TT_ratio', 'nonCDS_Gaps', 
             'nonCDS_Insertions', 'nonCDS_Deletions')



