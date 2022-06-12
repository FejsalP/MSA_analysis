library(shiny)
library(stringr) # needed to split columns in dataframe
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(tidyr)
################################################################################
############################# FUNCTIONS ########################################
################################################################################

# Creates sequence based on start and stop values
create_sequence <- function(start, end) {
  return(seq(start, end))
}

# Returns the number of the strand,
# e.g. virus_2020 returns 2020
get_ordered <- function(name) {
  name <- strsplit(name, split = "_")
  return(name[[1]][2])
}
# Trims leading white space
trim_leading <- function(x) sub("^\\s+", "", x)

# Add leading white space
add_leading_spaces <- function(x, y) {
  length_of_asterisk <- nchar(x)
  length_of_sequence <- y
  leading_spaces <- paste(replicate(length_of_sequence - length_of_asterisk, " "), collapse = "")
  new_string <- paste(leading_spaces, x)
  new_string <- substr(new_string, 2, 61)
  return(new_string)
}
# Intiailizes dataframe with given number of strands (rows) for NxN matrix
# and strand names (to change the names of rows)
initialize_df <- function(number_of_strands, strand_names) {
  df <- data.frame(matrix(
    ncol = number_of_strands,
    nrow = number_of_strands
  ),
  row.names = strand_names
  )
  colnames(df) <- strand_names
  return(df)
}
# Gets indices of columns of metrics for the entire sequence, CDS and nonCDS
# e.g. Similarity returns 1, 10, 19
get_selected_columns <- function(x) {
  if (x < 10) {
    y <- x + 9
    z <- y + 9
  } else if (x < 19) {
    x <- x - 9
    y <- x + 9
    z <- y + 9
  } else {
    x <- x - 9 - 9
    y <- x + 9
    z <- y + 9
  }
  return(c(x, y, z))
}
# Initializes x dataframes (empty dataframes)
# with dimension x dimension dimensions, with row_names as row names
initialize_dataframes <- function(dimension, row_names, x) {
  dataframes <- c()

  for (i in seq(1, x)) {
    dataframe <- list(initialize_df(dimension, row_names))
    dataframes <- c(dataframes, dataframe)
  }
  return(dataframes)
}
##### FILLS THE METRICS

fill_metrics <- function(metrics, df_clustal, indices_with_mutations, sequence) {
  number_of_strands <- nrow(metrics[[1]])
  for (i in seq(1, number_of_strands)) {
    sequence1 <- df_clustal[i, 2]
    length <- nchar(sequence1)
    for (j in seq(1, number_of_strands)) {
      sequence2 <- df_clustal[j, 2]
      mutations <- 0
      transitions <- 0
      transversions <- 0
      gaps <- 0
      insertions <- 0
      deletions <- 0
      Ns <- 0
      # iterate over all characters in a strand, all strands have same length
      for (k in sequence) {
        if (!(k %in% indices_with_mutations)) {
          next
        }
        if (substr(sequence1, k, k) == "N" || substr(sequence2, k, k) == "N") {
          Ns <- Ns + 1
        }
        if (substr(sequence1, k, k) != substr(sequence2, k, k)) {
          mutations <- mutations + 1
          if (substr(sequence1, k, k) == "-") {
            gaps <- gaps + 1
            insertions <- insertions + 1
          } else if (substr(sequence2, k, k) == "-") {
            gaps <- gaps + 1
            deletions <- deletions + 1
          } else if ((substr(sequence1, k, k) == "A" && substr(sequence2, k, k) == "G") ||
            (substr(sequence1, k, k) == "G" && substr(sequence2, k, k) == "A")) {
            transitions <- transitions + 1
          } else if ((substr(sequence1, k, k) == "T" && substr(sequence2, k, k) == "C") ||
            (substr(sequence1, k, k) == "C" && substr(sequence2, k, k) == "T")) {
            transitions <- transitions + 1
          } else {
            transversions <- transversions + 1
          }
        }
      }
      metrics[[1]][i, j] <- 1 - (mutations / length)
      metrics[[2]][i, j] <- mutations
      metrics[[3]][i, j] <- metrics[[2]][i, j] / length
      metrics[[4]][i, j] <- transitions
      metrics[[5]][i, j] <- transversions
      if (transversions != 0) {
        metrics[[6]][i, j] <- transitions / transversions
      } else {
        metrics[[6]][i, j] <- 0
      }
      metrics[[7]][i, j] <- gaps
      metrics[[8]][i, j] <- insertions
      metrics[[9]][i, j] <- deletions
    }
  }
  return(metrics)
}

# Creates and fills summary pairwise dataframe
fill_summary_pairwise <- function(MSA_metrics, CDS_metrics, nonCDS_metrics, number_of_strands, strand_names) {
  summary_pairwise <- data.frame(matrix(ncol = 27, nrow = number_of_strands - 1), # except first strand
    row.names = strand_names[2:length(strand_names)]
  )
  for (i in seq(1, length(strand_names) - 1)) {
    summary_pairwise[i, 1] <- MSA_metrics[[1]][i + 1, i]
    summary_pairwise[i, 2] <- MSA_metrics[[2]][i + 1, i]
    summary_pairwise[i, 3] <- MSA_metrics[[3]][i + 1, i]
    summary_pairwise[i, 4] <- MSA_metrics[[4]][i + 1, i]
    summary_pairwise[i, 5] <- MSA_metrics[[5]][i + 1, i]
    summary_pairwise[i, 6] <- MSA_metrics[[6]][i + 1, i]
    summary_pairwise[i, 7] <- MSA_metrics[[7]][i + 1, i]
    summary_pairwise[i, 8] <- MSA_metrics[[8]][i + 1, i]
    summary_pairwise[i, 9] <- MSA_metrics[[9]][i + 1, i]
    summary_pairwise[i, 10] <- CDS_metrics[[1]][i + 1, i]
    summary_pairwise[i, 11] <- CDS_metrics[[2]][i + 1, i]
    summary_pairwise[i, 12] <- CDS_metrics[[3]][i + 1, i]
    summary_pairwise[i, 13] <- CDS_metrics[[4]][i + 1, i]
    summary_pairwise[i, 14] <- CDS_metrics[[5]][i + 1, i]
    summary_pairwise[i, 15] <- CDS_metrics[[6]][i + 1, i]
    summary_pairwise[i, 16] <- CDS_metrics[[7]][i + 1, i]
    summary_pairwise[i, 17] <- CDS_metrics[[8]][i + 1, i]
    summary_pairwise[i, 18] <- CDS_metrics[[9]][i + 1, i]
    summary_pairwise[i, 19] <- nonCDS_metrics[[1]][i + 1, i]
    summary_pairwise[i, 20] <- nonCDS_metrics[[2]][i + 1, i]
    summary_pairwise[i, 21] <- nonCDS_metrics[[3]][i + 1, i]
    summary_pairwise[i, 22] <- nonCDS_metrics[[4]][i + 1, i]
    summary_pairwise[i, 23] <- nonCDS_metrics[[5]][i + 1, i]
    summary_pairwise[i, 24] <- nonCDS_metrics[[6]][i + 1, i]
    summary_pairwise[i, 25] <- nonCDS_metrics[[7]][i + 1, i]
    summary_pairwise[i, 26] <- nonCDS_metrics[[8]][i + 1, i]
    summary_pairwise[i, 27] <- nonCDS_metrics[[9]][i + 1, i]
  }
  return(summary_pairwise)
}

# Make a reference list, which will contain all mutation numbers at every point in the MSA. E.g., if in
# MSA index 5, 3 sequences have a mutation, index five in the reference list will contain the number 3.
# Ideally, five reference lists will be made, two containing transitions and transversions, one containing
# point mutations (sum of the previous two), one containing Gaps, and the fifth being the sum of the
# previous point mutations and gaps. This can be useful in order to observe different types of
# mutations individually, or all mutations together.
create_reference_list <- function(df, indices_with_mutations) {
  reference_list <- c()
  original_sequence <- df$sequence[1]
  reference_list <- rep(0, nchar(original_sequence))
  number_of_strands <- nrow(df)
  for (i in indices_with_mutations) {
    count <- 0
    for (j in seq(2, number_of_strands)) {
      if (substr(df$sequence[1], i, i) != substr(df$sequence[j], i, i)) {
        count <- count + 1
      }
    }
    reference_list[i] <- count
  }
  return(as.data.frame(reference_list))
}
# Reference list for transitions
# the cases where A > G, G > A, T > C, or C > T; > is symbol for mutation
# A > G
# G > A
# T > C
# C > T
create_reference_list_transitions <- function(df, indices_with_mutations) {
  reference_list <- c()
  original_sequence <- df$sequence[1]
  reference_list <- rep(0, nchar(original_sequence))
  number_of_strands <- nrow(df)

  for (i in indices_with_mutations) {
    count <- 0
    for (j in seq(2, number_of_strands)) {
      if ((substr(df$sequence[1], i, i) == "A" & substr(df$sequence[j], i, i) == "G") |
        (substr(df$sequence[1], i, i) == "G" & substr(df$sequence[j], i, i) == "A") |
        (substr(df$sequence[1], i, i) == "T" & substr(df$sequence[j], i, i) == "C") |
        (substr(df$sequence[1], i, i) == "C" & substr(df$sequence[j], i, i) == "T")
      ) {
        count <- count + 1
      }
    }
    reference_list[i] <- count
  }
  return(as.data.frame(reference_list))
}
# Reference list for transversions
# (the cases where A > T or C, G > T or C, T > A or G, or C > A or G)
# A > T or C
# G > T or C
# T > A or G
# C > A or G

create_reference_list_transversions <- function(df, indices_with_mutations) {
  reference_list <- c()
  original_sequence <- df$sequence[1]
  reference_list <- rep(0, nchar(original_sequence))
  number_of_strands <- nrow(df)

  for (i in indices_with_mutations) {
    count <- 0
    for (j in seq(2, number_of_strands)) {
      if ((substr(df$sequence[1], i, i) == "A" & (substr(df$sequence[j], i, i) == "T" | substr(df$sequence[j], i, i) == "C")) |
        (substr(df$sequence[1], i, i) == "G" & (substr(df$sequence[j], i, i) == "T" | substr(df$sequence[j], i, i) == "C")) |
        (substr(df$sequence[1], i, i) == "T" & (substr(df$sequence[j], i, i) == "A" | substr(df$sequence[j], i, i) == "G")) |
        (substr(df$sequence[1], i, i) == "C" & (substr(df$sequence[j], i, i) == "A" | substr(df$sequence[j], i, i) == "G"))
      ) {
        count <- count + 1
      }
    }
    reference_list[i] <- count
  }
  return(as.data.frame(reference_list))
}

# Reference list for gaps
# a gap existing in one letter, where there is a letter in another sequence

create_reference_list_gaps <- function(df, indices_with_mutations) {
  reference_list <- c()
  original_sequence <- df$sequence[1]
  reference_list <- rep(0, nchar(original_sequence))
  number_of_strands <- nrow(df)

  for (i in indices_with_mutations) {
    count <- 0
    for (j in seq(2, number_of_strands)) {
      if ((substr(df$sequence[1], i, i) == "-" & substr(df$sequence[j], i, i) != "-") | # insertion?
        (substr(df$sequence[1], i, i) != "-" & substr(df$sequence[j], i, i) == "-") # deletion
      ) {
        count <- count + 1
      }
    }
    reference_list[i] <- count
  }
  return(as.data.frame(reference_list))
}


create_reference_list_grouped <- function(reference_list, reference_groups) {
  reference_list$group <- c(0, rep(1:(nrow(reference_list) - 1) %/% reference_groups))

  grouped_reference_list <- group_by(reference_list, group) %>%
    summarise(number_of_mutations = sum(reference_list))

  return(grouped_reference_list)
}
# Given vector of numbers, return a list of consecutive numbers
# 1 2 3 4 9 10 15 16 - returns ((1,2,3,4), (9, 10), (15, 16))
create_df_ranges <- function(sequence) {
  list_of_sequences <- list()
  ranges <- list()
  for (i in seq(1, length(sequence) - 1)) {
    if ((sequence[i + 1] - sequence[i]) == 1) {
      ranges <- append(ranges, sequence[i])
      if (i != (length(sequence) - 1)) {
        next
      }
    }
    ranges <- append(ranges, sequence[i])
    if (i == (length(sequence) - 1)) {
      ranges <- append(ranges, sequence[i + 1])
    }
    list_of_sequences <- append(list_of_sequences, list(ranges))
    ranges <- list()
  }
  return(list_of_sequences)
}

get_added_number_of_mutations <- function(start, end, reference_list) {
  number_of_mutations <- c()
  for (i in seq(1, length(start))) {
    number_of_mutations <- append(number_of_mutations, sum(reference_list[start[i]:end[i]]))
  }
  return(number_of_mutations)
}

# Creates a reference list that has separated coding and noncoding sequence
# Given coding and noncoding sequences return a reference list
# E.g. nonCDS1 1:300 - 50 mutations, CDS1 301:500, nonCDS2 501:700, CDS2 701:900..
create_reference_list_separated <- function(coding_sequence_index, noncoding_sequence_index, reference_list) {
  coding_sequence_indices_df <- create_df_ranges(coding_sequence_index)
  noncoding_sequence_indices_df <- create_df_ranges(noncoding_sequence_index)

  # Get starting indices of coding and noncoding sequences
  coding_sequence_starting_indices <- lapply(coding_sequence_indices_df, function(x) {
    return(x[[1]])
  })
  noncoding_sequence_starting_indices <- lapply(noncoding_sequence_indices_df, function(x) {
    return(x[[1]])
  })
  # Get ending indices of coding and noncoding sequences
  coding_sequence_ending_indices <- lapply(coding_sequence_indices_df, function(x) {
    return(x[length(x)])
  })
  noncoding_sequence_ending_indices <- lapply(noncoding_sequence_indices_df, function(x) {
    return(x[length(x)])
  })

  # Store into a dataframe, also add Type column to reference the number of CDS/nonCDS portion
  coding_sequence_start_end_df <- data.frame(
    "Start" = unlist(coding_sequence_starting_indices),
    "End" = unlist(coding_sequence_ending_indices),
    "Type" = paste0(rep("CDS_", length(coding_sequence_starting_indices)), seq(1, length(coding_sequence_starting_indices)))
  )
  noncoding_sequence_start_end_df <- data.frame(
    "Start" = unlist(noncoding_sequence_starting_indices),
    "End" = unlist(noncoding_sequence_ending_indices),
    "Type" = paste0(rep("nonCDS_", length(noncoding_sequence_starting_indices)), seq(1, length(noncoding_sequence_starting_indices)))
  )

  combined_sequences <- rbind(coding_sequence_start_end_df, noncoding_sequence_start_end_df)
  combined_sequences <- combined_sequences[order(combined_sequences$Start), ]
  combined_sequences["Number_of_mutations"] <- get_added_number_of_mutations(combined_sequences$Start, combined_sequences$End, reference_list)
  # Used to reorder barplot
  combined_sequences["Count"] <- seq(1, nrow(combined_sequences))
  # Drops start and end column, redundant
  combined_sequences <- subset(combined_sequences, select = -c(1, 2))
  return(combined_sequences)
}
# create_bar_plot <- function(df, title, selected_columns, selected_names_of_columns){
#
#   ggplot(df,
#          aes(x = Strands,
#              y = df[,selected_columns],
#              fill= Strands)) +
#     geom_bar(stat = 'identity',
#              position = 'dodge') +
#     labs(title = title, x = 'Strands', y = colnames(df)[which(colnames(df) == selected_names_of_columns)]) +
#     theme(axis.ticks.x = element_blank(),
#           axis.text.x = element_blank()) +
#     scale_fill_brewer(palette="Paired")
# }
# Creates barplot for strands

create_bar_plot <- function(df,
                            selected_columns,
                            selected_names_of_columns,
                            title,
                            x_axis,
                            y_axis,
                            title_size,
                            x_axis_size,
                            type) {
  ggplot(
    df,
    aes(
      x = Strands,
      y = df[, selected_columns],
      fill = Strands
    )
  ) +
    geom_bar(
      stat = "identity",
      position = "dodge"
    ) +
    theme_light() +
    labs(title = title, x = x_axis, y = paste(y_axis, " - ", type)) + #  bilo  title = title, x = 'Strands', y = colnames(df)[which(colnames(df) == selected_names_of_columns)]
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = x_axis_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5)
    ) +
    # theme(axis.ticks.x = element_blank(),
    #       axis.text.x = element_blank()) +
    scale_fill_brewer(palette = "Paired")
}

# Creates histogram for reference lists, ungrouped

create_reference_list_histogram <- function(values,
                                            breaks,
                                            title,
                                            x_axis,
                                            x_axis_size,
                                            y_axis_size,
                                            title_size,
                                            y_axis_name) {
  reference_list_grouped <- create_reference_list_grouped(values, breaks)
  histogram <- ggplot(
    reference_list_grouped,
    aes(
      x = group,
      y = number_of_mutations,
      fill = "red"
    )
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_light() +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = y_axis_size, face = "bold"),
      axis.title.x = element_text(size = x_axis_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5)
    ) +
    labs(title = title, x = x_axis, y = paste("Number of", y_axis_name))


  return(histogram)
}

create_reference_list_separated_histogram <- function(values_df,
                                                      title,
                                                      x_axis,
                                                      x_axis_size,
                                                      y_axis_size,
                                                      title_size,
                                                      x_axis_ticks_size,
                                                      y_axis_name) {
  values_df <- transform(values_df, group = ifelse(substr(Type, 1, 3) == "non", "nonCDS", "CDS"))
  ggplot(
    values_df,
    aes(
      x = reorder(Type, Count),
      y = Number_of_mutations,
      fill = group
    )
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_light() +
    theme(
      axis.text.x = element_text(size = x_axis_ticks_size),
      axis.title.x = element_text(size = x_axis_size, face = "bold"),
      axis.title.y = element_text(size = y_axis_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    ) +
    labs(title = title, x = x_axis, y = paste("Number of", y_axis_name))
}

create_mutations_type_bar_plot <- function(transitions,
                                           transversions,
                                           gaps,
                                           title,
                                           x_axis,
                                           x_axis_size,
                                           y_axis_size,
                                           title_size,
                                           x_axis_ticks_size) {
  mutations_types <- cbind(transitions, transversions, gaps)

  # count <- subset(mutations_types, select = c(3))
  mutations_types <- subset(mutations_types, select = c(1, 2, 5, 8, 3))
  types <- mutations_types$Type
  mutations_types <- unite(mutations_types, Type, c(Type, Count), sep = "-")

  mutations_types <- subset(mutations_types, select = -c(5))


  colnames(mutations_types) <- c("Type", "Transitions", "Transversions", "Gaps")


  mutations_types <- mutations_types %>%
    pivot_longer(!Type, names_to = "Mutations", values_to = "Mutations_count")


  mutations_types <- separate(mutations_types, col = Type, into = c("Type", "Count"), sep = "-")
  mutations_types$Type <- as.factor(mutations_types$Type)
  levels(mutations_types$Type) <- types
  ggplot(mutations_types, aes(x = Type, y = Mutations_count, fill = Mutations)) +
    geom_bar(stat = "identity") +
    theme_light() +
    labs(title = title, x = x_axis, y = "Number of mutations") +
    theme(
      axis.text.x = element_text(size = x_axis_ticks_size),
      axis.title.x = element_text(size = x_axis_size, face = "bold"),
      axis.title.y = element_text(size = y_axis_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    )
}


new_CDI <- function(coding_indices, sequence) {

  # to be returned
  # stores gaps as 0, starts counting from the first nucleobase


  coding_sequence_indices_df <- create_df_ranges(coding_indices)
  # Get starting indices
  coding_sequence_starting_indices <- lapply(coding_sequence_indices_df, function(x) {
    return(x[[1]])
  })
  # Get ending indices of coding
  coding_sequence_ending_indices <- lapply(coding_sequence_indices_df, function(x) {
    return(x[length(x)])
  })

  # Store into a dataframe
  coding_sequence_start_end_df <- data.frame(
    "Start" = unlist(coding_sequence_starting_indices),
    "Stop" = unlist(coding_sequence_ending_indices)
  )

  new_coding_sequence_start_end <- data.frame(matrix(ncol = 2, nrow = ))
  coding_indices_2 <- c()
  count <- 1
  for (i in seq(1, nchar(sequence))) {
    if (substr(sequence, i, i) == "-") {
      coding_indices_2 <- append(coding_indices_2, 0)
    } else {
      coding_indices_2 <- append(coding_indices_2, count)
      count <- count + 1
    }
  }
  for (i in seq(1, nrow(coding_sequence_start_end_df))) {
    coding_sequence_start_end_df[i, 1] <- which(coding_sequence_start_end_df[i, 1] == coding_indices_2)[[1]]
    coding_sequence_start_end_df[i, 2] <- which(coding_sequence_start_end_df[i, 2] == coding_indices_2)[[1]]
  }

  # Creates Sequences column that has sequence from start to end for each row
  if (nrow(coding_sequence_start_end_df) > 1) {
    coding_sequence_start_end_df$Sequences <- mapply(create_sequence, coding_sequence_start_end_df$Start, coding_sequence_start_end_df$Stop)
  } else {
    coding_sequence_start_end_df$Sequences <- list(mapply(create_sequence, coding_sequence_start_end_df$Start, coding_sequence_start_end_df$Stop))
  }
  # Creates new coding sequence index
  coding_sequence_index <- c()

  for (i in seq(1, nrow(coding_sequence_start_end_df))) {
    coding_sequence_index <- append(coding_sequence_index, unlist(coding_sequence_start_end_df[i, ]$Sequences))
  }
  coding_sequence_index <- unlist(coding_sequence_index)
  coding_sequence_index <- unique(coding_sequence_index)

  return(coding_sequence_index)
}

##### VARIABLES

metrics <- c(
  "Similarity", "Mutations", "Mutation_Frequency", "Transitions",
  "Transversions", "TT_ratio", "Gaps", "Insertions", "Deletions",
  "CDS_Similarity", "CDS_Mutations", "CDS_Mutation_Frequency",
  "CDS_Transitions", "CDS_Transversions", "CDS_TT_ratio",
  "CDS_Gaps", "CDS_Insertions", "CDS_Deletions",
  "nonCDS_Similarity", "nonCDS_Mutations",
  "nonCDS_Mutation_Frequency", "nonCDS_Transitions",
  "nonCDS_Transversions", "nonCDS_TT_ratio", "nonCDS_Gaps",
  "nonCDS_Insertions", "nonCDS_Deletions"
)
