library(shiny)
library(stringr) #needed to split columns in dataframe
library(dplyr) 

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
initialize_df <- function (number_of_strands, strand_names){
  df <- data.frame(matrix(ncol = number_of_strands, 
                          nrow = number_of_strands), 
                   row.names = strand_names)
  colnames(df) <- strand_names
  return (df)
}

ui <- fluidPage(
  sliderInput(inputId='num',
              label = 'Choose a number',
              value = 25, min = 1, max = 100),
  fileInput(inputId = 'clustalFile', 
            label = 'Choose a .clustal file',
            accept = '.clustal'),
  fileInput(inputId = 'csvFile', 
            label = 'Choose a .csv file',
            accept = '.csv'),
  selectInput(inputId = 'metric', 
              label = 'Select the metric',
              choices =
                c('Similarity', 'Mutations')),
  actionButton(inputId ='update', label ='Process'),
  plotOutput("plot"),
  plotOutput("plot2")
  )

server <- function (input, output) {
  setwd('C:\\Users\\fejsa\\OneDrive\\Desktop\\Graduation project\\COMSAA\\input')
  #Load .csv file and store it into start_stop reactive expression
  
  # Processes the .clustal file and .csv file when button is clicked
  observeEvent(input$update, {
    start_stop <- reactive({
      
    })
    df_clustal <- reactive({
      
      
    })
    
  })
  
  #Load .clustal file and store it into start_stop reactive expression
  
  summary1 <- eventReactive(input$update, {
    #Load .csv file
    csvFile <- input$csvFile[1,1]
    start_stop <- read.csv(csvFile)
    #Load .clustal file
    clustalFile <- input$clustalFile[1,1]
    df_clustal <- read.delim(clustalFile) #first_df equivalent
    
    #Getting the name of the virus
    name_of_virus <- unlist(strsplit(csvFile, split='.', fixed=TRUE))[1]
    #Splitting rows on strands and sequences
    df_clustal[c('strand', 'sequence')] <- 
      str_split_fixed(df_clustal$CLUSTAL.O.1.2.4..multiple.sequence.alignment, ' ', 2)
    # Changing the name of the rows with only asterisks
    df_clustal$strand[df_clustal$strand == ''] <- 'asterisk'
    
    #drop first column that has everything combined
    df_clustal <- subset(df_clustal, select = -c(1))
    #clean leading spaces
    df_clustal$sequence <- lapply(df_clustal$sequence, trim_leading)
    
    # calculate number of characters for each line, 
    # all except the last one should be of length 60
    first_strand <- df_clustal[1,1]
    length_of_lines <- 
      lapply(df_clustal[df_clustal$strand==first_strand,2],nchar)
    length_of_sequence <- nchar(df_clustal$sequence[1])
    
    asterisk_only <- df_clustal$sequence[df_clustal$strand=='asterisk']
    
    # Fixing the length of rows with asterisks
    # All rows have the same length as the rows with sequences
    df_clustal$sequence[df_clustal$strand=='asterisk'] <- 
      mapply (add_leading_spaces, asterisk_only, length_of_lines)
    
    # groups sequences based on strand
    # strands are now in one line
    df_clustal <- df_clustal %>%
      group_by(strand) %>%
      summarise(sequence = paste(sequence, collapse = ""))
    
    #indices with mutations
    indices_with_mutations <- 
      which(strsplit(df_clustal$sequence[df_clustal$strand == 'asterisk'], "")[[1]] == " ") 
    #indices with no mutations
    indices_without_mutations <- 
      which(strsplit(df_clustal$sequence[df_clustal$strand == 'asterisk'], "")[[1]] == "*") 
    
    # extracting asterisk row
    asterisk_df <- strsplit(df_clustal$sequence[1], split='')
    df_clustal <- subset(df_clustal,strand != "asterisk" )
    
    # Get numbers of strands; #e.g. seq1997 => 1997
    mutation_numbers <- as.numeric(unlist(lapply(df_clustal$strand, get_ordered)))
    # Orders strands in ascending order
    df_clustal <- df_clustal[order(mutation_numbers),]
    
    # Get names and number of strands
    strand_names <- df_clustal$strand
    number_of_strands <- nrow(df_clustal)
    
    # Initializing dataframes for MSA
    # percent identity matrix
    MSA_similarity <- initialize_df(number_of_strands, strand_names)
    # mutation number matrix
    MSA_mutation_num <- initialize_df(number_of_strands, strand_names)
    # transitions number matrix
    MSA_transition_num <- initialize_df(number_of_strands, strand_names)
    # transversions number matrix
    MSA_transversion_num <- initialize_df(number_of_strands, strand_names)
    # transition/transvertion ratio matrix
    MSA_tt_ratio <- initialize_df(number_of_strands, strand_names)
    # number of gaps matrix
    MSA_gaps = initialize_df(number_of_strands, strand_names)
    # number of insertions matrix 
    MSA_insertions <- initialize_df(number_of_strands, strand_names)
    # number of deletions matrix
    MSA_deletions <- initialize_df(number_of_strands, strand_names)
    
    
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
        for (k in indices_with_mutations){
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
        MSA_similarity[i, j] <- 1 - (mutations/length)
        MSA_mutation_num[i, j] <- mutations
        MSA_transition_num[i, j] <- transitions
        MSA_transversion_num [i, j] <- transversions
        if(transversions != 0){
          MSA_tt_ratio[i, j] <- transitions/transversions
        }
        else{
          MSA_tt_ratio[i, j] <- 0
        }
        MSA_gaps [i, j] <- gaps
        MSA_insertions[i,j] <- insertions
        MSA_deletions[i, j] <- deletions
      }
    }
    ################################################################################
    ############################ CODING SEQUENCE ###################################
    ################################################################################
    
    start_stop <- subset(start_stop, select = c('Start', 'Stop'))
    start_stop <- start_stop[order(start_stop$Start),]
    
    #if there is only one row then mapply unlists the sequence
    if(nrow(start_stop) > 1){
      start_stop$Sequences <- mapply(create_sequence, start_stop$Start, start_stop$Stop)
    }
    else{
      start_stop$Sequences <- list(mapply(create_sequence, start_stop$Start, start_stop$Stop)) 
    }
    coding_sequence_index <- c()
    
    for (i in seq(1, nrow(start_stop))){
      coding_sequence_index <- append(coding_sequence_index, unlist(start_stop[i, ]$Sequences))
    }
    coding_sequence_index <- unlist(coding_sequence_index)
    coding_sequence_index <- unique(coding_sequence_index)
    
    noncoding_sequence_index <- setdiff(seq(1, nchar(df_clustal[1, 2])), coding_sequence_index)
    
    # percent identity matrix
    CDS_similarity <- initialize_df(number_of_strands, strand_names)
    # mutation number matrix
    CDS_mutation_num <- initialize_df(number_of_strands, strand_names)
    # transitions number matrix
    CDS_transition_num <- initialize_df(number_of_strands, strand_names)
    # transversions number matrix
    CDS_transversion_num <- initialize_df(number_of_strands, strand_names)
    # transition/transvertion ratio matrix
    CDS_tt_ratio <- initialize_df(number_of_strands, strand_names)
    # number of gaps matrix
    CDS_gaps = initialize_df(number_of_strands, strand_names)
    # number of insertions matrix 
    CDS_insertions <- initialize_df(number_of_strands, strand_names)
    # number of deletions matrix
    CDS_deletions <- initialize_df(number_of_strands, strand_names)
    length(coding_sequence_index)
    
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
        for (k in coding_sequence_index){
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
        CDS_similarity[i, j] <- 1 - (mutations/length)
        CDS_mutation_num[i, j] <- mutations
        CDS_transition_num[i, j] <- transitions
        CDS_transversion_num [i, j] <- transversions
        if(transversions != 0){
          CDS_tt_ratio[i, j] <- transitions/transversions
        }
        else{
          CDS_tt_ratio[i, j] <- 0
        }
        CDS_gaps [i, j] <- gaps
        CDS_insertions[i,j] <- insertions
        CDS_deletions[i, j] <- deletions
      }
    }
    
    ################################################################################
    ############################ NON CODING SEQUECE ################################
    ################################################################################
    
    # percent identity matrix
    nonCDS_similarity <- initialize_df(number_of_strands, strand_names)
    # mutation number matrix
    nonCDS_mutation_num <- initialize_df(number_of_strands, strand_names)
    # transitions number matrix
    nonCDS_transition_num <- initialize_df(number_of_strands, strand_names)
    # transversions number matrix
    nonCDS_transversion_num <- initialize_df(number_of_strands, strand_names)
    # transition/transvertion ratio matrix
    nonCDS_tt_ratio <- initialize_df(number_of_strands, strand_names)
    # number of gaps matrix
    nonCDS_gaps = initialize_df(number_of_strands, strand_names)
    # number of insertions matrix 
    nonCDS_insertions <- initialize_df(number_of_strands, strand_names)
    # number of deletions matrix
    nonCDS_deletions <- initialize_df(number_of_strands, strand_names)
    
    
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
        for (k in noncoding_sequence_index){
          if(!(k %in% indices_with_mutations))
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
        nonCDS_similarity[i, j] <- 1 - (mutations/length)
        nonCDS_mutation_num[i, j] <- mutations
        nonCDS_transition_num[i, j] <- transitions
        nonCDS_transversion_num [i, j] <- transversions
        if(transversions != 0){
          nonCDS_tt_ratio[i, j] <- transitions/transversions
        }
        else{
          nonCDS_tt_ratio[i, j] <- 0
        }
        nonCDS_gaps [i, j] <- gaps
        nonCDS_insertions[i,j] <- insertions
        nonCDS_deletions[i, j] <- deletions
      }
    }
    length_of_sequence <- length(indices_with_mutations) + 
      length(indices_without_mutations)
    
    #Getting Mutations frequency
    MSA_mutation_frequency <- MSA_mutation_num/length_of_sequence
    CDS_mutation_frequency <- CDS_mutation_num/length(coding_sequence_index)
    nonCDS_mutation_frequency <- nonCDS_mutation_num/length(noncoding_sequence_index)
    #Creating summary dataframe
    summary <- cbind(MSA_similarity[1], MSA_mutation_num[1], 
                     MSA_mutation_frequency[1],MSA_transition_num[1], 
                     MSA_transversion_num[1], MSA_tt_ratio[1], MSA_gaps[1], 
                     MSA_insertions[1], MSA_deletions[1],
                     CDS_similarity[1], CDS_mutation_num[1], 
                     CDS_mutation_frequency[1],CDS_transition_num[1], 
                     CDS_transversion_num[1], CDS_tt_ratio[1], CDS_gaps[1], 
                     CDS_insertions[1], CDS_deletions[1],
                     nonCDS_similarity[1], nonCDS_mutation_num[1], 
                     nonCDS_mutation_frequency[1],nonCDS_transition_num[1], 
                     nonCDS_transversion_num[1], nonCDS_tt_ratio[1], nonCDS_gaps[1], 
                     nonCDS_insertions[1], nonCDS_deletions[1])
    summary_column_names <- c('Similarity', "Mutations", "Mutation_Frequency", 
                              "Transitions", "Transversions", "TT_ratio", "Gaps", 
                              "Insertions", "Deletions",
                              'CDS_Similarity', "CDS_Mutations", "CDS_Mutation_Frequency", 
                              "CDS_Transitions", "CDS_Transversions", "CDS_TT_ratio", 
                              "CDS_Gaps", "CDS_Insertions", "CDS_Deletions",
                              'nonCDS_Similarity', "nonCDS_Mutations", 
                              "nonCDS_Mutation_Frequency", "nonCDS_Transitions", 
                              "nonCDS_Transversions", "nonCDS_TT_ratio", "nonCDS_Gaps", 
                              "nonCDS_Insertions", "nonCDS_Deletions") 
    colnames(summary) <- summary_column_names
    summary_pairwise <- data.frame(matrix(ncol = 27, nrow = number_of_strands-1), #except first strand 
                                   row.names = strand_names[2:length(strand_names)])
    colnames(summary_pairwise) <- summary_column_names
    ## Filling the summary_pairwise dataframe
    for (i in seq(1,length(strand_names)-1)){
      summary_pairwise[i, 1] = MSA_similarity[i+1, i]
      summary_pairwise[i, 2] = MSA_mutation_num[i+1, i]
      summary_pairwise[i, 3] = MSA_mutation_frequency[i+1, i]
      summary_pairwise[i, 4] = MSA_transition_num[i+1, i]
      summary_pairwise[i, 5] = MSA_transversion_num[i+1, i]
      summary_pairwise[i, 6] = MSA_tt_ratio[i+1, i]
      summary_pairwise[i, 7] = MSA_gaps[i+1, i]
      summary_pairwise[i, 8] = MSA_insertions[i+1, i]
      summary_pairwise[i, 9] = MSA_deletions[i+1, i]
      summary_pairwise[i, 10] = CDS_similarity[i+1, i]
      summary_pairwise[i, 11] = CDS_mutation_num[i+1, i]
      summary_pairwise[i, 12] = CDS_mutation_frequency[i+1, i]
      summary_pairwise[i, 13] = CDS_transition_num[i+1, i]
      summary_pairwise[i, 14] = CDS_transversion_num[i+1, i]
      summary_pairwise[i, 15] = CDS_tt_ratio[i+1, i]
      summary_pairwise[i, 16] = CDS_gaps[i+1, i]
      summary_pairwise[i, 17] = CDS_insertions[i+1, i]
      summary_pairwise[i, 18] = CDS_deletions[i+1, i]
      summary_pairwise[i, 19] = nonCDS_similarity[i+1, i]
      summary_pairwise[i, 20] = nonCDS_mutation_num[i+1, i]
      summary_pairwise[i, 21] = nonCDS_mutation_frequency[i+1, i]
      summary_pairwise[i, 22] = nonCDS_transition_num[i+1, i]
      summary_pairwise[i, 23] = nonCDS_transversion_num[i+1, i]
      summary_pairwise[i, 24] = nonCDS_tt_ratio[i+1, i]
      summary_pairwise[i, 25] = nonCDS_gaps[i+1, i]
      summary_pairwise[i, 26] = nonCDS_insertions[i+1, i]
      summary_pairwise[i, 27] = nonCDS_deletions[i+1, i]
    }
    View(summary)
    summary
  })
  
  output$plot <- renderPlot({
    metric <- input$metric
    sum <- summary1()
    hist(sum[,metric], breaks = input$num)
  })
}


shinyApp(ui = ui, server = server)