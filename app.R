source('variables_and_functions_and_libraries.R')

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
              choices = metrics),
  actionButton(inputId ='update', label ='Process'),
  textInput(inputId='title', label = 'Title'),
  radioButtons(inputId='plot_type', label = 'Select the type of plot',
              choices=c('Histogram', 'Bar chart')), #violin, boxplot, dotplot, stem and leaf?
                                                   #scatter plot?
  plotOutput('plot', width='1200px', height='500px'),
  )

server <- function (input, output) {  #Load .clustal file and store it into start_stop reactive expression
  
  summary1 <- eventReactive(input$update, {
    #Load .csv file
    start_stop <- read.csv(input$csvFile$datapath)
    #Load .clustal file
    df_clustal <- read.delim(input$clustalFile$datapath) #first_df equivalent
    
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
      summarise(sequence = paste(sequence, collapse = ''))
    
    #indices with mutations
    indices_with_mutations <- 
      which(strsplit(df_clustal$sequence[df_clustal$strand == 'asterisk'], '')[[1]] == ' ') 
    #indices with no mutations
    indices_without_mutations <- 
      which(strsplit(df_clustal$sequence[df_clustal$strand == 'asterisk'], '')[[1]] == '*') 
    
    # extracting asterisk row
    asterisk_df <- strsplit(df_clustal$sequence[1], split='')
    df_clustal <- subset(df_clustal,strand != 'asterisk' )
    
    # Get numbers of strands; #e.g. seq1997 => 1997
    mutation_numbers <- as.numeric(unlist(lapply(df_clustal$strand, get_ordered)))
    # Orders strands in ascending order
    df_clustal <- df_clustal[order(mutation_numbers),]
    
    # Get names and number of strands
    strand_names <- df_clustal$strand
    number_of_strands <- nrow(df_clustal)
    
    length_of_sequence <- length(indices_with_mutations) + 
      length(indices_without_mutations)
    MSA_metrics <- initialize_dataframes(number_of_strands, strand_names, 9)
    MSA_metrics <- fill_metrics (MSA_metrics, df_clustal, indices_with_mutations, seq(1, length_of_sequence))
    
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
    
    CDS_metrics <- initialize_dataframes(number_of_strands, strand_names, 9)
    CDS_metrics <- fill_metrics (MSA_metrics, df_clustal, indices_with_mutations, 
                                 coding_sequence_index)
    
    ################################################################################
    ############################ NON CODING SEQUECE ################################
    ################################################################################
    
    nonCDS_metrics <- initialize_dataframes(number_of_strands, strand_names, 9)
    nonCDS_metrics <- fill_metrics (MSA_metrics, df_clustal, 
                                    indices_with_mutations, noncoding_sequence_index)
    
    #Creating summary dataframe
    summary <- cbind(MSA_metrics[[1]][,1], MSA_metrics[[2]][,1], 
                     MSA_metrics[[3]][,1], MSA_metrics[[4]][,1], 
                     MSA_metrics[[5]][,1], MSA_metrics[[6]][,1],
                     MSA_metrics[[7]][,1], MSA_metrics[[8]][,1], 
                     MSA_metrics[[9]][,1],
                     CDS_metrics[[1]][,1], CDS_metrics[[2]][,1], 
                     CDS_metrics[[3]][,1], CDS_metrics[[4]][,1], 
                     CDS_metrics[[5]][,1], CDS_metrics[[6]][,1],
                     CDS_metrics[[7]][,1], CDS_metrics[[8]][,1], 
                     CDS_metrics[[9]][,1],
                     nonCDS_metrics[[1]][,1], nonCDS_metrics[[2]][,1], 
                     nonCDS_metrics[[3]][,1], nonCDS_metrics[[4]][,1], 
                     nonCDS_metrics[[5]][,1], nonCDS_metrics[[6]][,1],
                     nonCDS_metrics[[7]][,1], nonCDS_metrics[[8]][,1], 
                     nonCDS_metrics[[9]][,1])
    summary_column_names <- metrics
    colnames(summary) <- summary_column_names
    
    ## Create and fill the summary_pairwise dataframe
    summary_pairwise <- fill_summary_pairwise(MSA_metrics, CDS_metrics, nonCDS_metrics, number_of_strands, strand_names)
    
    summary <- cbind(summary, rownames(summary))
    colnames(summary)[colnames(summary) == 'rownames(summary)'] <- 'strands'
    summary
  })
  output$plot <- renderPlot({
    metric <- input$metric
    sum <- summary1()
    if (input$plot_type == 'Histogram'){
      hist(sum[,metric], 
           breaks = input$num, 
           main = isolate(input$title))
    }
    else {
      selected_column <- which(colnames(sum) == metric)
      selected_columns <- get_selected_columns(selected_column)
      selected_names_of_columns <- colnames(sum)[selected_columns]
      new_sum <- as.matrix(sum[,selected_columns[c(1,2,3)]])
      colors <- rainbow(nrow(sum))

      entire_sequence <- 
        ggplot(sum,                                      # Grouped barplot using ggplot2
               aes(x = strands,
                   y = sum[,selected_columns[1]],
                   fill= strands)) +
        geom_bar(stat = 'identity',
                 position = 'dodge') + 
        labs(title = 'Title', x = 'Strands', y = colnames(sum)[which(colnames(sum) == selected_names_of_columns[1])]) +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) + 
        scale_fill_brewer(palette="Paired")
      coding_sequence <- 
        ggplot(sum,                                      # Grouped barplot using ggplot2
               aes(x = strands,
                   y = sum[,selected_columns[2]],
                   fill= strands)) +
        geom_bar(stat = 'identity',
                 position = 'dodge') + 
        labs(title = 'Title', x = 'Strands', y = colnames(sum)[which(colnames(sum) == selected_names_of_columns[2])]) +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) +
        scale_fill_brewer(palette="Paired")
      noncoding_sequence <- 
        ggplot(sum,                                      # Grouped barplot using ggplot2
               aes(x = strands,
                   y = sum[,selected_columns[3]],
                   fill= strands)) +
        geom_bar(stat = 'identity',
                 position = 'dodge') + 
        labs(title = 'Title', x = 'Strands', y = colnames(sum)[which(colnames(sum) == selected_names_of_columns[3])]) +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) +
        scale_fill_brewer(palette="Paired")
      
      ggarrange(entire_sequence, coding_sequence, noncoding_sequence, nrow = 1)
      # legend ('topleft',rownames(new_sum),
      #         cex=1.0,
      #         fill=colors)
      # 
      
    }
  })
}


shinyApp(ui = ui, server = server)