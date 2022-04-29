source('variables_and_functions_and_libraries.R')

#set max size limit to 30MB
options(shiny.maxRequestSize = 30*1024^2)
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId='num',
                  label = 'Choose the number of bins for reference list',
                  value = 25, min = 1, max = 100),
      sliderInput(inputId='reference_list_bins',
                  label = 'Choose the number of genes per bin for reference list',
                  value = 25, min = 1, max = 1000),
      sliderInput(inputId='reference_list_separated_bins', 
                  label = 'Choose the number of genes per bin for separated reference list', 
                  value = 30, min = 1, max = 1000),
      fileInput(inputId = 'clustalFile', 
                label = 'Choose a .clustal file',
                accept = '.clustal'), #can be used to accept two files instead of two separate fileInputs
      fileInput(inputId = 'csvFile', 
                label = 'Choose a .csv file',
                accept = '.csv'),
      selectInput(inputId = 'metric', 
                  label = 'Select the metric',
                  choices = metrics),
      actionButton(inputId ='processButton', label ='Process'),
      textInput(inputId='title', label = 'Title'),
      radioButtons(inputId='plot_type', label = 'Select the type of plot',
                   choices=c('Histogram', 'Bar chart')), #violin, boxplot, dotplot, stem and leaf?
      #scatter plot?
      selectInput(inputId= 'reference_list_type',
                  label = 'Select type of reference list',
                  choices = c('Transition', 'Transversion',
                              'Point mutations', 'Gaps',
                              'Mutations')),
      radioButtons(inputId='reference_plot', label = 'Select the type of reference list plot',
                   choices=c('Normal', 'Grouped')),
      "Click here to download summary : ",
      downloadButton("downloadData", "Download"),
    ),
  # Show a plot and dataframe
    mainPanel(
      plotOutput('plot', width='1200px', height='500px'),
      tableOutput('summary_dataframe'),
      tableOutput('df_clustal'),
      plotOutput('reference_list_plot', width='1200px', height='500px'),
      plotOutput('reference_list_separated_plot', width='1200px', height='500px'),
      
  )
  )
)
server <- function (input, output) {  #Load .clustal file and store it into start_stop reactive expression

  variables <- reactiveValues(clustal_file = NULL, csv_file = NULL, 
                              df_clustal_with_asterisk = NULL, 
                              indices_with_mutations = '', 
                              indices_without_mutations = '',
                              reference_list = NULL, #for the entire list (mutations)
                              reference_list_separated = NULL, #for separated NONCDS and CDS (mutations)
                              chosenMetric = 'Similarity')
  reference_lists <- reactiveValues(transitions = NULL,
                                    transitions_separated = NULL,
                                    transversions = NULL,
                                    transversions_separated = NULL,
                                    point_mutations = NULL,
                                    point_mutations_separated = NULL,
                                    gaps = NULL,
                                    gaps_separated = NULL,
                                    mutations = NULL,
                                    mutations_separated = NULL)
  
  ##### INPUTS #####
  # clustal file, csv file, chosen metric
  observe({
    variables$clustal_file <- input$clustalFile$datapath
  })
  observe({
    variables$csv_file <- input$csvFile$datapath
  })
  observe({
    variables$chosenMetric <- input$metric
  })
  
  #####
  
  df_clustal <- reactive({
    validate(
      need(variables$clustal_file != "", "Please select a data set123")
    )
    
    #Load .clustal file
    df_clustal <- read.delim(variables$clustal_file)
    
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
    variables$indices_with_mutations <-
      which(strsplit(df_clustal$sequence[df_clustal$strand == 'asterisk'], '')[[1]] == ' ')
    #indices with no mutations
    variables$indices_without_mutations <-
      which(strsplit(df_clustal$sequence[df_clustal$strand == 'asterisk'], '')[[1]] == '*')
    
    # extracting asterisk row
    df_clustal <- subset(df_clustal,strand != 'asterisk' )
    
    # Get numbers of strands; #e.g. seq1997 => 1997
    mutation_numbers <- as.numeric(unlist(lapply(df_clustal$strand, get_ordered)))
    # Orders strands in ascending order
    df_clustal <- df_clustal[order(mutation_numbers),]
    
    
    # Creates 5 reference lists 
    # TODO
    reference_lists$transitions <- create_reference_list_transitions(df_clustal, variables$indices_with_mutations)
    reference_lists$transversions <- create_reference_list_transversions(df_clustal, variables$indices_with_mutations)
    reference_lists$point_mutations <- reference_lists$transitions + reference_lists$transversions
    reference_lists$gaps <- create_reference_list_gaps(df_clustal, variables$indices_with_mutations)
    reference_lists$mutations <- reference_lists$gaps + reference_lists$point_mutations
    
    variables$reference_list <- create_reference_list(df_clustal, 
                                                      variables$indices_with_mutations)
    df_clustal
    
  })
  
  
  summary1 <- eventReactive(input$processButton, {
    df_clustal <- df_clustal()
    # Get names and number of strands
    strand_names <- df_clustal$strand
    number_of_strands <- nrow(df_clustal)

    length_of_sequence <- length(variables$indices_with_mutations) +
      length(variables$indices_without_mutations)
    MSA_metrics <- initialize_dataframes(number_of_strands, strand_names, 9)
    MSA_metrics <- fill_metrics (MSA_metrics, df_clustal, variables$indices_with_mutations, seq(1, length_of_sequence))

    ################################################################################
    ############################ CODING SEQUENCE ###################################
    ################################################################################

    #Load .csv file
    start_stop <- read.csv(variables$csv_file)

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

    # Creating reference lists that separate based on coding and non-coding regions
    reference_lists$transitions_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$transitions$reference_list)
    reference_lists$transversions_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$transversions$reference_list)
    reference_lists$point_mutations_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$point_mutations$reference_list)
    reference_lists$gaps_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$gaps$reference_list)
    reference_lists$mutations_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$mutations$reference_list)
    
    CDS_metrics <- initialize_dataframes(number_of_strands, strand_names, 9)
    CDS_metrics <- fill_metrics (MSA_metrics, df_clustal, variables$indices_with_mutations,
                                 coding_sequence_index)

    ################################################################################
    ############################ NON CODING SEQUECE ################################
    ################################################################################

    nonCDS_metrics <- initialize_dataframes(number_of_strands, strand_names, 9)
    nonCDS_metrics <- fill_metrics (MSA_metrics, df_clustal,
                                    variables$indices_with_mutations, noncoding_sequence_index)
  
    #Creating summary dataframe
    summary <- data.frame(cbind(MSA_metrics[[1]][,1], MSA_metrics[[2]][,1],
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
                     nonCDS_metrics[[9]][,1]))
    summary_column_names <- metrics
    colnames(summary) <- summary_column_names

    ## Create and fill the summary_pairwise dataframe
    summary_pairwise <- fill_summary_pairwise(MSA_metrics, CDS_metrics, nonCDS_metrics, number_of_strands, strand_names)
    colnames(summary_pairwise) <- summary_column_names
    summary <- cbind(summary, rownames(summary))
    colnames(summary)[colnames(summary) == 'rownames(summary)'] <- 'strands'
    summary
  })
  

  # OUTPUTS
  #####
  output$df_clustal <- renderTable({
    df_clustal()
  })
  output$summary_dataframe <- renderTable({
    df <- summary1()
  })
  output$plot <- renderPlot({
    metric <- variables$chosenMetric
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
    }
  })

  output$reference_list_plot <- renderPlot({
    validate(
      need(reference_lists$transitions != "" | reference_lists$transversions != "" | 
           reference_lists$point_mutations != "" | reference_lists$gaps != "" | 
           reference_lists$mutations != "", "Please select a data set"),
    )
    if(input$reference_list_type == 'Transition'){
      reference_list_histogram <- create_reference_list_histogram(reference_lists$transitions, 
                                                                  input$reference_list_bins, 
                                                                  isolate(input$title))
    }
    else if(input$reference_list_type == 'Transversion'){
      reference_list_histogram <- create_reference_list_histogram(reference_lists$transversions, 
                                                                  input$reference_list_bins, 
                                                                  isolate(input$title))
    }
    else if(input$reference_list_type == 'Point mutations'){
      reference_list_histogram <- create_reference_list_histogram(reference_lists$point_mutations, 
                                                                  input$reference_list_bins, 
                                                                  isolate(input$title))
    }
    else if(input$reference_list_type == 'Gaps'){
      reference_list_histogram <- create_reference_list_histogram(reference_lists$gaps, 
                                                                  input$reference_list_bins, 
                                                                  isolate(input$title))
    }
    else if(input$reference_list_type == 'Mutations'){
      reference_list_histogram <- create_reference_list_histogram(reference_lists$mutations, 
                                                                  input$reference_list_bins, 
                                                                  isolate(input$title))
    }
    plot(reference_list_histogram)
  })
  
  output$reference_list_separated_plot <- renderPlot({
    validate(
      need(reference_lists$transitions_separated != "" | reference_lists$transversions_separated != "" | 
             reference_lists$point_mutations_separated != "" | reference_lists$gaps_separated != "" | 
             reference_lists$mutations_separated != "", "Please select a data set"),
    )
    if(input$reference_list_type == 'Transition'){
      reference_list_separated_histogram <- create_reference_list_separated_histogram(reference_lists$transitions_separated, 
                                                                  isolate(input$title))
    }
    else if(input$reference_list_type == 'Transversion'){
      reference_list_separated_histogram <- create_reference_list_separated_histogram(reference_lists$transversions_separated, 
                                                                  isolate(input$title))
    }
    else if(input$reference_list_type == 'Point mutations'){
      reference_list_separated_histogram <- create_reference_list_separated_histogram(reference_lists$point_mutations_separated, 
                                                                  isolate(input$title))
    }
    else if(input$reference_list_type == 'Gaps'){
      reference_list_separated_histogram <- create_reference_list_separated_histogram(reference_lists$gaps_separated, 
                                                                  isolate(input$title))
    }
    else if(input$reference_list_type == 'Mutations'){
      reference_list_separated_histogram <- create_reference_list_separated_histogram(reference_lists$mutations_separated, 
                                                                  isolate(input$title))
    }
    plot(reference_list_separated_histogram)
  })
  

  output$reference_list_grouped_plot <- renderPlot({
    validate(
      need(variables$reference_list != "", "Please select a data set - dotplot")
    )
    reference_list_grouped <- create_reference_list_grouped(variables$reference_list, input$reference_groups)
    ggplot(reference_list_grouped,
           aes(x = group,
               y = number_of_mutations,
               fill = 'red'))+
      geom_bar(stat='identity', position ='dodge')+
      theme_dark()
  })

  # output$reference_list_separated_plot <- renderPlot({
  #   validate(
  #     need(variables$reference_list_separated != "", "Needs to be processed first")
  #   )
  #   ggplot(variables$reference_list_separated,
  #          aes(x = reorder(Type, Count),
  #              y = Number_of_mutations,
  #              fill = 'red'))+
  #     geom_bar(stat='identity', position ='dodge')+
  #     theme_dark()
  # })
}


  

shinyApp(ui = ui, server = server)