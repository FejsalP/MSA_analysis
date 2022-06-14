source("variables_and_functions_and_libraries.R")

# set max size limit to 30MB
options(shiny.maxRequestSize = 30 * 1024^2)
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      # FILE INPUTS
      fileInput(
        inputId = "clustalFile",
        label = "Choose a .clustal file",
        accept = ".clustal"
      ),
      fileInput(
        inputId = "csvFile",
        label = "Choose a .csv file",
        accept = ".csv"
      ),
      actionButton(inputId = "processButton", label = "Process"),
      
      textOutput('test'),
      
      ## GRAPHS

      # GRAPH 1 OPTIONS - METRIC GRAPH
      tags$h4("Metric graph"),
      selectInput(
        inputId = "metric",
        label = "Select the metric",
        choices = metrics
      ),
      # checkboxGroupInput(
      #   "chosen_graphs", "Graphs to show:", # CHOOSE how many plots, 1, 2, or 3
      #   c(
      #     "Entire" = "entire",
      #     "Coding" = "CDS",
      #     "Non-Coding" = "nonCDS"
      #   )
      # ),
      # # # names and sizes
      textInput(inputId = "metric_graph_title", label = "Title"),
      textInput(inputId = "metric_graph_x_axis", label = "X axis name"),
      textInput(inputId = "metric_graph_y_axis", label = "Y axis name"),
      sliderInput(
        inputId = "metric_graph_x_axis_text_size",
        label = "Choose the size of text on X axis",
        value = 12,
        min = 1,
        max = 20
      ),
      sliderInput(
        inputId = "metric_graph_title_text_size",
        label = "Choose the size of text on title",
        value = 12,
        min = 1,
        max = 35
      ),

      # # GRAPH 2 OPTIONS - REFERENCE LIST
      tags$h4("Reference list graph"),
      selectInput(
        inputId = "reference_list_type",
        label = "Select type of reference list",
        choices = c(
          "Transition",
          "Transversion",
          "Point mutations",
          "Gaps",
          "Mutations"
        )
      ),
      sliderInput(
        inputId = "reference_list_bins",
        label = "Choose the number of nucleotides per bin for reference list",
        value = 25,
        min = 1,
        max = 1000
      ),

      # names and sizes
      textInput(inputId = "reference_list_graph_title", label = "Title"),
      textInput(inputId = "reference_list_graph_x_axis", label = "X axis name"),
      sliderInput(
        inputId = "reference_list_graph_x_axis_text_size",
        label = "Choose the size of text on X axis",
        value = 12,
        min = 1,
        max = 20
      ),
      sliderInput(
        inputId = "reference_list_graph_y_axis_text_size",
        label = "Choose the size of text on Y axis",
        value = 12,
        min = 1,
        max = 20
      ),
      sliderInput(
        inputId = "reference_list_graph_title_text_size",
        label = "Choose the size of text on title",
        value = 12,
        min = 1,
        max = 35
      ),

      # GRAPH 3 OPTIONS - REFERENCE LIST GROUPED
      tags$h4("Reference list grouped graph"),
      selectInput(
        inputId = "reference_list_grouped_type",
        label = "Select type of reference list",
        choices = c(
          "Transition",
          "Transversion",
          "Point mutations",
          "Gaps",
          "Mutations"
        )
      ),
      textInput(inputId = "reference_list_grouped_graph_title", label = "Title"),
      textInput(inputId = "reference_list_grouped_graph_x_axis", label = "X axis name"),
      sliderInput(
        inputId = "reference_list_grouped_graph_x_axis_text_size",
        label = "Choose the size of text on X axis",
        value = 12,
        min = 1,
        max = 20
      ),
      sliderInput(
        inputId = "reference_list_grouped_graph_y_axis_text_size",
        label = "Choose the size of text on Y axis",
        value = 12,
        min = 1,
        max = 20
      ),
      sliderInput(
        inputId = "reference_list_grouped_graph_title_text_size",
        label = "Choose the size of text on title",
        value = 12,
        min = 1,
        max = 35
      ),
      sliderInput(
        inputId = "reference_list_grouped_graph_x_axis_text_ticks_size",
        label = "Choose the size of text on X axis ticks",
        value = 12,
        min = 1,
        max = 20
      ),
      # GRAPH 4 OPTIONS - REFERENCE LIST GROUPED - TYPES OF MUTATIONS
      tags$h4("Reference list grouped (types of mutations) graph"),
      textInput(inputId = "reference_list_grouped_mutations_type_graph_title", label = "Title"),
      textInput(inputId = "reference_list_grouped_mutations_type_graph_x_axis", label = "X axis name"),
      sliderInput(
        inputId = "reference_list_grouped_mutations_type_graph_x_axis_text_size",
        label = "Choose the size of text on X axis",
        value = 12,
        min = 1,
        max = 20
      ),
      sliderInput(
        inputId = "reference_list_grouped_mutations_type_graph_y_axis_text_size",
        label = "Choose the size of text on Y axis",
        value = 12,
        min = 1,
        max = 20
      ),
      sliderInput(
        inputId = "reference_list_grouped_mutations_type_graph_title_text_size",
        label = "Choose the size of text on title",
        value = 12,
        min = 1,
        max = 35
      ),
      sliderInput(
        inputId = "reference_list_grouped_mutations_type_graph_x_axis_text_ticks_size",
        label = "Choose the size of text on X axis",
        value = 12,
        min = 1,
        max = 20
      ),


      # # DOWNLOAD BUTTONS
      "Click here to download summary table: ",
      downloadButton("download_summary", "Summary"),
      br(),
      "Click here to download sequences: ",
      downloadButton("download_sequences", "Sequences"),
      width = 2,
    ),
    # Show a plot and dataframe
    mainPanel(
      tags$h1("Plots"),
      fluidRow(plotOutput("plot")),
      br(),
      fluidRow(plotOutput("reference_list_plot")),
      br(),
      fluidRow(
        column(
          6, plotOutput("reference_list_separated_plot")
        ),
        column(6, plotOutput(
          "mutations_type_plot"
        ))
      ),
      width = 10
      # tableOutput('summary_dataframe'),
      # tableOutput('df_clustal'),
    )
  )
)
server <- function(input, output) { # Load .clustal file and store it into start_stop reactive expression


  variables <- reactiveValues(
    clustal_file = NULL,
    csv_file = NULL,
    df_clustal_with_asterisk = NULL,
    indices_with_mutations = "",
    indices_without_mutations = "",
    reference_list = NULL,
    # for the entire list (mutations)
    reference_list_separated = NULL,
    # for separated NONCDS and CDS (mutations)
    chosenMetric = "Similarity",
    summary_pairwise = NULL,
    plotCheck = "",
    progress = "Start"
  )
  reference_lists <- reactiveValues(
    transitions = NULL,
    transitions_separated = NULL,
    transversions = NULL,
    transversions_separated = NULL,
    point_mutations = NULL,
    point_mutations_separated = NULL,
    gaps = NULL,
    gaps_separated = NULL,
    mutations = NULL,
    mutations_separated = NULL
  )

  output$test <- renderText({
    variables$progress
  })
  
  ##### INPUTS #####
  # clustal file, csv file, chosen metric
  observe({
    # "C:\\Users\\fejsa\\OneDrive\\Desktop\\Graduation project\\COMSAA\\input\\MERS MSA.clustal"
    # input$clustalFile$datapath
    variables$clustal_file <-
      input$clustalFile$datapath
  })
  observe({
    # "C:\\Users\\fejsa\\OneDrive\\Desktop\\Graduation project\\COMSAA\\input\\MERS MSA.csv"
    # input$csvFile$datapath
    variables$csv_file <-
      input$csvFile$datapath
  })
  observe({
    variables$chosenMetric <- input$metric
  })

  #####


  df_clustal <- reactive({
    validate(need(variables$clustal_file != "", "Please select a data set"))

    # Load .clustal file
    df_clustal <- read.delim(variables$clustal_file)

    # Splitting rows on strands and sequences
    df_clustal[c("strand", "sequence")] <-
      str_split_fixed(
        df_clustal$CLUSTAL.O.1.2.4..multiple.sequence.alignment,
        " ",
        2
      )
    # Changing the name of the rows with only asterisks
    df_clustal$strand[df_clustal$strand == ""] <- "asterisk"


    # drop first column that has everything combined
    df_clustal <- subset(df_clustal, select = -c(1))
    # clean leading spaces
    df_clustal$sequence <- lapply(df_clustal$sequence, trim_leading)

    # calculate number of characters for each line,
    # all except the last one should be of length 60
    first_strand <- df_clustal[1, 1]
    length_of_lines <-
      lapply(df_clustal[df_clustal$strand == first_strand, 2], nchar)
    length_of_sequence <- nchar(df_clustal$sequence[1])

    asterisk_only <-
      df_clustal$sequence[df_clustal$strand == "asterisk"]


    # Fixing the length of rows with asterisks
    # All rows have the same length as the rows with sequences
    df_clustal$sequence[df_clustal$strand == "asterisk"] <-
      mapply(add_leading_spaces, asterisk_only, length_of_lines)

    # groups sequences based on strand
    # strands are now in one line
    df_clustal <- df_clustal %>%
      group_by(strand) %>%
      summarise(sequence = paste(sequence, collapse = ""))

    # indices with mutations

    variables$indices_with_mutations <-
      which(strsplit(df_clustal$sequence[df_clustal$strand == "asterisk"], "")[[1]] == " ")
    # indices with no mutations
    variables$indices_without_mutations <-
      which(strsplit(df_clustal$sequence[df_clustal$strand == "asterisk"], "")[[1]] == "*")

    # extracting asterisk row
    df_clustal <- subset(df_clustal, strand != "asterisk")


    # Get numbers of strands; #e.g. seq1997 => 1997
    mutation_numbers <-
      as.numeric(unlist(lapply(df_clustal$strand, get_ordered)))
    # Orders strands in ascending order
    df_clustal <- df_clustal[order(mutation_numbers), ]


    # Creates 5 reference lists

    reference_lists$transitions <- create_reference_list_transitions(df_clustal, variables$indices_with_mutations)
    reference_lists$transversions <- create_reference_list_transversions(df_clustal, variables$indices_with_mutations)
    reference_lists$point_mutations <- reference_lists$transitions + reference_lists$transversions
    reference_lists$gaps <- create_reference_list_gaps(df_clustal, variables$indices_with_mutations)
    reference_lists$mutations <- reference_lists$gaps + reference_lists$point_mutations

    variables$reference_list <- create_reference_list(
      df_clustal,
      variables$indices_with_mutations
    )
    variables$progress <- "Created DF"
    df_clustal
  })


  summary1 <- eventReactive(input$processButton, {
    withProgress(message = 'Processing', value = 0, {
    incProgress(1/5, detail = paste("Finished creating df"))
    df_clustal <- df_clustal()
    # Get names and number of strands
    strand_names <- df_clustal$strand
    number_of_strands <- nrow(df_clustal)

    length_of_sequence <- length(variables$indices_with_mutations) +
      length(variables$indices_without_mutations)
    incProgress(1/5, detail = paste("Calculating entire sequence"))
    MSA_metrics <- initialize_dataframes(number_of_strands, strand_names, 9)
    MSA_metrics <- fill_metrics(MSA_metrics, df_clustal, variables$indices_with_mutations, seq(1, length_of_sequence))
    ################################################################################
    ############################ CODING SEQUENCE ###################################
    ################################################################################
    incProgress(1/5, detail = paste("Loading .csv file"))
    
    # Load .csv file
    start_stop <- read.csv(variables$csv_file)

    start_stop <- subset(start_stop, select = c("Start", "Stop"))
    start_stop <- start_stop[order(start_stop$Start), ]
    # if there is only one row then mapply unlists the sequence
    if (nrow(start_stop) > 1) {
      start_stop$Sequences <- mapply(create_sequence, start_stop$Start, start_stop$Stop)
    } else {
      start_stop$Sequences <- list(mapply(create_sequence, start_stop$Start, start_stop$Stop))
    }
    coding_sequence_index <- c()

    for (i in seq(1, nrow(start_stop))) {
      coding_sequence_index <- append(coding_sequence_index, unlist(start_stop[i, ]$Sequences))
    }
    coding_sequence_index <- unlist(coding_sequence_index)
    coding_sequence_index <- unique(coding_sequence_index)

    # change to how coding and noncoding indices are calculated
    # ignore gaps in the reference sequence, find new indices
    # CDI - coding sequence index
    reference_sequence <- as.character(df_clustal[1, 2])
    coding_sequence_index <- new_CDI(coding_sequence_index, reference_sequence)


    noncoding_sequence_index <- setdiff(seq(1, nchar(df_clustal[1, 2])), coding_sequence_index)


    # Creating reference lists that separate based on coding and non-coding regions
    reference_lists$transitions_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$transitions$reference_list)
    reference_lists$transversions_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$transversions$reference_list)
    reference_lists$point_mutations_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$point_mutations$reference_list)
    reference_lists$gaps_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$gaps$reference_list)
    reference_lists$mutations_separated <- create_reference_list_separated(coding_sequence_index, noncoding_sequence_index, reference_lists$mutations$reference_list)
    incProgress(1/5, detail = paste("Calculating CDS"))
    CDS_metrics <- initialize_dataframes(number_of_strands, strand_names, 9)
    CDS_metrics <- fill_metrics(
      MSA_metrics, df_clustal, variables$indices_with_mutations,
      coding_sequence_index
    )
    ################################################################################
    ############################ NON CODING SEQUECE ################################
    ################################################################################
    incProgress(1/5, detail = paste("Calculating nonCDS"))
    nonCDS_metrics <- initialize_dataframes(number_of_strands, strand_names, 9)
    nonCDS_metrics <- fill_metrics(
      MSA_metrics, df_clustal,
      variables$indices_with_mutations, noncoding_sequence_index
    )
    # Creating summary dataframe
    summary <- data.frame(cbind(
      MSA_metrics[[1]][, 1], MSA_metrics[[2]][, 1],
      MSA_metrics[[3]][, 1], MSA_metrics[[4]][, 1],
      MSA_metrics[[5]][, 1], MSA_metrics[[6]][, 1],
      MSA_metrics[[7]][, 1], MSA_metrics[[8]][, 1],
      MSA_metrics[[9]][, 1],
      CDS_metrics[[1]][, 1], CDS_metrics[[2]][, 1],
      CDS_metrics[[3]][, 1], CDS_metrics[[4]][, 1],
      CDS_metrics[[5]][, 1], CDS_metrics[[6]][, 1],
      CDS_metrics[[7]][, 1], CDS_metrics[[8]][, 1],
      CDS_metrics[[9]][, 1],
      nonCDS_metrics[[1]][, 1], nonCDS_metrics[[2]][, 1],
      nonCDS_metrics[[3]][, 1], nonCDS_metrics[[4]][, 1],
      nonCDS_metrics[[5]][, 1], nonCDS_metrics[[6]][, 1],
      nonCDS_metrics[[7]][, 1], nonCDS_metrics[[8]][, 1],
      nonCDS_metrics[[9]][, 1]
    ))
    summary_column_names <- metrics
    colnames(summary) <- summary_column_names

    ## Create and fill the summary_pairwise dataframe
    variables$summary_pairwise <- fill_summary_pairwise(MSA_metrics, CDS_metrics, nonCDS_metrics, number_of_strands, strand_names)
    colnames(variables$summary_pairwise) <- summary_column_names


    summary <- cbind(summary, strand_names)
    colnames(summary)[colnames(summary) == "strand_names"] <- "Strands"

    summary <- cbind(summary, rownames(summary))
    colnames(summary)[colnames(summary) == "rownames(summary)"] <- "Count"
    variables$plotCheck <- "Pass"
    variables$progress <- "Finished"
    })
    
    summary
  })

  # OUTPUTS
  #####
  output$plot <- renderPlot({
    #     validate(need(variables$plotCheck != "", "Please select a data set"),
    # )
    metric <- variables$chosenMetric
    sum <- summary1()

    selected_column <- which(colnames(sum) == metric)
    selected_columns <- get_selected_columns(selected_column)
    selected_names_of_columns <- colnames(sum)[selected_columns]
    entire_sequence <- create_bar_plot(
      sum,
      selected_columns[1],
      selected_names_of_columns[1],
      input$metric_graph_title,
      input$metric_graph_x_axis,
      input$metric_graph_y_axis,
      input$metric_graph_title_text_size,
      input$metric_graph_x_axis_text_size,
      ""
    )
    coding_sequence <- create_bar_plot(
      sum,
      selected_columns[2],
      selected_names_of_columns[2],
      input$metric_graph_title,
      input$metric_graph_x_axis,
      input$metric_graph_y_axis,
      input$metric_graph_title_text_size,
      input$metric_graph_x_axis_text_size,
      "CDS"
    )
    noncoding_sequence <- create_bar_plot(
      sum,
      selected_columns[3],
      selected_names_of_columns[3],
      input$metric_graph_title,
      input$metric_graph_x_axis,
      input$metric_graph_y_axis,
      input$metric_graph_title_text_size,
      input$metric_graph_x_axis_text_size,
      "nonCDS"
    )

    ggarrange(entire_sequence, coding_sequence, noncoding_sequence,
      nrow =
        1
    )
  })

  output$reference_list_plot <- renderPlot({
    validate(
      need(reference_lists$transitions != "" | reference_lists$transversions != "" |
        reference_lists$point_mutations != "" | reference_lists$gaps != "" |
        reference_lists$mutations != "", "Please select a data set"),
    )
    if (input$reference_list_type == "Transition") {
      reference_list_histogram <- create_reference_list_histogram(
        reference_lists$transitions,
        input$reference_list_bins,
        input$reference_list_graph_title,
        input$reference_list_graph_x_axis,
        input$reference_list_graph_x_axis_text_size,
        input$reference_list_graph_y_axis_text_size,
        input$reference_list_graph_title_text_size,
        "transitions"
      )
    } else if (input$reference_list_type == "Transversion") {
      reference_list_histogram <- create_reference_list_histogram(
        reference_lists$transversions,
        input$reference_list_bins,
        input$reference_list_graph_title,
        input$reference_list_graph_x_axis,
        input$reference_list_graph_x_axis_text_size,
        input$reference_list_graph_y_axis_text_size,
        input$reference_list_graph_title_text_size,
        "transversions"
      )
    } else if (input$reference_list_type == "Point mutations") {
      reference_list_histogram <- create_reference_list_histogram(
        reference_lists$point_mutations,
        input$reference_list_bins,
        input$reference_list_graph_title,
        input$reference_list_graph_x_axis,
        input$reference_list_graph_x_axis_text_size,
        input$reference_list_graph_y_axis_text_size,
        input$reference_list_graph_title_text_size,
        "point mutations"
      )
    } else if (input$reference_list_type == "Gaps") {
      reference_list_histogram <- create_reference_list_histogram(
        reference_lists$gaps,
        input$reference_list_bins,
        input$reference_list_graph_title,
        input$reference_list_graph_x_axis,
        input$reference_list_graph_x_axis_text_size,
        input$reference_list_graph_y_axis_text_size,
        input$reference_list_graph_title_text_size,
        "gaps"
      )
    } else if (input$reference_list_type == "Mutations") {
      reference_list_histogram <- create_reference_list_histogram(
        reference_lists$mutations,
        input$reference_list_bins,
        input$reference_list_graph_title,
        input$reference_list_graph_x_axis,
        input$reference_list_graph_x_axis_text_size,
        input$reference_list_graph_y_axis_text_size,
        input$reference_list_graph_title_text_size,
        "mutations"
      )
    }
    plot(reference_list_histogram)
  })

  output$reference_list_separated_plot <- renderPlot({
    validate(
      need(reference_lists$transitions_separated != "" | reference_lists$transversions_separated != "" |
        reference_lists$point_mutations_separated != "" | reference_lists$gaps_separated != "" |
        reference_lists$mutations_separated != "", "Please select a data set"),
    )
    if (input$reference_list_grouped_type == "Transition") {
      reference_list_separated_histogram <- create_reference_list_separated_histogram(
        reference_lists$transitions_separated,
        input$reference_list_grouped_graph_title,
        input$reference_list_grouped_graph_x_axis,
        input$reference_list_grouped_graph_x_axis_text_size,
        input$reference_list_grouped_graph_y_axis_text_size,
        input$reference_list_grouped_graph_title_text_size,
        input$reference_list_grouped_graph_x_axis_text_ticks_size,
        "transitions"
      )
    } else if (input$reference_list_grouped_type == "Transversion") {
      reference_list_separated_histogram <- create_reference_list_separated_histogram(
        reference_lists$transversions_separated,
        input$reference_list_grouped_graph_title,
        input$reference_list_grouped_graph_x_axis,
        input$reference_list_grouped_graph_x_axis_text_size,
        input$reference_list_grouped_graph_y_axis_text_size,
        input$reference_list_grouped_graph_title_text_size,
        input$reference_list_grouped_graph_x_axis_text_ticks_size,
        "transversions"
      )
    } else if (input$reference_list_grouped_type == "Point mutations") {
      reference_list_separated_histogram <- create_reference_list_separated_histogram(
        reference_lists$point_mutations_separated,
        input$reference_list_grouped_graph_title,
        input$reference_list_grouped_graph_x_axis,
        input$reference_list_grouped_graph_x_axis_text_size,
        input$reference_list_grouped_graph_y_axis_text_size,
        input$reference_list_grouped_graph_title_text_size,
        input$reference_list_grouped_graph_x_axis_text_ticks_size,
        "point mutations"
      )
    } else if (input$reference_list_grouped_type == "Gaps") {
      reference_list_separated_histogram <- create_reference_list_separated_histogram(
        reference_lists$gaps_separated,
        input$reference_list_grouped_graph_title,
        input$reference_list_grouped_graph_x_axis,
        input$reference_list_grouped_graph_x_axis_text_size,
        input$reference_list_grouped_graph_y_axis_text_size,
        input$reference_list_grouped_graph_title_text_size,
        input$reference_list_grouped_graph_x_axis_text_ticks_size,
        "gaps"
      )
    } else if (input$reference_list_grouped_type == "Mutations") {
      reference_list_separated_histogram <- create_reference_list_separated_histogram(
        reference_lists$mutations_separated,
        input$reference_list_grouped_graph_title,
        input$reference_list_grouped_graph_x_axis,
        input$reference_list_grouped_graph_x_axis_text_size,
        input$reference_list_grouped_graph_y_axis_text_size,
        input$reference_list_grouped_graph_title_text_size,
        input$reference_list_grouped_graph_x_axis_text_ticks_size,
        "mutations"
      )
    }
    plot(reference_list_separated_histogram)
  })

  output$mutations_type_plot <- renderPlot({
    validate(
      need(reference_lists$transitions_separated != "" | reference_lists$transversions_separated != "" |
        reference_lists$gaps_separated != "", "Please select a data set"),
    )
    mutations_types <- create_mutations_type_bar_plot(
      reference_lists$transitions_separated,
      reference_lists$transversions_separated,
      reference_lists$gaps_separated,
      input$reference_list_grouped_mutations_type_graph_title,
      input$reference_list_grouped_mutations_type_graph_x_axis,
      input$reference_list_grouped_mutations_type_graph_x_axis_text_size,
      input$reference_list_grouped_mutations_type_graph_x_axis_text_size,
      input$reference_list_grouped_mutations_type_graph_title_text_size,
      input$reference_list_grouped_mutations_type_graph_x_axis_text_ticks_size
    )
    plot(mutations_types)
  })
  output$download_summary <- downloadHandler(
    filename = function() {
      paste(".csv", sep = "")
    },
    content = function(file) {
      write.csv(summary1(), file)
    }
  )
  output$download_sequences <- downloadHandler(
    filename = function() {
      paste(".csv", sep = "")
    },
    content = function(file) {
      write.csv(df_clustal(), file)
    }
  )
}
shinyApp(ui = ui, server = server)
