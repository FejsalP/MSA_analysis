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

output$reference_list_histogram <- renderPlot({
  validate(
    need(variables$reference_list != "", "Please select a data set - histogram")
  )
  hist(variables$reference_list$reference_list,
       breaks = input$num,
       main = isolate(input$title))
})

output$reference_list_plot <- renderPlot({
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

output$reference_list_separated_plot <- renderPlot({
  validate(
    need(variables$reference_list_separated != "", "Needs to be processed first")
  )
  ggplot(variables$reference_list_separated,
         aes(x = reorder(Type, Count),
             y = Number_of_mutations,
             fill = 'red'))+
    geom_bar(stat='identity', position ='dodge')+
    theme_dark()
})