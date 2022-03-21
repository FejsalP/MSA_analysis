data <- data.frame(values = c(4, 1, 3, 6, 7, 3),  # Create example data
                   group = rep(c("group 1",
                                 "group 2",
                                 "group 3"),
                               each = 2),
                   subgroup = LETTERS[1:2])
data
data_base <- reshape(data,                        # Modify data for Base R barplot
                     idvar = "subgroup",
                     timevar = "group",
                     direction = "wide")
row.names(data_base) <- data_base$subgroup
data_base <- data_base[ , 2:ncol(data_base)]
colnames(data_base) <- c("group 1", "group 2", "group 3")
data_base <- as.matrix(data_base)
data_base

ggplot(data,                                      # Grouped barplot using ggplot2
       aes(x = group,
           y = values,
           fill = subgroup)) +
  geom_bar(stat = "identity",
           position = "dodge")

new_df <- data.frame(strands = c(1,2,4,1,1,1,4,3))
new_df <- cbind(new_df, rownames(new_df))
colnames(new_df)[colnames(new_df) == "rownames(new_df)"] <- "strands1"

ggplot(new_df,
       aes(x = strands, y = strands1)) +
  geom_bar(stat = "identity",
           position = "dodge")
