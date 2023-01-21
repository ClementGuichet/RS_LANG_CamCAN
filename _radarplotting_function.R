library(fmsb)
library(RColorBrewer)
# Helper functions for radar plot

# Dataframe must be organized as follow:

# rows = facet grid different radar plots
# cols = labels within each radar plot
# This function can handle multiple radar plotting with a different set of labels for each radar plot
#
# Currently allowing 12 colors at once

radarplotting <- function(dataframe, max, min, mfrow_x, mfrow_y, inverse_grad = FALSE, palette, alpha, label_size, round = TRUE) {
  # Set the max & min of each column
  if (inverse_grad == FALSE) {
    max_min <- as.data.frame(t(data.frame(max = rep(max, times = ncol(dataframe)), min = rep(min, times = ncol(dataframe)))))
  } else {
    max_min <- as.data.frame(t(data.frame(max = rep(min, times = ncol(dataframe)), min = rep(max, times = ncol(dataframe)))))
  }
  rownames(max_min) <- c("Max", "Min")
  names(max_min) <- names(dataframe)

  Radar_plot_1 <- rbind(max_min, dataframe)

  # Reduce plot margin using par()
  # Split the screen in 3 parts
  op <- par(mar = c(1, 2, 2, 2))
  par(mfrow = c(mfrow_x, mfrow_y))

  for (i in 1:nrow(dataframe)) {
    # When plotting a specific row, remove any columns that takes NA
    Radar_plot_subset <- Radar_plot_1[c(1, 2, i + 2), ] %>%
      select_if(~ !any(is.na(.)))
    # Define colors and titles
    colors <- palette
    if (round == TRUE) {
      display_max <- round(max(Radar_plot_1))
      display_min <- round(min(Radar_plot_1))
      display_min_2 <- round(min(Radar_plot_1)) + (display_max - display_min) / 4
      display_min_3 <- round(min(Radar_plot_1)) + 2 * (display_max - display_min) / 4
      display_min_4 <- round(min(Radar_plot_1)) + 3 * (display_max - display_min) / 4
    } else {
      display_max <- max(Radar_plot_1)
      display_min <- min(Radar_plot_1)
      display_min_2 <- min(Radar_plot_1) + (display_max - display_min) / 4
      display_min_3 <- min(Radar_plot_1) + 2 * (display_max - display_min) / 4
      display_min_4 <- min(Radar_plot_1) + 3 * (display_max - display_min) / 4
    }

    if (inverse_grad == FALSE) {
      create_beautiful_radarchart(
        data = Radar_plot_subset, caxislabels = c(display_min, display_min_2, display_min_3, display_min_4, display_max),
        color = colors[i], title = row.names(Radar_plot_subset)[3], alpha = , label_size = label_size
      )
    } else {
      create_beautiful_radarchart(
        data = Radar_plot_subset, caxislabels = c(display_max, display_min_4, display_min_3, display_min_2, display_min),
        color = colors[i], title = row.names(Radar_plot_subset)[3], alpha = alpha, label_size = label_size
      )
    }
  }
}


radarplotting_overlap <- function(dataframe, max, min, mfrow_x, mfrow_y, title_fill, inverse_grad = FALSE, round = TRUE, palette, alpha, label_size) {
  # Set the max & min of each column
  if (inverse_grad == FALSE) {
    max_min <- as.data.frame(t(data.frame(max = rep(max, times = ncol(dataframe)), min = rep(min, times = ncol(dataframe)))))
  } else {
    max_min <- as.data.frame(t(data.frame(max = rep(min, times = ncol(dataframe)), min = rep(max, times = ncol(dataframe)))))
  }
  rownames(max_min) <- c("Max", "Min")
  names(max_min) <- names(dataframe)

  Radar_plot_1 <- rbind(max_min, dataframe) %>% mutate_all(., ~ replace(., is.na(.), 0))

  # Reduce plot margin using par()
  # Split the screen in 3 parts
  op <- par(mar = c(1, 2, 2, 2))
  par(mfrow = c(mfrow_x, mfrow_y))

  colors <- palette
  if (round == TRUE) {
    display_max <- round(max(Radar_plot_1))
    display_min <- round(min(Radar_plot_1))
    display_min_2 <- round(min(Radar_plot_1)) + (display_max - display_min) / 4
    display_min_3 <- round(min(Radar_plot_1)) + 2 * (display_max - display_min) / 4
    display_min_4 <- round(min(Radar_plot_1)) + 3 * (display_max - display_min) / 4
  } else {
    display_max <- max(Radar_plot_1)
    display_min <- min(Radar_plot_1)
    display_min_2 <- min(Radar_plot_1) + (display_max - display_min) / 4
    display_min_3 <- min(Radar_plot_1) + 2 * (display_max - display_min) / 4
    display_min_4 <- min(Radar_plot_1) + 3 * (display_max - display_min) / 4
  }


  if (inverse_grad == FALSE) {
    create_beautiful_radarchart(
      data = Radar_plot_1, caxislabels = c(display_min, display_min_2, display_min_3, display_min_4, display_max),
      color = colors, title = title_fill, alpha = alpha, label_size = label_size
    )
  } else {
    create_beautiful_radarchart(
      data = Radar_plot_1, caxislabels = c(display_max, display_min_4, display_min_3, display_min_2, display_min),
      color = colors, title = title_fill, alpha = alpha, label_size = label_size
    )
  }
}


create_beautiful_radarchart <- function(data, color = "#00AFBB",
                                        vlabels = colnames(data), label_size = label_size,
                                        caxislabels = NULL, title = NULL, alpha = 0.2, ...) {
  radarchart(
    data,
    axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, alpha), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 3, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black",
    # Variable labels
    vlcex = label_size, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}
