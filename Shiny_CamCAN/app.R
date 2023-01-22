library(shiny)
library(tidyverse)
library(data.table)
library(ggpubr)
library(bslib)

light <- bs_theme(bootswatch = "minty")
dark <- bslib::bs_theme(
  bg = "#101010",
  fg = "#FDF7F7",
  primary = "#ED79F9",
  secondary = "#1e90ff"
)

ui <- fluidPage(
  theme = light,
  titlePanel(
    h1("Topological trajectory across lifespan"),
    windowTitle = "Topological Trajectory - Shiny App"
  ),
  p(""),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      p("Explore the most likely topological reconfiguration trajectory"),
      checkboxInput("dark_mode", "Dark mode", value = FALSE),
      p("Visit", em(a("github.com/ClementGuichet/RS_LANG_CamCAN", href = "https://github.com/ClementGuichet/RS_LANG_CamCAN")), "for more information"),
      checkboxGroupInput("RSN_choice",
        h4("Pick resting-state networks"),
        choices = list(
          "Auditory" = "Auditory",
          "CON" = "CON",
          "DAN" = "DAN",
          "DMN" = "DMN",
          "FPN" = "FPN",
          "Language" = "Language",
          "SMN" = "SMN",
          "PMM" = "PMM",
          "VMM" = "VMM",
          "Visual_1" = "Visual_1",
          "Visual_2" = "Visual_2"
        )
      ),
      checkboxInput("all", "Plot all regions"),
      selectInput("functional_role",
        h4("Explore modular or interareal functional role reconfiguration"),
        choices = list(
          "Modular" = "Modular",
          "Interareal" = "Interareal"
        )
      ),
      # Input: Decimal interval with step value ----
      sliderInput("decimal", "Z-scored probability threshold for each region's trajectories:",
        min = -3, max = 3,
        value = 1.5, step = 0.25
      ),
      actionButton("button", "Plot trajectory", width = 460, style = "font-size: 20px"),
    ),
    mainPanel(
      h3(textOutput("select_var"), align = "center"),
      shinycssloaders::withSpinner(
        plotOutput("plot", height = 650),
        hide.ui = FALSE
      )
    )
  )
)


# Define server logic ----
server <- function(input, output, session) {
  observe(session$setCurrentTheme(
    if (isTRUE(input$dark_mode)) dark else light
  ))

  observe({
    x <- input$all

    if (!(is.null(x))) {
      y <- character(1)

      updateCheckboxGroupInput(session, "RSN_choice",
        choices = list(
          "Auditory" = "Auditory",
          "CON" = "CON",
          "DAN" = "DAN",
          "DMN" = "DMN",
          "FPN" = "FPN",
          "Language" = "Language",
          "SMN" = "SMN",
          "PMM" = "PMM",
          "VMM" = "VMM",
          "Visual_1" = "Visual_1",
          "Visual_2" = "Visual_2"
        ),
        selected = y
      )
    }
  })

  observe({
    x <- input$RSN_choice

    if (!(is.null(x))) {
      y <- character(1)

      updateCheckboxInput(session, "all", "Plot all regions", value = y)
    }
  })

  plot_button <- eventReactive(input$button, {
    input$RSN_choice
  })

  plot_button_bis <- eventReactive(input$button, {
    input$functional_role
  })

  sliderValue <- eventReactive(input$button, {
    input$decimal
  })

  output$select_var <- renderText({
    if (length(plot_button()) > 1) {
      list_input <- list()
      for (i in 1:length(plot_button())) {
        list_input[[i]] <- plot_button()[i]
        tmp_bis <- rbindlist(lapply(list_input, as.data.table))
        RSN_text <- capture.output(cat(tmp_bis %>% as.matrix(), sep = "|"))
      }
      paste("This is the most probable trajectory for", RSN_text, "regions combined.")
    } else if (length(plot_button()) == 1) {
      list_input <- list()
      for (i in 1:length(plot_button())) {
        list_input[[i]] <- plot_button()[i]
        tmp_bis <- rbindlist(lapply(list_input, as.data.table))
        RSN_text <- capture.output(cat(tmp_bis %>% as.matrix(), sep = "|"))
      }
      paste("This is the most probable trajectory for", RSN_text, "regions.")
    } else {
      paste("This is the most probable trajectory for all regions combined.")
    }
  })

  output$plot <- renderPlot({
    modular <- read.csv("modular_CamCAN.csv") %>%
      dplyr::select(-X) %>%
      plyr::rename(c("X1st_network" = "1st_network"))
    interareal <- read.csv("interareal_CamCAN.csv") %>%
      dplyr::select(-X) %>%
      plyr::rename(c("X1st_network" = "1st_network"))

    ############################################################################
    ############################################################################
    ############################################################################

    trajectory <- function(composition, list_RSN, threshold) {
      if (colnames(composition)[4] == "Hub_consensus") {
        modular$Age_group <- factor(modular$Age_group, levels = c("Young", "Middle", "Old"))
        if (list_RSN != "All") {
          Cond_PMF <- modular %>%
            filter(grepl(list_RSN, `1st_network`)) %>%
            spread(Age_group, freq) %>%
            arrange(Hub_consensus) %>%
            group_by(Region) %>%
            group_split()
        } else {
          Cond_PMF <- modular %>%
            spread(Age_group, freq) %>%
            arrange(Hub_consensus) %>%
            group_by(Region) %>%
            group_split()
        }
      } else {
        interareal$Age_group <- factor(interareal$Age_group, levels = c("Young", "Middle", "Old"))

        if (list_RSN != "All") {
          Cond_PMF <- interareal %>%
            filter(grepl(list_RSN, `1st_network`)) %>%
            spread(Age_group, freq) %>%
            arrange(Bridgeness) %>%
            group_by(Region) %>%
            group_split()
        } else {
          Cond_PMF <- interareal %>%
            spread(Age_group, freq) %>%
            arrange(Bridgeness) %>%
            group_by(Region) %>%
            group_split()
        }
      }

      # Convert an adjacency dataframe to a 2-column dataframe
      adjacency_to_2col <- function(data) {
        crossdata <- lapply(rownames(data), function(x) sapply(colnames(data), function(y) list(x, y, data[x, y])))
        crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
        crossdatamat <- t(crossdatatmp)
        if (colnames(data)[5] == "YM") {
          colnames(crossdatamat) <- c("Young", "Middle", "Value")
        } else {
          colnames(crossdatamat) <- c("Middle", "Old", "Value")
        }
        crossdatadf <- as.data.frame(crossdatamat, stringsAsFactors = F)
        crossdatadf[, 3] <- as.numeric(crossdatadf[, 3] %>% na.omit())
        return(crossdatadf %>% na.omit())
      }

      # For each region, find the most probable trajectory
      Cond_PMF_list <- list()
      for (i in 1:length(Cond_PMF)) {
        tmp <<- rbindlist(Cond_PMF[i])
        # Beware of the order here, first young then old
        outer_young_to_middle <- outer(tmp[, 4] %>% as.matrix(), tmp[, 5] %>% as.matrix()) %>% as.data.frame()
        if (colnames(composition)[4] == "Hub_consensus") {
          rownames(outer_young_to_middle) <- unlist(tmp$Hub_consensus)
          colnames(outer_young_to_middle) <- unlist(tmp$Hub_consensus)
        } else {
          rownames(outer_young_to_middle) <- unlist(tmp$Bridgeness)
          colnames(outer_young_to_middle) <- unlist(tmp$Bridgeness)
        }

        outer_young_to_middle <- outer_young_to_middle %>% mutate(YM = "YM")
        crossdatadf <- adjacency_to_2col(outer_young_to_middle)
        outer_young_to_middle_2 <- crossdatadf

        tmp_young_to_middle <- cbind(tmp %>% slice(1:nrow(outer_young_to_middle_2)) %>%
          dplyr::select(`1st_network`, Region), outer_young_to_middle_2) %>%
          plyr::rename(c("Value" = "Value_YM"))

        outer_middle_to_old <- outer(tmp_young_to_middle[, 5] %>% as.matrix(), tmp[, 6] %>% as.matrix()) %>% as.data.frame()
        if (colnames(composition)[4] == "Hub_consensus") {
          # Create a vector to bypass duplicate rows not allowed
          bypass_duplicate_row <- tmp_young_to_middle %>%
            mutate(helper_vector = rep(seq(1:nrow(.)))) %>%
            unite(., Young, "Young", "helper_vector", remove = FALSE) %>%
            unite(., Middle, "Middle", "helper_vector", remove = FALSE) %>%
            dplyr::select(-helper_vector)

          rownames(outer_middle_to_old) <- unlist(bypass_duplicate_row$Middle)
          colnames(outer_middle_to_old) <- unlist(tmp$Hub_consensus)
        } else {
          # Create a vector to bypass duplicate rows not allowed
          bypass_duplicate_row <- tmp_young_to_middle %>%
            mutate(helper_vector = rep(seq(1:nrow(.)))) %>%
            unite(., Young, "Young", "helper_vector", remove = FALSE) %>%
            unite(., Middle, "Middle", "helper_vector", remove = FALSE) %>%
            dplyr::select(-helper_vector)

          rownames(outer_middle_to_old) <- unlist(bypass_duplicate_row$Middle)
          colnames(outer_middle_to_old) <- unlist(tmp$Bridgeness)
        }

        outer_middle_to_old <- outer_middle_to_old %>% mutate(MO = "MO")
        crossdatadf <- adjacency_to_2col(outer_middle_to_old) %>% mutate(Value_norm = as.numeric(scale(Value)))

        # Probabilistic Threshold selection
        max <- crossdatadf %>% dplyr::slice_max(Value_norm, n = 1)
        if (max$Value_norm[1] >= threshold) {
          outer_middle_to_old_2 <- crossdatadf %>% filter(Value_norm >= threshold)
        } else {
          outer_middle_to_old_2 <- crossdatadf %>% dplyr::slice_max(Value_norm, n = 1)
        }

        tmp_middle_to_old <- cbind(tmp %>% slice(1:nrow(outer_middle_to_old_2)) %>%
          dplyr::select(`1st_network`, Region), outer_middle_to_old_2) %>%
          plyr::rename(c("Value" = "Value_MO"))

        full_trajectory <- tmp_middle_to_old %>%
          merge(., bypass_duplicate_row %>% dplyr::select(c("Young", "Middle", "Value_YM")), by = c("Middle"))

        Cond_PMF_list[[i]] <- full_trajectory
      }

      if (colnames(composition)[4] == "Hub_consensus") {
        Cond_PMF_final <<- rbindlist(Cond_PMF_list) %>%
          # Giving the original labels back
          mutate(Young = ifelse(base::startsWith(Young, "Connector"), "Connector",
            ifelse(base::startsWith(Young, "Provincial"), "Provincial",
              ifelse(base::startsWith(Young, "Peripheral"), "Peripheral",
                ifelse(base::startsWith(Young, "Satellite"), "Satellite", Young)
              )
            )
          )) %>%
          mutate(Middle = ifelse(base::startsWith(Middle, "Connector"), "Connector",
            ifelse(base::startsWith(Middle, "Provincial"), "Provincial",
              ifelse(base::startsWith(Middle, "Peripheral"), "Peripheral",
                ifelse(base::startsWith(Middle, "Satellite"), "Satellite", Middle)
              )
            )
          )) %>%
          # Identify each region with a unique label
          mutate(helper_vector = rep(seq(nrow(.)))) %>%
          unite(., Region, c("Region", "helper_vector"), remove = FALSE) %>%
          # Transform  back to an alluvial format
          pivot_longer(
            cols = c("Young", "Middle", "Old"),
            names_to = "Age_group",
            values_to = "Hub_consensus"
          )

        Cond_PMF_final$Age_group <- factor(Cond_PMF_final$Age_group, levels = c("Young", "Middle", "Old"))
      } else {
        Cond_PMF_final <<- rbindlist(Cond_PMF_list) %>%
          mutate(Young = ifelse(base::startsWith(Young, "Global_Bridge"), "Global_Bridge",
            ifelse(base::startsWith(Young, "Local_Bridge"), "Local_Bridge",
              ifelse(base::startsWith(Young, "Super_Bridge"), "Super_Bridge",
                ifelse(base::startsWith(Young, "Not_a_Bridge"), "Not_a_Bridge", Young)
              )
            )
          )) %>%
          mutate(Middle = ifelse(base::startsWith(Middle, "Global_Bridge"), "Global_Bridge",
            ifelse(base::startsWith(Middle, "Local_Bridge"), "Local_Bridge",
              ifelse(base::startsWith(Middle, "Super_Bridge"), "Super_Bridge",
                ifelse(base::startsWith(Middle, "Not_a_Bridge"), "Not_a_Bridge", Middle)
              )
            )
          )) %>%
          # Identify each region with a unique label
          mutate(helper_vector = rep(seq(nrow(.)))) %>%
          unite(., Region, c("Region", "helper_vector"), remove = FALSE) %>%
          # Transform  back to an alluvial format
          pivot_longer(
            cols = c("Young", "Middle", "Old"),
            names_to = "Age_group",
            values_to = "Bridgeness"
          )

        Cond_PMF_final$Age_group <- factor(Cond_PMF_final$Age_group, levels = c("Young", "Middle", "Old"))
      }


      library(ggalluvial)
      if (colnames(composition)[4] == "Hub_consensus") {
        display_cluster <- Cond_PMF_final %>%
          group_by(Age_group, Hub_consensus) %>%
          summarize(s = n()) %>%
          arrange(Age_group, desc(Hub_consensus)) %>%
          .$Hub_consensus

        alluvial_cluster <- ggplot(
          Cond_PMF_final,
          aes(x = Age_group, stratum = Hub_consensus, alluvium = Region, fill = Hub_consensus)
        ) +
          geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
          geom_stratum(alpha = .8) +
          scale_x_discrete(expand = c(.1, .1)) +
          # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
          geom_text(stat = "stratum", label = display_cluster) +
          scale_fill_brewer(palette = "Oranges", direction = 1) +
          labs(title = paste0("z-scored Probability threshold: ", threshold)) +
          labs(x = "Age group") +
          theme_pubclean()

        alluvial_cluster
      } else {
        display_cluster <- Cond_PMF_final %>%
          group_by(Age_group, Bridgeness) %>%
          summarize(s = n()) %>%
          arrange(Age_group, desc(Bridgeness)) %>%
          .$Bridgeness

        alluvial_cluster <- ggplot(
          Cond_PMF_final,
          aes(x = Age_group, stratum = Bridgeness, alluvium = Region, fill = Bridgeness)
        ) +
          geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
          geom_stratum(alpha = .8) +
          scale_x_discrete(expand = c(.1, .1)) +
          # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
          geom_text(stat = "stratum", label = display_cluster) +
          scale_fill_brewer(palette = "Oranges", direction = 1) +
          labs(title = paste0("z-scored Probability threshold: ", threshold)) +
          labs(x = "Age group") +
          theme_pubclean()

        alluvial_cluster
      }
    }

    ############################################################################
    ############################################################################
    ############################################################################

    input_composition <- switch(plot_button_bis(),
      "Modular" = modular,
      "Interareal" = interareal
    )

    if (length(plot_button()) == 1) {
      RSN <- input$RSN_choice
      threshold <- sliderValue()
    } else if (length(plot_button()) > 1) {
      list_input <- list()
      for (i in 1:length(plot_button())) {
        list_input[[i]] <- plot_button()[i]
        tmp_bis <- rbindlist(lapply(list_input, as.data.table))
        RSN <- capture.output(cat(tmp_bis %>% as.matrix(), sep = "|"))
      }
      threshold <- sliderValue()
    } else if (length(plot_button()) < 1) {
      RSN <- "All"
      threshold <- sliderValue()
    }

    trajectory(input_composition, RSN, threshold)
  })
}

shinyApp(ui = ui, server = server)
