library(shiny)
library(ggplot2)
library(dplyr)
library(plotly)
library(tidyr)

# UI
ui <- fluidPage(
  titlePanel("Interactive 96-Well Plate Layout with Amplification Plots"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload qPCR Data (CSV)", accept = c(".csv")),
      textInput("alleleFAM", "Assign Allele Name for FAM (Blue):", value = "FAM"),
      textInput("alleleHEX", "Assign Allele Name for HEX (Red):", value = "HEX"),
      #checkboxInput("swapAxes", "Swap Discrimination Plot Axes", value = FALSE), # Add swap option
      actionButton("generate", "Generate Plate Layout"),
      br(),
      tags$h4("96-Well Plate Layout"),
      plotlyOutput("interactivePlateLayout", height = "800px"),
      width = 5
      #actionButton("generate", "Generate Plate Layout"),
      #br(),
      #tags$h4("96-Well Plate Layout"),
      #plotlyOutput("interactivePlateLayout", height = "450px"), # Interactive plate layout
      #width = 5 # Adjust sidebar width
    ),
    mainPanel(
      tags$h4("Amplification Plot"),
      plotOutput("amplificationPlot", height = "300px"), # Amplification plot
      tags$h4("Allelic Discrimination Plot"),
      plotlyOutput("allelicDiscriminationPlot", height = "400px"),
      width = 7 # Adjust main panel width
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive function to process the uploaded data
  process_data <- reactive({
    req(input$file1)
    data <- read.csv(input$file1$datapath)
    
    # Expand well positions into individual rows
    data <- data %>%
      rowwise() %>%
      mutate(WELL = strsplit(as.character(WELL.POSITION), "\\|")) %>%
      unnest(cols = c(WELL)) %>%
      mutate(
        ROW = substr(WELL, 1, 1), # Extract row (A-H)
        COL = as.numeric(substr(WELL, 2, nchar(WELL))) # Extract column (1-12)
      )
    data
  })
  
  safeAlleleFAM <- reactive({
    ifelse(input$alleleFAM == "", "FAM", input$alleleFAM)
  })
  # 
  safeAlleleHEX <- reactive({
    ifelse(input$alleleHEX == "", "HEX", input$alleleHEX)
  })
  
  # Generate interactive plate layout
  output$interactivePlateLayout <- renderPlotly({
    req(input$generate)
    data <- process_data()
    
    # Check if SAMPLE.ID exists in data
    if (!"SAMPLE.ID" %in% colnames(data)) {
      stop("The input data must contain a column named 'SAMPLE.ID'.")
    }
    
    # Create a blank 96-well plate structure
    plate <- expand.grid(
      ROW = LETTERS[1:8], # Rows A-H
      COL = 1:12          # Columns 1-12
    )
    
    # Merge with input data to map wells
    plate <- plate %>%
      left_join(data, by = c("ROW", "COL")) # Ensure SAMPLE.ID is merged
    
    # Parse fluorescence values for FAM and HEX
    plate <- plate %>%
      rowwise() %>%
      mutate(
        #Show by diff of max to min
        FAM = ifelse(DYE.NAME == "FAM", max(as.numeric(unlist(strsplit(as.character(Y), "\\|")))) - min(as.numeric(unlist(strsplit(as.character(Y), "\\|")))), NA),
        HEX = ifelse(DYE.NAME == "HEX", max(as.numeric(unlist(strsplit(as.character(Y), "\\|")))) - min(as.numeric(unlist(strsplit(as.character(Y), "\\|")))), NA)
      ) %>%
      group_by(ROW, COL, SAMPLE.ID) %>% # Keep SAMPLE.ID in the pipeline
      summarise(
        FAM = max(FAM, na.rm = TRUE),
        HEX = max(HEX, na.rm = TRUE)
      ) %>%

      mutate(
        Color = case_when(
          FAM >= 5.0 & HEX < 5.0 ~ paste0(safeAlleleFAM(), "/", safeAlleleFAM()),   # FAM significantly higher
          HEX >= 5.0 & FAM < 5.0 ~ paste0(safeAlleleHEX(), "/", safeAlleleHEX()),   # HEX significantly higher
          FAM >= 5.0 & HEX >= 5.0 ~ paste0(safeAlleleFAM(), "/", safeAlleleHEX()),   # Both increased
          FAM >= 0.0 & FAM < 5.0 & HEX >= 0 & HEX < 5.0 ~ "Undetermined",          # no significant increase
          FAM == -Inf | HEX == -Inf ~ "No Data"                                   # No Sample 
        )
      )
    
      color_mapping <- setNames(
        c("blue", "red", "green", "gray","white"),
        c(
          paste0(safeAlleleFAM(), "/", safeAlleleFAM()),
          paste0(safeAlleleHEX(), "/", safeAlleleHEX()),
          paste0(safeAlleleFAM(), "/", safeAlleleHEX()),
          "Undetermined",
          "No Data"
        )
      )
    # Generate ggplot with well information
    p <- ggplot(plate, aes(x = COL, y = ROW)) +
      geom_point(
        aes(fill = Color, text = paste0(
          "WELL: ", ROW, COL, "\n",
          "Sample ID: ", SAMPLE.ID, "\n",
          input$alleleFAM, ": ", FAM, "\n",
          input$alleleHEX, ": ", HEX
        )),
        shape = 21, size = 11, color = "dark gray" # shape line 
      ) +
      geom_text(
        aes(label = stringr::str_trunc(SAMPLE.ID, width = 8)), # Truncate long Sample IDs
        size = 2, color = "dark gray", na.rm = TRUE # Letter Color
      ) +
      
      scale_color_manual(values = color_mapping) +
      
      scale_fill_identity() + # Use calculated colors directly
      scale_y_discrete(limits = rev(LETTERS[1:8])) + # Reverse rows to match plate view
      scale_x_continuous(breaks = 1:12) +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        #legend.position = "bottom",
        plot.title = element_text(size = 14, hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      #labs(title = "96-Well Plate Layout with Sample IDs")
      labs()
    # Convert ggplot to plotly for interactivity
    p <- ggplotly(p, tooltip = "text")
    
    #Register the plotly_click event
    event_register(p, "plotly_click")
    
    return(p)
  })
  
  
  # Observe clicks on the interactive plot
  selected_amplification <- reactiveVal(NULL)
  observeEvent(event_data("plotly_click"), {
    click_data <- event_data("plotly_click")
    if (!is.null(click_data)) {
      row_clicked <- LETTERS[9 - round(click_data$y)] # Reverse mapping for rows
      col_clicked <- round(click_data$x)
      
      # Filter the data for the clicked well
      data <- process_data()
      selected <- data %>%
        filter(ROW == row_clicked, COL == col_clicked)
      
      selected_amplification(selected)
    }
  })
  
    # Generate the amplification plot for the selected well
  output$amplificationPlot <- renderPlot({
    req(selected_amplification())
    selected <- selected_amplification()
    
    # Ensure X and Y values exist for the plot
    if (!is.null(selected$Y)) {
      # Parse cycle order (X) and fluorescence values (Y)
      cycle <- 1:40
      
      # Parse fluorescence values from Y
      #fluorescence <- as.numeric(unlist(strsplit(as.character(selected$Y), "\\|")))
      y_values <- strsplit(as.character(selected$Y), "\\|")
      dyes <- selected$DYE.NAME
      
      # Combine into a data frame
      plot_data <- data.frame(
        #Cycle = cycle[1:length(fluorescence)], # Ensure cycle matches the length of Y
        Cycle = rep(cycle, times = length(dyes)),
        Fluorescence = as.numeric(unlist(y_values)),
        #Fluorescence = fluorescence,
        #Dye = selected$DYE.NAME
        Dye = rep(unlist(dyes), each = length(cycle))
      )
      
      # Generate amplification plot
      ggplot(plot_data, aes(x = Cycle, y = Fluorescence, color = Dye)) +
        geom_line(size = 1.2)+ 
        scale_color_manual(values = c("FAM" = "blue", "HEX" = "red")) +
        labs(
          title = paste("Amplification Plot for Well", selected$WELL),
          x = "Cycle",
          y = "Fluorescence Intensity",
          color = "Dye"
        ) +
        scale_x_continuous(breaks = seq(0, 40, 1), limits = c(1, 40)) + # Fixed cycle range
        #scale_x_continuous(breaks = 1:40, limits = c(1, 40)) + # Show every cycle
        scale_y_continuous(breaks = seq(0, 40, 2))+
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.position = "bottom"
        )
    } else {
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "No Amplification Data", size = 5, hjust = 0.5) +
        theme_void()
    }
  }) # Best Version 
  
  output$allelicDiscriminationPlot <- renderPlotly({
    req(input$generate)
    data <- process_data()
    
    allelic_data <- data %>%
      group_by(ROW, COL) %>%
      summarise(
        FAM = max(ifelse(DYE.NAME == "FAM", max(as.numeric(unlist(strsplit(as.character(Y), "\\|")))), 0), na.rm = TRUE),
        HEX = max(ifelse(DYE.NAME == "HEX", max(as.numeric(unlist(strsplit(as.character(Y), "\\|")))), 0), na.rm = TRUE)
      ) %>%
      mutate(
        Category = case_when(
          FAM > HEX & FAM > 1.5 * HEX ~ "FAM High",
          HEX > FAM & HEX > 1.5 * FAM ~ "HEX High",
          FAM > 1.0 & HEX > 1.0 ~ "Both High",
          TRUE ~ "Low Signal"
        )
      )
    
    # Determine axes based on user choice
    x_axis <- if (input$swapAxes) "HEX" else "FAM"
    y_axis <- if (input$swapAxes) "FAM" else "HEX"
    x_label <- if (input$swapAxes) input$alleleHEX else input$alleleFAM
    y_label <- if (input$swapAxes) input$alleleFAM else input$alleleHEX
    
    p <- ggplot(allelic_data, aes_string(x = x_axis, y = y_axis, color = "Category")) +
      geom_point(size = 3, alpha = 0.8) +
    
    #p <- ggplot(allelic_data, aes(x = FAM, y = HEX, color = Category)) +
     # geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(
        values = c("FAM High" = "blue", "HEX High" = "red", "Both High" = "green", "Low Signal" = "gray")
      ) +
      labs(
        title = "Allelic Discrimination Plot",
        #x = "FAM Fluorescence Intensity",
        #y = "HEX Fluorescence Intensity",
        x = paste(input$alleleFAM, "(FAM)"),
        y = paste(input$alleleHEX, "(HEX)"),
        color = "Category"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "bottom"
      )
    
    ggplotly(p, tooltip = c("x", "y", "color"))
  })
}

# Run the app
shinyApp(ui = ui, server = server)
