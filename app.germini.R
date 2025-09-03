# This is a self-contained R Shiny web application for uploading and
# previewing different file formats (CSV, XLSX, TSV).

# Load the required libraries
library(shiny)
library(readr)  # For reading CSV and TSV
library(readxl) # For reading XLSX
library(dplyr)  # For data manipulation
library(tidyr)  # For data transformation
library(ggplot2) # For plotting and visualizations

# Define the user interface (UI) with a navigation menu and sidebar layout
ui <- navbarPage(
  title = "qPCR Data Analyzer",
  
  # The first tab, dedicated to file uploading and data preview
  tabPanel("Load Input",
           sidebarLayout(
             sidebarPanel(
               h4("Upload Your Data"),
               p("Select a file type and then browse for your file."),
               
               # Radio buttons to select file format
               radioButtons("file_type", "Select File Type:",
                            choices = c("CSV" = "csv", "XLSX" = "xlsx", "TSV" = "tsv"),
                            selected = "csv",
                            inline = TRUE),
               
               # Single file input for browsing
               fileInput("data_file", "Browse for File:",
                         multiple = FALSE,
                         accept = c(".csv", ".xlsx", ".tsv")),
               
               # Action buttons for user flow
               actionButton("load_data", "Load Data", class = "btn-primary"),
               hr(),
               actionButton("analyze_data", "Analyze Data", class = "btn-success")
             ),
             
             # Main panel for displaying the uploaded data
             mainPanel(
               h4("Preview of Uploaded Data"),
               p("The content of the most recently uploaded file will be displayed here."),
               dataTableOutput("contents")
             )
           )
  ),
  
  # A tab for data analysis and visualizations
  tabPanel("Data Analysis",
           
           fluidRow(
             # Left column (50% width) for the plate layout
             column(6,
                    fluidRow(
                      column(6, h4(strong("Plate Layout:"))),
                      column(6, selectInput("plate_size", NULL,
                                            choices = c("96-Wells Plate" = "96", "384-Wells Plate" = "384"),
                                            selected = "96"))
                    ),
                    # Add a click event to the plot output
                    plotOutput("plate_layout_plot", click = "plate_click"),
                    
                    fluidRow(
                      column(12,
                             # Text input for searching Sample ID
                             textInput("search_sample", "Search Sample ID:"),
                             # Checkbox group for selecting wells
                             checkboxGroupInput("selected_wells", "Select Well Position or Sample:")
                      )
                    )
             ),
             # Right column (50% width) for the other plots
             column(6,
                    # Top-right plot (25% total width)
                    fluidRow(
                      column(12,
                             h4("ΔRn Scatter Plot"),
                             plotOutput("scatter_plot")
                      )
                    ),
                    # Bottom-right plot (25% total width)
                    fluidRow(
                      column(12,
                             h4("Amplification Plot (for selected samples)"),
                             p("This section will be developed in a future step.")
                      )
                    )
             )
           ),
           
           hr()
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # Reactive values to hold the data at different processing stages
  raw_data_rv <- reactiveVal(NULL)
  processed_data_rv <- reactiveVal(NULL)
  full_cycle_data_rv <- reactiveVal(NULL)
  
  # Reactive value to store the selected wells
  selected_wells_rv <- reactiveVal(character(0))
  
  # Observer for the "Load Data" button to read the file
  observeEvent(input$load_data, {
    req(input$data_file)
    
    file_type <- input$file_type
    
    # Read the file based on the selected type
    data <- switch(file_type,
                   "csv" = read_csv(input$data_file$datapath, show_col_types = FALSE),
                   "xlsx" = read_excel(input$data_file$datapath),
                   "tsv" = read_tsv(input$data_file$datapath, show_col_types = FALSE)
    )
    
    # Store the raw data for later processing
    raw_data_rv(data)
    
    # Display a preliminary preview of the raw data (without processing)
    output$contents <- renderDataTable({
      req(raw_data_rv())
      head(raw_data_rv()) # Show only the head to avoid performance issues
    }, options = list(
      rownames = FALSE,
      pageLength = 10
    ))
  })
  
  # Observer for the "Analyze Data" button to perform all calculations and plotting
  observeEvent(input$analyze_data, {
    req(raw_data_rv())
    
    # The data processing pipeline
    data <- raw_data_rv() %>%
      # Organize data to have one row per cycle and dye
      pivot_longer(cols = c(X, Y), names_to = "Measurement", values_to = "Value") %>%
      separate_rows(Value, sep = "\\|", convert = TRUE) %>%
      select(-Measurement) %>% # We don't need this column anymore
      # Add a cycle number
      group_by(`SAMPLE ID`, `WELL POSITION`, `DYE NAME`) %>%
      mutate(Cycle = row_number()) %>%
      ungroup() %>%
      # Calculate the baseline for each sample (mean of the first 7 cycles)
      group_by(`SAMPLE ID`, `WELL POSITION`, `DYE NAME`) %>%
      mutate(Baseline = mean(Value[Cycle <= 7], na.rm = TRUE)) %>%
      ungroup() %>%
      # Calculate Delta Rn (ΔRn)
      mutate(delta_Rn = Value - Baseline)
    
    # Store the full cycle data for the fluorescence plot
    full_cycle_data_rv(data)
    
    # Get the mean ΔRn for each sample and dye
    processed_data <- data %>%
      group_by(`SAMPLE ID`, `WELL POSITION`, `DYE NAME`) %>%
      summarise(Mean_Delta_Rn = mean(delta_Rn, na.rm = TRUE), .groups = 'drop') %>%
      # Pivot wider to get separate columns for FAM and HEX ΔRn
      pivot_wider(
        names_from = `DYE NAME`,
        values_from = Mean_Delta_Rn
      ) %>%
      # Rename the new columns to the desired format
      rename(`ΔRn (FAM)` = FAM, `ΔRn (HEX)` = HEX) %>%
      # Round the values in the ΔRn columns to 5 decimal places
      mutate(`ΔRn (FAM)` = round(`ΔRn (FAM)`, 5),
             `ΔRn (HEX)` = round(`ΔRn (HEX)`, 5)) %>%
      # Select and order the columns for the preview table
      select(`WELL POSITION`, `SAMPLE ID`, `ΔRn (FAM)`, `ΔRn (HEX)`)
    
    # Store the processed data in the reactive value
    processed_data_rv(processed_data)
    
    # Update the preview table with the final processed data
    output$contents <- renderDataTable({
      processed_data_rv()
    }, options = list(
      rownames = FALSE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'))
    ))
  })
  
  # Create a reactive expression for checkbox choices
  filtered_wells <- reactive({
    req(processed_data_rv())
    
    # Get all unique well positions and their corresponding Sample IDs
    all_wells <- processed_data_rv() %>%
      select(`WELL POSITION`, `SAMPLE ID`)
    
    # Filter based on search input
    if (input$search_sample != "") {
      all_wells <- all_wells %>%
        filter(grepl(input$search_sample, `SAMPLE ID`, ignore.case = TRUE))
    }
    
    # Create a named vector for the checkbox choices
    setNames(all_wells$`WELL POSITION`, paste(all_wells$`WELL POSITION`, all_wells$`SAMPLE ID`))
  })
  
  # Update the checkbox group when the filtered_wells reactive changes
  observe({
    updateCheckboxGroupInput(session, "selected_wells",
                             choices = filtered_wells()
    )
  })
  
  # Observer to handle clicks on the plate layout plot
  observeEvent(input$plate_click, {
    req(processed_data_rv())
    
    # Use nearPoints to find the closest well to the click
    clicked_point <- nearPoints(processed_data_rv(), input$plate_click, allRows = FALSE, threshold = 10, maxpoints = 1)
    
    # Update the reactive value with the selected well position
    if (nrow(clicked_point) > 0) {
      well <- clicked_point$`WELL POSITION`[1]
      current_selections <- input$selected_wells
      if (well %in% current_selections) {
        # Remove from selection
        selected_wells_rv(current_selections[current_selections != well])
      } else {
        # Add to selection
        selected_wells_rv(c(current_selections, well))
      }
    }
  })
  
  # Update the checkbox group input based on the reactive selected_wells_rv
  observe({
    updateCheckboxGroupInput(session, "selected_wells",
                             selected = selected_wells_rv()
    )
  })
  
  # Render the Plate Layout plot
  output$plate_layout_plot <- renderPlot({
    req(processed_data_rv())
    
    # Get selected plate size
    plate_size <- input$plate_size
    
    if (plate_size == "96") {
      # 96-well plate dimensions (8 rows x 12 columns)
      rows <- rev(LETTERS[1:8]) # Reverse for correct plotting order
      cols <- 1:12
    } else {
      # 384-well plate dimensions (16 rows x 24 columns)
      rows <- rev(LETTERS[1:16]) # Reverse for correct plotting order
      cols <- 1:24
    }
    
    # Create a full grid of well positions
    full_grid <- expand.grid(Row = rows, Col = cols) %>%
      mutate(Well = paste0(Row, Col))
    
    # Join with the processed data
    plot_data <- left_join(full_grid, processed_data_rv(), by = c("Well" = "WELL POSITION")) %>%
      # Create a combined label for text
      mutate(label_text = paste0("FAM:", `ΔRn (FAM)`, "\nHEX:", `ΔRn (HEX)`))
    
    # Create the base plate layout plot
    p <- ggplot(plot_data, aes(x = as.factor(Col), y = as.factor(Row))) +
      geom_tile(aes(fill = `ΔRn (FAM)`), color = "black", linewidth = 0.5) +
      geom_text(aes(label = label_text), size = 3, color = "black") + 
      scale_fill_gradient(low = "yellow", high = "red", na.value = "lightgray", name = "FAM ΔRn") +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(angle = 90, hjust = 1, vjust = 0.5),
        aspect.ratio = ifelse(plate_size == "96", 8/12, 16/24),
        panel.grid = element_blank()
      )
    
    # If a well is selected, highlight it
    if (!is.null(input$selected_wells)) {
      selected_data <- plot_data %>% filter(Well %in% input$selected_wells)
      p <- p + geom_tile(data = selected_data, fill = NA, color = "red", linewidth = 2)
    }
    
    p
  })
  
  # Render the Scatter plot
  output$scatter_plot <- renderPlot({
    req(processed_data_rv())
    
    # Filter the data based on selected wells
    if (is.null(input$selected_wells)) {
      plot_data <- processed_data_rv()
    } else {
      plot_data <- processed_data_rv() %>%
        filter(`WELL POSITION` %in% input$selected_wells)
    }
    
    ggplot(plot_data, aes(x = `ΔRn (HEX)`, y = `ΔRn (FAM)`)) +
      geom_point(aes(color = `SAMPLE ID`), size = 3) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      labs(title = "FAM vs. HEX ΔRn Values",
           x = "ΔRn (HEX)",
           y = "ΔRn (FAM)",
           color = "Sample ID") +
      theme_minimal()
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
