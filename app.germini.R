# This is a self-contained R Shiny web application for uploading and
# previewing different file formats (CSV, XLSX, TSV).

# Load the required libraries
library(shiny)
library(readr)  # For reading CSV and TSV
library(readxl) # For reading XLSX
library(dplyr)  # For data manipulation
library(tidyr)  # For separating rows

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
               
               # Load data button
               actionButton("load_data", "Load Data", class = "btn-primary")
             ),
             
             # Main panel for displaying the uploaded data
             mainPanel(
               h4("Preview of Processed Data"),
               p("The organized and calculated data will be displayed here."),
               dataTableOutput("contents")
             )
           )
  ),
  
  # A simple, empty tab for future use
  tabPanel("Data Analysis",
           p("This is where the analysis and visualizations will go in future steps.")
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # A reactive value to hold the most recently uploaded data
  uploaded_data_rv <- reactiveVal(NULL)
  
  # Observer for the "Load Data" button to read and process the file
  observeEvent(input$load_data, {
    req(input$data_file)
    
    file_type <- input$file_type
    
    # Read the file based on the selected type
    df <- switch(file_type,
                 "csv" = read_csv(input$data_file$datapath),
                 "xlsx" = read_excel(input$data_file$datapath),
                 "tsv" = read_tsv(input$data_file$datapath)
    )
    
    # --- Data Processing and Calculation ---
    
    # 1. Organize data to different line for hex and fam of same SAMPLE ID or WELL position.
    processed_data <- df %>%
      select(`DYE NAME`, `SAMPLE ID`, `WELL POSITION`, X, Y) %>%
      pivot_longer(cols = c(X, Y), names_to = "measurement", values_to = "value") %>%
      separate_rows(value, sep = "\\|", convert = TRUE) %>%
      group_by(`DYE NAME`, `SAMPLE ID`, `WELL POSITION`) %>%
      mutate(Cycle = 1:n()) %>%
      pivot_wider(names_from = "measurement", values_from = "value")
    
    # 2. Calculate baseline fluorescence of FAM and HEX based on first 7 cycles
    baseline_data <- processed_data %>%
      filter(Cycle <= 7) %>%
      group_by(`DYE NAME`, `SAMPLE ID`, `WELL POSITION`) %>%
      summarise(Baseline = mean(Y, na.rm = TRUE), .groups = 'drop')
    
    # 3. Calculate ΔRn (FAM) and ΔRn (HEX)
    final_data <- processed_data %>%
      left_join(baseline_data, by = c("DYE NAME", "SAMPLE ID", "WELL POSITION")) %>%
      mutate(delta_Rn = Y - Baseline) %>%
      # Pivot wider to get separate columns for FAM and HEX delta_Rn
      pivot_wider(id_cols = c(`SAMPLE ID`, `WELL POSITION`, Cycle),
                  names_from = `DYE NAME`,
                  values_from = delta_Rn,
                  names_prefix = "Delta_Rn_") %>%
      rename("ΔRn (FAM)" = Delta_Rn_FAM, "ΔRn (HEX)" = Delta_Rn_HEX) %>%
      select(`SAMPLE ID`, `WELL POSITION`, `ΔRn (FAM)`, `ΔRn (HEX)`) %>%
      distinct() # Ensure unique rows
    
    # Store the processed data in the reactive value
    uploaded_data_rv(final_data)
  })
  
  # Output the contents of the uploaded file to a table
  output$contents <- renderDataTable({
    uploaded_data_rv()
  })
  
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
