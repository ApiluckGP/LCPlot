# This is a self-contained R Shiny web application for uploading and
# previewing different file formats (CSV, XLSX, TSV).

# Load the required libraries
library(shiny)
library(readr)  # For reading CSV and TSV
library(readxl) # For reading XLSX
library(dplyr)  # For data manipulation
library(tidyr)  # For data transformation

# Define the user interface (UI) with a navigation menu and sidebar layout
ui <- navbarPage(
  title = "File Uploader App",
  
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
               h4("Preview of Uploaded Data"),
               p("The content of the most recently uploaded file will be displayed here."),
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
  processed_data_rv <- reactiveVal(NULL)
  
  # Observer for the "Load Data" button to read and process the file
  observeEvent(input$load_data, {
    req(input$data_file)
    
    file_type <- input$file_type
    
    # Read the file based on the selected type
    data <- switch(file_type,
                   "csv" = read_csv(input$data_file$datapath),
                   "xlsx" = read_excel(input$data_file$datapath),
                   "tsv" = read_tsv(input$data_file$datapath)
    )
    
    # --- Data Processing and Calculation Pipeline ---
    processed_data <- data %>%
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
      mutate(delta_Rn = Value - Baseline) %>%
      # Get the mean ΔRn for each sample and dye
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
  })
  
  # Output the contents of the uploaded file to a table
  output$contents <- renderDataTable({
    processed_data_rv()
  }, options = list(
    rownames = FALSE,
    columnDefs = list(list(className = 'dt-center', targets = '_all'))
  ))
  
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
