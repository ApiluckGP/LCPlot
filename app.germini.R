# This is a self-contained R Shiny web application for uploading and
# previewing different file formats (CSV, XLSX, TSV).

# Load the required libraries
library(shiny)
library(readr)  # For reading CSV and TSV
library(readxl) # For reading XLSX
library(dplyr)  # For data manipulation

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
  uploaded_data_rv <- reactiveVal(NULL)
  
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
    
    # Store the data in the reactive value
    uploaded_data_rv(data)
  })
  
  # Output the contents of the uploaded file to a table
  output$contents <- renderDataTable({
    uploaded_data_rv()
  })
  
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
