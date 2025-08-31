library(shiny)
library(ggplot2)
library(dplyr)

ui <- fluidPage(
  titlePanel("Allelic Discrimination Plot"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload qPCR Data (CSV)"),
      uiOutput("dyeSelection"),  # Dynamic UI for selecting dyes
      actionButton("plot", "Generate Plot")
    ),
    
    mainPanel(
      tableOutput("fileContents"),  # Table to display file contents
      plotOutput("allelicDiscriminationPlot")
    )
  )
)

server <- function(input, output) {
  
  # Reactive function to load and process the data
  data <- reactive({
    req(input$file1)  # Ensure a file is uploaded
    df <- tryCatch({
      read.csv(input$file1$datapath, stringsAsFactors = FALSE)
    }, error = function(e) {
      showNotification("Error loading file", type = "error")
      return(NULL)
    })
    
    # Print column names to console for debugging
    print(colnames(df))
    
    # Automatically clean column names to handle spaces and special characters
    colnames(df) <- make.names(colnames(df), unique = TRUE)
    
    # Ensure that cleaned column names exist
    if (!"SAMPLE.ID" %in% names(df)) {
      showNotification("Column 'SAMPLE ID' not found. Please check the file.", type = "error")
      return(NULL)
    }
    if (!"WELL.POSITION" %in% names(df)) {
      showNotification("Column 'WELL POSITION' not found. Please check the file.", type = "error")
      return(NULL)
    }
    if (!"DYE.NAME" %in% names(df)) {
      showNotification("Column 'DYE NAME' not found. Please check the file.", type = "error")
      return(NULL)
    }
    
    # Filter out rows with blank SAMPLE.ID
    df <- df %>% filter(!is.na(SAMPLE.ID))
    
    # Show the first 10 rows of the data after loading
    output$fileContents <- renderTable({
      if (!is.null(df)) {
        head(df, 10)
      }
    })
    
    return(df)
  })
  
  # Dynamically generate dropdowns for dye selection based on data
  output$dyeSelection <- renderUI({
    df <- data()
    if (is.null(df)) return(NULL)  # If no data, don't generate UI
    
    # Get unique dyes present in the dataset
    dyes <- unique(df$DYE.NAME)
    
    # Check if any dyes are available
    if (length(dyes) == 0) {
      showNotification("No dyes found in the data", type = "error")
      return(NULL)
    }
    
    # Only allow the selection of dyes that are present in the dataset
    tagList(
      selectInput("dye1", "Select Dye 1 (Allele 1)", choices = dyes, selected = dyes[1]),
      selectInput("dye2", "Select Dye 2 (Allele 2)", choices = dyes, selected = dyes[2])
    )
  })
  
  # Plot generation based on user input and selected dyes
  observeEvent(input$plot, {
    df <- data()
    req(df, input$dye1, input$dye2)  # Ensure file and dyes are selected
    
    # Check if selected dyes are available
    if (!(input$dye1 %in% df$DYE.NAME) | !(input$dye2 %in% df$DYE.NAME)) {
      showNotification("Selected dyes are not available in the data", type = "error")
      return(NULL)
    }
    
    # Extract data for the selected dyes
    dye1_data <- df %>% filter(DYE.NAME == input$dye1)
    dye2_data <- df %>% filter(DYE.NAME == input$dye2)
    
    # Extract the final cycle fluorescence (end-point)
    dye1_data$End_Fluorescence <- sapply(strsplit(dye1_data$Y, "\\|"), function(x) as.numeric(tail(x, 1)))
    dye2_data$End_Fluorescence <- sapply(strsplit(dye2_data$Y, "\\|"), function(x) as.numeric(tail(x, 1)))
    
    # Merge the data by WELL.POSITION for plotting
    merged_data <- merge(dye1_data[, c("SAMPLE.ID", "WELL.POSITION", "End_Fluorescence")], 
                         dye2_data[, c("SAMPLE.ID", "WELL.POSITION", "End_Fluorescence")], 
                         by = "WELL.POSITION", suffixes = c(paste0("_", input$dye1), paste0("_", input$dye2)))
    
    # Check if merged_data has the expected columns
    print(colnames(merged_data))
    
    # Generate the plot with backticks to handle periods in column names
    output$allelicDiscriminationPlot <- renderPlot({
      ggplot(merged_data, aes(x = merged_data[[paste0("End_Fluorescence_", input$dye1)]], 
                              y = merged_data[[paste0("End_Fluorescence_", input$dye2)]], 
                              color = merged_data[["SAMPLE.ID"]])) +  # Handle column name with period
        geom_point(size = 4) +
        labs(title = "Allelic Discrimination Plot", 
             x = paste(input$dye1, "Fluorescence"), 
             y = paste(input$dye2, "Fluorescence")) +
        theme_minimal()
    })
  })
}

shinyApp(ui = ui, server = server)
