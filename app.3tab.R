library(shiny)
library(ggplot2)
library(dplyr)
library(plotly)
library(qpcR)

ui <- fluidPage(
  titlePanel("Allelic Discrimination Plot with Clustering and Tabs"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload qPCR Data (CSV)"),
      fluidRow(
        #column(6, uiOutput("alleleAssignment1")),  # Allele 1 assignment (left)
        column(4, textInput("allele1", "Assign Allele Dye 1:", value = "")),  # Free text input for Allele 1
        column(4, uiOutput("dyeSelection1"))  # Dye 1 selection (right)
      ),
      fluidRow(
        #column(6, uiOutput("alleleAssignment2")),  # Allele 2 assignment (left)
        column(4, textInput("allele2", "Assign Allele Dye 2:", value = "")),   # Free text input for Allele 2
        column(4, uiOutput("dyeSelection2"))  # Dye 2 selection (right)
      ),
      actionButton("plot", "Generate Plot"),
      #column (6, column( 2, checkboxInput("cluster", "Apply K-Means Clustering (k=3)", value = TRUE))),  # Clustering checkbox enabled by default
      checkboxInput("cluster", "Apply K-Means Clustering (k=3)", value = TRUE)
      , width = 4),  # Reduced sidebar width (4 out of 12 columns)
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data Table", tableOutput("fileContents")),  # Table for CSV input
        tabPanel("Allelic Discrimination Plot", plotlyOutput("interactivePlot")),  # 2D plot with clustering and circles
        tabPanel("Amplification Plot",
                 selectInput("amplificationDye", "Select Dye for Amplification Plot", choices = NULL),  # Dye selection dropdown
                 plotlyOutput("amplificationPlot")  # New amplification plot tab
        )
      )
      , width = 8)  # Main panel takes the remaining space (8 out of 12 columns)
  )
)

server <- function(input, output, session) {
  
  # Reactive function to load and process the data
  data <- reactive({
    req(input$file1)  # Ensure a file is uploaded
    df <- tryCatch({
      read.csv(input$file1$datapath, stringsAsFactors = FALSE)
    }, error = function(e) {
      showNotification("Error loading file", type = "error")
      return(NULL)
    })
    
    # Automatically clean column names to handle spaces and special characters
    colnames(df) <- make.names(colnames(df), unique = TRUE)
    
    # Ensure that necessary columns exist
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
    
    
    # Show the first 10 rows of the data after loading
    output$fileContents <- renderTable({
      df <- data()
      
      # Rename the columns
      colnames(df) <- gsub("DYE.NAME", "DYE NAME", colnames(df))
      colnames(df) <- gsub("SAMPLE.ID", "SAMPLE ID", colnames(df))
      colnames(df) <- gsub("WELL.POSITION", "WELL", colnames(df))
      
      # Filter out rows where SAMPLE ID is missing
      df <- df[!is.na(df$`SAMPLE ID`) & df$`SAMPLE ID` != "", ]
      
      # Remove the X and Y columns
      #df <- df %>% select(-X, -Y)
      
      # Show the updated table with renamed columns
      head(df, 20)  # Display only the first 20 rows
    })
    
    # Update the dropdown for the amplification plot dye selection dynamically
    updateSelectInput(session, "amplificationDye", choices = unique(df$DYE.NAME))
    
    return(df)
  })
  
  
  # Dynamically generate dropdowns for dye selection (Dye 1)
  output$dyeSelection1 <- renderUI({
    df <- data()
    if (is.null(df)) return(NULL)  # If no data, don't generate UI
    
    # Get unique dyes present in the dataset
    dyes <- unique(df$DYE.NAME)
    
    # Ensure dyes exist before rendering the UI
    if (length(dyes) == 0) {
      showNotification("No dyes found in the data", type = "error")
      return(NULL)
    }
    
    # Dye 1 selection dropdown
    selectInput("dye1", "Dye 1", choices = dyes, selected = dyes[1])
  })
  
  # Dynamically generate dropdowns for dye selection (Dye 2)
  output$dyeSelection2 <- renderUI({
    df <- data()
    if (is.null(df)) return(NULL)  # If no data, don't generate UI
    
    # Get unique dyes present in the dataset
    dyes <- unique(df$DYE.NAME)
    
    # Ensure dyes exist before rendering the UI
    if (length(dyes) == 0) {
      showNotification("No dyes found in the data", type = "error")
      return(NULL)
    }
    
    # Dye 2 selection dropdown
    selectInput("dye2", "Dye 2", choices = dyes, selected = dyes[2])
  })
  
  # Dynamically generate dropdowns for allele assignment (Allele 1)
  output$alleleAssignment1 <- renderUI({
    # Define allele options
    alleles <- c("A", "T", "C", "G")
    
    # Allele assignment dropdown for Dye 1
    selectInput("allele1", "Assign Allele for Dye 1", choices = alleles, selected = "A")
  })
  
  # Dynamically generate dropdowns for allele assignment (Allele 2)
  output$alleleAssignment2 <- renderUI({
    # Define allele options
    alleles <- c("A", "T", "C", "G")
    
    # Allele assignment dropdown for Dye 2
    selectInput("allele2", "Assign Allele for Dye 2", choices = alleles, selected = "T")
  })
  
  # Plot generation based on user input and selected dyes and allele assignment
  observeEvent(input$plot, {
    df <- data()
    req(df, input$dye1, input$dye2, input$allele1, input$allele2)  # Ensure file, dyes, and alleles are selected
    
    # Check if selected dyes are available
    if (!(input$dye1 %in% df$DYE.NAME) | !(input$dye2 %in% df$DYE.NAME)) {
      showNotification("Selected dyes are not available in the data", type = "error")
      return(NULL)
    }
    
    # Extract data for the selected dyes
    dye1_data <- df %>% filter(DYE.NAME == input$dye1)
    dye2_data <- df %>% filter(DYE.NAME == input$dye2)
    
    # Extract the end-point fluorescence (last value) from the X and Y columns
    dye1_data$End_Fluorescence <- sapply(strsplit(dye1_data$Y, "\\|"), function(x) as.numeric(tail(x, 1)))
    dye2_data$End_Fluorescence <- sapply(strsplit(dye2_data$Y, "\\|"), function(x) as.numeric(tail(x, 1)))
    
    # Combine both dyes into a single data frame for the 2D plot
    combined_data <- data.frame(
      SAMPLE.ID = dye1_data$SAMPLE.ID,
      WELL.POSITION = dye1_data$WELL.POSITION,
      Dye1_Fluorescence = dye1_data$End_Fluorescence,
      Dye2_Fluorescence = dye2_data$End_Fluorescence
    )
    
    # Apply K-means clustering to classify as Homozygous (REF), Heterozygous, and Homozygous (Variant)
    if (input$cluster) {
      kmeans_result <- kmeans(combined_data[, c("Dye1_Fluorescence", "Dye2_Fluorescence")], centers = 3)
      combined_data$Cluster <- as.factor(kmeans_result$cluster)
      
      # Assign names to clusters based on fluorescence values
      combined_data$Cluster <- recode(combined_data$Cluster, 
                                      `1` = paste("Homozygous (", input$allele1, ")", sep = ""),
                                      `2` = paste("Heterozygous (", input$allele1, "/", input$allele2, ")", sep = ""),
                                      `3` = paste("Homozygous (", input$allele2, ")", sep = ""))
    }
    
    # Generate 2D plot using plotly with allele and dye name in axis labels
    output$interactivePlot <- renderPlotly({
      plot <- ggplot(combined_data, aes(x = Dye1_Fluorescence, y = Dye2_Fluorescence, 
                                        text = SAMPLE.ID, color = Cluster)) +
        geom_point(size = 4) +
        labs(title = "Allelic Discrimination Plot", 
             x = paste("Allele", input$allele1, "(", input$dye1, ")"),  # Show Allele and Dye 1
             y = paste("Allele", input$allele2, "(", input$dye2, ")")) +  # Show Allele and Dye 2
        theme_minimal()
      
      # Add circles (ellipse) around clusters with boundary
      plot <- plot + stat_ellipse(aes(color = Cluster), type = "norm", level = 0.95, linetype = "dashed")
      
      # Convert ggplot to plotly for interactivity (hover tooltip with SampleID)
      ggplotly(plot, tooltip = "text")
    })
    
    # Generate amplification plot using plotly
    output$amplificationPlot <- renderPlotly({
      # Ensure the selected dye data is available
      req(input$amplificationDye)
      
      # Filter data for the selected dye
      dye_data <- df %>% filter(DYE.NAME == input$amplificationDye)
      
      # Extract the cycle and fluorescence data for amplification from columns X and Y
      dye_cycles <- strsplit(dye_data$X, "\\|")  # Split the cycle data
      dye_fluorescence <- strsplit(dye_data$Y, "\\|")  # Split the fluorescence data
      
      # Ensure the length of cycles and fluorescence values match for each sample
      dye_cycles <- lapply(dye_cycles, as.numeric)
      dye_fluorescence <- lapply(dye_fluorescence, as.numeric)
      
      # Create a data frame for the amplification plot
      amplification_data <- data.frame(
        Cycles = unlist(dye_cycles),
        Fluorescence = unlist(dye_fluorescence),
        Sample = rep(dye_data$`SAMPLE.ID`, sapply(dye_cycles, length))  # Add the sample IDs for labeling
      )
      
      # Generate the amplification plot
      plot <- ggplot(amplification_data, aes(x = Cycles, y = Fluorescence, color = Sample)) +
        geom_line() +
        labs(title = paste("Amplification Plot for", input$amplificationDye),
             x = "Cycles", y = "Fluorescence") +
        theme_minimal()
      
      # Convert ggplot to plotly for interactivity
      ggplotly(plot)
    })
    
  })
}

shinyApp(ui = ui, server = server)
