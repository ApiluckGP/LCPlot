library(shiny)
library(ggplot2)
library(dplyr)
library(plotly)
library(tidyr)

# UI
ui <- fluidPage(
  titlePanel("Interactive 96-Well Plate Layout with Clustering"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload qPCR Data (CSV)", accept = c(".csv")),
      textInput("alleleFAM", "Assign Allele Name for FAM:", value = "Allele 1"),
      textInput("alleleHEX", "Assign Allele Name for HEX:", value = "Allele 2"),
      actionButton("generate", "Generate Plate Layout"),
      br(),
      plotlyOutput("interactivePlateLayout", height = "800px"),
      width = 5
    ),
    mainPanel(
      plotOutput("amplificationPlot", height = "300px"),
      plotlyOutput("allelicDiscriminationPlot", height = "400px"),
      width = 7
    )
  )
)

# Server
server <- function(input, output, session) {
  process_data <- reactive({
    req(input$file1)
    data <- read.csv(input$file1$datapath)
    data %>%
      rowwise() %>%
      mutate(WELL = strsplit(as.character(WELL.POSITION), "\\|")) %>%
      unnest(cols = c(WELL)) %>%
      mutate(
        ROW = substr(WELL, 1, 1),
        COL = as.numeric(substr(WELL, 2, nchar(WELL))),
        Y_values = strsplit(as.character(Y), "\\|")
      ) %>%
      mutate(
        End_Fluorescence = as.numeric(sapply(Y_values, function(x) if (length(x) > 0) tail(x, 1) else NA))
      ) %>%
      select(-Y_values)
  })
  
  output$allelicDiscriminationPlot <- renderPlotly({
    req(input$generate)
    data <- process_data()
    
    allelic_data <- data %>%
      filter(DYE.NAME %in% c("FAM", "HEX")) %>%
      group_by(ROW, COL, SAMPLE.ID, WELL) %>%
      summarise(
        FAM = max(ifelse(DYE.NAME == "FAM", End_Fluorescence, NA), na.rm = TRUE),
        HEX = max(ifelse(DYE.NAME == "HEX", End_Fluorescence, NA), na.rm = TRUE),
        .groups = 'drop'
      )
    
    p <- ggplot(allelic_data, aes(x = FAM, y = HEX, label = paste(SAMPLE.ID, WELL))) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text(aes(label = paste(SAMPLE.ID, WELL)), hjust = -0.2, vjust = -0.2, size = 3) +
      labs(title = "Allelic Discrimination Plot", x = input$alleleFAM, y = input$alleleHEX) +
      theme_minimal()
    
    ggplotly(p, tooltip = c("x", "y", "label"))
  })
}

shinyApp(ui = ui, server = server)
