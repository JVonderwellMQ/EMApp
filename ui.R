library(shiny)

fluidPage(
  
  titlePanel("EM Webapp"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("inputfile", label = h3("Upload your CSV file")),
      textOutput("filename"),
      
      hr(),
      numericInput("numModes", label = h3("Number of modes"), value = 1),
      hr(),
      
      h4("EM Results"),
      verbatimTextOutput("em_results")
    ),
    
    mainPanel(
      h3("Data Preview"),
      tableOutput("preview"),
      
      hr(),
      
      plotOutput("hist_em")
    )
  )
)
