library(shiny)

fluidPage(
  
  titlePanel("EM Webapp"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("inputfile", label = h3("Upload your CSV file")),
      textOutput("filename"),
      
      hr(),
      radioButtons("modeSelection", label = h3("How should the number of modes be selected"),
                   choices = list("AIC" = 1, "BIC" = 2, "Manual" = 3), 
                   selected = 1),
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
      
      plotOutput("hist_em"),
      
      hr(),
      
      plotOutput("ic_plot")
    )
  )
)
