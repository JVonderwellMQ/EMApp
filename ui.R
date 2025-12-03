library(shiny)

fluidPage(
  
  titlePanel("EM Webapp"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("inputfile", label = h3("Upload your CSV file")),
      textOutput("filename"),
      numericInput("column", label = h3("Which column to fit?"), value = 1),
      hr(),
      radioButtons("modeSelection", label = h3("How should the number of modes be selected"),
                   choices = list("AIC" = 1, "BIC" = 2, "Manual" = 3), 
                   selected = 1),
      numericInput("numModes", label = h3("Number of modes"), value = 1),
      hr(),
      
      sliderInput("iterSelect", "Iteration to view:",
                  min = 1, max = 1, value = 1, step = 1), # max updated dynamically,
      checkboxInput("nonParametric", label = "Nonparametric estimation", value = FALSE)
    ),
    
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Histogram", plotOutput("hist_em")),
        tabPanel("Parameters", tableOutput("parameter_table")),
        tabPanel("Input Data", tableOutput("data_preview")),
        tabPanel("AIC/BIC Plot", plotOutput("ic_plot"))
      )
    )
  )
)
