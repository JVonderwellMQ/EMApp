library(shiny)
library(ggplot2)

em_gmm <- function(x, modes = 3, max_iter = 200, tol = 1e-6) {
  x <- as.numeric(x)
  
  # Initialize estimates
  mu <- sample(x, modes)
  #quantiles <- 1:modes / (modes + 1)
  #mu <- quantile(x, probs = quantiles)
  
  sigma <- rep(sd(x), modes)
  
  ## probability of a point belonging to each cluster.
  pi <- rep(1/modes, modes)
  
  loglik_prev <- -Inf
  
  for (i in 1:max_iter) {
    
    # tau: n × K responsibility matrix
    # row i = x_i
    # col k = mixture component k
    tau <- sapply(1:modes, function(k)
      pi[k] * dnorm(x, mean = mu[k], sd = sigma[k])
    )
    
    # gamma: n x K matrix
    # gamma[i,k] = probability that x_i belongs to component k
    gamma <- tau / rowSums(tau)
    
    # pi: K x 1 vector with weights for each component k
    pi <- colMeans(gamma)
    
    mu <- colSums(gamma * x) / colSums(gamma)
    sigma <- sqrt(
      colSums(gamma * (x - matrix(mu, nrow = length(x), ncol = modes, byrow = TRUE))^2) / 
      colSums(gamma))
    
    loglik <- sum(log(rowSums(tau)))
    
    if (abs(loglik - loglik_prev) < tol) break
    loglik_prev <- loglik
  }
  
  list(
    pi = pi, 
    mu = mu, 
    sigma = sigma, 
    loglik = loglik)
}

# -------------------------------
# Shiny Server
# -------------------------------
function(input, output, session) {
  
  data <- reactive({
    validate(
      need(input$inputfile, "Please upload a file."),
      need(grepl("\\.csv$", input$inputfile$name, ignore.case = TRUE),
           "File must be a .csv")
    )
    read.csv(input$inputfile$datapath)
  })
  
  output$filename <- renderText({
    req(data())
    paste("File name:", input$inputfile$name)
  })
  
  output$preview <- renderTable({
    head(data())
  })
  
  # EM results based on first numeric column
  em_results <- reactive({
    df <- data()
    num_cols <- sapply(df, is.numeric)
    validate(need(any(num_cols), "No numeric columns found."))
    col <- df[[which(num_cols)[1]]]
    em_gmm(col)
  })
  
  output$em_results <- renderPrint({
    req(em_results())
    print(em_results())
  })
  
  # ---------------------------------------------
  # Histogram + EM mixture density curve
  # ---------------------------------------------
  output$hist_em <- renderPlot({
    req(data(), em_results())
    
    df <- data()
    
    # first numeric column
    num_cols <- sapply(df, is.numeric)
    col <- df[[which(num_cols)[1]]]
    colname <- names(df)[which(num_cols)[1]]
    
    # EM parameters
    em <- em_results()
    pi <- em$pi
    mu <- em$mu
    sigma <- em$sigma
    
    # density curve for plotting
    xgrid <- seq(min(col), max(col), length.out = 400)
    mix_density <- rowSums(mapply(function(p, m, s) p * dnorm(xgrid, m, s),
                                  pi, mu, sigma, SIMPLIFY = TRUE))
    
    ggplot(data.frame(x = col), aes(x)) +
      geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.4, fill = "blue") +
      geom_line(data = data.frame(x = xgrid, y = mix_density),
                aes(x, y), linewidth = 1.2) +
      labs(
        x = colname,
        y = "Density",
        title = "Histogram with Fitted EM Gaussian Mixture"
      ) +
      theme_minimal()
  })
}
