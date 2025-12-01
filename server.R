  library(shiny)
  library(ggplot2)
  
  MAX_NUM_MODES = 10
  MAX_NUM_ITERATIONS = 1000
  TOL = 1e-6
  
  ics <- reactiveValues(aics = NULL, bics = NULL)
  
  # Runs EM for a gaussian mixture model given data and a number of modes.
  # Outputs a list with length equal to the number of iterations run.
  # Each list contains the parameters of the model, each stored as a list of length
  #   equal to the number of modes, as well as the log likelihood.
  # Output params: pi, mu, sigma, loglik.
  em_gmm <- function(x, modes, max_iter = MAX_NUM_ITERATIONS, tol = TOL) {
    x <- as.numeric(x)
    
    quantiles <- sapply(1:modes, function(m) m/(modes+1))
    mu <- quantile(x, probs = quantiles)
    sigma <- rep(sd(x), modes)
    pi <- rep(1/modes, modes)
    
    loglik_prev <- -Inf
    history <- vector("list", max_iter)
    
    for(i in 1:max_iter){
      tau <- sapply(1:modes, function(k) pi[k] * dnorm(x, mu[k], sigma[k]))
      gamma <- tau / rowSums(tau)
      
      pi <- colMeans(gamma)
      mu <- colSums(gamma * x) / colSums(gamma)
      sigma <- sqrt(colSums(gamma*(x - matrix(mu, length(x), modes, TRUE))^2) /
                      colSums(gamma))
      
      loglik <- sum(log(rowSums(tau)))
      history[[i]] <- list(mu=mu, sigma=sigma, pi=pi, loglik=loglik)
      
      if(abs(loglik - loglik_prev) < tol) break
      loglik_prev <- loglik
    }
    history[1:i]
  }
  
  em_aic <- function(numModes, loglik) {
    # mu + sigma + (pi - 1)
    k <- 3 * numModes - 1
    
    2 * k - 2 * loglik
  }
  
  em_bic <- function(numModes, loglik, n) {
    # mu + sigma + (pi - 1)
    k <- 3 * numModes - 1
    
    log(n) * k - 2 * loglik
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
    
    column <- reactive({
      req(data())
      df <- data()
      validate(
        need(input$column > 0, "Column number must be > 0."),
        need(input$column <= ncol(df), "Column number not found in dataset.")
      )
      col <- df[[input$column]]
    })
    
    output$filename <- renderText({
      req(data())
      paste("File name:", input$inputfile$name)
    })
    
    output$parameter_table <- renderTable({
      history <- em_results()
      req(length(history) > 0)
      
      df <- do.call(rbind, lapply(seq_along(history), function(i) {
        data.frame(
          Iteration = i,
          Component = seq_along(history[[i]]$mu),
          Mu = round(history[[i]]$mu, 4),
          Sigma = round(history[[i]]$sigma, 4),
          Pi = round(history[[i]]$pi, 4),
          Loglik = round(history[[i]]$loglik, 4)
        )
      }))
      
      df
    })
    
    output$data_preview <- renderTable({
      req(data())
      head(data(), 20)   # show first 20 rows
    })
  
    # EM results based on the chosen column
    all_modes <- reactive({
      # Run EM for 1:MAX_NUM_MODES
      lapply(1:MAX_NUM_MODES,function(k) em_gmm(column() ,k))
    })
    
    num_modes <- reactive({
      all_modes_em <- all_modes()
      
      # Compute AIC/BIC from final iteration of each mode
      finals <- lapply(all_modes_em, function(h) h[[length(h)]])
      ics$aics <- sapply(finals, function(h) em_aic(length(h$mu), loglik=h$loglik))
      ics$bics <- sapply(finals, function(h) em_bic(length(h$mu), h$loglik, length(column())))
      
      # Determine chosen number of components
      switch(input$modeSelection,
             "1" = which.min(ics$aics),
             "2" = which.min(ics$bics),
             "3" = input$numModes)
    })
    
    # EM results for the chosen number of modes.
    em_results <- reactive({
      all_modes_em <- all_modes()
      num_modes_chosen <- num_modes()
      all_modes_em[[num_modes_chosen]]
    })
    
    current_iteration <- reactive({
      history <- em_results()
      req(input$iterSelect <= length(history))
      history[[input$iterSelect]]
    })
    
    observeEvent(em_results(), {
      history <- em_results()
      updateSliderInput(session, "iterSelect",
                        min = 1,
                        max = length(history),
                        value = length(history))
      updateNumericInput(session, "numModes", value = num_modes())
    })
    
    output$em_results <- renderPrint({
      req(em_results())
      print(em_results())
    })
    
    # ---------------------------------------------
    # Histogram + EM mixture density curve
    # ---------------------------------------------
    output$hist_em <- renderPlot({
      df <- data()
      col <- df[[which(sapply(df, is.numeric))[input$column]]]
      
      em <- current_iteration()
      
      xgrid <- seq(min(col), max(col), length.out = 400)
      mix <- rowSums(mapply(function(p, m, s) p * dnorm(xgrid, m, s), em$pi, em$mu, em$sigma))
      
      ggplot(data.frame(x=col), aes(x=x)) +
        geom_histogram(aes(y=..density..), bins=50, alpha=0.4, fill="blue") +
        geom_line(data=data.frame(x=xgrid, y=mix), aes(x=x, y=y), color="red", linewidth=1.2) +
        labs(title=paste("Iteration", input$iterSelect, "— Mixture Density"), y="Density") +
        theme_minimal()
    })
    
    # ---------------------------------------------
    # IC plot
    # ---------------------------------------------
    output$ic_plot <- renderPlot({
      req(ics$aics)
      req(length(ics$aics) > 0)
      req(ics$bics)
      req(length(ics$bics) > 0)
      
      modes <- 1:length(ics$aics)
      
      # Find minimums
      aic_min_index <- which.min(ics$aics)
      aic_min_value <- ics$aics[aic_min_index]
      bic_min_index <- which.min(ics$bics)
      bic_min_value <- ics$bics[bic_min_index]
      
      ggplot() +
        # AIC line and points
        geom_line(aes(x=factor(modes), y=ics$aics, color="AIC", group=1), size=1) +
        geom_point(aes(x=factor(modes), y=ics$aics, color="AIC"), size=3) +
        geom_point(aes(x=factor(aic_min_index), y=aic_min_value, color="AIC"), size=4, shape=18) +
        geom_vline(xintercept = aic_min_index, linetype = "dashed", color = "blue") +
        
        # BIC line and points
        geom_line(aes(x=factor(modes), y=ics$bics, color="BIC", group=1), size=1) +
        geom_point(aes(x=factor(modes), y=ics$bics, color="BIC"), size=3) +
        geom_point(aes(x=factor(bic_min_index), y=bic_min_value, color="BIC"), size=4, shape=18) +
        geom_vline(xintercept = bic_min_index, linetype = "dashed", color = "green") +
        
        # Legend
        scale_color_manual(name="Criterion", values=c("AIC"="blue","BIC"="green")) +
        
        labs(
          x = "Number of Modes",
          y = "Information Criterion",
          title = "AIC and BIC vs Number of Modes"
        ) +
        theme_minimal()
    })
  }
