library(shiny)
library(ggplot2)

MAX_NUM_MODES = 10

ics <- reactiveValues(aics = NULL, bics = NULL)

em_gmm <- function(x, modes, max_iter = 1000, tol = 1e-6, min_sigma = 0.1) {
  x <- as.numeric(x)
  min_sigma <- 0.05 * sd(x)
  
  # Initialize estimates
  #mu <- sample(x, modes)
  quantiles <- sapply(1:modes, function(mode) mode / (modes + 1))
  mu <- quantile(x, probs = quantiles)
  sigma <- rep(sd(x), modes)
  pi <- rep(1/modes, modes)
  loglik_prev <- -Inf
  
  history <- vector("list", max_iter)
  
  for (i in 1:max_iter) {
    tau <- sapply(1:modes, function(k)
      pi[k] * dnorm(x, mean = mu[k], sd = sigma[k])
    )
    
    gamma <- tau / rowSums(tau)
    
    pi <- colMeans(gamma)
    mu <- colSums(gamma * x) / colSums(gamma)
    
    # Standard EM update with minimum sigma clamp
    sigma <- sqrt(
      colSums(gamma * (x - matrix(mu, nrow=length(x), ncol=modes, byrow=TRUE))^2) / colSums(gamma)
    )
    sigma <- pmax(sigma, min_sigma)  # <- enforce minimum sigma
    
    loglik <- sum(log(rowSums(tau)))
    
    history[[i]] <- list(mu=mu, sigma=sigma, pi=pi, loglik=loglik)
    
    if (abs(loglik - loglik_prev) < tol) break
    loglik_prev <- loglik
  }
  
  history = history[1:i]
  
  list(
    num_modes = modes,
    pi = pi, 
    mu = mu, 
    sigma = sigma, 
    loglik = loglik,
    history = history
  )
}

em_aic <- function(emResults) {
  # mu + sigma + (pi - 1)
  k <- 3 * emResults$num_modes - 1
  
  2 * k - 2 * emResults$loglik
}

em_bic <- function(emResults, n) {
  # mu + sigma + (pi - 1)
  k <- 3 * emResults$num_modes - 1
  
  log(n) * k - 2 * emResults$loglik
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
    req(em_results())
    history <- em_results()$history
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

  # EM results based on first numeric column
  em_results <- reactive({
    df <- data()
    num_cols <- sapply(df, is.numeric)
    validate(need(any(num_cols), "No numeric columns found."))
    col <- df[[which(num_cols)[1]]]
    
    all_modes_results = lapply(1:MAX_NUM_MODES, function(n) em_gmm(col, modes = n))
    
    ics$aics <- sapply(all_modes_results, function(result) em_aic(result))
    ics$bics <- sapply(all_modes_results, function(result) em_bic(result, n = length(col)))
    
    req(input$modeSelection)
    mode_selection <- as.numeric(input$modeSelection)
  
    if (mode_selection == 1) {
      num_modes <- which.min(ics$aics)
    } else if (mode_selection == 2) {
      num_modes <- which.min(ics$bics)
    } else if (mode_selection == 3) {
      req(!is.null(input$numModes), !is.na(input$numModes))
      num_modes <- input$numModes
    }
    
    all_modes_results[[num_modes]]
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
  
  # ---------------------------------------------
  # AIC plot
  # ---------------------------------------------
  output$ic_plot <- renderPlot({
    req(ics$aics)
    req(length(ics$aics) > 0)
    req(ics$bics)
    req(length(ics$bics) > 0)
    
    df <- data.frame(
      modes = 1:length(ics$aics),
      AIC = ics$aics,
      BIC = ics$bics
    )
    
    # Find minimums
    aic_min_index <- which.min(ics$aics)
    aic_min_value <- ics$aics[aic_min_index]
    bic_min_index <- which.min(ics$bics)
    bic_min_value <- ics$bics[bic_min_index]
    
    ggplot(df, aes(x = factor(modes))) +
      geom_line(aes(y = AIC, group = 1), color = "blue", size = 1) +
      geom_point(aes(y = AIC), color = "blue", size = 3) +
      geom_point(data = data.frame(modes = aic_min_index, IC = aic_min_value),
                 aes(x = factor(modes), y = IC),
                 color = "blue", size = 4) +
      geom_vline(xintercept = aic_min_index, linetype = "dashed", color = "blue") +
      
      geom_line(aes(y = BIC, group = 1), color = "green", size = 1) +
      geom_point(aes(y = BIC), color = "green", size = 3) +
      geom_point(data = data.frame(modes = bic_min_index, IC = bic_min_value),
                 aes(x = factor(modes), y = IC),
                 color = "green", size = 4) +
      geom_vline(xintercept = bic_min_index, linetype = "dashed", color = "green") +
      
      labs(
        x = "Number of Modes",
        y = "Information Criterion",
        title = "AIC and BIC vs Number of Modes"
      ) +
      theme_minimal()
    
  })
  
}
