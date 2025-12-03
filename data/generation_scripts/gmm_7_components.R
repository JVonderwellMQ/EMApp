gmm_generation <- function(n, weights, means, sds, plot = TRUE) {
  set.seed(6986897)
  
  stopifnot(n > 0)
  stopifnot(sum(weights) == 1)
  stopifnot(length(weights) == length(means))
  stopifnot(length(weights) == length(sds))
  
  n_components <- length(weights)
  
  # draw component for each sample
  components <- sample(1:n_components, n, replace = TRUE, prob = weights)
  
  # generate mixture data
  data <- rnorm(n, mean = means[components], sd = sds[components])
  
  if (plot) {
    hist(data, breaks=40, freq=FALSE, col="lightgray")
    
    mixture_density <- function(x) {
      rowSums(mapply(function(w, m, s) w * dnorm(x, m, s),
                     weights, means, sds, SIMPLIFY = TRUE))
    }
    curve(mixture_density, add = TRUE, col = "red", lwd = 2)
  }
  return(data)
}

weights <- c(0.1, 0.2, 0.2, 0.1, 0.05, 0.25, 0.1)
means <- c(10, 150, 200, 400, 550, 800, 1000)
sds <- c(10, 5, 10, 12, 8, 16, 30)
data <- gmm_generation(1000, weights, means, sds)

filename <- paste0("gmm_", length(weights), "_components.csv")
filepath <- paste0("data/", filename)
write.csv(data, filepath, row.names = FALSE)