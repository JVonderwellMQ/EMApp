gmm_generation <- function(n, weights, means, sds, plot = TRUE) {
  set.seed(6986897)
  
  stopifnot(n > 0)
  stopifnot(abs(sum(weights) - 1) <= 0.0001 )
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

weights <- rep(0.1, 10)
means <- runif(10, 0, 1000)
sds <- runif(10, 5, 25)
data <- gmm_generation(1000, weights, means, sds)

filename <- paste0("gmm_", length(weights), "_components.csv")
filepath <- paste0("data/", filename)
write.csv(data, filepath, row.names = FALSE)