# Shiny app: Demonstrating Maximum Likelihood Estimation for multiple univariate distributions
# Saves as app.R and run with: shiny::runApp('path/to/this/file') or open in RStudio and click Run App

library(shiny)
library(ggplot2)
library(numDeriv) # for numeric Hessian

# Helper: negative log-likelihoods
nll_normal <- function(params, x){
  mu <- params[1]
  sigma <- params[2]
  if(sigma <= 0) return(1e10)
  -sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}

nll_exp <- function(params, x){
  rate <- params[1]
  if(rate <= 0) return(1e10)
  -sum(dexp(x, rate = rate, log = TRUE))
}

nll_pois <- function(params, x){
  lambda <- params[1]
  if(lambda <= 0) return(1e10)
  -sum(dpois(x, lambda = lambda, log = TRUE))
}

nll_binom <- function(params, x, size){
  p <- params[1]
  if(p <= 0 || p >= 1) return(1e10)
  -sum(dbinom(x, size = size, prob = p, log = TRUE))
}

nll_gamma <- function(params, x){
  shape <- params[1]
  scale <- params[2]
  if(shape <= 0 || scale <= 0) return(1e10)
  -sum(dgamma(x, shape = shape, scale = scale, log = TRUE))
}

nll_beta <- function(params, x){
  a <- params[1]
  b <- params[2]
  if(a <= 0 || b <= 0) return(1e10)
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}

ui <- fluidPage(
  titlePanel("MLE explorer: 6 univariate distributions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Distribution:",
                  choices = c("Normal", "Exponential", "Poisson", "Binomial", "Gamma", "Beta")),
      numericInput("n", "Sample size:", value = 100, min = 5, step = 1),
      numericInput("seed", "Random seed (0 = random):", value = 0, step = 1),
      hr(),
      # True params controls (conditional UI in server)
      uiOutput("params_ui"),
      actionButton("gen", "Generate sample / compute MLE"),
      hr(),
      checkboxInput("show_ll", "Show log-likelihood curve(s)", value = TRUE),
      checkboxInput("show_contour", "Show 2D log-likelihood contour (for 2-parameter models)", value = TRUE),
      hr(),
      downloadButton("downloadData", "Download sample (CSV)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("mainPlot", height = "480px")),
        tabPanel("Estimates", verbatimTextOutput("estimates")),
        tabPanel("Log-likelihood", plotOutput("llPlot", height = "480px")),
        tabPanel("Notes",
                 h4("How to use"),
                 tags$ul(
                   tags$li("Select distribution and set sample size and true parameters (for simulation)."),
                   tags$li("Click 'Generate sample / compute MLE' to simulate and compute estimates."),
                   tags$li("For 1-parameter models you'll see a log-likelihood curve; for 2-parameter models a contour will be shown if selected."),
                   tags$li("Standard errors are computed using the observed information (inverse Hessian) when possible.")
                 )
        )
      )
    )
  )
)

server <- function(input, output, session){
  
  # Dynamic UI for true params depending on distribution
  output$params_ui <- renderUI({
    switch(input$dist,
           "Normal" = tagList(
             numericInput("true_mu", "True mu:", value = 0),
             numericInput("true_sigma", "True sigma:", value = 1, min = 1e-6)
           ),
           "Exponential" = numericInput("true_rate", "True rate (lambda):", value = 1, min = 1e-6),
           "Poisson" = numericInput("true_lambda", "True lambda:", value = 2, min = 0),
           "Binomial" = tagList(
             numericInput("true_size", "Binomial size (n):", value = 10, min = 1, step = 1),
             numericInput("true_p", "True p:", value = 0.3, min = 0, max = 1, step = 0.01)
           ),
           "Gamma" = tagList(
             numericInput("true_shape", "True shape (k):", value = 2, min = 1e-6),
             numericInput("true_scale", "True scale (theta):", value = 1, min = 1e-6)
           ),
           "Beta" = tagList(
             numericInput("true_a", "True alpha:", value = 2, min = 1e-6),
             numericInput("true_b", "True beta:", value = 5, min = 1e-6)
           )
    )
  })
  
  # Reactive storage for data and results
  rv <- reactiveValues(data = NULL, results = NULL)
  
  observeEvent(input$gen, {
    # seed
    if(input$seed != 0) set.seed(input$seed)
    
    n <- input$n
    dist <- input$dist
    x <- NULL
    res <- list()
    
    if(dist == "Normal"){
      x <- rnorm(n, mean = input$true_mu, sd = input$true_sigma)
      # analytic MLEs
      mu_hat <- mean(x)
      sigma_hat <- sqrt(mean((x - mu_hat)^2))
      # numeric check via optim
      opt <- optim(par = c(mu_hat, sigma_hat), fn = nll_normal, x = x, hessian = TRUE)
      hess <- tryCatch(numDeriv::hessian(func = function(p) nll_normal(p, x), x = opt$par), error = function(e) NULL)
      se <- NA
      if(!is.null(hess) && all(is.finite(hess))){
        cov <- tryCatch(solve(hess), error = function(e) NULL)
        if(!is.null(cov)) se <- sqrt(diag(cov))
      }
      res$est <- c(mu = opt$par[1], sigma = opt$par[2])
      res$se <- se
      res$opt <- opt
    }
    
    if(dist == "Exponential"){
      x <- rexp(n, rate = input$true_rate)
      # analytic MLE: rate = 1/mean(x)
      rate_hat <- 1/mean(x)
      opt <- optim(par = rate_hat, fn = nll_exp, x = x, hessian = TRUE)
      hess <- tryCatch(numDeriv::hessian(func = function(p) nll_exp(p, x), x = opt$par), error = function(e) NULL)
      se <- NA
      if(!is.null(hess) && is.finite(hess)){
        cov <- tryCatch(1 / hess, error = function(e) NULL) # since 1x1
        if(!is.null(cov)) se <- sqrt(cov)
      }
      res$est <- c(rate = opt$par[1])
      res$se <- se
      res$opt <- opt
    }
    
    if(dist == "Poisson"){
      lambda <- input$true_lambda
      x <- rpois(n, lambda = lambda)
      lambda_hat <- mean(x)
      opt <- optim(par = lambda_hat, fn = nll_pois, x = x, hessian = TRUE)
      hess <- tryCatch(numDeriv::hessian(func = function(p) nll_pois(p, x), x = opt$par), error = function(e) NULL)
      se <- NA
      if(!is.null(hess) && is.finite(hess)){
        cov <- tryCatch(1 / hess, error = function(e) NULL)
        if(!is.null(cov)) se <- sqrt(cov)
      }
      res$est <- c(lambda = opt$par[1])
      res$se <- se
      res$opt <- opt
    }
    
    if(dist == "Binomial"){
      size <- input$true_size
      p <- input$true_p
      x <- rbinom(n, size = size, prob = p)
      p_hat <- mean(x) / size
      opt <- optim(par = p_hat, fn = function(par) nll_binom(par, x = x, size = size), hessian = TRUE, method = "L-BFGS-B", lower = 1e-8, upper = 1-1e-8)
      hess <- tryCatch(numDeriv::hessian(func = function(par) nll_binom(par, x = x, size = size), x = opt$par), error = function(e) NULL)
      se <- NA
      if(!is.null(hess) && is.finite(hess)){
        cov <- tryCatch(1 / hess, error = function(e) NULL)
        if(!is.null(cov)) se <- sqrt(cov)
      }
      res$est <- c(p = opt$par[1])
      res$se <- se
      res$opt <- opt
    }
    
    if(dist == "Gamma"){
      shape <- input$true_shape
      scale <- input$true_scale
      x <- rgamma(n, shape = shape, scale = scale)
      # use method of moments to initialize
      m <- mean(x); v <- var(x)
      shape_init <- m^2 / v
      scale_init <- v / m
      opt <- optim(par = c(shape_init, scale_init), fn = nll_gamma, x = x, hessian = TRUE, method = "L-BFGS-B", lower = c(1e-8,1e-8))
      hess <- tryCatch(numDeriv::hessian(func = function(p) nll_gamma(p, x), x = opt$par), error = function(e) NULL)
      se <- NA
      if(!is.null(hess) && all(is.finite(hess))){
        cov <- tryCatch(solve(hess), error = function(e) NULL)
        if(!is.null(cov)) se <- sqrt(diag(cov))
      }
      res$est <- c(shape = opt$par[1], scale = opt$par[2])
      res$se <- se
      res$opt <- opt
    }
    
    if(dist == "Beta"){
      a <- input$true_a
      b <- input$true_b
      # draw beta. Need to ensure values are in (0,1)
      x <- rbeta(n, shape1 = a, shape2 = b)
      # initialize with method of moments
      m <- mean(x); v <- var(x)
      # approximate alpha/beta
      tmp <- (m*(1-m)/v - 1)
      a_init <- max(1e-3, m*tmp)
      b_init <- max(1e-3, (1-m)*tmp)
      opt <- optim(par = c(a_init, b_init), fn = nll_beta, x = x, hessian = TRUE, method = "L-BFGS-B", lower = c(1e-8,1e-8))
      hess <- tryCatch(numDeriv::hessian(func = function(p) nll_beta(p, x), x = opt$par), error = function(e) NULL)
      se <- NA
      if(!is.null(hess) && all(is.finite(hess))){
        cov <- tryCatch(solve(hess), error = function(e) NULL)
        if(!is.null(cov)) se <- sqrt(diag(cov))
      }
      res$est <- c(alpha = opt$par[1], beta = opt$par[2])
      res$se <- se
      res$opt <- opt
    }
    
    rv$data <- x
    rv$results <- res
  })
  
  output$mainPlot <- renderPlot({
    req(rv$data)
    x <- rv$data
    dist <- input$dist
    df <- data.frame(x = x)
    
    p <- NULL
    if(dist %in% c("Normal", "Exponential", "Gamma", "Beta")){
      # continuous: histogram + fitted density
      p <- ggplot(df, aes(x = x)) + geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5)
      est <- rv$results$est
      if(dist == "Normal"){
        p <- p + stat_function(fun = dnorm, args = list(mean = est['mu'], sd = est['sigma']), size = 1)
      }
      if(dist == "Exponential"){
        p <- p + stat_function(fun = dexp, args = list(rate = est['rate']), size = 1)
      }
      if(dist == "Gamma"){
        p <- p + stat_function(fun = dgamma, args = list(shape = est['shape'], scale = est['scale']), size = 1)
      }
      if(dist == "Beta"){
        # Beta is on (0,1) — use density if data in (0,1)
        p <- p + stat_function(fun = dbeta, args = list(shape1 = est['alpha'], shape2 = est['beta']), size = 1)
      }
      p <- p + theme_minimal() + ggtitle(paste0("Histogram and fitted density: ", dist))
    } else {
      # discrete: bar plot of counts + fitted pmf
      tab <- as.data.frame(table(x))
      names(tab) <- c('x', 'freq')
      tab$x <- as.numeric(as.character(tab$x))
      p <- ggplot(tab, aes(x = x, y = freq)) + geom_col(alpha = 0.6)
      est <- rv$results$est
      if(dist == "Poisson"){
        xs <- min(tab$x):max(tab$x)
        pmf <- dpois(xs, lambda = est['lambda']) * length(x)
        p <- p + geom_point(data = data.frame(x = xs, freq = pmf), aes(x = x, y = freq), color = 'red', size = 2)
      }
      if(dist == "Binomial"){
        size <- input$true_size
        xs <- 0:size
        pmf <- dbinom(xs, size = size, prob = est['p']) * length(x)
        p <- p + geom_point(data = data.frame(x = xs, freq = pmf), aes(x = x, y = freq), color = 'red', size = 2)
      }
      p <- p + theme_minimal() + ggtitle(paste0("Counts and fitted pmf: ", dist))
    }
    print(p)
  })
  
  output$estimates <- renderPrint({
    req(rv$results)
    cat("Maximum likelihood estimates:\n")
    print(rv$results$est)
    cat("\nStandard errors (from observed information, if available):\n")
    print(rv$results$se)
    cat("\nOptimizer output (convergence = 0 means success):\n")
    print(rv$results$opt)
  })
  
  output$llPlot <- renderPlot({
    req(rv$data, rv$results)
    x <- rv$data
    dist <- input$dist
    est <- rv$results$est
    
    if(!input$show_ll) return(NULL)
    
    # 1-parameter likelihood curve
    if(dist %in% c("Exponential", "Poisson", "Binomial")){
      # treat as single param
      if(dist == "Exponential"){
        grid <- seq(max(1e-6, est['rate']*0.1), est['rate']*3, length.out = 200)
        ll <- sapply(grid, function(r) -nll_exp(r, x))
        df <- data.frame(param = grid, loglik = ll)
        ggplot(df, aes(x = param, y = loglik)) + geom_line() + geom_vline(xintercept = est['rate'], linetype = "dashed") + theme_minimal() + xlab('rate') + ggtitle('Log-likelihood for rate')
      } else if(dist == "Poisson"){
        grid <- seq(max(1e-6, est['lambda']*0.1), max(est['lambda']*3, mean(x)+5), length.out = 200)
        ll <- sapply(grid, function(l) -nll_pois(l, x))
        df <- data.frame(param = grid, loglik = ll)
        ggplot(df, aes(x = param, y = loglik)) + geom_line() + geom_vline(xintercept = est['lambda'], linetype = "dashed") + theme_minimal() + xlab('lambda') + ggtitle('Log-likelihood for lambda')
      } else if(dist == "Binomial"){
        size <- input$true_size
        grid <- seq(1e-6, 1-1e-6, length.out = 200)
        ll <- sapply(grid, function(p) -nll_binom(p, x, size))
        df <- data.frame(param = grid, loglik = ll)
        ggplot(df, aes(x = param, y = loglik)) + geom_line() + geom_vline(xintercept = est['p'], linetype = "dashed") + theme_minimal() + xlab('p') + ggtitle('Log-likelihood for p')
      }
    }
    
    # 2-parameter: Normal, Gamma, Beta
    else if(dist == "Normal"){
      mu_grid <- seq(est['mu'] - 3*est['sigma'], est['mu'] + 3*est['sigma'], length.out = 120)
      sigma_grid <- seq(max(1e-6, est['sigma']/5), est['sigma']*3, length.out = 120)
      if(input$show_contour){
        llmat <- outer(mu_grid, sigma_grid, Vectorize(function(mu, sigma) -nll_normal(c(mu, sigma), x)))
        filled.contour(mu_grid, sigma_grid, llmat, color.palette = terrain.colors, xlab = 'mu', ylab = 'sigma', main = 'Log-likelihood (Normal)')
      } else {
        # plot profile over mu
        prof_mu <- sapply(mu_grid, function(mu) -optimize(function(s) nll_normal(c(mu, s), x), interval = c(1e-6, max(10, est['sigma']*5)))$objective)
        df <- data.frame(mu = mu_grid, loglik = prof_mu)
        ggplot(df, aes(x = mu, y = loglik)) + geom_line() + geom_vline(xintercept = est['mu'], linetype = 'dashed') + theme_minimal() + ggtitle('Profile log-likelihood for mu (Normal)')
      }
    }
    else if(dist == "Gamma"){
      shape_grid <- seq(max(1e-6, est['shape']/5), est['shape']*3, length.out = 120)
      scale_grid <- seq(max(1e-6, est['scale']/5), est['scale']*3, length.out = 120)
      if(input$show_contour){
        llmat <- outer(shape_grid, scale_grid, Vectorize(function(a, b) -nll_gamma(c(a, b), x)))
        filled.contour(shape_grid, scale_grid, llmat, color.palette = heat.colors, xlab = 'shape', ylab = 'scale', main = 'Log-likelihood (Gamma)')
      } else {
        prof_shape <- sapply(shape_grid, function(s) -optimize(function(sc) nll_gamma(c(s, sc), x), interval = c(1e-6, max(est['scale']*5, 10)))$objective)
        df <- data.frame(shape = shape_grid, loglik = prof_shape)
        ggplot(df, aes(x = shape, y = loglik)) + geom_line() + geom_vline(xintercept = est['shape'], linetype = 'dashed') + theme_minimal() + ggtitle('Profile log-likelihood for shape (Gamma)')
      }
    }
    else if(dist == "Beta"){
      a_grid <- seq(max(1e-6, est['alpha']/5), est['alpha']*3, length.out = 120)
      b_grid <- seq(max(1e-6, est['beta']/5), est['beta']*3, length.out = 120)
      if(input$show_contour){
        llmat <- outer(a_grid, b_grid, Vectorize(function(a, b) -nll_beta(c(a, b), x)))
        filled.contour(a_grid, b_grid, llmat, color.palette = topo.colors, xlab = 'alpha', ylab = 'beta', main = 'Log-likelihood (Beta)')
      } else {
        prof_a <- sapply(a_grid, function(a) -optimize(function(b) nll_beta(c(a, b), x), interval = c(1e-6, max(est['beta']*5, 10)))$objective)
        df <- data.frame(alpha = a_grid, loglik = prof_a)
        ggplot(df, aes(x = alpha, y = loglik)) + geom_line() + geom_vline(xintercept = est['alpha'], linetype = 'dashed') + theme_minimal() + ggtitle('Profile log-likelihood for alpha (Beta)')
      }
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function(){ paste0('sample_', input$dist, '.csv') },
    content = function(file){ write.csv(data.frame(x = rv$data), file, row.names = FALSE) }
  )
}

shinyApp(ui, server)
