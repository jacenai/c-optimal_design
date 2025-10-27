# app.R
# Shiny app: c-optimal designs for the bivariate Emax model (arbitrary n-point design)
# Inputs: Smax/Emax, k2/k1, SD50/ED50, rho, sigma2^2/sigma1^2, dose bounds, number of design points
# Requires: shiny, ggplot2, pso, numDeriv

library(shiny)
library(ggplot2)
library(pso)
library(numDeriv)
library(bslib)

options(scipen = 6, digits = 6)




# Model primitives
Sigma_inv <- function(s1, s2, rho){
  det <- (s1^2)*(s2^2)*(1 - rho^2)
  matrix(c( s2^2, -rho*s1*s2,
            -rho*s1*s2,  s1^2), 2, 2, byrow=TRUE) / det
}

J_x <- function(x, ED50, Emax, SD50, Smax){
  x_ed <- x + ED50
  x_sd <- x + SD50
  if (x_ed <= 0 || x_sd <= 0) return(matrix(0, 2, 4))
  dmu1_dED50 <- -Emax * x / (x_ed^2)
  dmu1_dEmax <-  x / x_ed
  dmu1_dSD50 <-  0
  dmu1_dSmax <-  0
  
  dmu2_dED50 <-  0
  dmu2_dEmax <-  0
  dmu2_dSD50 <- -Smax * x / (x_sd^2)
  dmu2_dSmax <-  x / x_sd
  
  matrix(c(dmu1_dED50, dmu1_dEmax, dmu1_dSD50, dmu1_dSmax,
           dmu2_dED50, dmu2_dEmax, dmu2_dSD50, dmu2_dSmax),
         nrow=2, byrow=TRUE)
}

M_one_point <- function(x, pars){
  Sinv <- Sigma_inv(pars$sigma1, pars$sigma2, pars$rho)
  J <- J_x(x, pars$ED50, pars$Emax, pars$SD50, pars$Smax)
  t(J) %*% Sinv %*% J
}

M_design <- function(xs, ws, pars, ridge=1e-12){
  M <- matrix(0, 4, 4)
  for (i in seq_along(xs)) if (ws[i] > 0) M <- M + ws[i] * M_one_point(xs[i], pars)
  M + ridge * diag(4)
}

# g(theta) from Eq. (4)
g_theta <- function(theta, k1, k2){
  ED50 <- theta[1]; Emax <- theta[2]; SD50 <- theta[3]; Smax <- theta[4]
  num_sqrt <- sqrt(k1*ED50*Emax * k2*SD50*Smax)
  num <- num_sqrt*(ED50 - SD50) - ED50*SD50*(k1*Emax - k2*Smax)
  den <- (k1*ED50*Emax - k2*SD50*Smax)
  num / den
}

grad_g <- function(pars){
  theta <- c(pars$ED50, pars$Emax, pars$SD50, pars$Smax)
  as.numeric(numDeriv::grad(function(th) g_theta(th, pars$k1, pars$k2), theta))
}

psi_of_design <- function(xs, ws, pars){
  M <- M_design(xs, ws, pars)
  grad <- grad_g(pars)
  Minv <- tryCatch(solve(M), error=function(e) NULL)
  if (is.null(Minv)) return(1e50)
  as.numeric(t(grad) %*% Minv %*% grad)
}

GET_check <- function(xs, ws, pars, grid){
  M <- M_design(xs, ws, pars)
  Minv <- solve(M)
  grad <- grad_g(pars)
  rhs <- as.numeric(t(grad) %*% Minv %*% grad)
  vals <- vapply(grid, function(x){
    Mx <- M_one_point(x, pars)
    lhs <- as.numeric(t(grad) %*% Minv %*% Mx %*% Minv %*% grad)
    lhs - rhs
  }, numeric(1))
  at_design <- vapply(xs, function(x){
    Mx <- M_one_point(x, pars)
    lhs <- as.numeric(t(grad) %*% Minv %*% Mx %*% Minv %*% grad)
    lhs - rhs
  }, numeric(1))
  list(max_violation = max(vals), grid_x = grid, phi = vals, at_design = at_design)
}


# Optimization 
clean_weights <- function(w){
  w[ w < 0 & w > -1e-10 ] <- 0
  w <- pmin(pmax(w, 0), 1)
  s <- sum(w)
  if (!is.finite(s) || s <= 0) rep(1/length(w), length(w))
  else if (abs(s - 1) > 1e-10) w / s else w
}

distinct_penalty <- function(xs, min_gap=1e-3, strength=1e6){
  xs <- sort(xs)
  if (length(xs) < 2) return(0)
  gmin <- min(diff(xs))
  if (gmin >= min_gap) 0 else strength * (min_gap - gmin)^2
}

# Generic n-point objective:
# v = [x1..xn, w1..w_{n-1}], with wn := 1 - sum_{i=1}^{n-1} wi
obj_n_point <- function(v, pars, xmin, xmax, n){
  xs <- v[seq_len(n)]
  ws_head <- v[(n+1):(n + (n-1))]
  w_last <- 1 - sum(ws_head)
  ws <- c(ws_head, w_last)
  
  # if weight sum > 1 or last negative by a lot, penalize (soft constraint)
  penalty <- 0
  if (w_last < -1e-6) penalty <- penalty + 1e6 * (abs(w_last))^2
  
  # sort doses and reorder weights accordingly
  ord <- order(xs)
  xs <- xs[ord]; ws <- ws[ord]
  
  ws <- clean_weights(ws)
  if (any(xs < xmin) || any(xs > xmax)) return(1e9 + penalty)
  penalty <- penalty + distinct_penalty(xs)
  
  psi_of_design(xs, ws, pars) + penalty
}

# Local polish via L-BFGS-B for n-point
polish_n <- function(sol, pars, xmin, xmax, n){
  x0 <- c(sol$par[seq_len(n + (n-1))])
  fn <- function(v) obj_n_point(v, pars, xmin, xmax, n)
  lower <- c(rep(xmin, n), rep(0, n-1))
  upper <- c(rep(xmax, n), rep(1, n-1))
  out <- tryCatch(optim(x0, fn, method="L-BFGS-B", lower=lower, upper=upper),
                  error=function(e) NULL)
  if (is.null(out) || !is.list(out)) return(sol)
  if (out$value <= sol$value) list(par=out$par, value=out$value) else sol
}

# Helper: random Dirichlet weights (length n), return first n-1 variables used by optimizer
rand_w_head <- function(n){
  w <- rgamma(n, shape=1, rate=1); w <- w / sum(w)
  w[1:(n-1)]
}

# Decode solution vector to canonical (sorted) xs/ws
decode_solution <- function(par_vec, xmin, xmax, n){
  xs <- par_vec[seq_len(n)]
  ws_head <- par_vec[(n+1):(n + (n-1))]
  w_last <- 1 - sum(ws_head)
  ws <- c(ws_head, w_last)
  ord <- order(xs); xs <- xs[ord]; ws <- ws[ord]
  ws <- clean_weights(ws)
  list(xs=xs, ws=ws)
}


# aesthetics and plotting
theme_paper <- function(base_size = 13, base_family = "serif") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(),
      axis.text  = element_text(),
      panel.grid.major.y = element_line(size = 0.25, colour = "grey85"),
      panel.grid.minor   = element_blank()
    )
}

# placeholder_phi <- function() {
#   ggplot(data.frame(x = c(0,1), y = c(0,0)), aes(x, y)) +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     labs(
#       x = "Dose x",
#       y = expression(phi(x)),
#       title = "GET sensitivity \u03C6(x)  (\u2264 0 at optimum)"
#       # subtitle = "Tip: set parameters and click Run to compute an optimal design"
#     ) +
#     coord_cartesian(xlim = c(0,1), ylim = c(-1,1), expand = TRUE) +
#     theme_paper()
# }
# 
# placeholder_table <- function() {
#   data.frame(
#     dose       = "\u2014",
#     weight     = "\u2014",
#     phi_at_x   = "\u2014",
#     check.names = FALSE
#   )
# }
# 
# placeholder_text <- function() {
#   cat("No design yet.\n",
#       "• Set parameter ratios and rho, choose n, then click Run.\n",
#       "• The plot and table show the GET sensitivity curve and design doses.\n",
#       sep = "")
# }


# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(

  # 1) Title
  fluidRow(column(8, offset = 2, h1("Efficient c-Optimal Emax Design Using PSO"))),
  
  
  # 2) Intro/Methods section rendered from markdown with MathJax
  withMathJax(),
  div(
    class = "prose mb-4",
    includeHTML("about.html")  # write your paper-style description here
  ),
  hr(),
  
  # 3) A separator heading for the working/app section
  fluidRow(column(8, offset = 2, h2("Implementation"))),
  
  fluidRow(
    column(
      width = 8, offset = 2,   # centered, 10/12 of the page

  sidebarLayout(

    sidebarPanel(
      width = 4,
      
      
      h4("Dose bounds"),
      numericInput("xmin", "Lower bound x_min", value = 0.0, step = 0.1),
      numericInput("xmax", "Upper bound x_max", value = 500.0, min = 0.1, step = 1),
      
      h4("Design size"),
      numericInput("npoints", "Number of design points (n ≥ 2)", value = 2, min = 2, step = 1),
      
      
      h4("Ratios & correlation"),
      numericInput("ratio_SE", "Smax / Emax", value = 1.0, min = 0.05, step = 0.05),
      numericInput("ratio_k21","k2 / k1",     value = 1.0, min = 0.05, step = 0.05),
      numericInput("ratio_SD_ED","SD50 / ED50", value = 2.0, min = 0.05, step = 0.05),
      numericInput("rho", "Correlation (rho) in [-0.99, 0.99]", min=-0.99, max=0.99, value=0.0, step=0.05),
      numericInput("ratio_var", HTML("&sigma;<sub>2</sub><sup>2</sup> / &sigma;<sub>1</sub><sup>2</sup>"),
                   value = 1.0, min = 0.1, step = 0.1),

     
      h4("Optimization controls"),
      numericInput("swarm", "PSO swarm size", value = 40, min = 20, step = 10),
      numericInput("iters",  "PSO iterations", value = 50, min = 100, step = 50),
      checkboxInput("polish", "Local polish (L-BFGS-B)", value = TRUE),
      actionButton("runbtn", "Run optimization", class = "btn-primary"),

      # hr(),
      # helpText("Scale invariance: Emax=1, k1=1, ED50=1, sigma1=1."),
      # helpText("Reconstruction: Smax=Smax/Emax, k2=k2/k1, SD50=SD50/ED50, sigma2=sqrt(sigma2^2/sigma1^2).")
    ),

    mainPanel(
      h4("Criterion value"),
      p(class = "text-muted small",
        "Tip: Set dose bounds, choose n, set parameter ratios and ρ, then click Run."),
      verbatimTextOutput("resText"),
      
      h4("Design table"),
      p(class = "text-muted small",
        "Doses and weights appear here after running optimization."),
      tableOutput("atDesignTable"),
      
      h4("General Equivalence (GET) plot"),
      p(class = "text-muted small",
        "Diagnostic plot: At the optimal design, ϕ(x)≤0 across the dose range (curve should lie at or below zero), touching zero at the selected design doses."),
      plotOutput("phiPlot", height = 320),
      
      # --- Interpretation section ---
      h4("Interpretation"),
      p(
        "In practical dose-finding contexts, this design identifies doses that most efficiently characterize both benefit (efficacy) and risk (toxicity), helping to refine the dose range for Phase I, II, or III trials."
      ),
      p(
        "Researchers can adjust the input ratios and correlation to reflect different clinical trial scenarios. The resulting design suggests dose levels and corresponding weights ",
        HTML("&#123;x<sub>i</sub>, w<sub>i</sub>&#125;"), 
        " that maximize inferential efficiency in estimating the clinically optimal dose across development phases."
      ),
      p(
        "For example, if the app returns an optimal 3-point design with ",
        HTML("{x<sub>1</sub>=1, x<sub>2</sub>=9, x<sub>3</sub>=500}"), 
        " and weights ",
        HTML("{w<sub>1</sub>=0.44, w<sub>2</sub>=0.25, w<sub>3</sub>=0.31}"),
        " with GET maximum less than or equal to 0, this indicates that, under the specified model parameters, the most statistically efficient allocation is to assign approximately 44% of subjects to a low dose (~1), 25% to a mid dose (~9), and 31% to a high dose (~500)."
      )
      
    )

      )
    )
  ),
  
  # 5) Optional references at bottom (also from markdown)
  hr(),
  fluidRow(column(8, offset = 2, h2("Supplementary"))),
  div(
    class = "prose mb-4",
    includeHTML("theory.html")
  )

)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session){
  
  
  
  # observeEvent(input$rho, ignoreInit = TRUE, {
  #   val <- input$rho
  #   if (is.null(val) || !is.finite(val)) return()
  #   # snap to step and clamp to [min, max]
  #   step <- 0.05
  #   minv <- -0.99
  #   maxv <-  0.99
  #   snapped <- round(val / step) * step
  #   clamped <- max(minv, min(maxv, snapped))
  #   if (!isTRUE(all.equal(val, clamped))) {
  #     updateNumericInput(session, "rho", value = clamped)
  #     showNotification(sprintf("rho must be in [%.2f, %.2f] (step %.2f). Adjusted to %.2f.",
  #                              minv, maxv, step, clamped),
  #                      type = "message", duration = 3)
  #   }
  # })
  
  
  reconstruct_params <- reactive({
    r_SE   <- input$ratio_SE          # Smax/Emax
    r_k21  <- input$ratio_k21         # k2/k1
    r_SDED <- input$ratio_SD_ED       # SD50/ED50
    r_var  <- input$ratio_var         # sigma2^2 / sigma1^2
    rho    <- input$rho
    
    validate(
      need(is.finite(input$rho) && input$rho >= -0.99 && input$rho <= 0.99,
           "Correlation rho must be between -0.99 and 0.99."),
      need(r_SE  > 0, "Smax/Emax must be > 0"),
      need(r_k21 > 0, "k2/k1 must be > 0"),
      need(r_SDED> 0, "SD50/ED50 must be > 0"),
      need(r_var > 0, "sigma2^2/sigma1^2 must be > 0"),
      need(input$xmax > input$xmin, "x_max must be > x_min"),
      need(input$npoints >= 2, "Number of design points must be at least 2"),
    )
    
    list(
      ED50  = 1.0,
      Emax  = 1.0,
      SD50  = r_SDED,          # SD50/ED50 with ED50 fixed to 1
      Smax  = r_SE,            # Smax/Emax with Emax fixed to 1
      sigma1= 1.0,
      sigma2= sqrt(r_var),     # since r_var = sigma2^2 / sigma1^2 and sigma1=1
      rho   = rho,
      k1    = 1.0,
      k2    = r_k21            # k2/k1 with k1 fixed to 1
    )
  })
  
  run_once <- eventReactive(input$runbtn, {
    withProgress(message = "Optimizing design...", value = 0, {
      incProgress(0.05, detail = "Preparing parameters")
      pars <- reconstruct_params()
      xmin <- input$xmin; xmax <- input$xmax
      n    <- max(2L, as.integer(input$npoints))
      swarm <- input$swarm; iters <- input$iters
      polish <- isTRUE(input$polish)
      # set.seed(2025)
      
      # Build PSO problem of dimension n + (n-1)
      D <- n + (n-1)
      
      # Objective wrapper
      fn <- function(v) obj_n_point(v, pars, xmin, xmax, n)
      
      # Bounds
      lower <- c(rep(xmin, n), rep(0, n-1))
      upper <- c(rep(xmax, n), rep(1,   n-1))
      
      # Random start: doses uniform, weights ~ Dirichlet then drop last
      par0 <- c(runif(n, xmin, xmax), rand_w_head(n))
      
      incProgress(0.20, detail = sprintf("Running PSO (n = %d)", n))
      sol <- psoptim(par = par0, fn = fn, lower = lower, upper = upper,
                     control = list(maxit = iters, s = swarm, trace = 0))
      
      incProgress(0.45, detail = "Polishing")
      if (polish) sol <- polish_n(sol, pars, xmin, xmax, n)
      
      dec <- decode_solution(sol$par, xmin, xmax, n)
      xs <- dec$xs; ws <- dec$ws
      psi <- psi_of_design(xs, ws, pars)
      
      incProgress(0.25, detail = "Checking GET")
      grid <- seq(xmin, xmax, length.out = 10000)
      get <- GET_check(xs, ws, pars, grid)
      
      list(n=n, xs=xs, ws=ws, psi=psi, GET=get)
    })
  })
  
  output$resText <- renderText({
    t0  <- proc.time()
    res <- run_once()
    elapsed <- (proc.time() - t0)[["elapsed"]]
    if (is.null(res)) return(placeholder_text())
    paste0(
      sprintf("%d-point design\n", res$n),
      sprintf("\nCriterion  \u03A8 = %.8f\n", res$psi),
      sprintf("GET max violation: %.6e\n", res$GET$max_violation),
      sprintf("Runtime: %.3f s\n", elapsed)
    )
  })
  
  output$atDesignTable <- renderTable({
    res <- run_once()
    if (is.null(res)) return(placeholder_table())  # <-- add return()
    data.frame(
      dose     = res$xs,
      weight   = res$ws,
      phi_at_x = as.numeric(res$GET$at_design)
    )
  }, digits = 6)
  
  
  output$phiPlot <- renderPlot({
    res <- run_once()
    if (is.null(res)) return(placeholder_phi())
    
    df     <- data.frame(x = res$GET$grid_x, phi = res$GET$phi)
    df_all <- data.frame(x = res$xs)
    keep   <- which(res$ws > 0.01)
    
    p <- ggplot(df, aes(x, phi)) + 
      geom_line() + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      labs( 
        x = "Dose x", 
        y = expression(phi(x)), 
        title = element_blank() 
            # "GET sensitivity \u03C6(x) (\u2264 0 at optimum)" 
      ) +
      theme_paper() + 
      # theme_classic() + 
      geom_vline(data = df_all, aes(xintercept = x), 
                 linetype = "dotted", color = "grey60")
    
    
    if (length(keep) > 0) {
      i_max <- which.max(df$phi)
      df_v <- data.frame(x = res$xs[keep], key = "Dose weight > 0.01")
      p <- p +
        geom_vline(data = df_v,
                   aes(xintercept = x, linetype = key, color = key),
                   linewidth = 0.6) +
        scale_color_manual(name = "", values = c("Dose weight > 0.01" = "red")) +
        scale_linetype_manual(name = "", values = c("Dose weight > 0.01" = "dotted")) +
        # move the legend of the vertical lines to the bottom
        theme(legend.position = "bottom") +
        # guides(color = guide_legend(override.aes = list(linewidth = 1.2))) +
        # report the GET max violation
        # annotate("text", x = Inf, y = Inf,
        #          label = sprintf("GET max violation: %.2e", res$GET$max_violation),
        #          hjust = 3, vjust = 1, size = 4)
        annotate("label",
                 x = Inf, y = Inf,
                 label = sprintf("GET max violation: %.2e", res$GET$max_violation),
                 hjust = 1.02, vjust = 0.45, size = 4.3,
                 label.size = 0, fill = "white", alpha = 0) +
        coord_cartesian(clip = "off") +
        theme(plot.margin = margin(6, 20, 6, 6))   # add right margin so it’s visible
    }
    
    p
  })
  
  
}

shinyApp(ui, server)
