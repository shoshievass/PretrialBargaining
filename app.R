library(shiny)
library(ggplot2)
library(dplyr)
library(ggthemes)

ui <- fluidPage(

  titlePanel("Settlement Dynamics", windowTitle = "Primitive Parameters"),

  plotOutput("settlementPlot"),

  hr(),


  fluidRow(
    column(4,
           sliderInput("yInput", "Optimism (y)", 0, 1, 0.15, pre = ""),
           sliderInput("lamInput", "Lambda", 0, 1, 0.1, pre = ""),
           sliderInput("alphaInput", "Alpha", 0, 1, 0.6, pre = "")
    ),
    column(2,
           sliderInput("cpInput", "Defendant's Negotiation Fee", 0, 100, 50, pre = "$"),
           sliderInput("cdInput", "Plaintiff's Negotiation Fee", 0, 100, 50, pre = "$")
    ),
    column(2,
           sliderInput("kdInput", "Defendant's Court Fee", 0, 500, 250, pre = "$"),
           sliderInput("kpInput", "Plaintiff's Court Fee", 0, 500, 300, pre = "$")
    ),
    column(4,
           sliderInput("JInput", "Verdict Fee", 0, 5000, 2000, pre = "$"),
           numericInput("TInput","T", 1000, min = 0, max = 20000)
    )
  )

)

server <- function(input, output) {

  output$settlementPlot <- renderPlot({
    Tmax = input$TInput # length of negotiation period
    y = input$yInput
    lambda = input$lamInput# poisson rate

    c_p = input$cpInput
    c_d = input$cdInput
    k_p = input$kpInput
    k_d = input$kdInput
    J = input$JInput
    alpha = input$alphaInput

    get_alpha_star <- function(c_p, c_d){
      c <- c_p + c_d

      return( (c_p * 1.0 / c) )
    }

    get_S_t_of_L_point <- function(alpha, Tmax, J, c_p, c_d, k_p, k_d, t){
      c <- c_p + c_d
      k <- k_p + k_d

      S_plus = J + alpha*(c*(Tmax - t) + k) - (c_p*(Tmax - t) + k_p)
      return(max(S_plus, 0))
    }

    get_S_t_of_L <- Vectorize(get_S_t_of_L_point)

    S_t = function(t){get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t)}

    getConditionalProbFirstArrival <- function(lambda, t_not, t){
      prob <- lambda * exp(-lambda * (t-t_not))
      return(prob)
    }

    get_Expected_S_t_of_L_point <- function(lambda, alpha, J, c_p, c_d, k_p, k_d, t_endpoint, t){
      c <- c_p + c_d
      k <- k_p + k_d

      integrand <- function(x) {
        S_plus = J + alpha*(c*(t_endpoint - x) + k) - (c_p*(t_endpoint - x) + k_p)
        S = max(S_plus, 0.0)

        return(S * getConditionalProbFirstArrival(lambda, t, x))
        # return(getProbFirstArrival(lambda, t, x) * get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, x))
      }

      E_S_t <- 1.0/(1.0 - exp(-lambda * (t_endpoint - t))) * integrate(f = integrand, lower = t, upper = t_endpoint)$value

      return(E_S_t)
    }

    get_Expected_S_t_of_L <- Vectorize(get_Expected_S_t_of_L_point)

    get_t_star_analytic <- function(y, alpha, Tmax, J, c_p, c_d, k_p, k_d){
      c <- c_p + c_d
      k <- k_p + k_d
      alpha_star <- get_alpha_star(c_p,c_d)

      tstar <- Tmax - (c/(lambda*y) - (J + alpha*k - k_p) )/((alpha - alpha_star)*c)

      return(tstar)
    }

    get_t_star_solve <- function(y, alpha, Tmax, J, c_p, c_d, k_p, k_d){
      s_equation <- function(t){
        eq <- get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t) - get_s_star(lambda, y, c_p, c_d)

        return(eq)
        return(min(eq, Tmax))
      }

      unique_root <- uniroot(s_equation, c(0, Tmax),extendInt = "yes")$root
      return(unique_root)
    }

    get_t_star_star_solve <- function(y, alpha, Tmax, J, c_p, c_d, k_p, k_d){
      if(get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, Tmax-0.01) < get_s_star(lambda, y, c_p, c_d)){
        return(Tmax)
      }
      else{
        s_equation <- function(t){
          eq <- get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t) - get_s_star(lambda, y, c_p, c_d)

          return(eq)
        }

        unique_root <- uniroot(s_equation, c(0, Tmax-0.01),extendInt = "yes")$root
        return(unique_root)
      }
    }

    get_s_star <- function(lambda, y, c_p, c_d){
      c <- c_p + c_d
      sstar <- (c / lambda) * (1.0 / y)
      return(sstar)
    }

    t_star_comp <- get_t_star_solve(y, alpha, Tmax, J, c_p, c_d, k_p, k_d)

    s_star_t = function(t){ get_s_star(lambda, y, c_p, c_d) }

    ## check for strong plaintiff or defendent
    alpha_star <- get_alpha_star(c_p,c_d)
    if(alpha > alpha_star | alpha == alpha_star){

      agreement <- function(t) {
        y <- S_t(t)
        y[t < t_star_comp] <- NA
        return(y)
      }

      disagreement <- function(t) {
        y <- S_t(t)
        y[t > t_star_comp] <- NA
        return(y)
      }

      print("Strong plaintiff!")
      ggplot(data.frame(x = c(0,Tmax)), aes(x = x)) +
        stat_function(fun = S_t, color = "black") +
        stat_function(fun = s_star_t, color = "blue", linetype = 2) +
        geom_vline(xintercept = t_star, linetype = 3) +
        stat_function(fun=agreement, geom="area", aes(fill = "Agree", alpha=0.2)) +
        stat_function(fun=disagreement, geom="area", aes(fill = "Disagree", alpha=0.2)) +
        scale_fill_manual(values = c("#84CA72","grey"),
                          name="Settlement Regime",
                          breaks=c("Agree", "Disagree"),
                          labels=c("Agreement", "Disagreement")) +
        scale_alpha_continuous(guide = F) +
        theme_hc()
    }
    else{
      print("Strong defendent!")

      S_t = function(t){get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t)}
      E_S_t = function(t){get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t)}
      t_star_star_comp <- get_t_star_star_solve(y, alpha, Tmax, J, c_p, c_d, k_p, k_d)

      agreement <- function(t) {
        y <- S_t(t)
        y[t < t_star_star_comp] <- NA
        return(y)
      }

      disagreement <- function(t) {
        y <- S_t(t)
        y[t > t_star_star_comp] <- NA
        return(y)
      }


      ggplot(data.frame(x = c(0,Tmax)), aes(x = x)) +
        stat_function(fun = S_t, color = "black") +
        stat_function(fun = s_star_t, color = "blue", linetype = 2) +
        stat_function(fun = E_S_t, color = "red") +
        # geom_vline(xintercept = t_star_star_comp, linetype = 3) +
        stat_function(fun=agreement, geom="area", aes(fill = "Agree", alpha=0.2)) +
        stat_function(fun=disagreement, geom="area", aes(fill = "Disagree", alpha=0.2)) +
        scale_fill_manual(values = c("#84CA72","grey"),
                          name="Settlement Regime",
                          breaks=c("Agree", "Disagree"),
                          labels=c("Agreement", "Disagreement")) +
        scale_alpha_continuous(guide = F) +
        theme_hc()
    }

  })
}

shinyApp(ui = ui, server = server)
