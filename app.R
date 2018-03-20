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
           sliderInput("lamInput", "Lambda", 0.0000001, 1, 0.2, pre = ""),
           sliderInput("alphaInput", "Alpha", 0, 1, 0.6, pre = "")
    ),
    column(2,
           sliderInput("cpInput", "Defendant's Negotiation Fee", 0, 100, 30, pre = "$"),
           sliderInput("cdInput", "Plaintiff's Negotiation Fee", 0, 100, 50, pre = "$")
    ),
    column(2,
           sliderInput("kdInput", "Defendant's Court Fee", 0, 500, 250, pre = "$"),
           sliderInput("kpInput", "Plaintiff's Court Fee", 0, 500, 300, pre = "$")
    ),
    column(4,
           sliderInput("JInput", "Verdict Fee", 0, 5000, 2000, pre = "$"),
           numericInput("TInput","T", 100, min = 0, max = 20000)
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
        S = pmax(S_plus, 0.0)

        return(S * getConditionalProbFirstArrival(lambda, t, x))
      }

      E_S_t <- 1.0/(1.0 - exp(-lambda * (t_endpoint - t))) * integrate(f = integrand, lower = t, upper = t_endpoint, stop.on.error=FALSE)$value

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
      }

      unique_root <- uniroot(s_equation, c(0, Tmax-0.01),extendInt = "yes")$root
      return(unique_root)
    }

    get_t_star_star_solve <- function(y, alpha, Tmax, J, c_p, c_d, k_p, k_d){
      print(paste0("E[S_t(Tmax)] = ", get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, (Tmax-0.01))))
      print(paste0("S_star = ", get_s_star(lambda, y, c_p, c_d)))

      if(get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, (Tmax-0.01)) < get_s_star(lambda, y, c_p, c_d)){
        return(Tmax)
      }
      else{
        s_equation <- function(t){
          eq <- get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t) - get_s_star(lambda, y, c_p, c_d)

          return(eq)
        }

        unique_root <- uniroot(s_equation, c(0, (Tmax-0.01)),extendInt = "yes")$root
        return(unique_root)
      }
    }

    get_t_star_star_star_solve <- function(y, alpha, Tmax, J, c_p, c_d, k_p, k_d, lambda){
      c <- c_p + c_d
      k <- k_p + k_d
      alpha_star <- get_alpha_star(c_p,c_d)
      s_star <- get_s_star(lambda, y, c_p, c_d)

      get_s_star_star <- function(t){
        exp_term = exp(-lambda *(Tmax - t))
        s_star_ext <- s_star - (exp_term/(1-exp_term)) * (J - k/y)
        return(s_star_ext)
      }

      # if(get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, Tmax-0.01) > get_s_star_star(Tmax-0.01)){
      #   return(Tmax)
      # }
      # else{
        s_equation <- function(t){
          eq <- get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t) - get_s_star_star(t)

          return(eq)
        }

        unique_root <- uniroot(s_equation, c(0, Tmax-0.01),extendInt = "yes")$root
        return(unique_root)
      # }
    }

    get_s_star <- function(lambda, y, c_p, c_d){
      c <- c_p + c_d
      sstar <- (c / lambda) * (1.0 / y)
      return(sstar)
    }

    get_y_hat <- function(J, k_p, k_d){
      k <- k_p + k_d
      return(k/J)
    }

    get_y_star <- function(lambda, alpha, c_p, c_d, k_p, k_d){
      c <- c_p + c_d
      k <- k_p + k_d
      ystar <- (c / lambda) * (1.0 / (J + alpha*k - k_p))
      return(ystar)
    }

    y_star <- get_y_star(lambda, alpha, c_p, c_d, k_p, k_d)
    y_hat <- get_y_hat(J, k_p, k_d)

    t_star_analytic <- get_t_star_analytic(y, alpha, Tmax, J, c_p, c_d, k_p, k_d)

    s_star_t = function(t){ get_s_star(lambda, y, c_p, c_d) }
    s_star = get_s_star(lambda, y, c_p, c_d)

    max_y = pmax(s_star, get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, 0), get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, Tmax))

    deadline_color <- case_when(y <= y_hat ~ "#84CA72",
                                y > y_hat ~  "grey")

    print(paste0("y_hat is ", y_hat, " and y_star is ", y_star))

    ## check for strong plaintiff or defendent
    alpha_star <- get_alpha_star(c_p,c_d)
    if(alpha > alpha_star ){
      # print("hi")
      if(y > y_hat & y <= y_star){
        print("y_hat < y < y_star!")
        S_t = function(t){get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t)}
        E_S_t = function(t){get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t)}
        t_star_star_star <- get_t_star_star_star_solve(y, alpha, Tmax, J, c_p, c_d, k_p, k_d, lambda)

        print(paste0("t_star is ", t_star_analytic, " and t_star_star_star is ", t_star_star_star))

        ## TOCHECK
        agreement <- function(t) {
          y <- S_t(t)
          y[t > t_star_star_star] <- NA
          y[t < t_star_analytic ] <- NA

          return(y)
        }

        disagreement <- function(t) {
          y <- S_t(t)
          y[t > t_star_analytic & t < t_star_star_star] <- NA
          return(y)
        }

        E_S_t = function(t){get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t)}


        print("Strong plaintiff! y_hat < y_star")
        ggplot(data.frame(x = c(0,Tmax)), aes(x = x)) +
          stat_function(fun = S_t, color = "black") +
          annotate(geom="text", label=" -- Expected Settlement Amount: S_t(L)", x=(Tmax*(0.25)), y=(max_y*1.8), color = "black") +
          stat_function(fun = s_star_t, color = "blue", linetype = 2) +
          geom_vline(xintercept = t_star_analytic, linetype = 3) +
          annotate(geom="text", label="t*", x=t_star_analytic, y=-200, hjust=-0.5) +
          geom_vline(xintercept = t_star_star_star, linetype = 3) +
          annotate(geom="text", label="t**", x=t_star_star_star, y=-200, hjust=-0.5) +
          geom_segment(x = Tmax, xend=Tmax, y=0, yend=S_t(Tmax), color = deadline_color, size=4) +
          stat_function(fun = E_S_t, color = "red") +
          annotate(geom="text", label="-- E[S_t(L) | arrival in (t, T)]", x=(Tmax*(0.25)), y=(max_y*1.6), color = "red") +
          stat_function(fun=agreement, geom="area", aes(fill = "Agree", alpha=0.2)) +
          stat_function(fun=disagreement, geom="area", aes(fill = "Disagree", alpha=0.2)) +
          scale_fill_manual(values = c("#84CA72","grey"),
                            name="Settlement Regime",
                            breaks=c("Agree", "Disagree"),
                            labels=c("Agreement", "Disagreement")) +
          scale_alpha_continuous(guide = F) +
          theme_hc()+
          labs(
            x = "Time",
            y = "",
            title = "Strong Plaintiff"
          ) +
          scale_y_continuous(breaks = c(s_star), labels = c("s*")) +
          xlim(0, Tmax)
      }
      else{
        print("ystar < yhat")
        agreement <- function(t) {
          y <- S_t(t)
          y[t <= t_star_analytic] <- NA
          return(y)
        }

        disagreement <- function(t) {
          y <- S_t(t)
          y[t > t_star_analytic] <- NA
          return(y)
        }

        print("Strong plaintiff!")
        ggplot(data.frame(x = c(0,Tmax)), aes(x = x)) +
          stat_function(fun = S_t, color = "black") +
          annotate(geom="text", label=" -- Expected Settlement Amount: S_t(L)", x=(Tmax*(0.25)), y=(max_y*1.8), color = "black") +
          stat_function(fun = s_star_t, color = "blue", linetype = 2) +
          geom_vline(xintercept = t_star_analytic, linetype = 3) +
          annotate(geom="text", label="t*", x=t_star_analytic, y=-200, hjust=-0.5) +
          geom_segment(x = Tmax, xend=Tmax, y=0, yend=S_t(Tmax), color = deadline_color, size=4) +
          stat_function(fun=agreement, geom="area", aes(fill = "Agree", alpha=0.2)) +
          stat_function(fun=disagreement, geom="area", aes(fill = "Disagree", alpha=0.2)) +
          scale_fill_manual(values = c("#84CA72","grey"),
                            name="Settlement Regime",
                            breaks=c("Agree", "Disagree"),
                            labels=c("Agreement", "Disagreement")) +
          scale_alpha_continuous(guide = F) +
          theme_hc()+
          labs(
            x = "Time",
            y = "",
            title = "Strong Plaintiff"
          ) +
          scale_y_continuous(breaks = c(s_star), labels = c("s*")) +
          xlim(0, Tmax)
      }
    }
    else if(alpha < alpha_star ){
      print("Strong defendent!")

      S_t = function(t){get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t)}
      E_S_t = function(t){get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t)}

      t_star_star_comp <- get_t_star_star_solve(y, alpha, Tmax, J, c_p, c_d, k_p, k_d)
      t_star_star_star_comp <- get_t_star_star_star_solve(y, alpha, Tmax, J, c_p, c_d, k_p, k_d, lambda)
      print(paste0("t_star_star is ", t_star_star_comp, " and t_star_star_star is ", t_star_star_star_comp))

      if(y < y_hat){
        t_thresh <- t_star_star_star_comp
        t_thresh_label <- "t***"
      }
      else{
        t_thresh <- t_star_star_comp
        t_thresh_label <- "t**"
      }

      agreement <- function(t) {
        y <- S_t(t)
        y[t <= t_thresh] <- NA
        return(y)
      }

      disagreement <- function(t) {
        y <- S_t(t)
        y[t > t_thresh] <- NA
        return(y)
      }


      ggplot(data.frame(x = c(0,Tmax)), aes(x = x)) +
        stat_function(fun = S_t, color = "black") +
        annotate(geom="text", label=" -- Expected Settlement Amount: S_t(L)", x=(Tmax*(0.75)), y=(max_y*1.8), color = "black") +
        stat_function(fun = s_star_t, color = "blue", linetype = 2) +
        stat_function(fun = E_S_t, color = "red") +
        annotate(geom="text", label="-- E[S_t(L) | arrival in (t, T)]", x=(Tmax*(0.75)), y=(max_y*1.6), color = "red") +
        geom_vline(xintercept = t_thresh, linetype = 3) +
        annotate(geom="text", label=t_thresh_label, x=t_thresh, y=-200, hjust=-0.5) +
        # geom_vline(xintercept = t_star_star_star_comp, linetype = 3) +
        geom_segment(x = Tmax, xend=Tmax, y=0, yend=S_t(Tmax), color = deadline_color, size=4) +
        stat_function(fun=disagreement, geom="area", aes(fill = "Agree", alpha=0.2)) +
        stat_function(fun=agreement, geom="area", aes(fill = "Disagree", alpha=0.2)) +
        scale_fill_manual(values = c("#84CA72","grey"),
                          name="Settlement Regime",
                          breaks=c("Agree", "Disagree"),
                          labels=c("Agreement", "Disagreement")) +
        scale_alpha_continuous(guide = F) +
        theme_hc() +
        labs(
          x = "Time",
          y = "",
          title = "Strong Defendant"
        ) +
        scale_y_continuous(breaks = c(s_star), labels = c("s*")) +
        xlim(0, Tmax)
    }
    else{
      print("alpha == alpha_star")

      agreement <- function(t) {
        y <- S_t(t)
        y[t <= t_star_analytic] <- NA
        return(y)
      }

      disagreement <- function(t) {
        y <- S_t(t)
        y[t > t_star_analytic] <- NA
        return(y)
      }


      ggplot(data.frame(x = c(0,Tmax)), aes(x = x)) +
        stat_function(fun = S_t, color = "black") +
        annotate(geom="text", label="Expected Settlement Amount: S_t(L)", x=(Tmax*(0.45)), y= (S_t(Tmax*(0.45))*(1)), hjust=-0.5, color = "black") +
        # stat_function(fun = s_star_t, color = "blue", linetype = 2) +
        # geom_vline(xintercept = t_star, linetype = 3) +
        geom_segment(x = Tmax, xend=Tmax, y=0, yend=S_t(Tmax), color = deadline_color, size=4) +
        stat_function(fun=agreement, geom="area", aes(fill = "Agree", alpha=0.2)) +
        stat_function(fun=disagreement, geom="area", aes(fill = "Disagree", alpha=0.2)) +
        scale_fill_manual(values = c("#84CA72","grey"),
                          name="Settlement Regime",
                          breaks=c("Agree", "Disagree"),
                          labels=c("Agreement", "Disagreement")) +
        scale_alpha_continuous(guide = F) +
        theme_hc()+
        labs(
          x = "Time",
          y = "",
          title = "Neutral"
        ) +
        xlim(0, Tmax)
      }

  })
}

shinyApp(ui = ui, server = server)
