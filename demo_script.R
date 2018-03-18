library(rootSolve)
library(tidyverse)

## Define functions

get_y <- function(q_p, q_d){
  return(q_p - q_d)
}

get_alpha_star <- function(c_p, c_d){
  c <- c_p + c_d

  return( (c_p * 1.0 / c) )
}

get_S_t_of_L_point <- function(alpha, Tmax, J, c_p, c_d, k_p, k_d, t){
  c <- c_p + c_d
  k <- k_p + k_d

  S_plus = J + alpha*(c*(Tmax - t) + k) - (c_p*(Tmax - t) + k_p) + 0.0
  return(max(S_plus, 0))
}

get_S_t_of_L <- Vectorize(get_S_t_of_L_point)

S_t = function(t){get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t)}

getProbFirstArrival <- function(lambda, t_not, t){
  prob <- lambda * exp(-lambda * t) / (exp(-lambda * t_not))
  return(prob)
}

get_Expected_S_t_of_L_point <- function(lambda, alpha, J, c_p, c_d, k_p, k_d, t_endpoint, t){
  c <- c_p + c_d
  k <- k_p + k_d

  integrand <- function(x) {
    S_plus = J + alpha*(c*(Tmax - x) + k) - (c_p*(Tmax - x) + k_p)
    S = max(S_plus, 0.0)

    return(S * getProbFirstArrival(lambda, t, x))
    # return(getProbFirstArrival(lambda, t, x) * get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, x))
  }

  E_S_t <- 1.0/(1.0 - exp(-lambda * (t_endpoint - t))) * integrate(f = integrand, lower = t, upper = t_endpoint)$value

  return(E_S_t)
}

get_Expected_S_t_of_L <- Vectorize(get_Expected_S_t_of_L_point)

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

get_s_star <- function(lambda, y, c_p, c_d){
  c <- c_p + c_d
  sstar <- (c / lambda) * (1.0 / y)
  return(sstar)
}

get_t_star_solve <- function(y, alpha, Tmax, J, c_p, c_d, k_p, k_d){
  s_equation <- function(t){
    eq <- get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t) - get_s_star(lambda, y, c_p, c_d)

    return(eq)
  }

  unique_root <- uniroot(s_equation, c(0, Tmax),extendInt = "yes")$root
  return(unique_root)
}

get_t_star_analytic <- function(y, alpha, Tmax, J, c_p, c_d, k_p, k_d){
  c <- c_p + c_d
  k <- k_p + k_d
  alpha_star <- get_alpha_star(c_p,c_d)

  tstar <- Tmax - (c/(lambda*y) - (J + alpha*k - k_p) )/((alpha - alpha_star)*c)

  return(tstar)
}

get_power_case <- function(alpha, alpha_star){
  if (alpha < alpha_star) {
    return('Defendent')
  } else{
    return('Plaintiff')
  }
}

## Define model parameters
Tmax = 400 # length of negotiation period
q_p = .6 # plaintiff prob winning
q_d = .5 # defendant prob losing
lambda = 0.1 # poisson rate

c_p = 65
c_d = 50
k_p = 300
k_d = 250
J = 5000
alpha = 0.4

t_range = seq(from = 0, t = Tmax, length = 100)

y = get_y(q_p, q_d)
alpha_star = get_alpha_star(c_p,c_d)
power = get_power_case(alpha,alpha_star)
S_t = function(t){get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t)}
E_S_t = function(t){get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t)}

curve(S_t, from = 0, to = Tmax)
curve(E_S_t, from = 0, to = Tmax)


t_star_an <- get_t_star_analytic(y, alpha, Tmax, J, c_p, c_d, k_p, k_d)
t_star_comp <- get_t_star_solve(y, alpha, Tmax, J, c_p, c_d, k_p, k_d)
s_star_t = function(t){ get_s_star(lambda, y, c_p, c_d) }

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

ggplot(data.frame(x = c(0,Tmax)), aes(x = x)) +
  stat_function(fun = S_t, color = "black") +
  stat_function(fun = s_star_t, color = "blue", linetype = 2) +
  stat_function(fun = E_S_t, color = "red") +
  geom_vline(xintercept = t_star_comp, linetype = 3) +
  stat_function(fun=agreement, geom="area", aes(fill = "Agree", alpha=0.2)) +
  stat_function(fun=disagreement, geom="area", aes(fill = "Disagree", alpha=0.2)) +
  scale_fill_manual(values = c("#84CA72","grey"),
                    name="Settlement Regime",
                    breaks=c("Agree", "Disagree"),
                    labels=c("Agreement", "Disagreement")) +
  scale_alpha_continuous(guide = F) +
  theme_hc()

## Strong defendent

