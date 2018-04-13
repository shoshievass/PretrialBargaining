library(rootSolve)
library(tidyverse)
library(ggthemes)


## Define model parameters
Tmax = 100 # length of negotiation period
q_p = .65 # plaintiff prob winning
q_d = .5 # defendant prob losing
lambda = 0.5 # poisson rate

c_p = 50
c_d = 30
k_p = 300
k_d = 250
J = 2000
alpha = 0.3

t_range = seq(from = 0, t = Tmax, length = 100)


## Define functions

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


y = get_y(q_p, q_d)
alpha_star = get_alpha_star(c_p,c_d)
power = get_power_case(alpha,alpha_star)
S_t = function(t){get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t)}
E_S_t = function(t){get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t)}
s_star_star_t = function(t){get_s_star_star(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t)}
f = function(t){getConditionalProbFirstArrival(0.1,0,t)}

curve(S_t, from = 0, to = Tmax)
curve(E_S_t, from = 0, to = Tmax)
curve(s_star_star_t, from = 0, to = Tmax)

S_t = function(t){get_S_t_of_L(alpha, Tmax, J, c_p, c_d, k_p, k_d, t)}
E_S_t = function(t){get_Expected_S_t_of_L(lambda, alpha, J, c_p, c_d, k_p, k_d, Tmax, t)}

t_star_star_comp <- get_t_star_star_solve(y, alpha, Tmax, J, c_p, c_d, k_p, k_d)
t_star_star_star_comp <- get_t_star_star_star_solve(y, alpha, Tmax, J, c_p, c_d, k_p, k_d, lambda)
print(paste0("t_star_star is ", t_star_star_comp, " and t_star_star_star is ", t_star_star_star_comp))

if(y < y_hat){
  t_thresh <- t_star_star_star_comp
  t_thresh_label <- "t***"
}
if(y >= y_hat){
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
# print(paste0("y: ", y))


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
## Strong defendent

