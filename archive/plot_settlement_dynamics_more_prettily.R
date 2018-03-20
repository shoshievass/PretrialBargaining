library(tidyverse)
library(ggthemes)
library(latex2exp)

tmax_val = 50
J_val = 2500
k_val = 175
k_p_val = 40
lam_val = 0.01
c_val = 25
alpha_val = 0.55
alpha_st_val = 0.5

F_tau_tex <- TeX('$F_{\\tau^*}(t)$')
f_tau_tex <- TeX('$f_{\\tau^*}(t)$')
H_tau_tex <- TeX('$H_{\\tau^*}(t)$')

# F_tau_tex_big <- Tex($$'F_{\tau^{\ast}}(t) =
# \left\{
# 	\begin{array}{ll}
# 		1 + e^{-\lambda t)} (\beta(t) - 1) & \mbox{if } t < \bar{t} \\
# 		1 + e^{-\lambda t)} (1-\hat{y}) & \mbox{if } t = \bar{t}
# 	\end{array}
# \right.
# $$')
#
# f_tau_text_big <- TeX('$f_{\tau^{\ast}}(t) = \lambda e^{-\lambda t} \left[ A\beta(t)^2 - \beta(t) + 1 \right]$)
#
# H_tau_tex_big <- TeX('$H_{\tau^{\ast}}(t) = \lambda\left[ \frac{A\beta(t)^2}{1 - \beta(t)} + 1 \right]$')

t = seq(0, tmax_val, length.out = 1000)

beta_t <- function(t, c = c_val, lam = lam_val, alpha = alpha_val, alpha_st = alpha_st_val, t_bar = tmax_val, J = J_val, k = k_val, k_p = k_p_val){
  beta_den = (alpha - alpha_st)*(t_bar - t)*c + J + alpha*k - k_p
  beta = (c / lam) / beta_den
  return(beta)
}

# plot(t, beta_t(t))

F_tau_st <- function(t, lam = lam_val){
  out = 1 + exp(-lam*t)*(beta_t(t) - 1)
  return(out)
}

# plot(t, F_tau_st(t))

f_tau_st <- function(t, lam = lam_val,alpha = alpha_val, alpha_st = alpha_st_val){
  out = lam * exp(-lam * t) * ((alpha - alpha_st) * beta_t(t) - beta_t(t) + 1)
  return(out)
}

# plot(t, f_tau_st(t))

H_tau_st <- function(t, lam = lam_val,alpha = alpha_val, alpha_st = alpha_st_val){
  out = lam * ( ((alpha - alpha_st) * (beta_t(t)^2)) / (1 - beta_t(t)) + 1)
  return(out)
}

plot(t, H_tau_st(t))

y_hat <- function(k = k_val, J = J_val){
  return(k/J)
}

y_star <- function(c = c_val, lam = lam_val, alpha = alpha_val, J = J_val, k = k_val, k_p = k_p_val){
  denom = J + alpha*k - k_p
  out = (c / lam) / denom
  return(out)
}

F_tau_deadline <- function(t, lam = lam_val, k =k_val, J=J_val, tmax = tmax_val){
  y_hat <- k/J
  deadline_val <- 1 + exp(-lam *t)*(1-y_hat)
  return(deadline_val)
}

f_tau_deadline <- function(t, c = c_val, lam = lam_val, alpha = alpha_val, J = J_val, k = k_val, k_p = k_p_val,tmax = tmax_val){
  y_hat <- k/J
  denom <- J + alpha*k - k_p
  y_star <- (c / lam) / denom

  deadline_val <- 1 + exp(-lam *t)*(y_hat-y_star)
  return(deadline_val)
}


F_tau_curve_deadline_val <- F_tau_st(tmax_val)
F_tau_deadline_val <- F_tau_deadline(tmax_val)

f_tau_curve_deadline_val <- f_tau_st(tmax_val)
f_tau_deadline_val <- f_tau_deadline(tmax_val)

H_tau_curve_deadline_val <- H_tau_st(tmax_val)


# plot(t,f_tau_st(t))


F_tau_df <- data.frame( t = t,
                        F_tau = F_tau_st(t),
                        # F_tau_deadline = c(rep(NA, (tmax_val - 1)), min(1, F_tau_deadline_val)),
                        f_tau = f_tau_st(t),
                        H_tau = H_tau_st(t)
)

F_tau_df_plot <- F_tau_df %>%
  ggplot(aes(x = t, y = F_tau)) +
  geom_line() +
  geom_segment(arrow = arrow(length=unit(0.2, "cm")),aes(x = tmax_val, xend=tmax_val, y = F_tau_curve_deadline_val, yend = min(F_tau_deadline_val, .99), color = "red")) +
  # geom_point(aes(x=tmax_val,y=F_tau_deadline, color = "red")) +
  labs(
    y = F_tau_tex,
    t = ""
  ) + guides(color=FALSE) +
  theme_tufte() +
  geom_rangeframe() +
  scale_x_continuous(breaks = c()) +
  scale_y_continuous(breaks = c())

ggsave(F_tau_df_plot, file = "F_tau_df_plot.jpg")

f_tau_df_plot <- F_tau_df %>%
  ggplot(aes(x = t, y = f_tau)) +
  geom_line() +
  geom_segment(arrow = arrow(length=unit(0.2, "cm")),aes(x = tmax_val, xend=tmax_val, y = f_tau_curve_deadline_val, yend = min(f_tau_deadline_val, f_tau_curve_deadline_val+0.00007), color = "red")) +
  # geom_point(aes(x=tmax_val,y=F_tau_deadline, color = "red")) +
  labs(
    y = f_tau_tex,
    t = ""
  ) + guides(color=FALSE) +
  theme_tufte() +
  geom_rangeframe() +
  scale_x_continuous(breaks = c()) +
  scale_y_continuous(breaks = c())

ggsave(f_tau_df_plot, file = "flit_tau_df_plot.jpg")


H_tau_df_plot <- F_tau_df %>%
  ggplot(aes(x = t, y = H_tau)) +
  geom_line() +
  geom_segment(arrow = arrow(length=unit(0.2, "cm")),aes(x = tmax_val, xend=tmax_val, y = H_tau_curve_deadline_val, yend = H_tau_curve_deadline_val+0.001, color = "red")) +
  # geom_point(aes(x=tmax_val,y=F_tau_deadline, color = "red")) +
  labs(
    y = H_tau_tex,
    t = ""
  ) + guides(color=FALSE) +
  theme_tufte() +
  geom_rangeframe() +
  scale_x_continuous(breaks = c()) +
  scale_y_continuous(breaks = c())

F_tau_df_plot
f_tau_df_plot
H_tau_df_plot

ggsave(H_tau_df_plot, file = "H_tau_df_plot.jpg")

