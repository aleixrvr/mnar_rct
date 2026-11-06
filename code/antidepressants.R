library(data.table)
d <- fread("antidepressant_data.csv"); n <- nrow(d)

# naive (complete cases)
cc <- d[a == 1]
p <- cc[, .(p = mean(y)), by = t][order(t)]$p
c(naive_ate    <- p[2] - p[1], naive_calor  <- qlogis(p[2]) - qlogis(p[1]))


# proposed plug-in estimators, W = (x, g)  [Algorithm, Section 3]
st <- d[, {
  cs  <- .SD[a == 1]
  pi  <- mean(t == 1)
  pO1 <- mean(cs[y == 1]$t); pO0 <- mean(cs[y == 0]$t)
  r0  <- (pi - pO0) / (pO1 - pO0)
  q1  <- mean(cs[t == 1]$y); q0 <- mean(cs[t == 0]$y)
  .(n_w = .N, n1 = sum(t == 1), n0 = sum(t == 0),
    rho1 = pO1 / pi * r0, rho0 = (1 - pO1) / (1 - pi) * r0,
    logOR = qlogis(q1) - qlogis(q0))
}, by = .(x, g)]

Phi1 <- st[, sum(rho1 * n1 / sum(n1))]
Phi0 <- st[, sum(rho0 * n0 / sum(n0))]

c(ATE = Phi1 - Phi0, CALOR = st[, sum(n_w / n * logOR)])
