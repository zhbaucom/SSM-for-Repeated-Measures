Sym.Enf <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}
