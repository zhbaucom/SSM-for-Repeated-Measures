library(stringr)

bdir <- function(ndir = "State-Space-Methods"){
  cdir <- str_split(getwd(), "/")[[1]]
  udir <- cdir[1:which(cdir == "State-Space-Methods")]
  setwd(paste(udir, collapse = "/"))
}



