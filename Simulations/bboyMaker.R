
# knitr::opts_chunk$set(eval = FALSE)
SimPrefix <- "ItoutNACC4 "

flist <- list.files("Simulations/SimOutReform")
flist <- flist[grepl(SimPrefix,flist)][c(2,3,4,5,6,7,10,11,9)]
B <- readRDS("NACC/B.RDS")
B <- c(B[-1], B[1])
SimList <- map(flist, ~readRDS(paste("Simulations/SimOutReform/",.x, sep = "")))
# map(Itout, "Fails")
names(SimList) <- flist


bboy <- map(names(SimList),function(s2){
  sl <- SimList[[s2]]
  OB <- map(c("Estimate", "UCL", "LCL"), function(ll){
    out <- plyr::adply(sl[[ll]], 1:3)
    out$Variable <- ll
    out
  }) %>%
    do.call("rbind", .)
  OB$Simulation <- s2
  OB
}) %>%
  do.call("rbind", .) %>%
  spread("Variable", "V1") %>%
  rename(Parameter = X1, Method = X2, Iteration = X3) %>%
  mutate(
    Method = factor(c("LME", "AR(1)", "LLT", "Bayes", "Part2", "Part4", "Part10")[Method], levels = c("LME", "AR(1)", "LLT", "Bayes", "Part2", "Part4", "Part10")),
    TrueParam = B[Parameter],
    CI.length = UCL - LCL,
    Covered = LCL < TrueParam & UCL > TrueParam,
    Bias = Estimate - TrueParam
  ) %>%
  mutate(
    GR = str_split(Simulation, "[.]"),
    GR = map_chr(GR, 1),
    GR = str_split(GR, " "),
    GR = map(GR, ~.x[-1]),
    V1 = map_chr(GR, 1),
    V2 = map_chr(GR, 2),
    plabel = map_chr(GR, function(x){
      if(any(grepl("AR1", x))){
        case_when(
          grepl("None", x[2]) ~ "$\\sigma^2 = 1$, $\\rho = 0$",
          grepl("Small", x[2]) ~ "$\\sigma^2 = 1$, $\\rho = 0.1$",
          grepl("Medium", x[2]) ~ "$\\sigma^2 = 1$, $\\rho = 0.5$",
          grepl("Large", x[2]) ~ "$\\sigma^2 = 1$, $\\rho = 0.9$"
        )
      }else{
        paste("$\\sigma^2_{\\epsilon} = ", x[1], ", \ \\sigma^2_{\\eta} = ", x[2], "$", sep = "")
      }
    }),
    plabel = factor(plabel)
  ) %>%
  select(-GR) %>%
  mutate(
    sigeps = case_when(
      V1 == "AR1" ~ "$\\sigma^1 = 1$",
      TRUE ~ paste("$\\sigma^2_\\varepsilon = ", V1, "$", sep = "")
    ),
    sigeta = case_when(
      V1 == "AR1" & grepl("None", V2) ~ "$\\rho = 0$",
      V1 == "AR1" & grepl("Small", V2) ~ "$\\rho = 0.1$",
      V1 == "AR1" & grepl("Medium", V2) ~ "$\\rho = 0.5$",
      V1 == "AR1" & grepl("Large", V2) ~ "$\\rho = 0.9$",
      TRUE ~ paste("$\\sigma^2_\\eta =", V2, "$", sep = "")
    )
  )

# 
# bboy %>%
#   filter(Method == "LME", Simulation == "ItoutNACC4 3 0.RDS") %>%
#   tibble() %>%
#   group_by(Parameter) %>%
#   summarise(Mean = mean(Estimate))
# 
# B
# 
# unique(bboy$Simulation)




levels(bboy$plabel) <- TeX(levels(bboy$plabel))
bboy$plabel <- factor(bboy$plabel, levels = levels(bboy$plabel)[c(
  which(!grepl("rho", levels(bboy$plabel))),
  which(grepl("rho", levels(bboy$plabel)))
)])

saveRDS(bboy, "Simulations/bboy.RDS")






