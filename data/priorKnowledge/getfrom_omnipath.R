#! /usr/bin/env Rscript
#args <- commandArgs(trailingOnly = TRUE)

#install.packages("BiocManager",  repos = "http://cran.us.r-project.org")
#remotes::install_github('saezlab/omnipathR')
#BiocManager::install("saezlab/decoupleR")
library(decoupleR)

sapply(c("mouse", "human"), FUN = function(organism){
  # progeny
  sapply(c(300, 500), FUN = function(top){
    res <- get_progeny(organism, top=top)
    write.csv(res, paste0("progeny", '_', organism, '_', top, '.csv'))
  })
  # collectri
  res <- get_collectri(organism)
  write.csv(res, paste0("collectri", '_', organism,'.csv'))
})

# ksn
res <- get_ksn_omnipath()
write.csv(res, paste0("ksn", '_', 'human', '.csv'))



#resource <- args[1]
#methodname = paste0("get_", resource)

#organism <- args[2]

#print(args)

#params <- list(args[3])
#sapply(c(500, 300), FUN = function(top){
#res <- do.call(methodname, c(organism, params))
#print(res)
#write.csv(res, paste0(resource, '_', organism, '_', sub("=", "", params), '.csv'))
#})

# sub("=", "", "top=300")
# sub(".*=(.*)", "\\1", "top=300")
