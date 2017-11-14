e <- new.env()

local({
   path <- paste0(Sys.getenv("HOME"), "/Dropbox/Rworld/Packages/pine/R/")
   source(paste0(path, "pine.R"), local=TRUE)
   source(paste0(path, "pine12.R"), local=TRUE)
   source(paste0(path, "pine36.R"), local=TRUE)
   gens       <- c(2, 10, 100, 101, 1000, 1001)
   p12        <- lapply(gens, pine12, plot=FALSE)
   names(p12) <- paste0("Gen", gens)
   p36        <- lapply(gens, pine36, plot=FALSE)
   names(p36) <- paste0("Gen", gens)
}, envir=e)

pine12_solutions <- e$p12
pine36_solutions <- e$p36
rm(e) # cleanup
