#This file contains the extra features for the initial load

# Welcome message
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste("---------------------",
                                "                    ",
                                "Welcome to volcano3D",
                                "                    ",
                                "---------------------", sep="\n"))
}


.onLoad <- function(libname, pkgname) {
    repos <- getOption("repos")
    repos["volcano3Ddata"] <- "http://katrionagoldmann.github.io/volcano3Ddata"
    options(repos = repos)
    invisible(repos)
}