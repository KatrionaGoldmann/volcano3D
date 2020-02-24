#This file contains the extra features for the initial load

# Welcome message
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste("---------------------",
                                "                    ",
                                "Welcome to volcano3D",
                                "                    ",
                                "---------------------", sep="\n"))
}


# Add the volcano3Ddata source to list of repos install loads from
.onLoad <- function(libname, pkgname) {
    repos = getOption("repos")
    repos["volcano3Ddata"] = "http://KatrionaGoldmann.github.io/volcano3Ddata"
    options(repos = repos)
    invisible(repos)
}