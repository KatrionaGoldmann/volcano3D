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
    repos["volcano3Ddata"] = "/home/katrionagoldmann/Documents/Analyses/volcano_package/volcano3Ddata/volcano3Ddata/"
    # when public: https://github.com/KatrionaGoldmann/volcano3Ddata
    options(repos = repos)
    invisible(repos)
}