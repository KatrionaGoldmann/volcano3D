#This file contains the extra features for the initial load

# Welcome message
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste("-----------------------",
                                "                      ",
                                "Welcome to volcano3D!",
                                "                      ",
                                "-----------------------", sep="\n"))
}


