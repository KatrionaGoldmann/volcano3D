# This function adds the volcano3Ddata source to list of libraries that will be 
# searched to packages when installing
.onLoad <- function(libname, pkgname) {
    repos = getOption("repos")
    repos["volcano3Ddata"] = "https://github.com/KatrionaGoldmann/volcano3Ddata"
    options(repos = repos)
    invisible(repos)
}
