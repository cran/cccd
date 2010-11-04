
.First.lib <- function(lib, pkg) {
	library.dynam("cccd",pkg,lib)
    gver <- read.dcf(file=system.file("DESCRIPTION", package=pkg), 
                      fields="Version")
    gdate <- read.dcf(file=system.file("DESCRIPTION", package=pkg), 
                      fields="Date")
    gpack <- read.dcf(file=system.file("DESCRIPTION", package=pkg), 
                      fields="Packaged")
	 gpack <- strsplit(gpack,split=";")[[1]][1]
	 cat(pkg,"version",gver,"\nDate:",gdate,"\nPackaged:",gpack,"\n") 
}
