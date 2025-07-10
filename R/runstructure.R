#' @title Runs a STRUCTURE analysis using a genlight object
#'
#' @description
#' This function takes a genlight object and runs a STRUCTURE analysis based on
#' functions from \code{strataG}
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ... Parameters to specify the STRUCTURE run (check \code{structureRun}
#'  within strataG.
#' for more details). Parameters are passed to the \code{structureRun} function.
#' For example you need to set the k.range and the type of model you would like
#' to run (noadmix, locprior) etc. If those parameter names do not tell you
#' anything, please make sure you familiarize with the STRUCTURE program
#' (Pritchard 2000).
#' @param exec Full path and name+extension where the structure executable is
#' located. E.g. \code{'c:/structure/structure.exe'} under Windows. For Mac and
#' Linux it might be something like \code{'./structure/structure'} if the
#' executable is in a subfolder 'structure' in your home directory
#' [default working directory "."].
#' @param plot.out Create an Evanno plot once finished. Be aware k.range needs
#' to be at least three different k steps [default TRUE].
#' @param plot_theme Theme for the plot. See details for options
#' [default theme_dartR()].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Set verbosity for this function (though structure output
#' cannot be switched off currently) [default NULL]
#' @details The function is basically a convenient wrapper around the beautiful
#' strataG function \code{structureRun} (Archer et al. 2016). For a detailed
#' description please refer to this package (see references below).
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @export
gl.run.structure <- function(x,
                             ...,
                            # exec = ".",
                             plot.out = TRUE,
                             plot_theme = theme_dartR(),
                             save2tmp = FALSE,
                             verbose = NULL) {
    
    pkg <- "purrr"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }

        # check that Structure is installed
     #   structure <- file.exists(exec)

        structure <- Sys.which("structure") ############ FOR CHECKING ###################
        exec <- structure
  
        # SET VERBOSITY
        verbose <- dartR.base::gl.check.verbosity(verbose)
        
        # FLAG SCRIPT START
        funname <- match.call()[[1]]
        dartR.base::utils.flag.start(func = funname,
                         build = "Jody",
                         verbose = verbose)
        
        # CHECK DATATYPE
   #     datatype <- utils.check.datatype(x, verbose = verbose)
        
    #    if (datatype != "SNP") {
     #       stop(error(
      #          "You need to provide a SNP genlight object (ploidy=2)!"
       #     ))
        #}
############################################

   if (tools::file_ext(x) == "csv") {
      file <- readr::read_csv(x)
   } else if (tools::file_ext(x) == "xlsx") {
      file <- readxl::read_excel(x)
   } else {
      stop("Input file should be in csv or xlsx format.")
   }


   file <- lapply(file, function(x) gsub("|", "/", x, fixed = TRUE))
   file <- as.data.frame(file)
   
   file[is.na(file)] <- "N"
   file <- file %>%
      mutate(across(everything(), ~ case_when(
         . == "N/A" ~ "N", 
         . == "NA" ~ "N",
         TRUE ~ .x))) %>%
      rename(Ind = 1, Pop = 2)

   ind <- as.character(file$Ind)
   pop <- as.character(file$Pop)
   fsnps_geno <- file[, 3:ncol(file)]
   
   fsnps_gen <- adegenet::df2genind(fsnps_geno, ind.names = ind, pop = pop, sep = "/", NA.char = "N", ploidy = 2, type = "codom")
   fsnps_gen@pop <- as.factor(file$Pop)

 # gg <- dartR::gi2gl(fsnps_gen, verbose = 0)
  gg <- fsnps_gen

  ############################################
        # DO THE JOB
      #  gg <- utils.structure.genind2gtypes(gl2gi(x, verbose = 0))
        
        sr <- dartR.base::utils.structure.run(gg, exec = exec, ...)
        
        ev <- dartR.base::utils.structure.evanno(sr)
        
        pa <- ((ev$plots$mean.ln.k + ev$plots$mean.ln.k) / 
                 (ev$plots$ln.ppk + ev$plots$delta.k)) + plot_theme
        
        # PRINTING OUTPUTS
        if (plot.out) {
            suppressMessages(print(pa))
        }
        
        # SAVE INTERMEDIATES TO TEMPDIR
        if (save2tmp & plot.out) {
            # check for '/' in match.call
            mc <- gsub("/", ":", as.character(funname))
            mc <- gsub(":", "-", mc)
            nmc <- gsub("/", "_over_", names(funname))
            nmc <- gsub(":", "-", nmc)
            
            # creating temp file names
            temp_plot <-
                tempfile(pattern = paste0("Plot", paste0(nmc, "_", mc,
                                                         collapse = "_")))
            
            # saving to tempdir
            saveRDS(pa, file = temp_plot)
            if (verbose >= 2) {
                cat(
                    report(
                        "  Saving the plot in ggplot format to the tempfile as",
                        temp_plot,
                        "using saveRDS\n"
                    )
                )
                cat(
                    report(
                        "  NOTE: Retrieve output files from tempdir using 
                        gl.list.reports() and gl.print.reports()\n"
                    )
                )
            }
        }
        
        # FLAG SCRIPT END
        if (verbose >= 1) {
            cat(report("Completed:", funname, "\n\n"))
        }
        
        # RETURN
        return(sr)
    
}
