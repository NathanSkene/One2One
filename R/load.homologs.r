#' Load the MGI homolog table
#'
#' Downloads and formats the MGI homolog table
#'
#' @return The full MGI homolog table
#'
#' @examples
#' allHomologs = load.homologs()
#' @importFrom utils download.file unzip read.table
#' @export
load.homologs <- function(){
    # Load the homolog data from MGI
    #read.table("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",sep="\t",stringsAsFactors = FALSE,quote="")
    #zipped file, need to unzip in temp directory
    temp <- tempdir()
    download.file("https://github.com/NathanSkene/One2One/files/6800774/HOM_AllOrganism.rpt.txt.zip",
                    paste0(temp,"/hom_vert.zip"),mode="wb")
    unzip(paste0(temp,"/hom_vert.zip"),"HOM_AllOrganism.rpt.txt",exdir=temp)
    hom_vert <- read.table(paste0(temp,"/HOM_AllOrganism.rpt.txt"),
                                    sep="\t",fill=TRUE,quote="",
                                    stringsAsFactors = FALSE)

    colnames(hom_vert) = hom_vert[1,]
    hom_vert = hom_vert[-1,]

    # The table is badly formatted, so drop rubbish
    hom_vert = suppressWarnings(hom_vert[!is.na(as.numeric(hom_vert$`HomoloGene ID`)),])

    # Provide some output to help the user
    print("Species for which homology classes are available:")
    print(unique(hom_vert$`Common Organism Name`))
    return(hom_vert)
}
