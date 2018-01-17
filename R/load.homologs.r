#' Load the MGI homolog table
#'
#' Downloads and formats the MGI homolog table
#'
#' @return The full MGI homolog table
#'
#' @examples
#' allHomologs = load.homologs()
#'
#' @export
load.homologs <- function(){
    # Load the homolog data from MGI
    #hom_vert = suppressWarnings(read.table("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",sep="\t",fill=TRUE,stringsAsFactors = FALSE))
    hom_vert = read.table("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",sep="\t",stringsAsFactors = FALSE,quote="")
    colnames(hom_vert) = hom_vert[1,]
    hom_vert = hom_vert[-1,]

    # The table is badly formatted, so drop rubbish
    hom_vert = suppressWarnings(hom_vert[!is.na(as.numeric(hom_vert$`HomoloGene ID`)),])

    # Provide some output to help the user
    print("Species for which homology classes are available:")
    print(unique(hom_vert$`Common Organism Name`))
    return(hom_vert)
}
