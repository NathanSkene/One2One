#' Get gene symbol from homologID
#'
#' Most outputs of this package are in the format of homolog IDs. To be useful they need to be converted to gene symbols. This functon does that.
#'
#' @param homoIDs The list of homolog IDs to be converted to gene symbols
#' @param species1 String with the name of the first species, i.e. "mouse"
#' @param species2 String with the name of the second species, i.e. "human"
#' @param ortholog_data Output of analyse.orthology() function
#' @param symbol_species The name of the species for which gene symbols should be returned. I.e. "mouse"
#'
#' @return List. The gene symbols equivilenet to homolog IDs in homoIDs
#'
#' @examples
#' allHomologs = load.homologs()
#' ortholog_data = analyse.orthology("human","mouse",allHomologs)
#' homoIDs_bg     = ortholog_data$species2_allGenes
#' symbols_bg = get.symbol.from.homologID(homoIDs_bg,"human","mouse",ortholog_data,"mouse")
#'
#' @export
get.symbol.from.homologID <- function(homoIDs,species1,species2,ortholog_data,symbol_species){
    # Which species do the symbols need to be returned from? I.e. MGI or HGNC symbols?
    whichSpecies = which(c(species1,species2)==symbol_species)+1
    if(length(whichSpecies)==0){stop("ERROR: symbol_species must be the same as either species1 or species2")}
    # These are the gene symbols:
    symbols = ortholog_data$orthologs_all[ortholog_data$orthologs_all$`HomoloGene ID` %in% homoIDs,whichSpecies]
    return(unique(symbols))
}
