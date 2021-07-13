#' Analyse orthology between two species
#'
#' Returns a table of 1:1 orthologs along with gene sets (in Homologene ID format) detailing how genes differ between the two species
#'
#' @param species1 String with the name of the first species, i.e. "mouse"
#' @param species2 String with the name of the second species, i.e. "human"
#' @param allHomologs Output of load.homologs() function
#'
#' @return List. The most important thing is orthologs_one2one.
#'
#' @examples
#' ortholog_data = analyse.orthology("human","mouse")
#' @import dplyr
#' @import magrittr
#' @export
analyse.orthology <- function(species1="mouse",species2="human",allHomologs=load.homologs()){
    #species1 = "mouse"
    #species2 = "human"

    species1 = check.species(requestedSpecies=species1,allSpecies=allHomologs$`Common Organism Name`)
    species2 = check.species(requestedSpecies=species2,allSpecies=allHomologs$`Common Organism Name`)

    species1_srt = gsub(",.*","",species1)
    species2_srt = gsub(",.*","",species2)

    hom_vert = allHomologs[allHomologs$`Common Organism Name` %in% c(species1,species2),c("HomoloGene ID","Common Organism Name","Symbol")] %>% unique
    species1_hom  = hom_vert[hom_vert$`Common Organism Name`==species1,] %>% dplyr::rename(species1.symbol=Symbol)
    species2_hom  = hom_vert[hom_vert$`Common Organism Name`==species2,] %>% dplyr::rename(species2.symbol=Symbol)
    colnames(species1_hom)[3]=sprintf("%s.symbol",species1_srt)
    colnames(species2_hom)[3]=sprintf("%s.symbol",species2_srt)
    allHomologsInSpecies = merge(species1_hom[,c(1,3)],species2_hom[,c(1,3)],by="HomoloGene ID",all=TRUE)

    species1_allGenes = unique(species1_hom$`HomoloGene ID`)
    species2_allGenes = unique(species2_hom$`HomoloGene ID`)
    total_numS1genes = length(species1_allGenes)
    total_numS2genes = length(species2_allGenes)
    print(sprintf("Full dataset contains %s genes from %s",total_numS1genes,species1_srt))
    print(sprintf("Full dataset contains %s genes from %s",total_numS2genes,species2_srt))

    # Keep only genes which are expressed in both species
    shared_genes = intersect(species1_hom$`HomoloGene ID`,species2_hom$`HomoloGene ID`)
    species1_sharedHom = species1_hom[species1_hom$`HomoloGene ID` %in% shared_genes,]
    species2_sharedHom = species2_hom[species2_hom$`HomoloGene ID` %in% shared_genes,]

    # Find genes which were deleted in one species
    species1_present_species2_deleted = setdiff(species1_allGenes,shared_genes)
    species2_present_species1_deleted = setdiff(species2_allGenes,shared_genes)
    print(sprintf("%s genes which are present in %s are deleted in %s",length(species1_present_species2_deleted),species1_srt,species2_srt))
    print(sprintf("%s genes which are present in %s are deleted in %s",length(species2_present_species1_deleted),species2_srt,species1_srt))

    # Get frequency of genes being duplicated
    species1_freq = data.frame(table(species1_hom$`HomoloGene ID`)) %>% dplyr::rename(HomoloGene.ID=Var1)
    species2_freq = data.frame(table(species2_hom$`HomoloGene ID`)) %>% dplyr::rename(HomoloGene.ID=Var1)
    species1_freq = species1_freq[order(species1_freq$Freq,decreasing = TRUE),]
    species2_freq = species2_freq[order(species2_freq$Freq,decreasing = TRUE),]
    species1_duplicated_homoloID = as.character(species1_freq[species1_freq$Freq>1,]$HomoloGene.ID)
    species2_duplicated_homoloID = as.character(species2_freq[species2_freq$Freq>1,]$HomoloGene.ID)
    print(sprintf("%s genes are duplicated in %s",length(species1_duplicated_homoloID),species1_srt))
    print(sprintf("%s genes are duplicated in %s",length(species2_duplicated_homoloID),species2_srt))

    # Drop genes that have more than one entry per species
    species1_onceOnly = species1_freq %>% dplyr::filter(Freq==1) %>% .[,"HomoloGene.ID"] %>% as.character()
    species2_onceOnly = species2_freq %>% dplyr::filter(Freq==1) %>% .[,"HomoloGene.ID"] %>% as.character()
    oncePerSpecies = intersect(species1_onceOnly,species2_onceOnly)
    species1_121 = species1_sharedHom[species1_sharedHom$`HomoloGene ID` %in% oncePerSpecies,]
    species2_121 = species2_sharedHom[species2_sharedHom$`HomoloGene ID` %in% oncePerSpecies,]
    colnames(species1_121)[3]=sprintf("%s.symbol",species1_srt)
    colnames(species2_121)[3]=sprintf("%s.symbol",species2_srt)

    # Which genes are duplicated in one species, but not the other
    species1_onceOnly_species2_dup = intersect(species1_onceOnly,species2_duplicated_homoloID)
    species2_onceOnly_species1_dup = intersect(species2_onceOnly,species1_duplicated_homoloID)
    species2_dup_species1_dup = intersect(species2_duplicated_homoloID,species1_duplicated_homoloID)

    # Get merged listing of 1:1 homologs
    merged_homologs = merge(species1_121[,c(1,3)],species2_121[,c(1,3)],by="HomoloGene ID")
    print(sprintf("%s are shared 1:1 between the two species",dim(merged_homologs)[1]))

    # Prepare results
    allRes = list()
    allRes$orthologs_all     = allHomologsInSpecies
    allRes$orthologs_one2one = merged_homologs
    allRes$species1_allGenes = species1_allGenes
    allRes$species2_allGenes = species2_allGenes
    allRes$shared_genes      = shared_genes
    allRes$species1_present_species2_deleted = species1_present_species2_deleted
    allRes$species2_present_species1_deleted = species2_present_species1_deleted
    allRes$species1_onceOnly_species2_dup    = species1_onceOnly_species2_dup
    allRes$species2_onceOnly_species1_dup    = species2_onceOnly_species1_dup
    allRes$species2_dup_species1_dup         = species2_dup_species1_dup
    return(allRes)
}
