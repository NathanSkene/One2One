#' Checks that the species names, i.e. "mouse" can be linked 1:1 to a species name from the MGI homolog table
#'
#' Returns the linked species name from the MGI homolog table, i.e. if you pass "mouse" it returns "mouse, laboratory".
#' Essentially it's a convenience function so you don't need to write out "mouse, laboratory". It's intended primarily
#' for use within the packages other functions and should not be generally called seperately.
#'
#' @param requestedSpecies String with the name of the first species, i.e. "mouse"
#' @param allSpecies String with the names of all species from the MGI homolog table
#'
#' @return The name of the species in the MGI homolog table.
#'
#' @examples
#' allHomologs = load.homologs()
#' species1 = check.species(requestedSpecies="mouse",allSpecies=allHomologs$`Common Organism Name`)
#'
#' @export
check.species <- function(requestedSpecies,allSpecies){
    allSpecies    = unique(allSpecies)
    simpleSpecies = gsub(",.*","",allSpecies)
    species = data.frame(all=allSpecies,simple=simpleSpecies)

    # Check species only refers to one of the species listed
    selectedSpecies = apply(species==requestedSpecies,1,sum)>0
    if(sum(selectedSpecies==1)){
        validSpecies=FALSE
        selSpec = allSpecies[selectedSpecies]
        print(sprintf("Selected species: %s",selSpec))
    }else{
        validSpecies=TRUE
        print("Available species options:")
        print(species)
        stop("ERROR: species does not match available options")
    }
    return(selSpec)
}
