#' Get MitoCarta Mouse Data
#'
#' This function downloads the MitoCarta mouse data from a specified URL, processes the data, and saves it to a local file.
#' It also returns the processed data.
#'
#' @return A data frame containing the processed MitoCarta mouse data.
#' 
#'
#' @examples
#' \dontrun{
#' get_MitoCarta_mouse()
#' }
get_MitoCarta_mouse <- function() {
    # URL of the Excel file
    url <-"https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls"

    # Temporary file to store the downloaded Excel file
    temp_file <- tempfile()

    # Download the file
    download.file(url, temp_file, mode = "wb")

    # Read the Excel file
    data <- readxl::read_excel(temp_file,sheet = "A Mouse MitoCarta3.0")

    # Don't forget to remove the temporary file when you're done with it
    unlink(temp_file)

    # Create a copy of the MitoCarta3.0_MitoPathways column
    data <- dplyr::mutate(data, MitoPathways_copy = MitoCarta3.0_MitoPathways)

    # Separate the gene sets into different rows
    data <- tidyr::separate_rows(data, MitoPathways_copy, sep = " \\| ")

    # Check if the ">" character exists in the MitoPathways_copy column
    data <- dplyr::mutate(data, GeneSet = ifelse(stringr::str_detect(MitoPathways_copy, ">"), NA, MitoPathways_copy))

    # Create another copy of the MitoPathways_copy column
    data <- dplyr::mutate(data, MitoPathways_copy2 = MitoPathways_copy)

    # Extract the gene set from the ontology
    data <- tidyr::separate(data, MitoPathways_copy, into = c("Group", "SubGroup", "GeneSet"), sep = " > ", extra = "merge", fill = "right") 
    data <- dplyr::mutate(data, GeneSet = dplyr::coalesce(GeneSet, stringr::str_extract(MitoPathways_copy2, "[^>]+$")))
    data <- dplyr::mutate(data, GeneSet = trimws(GeneSet)) 

    MitoPathways <- dplyr::filter(data, GeneSet != "0") %>% dplyr::select( "GeneSet", "Symbol", "EnsemblGeneID")
    
    save(MitoPathways, file = "./workflow/scripts/RNAscripts/data/mito_gs.RData",compress = T)

    # Print the data
    print(data)
    
    return(MitoPathways)
}
#' Get MitoCarta Mouse Data
#'
#' This function downloads the MitoCarta mouse data from a specified URL, processes the data, and saves it to a local file.
#' It also returns the processed data.
#'
#' @return A data frame containing the processed MitoCarta mouse data.
#' 
#'
#' @examples
#' \dontrun{
#' get_MitoCarta_mouse()
#' }
get_MitoCarta_human <- function() {
    # URL of the Excel file
    url <-"https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls"
    
    # Temporary file to store the downloaded Excel file
    temp_file <- tempfile()
    
    # Download the file
    download.file(url, temp_file, mode = "wb")
    
    # Read the Excel file
    data <- readxl::read_excel(temp_file,sheet = "A Human MitoCarta3.0")
    
    # Don't forget to remove the temporary file when you're done with it
    unlink(temp_file)
    
    # Create a copy of the MitoCarta3.0_MitoPathways column
    data <- dplyr::mutate(data, MitoPathways_copy = MitoCarta3.0_MitoPathways)
    
    # Separate the gene sets into different rows
    data <- tidyr::separate_rows(data, MitoPathways_copy, sep = " \\| ")
    
    # Check if the ">" character exists in the MitoPathways_copy column
    data <- dplyr::mutate(data, GeneSet = ifelse(stringr::str_detect(MitoPathways_copy, ">"), NA, MitoPathways_copy))
    
    # Create another copy of the MitoPathways_copy column
    data <- dplyr::mutate(data, MitoPathways_copy2 = MitoPathways_copy)
    
    # Extract the gene set from the ontology
    data <- tidyr::separate(data, MitoPathways_copy, into = c("Group", "SubGroup", "GeneSet"), sep = " > ", extra = "merge", fill = "right") 
    data <- dplyr::mutate(data, GeneSet = dplyr::coalesce(GeneSet, stringr::str_extract(MitoPathways_copy2, "[^>]+$")))
    data <- dplyr::mutate(data, GeneSet = trimws(GeneSet)) 
    
    MitoPathways <- dplyr::filter(data, GeneSet != "0") %>% dplyr::select( "GeneSet", "Symbol", "EnsemblGeneID_mapping_version_20200130")
    
    MitoPathways_hs <- MitoPathways |> dplyr::rename(EnsemblGeneID = "EnsemblGeneID_mapping_version_20200130" )
    
    save(MitoPathways_hs, file = "./workflow/scripts/RNAscripts/data/mito_gs_hs.RData",compress = T)
    
    # Print the data
    print(data)
    
    return(MitoPathways_hs)
}
