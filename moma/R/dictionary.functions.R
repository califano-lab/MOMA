#' Convert from entrez ids to hugo gene names
#' @importFrom utils data
#' @param entrez.ids : vector of entrez ids requires hugo2entrez to be loaded
#' @return : vector of hugo gene names
map_entrez <- function(entrez.ids) {
  utils::data("gene.map")
  # make sure no empty spaces around the entrez ids
  entrez.ids <- gsub(" ", "", entrez.ids)
  idx <- match(as.character(entrez.ids), as.character(gene.map$Entrez.IDs))
  mapped <- gene.map$Gene.Symbol[idx]
  
  if(sum(is.na(mapped)) > 0 ) {
    #print("Some entrez ids not mapped to genenames! Replaced with their original input")
    na.idx <- which(is.na(mapped))
    for (na in na.idx) {
      old.name <- entrez.ids[na]
      mapped[na] <- old.name
    }
  }
  
  return(as.character(mapped))
}



#' Convert from hugo gene names to entrez ids
#' @importFrom utils data
#' @param hugo.ids : vector of hugo gene names, requires hugo2entrez to be loaded
#' @return : vector of entrez ids 
map_hugo <- function(hugo.ids) {
  utils::data("gene.map")
  idx <- match(hugo.ids, gene.map$Gene.Symbol)
  mapped <- gene.map$Entrez.IDs[idx]
  
  if(sum(is.na(mapped)) > 0 ) {
    print("Some genenames not mapped to entrez ids! Replaced with their original input")
    na.idx <- which(is.na(mapped))
    for (na in na.idx) {
      old.name <- hugo.ids[na]
      mapped[na] <- old.name
    }
  }
  return(as.character(mapped))
}


