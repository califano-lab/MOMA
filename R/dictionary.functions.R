utils::globalVariables(c("gene.map"))

#' Convert from entrez ids to hugo gene names
#' @importFrom utils data
#' @param entrez.ids : vector of entrez ids requires hugo2entrez to be loaded
#' @examples
#' mapEntrez(c("29974", "5728"))
#' @seealso \code{\link[moma]{mapHugo}}
#' @return : vector of hugo gene names
#' @export
mapEntrez <- function(entrez.ids) {
  utils::data("gene.map")
  # make sure no empty spaces around the entrez ids
  entrez.ids <- gsub(" ", "", entrez.ids)
  idx <- match(as.character(entrez.ids), as.character(gene.map$Entrez.IDs))
  mapped <- gene.map$Gene.Symbol[idx]
  
  if(sum(is.na(mapped)) > 0 ) {
    print("Some entrez ids not mapped to genenames! Replaced with their original input")
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
#' @examples
#' mapHugo(c("A1CF","PTEN"))
#' @seealso \code{\link[moma]{mapEntrez}}
#' @return : vector of entrez ids 
#' @export
mapHugo <- function(hugo.ids) {
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


