#' Convert from entrez ids to hugo gene names
#' 
#' @param entrez.ids : vector of entrez ids requires hugo2entrez to be loaded
#' @return : vector of hugo gene names
#' @export
map.entrez <- function(entrez.ids) {
  mapped <- c()
  for (name in entrez.ids) {
    idx <- which(as.numeric(hugo2entrez[,2]) == as.numeric(name))
    if (length(idx) < 1) {
      mapped <- c(mapped, NA)
    } else {
      mapped <- c(mapped, hugo2entrez[idx, 1])
    }
  } 
  return (as.character(mapped))
}


#' Convert from hugo gene names to entrez ids
#' 
#' @param hugo.ids : vector of entrez ids, requires hugo2entrez to be loaded
#' @return : vector of entrez ids 
#' @export
map.hugo <- function(hugo.ids) {
  mapped <- c()
  for (name in hugo.ids) {
    idx <- which(hugo2entrez[,1] == name)
    if (length(idx) < 1) {
      mapped <- c(mapped, NA)
    } else {
      mapped <- c(mapped, as.numeric(hugo2entrez[idx, 2]))
    }
  } 
  return (as.character(mapped))
}