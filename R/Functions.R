#' CellMarker Enrichment
#'
#' Given a list of markers, indentifies the putative cell type using the comphensive catlog creatd in
#'
#' @param gene.list character vector with a list of \bold{Gene Symbols}
#' @param species Can either be \emph{human} or \emph{mouse} for the moment.
#'
#' @return  A \emph{list} two elememts. First element (enrichments) is a data.frame structured as follows:
#' \itemize{
#' \item {\bold{p.value}} {P-value of the enrichment, the lower the better}
#' \item {\bold{overlap}} {How many common genes between your costom list and the markers of the given cell type}
#' \item {\bold{signature}} {Total markers of the given cell type }
#' \item {\bold{genes}} {Common genes}
#' }
#' The second element is a list with the actual gene names for each enrichment.
#'
#' @examples
#' results=CMenrich(gene.list=c('NEUROD1','ARX','CD79A'),species='human') .# some gene names that I came up with ....
#'
#' @export
CMenrich = function (gene.list,species)
{

  if (species=='human')
    data('human.CMenrich')
  else
    data('mouse.CMenrich')

  #print(gene.list)

  gene.list=intersect(gene.list,all.markers)

  #print(gene.list)

  pop.size=length(all.markers)
  samp.size=length(gene.list)

  p.val=rep(0,length(markers))
  samp.hits=rep(0,length(markers))
  pop.hits=rep(0,length(markers))
  genes=list()

  for (k in 1:length(markers))
  {
    samp.hits[k]=length(intersect(gene.list,markers[[k]]))
    pop.hits[k]=length(markers[[k]])
    p.val[k]=phyper(samp.hits[k], pop.hits[k], pop.size-pop.hits[k], samp.size,lower.tail = FALSE)+dhyper(samp.hits[k], pop.hits[k], pop.size-pop.hits[k], samp.size)
    genes[[k]]=intersect(gene.list,markers[[k]])
  }

  ix=order(p.val)

  enrich=data.frame(Cell.type=cell.type, p.value=p.val,overlap=samp.hits,signature=pop.hits)
  enrich=enrich[ix,]

  result=list(enrichments=enrich,genes=genes[ix])

}
