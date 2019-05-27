cell_markers <- read.delim("C:/ncmls/varie/cell marker/Human_cell_markers.txt", stringsAsFactors=FALSE)
cell.type=c()
for (k in 1:nrow(cell_markers))
  cell.type[k]=sprintf('%s, %s, %s',cell_markers$tissueType[k],cell_markers$cancerType[k],cell_markers$cellName[k])

markers=list()
for (k in 1:nrow(cell_markers))
{
  dummy=unlist(strsplit(cell_markers$geneSymbol[k],','))
  for (j in 1:length(dummy))
    dummy[j]=gsub(" ", "", dummy[j], fixed = TRUE)
  markers[[k]]=dummy
}

all.markers=c()
for (k in 1:length(markers))
  all.markers=c(all.markers,markers[[k]])
all.markers=unique(all.markers)


rm(cell_markers,j,k,dummy)
