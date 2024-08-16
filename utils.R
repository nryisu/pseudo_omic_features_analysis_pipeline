#convert human gene to mouse 
convertHumanGeneList <- function(x){ 
        require("biomaRt")
        human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
        humanx <- unique(genesV2[, 2])
        return(humanx)
}

#based on 'Gapdh' to normalize 
GAPDHnorm <- function ( cts , gapdhThreshold = 0 ) {
  gapdhVector <- as.numeric ( cts [ "Gapdh" , ] )
  okay <- gapdhVector > gapdhThreshold
  cts <- cts [ , okay ]
  apply ( cts , 2, function ( x ) { x / x [ "Gapdh" ] } )
}


