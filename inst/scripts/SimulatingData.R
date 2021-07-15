## Simulation of HDX data
require(S4Vectors)
require(Spectra)

## Number of peptides
n <-  500

## min/max length
minL <- 5
maxL <- 25

## possible amino acids
AA <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
        "F", "P", "S", "T", "W", "Y", "V")

## Number of time points
Tp <- 3

## Number of replicates
R <- 2

## probability condition swapping
probz <- 0.05



simHDXdata <- function(n = n,
                       minL = minL,
                       maxL = maxL,
                       AA = AA,
                       Tp = Tp,
                       R = R,
                       probz = probz){
centroids <- list()
valuesz <- list()

for (j in 1:n){
    
    ## storage 
    centroids[[j]] <- matrix(NA, ncol = Tp + 1, nrow = R*2)
    
    # length of amino acid
    L <- sample(minL:maxL, size = 1)
    
    ## sample amino acids
    seq <- sample(AA, size = L, replace = TRUE)
    
    ## simulate unduerated spectra
    deut0 <- lapply(1:(2*R), function(z) isotopicDistributionHDX(sequence = paste0(seq, collapse = ""),
                                     incorp = 0,
                                     charge = 2))
    cent0 <- lapply(deut0, peptidemass)
    
    centroids[[j]][, 1] <- sapply(cent0, function(x) x[,1])
    
    incorp1 <- cumsum(rdirichlet(n = 1, alpha = 20/seq.int(Tp)))
    incorp1[length(incorp1)] <- incorp1[length(incorp1)] - 0.01 # stop numerical issues
    
    z <- sample(c(0, 1), size = 1, prob = c(1 - probz, probz))
    incorp2 <- cumsum(rdirichlet(n = 1, alpha = 20/seq.int(Tp)))
    incorp2[length(incorp2)] <- incorp2[length(incorp2)] - 0.01 # stop numerical issues
    
    for(i in seq.int(Tp)){

        deutn <- lapply(1:R, function(z) isotopicDistributionHDX(sequence = paste0(seq, collapse = ""),
                                                                   incorp = incorp1[i],
                                                                   charge = 2))
        centn <- lapply(deutn, peptidemass)
        
        # condition specific effect
        
        if (z == 1){
            deutc <- lapply(1:R, function(z) isotopicDistributionHDX(sequence = paste0(seq, collapse = ""),
                                                                     incorp = incorp2[i],
                                                                     charge = 2))
            centc <- lapply(deutc, peptidemass)
        } else{
            deutc <- lapply(1:R, function(z) isotopicDistributionHDX(sequence = paste0(seq, collapse = ""),
                                                                     incorp = incorp1[i],
                                                                     charge = 2))
            centc <- lapply(deutc, peptidemass)
        }

        centroids[[j]][, i + 1] <- sapply(c(centn, centc), function(z) z[,1]) 
    }
    valuesz[[j]] <- z
    
}

 return(res <- list(centroids = centroids, valuesz = valuesz))   
}
