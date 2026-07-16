MC <- function(k, notionals, EAD, LGD, PDs, one_minus_alpha, assets = NULL, Cholesky = NULL, DS = NULL, betas = NULL, sigma_M = NULL, sigmas = NULL, type = NULL, sim = 1){
    single_simulation <- function() {
        if(length(notionals) == 1) {
        notionals <- rep(notionals, length.out = assets)
        }
        if(length(PDs) == 1) {
            PDs <- rep(PDs, length.out = assets)
        }
        if(type == 1) {
            simulation <- replicate(k, {
                indfun <- ifelse(
                    runif(length(notionals)) <= PDs, 1, 0
                )
                sum(indfun*notionals*EAD*LGD)
            })
        } else if(type == 2) {
            simulation <- replicate(k, {
                indfun <- ifelse(
                    rnorm(length(notionals))%*%Cholesky <= qnorm(PDs), 1, 0
                )
                sum(indfun*notionals*EAD*LGD)
            })
        } else if(type == 3) {
            if(length(unique(notionals)) > 1) {
                warning("observe the notionals of the loans! the function is applicable only for assets with same notinal value. NA returned")
                return(NA)
            } else {
                simulation <- replicate(k, {
                    notionalsDS <- head(notionals, n = length(notionals)*DS)
                    indfun <- ifelse(
                        runif(length(notionalsDS)) <= head(PDs, n = length(notionalsDS)), 1, 0
                    )
                    sum(indfun*EAD*LGD*notionalsDS)/DS
                })
            }
        } else if(type == 4) {
            Dtd <- qnorm(PDs)
            W_i <- betas*sigma_M/sigmas
            simulation <- replicate(k, {
                Y <- rnorm(length(notionals)+1)
                indfun <- ifelse(
                    W_i*Y[length(notionals)+1]+(1-W_i)*Y[1:length(notionals)] <= Dtd, 1, 0
                )
                sum(indfun*notionals*EAD*LGD)
            })
        }
    Average_Portfolio_Loss <- mean(simulation)
    ELR <- Average_Portfolio_Loss/sum(notionals)
    VaR <- unname(quantile(simulation, probs = one_minus_alpha, type = 7))
    return(c(Average_Portfolio_Loss = Average_Portfolio_Loss, ELR = ELR, VaR = VaR))
    }
    if(sim == 1) {
        res <- single_simulation()
        cat("Average_Portfolio_Loss\t\tELR\t\tVaR\n", res["Average_Portfolio_Loss"], "\t\t\t\t", res["ELR"], "\t", res["VaR"], "\n")
        return(invisible(res))
    } else {
        results <- replicate(sim, single_simulation())
        VaRseries <- results["VaR", ]
        expVaR <- mean(VaRseries)
        MC_error <- sd(VaRseries)
        cat("Expected VaR:", expVaR, " | MC Error:", MC_error, "\n")
        return(invisible(VaRseries))
    }
}