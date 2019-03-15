#' @title visualize probabisitic mutaiton signature for the independent model
#' @description Generate visualization of mutation signatures for the model with
#'   substitution patterns and flanking bases represented by the indepenent
#'   representation.
#'
#' @param vF a matrix for mutation signature
#' @param numBases the number of flanking bases
#' @param baseCol the colour of the bases (A, C, G, T, plus strand, minus strand)
#' @param trDir the index whether the strand direction is plotted or not
#' @param charSize the size of the character
#' @param isScale the index whether the height of the flanking base is changed or not
#' @param alpha the parameter for the Renyi entropy (applicable only if the isScale is TRUE)
#' @param charLimit
#'
#'
visPMS <- function(vF, numBases, baseCol = NA, trDir, charSize = 1.2, scale = TRUE, alpha = 2,
    charLimit = 0.25) {

    if (is.na(baseCol)) {
        gg_color_hue6 <- hcl(h = seq(15, 375, length = 7), l = 65, c = 100)[1:6]
        baseCol <- c(gg_color_hue6[3], gg_color_hue6[5], gg_color_hue6[2], gg_color_hue6[1],
            gg_color_hue6[4], gg_color_hue6[6])
    }

    centerBase <- (1 + numBases)/2

    v1 <- vF[1, 1:6]
    V2 <- vF[2:(numBases), 1:4]
    A <- matrix(0, numBases, 4)
    B <- matrix(0, 4, 4)

    if (trDir == TRUE) {
        v3 <- vF[(numBases + 1), 1:2]
    }

    for (l in 1:numBases) {
        if (l < centerBase) {
            A[l, ] <- V2[l, ]
        } else if (l > centerBase) {
            A[l, ] <- V2[l - 1, ]
        }
    }
    A[centerBase, 2] <- sum(v1[1:3])
    A[centerBase, 4] <- sum(v1[4:6])

    B[2, c(1, 3, 4)] <- v1[1:3]/sum(v1[1:3])
    B[4, c(1, 2, 3)] <- v1[4:6]/sum(v1[4:6])

    num2base <- c("A", "C", "G", "T")

    renyi = function(p, tAlpha = alpha) {
        if (tAlpha == 1) {
            return(-sum(p * log2(p), na.rm = TRUE))
        } else {
            return(log(sum(p^tAlpha))/(1 - tAlpha))
        }
    }

    # sizes <- 0.5 * (2 - apply(A, MARGIN = 1, FUN = function(p, alpha = scale) {-sum(x *
    # log2(x), na.rm = TRUE)})); sizes <- sizes ** scale

    if (scale == FALSE) {
        sizes <- rep(1, numBases)
    } else {
        sizes <- 0.5 * (2 - apply(A, MARGIN = 1, FUN = renyi))
    }

    startx <- 0
    for (l in 1:numBases) {

        for (w in 1:4) {
            endx <- startx + A[l, w]
            polygon(c(startx, endx, endx, startx), c(0, 0, sizes[l], sizes[l]), col = baseCol[w],
                border = F)
            if (endx - startx > charLimit & sizes[l] > 0.5 & charSize > 0) {
                text(0.5 * (endx + startx), 0.5 * sizes[l], num2base[w], col = "white", cex = charSize)
            }
            startx <- endx
        }
        startx <- startx + 0.25
    }

    startx <- (centerBase - 1) * 1.25
    for (w in 1:4) {
        starty <- 2
        endx <- startx + A[centerBase, w]
        for (ww in 1:4) {
            endy <- starty + B[w, ww]
            polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col = baseCol[ww],
                border = F)
            if ((endy - starty > charLimit) & (endx - startx > charLimit) & charSize > 0) {
                text(0.5 * (endx + startx), 0.5 * (endy + starty), num2base[ww], col = "white",
                  cex = charSize)
            }
            starty <- endy
        }
        startx <- endx
        starty <- endy
    }

    # draw arrow
    xs <- c(1/3, 2/3, 2/3, 5/6, 1/2, 1/6, 1/3, 1/3) + (centerBase - 1) * 1.25
    ys <- c(1/4, 1/4, 1/2, 1/2, 3/4, 1/2, 1/2, 1/4) + 1

    polygon(xs, ys, col = 8, border = F)

    if (trDir == TRUE) {
        # draw direction bias
        startx <- (numBases - 1) * 1.25 + 0.24
        endx <- (numBases - 1) * 1.25 + 0.49
        starty <- 2
        endy <- starty + v3[1]
        polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col = baseCol[5],
            border = F)
        if (endy - starty > 1/8 & charSize > 0) {
            text(0.5 * (startx + endx), 0.5 * (starty + endy), "+", col = "white", cex = charSize)
        }

        startx <- (numBases - 1) * 1.25 + 0.51
        endx <- (numBases - 1) * 1.25 + 0.76
        starty <- 2
        endy <- starty + v3[2]
        polygon(c(startx, endx, endx, startx), c(starty, starty, endy, endy), col = baseCol[6],
            border = F)
        if (endy - starty > 1/8 & charSize > 0) {
            text(0.5 * (startx + endx), 0.5 * (starty + endy), "-", col = "white", cex = charSize)
        }

    }

}


#' Read the raw mutation data with the mutation feature vector format, estimate and plot both mutation signatures and their fractions
#'
#' @param numBases the number of flanking bases around the mutated position.
#' @param numSig the number of mutation signatures.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot
#' @export

hilda_plotSignature <- function(inputParam) {


    par(mar = c(0, 0, 0, 0))
    par(bg = rgb(0.9, 0.9, 0.9))
    par(xaxs = "i", yaxs = "i")
    par(mfrow = c(6, 1))


    for (i in 1:inputParam@signatureNum) {
        plot.new()
        plot.window(xlim = c(-0.3, 28.3), ylim = c(-0.3, 3.1))
        polygon(c(-0.25, 6.25, 6.25, -0.25), c(-0.25, -0.25, 3.09, 3.09), col = "white", border = "white")
        visPMS(inputParam@signatureFeatureDistribution[i, , ], numBases = inputParam@flankingBasesNum,
            trDir = inputParam@transcriptionDirection, charSize = 1)
    }

}


#' Estimate initial values for HiLDA test from the pmsignature result
#'
#' @param input a EstimatedParameters S4 class output by the pmsignature..
#' @param refGroup the indice indicating the samples in the reference group.
#' @param sigOrder the reference signature should be replaced on the first.
#' @param numSig the number of mutation signatures.
#'
#'
#' @importFrom abind abind
#' @export
#'


hilda_inits <- function(inputParam = Param, refGroup, sigOrder = NULL, ...) {

    n.sig <- inputParam@signatureNum

    if (is.null(sigOrder)) {
        sigOrder <- 1:n.sig
    }

    sig <- abind::abind(lapply(sigOrder, function(x) inputParam@signatureFeatureDistribution[x,
        , ]), along = 3)
    sampleNum <- length(inputParam@sampleList)
    fraction <- inputParam@sampleSignatureDistribution[, sigOrder]
    fraction[which(fraction == 0)] <- min(fraction[which(fraction > 0)]) * 0.01

    numSig <- inputParam@signatureNum
    caseGroup <- setdiff(1:sampleNum, refGroup)

    n.flanking <- inputParam@flankingBasesNum - 1
    n.feature <- inputParam@flankingBasesNum

    mat.alpha <- rbind(sirt::dirichlet.mle(fraction[refGroup, ])$alpha, sirt::dirichlet.mle(fraction[caseGroup,
        ])$alpha)

    inits <- list(list(p.states1 = array(sig[1, , sigOrder], dim = c(1, 6, n.sig)), p.states2 = array(sig[2:n.feature,
        1:4, sigOrder], dim = c(n.flanking, 4, n.sig)), alpha = mat.alpha), list(p.states1 = array(sig[1,
        , sigOrder], dim = c(1, 6, n.sig)), p.states2 = array(sig[2:n.feature, 1:4, sigOrder],
        dim = c(n.flanking, 4, n.sig)), alpha = mat.alpha))

    return(inits)
}


#' Apply HiLDA to statistically testing the global difference in burdens of mutation signatures between two groups
#'
#' @param inputG a MutationFeatureData S4 class output by the pmsignature.
#' @param inputParam a EstimatedParameters S4 class output by the pmsignature.
#' @param refGroup the indice indicating the samples in the reference group.
#' @param sigOrder the order of the mutational signatures.
#' @param n.iter number of total iterations per chain (default: 2000).
#' @param n.burnin length of burn (default: 0).
#' @param prob_M1 the probability of sampling the model 1 (no difference in the means)
#'
#'
#' @importFrom R2jags jags
#' @export


hilda_bayesfactor <- function(inputG = G, inputParam = Param, refGroup, sigOrder = NULL, n.iter = 2000,
    n.burnin = 0, prob_M1 = 0.5, ...) {
    n.sig <- inputParam@signatureNum

    if (is.null(sigOrder)) {
        sigOrder <- 1:n.sig
    }

    # reshape the counts of input data
    countwide <- as.data.frame(t(inputG@countData))
    colnames(countwide) <- c("type", "sample", "count")
    countlong <- reshape(countwide, idvar = "sample", timevar = "type", direction = "wide")
    countlong <- countlong[order(countlong[, 1]), ]
    countlong[is.na(countlong)] <- 0
    countlong <- countlong[, -1]

    # generate the known data for MCMC
    Num1 <- length(refGroup)
    Num2 <- length(inputG@sampleList) - Num1
    rownames(countlong) <- 1:(Num1 + Num2)
    mutationN <- as.vector(rowSums(countlong))
    caseGroup <- setdiff(1:(Num1 + Num2), refGroup)
    numBases <- length(inputG@possibleFeatures)

    X_G1 <- structure(.Data = rep(0, Num1 * max(mutationN[refGroup]) * numBases), .Dim = c(max(mutationN[refGroup]),
        numBases, Num1))
    X_G2 <- structure(.Data = rep(0, Num2 * max(mutationN[caseGroup]) * numBases), .Dim = c(max(mutationN[caseGroup]),
        numBases, Num2))

    for (i in refGroup) {
        X_G1[1:sum(countlong[i, ]), , which(refGroup == i)] <- t(inputG@featureVectorList[, rep(1:ncol(inputG@featureVectorList),
            countlong[i, ])])
    }

    for (i in caseGroup) {
        X_G2[1:sum(countlong[i, ]), , which(caseGroup == i)] <- t(inputG@featureVectorList[,
            rep(1:ncol(inputG@featureVectorList), countlong[i, ])])
    }

    # set up the MCMC
    N1 <- mutationN[refGroup]
    N2 <- mutationN[caseGroup]

    jdata <- list(I1 = Num1, I2 = Num2, K = n.sig, numStates = 6, numflank = numBases - 1, N1 = N1,
        N2 = N2, X_T = X_G1, X_B = X_G2, prob1 = prob_M1)

    # inits <- hilda_inits(inputParam, refGroup, c(2,3,1))

    var.s <- c("p.states1", "p.states2", "pM2", "alpha", "beta")

    model.fit <- R2jags::jags(model.file = system.file("models/bayesfactor.txt", package = "HiLDA"), data = jdata, parameters.to.save = var.s,
        n.chains = 2, n.iter = n.iter, n.burnin = n.burnin)

    return(model.fit)
}

#' Apply HiLDA to statistically testing the burdens of mutation signatures between two groups
#'
#' @param inputG a MutationFeatureData S4 class output by the pmsignature.
#' @param inputParam a EstimatedParameters S4 class output by the pmsignature.
#' @param refGroup the reference group
#' @param sigOrder the indice indicating the samples in the reference group.
#' @param n.iter number of total iterations per chain (default: 2000).
#' @param n.burnin length of burn (default: 0).
#'
#'
#' @importFrom R2jags jags
#' @export


hilda_test <- function(inputG = G, inputParam = Param, refGroup, sigOrder = NULL, n.iter = 2000,
    n.burnin = 0, ...) {
    n.sig <- inputParam@signatureNum

    if (is.null(sigOrder)) {
        sigOrder <- 1:n.sig
    }

    # reshape the counts of input data
    countwide <- as.data.frame(t(inputG@countData))
    colnames(countwide) <- c("type", "sample", "count")
    countlong <- reshape(countwide, idvar = "sample", timevar = "type", direction = "wide")
    countlong <- countlong[order(countlong[, 1]), ]
    countlong[is.na(countlong)] <- 0
    countlong <- countlong[, -1]

    # generate the known data for MCMC
    Num1 <- length(refGroup)
    Num2 <- length(inputG@sampleList) - Num1
    rownames(countlong) <- 1:(Num1 + Num2)
    mutationN <- as.vector(rowSums(countlong))
    caseGroup <- setdiff(1:(Num1 + Num2), refGroup)
    numBases <- length(inputG@possibleFeatures)

    X_G1 <- structure(.Data = rep(0, Num1 * max(mutationN[refGroup]) * numBases), .Dim = c(max(mutationN[refGroup]),
        numBases, Num1))
    X_G2 <- structure(.Data = rep(0, Num2 * max(mutationN[caseGroup]) * numBases), .Dim = c(max(mutationN[caseGroup]),
        numBases, Num2))

    for (i in refGroup) {
        X_G1[1:sum(countlong[i, ]), , which(refGroup == i)] <- t(inputG@featureVectorList[, rep(1:ncol(inputG@featureVectorList),
            countlong[i, ])])
    }

    for (i in caseGroup) {
        X_G2[1:sum(countlong[i, ]), , which(caseGroup == i)] <- t(inputG@featureVectorList[,
            rep(1:ncol(inputG@featureVectorList), countlong[i, ])])
    }

    # set up the MCMC
    N1 <- mutationN[refGroup]
    N2 <- mutationN[caseGroup]

    jdata <- list(I1 = Num1, I2 = Num2, K = n.sig, numStates = 6, numflank = numBases - 1, N1 = N1,
        N2 = N2, X_T = X_G1, X_B = X_G2)

    inits <- hilda_inits(inputParam, refGroup, sigOrder)

    var.s <- c("p.states1", "p.states2", "p", "alpha", "beta")

    model.fit <- R2jags::jags(model.file = system.file("models/hilda.txt", package = "HiLDA"), data = jdata, parameters.to.save = var.s,
        inits = inits, n.chains = 2, n.iter = n.iter, n.burnin = n.burnin)

    return(model.fit)
}

#' Output the maximum potential scale reduction statistic of all parameters estimated
#'
#' @param jags.output the output jags file generated by the jags function from the R2jags package.
#'
#'
#' @export
#

hilda_rhat <- function(jags.output) max(jags.output$BUGSoutput$summary[, "Rhat"])


#' Extract the posterior distributions of the mean differences in muational exposures
#'
#' @param jags.output the output jags file generated by the jags function from the R2jags package.
#'
#'
#' @export
#

hilda_posterior <- function(jags.output) {
    rownames <- rownames(jags.output$BUGSoutput$summary)
    beta.index <- stringr::str_detect(rownames, "beta")
    return(jags.output$BUGSoutput$summary[beta.index, ])
}

#' Compute the Bayes factor
#'
#' @param jags.output the output jags file generated by the jags function from the R2jags package.
#' @param prob_M1 the probability of sampling the model 1 (no difference in the means), default is 0.5
#'
#' @export
#

hilda_bayesfactor_result <- function(jags.output, prob_M1 = 0.5) {
    freq <- table(jags.output$BUGSoutput$sims.list$pM2)
    if (length(freq) == 1) {
        stop(paste("It got stuck in the model", as.numeric(names(freq)) + 1))
    }
    return(as.vector(freq)[2]/as.vector(freq)[1])
}
