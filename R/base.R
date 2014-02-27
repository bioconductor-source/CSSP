#'An S-4 class containing the model fit information for CSSP model.
#'
#'\describe{
#'\item{lambdax}{Sequencing depth of the input sample.}
#'\item{lambday}{Sequencing depth of the ChIP sample.}
#'\item{e0}{The normalization parameter for the ChIP sample.}
#'\item{pi0}{The pi_0 parameter of CSSP model, denoting the proportion of bins that are enriched.}
#'\item{mu.chip}{The vector of the estimated hyper means for the background model of the ChIP sample.}
#'\item{mu.input}{The vector of the estimated hyper means for the input sample.}
#'\item{mean.sig}{The vector of the hyper means for each signal component.}
#'\item{size.sig}{The vector of the size parameters for each signal component.}
#'\item{a}{The size parameter of the input sample model.}
#'\item{b}{The size parameter of the background model for the ChIP sample.}
#'\item{p.sig}{The vector of the proportions of enrichment as each signal component across all enrichment bins.}
#'\item{prob.zero}{The vector of the prior inflated probability at 0.}
#'\item{post.p.sig}{The matrix for the posterior probability of each bin being enriched as a signal component conditioning on the event that the bin is enriched. Each column corresponds to one signal component.}
#'\item{post.p.bind}{Posterior probability of each bin being enriched.}
#'\item{post.p.zero}{Posterior probability of the inflated probability at 0.}
#'\item{post.shape.sig}{The matrix for the shape parameters for the posterior gamma distributions of bin level poisson parameters, conditioning on the event that the bins are enriched as each signal component. Each column corresponds to one signal component.}
#'\item{post.scale.sig}{The matrix for the scale parameters of the posterior gamma distributions of bin level poisson parameters, conditioning on the event that the bins are enriched as each signal component. Each column corresponds to one signal component.}
#'\item{post.shape.back}{The shape parameters for the posterior gamma distributions of bin level poisson parameters, conditioning on each bin being enriched.}
#'\item{post.scale.back}{The scale parameters for the posterior gamma distributions of bin level poisson parameters, conditioning on each bin being unenriched.}
#'\item{n}{The number of mappable bins that are fitted by the model.}
#'\item{k}{The number of signal components.}
#'\item{map.id}{The indices for the mappable bins that are fitted by the model.}
#'\item{pvalue}{The continuously corrected p-values for a subset of ChIP sample bin counts against the background model.}
#'\item{cum.pval}{The cumulative distribution for p-values for a subset of ChIP sample bin counts against the background model.}
#'}
#'
#'@examples showClass("CSSPFit")
#'@name CSSPFit-class
#'@rdname CSSPFit-class
#'@exportClass CSSPFit
setClass("CSSPFit",
         representation=representation(
           lambdax="numeric",
           lambday="numeric",
           e0="numeric",
           pi0="numeric",
           a="numeric",
           b="numeric",
           mu.chip="numeric",
           mu.input="numeric",
           mean.sig="numeric",
           size.sig="numeric",
           p.sig="numeric",
           prob.zero="numeric",
           post.p.sig="matrix",
           post.p.bind="numeric",
           post.p.zero="numeric",
           post.shape.sig="matrix",
           post.shape.back="numeric",
           post.scale.sig="matrix",
           post.scale.back="numeric",
           n="numeric",
           k="numeric",
           map.id="numeric",
           pvalue="numeric",
           cum.pval="numeric"
           )
         )

#'An S-4 class containing the model fit information for a CSSP model.
#'
#'\describe{
#'\item{chrID}{The chromosome ID.}
#'\item{coord}{The genome coordinates for the starting positions of each bin.}
#'\item{tagCount}{The number of ChIP reads mapped to each bin.}
#'\item{mappability}{The mappability score of each bin.}
#'\item{gcContent}{The gc-content score of each bin.}
#'\item{input}{The number of input reads mapped to each bin.}
#'\item{dataType}{Either "unique" or "multi'.}
#'}
#'
#'@docType class
#'@name BinData-class
#'@rdname BinData-class
#'@exportClass BinData
setClass( Class="BinData",
         representation=representation(
           coord="numeric",
           tagCount="numeric",
           mappability="numeric",
           gcContent="numeric",
           input="numeric",
           dataType="character",
           chrID="character"
           )
         )

#' @name bindpos
#' @title An artificially constructed dataset containing enrichment positions on 5 chromosomes.
#' @description This data set contains artificially generated nucleotide-leve enrichment positions on a genome of 5 chromosomes.
#' @docType data
#' @usage example
#' @format A \link{list} containing the genome coordinates for enrichment sites on each of the 5 chromosomes.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name peakpos
#' @title An artificially generated dataset containing peak invervals on 5 chromosomes.
#' @description This data set contains the genome coordinates of artificially generated peak intervals on a genome of 5 chromosomes.
#' @docType data
#' @format a \link{list} of 2-column matrices. Each matrix contains the coordinates of the peak intervals for one chromosome.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name tagdat_chip
#' @title An artificially constructed dataset containing genome coordinates for aligned ChIP sample reads.
#' @description This dataset contains artificially generated genome coordinates for ChIP sample reads on a genome of 5 chromosomes. The sign of each read represents the strand direction, with 5' represented by positive numbers and 3' represented by negative numbers.
#' @docType data
#' @usage example
#' @format a \link{list} containing the reads coordinates on each of the 5 chromosomes.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name tagdat_input
#' @title An artificially constructed dataset containing genome coordinates for aligned input sample reads.
#' @description This dataset contains artificially generated genome coordinates for ChIP sample reads on a genome of 5 chromosomes. The sign of each read represents the strand direction, with 5' represented by positive numbers and 3' represented by negative numbers.
#' @docType data
#' @usage example
#' @format a \link{list} containing the reads coordinates on each of the 5 chromosomes.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name bin.data
#' @title An artificially constructed \link{BinData-class} class object.
#' @description This data set contains a typical example for a BinData class object,.
#' @docType data
#' @format a \link{BinData-class} class object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name bindata.chr1
#' @title An artificially constructed data.frame object that can be used by cssp.fit function.
#' @description This data set contains a typical example for a data.frame object that can be imported by \link{cssp.fit}.
#' @docType data
#' @format a \link{data.frame} class object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name sampleFit
#' @title A "CSSPFit" class object containing the fitted CSSP model for \link{bin.data}.
#' @description A \link{CSSPFit-class} class object constructed by fitting CSSP model on \link{bin.data}.
#' @docType data
#' @format a \link{CSSPFit-class} class object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

