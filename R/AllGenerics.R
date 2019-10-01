#Constructor
#' @export
#' @rdname SangerReadss-class
setGeneric("SangerReads", function(obj) standardGeneric("SangerReads"))


setMethod("SangerReads", "abif",
          function(obj) {
              res <- new("SangerReads")
          })


# setGeneric("forward", function(x) standardGeneric("forward"))
#
# setMethod("forward", signature("SangerReads"), function(object) {
#     3
# })

### sangerseqR package Try First
hetab1 <- read.abif(system.file("extdata", "heterozygous.ab1", package = "sangerseqR"))
str(hetab1, list.len = 20)

homoscf <- read.scf(system.file("extdata", "homozygous.scf", package = "sangerseqR"))
str(homoscf)
