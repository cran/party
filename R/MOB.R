## MOdel-Based partitioning
## requirements for model objects
##   - StatModel object
##   - for fitting: estfun, reweight, weights, some extractor for objective function
##   - inherited if available: print, predict, coef, summary, residuals, logLik

## generic function mob creates objects of class "mob"
setClass("mob", contains = "BinaryTree")

mob <- function(formula, weights, data = list(),
  model = glinearModel, control = mob_control(), ...)
{
  if(inherits(formula, "formula")) {
    ## convenience preprocessor for formula like
    ## y ~ x + z | a + b
    mobpp <- function(formula, data, model) {
      ff <- attr(modeltools:::ParseFormula(formula), "formula")
      ff$input[[3]] <- ff$input[[2]]
      ff$input[[2]] <- ff$response[[2]]
      ### <FIXME> we can't handle missing values properly yet
      cc <- complete.cases(data)
      if (!all(cc)) {
          warning(sum(cc), " missing values in ", 
                  sQuote("data"), " have been removed")
          data <- data[cc,]
      }
      ### </FIXME>
      dpp(model, as.formula(ff$input), other = list(part = as.formula(ff$blocks)), data = data)
    }
    formula <- mobpp(formula, data, model)
  }

  ## fit global model
  if(missing(weights)) weights <- rep(1, dimension(formula, "part")[1])
  fm <- fit(model, formula, ...)


  ## the main recursion function
  mob_fit <- function(obj, mf, weights, control) {
    ### fit a model for the current node
    obj <- reweight(obj, weights)

    ### set up node (empty if reweighting failed)
    if(inherits(obj, "try-error")) {
      node <- list(nodeID = NULL, weights = weights,
                   criterion = list(statistic = 0, criterion = 0, maxcriterion = 0),
                   terminal = TRUE, prediction = 0, model = obj)
      class(node) <- "TerminalNodeModel"
      node$nodeID <- as.integer(nodeid)
      nodeid <<- nodeid + 1
      return(node)
    }
    thisnode <- mob_fit_setupnode(obj, mf, weights, control)
    thisnode$nodeID <- as.integer(nodeid)
    nodeid <<- nodeid + 1
    thisnode$model <- obj

    ### split (if appropriate)
    if(!thisnode$terminal) {
      ## compute size of (potential) children
      childweights <- mob_fit_childweights(thisnode, mf, weights)

      ## stop or...
      if(any(sapply(childweights, sum) == 0)) {
        thisnode$terminal <- TRUE
        class(thisnode) <- "TerminalModelNode"
        return(thisnode)
      }
      ## ...recall for children
      thisnode$left <- mob_fit(obj, mf, weights = childweights$left, control)
      thisnode$right <- mob_fit(obj, mf, weights = childweights$right, control)
    }
    
    return(thisnode)
  }

  ## recursive partitioning  
  nodeid <- 1
  tr <- mob_fit(fm, formula, weights = weights, control = control)

  ## package into return object
  rval <- new("mob", tree = tr, 
    responses = initVariableFrame(formula@get("response"), trafo = NULL),
    data = formula)
  return(rval)
}

## control splitting parameters
mob_control <- function(alpha = 0.05, bonferroni = TRUE, minsplit = 20, objfun = deviance,
  breakties = FALSE, verbose = FALSE)
{
  rval <- list(alpha = alpha, bonferroni = bonferroni, minsplit = minsplit,
               objfun = objfun, breakties = breakties, verbose = verbose)
  class(rval) <- "mob_control"
  return(rval)
}


## S3 fitted model functions

print.mob <- function(x, ...) print(x@tree)

print.TerminalModelNode <- function (x, n = 1, ...) {
    print.TerminalNode(x, n = n, ...)
    if (!is.null(x$model))
        cat("Terminal node model\n")
        print(x$model)
        cat("\n")
}

predict.mob <- function(object, newdata = NULL, type = c("response", "node"), ...)
{
    if(is.null(newdata)) {
      newpart <- object@data@get("part")
      newinput <- object@data@get("input")
    } else {
      if(inherits(newdata, "ModelEnvFormula")) {
        newpart <- newdata@get("part")
        newinput <- newdata@get("input")
      } else {
        newpart <- object@data@get("part", data = newdata)
        newinput <- object@data@get("input", data = newdata)
      }
    }
    nobs <- NROW(newpart)
    newpart <- initVariableFrame(newpart, trafo = NULL)
    
    nodeIDs <- .Call("R_get_nodeID", object@tree, newpart, as.double(0.0), PACKAGE = "party")

    type <- match.arg(type)
    if(type == "response") {
        pred <- vector(mode = "list", length = nobs)
        for (n in unique(nodeIDs)) {
            node <- .Call("R_get_nodebynum", object@tree, as.integer(n), PACKAGE = "party")
            indx <- which(nodeIDs == n)
            pred[indx] <- predict(node$model, newdata = newinput[indx,,drop = FALSE], ...)
        }
	rval <- if(isTRUE(all.equal(sapply(pred, length), rep(1, nobs)))) unlist(pred) else pred
    } else {
      rval <- nodeIDs
    }
    
    return(rval)
}

residuals.mob <- function(object, ...)
{
    newpart <- object@data@get("part")
    newinput <- object@data@get("input")
    nobs <- NROW(newpart)
    newpart <- initVariableFrame(newpart, trafo = NULL)
    nodeIDs <- .Call("R_get_nodeID", object@tree, newpart, as.double(0.0), PACKAGE = "party")

    res <- vector(mode = "list", length = nobs)
    for (n in unique(nodeIDs)) {
    	node <- .Call("R_get_nodebynum", object@tree, as.integer(n), PACKAGE = "party")
    	indx <- which(nodeIDs == n)
    	res[indx] <- residuals(node$model, ...)[indx]
    }
    if(isTRUE(all.equal(sapply(res, length), rep(1, nobs)))) res <- unlist(res)    
    return(res)
}

fitted.mob <- function(object, ...)
  predict(object, ...)

coef.mob <- function(object, node = NULL, ...) {
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- sapply(nodes(object, node), function(z) coef(z$model))
  if(!is.null(dim(rval))) {
    rval <- t(rval)
    rownames(rval) <- node
  }
  return(rval)
}

summary.mob <- function(object, node = NULL, ...) {
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- lapply(nodes(object, node), function(z) summary(z$model))
  if(length(rval) < 2) rval <- rval[[1]]
    else {
    names(rval) <- node 
  }
  return(rval)
}

sctest.mob <- function(x, node = NULL, ...) {
  if(is.null(node)) node <- 1:max(terminal_nodeIDs(x@tree))
  rval <- lapply(nodes(x, node),
    function(z) rbind(statistic = z$criterion$statistic, p.value = 1-z$criterion$criterion))
  if(length(rval) < 2) rval <- rval[[1]]
    else {
    names(rval) <- node 
  }
  rval
}

logLik.mob <- function(object, node = NULL, ...) {
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- lapply(nodes(object, node), function(z) logLik(z$model))
  rval <- structure(sum(sapply(rval, head, 1)),
    #not supported by glinearModel# nall = sum(sapply(rval, function(z) attr(z, "nall"))),
    ## nobs = sum(sapply(rval, function(z) attr(z, "nobs"))),
    df = sum(sapply(rval, function(z) attr(z, "df"))) + nterminal(object@tree) - 1,
    class = "logLik")
  return(rval)
}

deviance.mob <- function(object, node = NULL, ...) {
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- sum(sapply(nodes(object, node), function(z) deviance(z$model)))
  return(rval)
}

weights.mob <- function(object, node = NULL, ...) {
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- rowSums(sapply(nodes(object, node), function(z) weights(z$model)))
  return(rval)
}