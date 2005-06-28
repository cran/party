
# $Id: Print.R,v 1.4 2005/06/28 15:40:16 hothorn Exp $

prettysplit <- function(x, inames = NULL, ilevels = NULL) {
    if (length(x) == 4)
        names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics")
    if (length(x) == 5)
        names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics",
                      "toleft")
    if (length(x) == 6)
        names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics",
                      "toleft", "table")
    if (x$ordered) {
        class(x) <- "orderedSplit"
    } else {
        class(x) <- "nominalSplit"
    }
    if (!is.null(ilevels)) {
        if (!is.null(ilevels[x[["variableID"]]]))
            attr(x$splitpoint, "levels") <- ilevels[[x[["variableID"]]]]
    }
    if (!is.null(inames)) x$variableName <- inames[x[["variableID"]]]
    return(x)
}

prettytree <- function(x, inames = NULL, ilevels = NULL) {
    names(x) <- c("nodeID", "weights", "criterion", "terminal",
                  "psplit", "ssplits", "prediction", "left", "right")
    if (is.null(inames) && extends(class(x), "BinaryTree"))
        inames <- x@inputnames
    names(x$criterion) <- c("criterion", "statistic", "maxcriterion")
    names(x$criterion$criterion) <- inames
    names(x$criterion$statistic) <- inames

    if (x$terminal) {
        class(x) <- "TerminalNode"
        return(x)
    }

    x$psplit <- prettysplit(x$psplit, inames = inames, ilevels = ilevels)
    if (length(x$ssplit) > 0)
        x$ssplit <- lapply(x$ssplit, prettysplit, inames = inames, 
                           ilevels = ilevels)

    class(x) <- "SplittingNode"
    x$left <- prettytree(x$left, inames = inames, ilevels = ilevels)   
    x$right <- prettytree(x$right, inames = inames, ilevels = ilevels)    
    return(x)
}
 
print.TerminalNode <- function(x, n = 1, ...) {
    cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ")* ", 
                    sep = "", collapse = ""),
        "weights =", sum(x$weights), "\n")
}
 
print.SplittingNode <- function(x, n = 1, ...) {
    cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ") ", sep=""))
    print(x$psplit, left = TRUE)
    cat(paste("; criterion = ", round(x$criterion[[3]], 3), 
              ", statistic = ", round(max(x$criterion[[1]]), 3), "\n", 
              collapse = "", sep = ""))
    print(x$left, n + 2)
    cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ") ", sep=""))
    print(x$psplit, left = FALSE)
    cat("\n")
    print(x$right, n + 2)
}

print.orderedSplit <- function(x, left = TRUE, ...) {
    if (!is.null(attr(x$splitpoint, "levels"))) {
        sp <- attr(x$splitpoint, "levels")[x$splitpoint]
    } else {
        sp <- x$splitpoint
    }
    if (!is.null(x$toleft)) left <- as.logical(x$toleft) == left
    if (left) {
        cat(x$variableName, "<=", sp)
    } else {
        cat(x$variableName, ">", sp)
    }
}

print.nominalSplit <- function(x, left = TRUE, ...) {

    levels <- attr(x$splitpoint, "levels")

    ### is > 0 for levels available in this node
    tab <- x$table

    if (left) {
        lev <- levels[as.logical(x$splitpoint) & (tab > 0)]
    } else {
        lev <- levels[!as.logical(x$splitpoint) & (tab > 0)]
    }

    txt <- paste("\{", paste(lev, collapse = ", "), "\}", collapse = "", sep = "")
    cat(x$variableName, "==", txt)
}


print.BinaryTreePartition <- function(x, ...)
    print(x@tree)

print.BinaryTree <- function(x, ...) {
    cat("\n")
    cat("\t Conditional tree with", length(unique(where(x))), 
        "terminal nodes\n\n")
    y <- x@responses
    if (y@is_censored) {
        cat("Response: ", names(y@variables), "(censored)\n")
    } else {
        if (y@ninputs > 1) {
            cat("Responses:", paste(names(y@variables), 
                                    collapse = ", "), "\n")
        }  else {
            cat("Response: ", names(y@variables), "\n")
        }
    }
    if (length(x@inputnames) > 1) {
        cat("Inputs: ", paste(x@inputnames, collapse = ", "), "\n")
    } else {
        cat("Input: ", x@inputnames, "\n")
    }
    cat("Number of observations: ", x@responses@nobs, "\n\n")
    print(x@tree)
}
