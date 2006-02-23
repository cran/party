
# $Id: Classes.R 2475 2006-02-23 15:11:21Z hothorn $

### Conditional Expectation and Covariance
setClass(Class = "ExpectCovar",
    representation = representation(
        expectation = "numeric",
        covariance  = "matrix",
        dimension   = "integer"
   )
)

### Expectation and Covariance of the influence function 
### (+ sum of weights)
setClass(Class = "ExpectCovarInfluence",
    representation = representation(
        sumweights = "numeric"
    ),
    contains = "ExpectCovar"
)

### Linear statistic with expectation and covariance
setClass(Class = "LinStatExpectCovar",
    representation = representation(
        linearstatistic = "numeric",
        expcovinf = "ExpectCovarInfluence"
    ),
    contains = "ExpectCovar"
)

### Memory for C_svd
setClass(Class = "svd_mem",
    representation = representation(
        svd    = "list",
        method = "character",
        jobu   = "character",
        jobv   = "character",
        u      = "matrix",
        v      = "matrix",
        s      = "numeric",
        p      = "integer"
    )
)

### with Moore-Penrose inverse of the covariance matrix
setClass(Class = "LinStatExpectCovarMPinv",
    representation = representation(
        MPinv  = "matrix",   
        rank   = "numeric",
        svdmem = "svd_mem"
    ), 
    contains = "LinStatExpectCovar"
)

################ Memory Classes #####################

setClass(Class = "TreeFitMemory",
    representation = representation(
        expcovinf         = "ExpectCovarInfluence",
        expcovinfss       = "ExpectCovarInfluence",
        linexpcov2sample  = "LinStatExpectCovar",
        weights           = "list",
        varmemory         = "list",
        varMmemory        = "list",
        Mscorematrices    = "list",
        dontuse           = "logical",
        dontusetmp        = "logical",
        splitstatistics   = "numeric"
    ), 
    validity = function(object) {
        ni <- length(dontuse)
        all(sapply(object@Mscorematrices, class) == "matrix") &&
        (length(varmemory) == ni && length(dontusetmp) == ni)
    }
)


##############  Tree Classes  ######################

setClassUnion("df_OR_list", c("data.frame", "list"))

setClass(Class = "VariableControl",
    representation = representation(
        teststattype = "factor",
        pvalue       = "logical",
        tol          = "numeric",
        maxpts       = "integer",
        abseps       = "numeric",
        releps       = "numeric"
    ),
    prototype = list(
        teststattype = factor("maxabs", levels = c("maxabs", "quadform")),
        pvalue       = TRUE,
        tol          = 1e-10,
        maxpts       = as.integer(25000),
        abseps       = 1e-4,
        releps       = 0.0
    )
)

setClass(Class = "SplitControl",
    representation = representation(
        minprob      = "numeric",
        minsplit     = "numeric",
        minbucket    = "numeric",
        tol          = "numeric",
        maxsurrogate = "integer"
    ),
    prototype = list(minprob = 0.01, minsplit = 20, 
                     minbucket = 7, tol = 1e-10, maxsurrogate = as.integer(0)
    )
)

setClass(Class = "GlobalTestControl",
    representation = representation(
        testtype     = "factor",
        nresample    = "integer",
        randomsplits = "logical",
        mtry         = "integer",
        mincriterion = "numeric"
    ),
    prototype = list(
        testtype = factor("Bonferroni", 
            levels = c("Bonferroni", "MonteCarlo", "Aggregated", 
                       "Univariate", "Teststatistic")),
        nresample = as.integer(9999),
        randomsplits = FALSE,
        mtry = as.integer(0),
        mincriterion = 0.95
    )
)

setClass(Class = "TreeGrowControl",
    representation = representation(
        stump           = "logical",
        maxdepth        = "integer",
        savesplitstats  = "logical"
    ),
    prototype = list(stump = FALSE, maxdepth = as.integer(0), savesplitstats = TRUE)
)

setClass(Class = "TreeControl",
    representation = representation(
        varctrl   = "VariableControl",
        splitctrl = "SplitControl",
        gtctrl    = "GlobalTestControl",
        tgctrl    = "TreeGrowControl"
    ),
    prototype = list(varctrl = new("VariableControl"),
                     splitctrl = new("SplitControl"),
                     gtctrl = new("GlobalTestControl"),
                     tgctrl = new("TreeGrowControl")
    )
)

setClass(Class = "VariableFrame",
    representation = representation(
        variables       = "df_OR_list", 
        transformations = "list", 
        is_nominal      = "logical", 
        is_ordinal      = "logical",
        is_censored     = "logical",
        ordering        = "list", 
        levels          = "list", 
        scores          = "list",
        has_missings    = "logical", 
        whichNA         = "list",
        jointtransf     = "matrix",
        nobs            = "integer",
        ninputs         = "integer")
)

setClass(Class = "LearningSample",
    representation = representation(
        responses = "VariableFrame",
        inputs    = "VariableFrame",
        weights   = "numeric",
        nobs      = "integer",
        ninputs   = "integer",
        menv      = "ModelEnv"
    )
)

### the tree structure itself is a list, 
### and we need to make sure that the tree slot excepts
### the S3 classes. 
setClass(Class = "SplittingNode", contains = "list")
setClass(Class = "TerminalNode", contains = "list")
setClass(Class = "orderedSplit", contains = "list")
setClass(Class = "nominalSplit", contains = "list")

### and we don't want to see warnings that class `Surv'
### (S3 method in `survival') is unknown
setClass(Class = "Surv", contains = "list")


### A class for partitions induced by recursive binary splits
setClass(Class = "BinaryTreePartition",
    representation = representation(
        tree     = "list",          # the basic tree structure as (named or
                                    # unnamed) list
        where    = "integer"        # the nodeID of the observations in the
                                    # learning sample
    ),
)

### A class for binary trees   
setClass(Class = "BinaryTree", 
    representation = representation(
        data                = "ModelEnv",
        responses           = "VariableFrame", # a list of response `variables'
                                               # for computing predictions
        cond_distr_response = "function",      # predict distribtion
        predict_response    = "function",      # predict responses
        prediction_weights  = "function",      # prediction weights
        get_where           = "function"       # node numbers
    ),
    contains = "BinaryTreePartition"
)

### A class for random forest  
setClass(Class = "RandomForest", 
    representation = representation(
        ensemble            = "list",
        data                = "ModelEnv",
        responses           = "VariableFrame", # a list of response `variables'
                                               # for computing predictions
        cond_distr_response = "function",      # predict distribtion
        predict_response    = "function",      # predict responses
        prediction_weights  = "function"      # prediction weights
    )
)

