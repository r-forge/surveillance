################################################################################
### Get a tidy() table of parameters from an "hhh4" fit
### (not yet sure if this is useful as multiple components are non-standard)
###
### Copyright (C) 2025 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################

## this could be registered (lazily) as an S3 method for generics::tidy()
tidy.hhh4 <- function (x, effects = "fixed",
                       exponentiate = FALSE, # or TRUE or an index vector
                       conf.int = FALSE, conf.level = 0.95,
                       ...) # (currently unused)
{
    allEffects <- c("fixed", "ran_pars", "ran_vals") # as in 'broom.mixed'
    effects <- match.arg(effects, c(allEffects, "random"), # same as "ran_vals"
                         several.ok = TRUE)
    idx2Exp <- if (isFALSE(exponentiate)) NULL else exponentiate

    estSE <- coef(x, se = TRUE, idx2Exp = idx2Exp,
                  amplitudeShift = TRUE) # we always want this
    colnames(estSE) <- c("estimate", "std.error") # 'broom' convention

    result_fixed <-
        if ("fixed" %in% effects) {
            estSE_fixed <- head(estSE, x$dim[1L])
            cbind(.tidy.hhh4_fixed(x),
                  "exp" = startsWith(rownames(estSE_fixed), "exp("),
                  estSE_fixed)
        }

    result_random <-
        if (any(c("random", "ran_vals") %in% effects)) {
            estSE_random <- tail(estSE, x$dim[2L])
            cbind(.tidy.hhh4_random(x),
                  "exp" = startsWith(rownames(estSE_random), "exp("),
                  estSE_random)
        }

    result <- rbind(result_fixed, result_random)

    if (!is.null(result) && conf.int) {
        cis <- confint(x, level = conf.level,
                       idx2Exp = idx2Exp, amplitudeShift = TRUE)
        cipars <- intersect(rownames(result), rownames(cis))
        result[cipars, c("conf.low", "conf.high")] <- cis[cipars,,drop=FALSE]
    }

    if ("ran_pars" %in% effects) { # TODO
        warning("\"ran_pars\" is not implemented yet, but see 'x$Sigma'")
    }

    ## drop superfluous metadata
    rownames(result) <- NULL
    if (all(is.na(result$unit)))
        result$unit <- NULL
    if (isFALSE(exponentiate))
        result$exp <- NULL
    if (length(effects) == 1L)
        result$effect <- NULL
    result
}

## create index columns for the fixed effects
.tidy.hhh4_fixed <- function (object)
{
    model <- terms(object)
    comps <- c("ar","ne","end")[unlist(model$terms["offsetComp",])]
    names_theta <- names(model$initialTheta)
    names_terms <- unlist(model$terms["name",])
    names_terms[names_terms %in% c("1", "ri(iid)", "ri(car)")] <- "(Intercept)"
    idxAS <- surveillance:::getCoefIdxRenamed(names_theta, amplitudeShift = TRUE)$AS
    names_terms[idxAS] <- paste0(
        rep_len(c("amplitude", "shift"), length(idxAS)), ".",
        as.integer(sub(".*\\(([0-9]+).*", "\\1", names_terms[idxAS])) / 2
    )
    namesW <- names_theta[model$nFE + seq_len(model$nd)]
    result <- data.frame(
        "component" = c(comps, rep("W", model$nd),
                        rep("dispersion", model$nOverdisp)),
        "effect" = "fixed",
        "unit" = NA_character_,
        "term" = c(names_terms, sub("^neweights\\.", "", namesW),
                   rep("overdisp", model$nOverdisp)),
        row.names = NULL, check.names = FALSE
    )
    if (model$nOverdisp > 1)
        result$unit[result$term == "overdisp"] <- levels(model$indexPsi)
    result
}

## create index columns for the random intercepts
.tidy.hhh4_random <- function (object)
{
    names_ri <- tail(names(object$coefficients), object$dim[2L])
    m <- regexec("^([^.]+)\\.([^.]+)\\.(.+)$", names_ri)
    tab <- do.call(rbind, regmatches(names_ri, m))[,-1L]
    result <- data.frame(
        "component" = tab[,1L], # = compsFE[model$indexRE]
        "effect" = "random",
        "unit" = tab[,3L],
        "term" = tab[,2L],
        row.names = NULL, check.names = FALSE
    )
    result
}
