#!/usr/bin/Rscript --vanilla
##
## Analyse the simulation outputs and produce summary files.
##
## It will read output files in `./data/model-output` and write summary files
## to `./data/model-analysis`.
##
##
## Usage:
##
##     ./analyse-output.R
##

SCRIPT <- 'analyse-output.R'


load_flu_data <- function() {
    data_file <- './data/influenza/influenza-case-distribution.ssv'
    df <- read.table(data_file, header = TRUE)
    df_med <- aggregate(pcnt_of_cases ~ SA3, data = df, FUN = median)
    names(df_med) <- c('SA3', 'median')
    df_sd <- aggregate(pcnt_of_cases ~ SA3, data = df, FUN = sd)
    names(df_sd) <- c('SA3', 'sd')
    df_stat <- merge(df_med, df_sd)
    df_popn <- unique(df[, c('SA3', 'pcnt_of_popn')])
    merge(df_stat, df_popn)
}


find_models <- function(path = NULL, pattern = NULL) {
    if (is.null(path)) {
        path <- './data/model-output'
    }
    if (is.null(pattern)) {
        pattern <- '.*-output.ssv'
    }
    list.files(path, pattern, full.names = TRUE)
}


find_stoch_models <- function(path = NULL) {
    find_models(path = path, pattern = '.*-output.ssv')
}


find_det_models <- function(path = NULL) {
    find_models(path = path, pattern = '.*-deterministic.ssv')
}


find_det_variable_frac_self_models <- function(path = NULL) {
    find_models(path = path, pattern = '.*-variable-frac-self.ssv')
}


calc_gaussian_error <- function(y, x_median, x_sd) {
    1 - dnorm(y, x_median, x_sd) / dnorm(x_median, x_median, x_sd)
}


calc_model_errors <- function(df, df_flu) {
    net_infs <- aggregate(cum_infs ~ sim, data = df, FUN = sum)
    names(net_infs) <- c('sim', 'final_size')
    df <- merge(df, net_infs)
    df$pcnt_of_cases <- 100 * df$cum_infs / df$final_size
    df <- merge(df, df_flu[, c('SA3', 'median', 'sd', 'pcnt_of_popn')],
                by.x = 'ix', by.y = 'SA3')

    df$mae <- abs(df$pcnt_of_cases - df$median)
    df$mae_baseline <- abs(df$pcnt_of_popn - df$median)
    df$mae_improved <- as.integer(df$mae < df$mae_baseline)

    df$gau <- calc_gaussian_error(df$pcnt_of_cases, df$median, df$sd)
    df$gau_baseline <- calc_gaussian_error(df$pcnt_of_popn, df$median, df$sd)
    df$gau_improved <- as.integer(df$gau < df$gau_baseline)

    df
}


calculate_errors <- function(df_flu, file_sa3, file_net) {
    if (file.exists(file_sa3) && file.exists(file_net)) {
        cat('Loading', file_sa3, '...\n')
        df_errors <- read.table(file_sa3, header = TRUE, as.is = TRUE)
        df_errors$SA3 <- factor(df_errors$SA3)
        cat('Loading', file_net, '...\n')
        df_net <- read.table(file_net, header = TRUE, as.is = TRUE)
        return(list(df_errors = df_errors, df_net = df_net))
    }

    df_errors <- NULL
    model_files <- find_det_models()
    n_files <- length(model_files)
    n_width <- nchar(n_files)
    i <- 0
    for (model_file in model_files) {
        i <- i + 1
        cat('Reading file ', formatC(i, width = n_width),
            ' of ', formatC(n_files, width = n_width),
            ': ', model_file, '...\n', sep = '')
        df <- read.table(model_file, header = TRUE)
        f_self_vals <- sort(unique(df$frac_self))
        f_cbd_vals <- sort(unique(df$frac_cbd))
        for (f_cbd in f_cbd_vals) {
            mask_cbd <-  df$frac_cbd == f_cbd
            for (f_self in f_self_vals) {
                mask <- mask_cbd & df$frac_self == f_self
                tbl <- df[mask, ]
                ## Calculate the MAEs and Gaussian errors.
                tbl_errors <- calc_model_errors(tbl, df_flu)
                tbl_errors$file <- sub('-deterministic.ssv', '',
                                       basename(model_file))
                tbl_errors$frac_cbd <- unique(tbl$frac_cbd)
                tbl_errors$frac_self <- unique(tbl$frac_self)
                df_errors <- rbind(df_errors, tbl_errors)
            }
        }
    }

    model_files <- find_det_variable_frac_self_models()
    n_files <- length(model_files)
    n_width <- nchar(n_files)
    i <- 0
    for (model_file in model_files) {
        i <- i + 1
        cat('Reading file ', formatC(i, width = n_width),
            ' of ', formatC(n_files, width = n_width),
            ': ', model_file, '...\n', sep = '')
        df <- read.table(model_file, header = TRUE)
        f_self_vals <- sort(unique(df$frac_self))
        f_cbd_vals <- sort(unique(df$frac_cbd))
        for (f_cbd in f_cbd_vals) {
            mask_cbd <-  df$frac_cbd == f_cbd
            for (f_self in f_self_vals) {
                mask <- mask_cbd & df$frac_self == f_self
                tbl <- df[mask, ]
                ## Calculate the MAEs and Gaussian errors.
                tbl_errors <- calc_model_errors(tbl, df_flu)
                tbl_errors$file <- sub('-variable-frac-self.ssv', '',
                                       basename(model_file))
                tbl_errors$frac_cbd <- unique(tbl$frac_cbd)
                tbl_errors$frac_self <- unique(tbl$frac_self)
                df_errors <- rbind(df_errors, tbl_errors)
            }
        }
    }

    cat('Aggregating errors ...\n')
    df_errors$cases <- df_errors$cum_infs
    df_errors$SA3 <- factor(df_errors$ix)
    df_errors <- df_errors[, c('file', 'frac_self', 'frac_cbd', 'SA3',
                               'cases', 'pcnt_of_cases', 'median', 'sd',
                               'popn', 'pcnt_of_popn',
                               'mae', 'mae_baseline', 'mae_improved',
                               'gau', 'gau_baseline', 'gau_improved')]

    ## Calculate flat sums of the MAE and Gaussian errors.
    df_net_mae <- aggregate(mae ~ frac_cbd + frac_self + file,
                            data = df_errors, FUN = sum)
    df_net_gau <- aggregate(gau ~ frac_cbd + frac_self + file,
                            data = df_errors, FUN = sum)

    ## Calculate population-weighted sums of the MAE and Gaussian errors.
    ## No need to store these values in df_errors, they only make sense
    ## when aggregated over the SA3s.
    df_errors$mae_ppn <- df_errors$mae * df_errors$pcnt_of_popn / 100
    df_errors$gau_ppn <- df_errors$gau * df_errors$pcnt_of_popn / 100
    df_ppn_mae <- aggregate(mae_ppn ~ frac_cbd + frac_self + file,
                            data = df_errors, FUN = sum)
    df_ppn_gau <- aggregate(gau_ppn ~ frac_cbd + frac_self + file,
                            data = df_errors, FUN = sum)
    df_errors$mae_ppn <- NULL
    df_errors$gau_ppn <- NULL

    df_net <- merge(df_net_mae, df_net_gau)
    df_net <- merge(df_net, df_ppn_mae)
    df_net <- merge(df_net, df_ppn_gau)

    ## Baseline errors (non-spatial model).
    net_mae_baseline <- sum(unique(df_errors[, c('SA3', 'mae_baseline')])
                            $mae_baseline)
    net_gau_baseline <- sum(unique(df_errors[, c('SA3', 'gau_baseline')])
                            $gau_baseline)
    df_net$mae_baseline <- net_mae_baseline
    df_net$mae_improved <- as.integer(df_net$mae < net_mae_baseline)
    df_net$gau_baseline <- net_gau_baseline
    df_net$gau_improved <- as.integer(df_net$gau < net_gau_baseline)

    ## Population-weighted baseline errors.
    mae_ppn_df <- unique(df_errors[, c('SA3', 'mae_baseline', 'pcnt_of_popn')])
    mae_ppn_df$mae_ppn <- mae_ppn_df$mae * mae_ppn_df$pcnt_of_popn / 100
    mae_ppn_baseline <- sum(mae_ppn_df$mae_ppn)
    gau_ppn_df <- unique(df_errors[, c('SA3', 'gau_baseline', 'pcnt_of_popn')])
    gau_ppn_df$gau_ppn <- gau_ppn_df$gau * gau_ppn_df$pcnt_of_popn / 100
    gau_ppn_baseline <- sum(gau_ppn_df$gau_ppn)
    df_net$mae_ppn_baseline <- mae_ppn_baseline
    df_net$mae_ppn_improved <- as.integer(df_net$mae_ppn < mae_ppn_baseline)
    df_net$gau_ppn_baseline <- gau_ppn_baseline
    df_net$gau_ppn_improved <- as.integer(df_net$gau_ppn < gau_ppn_baseline)

    ## Distinguish between ABS and Sygic data sets.
    df_errors$data_src <- 'ABS'
    df_errors$data_src[grepl("^sygic", df_errors$file)] <- 'Sygic'
    df_net$data_src <- 'ABS'
    df_net$data_src[grepl("^sygic", df_net$file)] <- 'Sygic'

    cat('Writing', file_sa3, '...\n')
    write.table(df_errors, file = file_sa3, row.names = FALSE, quote = FALSE)
    cat('Writing', file_net, '...\n')
    write.table(df_net, file = file_net, row.names = FALSE, quote = FALSE)

    return(list(df_errors = df_errors, df_net = df_net))
}


save_sa3_min_errors <- function(error_dfs, path_out) {
    df <- error_dfs$df_errors

    ## Restrict to two OD matrices.
    keep_files <- c('abs_all', 'abs_private', 'abs_public', 'sygic')
    df <- df[df$file %in% keep_files, ]

    ## Extract the smallest (in magnitude) error for each SA3.
    df$err <- df$pcnt_of_cases - df$median
    ## https://stackoverflow.com/a/40352069 ... handle duplicates???
    absmin <- function(x) { x[which.min( abs(x) )]}
    df <- aggregate(err ~ file + SA3, data = df, FUN = absmin)

    ## Only retain pertinent columns.
    keep_cols <- c('file', 'SA3', 'err')
    df <- df[, keep_cols]

    labels <- c('abs_all' = 'ABS (all modes)',
                'abs_private' = 'ABS (private)',
                'abs_public' = 'ABS (public)',
                'sygic' = 'Sygic')
    df$label <- labels[df$file]

    ## Write the entire table to a single file.
    out_file <- file.path(path_out, 'sa3-min-err.ssv')
    cat('Writing', out_file, '...\n')
    write.table(df, file = out_file, row.names = FALSE, quote = TRUE)

    ## For each OD matrix, write these smallest errors to disk.
    for (od_file in keep_files) {
        out_base <- paste0('sa3-min-err-', od_file, '.ssv')
        out_file <- file.path(path_out, out_base)

        cat('Writing', out_file, '...\n')
        write.table(df[df$file == od_file, keep_cols], file = out_file,
                    row.names = FALSE, quote = TRUE)
    }
}


save_sa3_deltas <- function(error_dfs, path_out) {
    df <- error_dfs$df_errors

    ## Restrict to a single value of frac_self, we only want to see the
    ## relative magnitudes of the change in each SA3.
    df <- df[df$frac_self == min(df$frac_self), ]

    ## Restrict to a single value of frac_cbd, it has only a minor effect on
    ## non-CBD SA3s.
    df <- df[df$frac_cbd == 0.33, ]

    ## Restrict to two OD matrices.
    keep_files <- c('abs_all', 'abs_private', 'abs_public', 'sygic')
    df <- df[df$file %in% keep_files, ]

    ## Calculate the change in burden, relative to the non-spatial model.
    df$delta <- df$pcnt_of_cases - df$pcnt_of_popn

    ## Convert to a discrete scale: increase or decrease, small or large.
    df$direction <- 'Increase'
    df$direction[df$delta < 0] <- 'Decrease'
    fixed_toln <- 0.1
    df$direction[abs(df$delta) <= fixed_toln] <- 'No change'
    df$magnitude <- '(< 1%)'
    df$magnitude[abs(df$delta) >= 1] <- '(> 1%)'
    df$magnitude[abs(df$delta) <= fixed_toln] <- '(< 0.1%)'
    df$delta <- paste(df$direction, df$magnitude)

    ## Only retain pertinent columns.
    keep_cols <- c('file', 'frac_self', 'frac_cbd', 'SA3', 'delta')
    df <- df[, keep_cols]

    labels <- c('abs_all' = 'ABS (all modes)',
                'abs_private' = 'ABS (private)',
                'abs_public' = 'ABS (public)',
                'sygic' = 'Sygic')
    df$label <- labels[df$file]

    ## Write the entire table to a single file.
    out_file <- file.path(path_out, 'sa3-deltas.ssv')
    cat('Writing', out_file, '...\n')
    write.table(df, file = out_file, row.names = FALSE, quote = TRUE)

    ## For each OD matrix, write these subsets to disk.
    for (od_file in keep_files) {
        out_base <- paste0('sa3-deltas-', od_file, '.ssv')
        out_file <- file.path(path_out, out_base)

        cat('Writing', out_file, '...\n')
        write.table(df[df$file == od_file, keep_cols], file = out_file,
                    row.names = FALSE, quote = TRUE)
    }
}


median_errors <- function(error_dfs, file_med, file_cls, file_sgn) {
    df_errors <- error_dfs$df_errors

    ## Calculate the median errors in each SA3, grouped by the value of the
    ## `frac_self` parameter (i.e., aggregating over OD files and `frac_cbd`).
    df_median <- aggregate(pcnt_of_cases ~ SA3 + frac_self, data = df_errors,
                           FUN = median)
    merge_cols <- c('SA3', 'median', 'sd', 'pcnt_of_popn')
    df_median <- merge(df_median, unique(df_errors[, merge_cols]))
    ## Calculate the change in burden, relative to the non-spatial model.
    df_median$delta <- df_median$pcnt_of_cases - df_median$pcnt_of_popn
    cat('Writing', file_med, '...\n')
    write.table(df_median, file = file_med,
                row.names = FALSE, quote = TRUE)

    ## Calculate the median errors in each SA3, grouped by the value of the
    ## `frac_cbd` parameter (i.e., aggregating over OD files and `frac_self`).
    file_cbd_med <- sub('\\.[^.]+$', '-cbd.ssv', file_med)
    df_cbd_med <- aggregate(pcnt_of_cases ~ SA3 + frac_cbd, data = df_errors,
                            FUN = median)
    merge_cols <- c('SA3', 'median', 'sd', 'pcnt_of_popn')
    df_cbd_med <- merge(df_cbd_med, unique(df_errors[, merge_cols]))
    ## Calculate the change in burden, relative to the non-spatial model.
    df_cbd_med$delta <- df_cbd_med$pcnt_of_cases - df_cbd_med$pcnt_of_popn
    cat('Writing', file_cbd_med, '...\n')
    write.table(df_cbd_med, file = file_cbd_med,
                row.names = FALSE, quote = TRUE)

    ## Classify the errors in each SA3.
    sa3_list <- sort(unique(df_median$SA3))
    sa3_type <- c()
    sa3_fself <- c()
    sa3_error <- c()
    sa3_bline <- c()
    sa3_better <- c()
    sa3_worse <- c()
    for (sa3 in sa3_list) {
        dfss <- df_median[df_median$SA3 == sa3, ]
        dfss <- dfss[order(dfss$frac_self), ]
        dfss$error <- dfss$pcnt_of_cases - dfss$median
        dfss$baseline <- dfss$pcnt_of_popn - dfss$median
        max_fs <- dfss[dfss$frac_self == max(dfss$frac_self), ]
        min_fs <- dfss[dfss$frac_self == min(dfss$frac_self), ]
        best_fs <- dfss[abs(dfss$error) == min(abs(dfss$error)), ]
        n <- nrow(dfss)
        n_worse <- sum(abs(dfss$error) > abs(dfss$baseline))
        n_better <- sum(abs(dfss$error) < abs(dfss$baseline))
        n_worse_pos <- sum((abs(dfss$error) > abs(dfss$baseline)) & dfss$error > 0)
        n_worse_neg <- sum((abs(dfss$error) > abs(dfss$baseline)) & dfss$error < 0)
        if (abs(n_worse - n_better) <= 2) {
            ## Identify SA3s where the errors are evenly split between being
            ## better or worse than the non spatial model (50%, +/- 1 SA3).
            sa3_type <- c(sa3_type, ' Even')
        } else if (sign(max_fs$error) != sign(min_fs$error)) {
            ## Identify SA3s where the errors range from negative to positive.
            sa3_type <- c(sa3_type, ' Spans +/-')
        }
        else if (n_worse > n_better + 2) {
            ## Identify SA3s where the errors are mostly/always worse than the
            ## non-spatial model.
            ## Also identify if the errors are exclusively under-estimating or
            ## over-estimating the burden of disease.
            if (n_worse == n_worse_neg) {
                sa3_type <- c(sa3_type, 'Worse, under-est')
            } else if (n_worse == n_worse_pos) {
                sa3_type <- c(sa3_type, 'Worse, over-est')
            } else {
                sa3_type <- c(sa3_type, 'Worse')
            }
        } else if (n_better > n_worse + 2) {
            ## Identify SA3s where the errors are mostly/always better than
            ## the non-spatial model.
            ## Also identify if the smallest errors are obtained for small or
            ## large values of the `frac_self` parameter.
            if (best_fs$frac_self < 0.5) {
                sa3_type <- c(sa3_type, 'Better, small')
            } else {
                sa3_type <- c(sa3_type, 'Better, large')
            }
        } else {
            ## Identify SA3s where the errors are evenly split between being
            ## better or worse than the non spatial model (50%, +/- 1 SA3).
            sa3_type <- c(sa3_type, 'Even')
        }
        sa3_fself <- c(sa3_fself, best_fs$frac_self)
        sa3_error <- c(sa3_error, best_fs$error)
        sa3_bline <- c(sa3_bline, best_fs$baseline)
        sa3_better <- c(sa3_better, n_better)
        sa3_worse <- c(sa3_worse, n_worse)
    }

    ## Assemble these classifications into a single data frame.
    df_classify <- data.frame(SA3 = sa3_list, scale = sa3_type,
                              frac_self = sa3_fself,
                              error = sa3_error, baseline = sa3_bline,
                              better = sa3_better, worse = sa3_worse)
    cat('Writing', file_cls, '...\n')
    write.table(df_classify, file = file_cls, row.names = FALSE, quote = TRUE)

    ## Calculate the mean error in each SA3 (i.e., the mean of the medians).
    df_sign <- aggregate(pcnt_of_cases - median ~ SA3, FUN = mean,
                         data = df_median)
    names(df_sign)[names(df_sign) != 'SA3'] <- 'scale'
    cat('Writing', file_sgn, '...\n')
    write.table(df_sign, file = file_sgn, row.names = FALSE, quote = TRUE)

    return(list(df_median = df_median, df_cbd_med = df_cbd_med,
                df_classify = df_classify, df_sign = df_sign))
}


analyse <- function() {
    df_flu <- load_flu_data()

    ## For each SA3, calculate the mean absolute error and the Gaussian error.
    file_sa3 <- './data/model-analysis/model-errors-sa3.ssv'
    ## Aggregate these errors as flat sums and as population-weighted sums.
    file_net <- './data/model-analysis/model-errors-net.ssv'
    ## Median errors for each SA3 and each value of frac_self.
    file_med <- './data/model-analysis/median-errors.ssv'
    ## Classify median errors, relative to the non-spatial model.
    file_cls <- './data/model-analysis/classify-sa3-errors.ssv'
    ## Classify median errors, relative to the influenza case data.
    file_sgn <- './data/model-analysis/classify-sa3-error-sign.ssv'
    ## Record the change in the distribution of infections over the SA3s.
    path_dta <- './data/model-analysis'

    ## Calculate errors for each mixing matrix; i.e., for each combination of:
    ## origin-destination file, `frac_self`, and `frac_cbd`.
    error_dfs <- calculate_errors(df_flu, file_sa3, file_net)

    ## Calculate the change in the spatial distribution of infection for
    ## several mixing matrices. and save them to disk.
    save_sa3_deltas(error_dfs, path_dta)

    save_sa3_min_errors(error_dfs, path_dta)

    ## Calculate median errors (and other summary statistics) for each SA3.
    median_dfs <- median_errors(error_dfs, file_med, file_cls, file_sgn)
}


main <- function(args = NULL) {
    analyse()
    invisible(0)
}


##
## Only call main() if this script is being run from the command-line.
##
if (! interactive()) {
    file.arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
    file.path <- substring(file.arg, 8)
    file.regex <- paste0('^(|.*/)', SCRIPT)
    if (grepl(file.regex, file.path)) {
        quit(status = main(args = commandArgs(TRUE)))
    }
}
