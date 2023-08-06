#!/usr/bin/Rscript --vanilla
##
## Plot the distribution of influenza cases across the SA3s.
##
##
## Usage:
##
##     ./plot-flu-distr.R input.ssv output.pdf
##


SCRIPT <- 'plot-flu-distr.R'

load_case_distribution <- function(data_file) {
    df <- read.table(data_file, header = TRUE, as.is = TRUE)

    SA4_names <- c('206' = 'Inner', '207' = 'Inner East',
                   '208' = 'Inner South', '209' = 'North East',
                   '210' = 'North West', '211' = 'Outer East',
                   '212' = 'South East', '213' = 'West')
    df$SA4 <- substring(df$SA3, 1, 3)
    df$SA4_name <- SA4_names[df$SA4]
    df$SA3_num <- as.integer(substring(df$SA3, 4, 5))
    df$SA3_label <- paste0('"', df$SA4_name, ' "*', df$SA3_num)
    df$SA3_label[df$SA3 == 20604] <- 'bold(CBD)'
    df$SA3_label <- factor(df$SA3_label, ordered = TRUE,
                           levels = unique(df$SA3_label))

    df$cbd <- df$SA3 == '20604'
    df$SA3_by_popn <- factor(
        df$SA3_label, ordered = TRUE,
        levels = unique(df$SA3_label[order(df$popn)]))

    return(df)
}

colour_SA3s <- function(df) {
    ## Colour every 7th SA3, and highlight the CBD.
    subset <- df$year == min(df$year)
    ordering <- order(df$pcnt_of_popn[subset])
    sa3_order <- df[subset, ]$SA3_label[ordering]
    sa3_coloured <- sa3_order[1:5 * 7]
    df$coloured <- 'any'
    df$coloured[df$cbd] <- 'cbd'
    df$coloured[df$SA3_label %in% sa3_coloured] <- 'coloured'
    sa3_cbd <- df[subset, ]$cbd[ordering]
    sa3_clr <- df[subset, ]$coloured[ordering]
    df_lbl <- data.frame(
        x = seq(0.1, 5.4, length.out = length(sa3_order)),
        y = 9.2,
        cbd = sa3_cbd,
        coloured = sa3_clr,
        label = sa3_order)
    return(list(df = df, df_lbl = df_lbl))
}

## Plot the distribution of influenza cases in each epidemic against the SA3
## and the SA3 resident population; in both cases, order the SA3s by their
## resident population size.
plot <- function(data_file, plot_file) {
    dfs <- colour_SA3s(load_case_distribution(data_file))
    df <- dfs$df
    df_lbl <- dfs$df_lbl

    sa3_labels <- levels(df$SA3_label)

    ## Plot the proportion of influenza cases in each epidemic that occurred
    ## in each SA3 against the SA3s (ordered by resident population).
    p_distr1 <- ggplot(df, aes(SA3_by_popn, pcnt_of_cases, colour = cbd)) +
        geom_boxplot() +
        ## Avoid creating duplicate labels for each SA3.
        geom_text(data = df[! duplicated(df$SA3_by_popn), ],
                  mapping = aes(SA3_by_popn, 9.2, label = SA3_label, colour = cbd),
                  parse = TRUE,
                  angle = 270, hjust = 0) +
        scale_x_discrete(
            labels = c(),
            ## Draw a grid line at every third SA3.
            breaks = levels(df$SA3_by_popn)[
                seq(1, length(levels(df$SA3_by_popn)), 3)]) +
        xlab('SA3') +
        ylab('% of cases') +
        scale_y_continuous(limits = c(0, 9.3), expand = c(0, 0),
                           breaks = c(0, 2, 4, 6, 8)) +
        scale_colour_manual(breaks = c(FALSE, TRUE),
                            values = c('black', '#d95f02')) +
        guides(colour = 'none') +
        theme_light()

    ## Plot the proportion of influenza cases in each epidemic that occurred
    ## in each SA3 against the resident population.
    p_distr2 <- ggplot(df, aes(x = pcnt_of_popn, y = pcnt_of_cases,
                               group = SA3, colour = coloured)) +
        geom_abline(intercept = 0, slope = 1, colour = '#afafaf') +
        geom_boxplot() +
        geom_text(data = df_lbl,
                  mapping = aes(x, y, label = label,
                                group = 1, colour = coloured),
                  parse = TRUE,
                  angle = 270, hjust = 0) +
        xlab("% of population") +
        ylab("% of cases") +
        scale_x_continuous(limits = c(0, 5.5), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 9.3), expand = c(0, 0),
                           breaks = c(0, 2, 4, 6, 8)) +
        scale_colour_manual(breaks = c('any', 'cbd', 'coloured'),
                            values = c('black', '#d95f02', '#1b9e77')) +
        guides(colour = 'none') +
        theme_light()

    ## Plot the proportion of influenza cases in each epidemic that occurred
    ## in each SA3 against the SA3s (ordered numerically), and show the
    ## resident population distribution for comparison.
    p_distr3 <- ggplot(df, aes(SA3_label, pcnt_of_cases, group = factor(SA3))) +
        geom_boxplot() +
        geom_point(mapping = aes(SA3_label, pcnt_of_popn),
                   shape = 4, colour = 'red') +
        xlab(NULL) +
        ylab('% of cases') +
        scale_x_discrete(labels = parse(text = sa3_labels)) +
        scale_y_continuous(limits = c(0, 8), expand = c(0, 0)) +
        guides(colour = 'none') +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    ## Plot the difference between the resident population distribution and
    ## the influenza case distribution (i.e., the non-spatial model error).
    df_median <- aggregate(pcnt_of_cases ~ SA3_label, data = df, FUN = median)
    df_median <- merge(df_median, unique(df[, c('SA3_label', 'pcnt_of_popn')]))
    p_distr4 <- ggplot(df_median, aes(SA3_label, pcnt_of_popn - pcnt_of_cases,
                                      group = SA3_label)) +
        geom_bar(stat = "identity") +
        xlab(NULL) +
        ylab('% of population - % of cases') +
        scale_x_discrete(labels = parse(text = sa3_labels)) +
        guides(colour = 'none') +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    p_delta <- ggplot(df, aes(SA3_label, pcnt_of_cases - pcnt_of_popn)) +
        geom_hline(yintercept = 0) +
        geom_point(colour = '#7f7f7f', shape = 1) +
        geom_point(data = df_median, colour = 'black', size = 2.5) +
        # geom_point(data = df_median, colour = 'black') +
        # geom_line(data = df_median, colour = 'black',
        #           mapping = aes(factor(SA3), pcnt_of_cases - pcnt_of_popn),
        #           group = 1) +
        xlab(NULL) +
        ylab('Change in notified cases (%)') +
        scale_x_discrete(labels = parse(text = sa3_labels)) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    cat ('Writing plots to', plot_file, '...\n')
    CairoPDF(plot_file, bg = 'transparent', w = 8, h = 6)
    print(p_distr1)
    print(p_distr2)
    print(p_distr3)
    print(p_distr4)
    print(p_delta)
    invisible(dev.off())
}


load_libs <- function() {
    ## Ensure that the required libraries are available.
    reqd_libs <- c('ggplot2', 'Cairo')
    failed <- c()
    for (lib in reqd_libs) {
        if (! require(lib, character.only = TRUE, quietly = TRUE)) {
            failed <- c(failed, lib)
        }
    }
    if (length(failed) > 0) {
        cat('\nCould not load:', paste(failed, collapse = ', '), '\n\n')
        quit(status = 10)
    }
}


main <- function(args = NULL) {
    if (is.null(args) || length(args) != 2) {
        cat('\nUsage:', SCRIPT, 'input.ssv', 'output.pdf\n\n')
        return(invisible(2))
    }

    load_libs()
    plot(args[1], args[2])

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
