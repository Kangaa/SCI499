#!/usr/bin/Rscript --vanilla
##
## Plot the changes in disease burden for each SA3.
##
##
## Usage:
##
##     ./plot-sa3-deltas.R input.ssv output.pdf
##


SCRIPT <- 'plot-sa3-deltas.R'


plot <- function(data_file, plot_file) {
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

    df$SA3 <- factor(df$SA3)

    ## Add points for the non-spatial model.
    df_1 <- df[df$frac_self == max(df$frac_self), ]
    df_1$frac_self <- 1
    df_1$pcnt_of_cases <- df_1$pcnt_of_popn
    df <- rbind(df, df_1)

    labels <- c('abs_all' = 'Journeys to work (all)',
                'abs_private' = 'Journeys to work (private)',
                'abs_public' = 'Journeys to work (public)',
                'sygic' = 'GPS data')
    df$modality <- labels[df$file]

    df_var_frac_self <- df[df$frac_self == 0.0, ]
    df <- df[df$frac_self != 0.0, ]

    p <- ggplot(mapping = aes(x = SA3_label, y = pcnt_of_cases - pcnt_of_popn,
                              colour = frac_self)) +
        geom_hline(yintercept = 0) +
        geom_point() +
        xlab(NULL) +
        ylab('Change in model infections (%)') +
        scale_colour_distiller(expression("\u03b4"[i]^H),
                               palette = 'Spectral',
                               limits = c(0.1, 1), breaks = c(0.1, 0.5, 1),
                               guide = guide_colourbar(direction = 'horizontal',
                                                       title.vjust = 1)) +
        scale_x_discrete(labels = parse(text = levels(df$SA3_label))) +
        facet_wrap(~ modality, scales = 'free_y') +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.position = 'top',
              legend.title = element_text(size = rel(1.25)),
              legend.text = element_text(size = rel(1.25)),
              strip.background = element_rect(size = 0, fill = 'transparent'),
              strip.text = element_text(colour = 'black', size = rel(1.5)),
              strip.text.x = element_text(margin = margin(b = 5)),
              strip.text.y = element_text(margin = margin(l = 10)),
              panel.spacing.x = unit(2, 'lines'),
              panel.spacing.y = unit(1, 'lines'),
              axis.title.y = element_text(size = rel(1.5)),
              plot.title = element_text(hjust = 0.5))

    p1 <- p %+% df[df$frac_cbd == 0.33, ]
    p2 <- p %+% df[df$frac_cbd == 0.67, ]

    p1a <- p1 +
        geom_point(data = df,
                   mapping = aes(x = SA3_label, y = median - pcnt_of_popn),
                   colour = 'black', shape = 4) +
        geom_point(data = df_var_frac_self[df_var_frac_self$frac_cbd == 0.33, ],
                   mapping = aes(x = SA3_label, y = pcnt_of_cases - pcnt_of_popn),
                   colour = 'black', shape = 0, size = 2)

    p2a <- p2 +
        geom_point(data = df,
                   mapping = aes(x = SA3_label, y = median - pcnt_of_popn),
                   colour = 'black', shape = 4) +
        geom_point(data = df_var_frac_self[df_var_frac_self$frac_cbd == 0.67, ],
                   mapping = aes(x = SA3_label, y = pcnt_of_cases - pcnt_of_popn),
                   colour = 'black', shape = 0, size = 2)

    ## Calculate the mean error for each combination of frac_self and frac_cbd.
    df_net <- rbind(df_var_frac_self, df)
    df_net$error <- df_net$pcnt_of_cases - df_net$median
    df_net$abs_error <- abs(df_net$error)
    df_err <- aggregate(abs_error ~ modality + frac_self + frac_cbd,
                        data = df_net, FUN = sum)
    df_err_var <- df_err[df_err$frac_self == 0, ]
    df_err <- df_err[df_err$frac_self > 0, ]

    p_error <- ggplot(df_err, aes(frac_self, abs_error,
                                  colour = factor(frac_cbd))) +
        geom_line() +
        geom_point(data = df_err_var,
                   mapping = aes(y = abs_error, colour = factor(frac_cbd)),
                   x = 1,
                   shape = 0, size = 3) +
        facet_wrap(~ modality, scale = 'free_y') +
        expand_limits(y = 0) +
        xlab(expression("\u03b4"[i]^H)) +
        ylab('Absolute error') +
        scale_colour_hue(expression("\u03b4"[C])) +
        theme_light() +
        theme(strip.background = element_rect(size = 0, fill = 'transparent'),
              strip.text = element_text(colour = 'black', size = rel(1.5)),
              strip.text.x = element_text(margin = margin(b = 5)),
              strip.text.y = element_text(margin = margin(l = 10)),
              panel.spacing.x = unit(2, 'lines'),
              panel.spacing.y = unit(1, 'lines'),
              axis.title.y = element_text(size = rel(1.5)),
              plot.title = element_text(hjust = 0.5))

    cat ('Writing plot to', plot_file, '...\n')
    CairoPDF(plot_file, bg = 'transparent', w = 10, h = 8)
    print(p1)
    print(p2)
    print(p1a)
    print(p2a)
    print(p_error)
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
