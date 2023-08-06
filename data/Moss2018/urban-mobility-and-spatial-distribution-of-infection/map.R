#!/usr/bin/Rscript --vanilla
##
## A set of convenience functions for plotting maps of metropolitan Melbourne
## and applying colour scales to each SA3.
##
##
## Interactive usage:
##
##     library(ggplot2)
##     library(viridis)
##     library(Cairo)
##
##     source('map.R')
##     map_data <- map.load_shapes()
##
##     sa3s <- c(20601, 20602, 20603, 20604, 20605, 20606, 20607, 20701,
##               20702, 20703, 20801, 20802, 20803, 20804, 20901, 20902,
##               20903, 20904, 21001, 21002, 21003, 21004, 21005, 21101,
##               21102, 21103, 21104, 21105, 21201, 21202, 21203, 21204,
##               21205, 21301, 21302, 21303, 21304, 21305)
##     df <- data.frame(SA3 = sa3s, scale = runif(length(sa3s)))
##
##     p <- map.plot(map_data, df) +
##         scale_fill_viridis("Random")
##
##     CairoPDF('output.pdf', w=6, h=4.9)
##     print(p)
##     dev.off()
##
##
## Command-line usage:
##
##     ./map.R input.ssv output.pdf
##


SCRIPT <- 'map.R'


map.load_libs <- function() {
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(RColorBrewer))
    suppressPackageStartupMessages(library(viridis))
    suppressPackageStartupMessages(library(Cairo))
}


map.load_shapes <- function(shape_dir = './data/shape-data') {
    file_vic <- file.path(shape_dir, 'victoria.ssv')
    file_mel <- file.path(shape_dir, 'melbourne.ssv')
    file_sa4 <- file.path(shape_dir, 'melbourne-sa4.ssv')
    vic.df <- read.table(file_vic, header = TRUE, as.is = TRUE)
    mel.df <- read.table(file_mel, header = TRUE, as.is = TRUE)
    sa4.df <- read.table(file_sa4, header = TRUE, as.is = TRUE)
    list(vic.df = vic.df, mel.df = mel.df, sa4.df = sa4.df)
}


map.plot <- function(map_data, sa3_scale, annotate = FALSE) {
    df <- merge(map_data$mel.df, sa3_scale, by.x = 'SA3_CODE', by.y = 'SA3')
    ## Draw SA4 borders in black on all maps except for the SA3 map, where we
    ## instead use a lighter colour so that the text annotation for each SA3
    ## stands out.
    if (annotate) {
        sa4_border <- '#5f5f5f'
    } else {
        sa4_border <- 'black'
    }
    p <- ggplot(df, aes(long, lat, group = group, fill = scale)) +
        geom_polygon(data = map_data$vic.df, fill = '#dadbdc') +
        geom_polygon() +
        geom_polygon(fill = 'transparent', colour = 'black',
                     size = 0.1, show.legend = FALSE) +
        geom_polygon(fill = 'transparent', colour = sa4_border,
                     size = 0.5, show.legend = FALSE,
                     data = map_data$sa4.df,
                     mapping = aes(long, lat, group = group)) +
        coord_equal(xlim = range(df$lon),
                    ylim = range(df$lat)) +
        scale_x_continuous(expand = c(0.025, 0)) +
        scale_y_continuous(expand = c(0.025, 0)) +
        theme(panel.background = element_rect(fill = 'white'),
              panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              plot.title = element_text(size = 16, hjust = 0.5))
    if (annotate) {
        ## Weight the contribution of each point by the distance between it
        ## and the next point.
        geom <- map_data$mel.df
        cycle_fn <- function(x) c(x[-1], x[1])
        geom$cyc_long <- ave(geom$long, geom$id, FUN = cycle_fn)
        geom$cyc_lat <- ave(geom$lat, geom$id, FUN = cycle_fn)
        geom$dist <- sqrt((geom$long - geom$cyc_long)^2 +
                          (geom$lat - geom$cyc_lat)^2)
        mean_longs <- sapply(split(geom, geom$id),
                             function(grp) weighted.mean(grp$long, grp$dist))
        mean_lats <- sapply(split(geom, geom$id),
                            function(grp) weighted.mean(grp$lat, grp$dist))
        geom$lat <- mean_lats[as.character(geom$id)]
        geom$long <- mean_longs[as.character(geom$id)]
        geom <- unique(geom[, c('SA3_CODE', 'long', 'lat')])
        geom$label <- geom$SA3_CODE %% 100
        geom$group <- 1
        geom$scale <- df$scale[1]
        p <- p + geom_text(data = geom,
                           mapping = aes(long, lat, label = label),
                           size = 3, fontface = 'bold')
    } else {
        ## Mark the CBD with a cross.
        cbd <- df[df$SA3_CODE == 20604, c('long', 'lat')]
        p <- p + annotate(geom = 'point',
                          x = mean(cbd$long), y = mean(cbd$lat),
                          shape = 4)
        ## When plotting the minimum model error, annotate those SA3s where
        ## the model predictions are never a good match.
        if ('err' %in% names(df) && all(df$scale == df$err)) {
            sa3_codes <- c(20607, 20703, 20801, 20802, 20901,
                           20903, 21103, 21201, 21202, 21301)
            geom <- df[df$SA3_CODE %in% sa3_codes, ]
            cycle_fn <- function(x) c(x[-1], x[1])
            geom$cyc_long <- ave(geom$long, geom$id, FUN = cycle_fn)
            geom$cyc_lat <- ave(geom$lat, geom$id, FUN = cycle_fn)
            geom$dist <- sqrt((geom$long - geom$cyc_long)^2 +
                              (geom$lat - geom$cyc_lat)^2)
            mean_longs <- sapply(split(geom, geom$id),
                                 function(grp) weighted.mean(grp$long, grp$dist))
            mean_lats <- sapply(split(geom, geom$id),
                                function(grp) weighted.mean(grp$lat, grp$dist))
            geom$lat <- mean_lats[as.character(geom$id)]
            geom$long <- mean_longs[as.character(geom$id)]
            geom <- unique(geom[, c('SA3_CODE', 'long', 'lat')])
            geom$label <- geom$SA3_CODE %% 100
            geom$group <- 1
            geom$scale <- df$scale[1]
            ## Manually tweak the position of this label, so that it is not
            ## obscured by the SA4 boundary.
            geom$long[geom$SA3_CODE == 20607] <- geom$long[geom$SA3_CODE == 20607] - 1250
            p <- p + geom_text(data = geom,
                               mapping = aes(long, lat, label = label),
                               size = 3, fontface = 'bold')
        }
    }
    p
}


map.default_opts <- function() {
    list(
        height = 4.5,
        width = 6.0,
        plot_type = 'continuous',
        diverge_at = 0,
        scale_min = NULL,
        scale_max = NULL,
        scale_log = FALSE,
        pad_cbar = 2,
        col_name = 'scale',
        facet_by <- c(),
        percentage = FALSE,
        annotate = FALSE,
        title = NULL)
}


map.parse_args <- function(args, reqd_args = c(1,2)) {
    opts <- map.default_opts()

    usage <- function(error = TRUE) {
           cat('\nUsage:', SCRIPT, '[options] input.ssv output.pdf\n\n')
           cat('Options:\n')
           cat('    -h, --help           Show this information\n')
           cat('\n')
           cat('    -c, --continuous     Data are continuous\n')
           cat('    -d, --discrete       Data are discrete\n')
           cat('    -v, --diverging      Data are diverging from zero\n')
           cat('\n')
           cat('    -n, --name COLUMN    The data column to plot\n')
           cat('    --diverge-at VAL     Data are diverging from VAL\n')
           cat('    --percentage         The data are percentages\n')
           cat('    --min MIN_VAL        The minimum value for the scale\n')
           cat('    --max MAX_VAL        The maximum value for the scale\n')
           cat('    --log                Use a logarithmic scale\n')
           cat('    -f, --facet COLUMN   Plot multiple facets\n')
           cat('\n')
           cat('    --annotate           Plot SA3 identifiers\n')
           cat('    -t, --title TITLE    Plot title\n')
           cat('\n')
           cat('    -w, --width LENGTH   Plot width (inches)\n')
           cat('    -g, --height LENGTH  Plot height (inches)\n')
           cat('    -p, --pad LENGTH     Colour bar padding (lines)\n')
           cat('                         Delete "map_shape_???.ssv" first\n')
           cat('\n')
           cat('Default values:\n')
           cat('    data:     ', opts$plot_type, '\n', sep='')
           cat('    name:     ', opts$col_name, '\n', sep='')
           cat('    title:    ', opts$title, '\n', sep='')
           cat('    width:    ', opts$width, '\n', sep='')
           cat('    height:   ', opts$height, '\n', sep='')
           cat('    pad:      ', opts$pad_cbar, '\n', sep='')
           cat('\n')

           if (error) {
               return(invisible(2))
           } else {
               return(invisible(0))
           }

    }

    ## Process leading command-line arguments beginning with '-'.
    while (length(args) >= 1 && substr(args[1], 1, 1) == '-') {
        if (args[1] == '--') {
            args <- args[-1]
            break
        } else if (args[1] %in% c('-h', '--help')) {
            return(usage(error = FALSE))
        } else if (args[1] %in% c('-c', '--continuous')) {
            opts$plot_type <- 'continuous'
            args <- args[-1]
        } else if (args[1] %in% c('-d', '--discrete')) {
            opts$plot_type <- 'discrete'
            args <- args[-1]
        } else if (args[1] %in% c('-v', '--diverging')) {
            opts$plot_type <- 'diverging'
            args <- args[-1]
        } else if (args[1] %in% c('-n', '--name')) {
            if (length(args) < 2) {
                return(usage())
            }
            opts$col_name <- args[2]
            args <- args[c(-1, -2)]
        } else if (args[1] %in% c('--diverge-at')) {
            if (length(args) < 2) {
                return(usage())
            }
            opts$diverge_at <- as.numeric(args[2])
            args <- args[c(-1, -2)]
        } else if (args[1] %in% c('--percentage')) {
            opts$percentage <- TRUE
            args <- args[-1]
        } else if (args[1] %in% c('--min')) {
            if (length(args) < 2) {
                return(usage())
            }
            opts$scale_min <- as.numeric(args[2])
            args <- args[c(-1, -2)]
        } else if (args[1] %in% c('--max')) {
            if (length(args) < 2) {
                return(usage())
            }
            opts$scale_max <- as.numeric(args[2])
            args <- args[c(-1, -2)]
        } else if (args[1] %in% c('--log')) {
            opts$scale_log <- TRUE
            args <- args[-1]
        } else if (args[1] %in% c('-f', '--facet')) {
            if (length(args) < 2) {
                return(usage())
            }
            opts$facet_by <- c(opts$facet_by, args[2])
            args <- args[c(-1, -2)]
            if (length(opts$facet_by) > 2) {
                cat('Error: too many facet columns\n')
                return(usage())
            }
        } else if (args[1] %in% c('--annotate')) {
            opts$annotate <- TRUE
            args <- args[-1]
        } else if (args[1] %in% c('-t', '--title')) {
            if (length(args) < 2) {
                return(usage())
            }
            opts$title <- args[2]
            args <- args[c(-1, -2)]
        } else if (args[1] %in% c('-g', '--height')) {
            if (length(args) < 2) {
                return(usage())
            }
            opts$height <- as.numeric(args[2])
            args <- args[c(-1, -2)]
        } else if (args[1] %in% c('-w', '--width')) {
            if (length(args) < 2) {
                return(usage())
            }
            opts$width <- as.numeric(args[2])
            args <- args[c(-1, -2)]
        } else if (args[1] %in% c('-p', '--pad')) {
            if (length(args) < 2) {
                return(usage())
            }
            opts$pad_cbar <- as.numeric(args[2])
            args <- args[c(-1, -2)]
        } else {
            cat(paste0('Error: unknown argument "', args[1], '"\n'))
            return(usage())
        }
    }

    ## Ensure any additional command-line arguments are expected.
    if (! length(args) %in% reqd_args) {
        return(usage())
    }

    opts$args <- args
    return(opts)
}


map.main <- function(args = NULL) {
    opts <- map.parse_args(args)
    if (! is.list(opts)) {
        return(invisible(opts))
    }

    in_file <- opts$args[1]
    if (length(opts$args) > 1) {
        plot_file <- opts$args[2]
    } else {
        re_ext <- '\\.[^.]+$'
        if (grepl(re_ext, in_file)) {
            ## Change the extension of the input file.
            plot_file <- sub(re_ext, '.pdf', in_file)
        } else {
            ## Append the PDF extension to the input file name.
            plot_file <- paste0(in_file, '.pdf')
        }
    }

    map.load_libs()
    map_data <- map.load_shapes()
    df <- read.table(in_file, header = TRUE, as.is = TRUE)
    if (! opts$col_name %in% names(df)) {
        cat(paste0('Error: missing data column "', opts$col_name, '"\n'))
        return(invisible(2))
    }
    if (opts$col_name == 'SA4') {
        SA4_names <- c('206' = 'Inner', '207' = 'Inner East',
                       '208' = 'Inner South', '209' = 'North East',
                       '210' = 'North West', '211' = 'Outer East',
                       '212' = 'South East', '213' = 'West')
        df[[opts$col_name]] <- paste0(
            SA4_names[as.character(df[[opts$col_name]])],
            ' (', df[[opts$col_name]], ')')
    }
    df$scale <- df[[opts$col_name]]

    plot_facet <- NULL
    if (length(opts$facet_by) == 1) {
        if (opts$facet_by[1] == 'label') {
            ## Use more descriptive names for each data source.
            new_labels <- c(
                'ABS (all modes)' = 'Journeys to work (all)',
                'ABS (private)' = 'Journeys to work (private)',
                'ABS (public)' = 'Journeys to work (public)',
                'Sygic' = 'GPS data'
            )
            df$label <- new_labels[df$label]
        }
        df$facet <- df[[opts$facet_by[1]]]
        plot_facet <- facet_wrap(~ facet)
    } else if (length(opts$facet_by) == 2) {
        df$facet_1 <- df[[opts$facet_by[1]]]
        df$facet_2 <- df[[opts$facet_by[2]]]
        plot_facet <- facet_wrap(facet_1 ~ facet_2)
    }

    if (opts$percentage) {
        scale_labels <- scales::unit_format(unit = '%', sep = '')
    } else {
        scale_labels <- waiver()
    }

    if (is.null(opts$scale_min) && is.null(opts$scale_max)) {
        scale_limits <- NULL
    } else if (is.null(opts$scale_min) || is.null(opts$scale_max)) {
        cat('Error: must specify minimum and maximum limits\n')
        return(invisible(2))
    } else {
        scale_limits <- c(opts$scale_min, opts$scale_max)
    }

    if (opts$plot_type == 'discrete') {
        df$scale <- factor(df$scale)
        p <- map.plot(map_data, df, annotate = opts$annotate) +
            scale_fill_brewer(opts$title, palette = 'Paired')
    } else if (opts$plot_type == 'diverging') {
        clrs <- brewer.pal(3, "BrBG")
        p <- map.plot(map_data, df, annotate = opts$annotate) +
            scale_fill_gradient2(opts$title,
                                 limits = scale_limits,
                                 labels = scale_labels,
                                 midpoint = opts$diverge_at,
                                 low = clrs[1],
                                 mid = clrs[2],
                                 high = clrs[3])
    } else if (opts$plot_type == 'continuous') {
        if (opts$scale_log) {
            scale_limits <- signif(range(df$scale), 2)
            p <- map.plot(map_data, df, annotate = opts$annotate) +
                scale_fill_viridis(opts$title,
                                   trans="log",
                                   breaks = scale_limits,
                                   limits = scale_limits,
                                   labels = scale_limits)
        } else {
            p <- map.plot(map_data, df, annotate = opts$annotate) +
                scale_fill_viridis(opts$title,
                                   limits = scale_limits,
                                   labels = scale_labels)
        }
    } else {
        cat('ERROR: unknown plot type "', opts$plot_type, '"\n', sep='')
        return(2)
    }

    ## Scale the colour bar to the height of the viewport ('npc'), and then
    ## subtract just more than 1 line of text ('lines') if there is a legend
    ## title.
    bar_height <- unit(1, 'npc') - unit(opts$pad_cbar, 'lines')
    if (! is.null(opts$title)) {
        bar_height <- bar_height - unit(1.10, 'lines')
    }

    if (opts$plot_type != 'discrete') {
        p <- p +
            guides(fill = guide_colorbar(barheight = bar_height))
    }

    if (! is.null(plot_facet)) {
        p <- p +
            plot_facet +
            theme(strip.background = element_rect(size = 0,
                                                  fill = 'transparent'),
                  strip.text = element_text(colour = 'black',
                                            size = rel(1.5)),
                  strip.text.x = element_text(margin = margin(b = 8, t = 4)),
                  panel.spacing.x = unit(2, 'lines'),
                  panel.spacing.y = unit(1, 'lines'))
    }

    p <- p +
        theme(plot.margin = margin(0, 0, 0, 0),
              legend.box.margin = margin(0, 0, 0, 0),
              panel.background = element_rect(size = 0),
              plot.background = element_rect(size = 0),
              legend.position = 'right',
              legend.justification = c(1, 0.5)
              )

    cat('Writing plot to', plot_file, '...\n')
    CairoPDF(plot_file, width = opts$width, height = opts$height,
             bg = 'transparent')
    print(p)
    dev.off()

    return(0)
}


##
## Only call main() if this script is being run from the command-line.
##
if (! interactive()) {
    file.arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
    file.path <- substring(file.arg, 8)
    file.regex <- paste0('^(|.*/)', SCRIPT)
    if (grepl(file.regex, file.path)) {
        quit(status = map.main(args = commandArgs(TRUE)))
    }
}
