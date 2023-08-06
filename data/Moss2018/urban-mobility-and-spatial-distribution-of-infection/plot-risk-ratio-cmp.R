#!/usr/bin/Rscript

library(ggplot2)
library(Cairo)

file_sim <- 'data/model-analysis/model-errors-sa3.ssv'
file_flu <- 'data/influenza/influenza-notif-ratio.ssv'
plot_file <- 'risk-ratio-cmp.pdf'

df_sim <- read.table(file_sim, header = TRUE, as.is = TRUE)

net_cases <- aggregate(cases ~ file + frac_self + frac_cbd, data = df_sim,
                       FUN = sum)
names(net_cases)[names(net_cases) == 'cases'] <- 'net_cases'
df_sim <- merge(df_sim, net_cases)
net_popn <- sum(unique(df_sim[, c('SA3', 'popn')])$popn)

df_sim$net_risk <- df_sim$net_cases / net_popn
df_sim$sa3_risk <- df_sim$cases / df_sim$popn
df_sim$rel_risk <- df_sim$sa3_risk / df_sim$net_risk

keep_cols <- c('file', 'frac_self', 'frac_cbd', 'SA3', 'popn',
               'net_risk', 'sa3_risk', 'rel_risk')
df_sim <- df_sim[, keep_cols]

df_flu <- read.table(file_flu, header = TRUE, as.is = TRUE)
df_sim <- merge(df_sim, unique(df_flu[, c('SA3', 'median_ratio')]))

## Calculate the interquartile range for influenza relative risks in each SA3.
df_lwr <- setNames(
    aggregate(ratio ~ SA3, data = df_flu,
              FUN = function(x) quantile(x, probs = 0.25, type = 1)),
    c('SA3', 'lwr_ratio'))
df_upr <- setNames(
    aggregate(ratio ~ SA3, data = df_flu,
              FUN = function(x) quantile(x, probs = 0.75, type = 1)),
    c('SA3', 'upr_ratio'))
df_flu <- merge(merge(df_flu, df_lwr), df_upr)

SA4_names <- c('206' = 'Inner', '207' = 'Inner East',
               '208' = 'Inner South', '209' = 'North East',
               '210' = 'North West', '211' = 'Outer East',
               '212' = 'South East', '213' = 'West')

df_flu$SA4 <- substring(df_flu$SA3, 1, 3)
df_flu$SA4_name <- SA4_names[df_flu$SA4]
df_flu$SA3_num <- as.integer(substring(df_flu$SA3, 4, 5))
df_flu$SA3_label <- paste0('"', df_flu$SA4_name, ' "*', df_flu$SA3_num)
df_flu$SA3_label[df_flu$SA3 == 20604] <- 'bold(CBD)'
df_flu$SA3_label <- factor(df_flu$SA3_label, ordered = TRUE,
                           levels = unique(df_flu$SA3_label))

df_flu_med <- unique(df_flu[, c('SA3_label', 'median_ratio',
                                'lwr_ratio', 'upr_ratio')])

labels <- c('abs_all' = 'Journeys to work (all)',
            'abs_private' = 'Journeys to work (private)',
            'abs_public' = 'Journeys to work (public)',
            'sygic' = 'GPS data')
df_sim$label <- labels[df_sim$file]

df_sim$SA4 <- substring(df_sim$SA3, 1, 3)
df_sim$SA4_name <- SA4_names[df_sim$SA4]
df_sim$SA3_num <- as.integer(substring(df_sim$SA3, 4, 5))
df_sim$SA3_label <- paste0('"', df_sim$SA4_name, ' "*', df_sim$SA3_num)
df_sim$SA3_label[df_sim$SA3 == 20604] <- 'bold(CBD)'
df_sim$SA3_label <- factor(df_sim$SA3_label, ordered = TRUE,
                           levels = unique(df_sim$SA3_label))

## Add points for the non-spatial model.
df_1 <- df_sim[df_sim$frac_self == max(df_sim$frac_self), ]
df_1$frac_self <- 1
df_1$rel_risk <- 1
df_sim <- rbind(df_sim, df_1)

plot_files <- c('abs_all', 'abs_private', 'abs_public', 'sygic')
p1 <-  ggplot(data = df_sim[df_sim$file %in% plot_files, ],
             aes(x = SA3_label, y = rel_risk, colour = frac_self)) +
    geom_hline(yintercept = 1) +
    geom_point() +
    geom_point(mapping = aes(x = factor(SA3_label), y = median_ratio),
               colour = 'black', shape = 4) +
    scale_colour_distiller(palette = 'Spectral') +
    facet_grid(frac_cbd ~ label)

## Measure the difference between the model risk ratios and the observed
## influenza case risk ratios, taking mean (mae) and population-weighted (wae)
## sums.
##
## NOTE: we ignore SA3 21002, it has an influenza case risk ratio of 0.00074.
df_sim$rr_diff <- abs(df_sim$rel_risk - df_sim$median_ratio)
df_sim$rr_wdiff <- df_sim$rr_diff * 37 * df_sim$popn / net_popn

df_mae <- aggregate(rr_diff ~ file + frac_self + frac_cbd,
                    data = df_sim[df_sim$SA3 != '21002', ], FUN = sum)
df_wae <- aggregate(rr_wdiff ~ file + frac_self + frac_cbd,
                    data = df_sim[df_sim$SA3 != '21002', ], FUN = sum)

df_mae <- df_mae[order(df_mae$rr_diff), ]
df_wae <- df_wae[order(df_wae$rr_wdiff), ]

## Values of frac_cbd: 0.20 0.25  0.33 0.50 0.67 0.75 0.80
df_plot <- df_sim[df_sim$file %in% plot_files
                  ## & df_sim$frac_self >= 0.75
                  & df_sim$frac_cbd == 0.20, ]

p0 <- ggplot(data = df_plot,
             aes(x = SA3_label, y = rel_risk, colour = frac_self,
                 group = interaction(frac_cbd, frac_self))) +
    geom_hline(yintercept = 1) +
    geom_point() +
    # geom_line() +
    xlab(NULL) +
    ylab("Relative Risk") +
    scale_colour_distiller(expression("\u03b4"[i]^H),
                           palette = 'Spectral',
                           breaks = range(df_plot$frac_self),
                           labels = range(df_plot$frac_self)) +
    scale_x_discrete(labels = parse(text = levels(df_sim$SA3_label))) +
    coord_cartesian(ylim = c(0.4, 2.1)) +
    facet_wrap(~ label) +
    theme_light() +
    guides(colour = guide_colourbar(direction = 'horizontal',
                                    title.vjust = 1)) +
    theme(strip.background = element_rect(size = 0,
                                          fill = 'transparent'),
          strip.text = element_text(colour = 'black',
                                    size = rel(1.25)),
          strip.text.x = element_text(margin = margin(b = 8, t = 4)),
          panel.spacing.x = unit(2, 'lines'),
          panel.spacing.y = unit(1, 'lines'),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = 'top',
          legend.title = element_text(size = rel(1.1)),
          legend.text = element_text(size = rel(1.1)))

p1 <- p0 +
    geom_point(mapping = aes(x = SA3_label, y = median_ratio),
               colour = 'black', shape = 4, size = 2)

p2 <- p0 +
    geom_pointrange(data = df_flu_med,
                    mapping = aes(x = SA3_label, y = median_ratio,
                                  ymin = lwr_ratio, ymax = upr_ratio,
                                  group = 1),
                    colour = 'black', shape = 4)

p3 <- ggplot(df_flu, aes(SA3_label, ratio)) +
    geom_hline(yintercept = 1) +
    geom_point(colour = '#7f7f7f', shape = 1) +
    geom_point(data = df_flu, mapping = aes(y = median_ratio),
               colour = 'black') +
    geom_line(data = df_flu, mapping = aes(y = median_ratio),
              colour = 'black',
              group = 1) +
    expand_limits(y = 0) +
    scale_x_discrete(labels = parse(text = levels(df_sim$SA3_label))) +
    xlab(NULL) +
    ylab('Relative Risk') +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


cat('Writing plot to', plot_file, '...\n')
CairoPDF(plot_file, width = 10, height = 6, bg = 'transparent')
print(p1)
print(p2)
print(p3)
invisible(dev.off())
