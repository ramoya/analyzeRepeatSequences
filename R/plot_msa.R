
#' Plot a multiple sequence alignment
#'
#' @description
#' Plots the alignment of multiple repetitive sequences.
#'
#'
#' @param unitorders_df         `data.frame` of a parsed multiple sequence
#'   alignment. Has columns `character`, `sample`, `count`, `seq`, `gp`. See
#'   dataset `unitorders_example` for format.
#' @param unit2color_df         `data.frame` of the conversion between repeat
#'   units and colors. Has columns `seq`, `color`. See dataset `unit2color` for
#'   format.
#' @param unitsToHighlight       Vector of repeat units to show in
#'   color. All other repeat units are shown in beige. Defaults to a zero length
#'   vector.
#' @param showLegend            `logical` indicating whether to plot the legend.
#'   Defaults to FALSE.
#' @param plotSubgroups         `logical` indicating whether `unitorders_df` has
#'   column `group` for clustering the sequences. Defaults to FALSE. Not
#'   compatible with neither plotVariabilityScore=T nor plotMinMax=T.
#' @param plotVariabilityScore  `logical` indicating whether to plot a score
#'   across the alignment representing repeat unit variability. See
#'   `calculate_variability_score` for details. Returns a list of two ggplot objects.
#' @param plotMinMax            `logical` indicating whether to show x-axis
#'   ticks at the first and last repeat units.
#'
#' @returns ggplot object showing the alignment of multiple repetitive sequences.
#' @examples
#' plot_msa(unitorders_example, unit2color)
#' plot_msa(unitorders_example, unit2color, unitsToHighlight = c("GATCCTGACCTTACTAGTTTACAATCACAG", "GACCCTGACCTGACTAGTTTACAATCACAT"))
#' plot_msa(unitorders_example, unit2color, showLegend = T)
#' plot_msa(unitorders_example, unit2color, showLegend = T) + theme(axis.text.y = element_text())
#' plot_msa(unitorders_withSubgroups_example, unit2color, plotSubgroups = T)
#' plot_msa(unitorders_example, unit2color, plotVariabilityScore = T)
#' plot_msa(unitorders_example, unit2color, plotMinMax = T)

#' @export
#' @import ggplot2
plot_msa <- function(unitorders_df,
                     unit2color_df,
                     unitsToHighlight = c(),
                     showLegend = F,
                     plotSubgroups = F,
                     plotVariabilityScore = F,
                     plotVariableRegion = F,
                     plotMinMax = F) {

  if(plotSubgroups == T & (plotVariabilityScore == T | plotMinMax == T)){
    stop("plotSubgroups=T not compatible with neither plotVariabilityScore=T nor plotMinMax=T")
  }

  l <- format_unitorders_to_plot(unitorders_df, unit2color_df,
                                 unitsToHighlight)
  unitorders_df <- l[[1]]
  ranked_units_in_plot <- l[[2]]
  ranked_colors_in_plot <- l[[3]]
  unit_counts_calculated_plus_colors_reordered <- l[[4]]

 #if alignment has gaps, don't show in legend
 if("-" %in% unit_counts_calculated_plus_colors_reordered$seq){
   leng = length(ranked_units_in_plot) - 1
 } else {
   leng = length(ranked_units_in_plot)
 }

  p <- ggplot(unitorders_df) +
    geom_bar(
      aes(
        x = sample,
        y = count,
        fill = seq,
        group = factor(gp)
      ),
      stat = "identity",
      width = 1,
      color = NA
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(
        family = "Courier",
        face = "bold",
        size = 8
      ),
      legend.text = element_text(family = "Courier", size = 8),
      panel.background = element_rect(fill = "black"),
      panel.grid = element_blank(),
      plot.background = element_blank(),
      axis.line = element_blank()
    ) +
    labs(fill = "Repeat units (frequency)") +
    guides(
      color = guide_legend(override.aes = list(size = 2)),
      fill = guide_legend(title.position = "top")
    ) +
    ##LEFT OFF HERE DO I NEED TO BE DROPPING THE LAST ROW??
    scale_fill_manual(
      breaks = ranked_units_in_plot[1:leng],
      values = ranked_colors_in_plot,
      labels = paste0(
        ranked_units_in_plot[1:leng],
        " (",
        round(
          unit_counts_calculated_plus_colors_reordered$count[1:leng] / sum(unit_counts_calculated_plus_colors_reordered$count[1:leng]),
          digits = 4
        ),
        ")"
      ),
      na.value = "black"
    ) +
    coord_flip()

  if (showLegend) {
    p <- p + theme(legend.position = "bottom")
  }

  if (plotSubgroups) {
    p <- p + facet_grid(
      rows = vars(group),
      scales = "free",
      space = "free"
    ) +
      theme(
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
      )
  }

  if (plotVariabilityScore) {
    # Define variability score
    variability_score_smoothed_df <- calculate_variability_score(unitorders_df)

    max_x_val <- max(variability_score_smoothed_df$x)

    # Line plot of variability score
    l <- ggplot(data = variability_score_smoothed_df) +
      geom_line(aes(x = x, y = y), size = 0.75) +
      scale_x_continuous(
        expand = c(0, 0),
        breaks = seq(0, max(variability_score_smoothed_df$x), by = 100)
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.4),
        expand = c(0, 0)
      ) +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), "pt")
      )

    # MSA
    p <- p +
      theme(
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "in"),
        panel.border = element_blank()
      )

    p <- list(l, p)
  }

  if (plotMinMax) {
    p <- p +
      scale_y_continuous(breaks = c(
        0,
        (
          unitorders_df %>% dplyr::group_by(sample) %>% dplyr::summarize(max = sum(count)) %>%
            dplyr::pull(max)
        )[1]),
        expand = c(0, 0)
      )
  }

  return(p)
}


#' Format `unitorders_df`
#'
#' @description Formats `unitorders_df` for plotting.
#'
#'
#' @param unitorders_df         `data.frame` of a parsed multiple sequence
#'   alignment. Has columns `character`, `sample`, `count`, `seq`, `gp`. See
#'   dataset `unitorders_example` for format.
#' @param unit2color_df         `data.frame` of the conversion between repeat
#'   units and colors. Has columns `seq`, `color`. See dataset `unit2color` for
#'   format.
#' @param unitsToHighlight       Vector of repeat units to show in
#'   color. All other repeat units are shown in beige. Defaults to a zero length
#'   vector.
#'
#' @returns a list of length 4 with the objects used to plot.


#' @export
format_unitorders_to_plot <- function(unitorders_df,
                                      unit2color_df,
                                      unitsToHighlight = c()) {

  unitorders_df$seq[unitorders_df$seq == ""] <- "-"

  # Keep color of units to highlight. Set all other colors to beige.
  if (length(unitsToHighlight) != 0) {
    unitorders_df[!(unitorders_df$seq %in% unitsToHighlight) &
                    unitorders_df$seq != "-", "seq"] <- "Other"
    new_color_df <- data.frame("#E9DAC4", "Other")
    names(new_color_df) <- c("color", "seq")
    unit2color_df <- rbind(unit2color_df, new_color_df)
  }

  # Calculate number of each repeat unit
  unit_counts_calculated <- unitorders_df %>%
    dplyr::group_by(seq) %>%
    dplyr::summarize(count = dplyr::n())

  print(unit_counts_calculated, n=50)

  # Add column for color
  if (length(unitsToHighlight) != 0) {
    unit_counts_calculated_plus_colors <- plyr::join(unit2color_df, unit_counts_calculated)
    unit_counts_calculated_plus_colors[is.na(unit_counts_calculated_plus_colors$count), "count"] <- 0
  } else {
    unit_counts_calculated_plus_colors <- plyr::join(unit_counts_calculated, unit2color_df)
  }

  # Order these three rows lowest, but only if they're in the df
  tmp_df <- unit_counts_calculated_plus_colors[order(unit_counts_calculated_plus_colors[, "count"], decreasing = T), ]

  if("30mer with frequency <= 0.0005" %in% tmp_df$seq){
    tmp_df <- rbind(tmp_df[!tmp_df$seq == "30mer with frequency <= 0.0005", ], tmp_df[tmp_df$seq == "30mer with frequency <= 0.0005", ])
  }
  if("Non-30mer (collective frequency = 0.00041)" %in% tmp_df$seq){
    tmp_df <- rbind(tmp_df[!tmp_df$seq == "Non-30mer (collective frequency = 0.00041)", ], tmp_df[tmp_df$seq == "Non-30mer (collective frequency = 0.00041)", ])
  }
  if("-" %in% tmp_df$seq){
    tmp_df <- rbind(tmp_df[!tmp_df$seq == "-", ], tmp_df[tmp_df$seq == "-", ])
  }
  unit_counts_calculated_plus_colors_reordered <- tmp_df



  print(unit_counts_calculated_plus_colors_reordered)
  ranked_colors_in_plot <- unit_counts_calculated_plus_colors_reordered$color
  ranked_units_in_plot <- unit_counts_calculated_plus_colors_reordered$seq

  # Assign levels to get colors and sample order right
  unitorders_df$seq <- factor(unitorders_df$seq, levels = ranked_units_in_plot)

  return(list(unitorders_df, ranked_units_in_plot, ranked_colors_in_plot, unit_counts_calculated_plus_colors_reordered))
}



