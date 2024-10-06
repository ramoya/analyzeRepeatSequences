# dfs beginning with unitorders are in a format for making plots & they have infrequent typical and atypical length units masked

#' Plot a multiple sequence alignment
#'
#' @description
#' `plot_msa` returns a ggplot object showing the alignment of multiple repetitive sequences.
#'
#'
#' @param unitorders_df         `data.frame` of a parsed multiple sequence
#'   alignment. Has columns `character`, `sample`, `count`, `seq`, `gp`. See
#'   dataset `unitorders_example` for format.
#' @param unit2ascii_df         `data.frame` of the conversion between repeat
#'   units and single ASCII characters. Has columns `seq`, `character`. See
#'   dataset `unit2ascii` for format.
#' @param unit2color_df         `data.frame` of the conversion between repeat
#'   units and colors. Has columns `seq`, `color`. See dataset `unit2color` for
#'   format.
#' @param unitsToHighlight       Vector of selected repeat units to show in
#'   color. All other repeat units are shown in beige. Defaults to zero length
#'   vector.
#' @param showLegend            `logical` indicating whether to plot the legend. Defaults to FALSE.
#' @param plotSubgroups         `logical` indicating whether `unitorders_df` has
#'   column `group` for clustering the sequences. Defaults to FALSE.
#' @param plotVariableRegion    `logical` indicating if input is . Defaults to FALSE.
#'
#' @param plotMinMax            `logical` indicating

#' @export
#' @import ggplot2
plot_msa <- function(unitorders_df,
                     unit2ascii_df,
                     unit2color_df,
                     unitsToHighlight = c(),
                     showLegend = F,
                     plotSubgroups = F,
                     plotVariabilityScore = F,
                     plotVariableRegion = F,
                     plotMinMax = F) {
  unitorders_df$seq[unitorders_df$seq == ""] <- "-"

  if (length(unitsToHighlight) != 0) {
    unitorders_df[!(unitorders_df$seq %in% unitsToHighlight) &
      unitorders_df$seq != "-", "seq"] <- "Other"
    # Consider all other repeat units as "other"
    new_color_df <- data.frame("#E9DAC4", "Other")
    names(new_color_df) <- c("color", "seq")
    unit2color_df <- rbind(unit2color_df, new_color_df)
  }

  # Calculate number of each repeat unit
  unit_counts_calculated <- unitorders_df %>%
    dplyr::group_by(seq) %>%
    dplyr::summarize(count = dplyr::n())

  # Add column for color
  if (length(unitsToHighlight) != 0) {
    unit_counts_calculated_plus_colors <- plyr::join(unit2color_df, unit_counts_calculated)
    unit_counts_calculated_plus_colors[is.na(unit_counts_calculated_plus_colors$count), "count"] <- 0
  } else {
    unit_counts_calculated_plus_colors <- plyr::join(unit_counts_calculated, unit2color_df)
  }

  # Order these three rows lowest
  tmp_df <- unit_counts_calculated_plus_colors[order(unit_counts_calculated_plus_colors[, "count"], decreasing = T), ]
  tmp_df_2 <- rbind(tmp_df[!tmp_df$seq == "30mer with frequency <= 0.0005", ], tmp_df[tmp_df$seq == "30mer with frequency <= 0.0005", ])
  tmp_df_3 <- rbind(tmp_df_2[!tmp_df_2$seq == "Non-30mer (collective frequency = 0.00041)", ], tmp_df_2[tmp_df_2$seq == "Non-30mer (collective frequency = 0.00041)", ])
  unit_counts_calculated_plus_colors_reordered <- rbind(tmp_df_3[!tmp_df_3$seq == "-", ], tmp_df_3[tmp_df_3$seq == "-", ])


  ranked_colors_in_plot <- unit_counts_calculated_plus_colors_reordered$color
  ranked_units_in_plot <- unit_counts_calculated_plus_colors_reordered$seq

  # Assign levels to get colors and sample order right
  unitorders_df$seq <- factor(unitorders_df$seq, levels = ranked_units_in_plot)

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
    scale_fill_manual(
      breaks = ranked_units_in_plot[1:length(ranked_units_in_plot) - 1],
      values = ranked_colors_in_plot,
      labels = paste0(
        ranked_units_in_plot[1:length(ranked_units_in_plot) - 1],
        " (",
        round(
          unit_counts_calculated_plus_colors_reordered$count[1:length(ranked_units_in_plot) - 1] / sum(unit_counts_calculated_plus_colors_reordered$count[1:length(ranked_units_in_plot) - 1]),
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



  if (plotVariableRegion) {
    p <- ggplot(
      unitorders_df,
      aes(
        x = allele_num,
        y = width,
        fill = factor(seq),
        group = factor(gp)
      )
    ) +
      geom_bar(
        stat = "identity",
        width = 0.8,
        color = NA
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(
          family = "Courier",
          face = "bold",
          size = 9
        ),
        legend.text = element_text(family = "Courier", size = 7)
      ) +
      labs(fill = "Repeat units (total count)") +
      guides(
        color = guide_legend(override.aes = list(size = 2)),
        fill = guide_legend(title.position = "top")
      ) +
      scale_fill_manual(
        breaks = ranked_units_in_plot,
        values = ranked_colors_in_plot,
        labels = paste0(
          ranked_units_in_plot,
          " (",
          unit_counts_calculated_plus_colors_reordered$count,
          ")"
        )
      ) +
      coord_flip()
  }

  if (plotMinMax) {
    p <- p +
      scale_y_continuous(breaks = c(
        0,
        (
          unitorders_df %>% dplyr::group_by(sample) %>% dplyr::summarize(max = sum(count)) %>%
            dplyr::pull(max)
        )[1],
        expand = c(0, 0)
      ))
  }

  return(p)
}

#' @export
calculate_consensus <- function(x, unit_frequency_df) {
  ordered_unit_frequencies <- sort(prop.table(table(x)), decreasing = TRUE)

  most_common_unit <- names(ordered_unit_frequencies[1])

  most_common_unit_frequency <- ordered_unit_frequencies[1]

  # If most common unit is a gap and it doesn't reach 65% frequency, choose most common non-gap unit
  if (most_common_unit == "-" &
    most_common_unit_frequency <= 0.65) {
    unit <- names(ordered_unit_frequencies[2])
    frequency <- ordered_unit_frequencies[2]
  } else {
    unit <- most_common_unit
    frequency <- most_common_unit_frequency
  }

  ordered_unit_frequencies_nogap <- ordered_unit_frequencies[names(ordered_unit_frequencies) != "-"]

  # If unit has equal frequency with another unit, choose the more abundant unit
  if (sum(frequency == ordered_unit_frequencies_nogap, na.rm = T) > 1) {
    # Merge unit_frequency_df with ASCII to create cipher
    unit_frequency_df <- merge(
      unit_frequency_df,
      unit_to_ascii_noduplicates,
      by.x = "seq",
      by.y = "seq"
    )
    # Sort unit_frequency_df according to frequency, decreasing
    unit_frequency_df <- unit_frequency_df[order(-unit_frequency_df$Frequency), ]
    characters_with_tied_frequencies <- names(ordered_unit_frequencies_nogap[frequency == ordered_unit_frequencies_nogap])

    unit <- unit_frequency_df[unit_frequency_df$character %in% characters_with_tied_frequencies, "character"][1]
  }

  return(c(unit, frequency))
}

#' @export
call_consensus_sequence <- function(unitorders_df, unit_frequency_df) {
  # Give each position a number in each sequence
  unitorders_df <- as.data.frame(unitorders_df %>% dplyr::group_by(sample) %>% dplyr::mutate(position = rev(1:dplyr::n())))
  # Mutate the df to have each position as a column
  unitorders_wide_df <- as.data.frame(tidyr::pivot_wider(
    unitorders_df[, c("sample", "position", "character")],
    names_from = position,
    values_from = character
  ))
  # Set the row names as the sample column
  row.names(unitorders_wide_df) <- unitorders_wide_df[, 1]
  unitorders_wide_df[, 1] <- NULL

  # For every column, call consensus unit
  consensus_seq_as_list <- apply(
    unitorders_wide_df,
    2,
    FUN = function(y) {
      calculate_consensus(y, unit_frequency_df = unit_frequency_df)[1]
    }
  )

  # Reverse the order of the units and collapse this list into a string
  consensus_seq_as_string <- paste(rev(consensus_seq_as_list), collapse = "")

  return(list(consensus_seq_as_string, unitorders_wide_df))
}

#' @export
make_unitorders_consensus <- function(consensus_seq_as_string, unit2ascii_df) {
  # Expand string to list of characters
  consensus_seq_as_list <- rev(unlist(strsplit(consensus_seq_as_string, "")))

  unitorders_consensus_df <- data.frame(
    character = consensus_seq_as_list,
    sample = rep("Consensus sequence", length(consensus_seq_as_list)),
    count = rep(2, length(consensus_seq_as_list))
  )
  unitorders_consensus_df$gp <- seq(nrow(unitorders_consensus_df))
  # Join with ASCII conversion table
  unitorders_consensus_df <- plyr::join(unitorders_consensus_df, unit2ascii_df)

  return(unitorders_consensus_df)
}

#' @export
calculate_fractional_abundance <- function(unitorders_group) {
  # At each position count fractional abundance of each unit
  # Give each position a number in each sequence
  unitorders_group_numbered <- as.data.frame(unitorders_group %>% dplyr::group_by(sample) %>% dplyr::mutate(position = rev(1:dplyr::n())))
  # Add dash to seq column for gaps
  unitorders_group_numbered$seq[unitorders_group_numbered$seq == ""] <- "-"
  # Calculate fractional abundance
  unit_occurrence_by_position <- unitorders_group_numbered %>%
    dplyr::group_by(position, seq) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::mutate(proportion = (n / sum(n)) * 100) %>%
    dplyr::arrange(seq, proportion, .by_group = T)

  return(unit_occurrence_by_position)
}

plot_fractional_abundance <- function(unitorders_group, unit2color_df) {
  # This function plots one group's fractional unit abundance
  unit_occurrence_by_position <- calculate_fractional_abundance(unitorders_group)
  # For plotting remove the gaps and make background black, fractional abundance not represented by
  unit_occurrence_by_position <- subset(unit_occurrence_by_position, seq != "-")

  # Get color
  unit_counts_calculated <- unitorders_group %>%
    dplyr::group_by(seq) %>%
    dplyr::summarize(count = dplyr::n())
  unit_counts_calculated_plus_colors <- plyr::join(unit2color_df, unit_counts_calculated)
  unit_counts_calculated_plus_colors[is.na(unit_counts_calculated_plus_colors$count), "count"] <- 0
  tmp_df <- unit_counts_calculated_plus_colors[order(unit_counts_calculated_plus_colors[, "count"], decreasing = T), ]
  tmp_df_2 <- rbind(tmp_df[!tmp_df$seq == "30mer with frequency <= 0.0005", ], tmp_df[tmp_df$seq == "30mer with frequency <= 0.0005", ])
  tmp_df_3 <- rbind(tmp_df_2[!tmp_df_2$seq == "Non-30mer (collective frequency = 0.00041)", ], tmp_df_2[tmp_df_2$seq == "Non-30mer (collective frequency = 0.00041)", ])
  unit_counts_calculated_plus_colors_reordered <- rbind(tmp_df_3[!tmp_df_3$seq == "-", ], tmp_df_3[tmp_df_3$seq == "-", ])
  ranked_colors_in_plot <- unit_counts_calculated_plus_colors_reordered$color
  ranked_units_in_plot <- unit_counts_calculated_plus_colors_reordered$seq

  # Plot as a barplot
  fractional_abundance_plot <- ggplot(unit_occurrence_by_position) +
    geom_bar(
      aes(
        x = position,
        y = proportion,
        fill = factor(seq),
        group = proportion
      ),
      stat = "identity",
      width = 1,
      position = "stack"
    ) +
    scale_fill_manual(
      breaks = ranked_units_in_plot,
      values = ranked_colors_in_plot,
      labels = paste0(
        ranked_units_in_plot,
        " (",
        unit_counts_calculated_plus_colors_reordered$count,
        ")"
      )
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = c(0, 0),
      name = "Percentage of repeat units"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      # axis.text.y = element_text(family = 'Courier', face = 'bold', size = 6),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "black")
    ) +
    ylab("Percentage of repeat units")

  return(fractional_abundance_plot)
} # useful for plotting fractional abundance of one consensus sequence

#' @export
make_repeat_summary_figure <- function(unitorders_group,
                                       ...,
                                       unit2color_df,
                                       unit2ascii_df,
                                       unit_frequency_df) {
  unitorders_dfs <- list(unitorders_group, ...)
  # Call consensus
  # Calculate fractional abundance and put in df unit_occurrences_by_position_dfs
  # unit_occurrences_by_position_dfs is a list of dataframes, each with the columns position, seq, n, proportion
  consensuses <- lapply(unitorders_dfs, function(x) {
    call_consensus_sequence(x, unit_frequency_df)
  })
  consensus_dfs <- lapply(
    seq(1:length(consensuses)),
    FUN = function(x) {
      data.frame(
        character = unlist(stringr::str_split(consensuses[[x]][[1]], pattern = "")),
        position = seq(1, length(unlist(
          stringr::str_split(consensuses[[x]][[1]], pattern = "")
        ))),
        source = x
      ) %>%
        magrittr::set_colnames(c("character", "position", "source"))
    }
  ) # use string of consensus sequence to get source (Type), position, and character

  # get fractional abundance at each position
  # VERSION WITH PROP.TABLE
  make_fractional_abundance_df <- function(x) {
    tmp <- do.call("rbind", as.list(unlist(
      apply(
        call_consensus_sequence(x, unit_frequency_df)[[2]],
        2,
        FUN = function(y) {
          sort(prop.table(table(y)), decreasing = TRUE)
        }
      ),
      recursive = F
    ))) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("to_split") %>%
      tidyr::separate_wider_delim(
        cols = "to_split",
        delim = ".",
        names = c("position", "character")
      ) %>%
      magrittr::set_colnames(c("position", "character", "proportion"))

    # convert character to seq
    tmp <- merge(tmp, unit2ascii_df, by = "character", all.x = T)
    tmp$position <- as.numeric(tmp$position)
    tmp$n <- tmp$proportion * (rep(dim(
      call_consensus_sequence(x, unit_frequency_df)[[2]]
    )[1], length(tmp$proportion)))
    tmp <- tmp %>%
      dplyr::ungroup() %>%
      dplyr::arrange(position, -proportion)
    return(tmp)
  }

  fractional_abundance_dfs <- lapply(unitorders_dfs, make_fractional_abundance_df)

  # combine consensus dfs
  consensus <- dplyr::bind_rows(consensus_dfs, .id = "source")
  consensus <- merge(consensus, unit2ascii_df, by = "character", all.x = T)
  consensus[is.na(consensus$seq), "seq"] <- "-"

  # combine fractional_abundance dfs
  unit_occurrence_by_position <- dplyr::bind_rows(fractional_abundance_dfs, .id = "source")
  unit_occurrence_by_position$source <- as.numeric(unit_occurrence_by_position$source)
  unit_occurrence_by_position[is.na(unit_occurrence_by_position$seq), "seq"] <- "-"

  # #previous
  # unit_occurrences_by_position_dfs <- lapply(unitorders_dfs, function(x) calculate_fractional_abundance(x))
  # #VERSION WITH n()
  # lapply(unitorders_dfs, function(x)
  #   y <- call_consensus_sequence(x, unit_counts_table_s5)[[2]] %>%
  #   tibble::rownames_to_column("sample") %>%
  #   tidyr::pivot_longer(cols = -c("sample"), names_to = "position", values_to = "character") %>%
  #   dplyr::group_by(position, character) %>% dplyr::summarize(n = dplyr::n()) %>%
  #   dplyr::mutate(proportion = (n / sum(n))*100) %>% dplyr::arrange(character, proportion, .by_group = T) %>%
  #   dplyr::select(position, character, proportion)
  #   # y$position <- as.numeric(y$position)
  #   # y <- y %>% dplyr::arrange(position, proportion)
  #   )
  # unit_occurrence_by_position <- dplyr::bind_rows(unit_occurrences_by_position_dfs, .id = 'source')
  # unit_occurrence_by_position$source <- as.numeric(unit_occurrence_by_position$source)
  # unit_occurrence_by_position_gapsdropped <- subset(unit_occurrence_by_position, !(seq == '-' & proportion <= 65))
  # consensus <- unit_occurrence_by_position_gapsdropped %>% dplyr::group_by(source, position) %>% dplyr::slice_max(proportion) %>% dplyr::ungroup()

  # Calculate new unit counts
  unit_counts_calculated <- consensus %>%
    dplyr::group_by(seq) %>%
    dplyr::summarize(count = dplyr::n())

  # Add colors to unit counts
  # The order here is important, default is left join
  # Need all units in legend for fractional abundance plot
  # However here I'm counting units in consensus
  # There's one unit in unit color df that isn't in any consensus
  # Still want this unit in the legend, give it a count of 0
  unit_counts_calculated_plus_colors <- plyr::join(unit2color_df, unit_counts_calculated)
  unit_counts_calculated_plus_colors[is.na(unit_counts_calculated_plus_colors$count), "count"] <- 0
  # Reorder samples
  # TODO Ideally these string values wouldn't be hard coded
  tmp_df <- unit_counts_calculated_plus_colors[order(unit_counts_calculated_plus_colors[, "count"], decreasing = T), ]
  tmp_df_2 <- rbind(tmp_df[!tmp_df$seq == "30mer with frequency <= 0.0005", ], tmp_df[tmp_df$seq == "30mer with frequency <= 0.0005", ])
  tmp_df_3 <- rbind(tmp_df_2[!tmp_df_2$seq == "Non-30mer (collective frequency = 0.00041)", ], tmp_df_2[tmp_df_2$seq == "Non-30mer (collective frequency = 0.00041)", ])
  unit_counts_calculated_plus_colors_reordered <- rbind(tmp_df_3[!tmp_df_3$seq == "-", ], tmp_df_3[tmp_df_3$seq == "-", ])

  ranked_colors_in_plot <- unit_counts_calculated_plus_colors_reordered$color
  ranked_units_in_plot <- unit_counts_calculated_plus_colors_reordered$seq

  # Right before plotting, remove gaps
  unit_occurrence_by_position <- subset(unit_occurrence_by_position, seq != "-")

  # Assign levels of unit sequences to get colors and sample order right
  unit_occurrence_by_position$seq <- factor(unit_occurrence_by_position$seq, levels = ranked_units_in_plot)

  rect_boundaries <- unit_occurrence_by_position %>%
    dplyr::group_by(source) %>%
    dplyr::summarize(
      xmax = max(position),
      xmin = 1,
      ymin = 0,
      ymax = Inf
    )

  # Use facetscales to customize the xaxis on every plot
  # I need a list of scales
  #
  # `4` = scale_y_continuous(limits = c(5, 25), breaks = seq(5, 25, 5)),
  # `f` = scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)),
  # `r` = scale_y_continuous(limits = c(10, 20), breaks = seq(10, 20, 2))
  #
  # # How to decide the interval for the breaks?
  # # Any given plot should only have 4 axis ticks
  #
  #
  #
  # scales_x <- unit_occurrence_by_position %>% group_by(source) %>% summarize(x_min=min(position),
  #                                                                x_max=floor(max(position)/50)*50,
  #                                                                breaks = floor((x_max/4)/50)*50) %>%
  #   stringr::str_glue_data(
  #     "`{source}` = scale_y_continuous(limits = c({x_min}, {x_max}), ",
  #     "breaks = seq({x_min}, {x_max}, {breaks}))") %>%
  #   stringr::str_flatten(", ") %>%
  #   stringr::str_c("list(", ., ")") %>%
  #   parse(text = .) %>%
  #   eval()


  seq_groups <- unique(unit_occurrence_by_position$source)
  max_x_all_groups <- max(unit_occurrence_by_position$position) + 1
  # Calculate relative widths of the groups
  relative_widths <- unlist(
    unit_occurrence_by_position %>% dplyr::group_by(source) %>% dplyr::summarize(x_max = max(position)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        max_x_max = max(x_max) * 1.2,
        rel_width = x_max / max_x_max
      ) %>% dplyr::select(rel_width)
  )

  plotOneGroup <- function(grp) {
    rng <- range(unit_occurrence_by_position[unit_occurrence_by_position$source == grp, "position"])
    o <- ggplot(unit_occurrence_by_position[unit_occurrence_by_position$source == grp, ]) +
      geom_rect(
        data = rect_boundaries[rect_boundaries$source == grp, ],
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax
        ),
        fill = "black"
      ) +
      geom_bar(
        aes(
          x = position,
          y = proportion,
          fill = seq,
          group = proportion
        ),
        stat = "identity",
        width = 1,
        position = "stack"
      ) +
      geom_tile(
        data = consensus[consensus$source == grp, ],
        aes(
          x = position,
          y = -1,
          fill = factor(seq)
        ),
        height = 1
      ) +
      # facet_grid(rows = vars(source), scales = 'free', space = 'free', drop = F) +
      # facetscales::facet_grid_sc(rows = vars(source), scales = list(x = scales_x)) +
      # facet_wrap(~source, scales = 'free') +
      # scale_x_discrete(drop = F) +
      scale_x_continuous(
        expand = c(0, 0),
        breaks = scales::pretty_breaks()(rng),
        limits = c(0, max_x_all_groups)
      ) +
      # tried automating the breaks plyr::round_any(max(unit_occurrence_by_position$position)/5, 10))
      # TODO automate the choice of distance between breakpoints
      scale_y_continuous(expand = c(0, 0)) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(3, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_text(
          family = "Courier",
          face = "bold",
          size = 9
        ),
        legend.text = element_text(family = "Courier", size = 9),
        legend.position = "none"
      ) +
      labs(fill = "Repeat units (count)") +
      scale_fill_manual(
        breaks = ranked_units_in_plot[1:length(ranked_units_in_plot)],
        values = ranked_colors_in_plot,
        labels = paste0(
          ranked_units_in_plot[1:length(ranked_units_in_plot)],
          " (",
          unit_counts_calculated_plus_colors_reordered$count[1:length(ranked_units_in_plot)],
          ")"
        ),
        drop = FALSE
      )
    # o <- o + patchwork::plot_spacer() + patchwork::plot_layout(widths= c(relative_widths[grp], 1-relative_widths[grp]))



    return(o)
  }

  plotList <- lapply(seq_groups, plotOneGroup)
  p <- patchwork::wrap_plots(plotList, ncol = 1)
  # ggplot(unit_occurrence_by_position) +
  #   geom_rect(data = rect_boundaries, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'black') +
  #   geom_bar(aes(x = position, y = proportion, fill = seq, group=proportion), stat='identity', width = 1, position = 'stack') +
  #   geom_tile(data = consensus, aes(x = position, y = -31, fill = factor(seq)), height = 30) +
  #   facet_grid(rows = vars(source), scales = 'free', space = 'free', drop = F) +
  #   #lemon::facet_rep_grid(rows = vars(source), scales = 'free', space = 'free', drop = F,
  #  #              repeat.tick.labels = TRUE) +
  #   #facetscales::facet_grid_sc(rows = vars(source), scales = list(x = scales_x)) +
  #   scale_fill_manual(breaks = ranked_units_in_plot[1:length(ranked_units_in_plot)],
  #                     values = ranked_colors_in_plot,
  #                     labels = paste0(ranked_units_in_plot[1:length(ranked_units_in_plot)], ' (', unit_counts_calculated_plus_colors_reordered$count[1:length(ranked_units_in_plot)], ')'),
  #                     drop = FALSE)

  return(p)
}

#' @export
call_consensus_frequency <- function(unitorders_df, unit_frequency_df) {
  # TODO call_consensus_sequence also uses this pivot and rename function
  # Give each position a number in each sequence
  unitorders_df <- as.data.frame(unitorders_df %>% dplyr::group_by(sample) %>% dplyr::mutate(position = rev(1:dplyr::n())))
  # Mutate the df to have each position as a column
  unitorders_wide_df <- as.data.frame(tidyr::pivot_wider(
    unitorders_df[, c("sample", "position", "character")],
    names_from = position,
    values_from = character
  ))
  # Set the row names as the sample column
  row.names(unitorders_wide_df) <- unitorders_wide_df[, 1]
  unitorders_wide_df[, 1] <- NULL

  # For every column, get the most common unit's frequency
  consensus_freq_as_list <- apply(
    unitorders_wide_df,
    2,
    FUN = function(y) {
      calculate_consensus(y, unit_frequency_df = unit_frequency_df)[2]
    }
  )

  return(consensus_freq_as_list)
}

#' @export
calculate_variability_score <- function(unitorders_df) {
  # Give each position a number in each sequence
  unitorders_with_position_df <- as.data.frame(
    unitorders_df %>% dplyr::group_by(sample) %>% dplyr::mutate(position = 1:dplyr::n())
  )
  # Mutate the df to have each position as a row
  unitorders_with_conservation_score_df <- as.data.frame(
    tidyr::pivot_wider(
      unitorders_with_position_df[, c("sample", "position", "character")],
      names_from = position,
      values_from = character
    )
  )
  row.names(unitorders_with_conservation_score_df) <- unitorders_with_conservation_score_df[, 1]
  unitorders_with_conservation_score_df[, 1] <- NULL
  # Estimate Shannon's entropy per position in MSA
  # where
  # N is number of samples in the alignment
  # K is number of unique units at each position in the alignment
  # x/N is fraction of units of each unit type
  N <- nrow(unitorders_with_conservation_score_df)
  K <- apply(unitorders_with_conservation_score_df, 2, function(x) {
    length(unique(x))
  })
  # #Instead of number of unique units at each position, do number of possible units
  S <- apply(unitorders_with_conservation_score_df, 2, function(y) {
    sum(unlist(lapply(table(y), function(x) {
      x / N * log2(x / N)
    }))) * -1
  })
  S_norm <- S / log2(K)
  # For the NaN values resulting from divide by 0, the calculation is 0/0 so set to 0
  S_norm[is.nan(S_norm)] <- 0

  # Format data
  S_norm_df <- data.frame(variability_score = S_norm)
  S_norm_df$x <- length(S_norm) - as.numeric(row.names(S_norm_df)) + 1
  # Smooth out the data a bit

  # Using a kernel filter
  x <- ksmooth(S_norm_df$x,
    S_norm_df$variability_score,
    "normal",
    bandwidth = 6
  )$x
  y <- ksmooth(S_norm_df$x,
    S_norm_df$variability_score,
    "normal",
    bandwidth = 6
  )$y
  S_norm_df_smoothed <- data.frame(y = y, x = x)

  return(S_norm_df_smoothed)
}

#' @export
summarize_variable_region <- function(vr_unitorders_df, positions_list) {
  # Define the length of the most common unit
  unit_length <- nchar(vr_unitorders_df[grepl("^[GATC]", vr_unitorders_df$seq), "seq"][1])
  # Replace non-30mers and infrequent 30mers with two different distinguishing sequences
  vr_unitorders_df[grepl("frequency", vr_unitorders_df$seq, fixed = T) &
    grepl("^N", vr_unitorders_df$seq), "seq"] <- strrep("N", unit_length)
  vr_unitorders_df[grepl("frequency", vr_unitorders_df$seq, fixed = T) &
    !grepl("^N", vr_unitorders_df$seq), "seq"] <- strrep("n", unit_length)
  # Reverse unit sequences
  vr_unitorders_df$seq_rev <- stringi::stri_reverse(vr_unitorders_df$seq)
  # Define function to reverse a each row string in a df col
  strReverse <- function(x) {
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
  }

  # First summarize alleles
  alleles_at_variable_region_df <- vr_unitorders_df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarize(
      allele_ascii_rev = paste(character, collapse = ""),
      allele_seq_rev = paste(seq_rev, collapse = ""),
      allele = strReverse(allele_ascii_rev),
      allele_num = factor(allele),
      allele_seq = strReverse(allele_seq_rev),
      group = dplyr::first(group)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(allele_num, allele, allele_seq) %>%
    dplyr::summarize(
      count = dplyr::n(),
      most_common_group = names(which.max(table(group))),
      samples = paste(sample, collapse = ",")
    )

  # Alleles by group
  alleles_by_group_df <- vr_unitorders_df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarize(
      allele_ascii_rev = paste(character, collapse = ""),
      allele_seq_rev = paste(seq_rev, collapse = ""),
      allele = strReverse(allele_ascii_rev),
      allele_num = factor(allele),
      allele_seq = strReverse(allele_seq_rev),
      group = dplyr::first(group)
    ) %>%
    ungroup() %>%
    dplyr::group_by(group, allele_num, allele, allele_seq) %>%
    dplyr::summarize(num_sequences = n_distinct(sample))

  # Define most common alleles
  n <- length(alleles_at_variable_region_df$count)
  num_alleles_tied_highest_count <- dim(alleles_at_variable_region_df[alleles_at_variable_region_df$count == sort(alleles_at_variable_region_df$count, partial = n - 1)[n], ])[1]
  if (num_alleles_tied_highest_count > 1) {
    most_common_allele_seq <- alleles_at_variable_region_df[alleles_at_variable_region_df$count == sort(alleles_at_variable_region_df$count, partial = n - 1)[n], ][1, "allele_seq"]$allele_seq
    second_most_common_allele_seq <- alleles_at_variable_region_df[alleles_at_variable_region_df$count == sort(alleles_at_variable_region_df$count, partial = n - 1)[n], ][2, "allele_seq"]$allele_seq

    most_common_allele <- alleles_at_variable_region_df[alleles_at_variable_region_df$count == sort(alleles_at_variable_region_df$count, partial = n - 1)[n], ][1, "allele"]$allele
    second_most_common_allele <- alleles_at_variable_region_df[alleles_at_variable_region_df$count == sort(alleles_at_variable_region_df$count, partial = n - 1)[n], ][2, "allele"]$allele

    warning(
      paste(
        "There are",
        num_alleles_tied_highest_count,
        "alleles tied for highest count. Order of ties depends on their original ordering, which is none (random). Choosing",
        "first as most common allele and second as second as second most common allele"
      )
    )
  } else if (num_alleles_tied_highest_count == 1) {
    most_common_allele_seq <- alleles_at_variable_region_df[alleles_at_variable_region_df$count
    == sort(alleles_at_variable_region_df$count, partial = n - 1)[n], "allele_seq"]$allele_seq
    second_most_common_allele_seq <- alleles_at_variable_region_df[alleles_at_variable_region_df$count
    == sort(alleles_at_variable_region_df$count, partial = n - 1)[n - 1], "allele_seq"]$allele_seq
    most_common_allele <- alleles_at_variable_region_df[alleles_at_variable_region_df$count == sort(alleles_at_variable_region_df$count, partial = n - 1)[n], "allele"]$allele
    second_most_common_allele <- alleles_at_variable_region_df[alleles_at_variable_region_df$count == sort(alleles_at_variable_region_df$count, partial = n - 1)[n - 1], "allele"]$allele
  }

  # Breakout most common alleles into 30bp segments in unitorders format by position
  most_common_allele_unitorders_df <-
    data.frame(
      seq = substring(
        most_common_allele_seq,
        seq(1, nchar(most_common_allele_seq) - 1, unit_length),
        seq(unit_length, nchar(most_common_allele_seq), unit_length)
      ),
      position = positions_list
    )
  second_most_common_allele_unitorders_df <-
    data.frame(
      seq = substring(
        second_most_common_allele_seq,
        seq(1, nchar(second_most_common_allele_seq) - 1, unit_length),
        seq(
          unit_length,
          nchar(second_most_common_allele_seq),
          unit_length
        )
      ),
      position = positions_list
    )

  dflist <- list()

  # Calculate edit distance in nt
  for (i in 1:length(positions_list)) {
    pos <- positions_list[i]
    most_common_unit_at_position <- most_common_allele_unitorders_df[most_common_allele_unitorders_df$position == pos, "seq"]
    second_most_common_unit_at_position <- second_most_common_allele_unitorders_df[second_most_common_allele_unitorders_df$position == pos, "seq"]
    vr_unitorders_edit_dist_single_position_df <-
      vr_unitorders_df %>%
      dplyr::filter(position == pos) %>%
      dplyr::mutate(
        edit_dist_most_common_allele_nt = purrr::map_int(
          seq,
          ~ DescTools::StrDist(., most_common_unit_at_position, method = "hamming")
        ),
        edit_dist_second_most_common_allele_nt = purrr::map_int(
          seq,
          ~ DescTools::StrDist(., second_most_common_unit_at_position, method = "hamming")
        )
      )

    # # If either most common or second common unit is a gap, we don't need to include it's edit dist in the sum
    # if(most_common_unit_at_position == strrep('-', unit_length)){
    #   vr_unitorders_edit_dist_single_position_df$edit_dist_most_common_allele_nt <- 0
    # }
    # if(second_most_common_unit_at_position == strrep('-', unit_length)){
    #   vr_unitorders_edit_dist_single_position_df$edit_dist_second_most_common_allele_nt <- 0
    # }
    # TODO change variable names because it seems like I'm commenting on the frequency of the unit
    # I intended it to read like "most common, unit at position" but it is unclear

    dflist[[i]] <- vr_unitorders_edit_dist_single_position_df
  }

  vr_unitorders_edit_dist_df <- do.call(rbind, dflist)
  # The for loop above and df concatenation changes the default order of the rows
  vr_unitorders_edit_dist_df <- vr_unitorders_edit_dist_df %>%
    group_by(sample) %>%
    arrange(desc(position)) %>%
    dplyr::ungroup()

  # If unit is an infrequent unit, add fixed value to edit distance
  vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep("N", unit_length), "edit_dist_most_common_allele_nt"] <- 3
  vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep("N", unit_length), "edit_dist_second_most_common_allele_nt"] <- 3
  vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep("n", unit_length), "edit_dist_most_common_allele_nt"] <- 2
  vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep("n", unit_length), "edit_dist_second_most_common_allele_nt"] <- 2
  # We only want to compare positions where both the most common allele and each distinct allele have sequence, not gaps
  # If one or the other has a gap, we can see that in the allele graphic
  # If both have gaps, the edit dist will be 0
  # vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep('-', unit_length), 'edit_dist_most_common_allele_nt'] <- 0
  # vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep('-', unit_length), 'edit_dist_second_most_common_allele_nt'] <- 0

  # Summarize alleles again and add up edit dist cols
  # Important: after the
  alleles_by_sample_df <- vr_unitorders_edit_dist_df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarize(
      allele_ascii_rev = paste(character, collapse = ""),
      allele_seq_rev = paste(seq_rev, collapse = ""),
      allele = strReverse(allele_ascii_rev),
      allele_num = factor(allele),
      allele_seq = strReverse(allele_seq_rev),
      group = dplyr::first(group),
      edit_dist_most_common_allele_nt_sum = sum(edit_dist_most_common_allele_nt),
      edit_dist_second_most_common_allele_nt_sum = sum(edit_dist_second_most_common_allele_nt)
    ) %>%
    dplyr::ungroup()

  alleles_at_variable_region_edit_dist_df <- alleles_by_sample_df %>%
    dplyr::group_by(allele_num, allele, allele_ascii_rev) %>%
    dplyr::summarize(
      count = dplyr::n(),
      most_common_group = names(which.max(table(group))),
      samples = paste(sample, collapse = ","),
      edit_dist_most_common_allele_nt = dplyr::first(edit_dist_most_common_allele_nt_sum),
      edit_dist_second_most_common_allele_nt = dplyr::first(edit_dist_second_most_common_allele_nt_sum)
    )

  # Also calculate edit distance between alleles in units
  alleles_at_variable_region_edit_dist_df$edit_dist_most_common_allele_units <- lapply(alleles_at_variable_region_edit_dist_df$allele, function(x) {
    DescTools::StrDist(x, most_common_allele, method = "hamming")
  })
  alleles_at_variable_region_edit_dist_df$edit_dist_most_common_allele_units <- as.numeric(alleles_at_variable_region_edit_dist_df$edit_dist_most_common_allele_units)
  alleles_at_variable_region_edit_dist_df$edit_dist_second_most_common_allele_units <- lapply(alleles_at_variable_region_edit_dist_df$allele, function(x) {
    DescTools::StrDist(x, second_most_common_allele, method = "hamming")
  })
  alleles_at_variable_region_edit_dist_df$edit_dist_second_most_common_allele_units <- as.numeric(
    alleles_at_variable_region_edit_dist_df$edit_dist_second_most_common_allele_units
  )

  # Report edit distance relative to allele length
  allele_length_unit <- nchar(most_common_allele)
  allele_length_nt <- nchar(most_common_allele_seq)
  alleles_at_variable_region_edit_dist_df$edit_dist_most_common_allele_units_frequency <- alleles_at_variable_region_edit_dist_df$edit_dist_most_common_allele_units / allele_length_unit
  alleles_at_variable_region_edit_dist_df$edit_dist_most_common_allele_nt_frequency <- alleles_at_variable_region_edit_dist_df$edit_dist_most_common_allele_nt / allele_length_nt
  alleles_at_variable_region_edit_dist_df$edit_dist_second_most_common_allele_units_frequency <- alleles_at_variable_region_edit_dist_df$edit_dist_second_most_common_allele_units / allele_length_unit
  alleles_at_variable_region_edit_dist_df$edit_dist_second_most_common_allele_nt_frequency <- alleles_at_variable_region_edit_dist_df$edit_dist_second_most_common_allele_nt / allele_length_nt

  return(list(
    alleles_at_variable_region_edit_dist_df,
    alleles_by_group_df
  ))
}

#' @export
summarize_variable_regions <- function(unitorders_df, threshold) {
  threshold_variability <- threshold
  unitorders_df$seq[unitorders_df$seq == ""] <- strrep("-", 30)
  variability_score_df <- calculate_variability_score(unitorders_df)
  print(variability_score_df)
  # Use a copy of the x position column to mark where variability score is larger than threshold, mark with 1
  variability_score_df$detect_runs <- variability_score_df$x
  variability_score_df[variability_score_df$y > threshold_variability, "detect_runs"] <- 1
  # Detect length of runs, doesn't specifically know to detect 1s, this is relying on the assumption that two neighboring positions will not have the exact same variability score
  runs <- rle(variability_score_df$detect_runs)
  # Summarize lengths of runs in df
  end <- cumsum(runs$lengths)
  print("END: ")
  print(end)
  start <- c(1, dplyr::lag(end)[-1] + 1)
  print("START: ")
  print(start)
  run_lengths_df <- data.frame(start, end)
  print("RUN LENGTHS: ")
  print(run_lengths_df)
  variable_regions_boundaries_df <- subset(run_lengths_df, end - start >= 2)
  print("Variable regions boundaries df")
  print(variable_regions_boundaries_df)
  variable_region_positions <- apply(variable_regions_boundaries_df, 1, function(x) {
    seq(x["start"], x["end"])
  }) # the names of this nested list are the row numbers where the range is defined
  # Give each position in the sequence a number
  print(variable_region_positions)
  unitorders_with_position_df <- as.data.frame(unitorders_df %>% dplyr::group_by(sample) %>% dplyr::mutate(position = rev(1:dplyr::n())))
  variable_region_unitorders_df <- lapply(seq(variable_region_positions), function(x) {
    unitorders_with_position_df[unitorders_with_position_df$position %in% variable_region_positions[[x]], ]
  })
  # This is a nested list of dfs and there are several modifications I want to make to every df in the list
  # Define a function to do all these modifications
  alleles_at_variable_regions_df <- mapply(
    summarize_variable_region,
    variable_region_unitorders_df,
    positions_list = variable_region_positions,
    SIMPLIFY = FALSE
  )

  # return(alleles_at_variable_regions_df)
  return(
    list("allele_dfs" = alleles_at_variable_regions_df, "position_lists" = variable_region_positions)
  )
}

#' @export
plot_alleles_at_variable_regions <- function(allele_at_variable_regions_df,
                                             unit2ascii_df,
                                             unit2color_df) {
  # #Order by edit dist to most common, this is commented out because df is ordered by allele_num
  # allele_at_variable_regions_df <- allele_at_variable_regions_df[order(allele_at_variable_regions_df$edit_dist_most_common_allele_units),]
  # allele_at_variable_regions_df$allele_num <- forcats::fct_reorder(allele_at_variable_regions_df$allele_num, allele_at_variable_regions_df$edit_dist_most_common_allele_units)

  # Columns for plotting format: character, seq, gp, allele
  # allele_num is a factored version of allele
  # Drop samples column
  unitorders_variable_region_df <- subset(allele_at_variable_regions_df, select = -c(samples)) %>%
    tidyr::separate_rows(allele_ascii_rev, sep = "(?!^)")
  # Drop any row that has a blank in character column (rename character column)
  colnames(unitorders_variable_region_df) <- c(
    "allele_num",
    "allele",
    "character",
    "count",
    "most_common_group",
    "edit_dist_most_common_allele_nt",
    "edit_dist_second_most_common_allele_nt",
    "edit_dist_most_common_allele_units",
    "edit_dist_second_most_common_allele_units"
  )
  # Merge df with unit sequences df
  unitorders_variable_region_df <- plyr::join(unitorders_variable_region_df, unit2ascii_df)
  # Add width of each unit
  unitorders_variable_region_df$width <- 2
  unitorders_variable_region_df$gp <- seq(nrow(unitorders_variable_region_df))
  unitorders_variable_region_df[is.na(unitorders_variable_region_df$seq), "seq"] <- "-"

  unitorders_variable_region_df$allele_num <- factor(unitorders_variable_region_df$allele_num,
    levels = rev(levels(
      allele_at_variable_regions_df$allele_num
    ))
  )

  vr_plot <- plot_msa(
    unitorders_variable_region_df,
    unit2ascii_df,
    unit2color_df,
    plotVariableRegion = T
  )

  return(vr_plot)
}
