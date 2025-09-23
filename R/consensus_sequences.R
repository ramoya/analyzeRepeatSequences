# For each position, handle infrequent alignment gaps and repeat units with equal frequency.
calculate_consensus <- function(x, unit_frequency_df, unit2ascii_df) {
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
      unit2ascii_df,
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

#' Define consensus sequence from multiple sequence alignment `data.frame`.
#'
#' @description
#' `call_consensus_sequence` defines a consensus sequence from a multiple
#'                           sequence alignment represented by `unitorders_df`.
#' @param unitorders_df         `data.frame` of a parsed multiple sequence
#'   alignment. Has columns `character`, `sample`, `count`, `seq`, `gp`. See
#'   dataset `unitorders_example` for format.
#' @param unit_frequency_df     `data.frame`cataloging all repeat units observed in
#'   the dataset and their frequency. Has columns `seq`, `length`, `count`,
#'                                                `Frequency`, `samples`,
#'                                                `num_samples`, `frac_samples`.
#'   See dataset `unit_frequency_example` for format.
#' @param unit2ascii_df         `data.frame` of the conversion
#'   between repeat units and single ASCII characters. Has columns `seq`,
#'   `character`. See dataset `unit2ascii` for format.
#'
#' @returns a list.
#'          First element: the consensus sequence with alignment gaps.
#'          Second element: the supporting multiple sequence alignment
#'                          (rows: individual sequences, columns: alignment position)
#' @examples
#' call_consensus_sequence(unitorders_example, unit_frequency_example, unit2ascii)
#' @export
call_consensus_sequence <- function(unitorders_df, unit_frequency_df, unit2ascii_df) {
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
      calculate_consensus(y, unit_frequency_df = unit_frequency_df,
                          unit2ascii_df = unit2ascii_df)[1]
    }
  )

  # Reverse the order of the units and collapse this list into a string
  consensus_seq_as_string <- paste(rev(consensus_seq_as_list), collapse = "")

  return(list(consensus_seq_as_string, unitorders_wide_df))
}

# Convert a consensus sequence into `unitorders_df` format.
#' @export
#' @keywords internal
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

# Calculate the fractional contribution of each repeat unit at each alignment position.
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

# Plot fractional abundance of a single multiple sequence alignment.
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
}

#' Make summary plot
#'
#' @description
#' `make_repeat_summary_figure` is a plotting function to summarize more than one
#'                              multiple sequence alignment and their consensus sequences.
#'
#'
#' @param unitorders_group         `data.frame` of a parsed multiple sequence
#'   alignment for a group of sequences. Has columns `character`, `sample`, `count`, `seq`, `gp`. See
#'   dataset `unitorders_group1` for format.
#' @param ...                   Optional additional `unitorders_group` `data.frame`s.
#' @param unit2color_df         `data.frame` of the conversion between repeat
#'   units and colors. Has columns `seq`, `color`. See dataset `unit2color` for
#'   format.
#' @param unit2ascii_df         `data.frame` of the conversion
#'   between repeat units and single ASCII characters. Has columns `seq`,
#'   `character`. See dataset `unit2ascii` for format.
#' @param unit_frequency_df     `data.frame`cataloging all repeat units observed in
#'   the dataset and their frequency. Has columns `seq`, `count`, `samples`, `no_samples`, `unit_length`.
#'   See dataset `unit_frequency_example` for format.
#'
#' @returns ggplot object of several multiple sequence alignments on the same x-axis scale, arranged vertically.
#' @examples
#' make_repeat_summary_figure(unitorders_group1_example,
#'                            unitorders_group2_example,
#'                            unit2color_df = unit2color,
#'                            unit2ascii_df = unit2ascii,
#'                            unit_frequency_df = unit_frequency_example)
#' @export
#' @import ggplot2
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
    call_consensus_sequence(x, unit_frequency_df, unit2ascii_df)
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
  )

  # get fractional abundance at each position
  make_fractional_abundance_df <- function(x) {
    tmp <- do.call("rbind", as.list(unlist(
      apply(
        call_consensus_sequence(x, unit_frequency_df, unit2ascii_df)[[2]],
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
      call_consensus_sequence(x, unit_frequency_df, unit2ascii_df)[[2]]
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

  # Calculate new unit counts
  unit_counts_calculated <- consensus %>%
    dplyr::group_by(seq) %>%
    dplyr::summarize(count = dplyr::n())

  # Add colors to unit counts
  unit_counts_calculated_plus_colors <- plyr::join(unit2color_df, unit_counts_calculated) #default left join to retain order
  unit_counts_calculated_plus_colors[is.na(unit_counts_calculated_plus_colors$count), "count"] <- 0
  # Reorder samples
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

    return(o)
  }

  plotList <- lapply(seq_groups, plotOneGroup)
  p <- patchwork::wrap_plots(plotList, ncol = 1)

  return(p)
}

