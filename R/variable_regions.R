#' Calculate variability score
#'
#' @description Calculates a Shannon's uncertainty score for each position in
#' the alignment to quantify variability. Variability is computed as:
#'
#' {\eqn{H(X)=-\sum_{i=1}^{K}p(x)_{i}log_{2}(p(x)_{i})}}
#'
#' Where:
#'
#' {\eqn{p(x)_{i}}} is the fraction of repeat units of unit type {\eqn{i}}.
#'
#' {\eqn{K}} is the number of different repeat units at position {\eqn{X}}.
#'
#' Gaps are included in the calculation.
#'
#' Last, the output value is normalized using {\eqn{H(X)_{normalized}=\frac{H(X)}{log_{2}(K)}}}
#' then a smoothing filter is applied using the R package ksmooth.
#'
#' @param unitorders_df         `data.frame` of a parsed multiple sequence
#'   alignment. Has columns `character`, `sample`, `count`, `seq`, `gp`. See
#'   dataset `unitorders_example` for format.
#'
#' @returns a `data.frame` summarizing the distinct sequences with columns `y` (variability scores) and
#'   `x` (alignment positions)
#'
#' @examples
#' calculate_variability_score(unitorders_example)
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
  # K is number of possible units at each position in the alignment
  # x/N is fraction of units of each unit type
  N <- nrow(unitorders_with_conservation_score_df)
  K <- apply(unitorders_with_conservation_score_df, 2, function(x) {
    length(unique(x))
  })
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

  # Smooth out the data using a kernel filter
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


#' Summarize alleles at variable region
#'
#' @description Called by `summarize_variable_regions`.
#'   Defines the two most common sequences as two alleles within a
#'   variable region. Calculates the edit distance between these alleles and the
#'   other sequences within the variable region. With this, users can set edit
#'   distance thresholds to assign each sequence to an allele.
#'
#'   Note: This function defines exactly two alleles. There may be more
#'   alleles that have high edit distances to both.
#'
#'
#' @param vr_unitorders_df         `data.frame` of a parsed multiple sequence
#'   alignment within the positions defined as a variable region by
#'   `summarize_variable_regions`. Has columns `character`, `sample`, `count`,
#'   `seq`, `gp`, `group`, `position`. See dataset
#'   `unitorders_variableregion_example` for format.
#' @param positions_list            Vector of positions defined as a variable region. See dataset `positions_list` for format.
#'
#' @returns a `data.frame` summarizing the unique sequences of a variable region.
#'
#' @examples
#' summarize_variable_region(unitorders_variableregion_example, positions_list)
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
  # Define function to reverse each string in a column
  strReverse <- function(x) {
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
  }

  if(!"group" %in% names(vr_unitorders_df)){
    vr_unitorders_df$group <- 1
  }

  # Summarize alleles
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

  # Define two most common alleles
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

  # Split most common allele sequences into repeat unit-length segments. Format by position like unitorders_df.
  most_common_allele_unitorders_df <-
    data.frame(seq = substring(
      most_common_allele_seq,
      seq(1, nchar(most_common_allele_seq) - 1, unit_length),
      seq(unit_length, nchar(most_common_allele_seq), unit_length)
    ),
    position = positions_list)
  second_most_common_allele_unitorders_df <-
    data.frame(seq = substring(
      second_most_common_allele_seq,
      seq(1, nchar(second_most_common_allele_seq) - 1, unit_length),
      seq(
        unit_length,
        nchar(second_most_common_allele_seq),
        unit_length
      )
    ),
    position = positions_list)

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

    dflist[[i]] <- vr_unitorders_edit_dist_single_position_df
  }

  vr_unitorders_edit_dist_df <- do.call(rbind, dflist)
  # The for loop above and rbind change the order of the rows. Reset it.
  vr_unitorders_edit_dist_df <- vr_unitorders_edit_dist_df %>%
    dplyr::group_by(sample) %>%
    dplyr::arrange(desc(position)) %>%
    dplyr::ungroup()

  # If unit is an infrequent unit, add fixed value to edit distance
  vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep("N", unit_length), "edit_dist_most_common_allele_nt"] <- 3
  vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep("N", unit_length), "edit_dist_second_most_common_allele_nt"] <- 3
  vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep("n", unit_length), "edit_dist_most_common_allele_nt"] <- 2
  vr_unitorders_edit_dist_df[vr_unitorders_edit_dist_df$seq == strrep("n", unit_length), "edit_dist_second_most_common_allele_nt"] <- 2

  # Summarize alleles again and sum edit distances
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

  return(alleles_at_variable_region_edit_dist_df)
}

#' Summarize alleles at all variable regions within a multiple sequence alignment
#'
#' @description Defines variable regions by variability score. Formats the input
#'   variable region sequences to `summarize_variable_region`.
#'
#'
#' @param unitorders_df         `data.frame` of a parsed multiple sequence
#'   alignment. Has columns `character`, `sample`, `count`, `seq`, `gp`. See
#'   dataset `unitorders_example` for format.
#' @param threshold            Threshold for variability score (range 0-1) to use for defining variable regions.
#'
#' @returns Named list with two elements.
#'   allele_dfs: List of `data.frame`, one per variable region. Each row is a unique sequence of a variable region.
#'   positions_lists: List of vectors defining the positions of the multiple sequences alignment that are variable regions.
#' @examples
#' summarize_variable_regions(unitorders_example, threshold = 0.25)
#' @export
summarize_variable_regions <- function(unitorders_df, threshold) {
  threshold_variability <- threshold

  unitorders_df$seq[unitorders_df$seq == ""] <- strrep("-", 30)

  variability_score_df <- calculate_variability_score(unitorders_df)

  # Define variable regions
  # Use a copy of the existing position column to mark where variability score is larger than threshold. Mark with 1.
  variability_score_df$detect_runs <- variability_score_df$x
  variability_score_df[variability_score_df$y > threshold_variability, "detect_runs"] <- 1

  # Detect length of runs.
  runs <- rle(variability_score_df$detect_runs)
  # Summarize lengths of runs.
  end <- cumsum(runs$lengths)
  start <- c(1, dplyr::lag(end)[-1] + 1)
  run_lengths_df <- data.frame(start, end)
  variable_regions_boundaries_df <- subset(run_lengths_df, end - start >= 2)
  variable_region_positions <- apply(variable_regions_boundaries_df, 1, function(x) {
    seq(x["start"], x["end"])
  }) # the names of this nested list are the row numbers where the range is defined

  # Assign position
  unitorders_with_position_df <- as.data.frame(unitorders_df %>% dplyr::group_by(sample) %>% dplyr::mutate(position = rev(1:dplyr::n())))
  #Subset to variable region positions
  variable_region_unitorders_df <- lapply(seq(variable_region_positions), function(x) {
    unitorders_with_position_df[unitorders_with_position_df$position %in% variable_region_positions[[x]], ]
  }) #This is a nested list of dfs

  #For each variable region, summarize its alleles
  alleles_at_variable_regions_df <- mapply(
    summarize_variable_region,
    variable_region_unitorders_df,
    positions_list = variable_region_positions,
    SIMPLIFY = FALSE
  )

  return(
    list("allele_dfs" = alleles_at_variable_regions_df, "position_lists" = variable_region_positions)
  )
}


#' Plot variable regions
#'
#' @description
#' For visualizing the unique sequences within variable regions identified by
#' `summarize_variable_regions`.
#'
#'
#' @param allele_at_variable_regions_df         `data.frame` of a parsed
#'   multiple sequence alignment. Has columns `character`, `sample`, `count`,
#'   `seq`, `gp`. See dataset `unitorders_example` for format.
#' @param unit2ascii_df                         `data.frame` of the conversion
#'   between repeat units and single ASCII characters. Has columns `seq`,
#'   `character`. See dataset `unit2ascii` for format.
#' @param unit2color_df                         `data.frame` of the conversion
#'   between repeat units and colors. Has columns `seq`, `color`. See dataset
#'   `unit2color` for format.
#'
#' @returns ggplot object showing the alignment of each unique sequence within a
#'   variable region.
#' @examples
#' plot_alleles_at_variable_regions(vr_alleles_summarized_example,
#'                                  unit2ascii_df = unit2ascii,
#'                                  unit2color_df = unit2color)

#' @export
#' @import ggplot2
plot_alleles_at_variable_regions <- function(allele_at_variable_regions_df,
                                             unit2ascii_df,
                                             unit2color_df) {
  # #Order by edit dist to most common, this is commented out because df is ordered by allele_num
  # allele_at_variable_regions_df <- allele_at_variable_regions_df[order(allele_at_variable_regions_df$edit_dist_most_common_allele_units),]
  # allele_at_variable_regions_df$allele_num <- forcats::fct_reorder(allele_at_variable_regions_df$allele_num, allele_at_variable_regions_df$edit_dist_most_common_allele_units)

  # Columns for plotting format: character, seq, gp, allele
  # allele_num is a factored version of allele
  # Drop samples column

  unitorders_variable_region_df <- allele_at_variable_regions_df %>% dplyr::select(!c(samples)) %>%
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

  l <- format_unitorders_to_plot(unitorders_variable_region_df, unit2color_df)
  unitorders_variable_region_df <- l[[1]]
  ranked_units_in_plot <- l[[2]]
  ranked_colors_in_plot <- l[[3]]
  unit_counts_calculated_plus_colors_reordered <- l[[4]]

  vr_plot <- ggplot(
    unitorders_variable_region_df,
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

  return(vr_plot)
}
