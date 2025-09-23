#' positions_list
#'
#' Example of multiple sequence alignment positions defining a variable region.
#'
#'
#' @format A vector of sequential integers.
#'
#' @source Second element of `summarize_variable_regions` output..
"positions_list"

#' unit_frequency_example
#'
#' Table of repeat unit frequencies with no grouping of infrequent repeat units.
#'
#' A `data.frame` describing the frequency of each repeat unit in the dataset.
#'
#'
#' @format A data frame with 158 rows and 7 variables.
#'
#' \describe{
#' \item{seq}{Sequence of repeat unit (character)}
#' \item{length}{Length in nucleotides of repeat unit (integer)}
#' \item{count}{Count of repeat unit observed in all sequences (integer)}
#' \item{Frequency}{Frequency of repeat unit abundance among all repeat units (numeric)}
#' \item{samples}{Comma-separated list of samples with at least one copy of
#'                repeat unit (character)}
#' \item{num_samples}{Count of samples with at least one copy of repeat unit (integer)}
#' \item{frac_samples}{Fraction of total samples with at least one copy of the
#'                     repeat unit (numeric)}
#' }
#'
#' @source Custom .py script to parse repeat units found in dataset.
"unit_frequency_example"

#' unit2ascii
#'
#' Conversion table for repeat units represented by single characters.
#'
#' A `data.frame` describing the conversion of a repeat unit to a single character
#' for alignment.
#'
#' @format A data frame with 37 rows and 2 variables.
#'
#' \describe{
#' \item{seq}{Sequence of repeat unit (character)}
#' \item{character}{Repeat unit encoded as a single ASCII character (character)}
#' }
#'
#' @source Custom .py script to parse repeat units found in dataset.
"unit2ascii"

#' unit2color
#'
#' Conversion table for repeat units represented by unique color.
#'
#' A `data.frame` describing the color to use for each repeat unit in plots.
#' Black is used for an alignment gap ("-").
#'
#' @format A data frame with 38 rows and 2 variables.
#'
#' \describe{
#' \item{seq}{Sequence of repeat unit (character)}
#' \item{color}{Hex color (character)}
#' }
#'
#' @source Custom .py script to parse repeat units found in dataset.
"unit2color"

#' unitorders_example
#'
#' Example of a parsed multiple sequence alignment.
#'
#' A `data.frame` with the results of a multiple sequence alignment of 3 repetitive sequences.
#' Each row is one repeat unit. Within a sample, rows are ordered from bottom to top in the 5' to 3' direction.
#'
#'
#' @format A data frame with 600 rows and 6 variables.
#'
#' \describe{
#' \item{character}{Repeat unit encoded as a single ASCII character (character)}
#' \item{sample}{Sequence ID (character)}
#' \item{count}{Width of each repeat unit on the plot (numeric)}
#' \item{seq}{Sequence of repeat unit (character)}
#' \item{gp}{Integer that increments one for every row in the ungrouped `data.frame`
#'           This gives every repeat unit in a sequence its own position, rather than clustering identical units together (integer)}
#' }
#'
#' @source Custom .py script to parse FASTA file of multiple sequence alignment.
"unitorders_example"

#' unitorders_group1_example
#'
#' Example of a parsed multiple sequence alignment for one group of sequences.
#'
#' A `data.frame` with the results of a multiple sequence alignment of 134 repetitive sequences.
#' Each row is one repeat unit. Within a sample, rows are ordered from bottom to top in the 5' to 3' direction.
#'
#'
#' @format A data frame with 33634 rows and 5 variables.
#'
#' \describe{
#' \item{character}{Repeat unit encoded as a single ASCII character (character)}
#' \item{sample}{Sequence ID (character)}
#' \item{count}{Width of each repeat unit on the plot (numeric)}
#' \item{seq}{Sequence of repeat unit (character)}
#' \item{gp}{Integer that increments one for every row in the ungrouped `data.frame`
#'           This gives every repeat unit in a sequence its own position,
#'           rather than clustering identical units together (integer)}
#' }
#'
#' @source Custom .py script to parse FASTA file of multiple sequence alignment.
"unitorders_group1_example"

#' unitorders_group2_example
#'
#' Example of a parsed multiple sequence alignment for one group of sequences.
#'
#' A `data.frame` with the results of a multiple sequence alignment of 8 repetitive sequences.
#' Each row is one repeat unit. Within a sample, rows are ordered from bottom to top in the 5' to 3' direction.
#'
#'
#' @format A data frame with 2712 rows and 5 variables.
#'
#' \describe{
#' \item{character}{Repeat unit encoded as a single ASCII character (character)}
#' \item{sample}{Sequence ID (character)}
#' \item{count}{Width of each repeat unit on the plot (numeric)}
#' \item{seq}{Sequence of repeat unit (character)}
#' \item{gp}{Integer that increments one for every row in the ungrouped `data.frame`
#'           This gives every repeat unit in a sequence its own position,
#'           rather than clustering identical units together (integer)}
#' }
#'
#' @source Custom .py script to parse FASTA file of multiple sequence alignment.
"unitorders_group2_example"

#' unitorders_withSubgroups_example
#'
#' Example of a parsed multiple sequence alignment.
#'
#' A `data.frame` with the results of a multiple sequence alignment of 3 repetitive sequences.
#' Each row is one repeat unit. Within a sample, rows are ordered from bottom to top in the 5' to 3' direction.
#'
#'
#' @format A data frame with 47171 rows and 7 variables.
#'
#' \describe{
#' \item{character}{Repeat unit encoded as a single ASCII character (character)}
#' \item{sample}{Sequence ID (character)}
#' \item{count}{Width of each repeat unit on the plot (numeric)}
#' \item{seq}{Sequence of repeat unit (character)}
#' \item{gp}{Integer that increments one for every row in the ungrouped `data.frame`
#'           This gives every repeat unit in a sequence its own position,
#'           rather than clustering identical units together (integer)}
#' \item{group}{Used with plotSubgroups = T to cluster sequences in the plot.
#'              Integer indicating multiple sequence alignment membership.}
#' \item{position}{Position in multiple sequence alignment.}
#' }
#'
#' @source Custom .py script to parse FASTA file of multiple sequence alignment.
"unitorders_withSubgroups_example"

#' unitorders_variableregion_example
#'
#' Subset of a parsed multiple sequence alignment at identified variable regions.
#'
#' A `data.frame` with 134 aligned sequences at one variable region.
#' Each row is one repeat unit. Within a sample, rows are ordered from bottom to top in the 5' to 3' direction.
#'
#'
#' @format A data frame with 2278 rows and 7 variables.
#'
#' \describe{
#' \item{character}{Repeat unit encoded as a single ASCII character (character)}
#' \item{sample}{Sequence ID (character)}
#' \item{count}{Width of each repeat unit on the plot (numeric)}
#' \item{seq}{Sequence of repeat unit (character)}
#' \item{gp}{Integer that increments one for every row in the ungrouped `data.frame`
#'           This gives every repeat unit in a sequence its own position, rather than clustering identical units together (integer)}
#' \item{group}{Integer indicating multiple sequence alignment membership.}
#' \item{position}{Position within multiple sequence alignment.}
#' }
#'
#' @source Custom .py script to parse FASTA file of multiple sequence alignment.
"unitorders_variableregion_example"


#' vr_alleles_summarized_example
#'
#' Subset of a parsed multiple sequence alignment at identified variable regions.
#'
#' A `data.frame` with counts of the unique sequences found at one variable region.
#' Each row is one variable region sequence. The sequences can be clustered into
#' two alleles by edit distances.
#'
#'
#' @format A data frame with 19 rows and 14 variables.
#'
#' \describe{
#' \item{allele_num}{Variable region sequence (factor)}
#' \item{allele}{Variable region sequence (character)}
#' \item{allele_ascii_rev}{Variable region sequence, reversed (character)}
#' \item{count}{Number of sequences with variable region sequence (numeric)}
#' \item{most_common_group}{Group of sequences most commonly identified in (integer)}
#' \item{samples}{Sequence IDs with variable region sequence (character)}
#' \item{edit_dist_most_common_allele_nt}{Edit distance in nucleotides to
#'                                        most common variable region sequence}
#' \item{edit_dist_second_most_common_allele_nt}{Edit distance in nucleotides to
#'                                               second most common variable region sequence}
#' \item{edit_dist_most_common_allele_units}{Edit distance in repeat units to
#'                                            most common variable region sequence}
#' \item{edit_dist_second_most_common_allele_units}{Edit distance in repeat units to
#'                                                 second most common variable region sequence}
#' \item{edit_dist_most_common_allele_units_frequency}{Edit distance in repeat units per
#'                                                     length of variable region sequence}
#' \item{edit_dist_most_common_allele_nt_frequency}{Edit distance in nucleotides per
#'                                                  length of variable region sequence}
#' \item{edit_dist_second_most_common_allele_units_frequency}{Edit distance in repeat units per
#'                                                            length of variable region sequence}
#' \item{edit_dist_second_most_common_allele_nt_frequency}{Edit distance in nucleotides per
#'                                                         length of variable region sequence}
#' }
#'
#' @source First element of `summarize_variable_regions` output.
"vr_alleles_summarized_example"
