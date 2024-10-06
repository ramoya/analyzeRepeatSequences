#' Example of a parsed multiple sequence alignment.
#'
#' A `data.frame` with the results of a multiple sequence alignment of 3 repetitive sequences.
#' Each row is one repeat unit. Within a sample, rows are ordered from bottom to top in the 5' to 3' direction.
#'
#'
#' @format Column name (type)
#'
#' \describe{
#' \item{character}{Repeat unit encoded as a single ASCII character (character)}
#' \item{sample}{Sequence ID (character)}
#' \item{count}{Width of each repeat unit on the plot (numeric)}
#' \item{seq}{Sequence of repeat unit (character)}
#' \item{gp}{Integer that increments one for every row in the ungrouped `data.frame`
#'           This gives every repeat unit in a sequence its own position, rather than clustering identical units together (integer)}
#' \item{group}{Used with plotSubgroups = T to cluster sequences in the plot}
#' }
"unitorders_example"
