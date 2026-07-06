mutateToA <- function(seq, n) {
  stopifnot(
    is.character(seq),
    length(seq) == 1,
    is.numeric(n),
    n >= 1
  )

  chars <- strsplit(seq, "")[[1]]
  is_upper <- chars == toupper(chars) & chars != tolower(chars)

  # Identify runs of contiguous uppercase letters
  r <- rle(is_upper)
  upper_runs <- which(r$values)

  # NEW: return empty string if n exceeds number of chunks
  if (n > length(upper_runs)) {
    return("")
  }

  # Compute start/end indices of each uppercase run
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1

  run_idx <- upper_runs[n]
  idx <- starts[run_idx]:ends[run_idx]

  # Replace odd-numbered positions (1-based within the chunk)
  odd_positions <- idx[seq_along(idx) %% 2 == 1]
  chars[odd_positions] <- "A"

  paste(chars, collapse = "")
}
