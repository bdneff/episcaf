mutate6 <- function(seq, replacement, n = 1, max_tries = 5, fail = c("error", "na")) {

  fail <- match.arg(fail)

  stopifnot(
    is.character(seq),
    length(seq) == 1,
    is.character(replacement),
    nchar(replacement) == 6,
    replacement == toupper(replacement),
    is.numeric(n),
    n >= 1,
    is.numeric(max_tries),
    max_tries >= 1
  )

  for (iter in seq_len(n)) {

    success <- FALSE

    for (attempt in seq_len(max_tries)) {

      chars <- strsplit(seq, "")[[1]]
      len <- length(chars)

      is_upper <- chars == toupper(chars) & chars != tolower(chars)
      is_lower <- chars == tolower(chars) & chars != toupper(chars)

      capital_idx <- which(is_upper)
      starts <- 1:(len - 5)

      valid_starts <- starts[vapply(starts, function(i) {
        six_idx <- i:(i + 5)

        # Condition 1: all six must be lowercase
        if (!all(is_lower[six_idx])) return(FALSE)

        # Condition 2: every position of the 6-mer must be >4 away from any capital
        if (length(capital_idx) > 0) {
          if (any(abs(outer(six_idx, capital_idx, "-")) <= 4)) {
            return(FALSE)
          }
        }

        TRUE
      }, logical(1))]

      if (length(valid_starts) > 0) {
        start <- sample(valid_starts, 1)
        chars[start:(start + 5)] <- strsplit(replacement, "")[[1]]
        seq <- paste(chars, collapse = "")
        success <- TRUE
        break
      }
    }

    if (!success) {
      msg <- sprintf(
        "Iteration %d failed after %d attempts: no valid lowercase 6-mer found.",
        iter, max_tries
      )
      if (fail == "error") stop(msg) else return(NA_character_)
    }
  }

  seq
}
