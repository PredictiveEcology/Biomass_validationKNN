## function to replace NAs in data.table columns
replaceNAs <- function(x, val = 0) {
  x[is.na(x)] <- val
  x
}
