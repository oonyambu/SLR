

largest_submatrix <- function(arr){
  # """
  # Obtain the largest submatrix
  # """
  
  find_index <- function(arr, idx){
    # """
    # Obtain the row and col slices containing non-na values
    # """
  n <- nrow(arr)
  m <- ncol(arr)
  if (length(idx) == 2){
    rows <- if (idx[1] < n) c(idx[1], idx[1] + 1) else c(idx[1] - 1, idx[1])
    cols <- if (idx[2] < m) c(idx[2], idx[2] + 1) else c(idx[2] - 1, idx[2])
    idx  <- c(rows, cols, 0)
  }
    row1 <- idx[1]
    row2 <- idx[2]
    col1 <- idx[3]
    col2 <- idx[4]
  # move up
  if (row1 > 1 && all(arr[(row1 - 1):row2, col1:col2]))
    row1 <-  row1 - 1 
  # move down
  if (row2 < n && all(arr[row1:(row2+1), col1:col2]))
    row2 <-  row2 + 1 
  # move left
  if (col1 > 1 && all(arr[row1:row2, (col1-1):col2]))
    col1 <- col1 - 1  
  # move right
  if (col2 < m && all(arr[row1:row2, col1:(col2+1)]))
    col2 <-  col2 + 1 
  # any movement? if yes, repeat the process above.
  size = (row1 - row2 + 1) * (col1 - col2 + 1)
  if (any(c(row1, row2, col1, col2) != idx[1:4]))
    find_index(arr, c(row1, row2, col1, col2, size))
  else c(row1, row2, col1, col2, size)
  }
  
  sub_indices <- function(arr_bool, index, vals = NULL){
    # """
    # Compare the slices and obtain the slices with the most values
    # """
  #print(index)
  row = index[,1]
  col = index[,2]
  indices = find_index(arr_bool,c(row[1], col[1]))
  # slices with max region
  vals = if(is.null(vals) || vals[5] < indices[5]) indices else vals
  
  # any non-nan value not within the have not been traversed? Repeat
  if (any(i2 <- !((indices[1] <= row) & (row <= indices[2]) &
                  (indices[3] <= col) & (col <= indices[4]))))
    sub_indices(arr_bool, cbind(row[i2], col[i2]), vals)
  else vals
  }
  
  arr_bool = !is.na(arr)
  vals = sub_indices(arr_bool, which(arr_bool, TRUE))
  arr[vals[1]:vals[2], vals[3]:vals[4]]
}
  
  
  
