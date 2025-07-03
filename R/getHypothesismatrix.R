# Function: get_extended_matrix
# Purpose: Extend a matrix V (with fewer than p columns) to full rank by adding columns from the identity matrix
get_extended_matrix <- function(V) {
  p <- nrow(V)  # Number of rows of V (assumed to be the dimension of the parameter space)
  
  # If V is square, it's already p x p, and extending it makes no sense in the model
  if(p == ncol(V)){
    stop("This is not part of the model")  # Error: V should not already be square
  }
  
  unit_vec <- diag(p)  # Identity matrix of size p (used to extract standard basis vectors)
  i <- 1  # Index for standard basis vectors
  
  rk <- Matrix::rankMatrix(V)  # Compute the current rank of V
  
  # Keep adding standard basis vectors as columns until full rank is reached
  while(rk < p){
    V_new <- cbind(V, unit_vec[,i])  # Add the i-th standard basis vector to V
    
    # Only accept the new column if it increases the rank
    if(Matrix::rankMatrix(V_new) > rk){
      V <- V_new           # Update V
      rk <- Matrix::rankMatrix(V_new)  # Update rank
    }
    
    i <- i + 1  # Move to the next basis vector
  }
  
  return(V)  # Return the extended matrix
}


# Function: get_hypothesis
# Purpose: Construct a hypothesis matrix C and a hypothesis vector zeta from a given vector v0 and matrix V,
# containing the vectorized matrices from the linear covariance structure modell
get_hypothesis <- function(v0, V){
  p <- nrow(V)  # Number of parameters (rows in V)
  q <- ncol(V)  # Number of constraints or model dimensions (columns in V)
  
  V_ext <- get_extended_matrix(V)  # Extend V to full rank (if not already)
  
  # Select the rows (q+1 to p) of the identity matrix to build the contrast matrix E
  # This corresponds to the "unexplained" part of the parameter space
  E <- diag(1, p, p)[(q + 1):p, ]
  
  # Compute the hypothesis matrix C by transforming E with the inverse of the extended V
  C <-  E %*% solve(V_ext)
  
  # Error handling: v0 must be of the same dimension as the number of parameters
  if(length(v0) != p){
    stop("v0 must have the same length as p")
  }
  
  zeta <- C %*% v0  # Compute the transformed (reduced) hypothesis vector
  
  # Return the hypothesis matrix and vector
  return(list("Hypothesenmatrix" = C,  # "C" matrix in the hypothesis H0: Cθ = zeta
              "Hypothesenvector" = zeta))
}