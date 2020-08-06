# This program plots a numerical histogram of the eigenvalue density of
# Gaussian random matrix ensembles (Orthogonal, Unitary and Symplectic).
# You will be asked to select the ensemble by providing the value of the
# Dyson beta index (1 for GOE, 2 for GUE, 4 for GSE), to provide the
# matrix size N, and to choose the number of matrices to be diagonalized.

# Reads the value of the Dyson beta index from the Command Window
prompt = "Choose value of beta (1 for GOE, 2 for GUE, 4 for GSE): "
beta = 1 # strtoi(readline(prompt))

# Reads the matrix size from the Command Window
prompt = "Choose matrix size: "
N = 8 # strtoi(readline(prompt))

# Reads the number of matrices to be diagonalized from the Command Window
prompt = "Choose number of matrices to be diagonalized: "
Nmatr = 1000 # strtoi(readline(prompt))

# x is an empty vector that will be used to collect all eigenvalues
x = c()

# The following if conditions select the proper ensemble and call the
# corresponding function Nmatr times in a loop. Matrices are diagonalized
# and the eigenvalues are collected in the vector x

rnorm_matrix <- function(N) {
    matrix(rnorm(N*N), nrow=N, ncol=N, byrow=TRUE)
}

if (beta == 1) { # Gaussian Orthogonal Ensemble
    bnd = sqrt(2*N)
    for (nm in 1:Nmatr) {
        M = rnorm_matrix(N)
        M = (M + t(M))/2
        y = eigen(M, symmetric=TRUE, only.values=TRUE)$values
        x = c(x, y)
    }
} else if (beta == 2) { # Gaussian Unitary Ensemble
    bnd = sqrt(4*N)
    for (nm in 1:Nmatr) {
        M = rnorm_matrix(N) + 1i*rnorm_matrix(N)
        M = (M + t(M))/2
        y = eigen(M, symmetric=TRUE, only.values=TRUE)$values
        x = c(x, y)
    }
} else if (beta == 4) { # Gaussian Symplectic Ensemble
    bnd = sqrt(8*N)
    for (nm in 1:Nmatr) {
        A = rnorm_matrix(N) + 1i*rnorm_matrix(N)
        B = rnorm_matrix(N) + 1i*rnorm_matrix(N)
        M = rbind(cbind(A, B), cbind(-Conj(B), Conj(A)))
        M = (M + t(M))/2
        y = unique(eigen(M, symmetric=TRUE, only.values=TRUE)$values)
        # The unique function gets rid of the double eigenvalues
        x = c(x, y)
    }
} else { # An error message is printed to screen if beta is different from 1, 2 or 4
    stop("Error: beta has to be equal to 1, 2, or 4")
}

# ------------------------------------------   
par(mfrow=c(1,1));
# 
hist(x, 50, freq=FALSE, 
     main="Normalized eigenvalue histogram")
abline(v=c(1,-1)*bnd, col="red", lty=2)

x2 <- seq(-1,1,length.out=100)*bnd
y2 <- sqrt(1-(x2/bnd)^2)/(2*pi)/sqrt(beta)
lines(x2, y2, col="red")


