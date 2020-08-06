# This program plots a histogram of the eigenvalue densities obtained via
# numerical diagonalization of Wishart-Laguerre random matrices (for 
# beta = 1,2,4) and compares them with the Marcenko-Pastur density for
# the corresponding matrix size. You will be asked to provide the sizes
# of the rectangular NxM matrices H used to build Wishart-Laguerre
# matrices of the type W = X'*X, and the number of matrices to be 
# diagonalized.

# Reads the matrix size from the Command Window, M >= N
prompt = "Choose first matrix size: "
N = 5 # strtoi(readline(prompt))
prompt = "Choose second matrix size: "
M = 20 # strtoi(readline(prompt))

if (N >= M) stop("ERROR: N must be larger than M") 

# Defining the Marcenko-Pastur density function
c = N/M
xmin = (1-1/sqrt(c))^2
xmax = (1+1/sqrt(c))^2
rho <- function(x) sqrt((x-xmin)*(xmax-x))/(2*pi*x)

# Reads the number of matrices to be diagonalized from the Command Window
prompt = "Choose number of matrices to be diagonalized: "
Nmatr = 10000 # strtoi(readline(prompt))

# These vectors will be used to collect all eigenvalues
x1 = c()
x2 = c()
x4 = c()

rnorm_matrix2 <- function(N,M) {
    matrix(rnorm(N*M), nrow=N, ncol=M, byrow=TRUE)
}

for (nm in 1:Nmatr) {
    # beta = 1
    beta = 1
    H = rnorm_matrix2(N,M)
    W = H %*% t(H)
    y1 = eigen(W, symmetric=TRUE, only.values=TRUE)$values
    x1 = c(x1, y1/(beta*N)) # Notice the rescaling of the eigenvalues
    
    # beta = 2
    beta = 2
    H = rnorm_matrix2(N,M) + 1i*rnorm_matrix2(N,M)
    W = H %*% t(Conj(H)) 
    y2 = eigen(W, symmetric=TRUE, only.values=TRUE)$values
    x2 = c(x2, y2/(beta*N)) # Notice the rescaling of the eigenvalues
    
    # beta = 4
    beta = 4;
    A = rnorm_matrix2(N,M) + 1i*rnorm_matrix2(N,M)
    B = rnorm_matrix2(N,M) + 1i*rnorm_matrix2(N,M)
    H = rbind(cbind(A, B), cbind(-Conj(B), Conj(A)))
    W = H %*% t(Conj(H))  
    y4 = unique(round(eigen(W, symmetric=TRUE, only.values=TRUE)$values,5))
    x4 = c(x4, y4/(beta*N)) # Notice the rescaling of the eigenvalues
    
}

# ------------------------------------------   
par(mfrow=c(2,2));

# Plotting the Marcenko-Pastur density
rx <- seq(xmin, xmax, length.out=100)
plot(rx, rho(rx), type="l", col="blue", main="Marcenko-Pastur density")

# Plotting the histograms for the three numerical eigenvalue densities
hist(x1, 30, freq=FALSE, col="blue", main="MP density, beta = 1")
lines(rx, rho(rx), col="red", lwd=2)

hist(x2, 30, freq=FALSE, col="cyan", main="MP density, beta = 2")
lines(rx, rho(rx), col="red", lwd=2)

hist(x4, 30, freq=FALSE, col="orange", main="MP density, beta = 4")
lines(rx, rho(rx), col="red", lwd=2)

