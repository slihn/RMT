# This program plots a numerical histogram of the eigenvalue spacing
# distribution for 2x2 GOE random matrices, and compares it with the
# Wigner surmise. You will be asked to choose the number of matrices to 
# be diagonalized.

# Definition of the Wigner surmise function
surmise <- function(s) s*exp(-s^2/4)/2

# Reads the number of matrices to be diagonalized from the Command Window
prompt = "\n Choose number of matrices to be diagonalized: "
Nmatr = 50000 # strtoi(readline(prompt))

# x is an empty vector that will be used to collect all spacings
x = c()

for (nm in 1:Nmatr) {
    x1 = rnorm(1) 
    x2 = rnorm(1)
    x3 = rnorm(1)/sqrt(2)
    M = matrix(c(x1, x3, x3, x2), nrow=2, ncol=2, byrow=TRUE)
    y = eigen(M, symmetric=TRUE, only.values=TRUE)$values
    x = c(x, abs(diff(y)))
}

# ------------------------------------------   
par(mfrow=c(1,1));
# 
hist(x, 50,freq=FALSE, 
     main="Normalized spacing histogram")

x2 <- seq(0, max(x),length.out=100)
y2 <- surmise(x2)
lines(x2, y2, col="red")

