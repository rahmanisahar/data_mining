cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}


library(ggplot2)
library(corrplot)
pdf(outfile)
m31 <- read.csv(infile,heade=TRUE, sep=',')
#small_m31=m31[c(2,3,4,5,7,8),]
small_m31 = m31
small_m31_1=small_m31[c(4:ncol(m31))]
x_cor=cor(small_m31_1, use='pairwise',method="pearson")
res1 <- cor.mtest(x_cor, 0.95)
res2 <- cor.mtest(x_cor, 0.99)
corrplot(x_cor, method = "circle", addrect = 3,tl.cex = 0.8,tl.col='black', insig = "blank",p.mat = res2[[1]],sig.level = 0.05)
dev.off()