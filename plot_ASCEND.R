#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if ( length(args)==0 ) {
  stop("At least one argument must be supplied", call.=FALSE)
}

out = args[1]
fit = args[2]
output = args[3]
substracted = as.logical(args[4])
NRMSD_threshold = as.numeric(args[5])

print(out)
print(fit)
print(output)
print(substracted)
print(NRMSD_threshold)

#########

X = read.table(out, header=T)
fit = read.table(fit, header=T)

if (substracted) {
  xy = aggregate(X[,c("cor.substracted", "n.pairs")], list(X$bin.left.bound), sum)
} else {
  xy = aggregate(X[,c("cor.pop", "n.pairs")], list(X$bin.left.bound), sum)
}
q = X$chrom[1]
xy = data.frame(x=subset(X, chrom==q)$bin.left.bound, y=xy[,2]/xy[,3])

par(mfrow=c(1,1))

POP = as.character(fit$pop)[1]

x = data.frame(
  Tf = 100 * subset(fit, chromosome=='JK.MEAN')$t,
  Tf.LOW = 100 * (subset(fit, chromosome=='JK.MEAN')$t - 1.96 * subset(fit, chromosome=='JK.SE')$t),
  Tf.UP = 100 * (subset(fit, chromosome=='JK.MEAN')$t + 1.96 * subset(fit, chromosome=='JK.SE')$t),
  If = exp(1) * subset(fit, chromosome=='JK.MEAN')$A,
  If.LOW = (exp(1) * subset(fit, chromosome=='JK.MEAN')$A) - 1.96 * (exp(1) * subset(fit, chromosome=='JK.SE')$A),
  If.UP = (exp(1) * subset(fit, chromosome=='JK.MEAN')$A) + 1.96 * (exp(1) * subset(fit, chromosome=='JK.SE')$A),
  NRMSD = subset(fit, chromosome=='MEAN')$NRMSD)

if (nrow(x)==0) {
  print('Fitting failed')
  y = as.data.frame(t(rep(NA, 7)))
  names(y) = names(x)
  x = y
}

if (!is.na(x['Tf']) & x['Tf'] < 0) x = NA * x
if (!is.na(x['If']) & x['If'] < 0) x = NA * x
if (!is.na(x['Tf.LOW']) & x['Tf.LOW'] < 0) x = NA * x
if (!is.na(x['If.LOW']) & x['If.LOW'] < 0) x = NA * x

yfit = subset(fit, chromosome=='JK.MEAN')$A * exp(-2*xy[,1]*subset(fit, chromosome=='JK.MEAN')$t) + subset(fit, chromosome=='JK.MEAN')$c

png(filename = output)

  a = adjustcolor("gray", alpha.f = 0.2)

  plot(xy[,1:2], main=paste0(POP,'\nSubtracted: ',substracted), 
       cex.main = 1.3, xlab='Genetic distance (cM)', ylab='Allele sharing correlation', pch=19, col='#404040')
  if (is.na(x$Tf)|x$NRMSD>NRMSD_threshold) rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 'gray')
  if (is.na(x$Tf)|x$NRMSD>NRMSD_threshold) points(xy[,1:2], pch=19, col='#404040')
  
  if (!is.na(x$Tf) && x$NRMSD<NRMSD_threshold) {
    xy = data.frame(xy, yfit = yfit)
    lines(xy[,c(1,3)], col='red', lwd=2)
    leg = paste0('Tf = ',round(x$Tf, 0),' [',round(x$Tf.LOW),'; ',round(x$Tf.UP),']\n')
    leg = paste0(leg, 'If = ',round(100*x$If, 2),'% [',round(100*x$If.LOW, 2),'; ',round(100*x$If.UP, 2),']\n')
    leg = paste0(leg, 'NRMSD = ',round(x$NRMSD, 3),'\n')
    legend('topright', legend = leg, col = 'red', cex = 1.3)
  } else {
    if (is.na(x$Tf)) {
    } else {
      leg = paste0('\n')
      leg = paste0(x$excluded, '\n')
      leg = paste0(leg, 'NRMSD = ',round(x$NRMSD, 3),'\n')
      legend('topright', legend = leg, col = 'red', cex = 1.3)
    }
  }

dev.off()



