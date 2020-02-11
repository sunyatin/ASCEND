# 14 may 2019

global_decay_curve = function(X, minD = 0.1, maxD = 20, substracted = TRUE, return_all = FALSE) {
  if (return_all) {
    xy <- aggregate(X[,c("cov.pop", "cov.bg", "cov.substracted", "n.pairs")], list(X$bin.left.bound), sum)
    q = X$chrom[1]
    xy <- data.frame(x = subset(X, chrom==q)$bin.left.bound, 
                     cov.pop = xy[,2]/xy[,5], 
                     cov.bg = xy[,3]/xy[,5], 
                     cov.substracted = xy[,4]/xy[,5])
    xy <- subset(xy, x>=minD & x<=maxD)
    return(xy)
  } else {
    if (substracted) {
      xy <- aggregate(X[,c("cov.substracted", "n.pairs")], list(X$bin.left.bound), sum)
    } else {
      xy <- aggregate(X[,c("cov.pop", "n.pairs")], list(X$bin.left.bound), sum)
    }
    q = X$chrom[1]
    xy <- data.frame(x=subset(X, chrom==q)$bin.left.bound, y=xy[,2]/xy[,3])
    xy <- subset(xy, x>=minD & x<=maxD)
    return(xy)
  }
}


fit_decay_curve = function(X, minD = 0.1, maxD = 20) {
  X = subset(X, x>=minD & x<=maxD)
  
  options(show.error.messages = FALSE)
  coef = c('A'=NA, 't'=NA, 'c'=NA)
  fit = try( nls(y~A*exp(-2*x*t)+c, data=X, start=list(A=.1, t=.1, c=1e-3)), silent = TRUE )
  if (class(fit) != "try-error") coef = coef(fit)
  options(show.error.messages = TRUE)
  
  return(coef)
}


get_trues = function(x) {
  # get the true values in a file name
  library(stringr)
  y = unlist(strsplit(x, '_'))
  y.num <- as.numeric(str_extract(y, "[0-9]+"))
  return( data.frame("model"=y[1], "Tf.true"=y.num[2], 'Nf.true'=y.num[3], 'If.true'=y.num[4]/(2*y.num[3]), 'D'=y.num[4], 'rep'=y.num[5]) )
}

get_fit_beforeAug8 = function(X, stat, 
                   Nf.through.If = F,
                   method = 'jacquelin', 
                   check.chrom = TRUE, 
                   n.chrom = 22, 
                   do2D = FALSE,
                   NA.if.CI.incl.0 = TRUE) {
  
  meth = method
  x = subset(X, method == meth)
  
  if (check.chrom) if (nrow(x)-3 < n.chrom) return(data.frame(pop = x[1,'pop'], term = 1,
                                                          A = NA, t = NA, c = NA,
                                                          Tf = NA, Tf.LOW = NA, Tf.UP = NA,
                                                          Nf = NA, Nf.LOW = NA, Nf.UP = NA,
                                                          If = NA, If.LOW = NA, If.UP = NA))
  
  mean = subset(x, chromosome == 'MEAN')
  mjk = subset(x, chromosome == 'JK.MEAN')
  sde = subset(x, chromosome == 'JK.SD')
  
  if (stat == 'mean' | stat == 'var') {
    if (Nf.through.If) {
      y = c('A' = mjk$A,
            't' = mjk$t,
            'c' = mjk$c,
            'Tf' = 100*mjk$t,
            'Tf.LOW' = 100*mjk$t - 1.96 * 100*sde$t,
            'Tf.UP' = 100*mjk$t + 1.96 * 100*sde$t,
            'Nf' = mjk$X.1.Alog.abs.c..,
            'Nf.LOW' = mjk$X.1.Alog.abs.c.. - 1.96 * sde$X.1.Alog.abs.c..,
            'Nf.UP' = mjk$X.1.Alog.abs.c.. + 1.96 * sde$X.1.Alog.abs.c..,
            'If' = 1/2*mjk$X.Alog.abs.c..,
            'If.LOW' = 1/2*mjk$X.Alog.abs.c.. - 1.96 * 1/2*sde$X.Alog.abs.c..,
            'If.UP' = 1/2*mjk$X.Alog.abs.c.. + 1.96 * 1/2*sde$X.Alog.abs.c.. )
    } else {
      y = c('A' = mjk$A,
            't' = mjk$t,
            'c' = mjk$c,
            'Tf' = 100*mjk$t,
            'Tf.LOW' = 100*mjk$t - 1.96 * 100*sde$t,
            'Tf.UP' = 100*mjk$t + 1.96 * 100*sde$t,
            'Nf' = mjk$X1.A,
            'Nf.LOW' = mjk$X1.A - 1.96 * sde$X1.A,
            'Nf.UP' = mjk$X1.A + 1.96 * sde$X1.A,
            'If' = 1/2*mjk$X.Alog.abs.c..,
            'If.LOW' = 1/2*mjk$X.Alog.abs.c.. - 1.96 * 1/2*sde$X.Alog.abs.c..,
            'If.UP' = 1/2*mjk$X.Alog.abs.c.. + 1.96 * 1/2*sde$X.Alog.abs.c.. )
    }
    
    if (F) {
      y = c('A' = mean$A,
            't' = mean$t,
            'c' = mean$c,
            'Tf' = 100*mean$t,
            'Tf.LOW' = 100*mean$t - 1.96 * 100*sde$t,
            'Tf.UP' = 100*mean$t + 1.96 * 100*sde$t,
            'Nf' = mean$X1.A,
            'Nf.LOW' = mean$X1.A - 1.96 * sde$X1.A,
            'Nf.UP' = mean$X1.A + 1.96 * sde$X1.A )
    }
  }
  
  if (stat=='mean.jk') {
    if (Nf.through.If) {
      y = c('A' = mean$A,
            't' = mean$t,
            'c' = mean$c,
            'Tf' = 100*mean$t,
            'Tf.LOW' = 100*mean$t - 1.96 * 100*sde$t,
            'Tf.UP' = 100*mean$t + 1.96 * 100*sde$t,
            'Nf' = mean$X.1.Alog.abs.c..,
            'Nf.LOW' = mean$X.1.Alog.abs.c.. - 1.96 * sde$X.1.Alog.abs.c..,
            'Nf.UP' = mean$X.1.Alog.abs.c.. + 1.96 * sde$X.1.Alog.abs.c..,
            'If' = 1/2*mean$X.Alog.abs.c..,
            'If.LOW' = 1/2*mean$X.Alog.abs.c.. - 1.96 * 1/2*sde$X.Alog.abs.c..,
            'If.UP' = 1/2*mean$X.Alog.abs.c.. + 1.96 * 1/2*sde$X.Alog.abs.c.. )
    } else {
      y = c('A' = mean$A,
            't' = mean$t,
            'c' = mean$c,
            'Tf' = 100*mean$t,
            'Tf.LOW' = 100*mean$t - 1.96 * 100*sde$t,
            'Tf.UP' = 100*mean$t + 1.96 * 100*sde$t,
            'Nf' = mean$X1.A,
            'Nf.LOW' = mean$X1.A - 1.96 * sde$X1.A,
            'Nf.UP' = mean$X1.A + 1.96 * sde$X1.A,
            'If' = 1/2*mean$X.Alog.abs.c..,
            'If.LOW' = 1/2*mean$X.Alog.abs.c.. - 1.96 * 1/2*sde$X.Alog.abs.c..,
            'If.UP' = 1/2*mean$X.Alog.abs.c.. + 1.96 * 1/2*sde$X.Alog.abs.c.. )
    }
  }
  
  if (NA.if.CI.incl.0) {
    if (!is.na(y['Tf.LOW']) & y['Tf.LOW'] < 0) y = NA * y
    if (!is.na(y['Nf.LOW']) & y['Nf.LOW'] < 0) y = NA * y
  }
  
  y = data.frame(pop = x[1,'pop'], term = 1, t(y))
  
  if (do2D) {
    stop('error')
    Y = y
    
    if (stat == 'mean') {
      y = c('A' = mean$k2,
            't' = mean$t2,
            'c' = mean$c,
            'Tf' = 100*mean$t2,
            'Tf.LOW' = 100*mean$t2 - 1.96 * 100*sde$t2,
            'Tf.UP' = 100*mean$t2 + 1.96 * 100*sde$t2,
            'Nf' = 1/mean$k2,
            'Nf.LOW' = 1/mean$k2 - 1.96 * 1/sde$k2,
            'Nf.UP' = 1/mean$k2 + 1.96 * 1/sde$k2 )
    } else if (stat == 'var') {
      y = c('A' = mean$k2,
            't' = mean$t2,
            'c' = mean$c,
            'Tf' = 100*mean$t2,
            'Tf.LOW' = 100*mean$t2 - 1.96 * 100*sde$t2,
            'Tf.UP' = 100*mean$t2 + 1.96 * 100*sde$t2,
            'Nf' = 1/mean$k2,
            'Nf.LOW' = 1/mean$k2 - 1.96 * 1/sde$k2,
            'Nf.UP' = 1/mean$k2 + 1.96 * 1/sde$k2 )
    }
    
    if (!is.na(y['Tf.LOW']) & y['Tf.LOW'] < 0) y = NA * y
    # cannot filter on Nf because the ste is poorly designed
    #if (!is.na(y['Nf.LOW']) & y['Nf.LOW'] < 0) y = NA * y
    
    y = data.frame(pop = x[1,'pop'], term = 2, t(y))
    y = rbind(Y, y)
    return(y)
    if (all(is.na(y$Nf))) return(y[1,])
    if (any(is.na(y$Nf))) return(y[!is.na(y$Nf),])
    if (any(is.na(y$Nf))) return(y[!is.na(y$Nf),])
    y = y[which.min(y$Nf),]
  }
  
  return(y)
}


get_fit = function(X, check.chrom = TRUE, n.chrom = 22, set.NA.if.weird = TRUE) {

  x = X
  
  if (check.chrom) if (nrow(x)-3 < n.chrom) return(data.frame(pop = x[1,'pop'], term = 1,
                                                              A = NA, t = NA, c = NA,
                                                              Tf = NA, Tf.LOW = NA, Tf.UP = NA,
                                                              Nf = NA, Nf.LOW = NA, Nf.UP = NA,
                                                              If = NA, If.LOW = NA, If.UP = NA))
  
  mean = subset(x, chromosome == 'MEAN')
  mjk = subset(x, chromosome == 'JK.MEAN')
  sde = subset(x, chromosome == 'JK.SD')

  if (nrow(X)==0) {
     return(data.frame(pop = x[1,'pop'], term = 1,
              A = NA, t = NA, c = NA,
              Tf = NA, Tf.LOW = NA, Tf.UP = NA,
              Nf = NA, Nf.LOW = NA, Nf.UP = NA,
              If = NA, If.LOW = NA, If.UP = NA))
  }
  
  If.factor = exp(1)
  
  y = c('A' = mjk$A,
        't' = mjk$t,
        'c' = mjk$c,
        'Tf' = 100*mjk$t,
        'Tf.LOW' = 100*mjk$t - 1.96 * 100*sde$t,
        'Tf.UP' = 100*mjk$t + 1.96 * 100*sde$t,
        'Nf' = mjk$X1.A,
        'Nf.LOW' = mjk$X1.A - 1.96 * sde$X1.A,
        'Nf.UP' = mjk$X1.A + 1.96 * sde$X1.A,
        'If' = If.factor*mjk$A,
        'If.LOW' = If.factor*mjk$A - 1.96 * If.factor*sde$A,
        'If.UP' = If.factor*mjk$A + 1.96 * If.factor*sde$A )

  if (set.NA.if.weird) {
    if (!is.na(y['Tf']) & y['Tf'] < 0) y = NA * y
    if (!is.na(y['If']) & y['If'] < 0) y = NA * y
    if (!is.na(y['Tf.LOW']) & y['Tf.LOW'] < 0) y = NA * y
    if (!is.na(y['If.LOW']) & y['If.LOW'] < 0) y = NA * y
  }
  
  y = data.frame(pop = x[1,'pop'], term = 1, t(y))
  
  return(y)
}


merge_curve_fit = function(xy, fit, n_exp) {
  if (!is.data.frame(fit)) fit = as.data.frame(t(fit))
  if (n_exp==1) {
    xy = data.frame(xy, y.hat = as.numeric(NA))
    if (is.na(fit$A)) return(xy)
    xy$y.hat = sapply(xy$x, function(x) fit$A*exp(-2*x*fit$t)+fit$c)
    return(xy)
  }
  if (n_exp==2) {
    xy = data.frame(xy, y.hat = as.numeric(NA))
    if (is.na(fit$A)) return(xy)
    xy$y.hat = sapply(xy$x, function(x) fit$A*exp(-2*x*fit$t)+fit$A2*exp(-2*x*fit$t2)+fit$c)
    return(xy)
  }
}

from_fit_generate_per_chr_curves = function(fit, method = 'jacquelin', minD = .1, maxD = 20, stepD = .1) {
  require(reshape2)
  check.chrom = TRUE
  meth = method
  X = subset(fit, method == meth & chromosome!='MEAN' & chromosome!='JK.MEAN' & chromosome!='JK.SD')
  if (check.chrom) if (nrow(X) != 22) return(NA)
  x = seq(minD, maxD, by = .1)
  Y = apply(X[,4:6], 1, function(z) z[1]*exp(-2*x*z[2])+z[3])
  rownames(Y) = x
  Y = reshape2::melt(Y)
  names(Y) = c('x', 'chrom', 'y.hat')
  return(Y)
}



pairwise_Kelley = function(hapfile, n_samples, breakpoints, chr_length) {
  dat = fread(hapfile, select = c(3, 6:(5+n_samples)))
  dat = dat[order(dat[,1]),]
  pos = c(unlist(dat[,1]))
  pos = unname(pos)
  dat = as.matrix(dat[,-1])
  
  n = 1
  for ( i in  1:(n_samples-1) ) {
    for ( j in (i+1):n_samples ) {
      xp = pos[rowSums(dat[,c(i,j)])==1]
      xp = diff(xp)
      xp = hist(xp, breaks = breakpoints, plot = FALSE)
      if (n==1) XP = data.frame(x=xp$mids, y=xp$counts) else XP$y = XP$y + xp$counts
      n = n + 1
    }
  }
  XP$y = XP$y / (n-1)
  XP$y = XP$y / (2*chr_length)
  
  return(XP)
}

population_Kelley = function(hapfile, n_samples, breakpoints, chr_length) {
  dat = fread(hapfile, select = c(3, 6:(5+n_samples)))
  dat = dat[order(dat[,1]),]
  pos = c(unlist(dat[,1]))
  pos = unname(pos)
  dat = as.matrix(dat[,-1])
  
  xp = rowSums(dat)
  xp = (xp!=0) & (xp!=ncol(dat))
  xp = diff(pos[xp])
  xp = hist(xp, breaks = breakpoints, plot = FALSE)
  XP = data.frame(x=xp$mids, y=xp$counts)
  
  XP$y = XP$y / (n_samples*chr_length)
  
  return(XP)
}


fit_decay_curve_2D = function(X, minD = 0.1, maxD = 20) {
  # a + b*exp(p*x) + c*exp(q*x)
  
  M = subset(X, x>=minD & x<=maxD)
  
  x <- M[,1]
  y <- M[,2]
  
  n <- nrow(M)
  
  S <- 0
  for (k in 2:n) {
    S <- c(S, S[length(S)] + 1/2 * (y[k]+y[k-1])*(x[k]-x[k-1]))
  }
  
  SS <- 0
  for (k in 2:n) {
    SS <- c(SS, SS[length(SS)] + 1/2 * (S[k]+S[k-1])*(x[k]-x[k-1]))
  }
  
  M1 <- rbind(
    c(
      sum(SS**2),
      sum(SS*S),
      sum(SS*x**2),
      sum(SS*x),
      sum(SS)
    ),
    c(
      sum(SS*S),
      sum(S**2),
      sum(S*x**2),
      sum(S*x),
      sum(S)
    ),
    c(
      sum(SS*x**2),
      sum(S*x**2),
      sum(x**4),
      sum(x**3),
      sum(x**2)
    ),
    c(
      sum(SS*x),
      sum(S*x),
      sum(x**3),
      sum(x**2),
      sum(x)
    ),
    c(
      sum(SS),
      sum(S),
      sum(x**2),
      sum(x),
      n
    ))
  
  M2 <- c(
    sum(SS*y),
    sum(S*y),
    sum(x**2*y),
    sum(x*y),
    sum(y)
  )
  
  K <- solve(M1) %*% M2
  
  p <- 1/2*(K[2]+sqrt(K[2]**2+4*K[1]))
  q <- 1/2*(K[2]-sqrt(K[2]**2+4*K[1]))
  
  N1 <- rbind(
    c(
      n,
      sum(exp(p*x)),
      sum(exp(q*x))
    ),
    c(
      sum(exp(p*x)),
      sum(exp(2*p*x)),
      sum(exp((p+q)*x))
    ),
    c(
      sum(exp(q*x)),
      sum(exp((p+q)*x)),
      sum(exp(2*q*x))
    )
  )
  
  N2 <- c(sum(y), sum(y*exp(p*x)), sum(y*exp(q*x)))
  
  H <- solve(N1) %*% N2
  
  a <- H[1]
  b <- H[2]
  c <- H[3]
  
  #mod = nls(y ~ a + b*exp(p*x) + c*exp(q*x), start = list(a=a, b=b, c=c, p=p, q=q))
  #coef = coef(mod)
  #a=coef[1]; b=coef[2]; c=coef[3]; p=coef[4]; q=coef[5]
  
  return( c("A"=b, "t"=-p/2, "c"=a, "A2"=c, "t2"=-q/2) )
}


fit_decay_curve_1D <- function(X, minD = 0.1, maxD = 20) {
  
  M = subset(X, x>=minD & x<=maxD)
  
  x <- M[,1]
  y <- M[,2]

  n= length(x)
  
  y = y[order(x)]
  x = x[order(x)]
  
  dSk = diff(x) * (y[2:n] + y[1:(n-1)]) /2
  Sk = diffinv(dSk)
  
  dx = x - x[1]
  dy = y - y[1]
  
  m1 = matrix(NA,2,2)
  m1[1,1] = sum(dx**2)
  m1[1,2] = m1[2,1] = sum(dx*Sk)
  m1[2,2] = sum(Sk**2)
  
  m2 = c( sum(dx*dy), sum(dy*Sk))
  
  nc = solve(m1) %*% m2
  
  c = nc[2,1]
  
  m3 = matrix(NA,2,2)
  m3[1,1] = n
  m3[1,2] = m3[2,1] = sum(exp(c*x))
  m3[2,2] = sum(exp(2*c*x))
  
  m4 = c(sum(y), sum(y*exp(c*x)))
  
  ab = solve(m3) %*% m4
  
  return (c('A'=ab[2,1], 't'=-1/2*c, 'c'=ab[1,1]))  
}



# sept 27th 2019 -- to check
weighted_jackknife = function(chrom, values, chrom_weights_matrix) {
  # Busing F, Meijer E, Leeden R (1999) Delete-m Jackknife for Unequal m. Statistics and Computing 9: 3–8.
  
  if (length(chrom)!=length(values)) stop('chrom vector should be the same size as values vector')
  if (ncol(chrom_weights_matrix)!=2) stop('chrom_weights_matrix should be a matrix with two columns: chrom ID, chrom weight')
  
  chrom = as.character(chrom)
  general_mean = mean(values, na.rm=T)
  chrom_means = c(by(values, chrom, FUN=mean, na.rm=T))
  chrom_weights_matrix = data.frame(q = as.character(chrom_weights_matrix[,1]),
                                    w = as.numeric(as.character(chrom_weights_matrix[,2])))
  
  CHR = data.frame(q=names(chrom_means), mean=chrom_means)
  CHR = merge(CHR, chrom_weights_matrix, by='q', sort=F, all.x=TRUE, all.y=FALSE)
  
  mrem = sapply(seq_along(CHR$mean), function(i) mean(CHR$mean[-i], na.rm=T))
  CHR$mean = mrem
  
  n = sum(CHR$w)
  g = nrow(CHR)

  mean = g*general_mean - sum( (1-CHR$w/n)*CHR$mean  )
  
  H = n/CHR$w
  term1 = sum( (1-CHR$w/n)*CHR$mean )
  term2 = H*general_mean - (H-1)*CHR$mean - g*general_mean + term1
  var = 1/g * sum( 1/(H-1)*term2**2 )
  sd = sqrt(var)
  
  return(c('general_mean'=general_mean,
           'jackknife_mean'=mean,
           'jackknife_sde'=sd))
}




