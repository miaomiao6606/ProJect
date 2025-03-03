library(openxlsx)
options(stringsAsFactors = F)
options(java.parameters = '-Xmx8000m')
library(nortest)
library(ggplot2)
library(vegan)
library(readr)
#library(rJava)
library("viper")
library(magrittr)
library(pheatmap)
library(mvtnorm)
library(swamp)#ComBat
#library('gPCA')#gpca
library('bapred')
library(Harman)#harman
library(sva)#SVA
library(grid)
library(gridExtra)
library(ggrepel)
library(ggfortify)
library(magrittr)
library(dplyr)
library(missForest)
library(missMDA)
library(pcaMethods)
library(PEMM)
library(viper)
library(stats)
library(reshape2)
library(bladderbatch)
library('imputeLCMD')
library("hydroGOF")
library("pls")
library(sn)
library("Amelia")
library(mice)

svd = function (x, nu = min(n, p), nv = min(n, p), LINPACK = TRUE) 
{
  # print("LINPACK:"); print(LINPACK) 
  if (!missing(LINPACK)) 
    warning("the LINPACK argument has been defunct since R 3.1.0")
  x <- as.matrix(x)
  if (any(!is.finite(x))) 
    stop("infinite or missing values in 'x'")
  dx <- dim(x)
  n <- dx[1L]
  p <- dx[2L]
  if (!n || !p) 
    stop("a dimension is zero")
  La.res <- La.svd(x, nu, nv)
  res <- list(d = La.res$d)
  if (nu) 
    res$u <- La.res$u
  if (nu) 
    res$u <- La.res$u
  if (nv) {
    if (is.complex(x)) 
      res$v <- Conj(t(La.res$vt))
    else res$v <- t(La.res$vt)
  }
  res
}
assignInNamespace("svd", svd, "base")


predi <- function(object, newdata, ncomp = 1:object$ncomp, comps,
                  type = c("response", "scores"),
                  na.action = na.pass, ...)
{
  if (missing(newdata) || is.null(newdata))
    newX <- model.matrix(object)
  else if (is.matrix(newdata)) {
    ## For matrices, simply check dimension:
    if (ncol(newdata) != length(object$Xmeans))
      stop("'newdata' does not have the correct number of columns")
    newX <- newdata
  } else {
    Terms <- delete.response(terms(object))
    m <- model.frame(Terms, newdata, na.action = na.action)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, m)
    mm <- model.matrix(Terms, m)
    saveattr <- attributes(mm)
    #print(mm)
    intercept <- which(saveattr$assign == 0)
    if (length(intercept))
    {
      mm <- mm[,-intercept, drop=FALSE]
      saveattr$dim <- dim(mm)
      saveattr$dimnames <- dimnames(mm)
      saveattr$assign <- saveattr$assign[-intercept]
      attributes(mm) <- saveattr
    }
    newX <- mm
  }
  
  nobs <- dim(newX)[1]
  
  ## Perform any scaling:
  if (!is.null(object$scale)) newX <- newX / rep(object$scale, each = nobs)
  type <- match.arg(type)
  if (type == "response") {
    if (missing(comps) || is.null(comps)) {
      ## Predict with models containing ncomp[1] components,
      ## ncomp[2] components, etc.
      if (missing(newdata)) return(fitted(object)[,,ncomp, drop=FALSE])
      B <- coef(object, ncomp = ncomp, intercept = TRUE)
      dPred <- dim(B)
      dPred[1] <- dim(newX)[1]
      dnPred <- dimnames(B)
      dnPred[1] <-
        if(is.null(dimnames(newX))) list(NULL) else dimnames(newX)[1]
      pred <- array(dim = dPred, dimnames = dnPred)
      for (i in seq(along = ncomp))
        pred[,,i] <- newX %*% B[-1,,i] + rep(B[1,,i], each = nobs)
      return(pred)
    } else {
      ## Predict with a model containing the components `comps'
      B <- rowSums(coef(object, comps = comps), dims = 2)
      B0 <- object$Ymeans - object$Xmeans %*% B
      pred <- newX %*% B + rep(B0, each = nobs)
      if (missing(newdata) && !is.null(object$na.action))
        pred <- napredict(object$na.action, pred)
      return(pred)
    }
  } else {
    if (missing(comps) || is.null(comps)) comps <- ncomp
    if (missing(newdata)) {
      TT <- object$scores[,comps]
      if (!is.null(object$na.action))
        TT <- napredict(object$na.action, TT)
    } else {
      if (is.null(object$projection))
        stop("`object' has no `projection' component.'")
      TT <- (newX - rep(object$Xmeans, each = nobs)) %*%
        object$projection[,comps]
    }
    return(TT)
  }
}
pfit <- function (formula, ncomp, data,  scale = TRUE, center = TRUE,...) 
{
  mf <- model.frame(formula,data)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  if (is.matrix(Y)) {
    if (is.null(colnames(Y))) 
      colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
  }
  else {
    Y <- as.matrix(Y)
    colnames(Y) <- deparse(formula[[2]])
  }
  mm <- model.matrix(mt, mf)
  saveattr <- attributes(mm)
  #print(mm)
  intercept <- which(saveattr$assign == 0)
  if (length(intercept))
  {
    mm <- mm[,-intercept, drop=FALSE]
    saveattr$dim <- dim(mm)
    saveattr$dimnames <- dimnames(mm)
    saveattr$assign <- saveattr$assign[-intercept]
    attributes(mm) <- saveattr
  }
  X <- mm
  nobj <- dim(X)[1]
  npred <- dim(X)[2]
  if (length(attr(mt, "term.labels")) == 1 && !is.null(colnames(mf[[attr(mt, 
                                                                         "term.labels")]]))) 
  { colnames(X) <- sub(attr(mt, "term.labels"), "", colnames(X))
  print("yes")}
  ncomp <- min(nobj - 1, npred)
  sdscale <- isTRUE(scale)
  val <- NULL
  if (sdscale) {
    scale <- sqrt(colSums((X - rep(colMeans(X), each = nobj))^2)/(nobj - 
                                                                    1))
    X <- X/rep(scale, each = nobj)
  }
  
  Y <- as.matrix(Y)
  dnX <- dimnames(X)
  dnY <- dimnames(Y)
  dimnames(X) <- dimnames(Y) <- NULL
  nobj  <- dim(X)[1]
  npred <- dim(X)[2]
  nresp <- dim(Y)[2]
  Xmeans <- colMeans(X)
  X <- X - rep(Xmeans, each = nobj)
  Ymeans <- colMeans(Y)
  Y <- Y - rep(Ymeans, each = nobj)
  R <- P <- matrix(0, ncol = ncomp, nrow = npred)
  tQ <- matrix(0, ncol = nresp, nrow = ncomp)
  B <- array(0, c(npred, nresp, ncomp))
  W <- P                       
  U <- TT <- matrix(0, ncol = ncomp, nrow = nobj)
  tsqs <- rep.int(1, ncomp)       
  fitted <- array(0, c(nobj, nresp, ncomp))
  XtY <- crossprod(X, Y)
  for (a in 1:ncomp) {
    if (nresp == 1) {
      w.a <- XtY / sqrt(c(crossprod(XtY)))
    } else {
      if (nresp < npred) {
        q <- eigen(crossprod(XtY), symmetric = TRUE)$vectors[,1]
        w.a <- XtY %*% q
        w.a <- w.a / sqrt(c(crossprod(w.a)))
      } else {
        w.a <- eigen(XtY %*% t(XtY), symmetric = TRUE)$vectors[,1]
      }
    }
    r.a <- w.a
    if (a > 5) {
      r.a <- r.a - colSums(crossprod(w.a, P[,1:(a-1), drop=FALSE]) %*% t(R[,1:(a-1), drop=FALSE]))
    } else if (a > 1) {
      for (j in 1:(a - 1))
        r.a <- r.a - c(P[,j] %*% w.a) * R[,j]
    }
    t.a <- X %*% r.a
    tsq <- c(crossprod(t.a))
    p.a <- crossprod(X, t.a) / tsq
    q.a <- crossprod(XtY, r.a) / tsq
    XtY <- XtY - (tsq * p.a) %*% t(q.a)
    R[,a]  <- r.a
    P[,a]  <- p.a
    tQ[a,] <- q.a
    B[,,a] <- R[,1:a, drop=FALSE] %*% tQ[1:a,, drop=FALSE]
    tsqs[a] <- tsq
    u.a <- Y %*% q.a / c(crossprod(q.a))
    if (a > 1) u.a <- u.a - TT %*% (crossprod(TT, u.a) / tsqs)
    U[,a]  <- u.a
    TT[,a] <- t.a
    W[,a]  <- w.a
    fitted[,,a] <- TT[,1:a] %*% tQ[1:a,, drop=FALSE]
  }
  residuals <- - fitted + c(Y)
  fitted <- fitted + rep(Ymeans, each = nobj) # Add mean
  objnames <- dnX[[1]]
  if (is.null(objnames)) objnames <- dnY[[1]]
  prednames <- dnX[[2]]
  respnames <- dnY[[2]]
  compnames <- paste("Comp", 1:ncomp)
  nCompnames <- paste(1:ncomp, "comps")
  dimnames(TT) <- dimnames(U) <- list(objnames, compnames)
  dimnames(R) <- dimnames(W) <- dimnames(P) <-
    list(prednames, compnames)
  dimnames(tQ) <- list(compnames, respnames)
  dimnames(B) <- list(prednames, respnames, nCompnames)
  dimnames(fitted) <- dimnames(residuals) <-
    list(objnames, respnames, nCompnames)
  class(TT) <- class(U) <- "scores"
  class(P) <- class(W) <- class(tQ) <- "loadings"
  z <- list(coefficients = B,
            scores = TT, loadings = P,
            loading.weights = W,
            Yscores = U, Yloadings = t(tQ),
            projection = R,
            Xmeans = Xmeans, Ymeans = Ymeans,
            fitted.values = fitted, residuals = residuals,
            Xvar = colSums(P * P) * tsqs,
            Xtotvar = sum(X * X))
  z$ncomp <- ncomp
  z$scale <- scale
  z$center <- TRUE
  z$validation <- val
  z$terms <- mt
  z$model <- mf
  z
}
initial <-function(data2,subdata2,flag=1)# use -1 when MNAR
{
  d=data2[!is.na(data2)]
  cp.est <- sn.mple(y=d,opt.method = "nlminb")$cp
  data=(data2-cp.est[1])/cp.est[2]
  meand=cp.est[1]
  sdd=cp.est[2]
  subdata=c()
  submean=c()
  subsd=c()
  trash=0
  for(i in 1:ncol(subdata2))
  {
    tmpre=subdata2[,i]
    d=tmpre[!is.na(tmpre)]
    cp.est <- sn.mple(y=d,opt.method = "nlminb")$cp
    #print(cp.est)
    submean=rbind(submean, cp.est[1])
    subsd=rbind(subsd,cp.est[2])
    kk=(tmpre-cp.est[1])/cp.est[2]
    subdata=cbind(subdata,kk)
  }
  colnames(subdata)=colnames(subdata2)
  rownames(submean)=colnames(subdata2)
  rownames(subsd)=colnames(subdata2)
  #subdata=scale(subdata2)
  lis=which(data %in% NA)
  impdata=data
  
  t=cbind(data,subdata)
  tuse=c()
  for(i in 1:nrow(t))
  {
    if(!is.na(sum(t[i,]))) tuse=rbind(tuse,t[i,])
  }
  
  if(length(tuse)==0 || nrow(tuse)<=ncol(tuse)*10)
  {
    for(k in 1:nrow(t))
    {
      #print(k)
      tmp=t[k,]
      if(length(which(is.na(tmp)))<length(tmp)*0.2)
      {
        tti=which(is.na(tmp))
        tmpmean=mean(tmp,na.rm=T)
        for(l in tti)
        {tmp[l]=tmpmean+runif(1,-1e-4,1e-4)}
        b=cbind(t(tmp))
        tuse=rbind(tuse,b)
        #}
      }
    }
  }
  
  model1res=c()
  model2res=c()
  idr=c()
  cnt=0
  if(length(tuse)==0 || nrow(tuse)<=10 || nrow(tuse)<ncol(tuse)*5)
  {
    reccor=c()
    for(j in 2:ncol(t))
    {
      tmpcor=c()
      for(i in 1:nrow(t))
      {
        if(!is.na(sum(t[i,c(1,j)]))) tmpcor=rbind(tmpcor,t[i,c(1,j)])
      }
      tmp=cor(tmpcor[,1],tmpcor[,2],method="pearson")
      reccor=rbind(reccor,cbind(j,tmp))
    }
    filid=reccor[order(reccor[,2]),1]
    trash=1
    tmpuse=tuse
    while(length(filid)>4 && trash<=length(filid)/2 && (length(tmpuse)==0 || nrow(tmpuse)<ncol(tmpuse)*5))
    {
      fil=filid[-c(1:trash)]
      fil=sort(fil)
      tt=t[,c(1,fil)]
      tmpuse=c()
      for(k in 1:nrow(tt))
      {
        #print(k)
        tmp=tt[k,]
        if(length(which(is.na(tmp)))<length(tmp)*0.2)
        {
          tti=which(is.na(tmp))
          tmpmean=mean(tmp,na.rm=T)
          for(l in tti)
          {tmp[l]=tmpmean+runif(1,-1e-4,1e-4)}
          b=cbind(t(tmp))
          tmpuse=rbind(tmpuse,b)
          #}
        }
      }
      trash=trash+1
    }
    tuse=tmpuse
  }
  ###
  if(nrow(tuse)>=10 && nrow(tuse)>ncol(tuse)*5)
  {
    for(i in 1:nrow(tuse))# 1:nrow(tuse)
    {
      
      corr=c()
      a=1:nrow(tuse)
      a=a[-i]
      df=tuse[i,2:ncol(tuse)]
      for(j in a)
      {
        corr=rbind(corr,cbind(j,cor(df,tuse[j,2:ncol(tuse)],method="pearson")))
      }
      subl=which(abs(corr[,2])>0.9)
      subl=corr[subl,1]
      if(length(subl)<ncol(tuse)*10 && nrow(tuse)>=ncol(tuse)*10)
      {
        corrtmp=corr[order(abs(corr[,2]),decreasing = T),]
        subl=corrtmp[1:(ncol(tuse)*10),1]
      }
      if(nrow(tuse)<ncol(tuse)*10 && nrow(tuse)>=ncol(tuse)*5)
      {
        subl=c(1:nrow(tuse))
      }
      subl=sort(subl)
      all=tuse[subl,]
      all=data.frame(all) 
      colnames(all)[1]="y2"
      model <- pfit(y2~., data=all)
      class(model)="mvr"
      distance=cor(all,all,method="pearson")
      distance=distance[1,2:ncol(distance)]
      
      tmpall=as.numeric(model$Xvar)/model$Xtotvar
      tmpall2=tmpall
      tmpid=model$ncomp
      if(tmpid>2)
      {
        for(ttt in 2:model$ncomp)
        {
          tmpall2[ttt]=tmpall[ttt]+tmpall2[ttt-1]
          if(is.na(tmpall2[ttt]))
          {break}
          if(tmpall2[ttt]>=0.9) {
            tmpid=ttt
            break}
          
        }
      }
      
      if(tmpid<=2) tmpid=2
      
      #jsid=5
      jsid=tmpid
      #jsid=10
      tmpsort=sort(abs(distance),decreasing=T)[jsid]
      disname=colnames(all)[2:ncol(all)]
      trashdisname=disname[abs(distance)<tmpsort]
      distance[abs(distance)<tmpsort]=0
      distance=abs(distance)/sum(abs(distance))
      #  pcr_pred=(sum(df*distance)+pcr_pred2)/2
      notall=all[,!colnames(all)%in%trashdisname]
      model2=pfit(y2~., data=notall)
      class(model2)="mvr"
      df=t(df)
      notdf=df[,!disname%in%trashdisname]
      notdf=data.frame(notdf)
      notdf=t(notdf)
      pcr_pred1=sum(df*distance)
      tmpall=as.numeric(model2$Xvar)/model2$Xtotvar
      tmpid2=model2$ncomp
      tmpall2=tmpall
      if(tmpid2>2)
      {
        for(ttt in 2:model2$ncomp)
        {
          tmpall2[ttt]=tmpall[ttt]+tmpall2[ttt-1]
          if(is.na(tmpall2[ttt]))
          {
            tmpid2=ttt-1
            break}
          if(tmpall2[ttt]>=1) {
            tmpid2=ttt
            break}
        }
      }
      if(tmpid2<model2$ncomp) cnt=cnt+1
      pcr_pred2=predi(model2, notdf,ncomp=tmpid2)
      if(is.na(pcr_pred2) ||abs(pcr_pred2)>abs(pcr_pred1*3) || abs(pcr_pred2*3)<abs(pcr_pred1)) {pcr_pred2=predi(model,df,ncomp=tmpid)
      cnt=cnt+1}
      if(is.na(pcr_pred2)) {pcr_pred2=pcr_pred1}
      if(abs(pcr_pred2)>abs(pcr_pred1*3) || abs(pcr_pred2*3)<abs(pcr_pred1)) {
        pcr_pred2=pcr_pred1
      }
      model1res=append(model1res,pcr_pred1)
      model2res=append(model2res,pcr_pred2)
      idr=append(idr,i)
    }
    real=tuse[idr,1]
    low =0
    high =1
    eps = 1e6
    reclow=rechigh=0.5
    if(length(idr)>10)
    {
      a1=quantile(real,0.9)
      ids=which(real<=a1)
      model1res=model1res[ids]
      model2res=model2res[ids]
      real=real[ids]
      while (low<=1)
      {
        combined=model1res*low+model2res*high
        x=sum((combined-real)^2)/length(real)
        if (x <eps)
        { eps=x  
        reclow=low
        rechigh=high
        }
        low=low+0.01
        high=high-0.01
      }
      reclow=reclow+cnt/length(real)+trash/(ncol(t)-1)
      if(reclow>=1) reclow=1
      rechigh=1-reclow
    }
  }else{
    reclow=rechigh=0.5
  }
  
  #####
  #reclow=0.67
  #rechigh=0.33
  if(trash>0)
  {fil2=fil-1
  subdata=subdata[,fil2]
  }
  for(i in 1:length(lis))
  {
    id=lis[i]
    dist=c()
    recordid=c()
    for(j in 1:ncol(subdata))
    {
      # tmp=discal(id,subdata[,j])
      tmp=subdata[,j]
      if(length(tmp[is.na(tmp)])<length(tmp) && !is.na(tmp[id]))
      {
        dist=cbind(dist,tmp)
        recordid=append(recordid,j)
      }
    }
    if(length(dist)==0)
    {  pcr_pred=mean(as.numeric(subdata2[id,]),na.rm=T)
    pcr_pred=(pcr_pred-meand)/sdd
    impdata[id]=pcr_pred
    }
    if(length(dist)>0)
    {
      selectrows=c()
      for(k in 1:nrow(dist))
      {
        #print(k)
        tmp=dist[k,]
        if(!is.na(sum(tmp)) )
        {
          # a=allright(tmp)
          # if(a==1){ b=cbind(k,t(tmp))
          b=cbind(k,t(tmp))
          selectrows=rbind(selectrows,b)
          #}
        }
      }
      if(nrow(selectrows)<=ncol(subdata)*10)
      {
        for(k in 1:nrow(dist))
        {
          #print(k)
          tmp=dist[k,]
          if(length(which(is.na(tmp)))<length(tmp)*0.2)
          {
            # a=allright(tmp)
            # if(a==1){ b=cbind(k,t(tmp))
            tti=which(is.na(tmp))
            tmpmean=mean(tmp,na.rm=T)
            for(l in tti)
            {tmp[l]=tmpmean+runif(1,-1e-4,1e-4)}
            b=cbind(k,t(tmp))
            selectrows=rbind(selectrows,b)
            #}
          }
        }
      }
      
      y=data[selectrows[,1]]
      js=which(y %in% NA)
      x=data.frame(selectrows[,2:ncol(selectrows)])
      colnames(x)=colnames(subdata)[recordid]
      if(length(js)>0){
        x2=x[-js,]
        y2=y[-js]
      }else{
        x2=x
        y2=y
      }
      all=cbind(y2,x2)
      all=data.frame(all)
      df=data.frame(t(subdata[id,recordid]))
      
      corr=c()
      a=1:nrow(all)
      for(j in a)
      {
        tmpdf=as.numeric(df)
        tmp=as.numeric(all[j,2:ncol(all)])
        corr=rbind(corr,cbind(j,cor(tmpdf,tmp,method="pearson")))
      }
      subl=which(abs(corr[,2])>0.9)
      subl=corr[subl,1]
      if(length(subl)>0 && length(subl)<ncol(all)*10)
      {
        corrtmp=corr[order(abs(corr[,2]),decreasing = T),]
        if(nrow(corrtmp)>ncol(all)*10)
        {subl=corrtmp[1:(ncol(all)*10),1]
        }else{
          subl=corrtmp[,1]
        }
      }else
      {
        subl=corr[,1]
      }
      subl=sort(subl)
      all=all[subl,]
      all=data.frame(all)
      colnames(all)[1]="y2"
      
      
      #df=data.frame(t(subdata[id,(jjs-1)[2:6]]))
      if(ncol(df)==1) colnames(df)="x2"
      if(ncol(all)<=3 || nrow(all)<=ncol(all)*5){
        pcr_pred=mean(as.numeric(subdata2[id,]),na.rm=T)
        pcr_pred=(pcr_pred-meand)/sdd
      }else{
        
        distance=cor(all,all,method="pearson")
        distance=distance[1,2:ncol(distance)]
        tmpmax=max(df)
        tmpmin=min(df)
        if(ncol(all)>3 && tmpmin<=max(all)*2 && tmpmax>=min(all)*0.5 && nrow(all)>ncol(all)*5){
          #jsid=ncol(all)/2
          #model <- plsr(y2~., data=all,scale=T, validation="none")
          model <- pfit(y2~., data=all)
          class(model)="mvr"
          tmpall=as.numeric(model$Xvar)/model$Xtotvar
          tmpall2=tmpall
          tmpid=model$ncomp
          if(tmpid>2)
          {
            for(ttt in 2:model$ncomp)
            {
              if(is.na(tmpall2[ttt]))
              {break}
              tmpall2[ttt]=tmpall[ttt]+tmpall2[ttt-1]
              if(tmpall2[ttt]>=0.9) {
                tmpid=ttt
                break}
            }
          }
          
          if(tmpid<=2) tmpid=2
          pcr_pred2=10000
          
          #jsid=5
          jsid=tmpid
          #jsid=10
          tmpsort=sort(abs(distance),decreasing=T)[jsid]
          disname=colnames(all)[2:ncol(all)]
          trashdisname=disname[abs(distance)<tmpsort]
          distance[abs(distance)<tmpsort]=0
          distance=abs(distance)/sum(abs(distance))
          if(tmpid==2 && max(distance)>=0.8)
          {
            distance[distance<max(distance)]=0
            distance[distance>0]=1
            trashdisname=disname[abs(distance)<1]
          }
          #  pcr_pred=(sum(df*distance)+pcr_pred2)/2
          if(nrow(all)>=model$ncomp*ncol(all))
          {
            notall=all[,!colnames(all)%in%trashdisname]
            #model2=plsr(y2~., data=notall,scale=T, validation="none")
            model2=pfit(y2~., data=notall)
            class(model2)="mvr"
            notdf=df[,!disname%in%trashdisname]
            if(length(notdf)==1)
            {
              notdf=data.frame(notdf)
              tmpname=disname[!disname%in%trashdisname]
              colnames(notdf)=tmpname
            }
            tmpid2=model2$ncomp
            if(model2$ncomp>2)
            {
              tmpall=as.numeric(model2$Xvar)/model2$Xtotvar
              
              tmpall2=tmpall
              for(ttt in 2:model2$ncomp)
              {
                tmpall2[ttt]=tmpall[ttt]+tmpall2[ttt-1]
                if(is.na(tmpall2[ttt]))
                {
                  tmpid2=ttt-1
                  break}
                if(tmpall2[ttt]>=1) {
                  tmpid2=ttt
                  break}
              }
            }
            pcr_pred2 <- predi(model2, notdf,ncomp=tmpid2)
            tmp=sum(df*distance)
            if(is.na(pcr_pred2)) pcr_pred2=10000
            if(abs(pcr_pred2)>abs(tmp*2) || abs(pcr_pred2*2)<abs(tmp)) pcr_pred2=10000
          }
          if(pcr_pred2>=10000 && ncol(all)>=max(ncol(subdata2)/3,5))
          { pcr_pred=sum(df*distance)
          #model <- plsr(y2~., data=all,scale=T, validation="none")
          model <- pfit(y2~., data=all)
          class(model)="mvr"
          tmpall=as.numeric(model$Xvar)/model$Xtotvar
          tmpid=model$ncomp
          if(tmpid>2)
          {
            tmpall2=tmpall
            for(ttt in 2:model$ncomp)
            {
              tmpall2[ttt]=tmpall[ttt]+tmpall2[ttt-1]
              if(tmpall2[ttt]>=1) {
                tmpid=ttt
                break}
            }
          }
          pcr_pred2 <- predi(model, df,ncomp=tmpid)
          tmp=sum(df*distance)
          if(is.na(pcr_pred2) || pcr_pred2>1e3) {pcr_pred2=pcr_pred
          }else if(abs(pcr_pred2)>abs(tmp*3) || abs(pcr_pred2*3)<abs(tmp)) {pcr_pred2=pcr_pred
          }else {pcr_pred=(sum(df*distance)*reclow+pcr_pred2*rechigh)}
          }else if(pcr_pred2>=10000 && ncol(all)<max(ncol(subdata2)/3,5)){
            pcr_pred=sum(df*distance)
          }else{
            pcr_pred=(sum(df*distance)*reclow+pcr_pred2*rechigh)
          }
          
        }else{
          distance=abs(distance)/sum(abs(distance))
          pcr_pred=sum(df*distance)
        }
      }
      # relation=lm(y2~.,data=all)
      impdata[id]=pcr_pred
    }
  }
  impdata2=impdata
  if(flag==-1)
  {
    tmpmin=min(data,na.rm=T)
    for(i in 1:length(lis))
    {
      id=lis[i]
      por=length(impdata[impdata<=impdata[id]])/length(impdata)
      tpor=length(impdata[impdata<=tmpmin])/length(impdata)
      allmin=min(c(subdata,data),na.rm=T)
      if(por >tpor)
      {
        impdata2[id]=runif(1,allmin,tmpmin)
      }
      
    }
  }
  return(list(imp.data=impdata2,mislis=lis,meand=meand,sdd=sdd,restsubmeta1=subdata))
}


#Usage
# function: initial(data2,subdata2,flag)
# Explain:
#   data2: the col of data that u want to impute. We recommend to do log transformation before imputation.
#   subdata2: the rest of the data without the col you want to impute. We recommend to do log transformation before imputation.
#     *Note: if you want to impute the whole dataset, you may need to do a for loop, and impute col by col.
#   flag: -1 means pure MNAR imputation. Otherwise is set to mixed imputation.
# 
# Example:
# We are trying to impute the first col of a dataset named RC:
# 
# RCN1=log(RC[,1]) #we take the first col
# restRC1=log(RC[,-1]) #rest of the data
# impRC=initial(RCN1, restRC1,-1) #we start the imputation
# impdata=impRC$imp.data #this is the normalized output
# meand=impRC$meand # this is the mean value of the output
# sdd=impRC$sdd # this is the sd of the output
# imptmp=exp(impdata * sdd + meand) # this is the final first col after the imputation.



