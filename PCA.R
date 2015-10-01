library(scatterplot3d)
library(rgl)
library(RColorBrewer)


runPCAexample = function (){
  #everything until the last line is for simulating data
  #set the number of experimental conditions, number of replicates, and number of genes
  nConditions = 3
  nReps = 3
  nGenes = 10000
  #spikes are signal, spikeCount number of random genes will be shifted by spikeSize
  spikeCount = 400
  spikeSize = 2
  #conditions and reps are arrays of strings describing each column of data table
  conditions = getFakeConditions(nConditions = nConditions, nReps = nReps)
  reps = getFakeReps(nConditions = nConditions, nReps = nReps)
  #randomized dataset with genes as rows and condition+replicates as columns
  #replace with a similar matrix object
  data = getFakeData( conditions, reps, nGenes = nGenes,spikeCount = spikeCount, spikeSize = spikeSize)
  
  #draws a scatterplot table of top n principle components
  #you must provide a condition for each column in data and a replicate
  #examples would be c('ctrl','ctrl','exp','exp') and c('r1','r2','r1','r2') for a simple experiment versus control with 2 replicates of each
  #you can provide a title or suppress the key if you want to
  #if you're familiar with basic plots in R you can fiddle with the point shape/weight/size/color in the plotPCA function
  geneList = plotPCA(data, conditions, reps,n = 6, title = "PCA Example")
}

runAffyPCAexample = function (){
  #for this example i have already removed the header lines from the summary file
  #for time purposes, i am only reading 20,000 lines of the summary file, remove the nrows argument to load the full file
  data=as.matrix(read.table(file = "H:/R_workspace/example_summary.txt",row.names = 1,stringsAsFactors = F))
  labels=as.character(data[1,1:ncol(data)])
  data=data[2:nrow(data),]
  rnames=rownames(data)
  data=matrix(as.numeric(data),ncol = length(labels),byrow = T)
  rownames(data)=rnames
  labels=strsplit(labels,"[ _]")
  labels=matrix(unlist(labels), nrow = length(labels),byrow = T)
  #IMPORTANT - i renamed your .CEL files for my analysis and you will have to change which columns of labels to use as appropriate for you sample names
  conditionLabels=labels[,3]
  replicateLabels=labels[,5]
  
  #draws a scatterplot table of top n principle components
  #you must provide a condition for each column in data and a replicate
  #examples would be c('ctrl','ctrl','exp','exp') and c('r1','r2','r1','r2') for a simple experiment versus control with 2 replicates of each
  #you can provide a title or suppress the key if you want to
  #if you're familiar with basic plots in R you can fiddle with the point shape/weight/size/color in the plotPCA function
  geneList = plotPCA(data, conditionLabels, replicateLabels,n = 3, title = "PCA Example")
}

plotPCA=function(data, conditionLabels, replicateLabels, secondaryConditionLabels = NA, n = 6, title = "", key = TRUE, displayReps = F, do3D = F, doVarPlot = F){
  ## performs PCA and constructs a table of scatterplots representing results
  #  args
  #  data, a matrix like object with genes are rows and condition/replicates as columns
  #  conditionLabels, length must equal ncol of data.  describes condition being varied across data columns ex: wt,wt,exp,exp
  #  repliateLabels, length must equal ncol of data.  ex: r1,r2,r1,r2 or a,b,a,b
  #  n, the number of principle components to draw.
  #  title, a descriptive title to draw at top.  will also add a subtitle describing dataset characteristics
  #  key, whether the key should be drawn
  #
  # returns list of gene names sorted by their influence on the PCA
  
  #how many conditional attributes there are
  labelComplexity = 2
  metaConditionLabels = conditionLabels
  if(is.na(secondaryConditionLabels[1]) || secondaryConditionLabels == conditionLabels){
    secondaryConditionLabels = conditionLabels
    labelComplexity = 1
  } else{
    metaConditionLabels = paste(conditionLabels, secondaryConditionLabels)
  }
  
  data=as.matrix(data)
  uniqueGroups=unique(conditionLabels)
  uniqueSecondary=unique(secondaryConditionLabels)
  uniqueMeta=unique(metaConditionLabels)
  nCond1 = length(uniqueGroups)
  nCond2 = length(uniqueSecondary)
  nConditions = length(uniqueMeta)
  nReps = length(unique(replicateLabels))
  #modify how color groups are set if you don't like the colors
  pointGroups=1:length(uniqueSecondary)
  #   clrRamp = colorRamp(colors = c('darkorange','red','blue')) #replaced color ramp with color brewer
  #   colorGroups=rgb(clrRamp((1:nCond1-1)/(nCond1-1))/255)
  colorGroups = RColorBrewer::brewer.pal(n = length(uniqueGroups), name = 'Dark2')
  names(colorGroups)=uniqueGroups
  names(pointGroups)=uniqueSecondary
  condPrin=prcomp(t(data))
  
  #print("Variance explained by PC:")
  cumVar = paste(format(cumsum((condPrin$sdev)^2) / sum(condPrin$sdev^2),digits = 3,scientific = F ), collapse = '\n')
  indiVar = paste(substr(format(((condPrin$sdev)^2) / sum(condPrin$sdev^2),digits = 3,scientific = F ),start = 1,stop = 5), collapse = '\n')
  if(doVarPlot){
    layout(1)
    plot(x = 0:1, y = 0:1, type="n", axes=F, xlab="", ylab="")
    text("Variance explained by PC:",x = .2,y = .9, adj = c(0,1))
    text(paste(1:length(condPrin$sdev), collapse = '\n'),x = .1,y = .7, adj = c(0,1))
    text(indiVar,x = .2,y = .7, adj = c(0,1))
    text("Cumulative",x = .35,y = .8, adj = c(0,1))
    text(cumVar,x = .35,y = .7,adj = c(0,1))
  }
  
  geneImpact=condPrin$rotation
  if(do3D){
    layout(1:2, heights = c(.1,1))
    #dev.off()
    par(mai = rep(0,4))
    plot(c(0,1), type = 'n', axes = F, xlab = '', ylab = '')
    legend(x = 'center', legend = uniqueMeta,, col = colorGroups[uniqueMeta],pch = pointGroups[uniqueMeta], horiz = T)
    scatterplot3d(condPrin$x[,1],condPrin$x[,2],condPrin$x[,3], 
                  pch = pointGroups[secondaryConditionLabels],
                  color = colorGroups[conditionLabels],
                  type='h',
                  xlab = 'PC1',ylab = 'PC2',zlab = 'PC3',
                  lwd=2,axis=T)
  }
  #plot3d(condPrin$x[,1],condPrin$x[,2],condPrin$x[,3],pch =  pointGroups[secondaryConditionLabels],col=colorGroups[conditionLabels],size = 10,xlab = 'PC1',ylab = 'PC2',zlab = 'PC3',axes=T)
  #   
  
  subTitle=paste(nrow(data), 'genes,', nConditions, 'conditions,', nReps, 'replicates')
  plotsArea=mat=matrix(data=1:n^2,nrow=n,ncol=n)
  if(title!="" || key){
    #add a top row for title
    a=matrix(0,nrow=1 + nrow(plotsArea), ncol=ncol(plotsArea))
    a[1,]=1
    a[2:nrow(a),]=plotsArea+1
    plotsArea=a
  }
  
  nf=layout(plotsArea)
  #layout.show(nf)
  par(mar=rep(0,4))
  if(key || title!=""){
    frame()
  }
  if(key){
    #draw the key
    legend(x = .01,xjust = 0, y = .9, yjust = .5, legend = uniqueMeta, col = colorGroups[sort(rep(1:nCond1,nConditions/nCond1))],pch = rep(pointGroups[1:nCond2],nConditions/nCond2),cex = .8,horiz = T)
  }
  
  if(title!=""){
    #draw the title and subtitle
    text(x = .5, y = .5,title,cex = 3)
    #text(x = .5, y = .2,subTitle,cex = 1)
  }
  
  for(c in 1:n){
    for(r in 1:n){
      if(c != r){
        #modify this line to change how the scatterplots are drawn
        if(displayReps){#draw points as colored text
          plot(condPrin$x[,c],condPrin$x[,r], type = 'n', xlab="",ylab="",xaxt='n',yaxt='n')  
          for(i in 1:length(replicateLabels)){
            text(condPrin$x[i,c], condPrin$x[i,r], replicateLabels[i], col=colorGroups[conditionLabels[i]] )
          }
          
        }else{#draw points as colored shapes
          plot(condPrin$x[,c],condPrin$x[,r],col=colorGroups[conditionLabels],pch=pointGroups[secondaryConditionLabels],lwd=2,xlab="",ylab="",xaxt='n',yaxt='n')  
        }
        
      }
      else{
        boxPlotData=list()#matrix(0,ncol = nConditions, nrow = nReps)
        for(i in 1:length(uniqueMeta)){
          k = uniqueMeta[i]
          subset=metaConditionLabels==k
          boxPlotData[[i]] = condPrin$x[subset,c]
        }
        sets = matrix(0,nrow = nrow(data), ncol = 1)
        boxplot(boxPlotData,yaxt='n',xaxt='n',col = colorGroups[uniqueGroups], lwd = .2, border = colorGroups[uniqueGroups])
        
      }
    }
  }
  
  a= matrix(0,nrow = nrow(data), ncol=2)
  a[,1] = rownames(data)
  a[,2] = apply(geneImpact, FUN = function(x) max(abs(x)), MARGIN = 1)
  a = a[order(a[,2],decreasing = TRUE)]
  return(list(a, geneImpact))
}

getFakeConditions = function (nConditions=3, nReps=2){
  possibleText = c("Alpha", "Beta", "Gamma", "Delta", "Epsilon","Omega","Apple","Banana","Candy")
  conditionLabels=sort(rep(possibleText[1:nConditions],nReps))
  return(conditionLabels)
}

getFakeReps = function(nConditions=3, nReps=2){
  repLabels= rep(c(1:nReps), nConditions)
  return(repLabels)
}

getFakeData= function( conditionLabels, repLabels,nGenes=10000, seed=1,spikeCount=100,spikeSize=10){
  #create some data
  set.seed(seed)
  nConditions = length(unique(conditionLabels))
  nReps = length(unique(repLabels))
  a=rnorm(n = nGenes*nConditions*nReps,mean = 50,sd = 2 )
  
  
  data=matrix(a,nrow = nGenes,ncol = nConditions * nReps)
  rownames(data)=paste('gene',1:nGenes)
  colnames(data) = conditionLabels
  uniqueGroups=unique(conditionLabels)
  alternate=TRUE
  for(i in 1:spikeCount){
    randCondition=uniqueGroups[sample.int(n = nConditions,size = 1)]
    subset=conditionLabels==randCondition
    #3 genes are dynamic in condition 2 and 3
    randGene=sample.int(n = nGenes,size = 1)
    s=spikeSize
    if(alternate){
      alternate = FALSE
      s=-spikeSize
    }
    else{
      alternate = TRUE
    }
    data[randGene,subset] = data[randGene,subset]-s
    
  }
  return(data)
}

scatter3D = function (x, y, z, ..., colvar = z, phi = 40, theta = 40, col = NULL, 
                      NAcol = "white", colkey = NULL, panel.first = NULL, clim = NULL, 
                      clab = NULL, bty = "b", CI = NULL, surf = NULL, add = FALSE, 
                      plot = TRUE) 
{
  plist <- initplist(add)
  dot <- splitdotpersp(list(...), bty, NULL, x, y, z, plist = plist)
  x <- as.vector(x)
  y <- as.vector(y)
  z <- as.vector(z)
  len <- length(x)
  if (length(y) != len) 
    stop("'y' should have same length as 'x'")
  if (length(z) != len) 
    stop("'z' should have same length as 'x'")
  if (len > 1 & ispresent(colvar)) {
    if (length(colvar) != len) 
      stop("'colvar' should have same length as 'x', 'y' and 'z'")
    colvar <- as.vector(colvar)
    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
    if (dot$clog) {
      colvar <- log(colvar)
      clim <- log(clim)
    }
    if (is.null(col)) 
      col <- jet.col(100)
    iscolkey <- is.colkey(colkey, col)
    if (iscolkey) 
      colkey <- check.colkey(colkey)
    if (!is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
    Col <- variablecol(colvar, col, NAcol, clim)
    if (length(Col) == 1) 
      Col <- rep(Col, length.out = len)
  }
  else {
    if (is.null(col)) 
      col <- "black"
    if (!is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
    Col <- rep(col, length.out = len)
    iscolkey <- FALSE
  }
  if (is.null(plist)) {
    do.call("perspbox", c(alist(x = range(x), y = range(y), 
                                z = range(z, na.rm = TRUE), phi = phi, theta = theta, 
                                plot = plot, colkey = colkey, col = col), dot$persp))
    plist <- getplist()
  }
  fit <- NULL
  if (!is.null(surf)) {
    if (!is.list(surf)) 
      stop("'surf' should be a 'list' or 'NULL'")
    fit <- surf$fit
    surf$fit <- NULL
  }
  if (!is.null(fit)) {
    if (!is.null(CI)) {
      if (!is.null(CI$z)) 
        stop("cannot combine a confidence interval (CI) on 'z' with 'fit' in 'surf'")
    }
    else CI <- list()
    if (length(fit) != length(z)) 
      stop("'fit', argument of 'surf'  should be of equal size of 'z'")
    disttoz <- fit - z
    CIz <- matrix(ncol = 2, data = c(-disttoz, disttoz))
    CIz[CIz > 0] <- 0
    CI$z = CIz
    CI$alen = 0
  }
  isCI <- is.list(CI)
  if (isCI) 
    CI <- check.CI(CI, len, 3)
  if (is.null(CI) & any(Col == "transparent")) {
    ii <- Col != "transparent"
    x <- x[ii]
    y <- y[ii]
    z <- z[ii]
    Col <- Col[ii]
    len <- length(x)
  }
  if (is.function(panel.first)) 
    panel.first(plist$mat)
  if (!is.null(surf)) {
    if (is.null(surf$colvar)) 
      surf$colvar <- surf$z
    if (is.null(surf[["col"]])) {
      surf$col <- col
      if (is.null(surf$clim)) 
        surf$clim <- clim
    }
    if (is.null(surf$clim)) 
      surf$clim <- range(surf$colvar)
    surf$colvar[surf$colvar < min(surf$clim)] <- NA
    surf$colvar[surf$colvar > max(surf$clim)] <- NA
    surf$z[surf$z < dot$persp$zlim[1]] <- NA
    surf$z[surf$z > dot$persp$zlim[2]] <- NA
    spoly <- do.call("addimg", c(alist(poly = NULL, plist = plist), 
                                 surf))
  }
  else spoly <- NULL
  sseg <- NULL
  dtype <- dot$points$type
  dot$points$type <- NULL
  lwd <- dot$points$lwd
  if (is.null(lwd)) 
    lwd <- 1
  lty <- dot$points$lty
  if (is.null(lty)) 
    lty <- 1
  alpha <- dot$alpha
  if (is.null(alpha)) 
    alpha <- NA
  alpha <- rep(alpha, length.out = len)
  if (is.null(dtype)) 
    dtype <- "p"
  if (dtype == "h") {
    zmin <- dot$persp$zlim[1]
    Proj <- project(x, y, 0.5 * (z + zmin), plist, FALSE)
    sseg <- list(x.from = x, x.to = x, y.from = y, y.to = y, 
                 z.from = rep(zmin, len), z.to = z, col = Col, lwd = rep(lwd, 
                                                                         length.out = len), lty = rep(lty, length.out = len), 
                 alpha = alpha, proj = Proj)
    class(sseg) <- "segments"
  }
  else if (dtype %in% c("b", "l", "o")) {
    LCol <- Col
    if (length(LCol) > 1) {
      LCol <- cbind(Col[-1], Col[-len])
      LCol <- apply(LCol, MARGIN = 1, FUN = MeanColors)
      if (!is.null(dot$alpha)) 
        LCol <- setalpha(LCol, dot$alpha)
    }
    Proj <- project(0.5 * (x[-len] + x[-1]), 0.5 * (y[-len] + 
                                                      y[-1]), 0.5 * (z[-len] + z[-1]), plist)
    sseg <- list(x.from = c(sseg$x.from, x[-len]), x.to = c(sseg$x.to, 
                                                            x[-1]), y.from = c(sseg$y.from, y[-len]), y.to = c(sseg$y.to, 
                                                                                                               y[-1]), z.from = c(sseg$z.from, z[-len]), z.to = c(sseg$z.to, 
                                                                                                                                                                  z[-1]), col = c(sseg$col, LCol), lwd = c(sseg$lwd, 
                                                                                                                                                                                                           rep(lwd, length.out = len - 1)), lty = c(sseg$lty, 
                                                                                                                                                                                                                                                    rep(lty, length.out = len - 1)), alpha = c(sseg$alpha, 
                                                                                                                                                                                                                                                                                               alpha), proj = c(sseg$proj, Proj))
    class(sseg) <- "segments"
  }
  pch <- dot$points$pch
  if (is.null(pch)) 
    pch <- 1
  bg <- dot$points$bg
  if (is.null(bg)) 
    bg <- 1
  cex <- dot$points$cex
  if (is.null(cex)) 
    cex <- 1
  if (dtype == "l") 
    dopoints <- FALSE
  else dopoints <- TRUE
  CIpt <- NULL
  pt <- NULL
  if (!is.null(CI)) {
    CIpt <- list(x.from = NULL, y.from = NULL, z.from = NULL, 
                 x.to = NULL, y.to = NULL, z.to = NULL, x.mid = x, 
                 y.mid = y, z.mid = z, col = Col, pch = rep(pch, length.out = len), 
                 cex = rep(cex, length.out = len), bg = rep(bg, length.out = len), 
                 alpha = alpha)
    class(CIpt) <- "CIpt"
    CIpt$CIpar <- CI
    CIpt$CIpar$alen <- rep(mean(par("fin")) * CI$alen, len)
    if (is.null(CIpt$CIpar$col)) 
      CIpt$CIpar$col <- Col
    if (is.null(CIpt$CIpar$lwd)) 
      CIpt$CIpar$lwd <- 1
    if (is.null(CIpt$CIpar$lty)) 
      CIpt$CIpar$lty <- 1
    CIpt$CIpar <- setdots(CIpt$CIpar, len)
    addCI <- function(x.from, y.from, z.from, x.to, y.to, 
                      z.to, tr) {
      tr$x.from <- cbind(tr$x.from, x.from)
      tr$y.from <- cbind(tr$y.from, y.from)
      tr$z.from <- cbind(tr$z.from, z.from)
      tr$x.to <- cbind(tr$x.to, x.to)
      tr$y.to <- cbind(tr$y.to, y.to)
      tr$z.to <- cbind(tr$z.to, z.to)
      tr
    }
    if (!is.null(CI$x)) 
      CIpt <- addCI(x - CI$x[, 1], y, z, x + CI$x[, 2], 
                    y, z, CIpt)
    if (!is.null(CI$y)) 
      CIpt <- addCI(x, y - CI$y[, 1], z, x, y + CI$y[, 
                                                     2], z, CIpt)
    if (!is.null(CI$z)) 
      CIpt <- addCI(x, y, z - CI$z[, 1], x, y, z + CI$z[, 
                                                        2], CIpt)
    nCI <- ncol(CIpt$x.from)
    CIpt$nCI <- rep(nCI, len)
    if (nCI < 3) {
      n0 <- matrix(nrow = len, ncol = 3 - nCI, data = NA)
      CIpt$x.from <- cbind(CIpt$x.from, n0)
      CIpt$y.from <- cbind(CIpt$y.from, n0)
      CIpt$z.from <- cbind(CIpt$z.from, n0)
      CIpt$x.to <- cbind(CIpt$x.to, n0)
      CIpt$y.to <- cbind(CIpt$y.to, n0)
      CIpt$z.to <- cbind(CIpt$z.to, n0)
    }
    CIpt$dopoints <- TRUE
  }
  else if (dopoints) {
    ii <- which(!is.na(x) & !is.na(y) & !is.na(z))
    pt <- list(x.mid = x[ii], y.mid = y[ii], z.mid = z[ii], 
               col = Col[ii], pch = rep(pch, length.out = length(ii)), 
               cex = rep(cex, length.out = length(ii)), bg = rep(bg, 
                                                                 length.out = length(ii)), alpha = alpha)
    class(pt) <- "pt"
  }
  Proj <- project(x, y, z, plist)
  if (!is.null(CI)) 
    CIpt$proj <- Proj
  else if (!is.null(pt)) 
    pt$proj <- Proj
  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
                         dot$clog, type = "scatter3D")
  plist <- plot.struct.3D(plist, pt = pt, CIpt = CIpt, poly = spoly, 
                          segm = sseg, plot = plot)
  setplist(plist)
  invisible(plist$mat)
}

