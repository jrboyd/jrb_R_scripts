textplot.matrix <-
function(object,
                            halign=c("center","left","right"),
                            valign=c("center","top","bottom"),
                            cex, cmar=2, rmar=0.5,
                            show.rownames=TRUE, show.colnames=TRUE,
                            hadj=1,
                            vadj=1,
                            mar= c(1,1,4,1)+0.1,
                            col.data=par("col"),
                            col.rownames=par("col"),
                            col.colnames=par("col"),
                            ... )
{
  
  if(is.vector(object))
    object <- t(as.matrix(object))
  else
    object <- as.matrix(object)
  
  # check dimensions of col.data, col.rownames, col.colnames
  if(length(col.data)==1)
    col.data <- matrix(col.data, nrow=nrow(object), ncol=ncol(object))
  else
    if( nrow(col.data)!=nrow(object) || ncol(col.data)!=ncol(object) )
      stop("Dimensions of 'col.data' do not match dimensions of 'object'.")
  
  if(length(col.rownames)==1)
    col.rownames <- rep(col.rownames, nrow(object))      
  
  if(length(col.colnames)==1)
    if(show.rownames)
      col.colnames <- rep(col.colnames, ncol(object)+1)
  else
    col.colnames <- rep(col.colnames, ncol(object))
  
  halign=match.arg(halign)
  valign=match.arg(valign)
  
  opar <- par()[c("mar","xpd","cex")]
  on.exit( par(opar) )
  par(mar=mar, xpd=FALSE )
  
  # setup plot area
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1), log = "", asp=NA)
  
  
  
  # add 'r-style' row and column labels if not present
  if( is.null(colnames(object) ) )
    colnames(object) <- paste( "[,", 1:ncol(object), "]", sep="" )
  if( is.null(rownames(object)) )
    rownames(object) <- paste( "[", 1:nrow(object), ",]", sep="")
  
  
  # extend the matrix to include row and column labels
  if( show.rownames )
  {
    object <- cbind( rownames(object), object )
    col.data <- cbind( col.rownames, col.data )
    
  }
  if( show.colnames )
  {
    object <- rbind( colnames(object), object )
    col.data <- rbind( col.colnames, col.data )
  }
  
  # set the character size
  if( missing(cex) )
  {
    cex <- 1.0
    lastloop <- FALSE
  }
  else
  {
    lastloop <- TRUE
  }
  
  for (i in 1:20)
  {
    oldcex <- cex
    
    width  <- sum(
      apply( object, 2,
             function(x) max(strwidth(x,cex=cex) ) )
    ) +
      strwidth('M', cex=cex) * cmar * ncol(object)
    
    height <- strheight('M', cex=cex) * nrow(object) * (1 + rmar)
    
    if(lastloop) break
    
    cex <- cex / max(width,height)
    
    if (abs(oldcex - cex) < 0.001)
    {
      lastloop <- TRUE
    }
  }
  
  
  # compute the individual row and column heights
  rowheight<-strheight("W",cex=cex) * (1 + rmar)
  colwidth<- apply( object, 2, function(XX) max(strwidth(XX, cex=cex)) ) +
    strwidth("W")*cmar
  
  
  width  <- sum(colwidth)
  height <- rowheight*nrow(object)
  
  # setup x alignment
  if(halign=="left")
    xpos <- 0
  else if(halign=="center")
    xpos <- 0 + (1-width)/2
  else #if(halign=="right")
    xpos <- 0 + (1-width)
  
  # setup y alignment
  if(valign=="top")
    ypos <- 1
  else if (valign=="center")
    ypos <- 1 - (1-height)/2
  else #if (valign=="bottom")
    ypos <- 0 + height
  #xpos = xpos - .1
  x <- xpos
  y <- ypos
  
  # iterate across elements, plotting them
  xpos<-x
  for(i in 1:ncol(object)) {
    xpos <- xpos + hadj*colwidth[i]
    for(j in 1:nrow(object)) {
      ypos<-y-(j-1)*rowheight
      if( (show.rownames && i==1) || (show.colnames && j==1) )
        text(xpos, ypos, object[j,i], adj=c(hadj,vadj), cex=cex, font=2,
             col=col.data[j,i], ... )
      else
        text(xpos, ypos, object[j,i], adj=c(hadj,vadj), cex=cex, font=1,
             col=col.data[j,i], ... )
    }
    xpos <- xpos + (1-hadj)*colwidth[i]
  }
  
  par(opar)
}
