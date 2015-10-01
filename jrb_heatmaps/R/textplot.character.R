textplot.character <-
function (object,
                                halign = c("center", "left", "right"),
                                valign = c("center", "top", "bottom"),
                                cex, fixed.width=TRUE,
                                cspace=1,
                                lspace=1,
                                mar=c(0,0,3,0)+0.1,
                                tab.width=8,
                                ...)
{
  object <- paste(object,collapse="\n",sep="")
  object <- replaceTabs(object, width=tab.width)
  
  halign = match.arg(halign)
  valign = match.arg(valign)
  plot.new()
  
  opar <- par()[c("mar","xpd","cex","family")]
  on.exit( par(opar) )
  
  par(mar=mar,xpd=FALSE )
  if(fixed.width)
    par(family="mono")
  
  plot.window(xlim = c(0, 1), ylim = c(0, 1), log = "", asp = NA)
  
  slist   <- unlist(lapply(object, function(x) strsplit(x,'\n')))
  slist   <- lapply(slist, function(x) unlist(strsplit(x,'')))
  
  slen    <- sapply(slist, length)
  slines  <- length(slist)
  
  if (missing(cex))
  {
    lastloop <- FALSE
    cex <- 1
  }
  else
    lastloop <- TRUE
  
  
  for (i in 1:20)
  {
    oldcex <- cex
    #cat("cex=",cex,"\n")
    #cat("i=",i,"\n")
    #cat("calculating width...")
    cwidth  <- max(sapply(unlist(slist), strwidth,  cex=cex)) * cspace
    #cat("done.\n")
    #cat("calculating height...")
    cheight <- max(sapply(unlist(slist), strheight, cex=cex)) * ( lspace + 0.5 )
    #cat("done.\n")
    
    width <- strwidth(object, cex=cex)
    height <- strheight(object, cex=cex)
    
    if(lastloop) break
    
    cex <- cex  / max(width, height)
    
    if (abs(oldcex - cex) < 0.001)
    {
      lastloop <- TRUE
    }
    
  }
  
  if (halign == "left")
    xpos <- 0
  else if (halign == "center")
    xpos <- 0 + (1 - width)/2
  else xpos <- 0 + (1 - width)
  
  if (valign == "top")
    ypos <- 1
  else if (valign == "center")
    ypos <- 1 - (1 - height)/2
  else ypos <- 1 - (1 - height)
  
  text(x=xpos, y=ypos, labels=object, adj=c(0,1),
       cex=cex, ...)
  
  par(opar)
  invisible(cex)
}
