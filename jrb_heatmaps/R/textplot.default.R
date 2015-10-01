textplot.default <-
function(object,
                             halign=c("center","left","right"),
                             valign=c("center","top","bottom"),
                             cex, ... )
{
  
  if (is.matrix(object) || (is.vector(object) && length(object)>1) )
    return(textplot.matrix(object, halign, valign, cex, ... ))
  
  halign <- match.arg(halign)
  valign <- match.arg(valign)
  
  textplot.character(object, halign,  valign, cex, ...)
}
