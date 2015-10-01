replaceTabs <-
function(text, width=8)
{
  text <- as.character(text)
  retval <- sapply(text, replaceTabs.inner) 
  names(retval) <- names(text)
  retval
}
