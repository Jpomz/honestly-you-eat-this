# Functions for synthesising food webs from a registry of 
# known trophic links
#Gray et al. 2015
.StripWhitespace <- function(v)
{
  # Returns v with leading and trailing whitespace removed
  return (gsub('^\\s+|\\s+$', '', v, perl=TRUE))
}

.AddRegistryLinks <- function(synthesized, registry, res.col, con.col, 
                              res.method, con.method, reference.cols, 
                              debug)
{
  stopifnot(res.col %in% colnames(synthesized))
  stopifnot(res.col %in% colnames(registry))
  stopifnot(con.col %in% colnames(synthesized))
  stopifnot(con.col %in% colnames(registry))
  
  if(any(!reference.cols %in% colnames(registry)))
  {
    missing <- setdiff(reference.cols, colnames(registry))
    stop('Missing reference columns: ', paste(missing, collapse=','))
  }
  
  # Match on multiple columns by pasting
  # http://www.r-bloggers.com/identifying-records-in-data-frame-a-that-are-not-contained-in-data-frame-b-%E2%80%93-a-comparison/
  # Strip whitespace and convert to lowercase
  a <- lapply(synthesized[,c(res.col,con.col)], .StripWhitespace)
  a <- lapply(a, tolower)
  a <- do.call('paste', a)
  
  b <- lapply(registry[,c(res.col,con.col)], .StripWhitespace)
  b <- lapply(b, tolower)
  b <- do.call('paste', b)
  
  Debug <- function(...) if(debug) cat(...)
  Debug('Matching on resource', res.method, 'consumer', con.method, '\n')
  
  # Match where both resource and consumer columns are non-empty and a link 
  # does not already exist
  to.match <- ''!=synthesized[,res.col] & ''!=synthesized[,con.col] & 
    is.na(synthesized$N.records)
  for(row in which(to.match))
  {
    N <- sum(a[row]==b)
    if(N>0)
    {
      synthesized$res.method[row] <- res.method
      synthesized$con.method[row] <- con.method
      synthesized$N.records[row] <- N
      for(col in reference.cols)
      {
        synthesized[,col] <- paste(sort(unique(registry[a[row]==b,col])), collapse=',')
      }
    }
  }
  
  return (synthesized)
}

.AllPossible <- function(nodes, exclude.cols=NULL)
{
  # A data.frame of all possible links
  resource <- consumer <- nodes[,setdiff(colnames(nodes), exclude.cols),drop=FALSE]
  colnames(resource) <- paste0('res.', tolower(colnames(resource)))
  colnames(resource)['res.node'==colnames(resource)] <- 'resource'
  colnames(consumer) <- paste0('con.', tolower(colnames(consumer)))
  colnames(consumer)['con.node'==colnames(consumer)] <- 'consumer'
  
  synthesized <- data.frame(resource[rep(1:nrow(resource), each=nrow(resource)),,drop=FALSE], 
                            consumer[rep(1:nrow(resource), times=nrow(resource)),,drop=FALSE])
  return (synthesized)
}
.ARWorker <- function(synthesized, registry, methods, reference.cols, debug)
{
  # Add links for each level in methods
  for(outer in methods)
  {
    for(inner in methods[1:which(methods==outer)])
    {
      synthesized <- .AddRegistryLinks(synthesized, registry, 
                                       ifelse('exact'==outer, 'resource', paste0('res.', outer)), 
                                       ifelse('exact'==inner, 'consumer', paste0('con.', inner)), 
                                       outer, inner, reference.cols, debug)
      
      if(inner!=outer)
      {
        synthesized <- .AddRegistryLinks(synthesized, registry, 
                                         ifelse('exact'==inner, 'resource', paste0('res.', inner)), 
                                         ifelse('exact'==outer, 'consumer', paste0('con.', outer)), 
                                         inner, outer, reference.cols, debug)
      }
    }
  }
  
  synthesized <- droplevels(synthesized[!is.na(synthesized$N.records),])
  rownames(synthesized) <- NULL
  return (synthesized)
}

WebBuilder <- function(nodes, registry, 
                       methods=c('exact','genus','subfamily','family','order','class'),
                       reference.cols=c('linkevidence', 'source.id'), debug=FALSE)
{
  # nodes - a data.frame containing columns 'node' (character - unique node
  #         names) and either 'minimum.method' the minumum method for that
  #         taxon or 'minimum.res.method' and 'minimum.con.method' - the
  #         minimum methods for that taxon as a resource and a consumer
  #         respectively. Values in 'minimum.method', 'minimum.res.method' and
  #         'minimum.con.method' must be in the 'methods' argument.
  # registry - a data.frame of known trophic interactions
  # methods - character vector of methods, in order
  # reference.cols - names of columns in registry that will be included in the
  #         returned data.frame
  
  # Checks
  stopifnot(0<nrow(nodes))
  stopifnot(nrow(nodes)==length(unique(nodes$node)))
  stopifnot(all(c('resource','consumer') %in% colnames(registry)))
  
  error.msg <- paste('nodes should contain either "minimum.method" or', 
                     '"minimum.res.method" and "minimum.con.method"')
  if('minimum.method' %in% colnames(nodes))
  {
    if(any(c('minimum.res.method', 'minimum.con.method') %in% colnames(nodes)))
    {
      stop(error.msg)
    }
    
    # Use minimum.method for both ends of the links
    nodes$minimum.res.method <- nodes$minimum.con.method <- nodes$minimum.method
  }
  else
  {
    if(!all(all(c('minimum.res.method', 'minimum.con.method') %in% colnames(nodes))))
    {
      stop(error.msg)
    }
    # Use the provided minimum.res.method and minimum.con.method
  }
  
  stopifnot(all(nodes$minimum.res.method %in% methods))
  stopifnot(all(nodes$minimum.con.method %in% methods))
  
  synthesized <- .AllPossible(nodes,
                              c('minimum.res.method','minimum.con.method','minimum.method'))
  
  synthesized[,c('res.method','con.method')] <- ''
  synthesized$N.records <- NA
  synthesized[,reference.cols] <- ''
  
  links <- .ARWorker(synthesized, registry, methods, reference.cols, debug)
  
  # Set factor levels
  # What client has specified
  nodes$minimum.res.method <- factor(nodes$minimum.res.method, levels=methods)
  nodes$minimum.con.method <- factor(nodes$minimum.con.method, levels=methods)
  
  # What RegistryLinks found
  links$res.method <- factor(links$res.method, levels=methods)
  links$con.method <- factor(links$con.method, levels=methods)
  
  links <- links[as.integer(links$res.method)<=as.integer(nodes$minimum.res.method)[match(links$resource, nodes$node)] &
                   as.integer(links$con.method)<=as.integer(nodes$minimum.con.method)[match(links$consumer, nodes$node)],]
  return(droplevels(links))
}
