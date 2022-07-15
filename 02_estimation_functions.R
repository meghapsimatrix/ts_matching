
# simple matching ---------------------------------------------------------

match_them <- function(dat, 
                       ps_method = "nearest",
                       exact = NULL, 
                       ps = dat$ps_unit){
  
  
  m_out <- matchit(D ~ ps,  # the rhs doesn't matter I think 
                   caliper = .25,
                   method = ps_method,
                   exact = exact,
                   distance = ps,
                   data = dat)
  
  match_dat <- match.data(m_out)
  
  return(match_dat)
  
}