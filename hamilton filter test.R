#---------------------------------------------------------------------
#use this function to get the same data that Matt Brigida uses
require(XML)
require(xts)
library(EIAdata)

getMonEIA <- function(ID, key) {
  
  ID <- unlist(strsplit(ID, ";"))
  key <- unlist(strsplit(key, ";"))
  
  url <- paste("http://api.eia.gov/series?series_id=", ID, "&api_key=", key,
               "&out=xml", sep = "")
  
  doc <- xmlParse(file = url, isURL = TRUE)
  
  df <- xmlToDataFrame(nodes = getNodeSet(doc, "//data/row"))
  
  df <- plyr::arrange(df, df$date)
  
  date <- as.Date(paste(as.character(levels(df[, 1]))[df[, 1]], "01", sep = ""),
                  "%Y%m%d")
  values <- as.numeric(levels(df[, -1]))[df[, -1]]
  
  tmp <- data.frame(date=date,x=values)
  names(tmp) <- c('date',paste(ID))
  return(tmp)
  #xts_data <- xts(values, order.by = date)
  #names(xts_data) <- sapply(strsplit(ID, "-"), paste, collapse = ".")
  
  #assign(sapply(strsplit(ID, "-"), paste, collapse = "."), xts_data, envir = .GlobalEnv)
}
#----------------------------------------------------------------------#


#-----------------------------------------------------------------------------------
#in order to pull the data we will need to provide a key.  The EIA will send you a key for
# accessing their open data API....details can be found here:
#http://complete-markets.com/2014/01/r-functions-to-interact-with-the-eias-application-programming-interface-api/

key <- '80DAEF6FDCF88254A3D4E2F5F4CBF99C'
#-----------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
# use the function above to get monthly data from the EIA on natural gas prices and West
# Texas Crude Oil...then filter the series so they are the same length

lng <- getMonEIA(ID='NG.RNGWHHD.M',key=key)
names(lng) <- c('date','lng')
oil <- getMonEIA(ID='PET.RWTC.M',key=key)
names(oil) <- c('date','oil')

start <- max(c(min(lng$date),min(oil$date)))
end <- min(max(lng$date),max(oil$date))

lnoil <- log(as.vector(oil$oil[which(oil$date>=start & oil$date<=end)]))
lnng <- log(as.vector(lng$lng[which(lng$date>=start & lng$date<=end)]))
#-----------------------------------------------------------------------------------


################################################################################
################################################################################
################################################################################
# The Hamilton Filter

#need to write a function to do the calculation then optimize that function

#NOTE: this is not all that general.  It is set up to deal with a simple univariate
# regression model where:

# y(t) = alpha1 + beta1*x(t) if we are in regime 1
# y(t) = alpha2 + beta2*x(t) if we are in regime 2

# we need to generalize this later to be amenable to different regression forms

mrs.est <- function(theta,x,y){
  
  
  y=lnng
  x=lnoil
  alpha1 <- 1
  alpha2 <- 2
  alpha3 <- 3
  alpha4 <- 4
  alpha5 <- 5
  alpha6 <- 6
  
  p11 <- 1/(1+exp(-0.5))
  p22 <- 1/(1+exp(-1))
  
  #in order to make inference about what state we are in in period t we need the conditional
  # densities given the information set through t-1
  f1 <- (1/(alpha5*sqrt(2*pi)))*exp(-((y-alpha1-(alpha3*x))^2)/(2*(alpha5^2)))
  f2 <- (1/(alpha6*sqrt(2*pi)))*exp(-((y-alpha2-(alpha4*x))^2)/(2*(alpha6^2)))
 
  #p_s
  f <- matrix(cbind(f1,f2),nc=2)
  
  #S.forecast is the state value looking forward conditional on info up to time t
  #S.inf is the updated state value
  S.forecast <- rep(0,2*length(y))
  S.forecast <- matrix(S.forecast,nrow=(length(y)),ncol=2)
  
  S.inf <- S.forecast
  o.v <- c(1,1)
  
  #p_miu
  P<-matrix(c(p11,(1-p11),(1-p22),p22),nr=2,nc=2)
 
  model.lik <- rep(0,length(y))
  
  S.inf[1,] <- (c(p11,p22)*f[1,])/as.numeric(o.v %*% (c(p11,p22)*f[1,]))
  
  for(t in 1:(length(y)-1)){
    #in time t we first make our forecast of the state in t+1 based on the
    # data up to time t, then we update that forecast based on the data
    # available in t+1
    S.forecast[t+1,] <- P%*%S.inf[t,]
    S.inf[t+1,] <- (S.forecast[t+1,]*f[t+1,])/(S.forecast[t+1,] %*% f[t+1,])
    model.lik[t+1] <- o.v%*%(S.forecast[t+1,]*f[t+1,])
  }
  
  logl <- sum(log(model.lik[2:length(model.lik)]))
  return(-logl)
}
################################################################################
################################################################################
################################################################################



