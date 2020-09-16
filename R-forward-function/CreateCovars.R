
create_covars <- function(endTime="2020-09-07",predTime=NULL) {
  
  t0 <- as.Date("2019-12-31")
  tf <- as.Date(endTime)

  dat <- read.csv("covidtesting.csv", header=TRUE)
  dat$Reported.Date <- as.Date(dat$Reported.Date)
  dat <- subset(dat, Reported.Date>as.Date("2020-02-03") & Reported.Date<=tf)
  dat$Reported.Date <- as.numeric(as.Date(dat$Reported.Date))-as.numeric(t0)
  dat <- subset(dat, select=c("Reported.Date",
                              "Total.patients.approved.for.testing.as.of.Reporting.Date",
                              "Under.Investigation"))
  colnames(dat) <- c("Date","TotalTests","Pending")

  time <- seq(1,dat$Date[length(dat$Date)])
  total_tests <- approx(dat$Date,dat$TotalTests,xout=time,method="linear")$y
  total_tests <- ceiling(total_tests)
  
  
  pending <- rep(0,length(time))
  pending[dat$Date] <- dat$Pending 
  
  tests <- rep(NA,length(time))
  
  for (i in 2:length(time)) {
    tests[i] <- max((total_tests[i]) - (total_tests[i-1]),0)
  }
  tests <- ceiling(filter(tests, rep(1/2,2)) )
  tests[length(tests)] <- tests[length(tests)-1]
  
  
  
  data <- data.frame(time,ceiling(as.numeric(tests)))
  colnames(data) <- c("time","tests")
  
  
  
  
  
  if (!is.null(predTime)) {
    predTime <- as.Date(predTime)
    add_time <- seq(as.numeric(tf)+1,as.numeric(as.Date(predTime)))-as.numeric(t0)
    add_tests <- rep(10e3,length(add_time))
    add_data <- cbind(add_time,add_tests)
    colnames(add_data) <- c("time","tests")
    data <- rbind(data,add_data)
    data <- as.data.frame(data)
    
  }
  
  return(data)
}
