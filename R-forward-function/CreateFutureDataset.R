create_dataset <- function(path_to_data, endTime="2020-09-07",predTime=NULL) {
  t0 <- as.Date("2019-12-31")
  tf <- as.Date(endTime)
  data <-  as.data.frame(read.csv(file.path(path_to_data)))
  
  if (!is.null(predTime)) {
    predTime <- as.Date(predTime)
    add_time <- seq(as.numeric(tf)+1,as.numeric(as.Date(predTime)))-as.numeric(t0)
    add_reports <- rep(NA,length(add_time))
    add_deaths <- rep(NA,length(add_time))
    add_hospital <- rep(NA,length(add_time))
    add_ICU <- rep(NA,length(add_time))
    add_ventilator <- rep(NA,length(add_time))
    add_data <- cbind(add_time,add_reports,add_deaths,add_hospital,add_ICU,add_ventilator)
    colnames(add_data) <- c("time","reports","deaths","hospital","ICU","ventilator")
    data <- rbind(data,add_data)
    data <- as.data.frame(data)
    
  }
  
  return(data)
}