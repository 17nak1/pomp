
create_dataset <- function(endTime="2020-05-04",predTime=NULL) {

  t0 <- as.Date("2019-12-31")
  tf <- as.Date(endTime)
  dat <- read.csv("../samples/ON.csv", header=TRUE)
  dat$Date <- as.Date(dat$Date)

  delay <- 0
  dat <- subset(dat, Date>as.Date("2020-01-01")&Date<=tf-delay)
  dat$Date <- as.numeric(dat$Dat)
  dat$Date <- dat$Date - as.numeric(t0)

  time <- seq(1,as.numeric(tf)-as.numeric(t0),by=1)
  reports <- rep(0,length(time))
  reports[dat$Date] <- dat$Cases
  reports[(length(time)-delay):length(time)] <- NA


  dat <- read.csv("../samples/covidtesting.csv", header=TRUE)
  dat$Reported.Date <- as.Date(dat$Reported.Date)
  dat <- subset(dat, Reported.Date>as.Date("2020-03-17") & Reported.Date<=tf)
  dat$Reported.Date <- as.numeric(as.Date(dat$Reported.Date))-as.numeric(t0)
  dat <- subset(dat, select=c("Reported.Date","Deaths",
                              "Number.of.patients.hospitalized.with.COVID.19",
                              "Number.of.patients.in.ICU.with.COVID.19",
                              "Number.of.patients.in.ICU.on.a.ventilator.with.COVID.19"))
  colnames(dat) <- c("Date","Deaths","Hospital","ICU","Ventilator")

  total_deaths <- rep(0,length(time))
  total_deaths[dat$Date] <- dat$Deaths
  deaths <- rep(0,length(time))

  for (i in 2:length(time)) {
    deaths[i] <- total_deaths[i]-total_deaths[i-1]
  }

  hospital <- rep(NA,length(time))
  hospital[dat$Date] <- dat$Hospital

  ICU <- rep(NA,length(time))
  ICU[dat$Date] <- dat$ICU

  ventilator <- rep(NA,length(time))
  ventilator[dat$Date] <- dat$Ventilator

  data <- data.frame(time,reports,deaths,hospital,ICU,ventilator)
  colnames(data) <- c("time","reports","deaths","hospital","ICU","ventilator")

  ind <- which(data$time < 23)
  data[ind,c("hospital","ICU","ventilator")] <- 0

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
