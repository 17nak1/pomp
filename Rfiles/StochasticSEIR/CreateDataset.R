

Biweekly=function(Data){
  n=nrow(Data)/2
  m=ncol(Data)
  mat = matrix(,n,m)
  for (i in 0:n - 1 ){
    mat[i + 1,]= rep(0, m)
    for (j in 1:2) {
      x = (2*i)+j
      mat[i+1,] = c(mat[i+1,]) + c(Data[x,])
    }
  }
  return(mat)
}

"Loading Datasets"
#daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
#datfile <- file.path(tempdir(),"twentycities.rda")
#download.file(daturl,destfile=datfile,mode="wb")
datfile <- "twentycities.rda"
load(datfile)

demog$town = factor(demog$town)
measles$town = factor(measles$town)

"creating City datasets"
for (names in c("London")) {
  tmp<- subset(demog, town == names)
  tmp<-tmp[,-1]
  tmp %>% subset(year>=1944 & year<1964) %>%
    summarize(
      time=seq(from=min(year),to=max(year),by=1/12),
      pop=predict(smooth.spline(x=year,y=pop),x=time)$y,
      birthrate=predict(smooth.spline(x=year+0.5,y=births),x=time-4)$y
    ) -> covar
  
  assign( paste0(names,"_covar"),covar)
}


"Cases"
for (names in c("London")) {
  tmp <- subset(measles, town == names)
  tmp %>%
    dcast(date~"cases", fun.aggregate = sum) %>%
    mutate(year=as.integer(format(date,"%Y"))) %>%
    subset(year>=1944 & year<1965) %>%
    mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
    subset(time>1944 & time<1965, select=c(time,cases)) -> tmp #Weekly incidence data
  
  ################## Removing the column with town names
  covar<- subset(demog, town == names)
  covar<-covar[,-1]
  ##################
  
  covar %>% subset(year>=1944 & year<1964) %>%
    summarize(
      time= tmp$time,
      
      birthrate=(predict(smooth.spline(x=year,y=births),x=time-4)$y)/52,
      
      pop=(predict(smooth.spline(x=year,y=pop),x=time)$y)/52 ) -> covar #weekly birth and pop data, birth adjusted for marternal immunity
  
  merge(tmp,covar, by = "time")->x
  
  #Converting to Bi-weeks
  x<- as.matrix(x)
  Bdat = Biweekly(x)# Biweekly data
  
  #here
  result = cbind(time=seq(from = 1944+14/365.25 , by = (14/365.25), length.out = 547),cases=Bdat[,2])
  result = as.data.frame(result)
  
  
  
  assign( paste0(names,"_BiData"),result)
  
}

London_BiData <- rbind(c(1944,NA),London_BiData)
setwd("~/Documents/R/StochasticSEIR")
write.csv(London_covar, file="London_covar.csv", row.names=F)
write.csv(London_BiData, "London_BiData.csv",row.names=F)

# rm(Bdat,coord,covar,demog,measles,result,tmp,x,datfile)