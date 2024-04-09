#' Function to estimate occupancy trends
#'@param z Data frame output from \code{fit_occ}, containing estimated occupancy index and standard errors.
#'@param startyear Optional. First year to calculate the trend from. Default will select the first year of the index.
#'@param endyear Optional. Year to calculate the trend up to. Default will select the last year of the index.
#'@param predictions Option to include predicted values in the output e.g. for plotting.
#'@export

# Estimate trends
fit_trend <- function(z,
                      startyear=NULL,
                      endyear=NULL,
                      predictions = FALSE,
                      regionalise = NULL){

  if(quantile(z$psiA, .9, na.rm=TRUE) < .9){
    z <- subset(z, psiA < .99)
  }

  if(is.null(startyear)) startyear <- min(z$Year)
  if(is.null(endyear)) endyear <- max(z$Year)

  z <- subset(z, Year %in% startyear:endyear)
  trds <- NULL

  if(!is.null(regionalise)){
    for(r in unique(z[,get(regionalise)])){
      zt <- z[get(regionalise) == r]
      if(nrow(zt)>=3 && (min(zt$Year)-startyear) <= 5){
        fit1 <- try(glm(psiA ~ Year, data = zt, weights=1/psiA_SD,
                        family=binomial(link="logit")), silent=TRUE)
        if(class(fit1)[1] != "try-error"){
          pred <- data.frame(Year=startyear:endyear)
          fit1.pred <- predict(fit1,
                               newdata = pred,
                               weights = 1/zt$psiA_SD,
                               se.fit = TRUE)
          pred$pred <- expit(fit1.pred$fit)
          pred$upp <- expit(fit1.pred$fit + 1.96 * fit1.pred$se.fit)
          pred$low <- expit(fit1.pred$fit - 1.96 * fit1.pred$se.fit)
          pred[, regionalise] <- r

          trd_ci1 <- round(pcfunc2(head(pred$upp,1),tail(pred$low,1)),3)
          trd_ci2 <- round(pcfunc2(head(pred$low,1),tail(pred$upp,1)),3)
          trd_ci <- sort(c(trd_ci1,trd_ci2))
          trd_se <- (trd_ci[2]-trd_ci[1])/1.96/2

          trd <- data.frame(Startyear = startyear,
                            Endyear = endyear,
                            nyears = nrow(zt),
                            minyear = min(zt$Year),
                            trd = round(pcfunc(pred$pred),3),
                            trd_se = trd_se,
                            trd_ci1 = trd_ci1,
                            trd_ci2 = trd_ci2,
                            coef = coef(summary(fit1))[2,1],
                            SE = coef(summary(fit1))[2,2],
                            pval = coef(summary(fit1))[2,4])

          trd$sig <- if(!is.nan(trd$pval)){sigfunc(trd$pval)} else {NA}

          trd[, regionalise] <- r

          if(predictions){
            if(!is.null(trds)){
              trds <- list(trd = rbind(trds$trd, trd),
                           predictions = rbind(trds$predictions,
                                               pred))
            } else {
              trds <- list(trd = trd,
                           predictions = pred)
            }
          } else {
            trds <- rbind(trds, trd)
          }

        }}

    }

  } else {

    if(nrow(z)>=3 && (min(z$Year)-startyear) <= 5){
      fit1 <- try(glm(psiA ~ Year, data = z, weights=1/psiA_SD,
                      family=binomial(link="logit")), silent=TRUE)
      if(class(fit1)[1] != "try-error"){
        fit1.pred <- predict(fit1,
                             newdata=data.frame(Year=startyear:endyear),
                             weights=1/z$psiA_SD,
                             se.fit=TRUE)
        upp <- expit(fit1.pred$fit + 1.96 * fit1.pred$se.fit)
        pred <- expit(fit1.pred$fit)
        low <- expit(fit1.pred$fit - 1.96 * fit1.pred$se.fit)
        trd_ci1 <- round(pcfunc2(head(upp,1),tail(low,1)),3)
        trd_ci2 <- round(pcfunc2(head(low,1),tail(upp,1)),3)
        trd_ci <- sort(c(trd_ci1,trd_ci2))
        trd_se <- (trd_ci[2]-trd_ci[1])/1.96/2

        trds <- data.frame(Startyear=startyear,
                           Endyear=endyear,
                           nyears=nrow(z),
                           minyear=min(z$Year),
                           trd = round(pcfunc(pred),3),
                           trd_se = trd_se,
                           trd_ci1 = trd_ci1,
                           trd_ci2 = trd_ci2,
                           coef = coef(summary(fit1))[2,1],
                           SE = coef(summary(fit1))[2,2],
                           pval = coef(summary(fit1))[2,4])

        trds$sig <- if(!is.nan(trds$pval)){sigfunc(trds$pval)} else {NA}

        if(predictions){
          trds <- list(trd = trds,
                       predictions = data.frame(Year=startyear:endyear,
                                                pred = pred,
                                                low = low,
                                                upp = upp))
        }

      }}
  }
  return(trds)
}


# Trend utility functions
expit <- function(x){ exp(x)/(1+exp(x)) }
pcfunc <- function(pred){(tail(pred,1)-head(pred,1))/head(pred,1)*100}
pcfunc2 <- function(pred1,pred2){(pred2-pred1)/pred1*100}

sigfunc <- function(x){y <- ""; if(x <= 0.05) y <- "*";
  if(x <= 0.01) y <- "**";
  if(x <= 0.001) y <- "***";
  return(y)}

