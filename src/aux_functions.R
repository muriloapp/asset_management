library("here")
library("xts")
library("PerformanceAnalytics")
library("qrmdata")


download_and_save_data <- function(){
  data("FTSE_const"); prices <- FTSE_const
  prices    <- xts(prices, order.by = index(prices))
  prices    <- prices["2012/2015"]
  endofweek <- xts::endpoints(prices, on = "weeks")
  prices    <- prices[endofweek, 1:52] # Take the first 52 stocks
  rets      <- PerformanceAnalytics::CalculateReturns(prices, method = "discrete")
  rets      <- rets[-1,]
  rets$DLG.L <- NULL # Trading begin after initial period
  rets$III.L <- NULL # Data issues
  save(rets, file = here("data", "FTSE_const_rets.rda"))
}