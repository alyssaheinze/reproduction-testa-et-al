### Generate weights



get_acte <- function(df = acte, z = c("Selection", "Control")){
  
  #include weights
  
  acte_weights <- df |>
    mutate(weights = case_when(
      C=="Choice" & avoid01 == 0 & D_ch == "Control"~ weight$select_weights,
      C=="Choice" & avoid01 == 1 & D_ch == "Control"~ weight$avoid_weights,
      TRUE ~ 1)
    )
  
  ### Get weighted means and weighted SEs 
  weighted_means_sds <- acte_weights |> 
    filter(!is.na(treatment),
           !is.na(dv_pca_metoo)) |>
    group_by(treatment) |>
    summarise(avg = sum(weights*dv_pca_metoo, 
                        na.rm = T)/sum(weights, na.rm =T),
              se = weighted.var(dv_pca_metoo, 
                                weights, na.rm = T))
  
  ## Get the first part of the expression and its SE 
  ## This is just the mean and se for the select01 group
  y_select <- df |>
    filter(C == "Choice") |>
    summarise(avg = mean(select01),
              se = sd(select01)/sqrt(n()))
  
  ### Delta Method
  avg_x <- weighted_means_sds %>%
    # filter(treatment != "Treatment")%>%
    # summarise(diff = sum(abs(avg)))%>%
    summarise(diff = avg[treatment == z[1]] - avg[treatment == z[2]]) %>%
    pull()
  
  x_se <- weighted_means_sds %>% 
    # filter(treatment != "Treatment")%>%
    # summarise(se = sqrt(sum(se)))%>%
    summarise(se = sqrt(sum(se[treatment == z[1]], se[treatment == z[2]]))) %>%
    pull()
  
  denom <- ifelse(z[1] == "Selection", y_select$avg, 1-y_select$avg)
  
  # estimate <- avg_x /y_select$avg
  estimate <- avg_x / denom
  
  ## function that we seek to estimate 
  ## pass this to the jacobian function which calculates the gradient 
  ratio <- function(x){
    x[1]/x[2]
  }
  
  ## Variance covariance matrix 
  vcov <- diag(c(x_se,y_select$se)^2)
  
  ## this function by default will output this in 
  ## transpose form, so take the transpose to match 
  ## discussion above
  # grad_g <- t(numDeriv::jacobian(func = ratio, c(avg_x, y_select$avg)))
  grad_g <- t(numDeriv::jacobian(func = ratio, c(avg_x, denom)))
  se_b <- sqrt(t(grad_g) %*% vcov %*%grad_g)
  
  lwr_b = estimate - 1.96*se_b 
  upp_b = estimate + 1.96*se_b
  
  ### SIDE: ALTERNATIVE FUNCTION DELTA METHOD
  # mvec <- c(x=avg_x, y=denom)
  # est <- car::deltaMethod(mvec,"x/y",vcov,level=.95)
  
  round(c(estimate = estimate, lwr_b = lwr_b, upp_b = upp_b),2)
  
}

# test function

get_acte(df = acte_weights, z = c("Selection", "Control")) # ACTE select
get_acte(df = acte_weights, z = c("Treatment", "Selection")) # ACTE avoid
get_acte(df = acte[acte$Partisanship=="Republicans",], # ACTE avoid - republican
         z = c("Selection", "Control")) 
get_acte(df = acte_weights[acte_weights$Partisanship=="Republicans",], # ACTE avoid - republican
         z = c("Treatment", "Selection")) 



