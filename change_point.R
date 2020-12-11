library(ISLR)
library(splines)
library(SplinesUtils)
library(tidyverse)


################################################################################
# Reading in data
################################################################################

file <- file.choose()
df <- read_csv(file)

# us-counties.csv
df.cases <- df %>% 
  select(date, county, cases) %>%
  group_by(date) %>%
  summarise(cases = sum(cases)) %>%
  mutate(daily_cases = cases  - dplyr::lag(cases, order_by = date))

# us.csv
df.cases <- df %>% 
  mutate(daily_cases = cases  - dplyr::lag(cases, order_by = date, default = first(cases)))

################################################################################
# First Spline Fitting Attempts
################################################################################

# spline fit
knot_1 <- 188425
degree_1 <- 1
splineTerm_1 <- "bs(df.cases$cases, degree = 1, knots = 188425)"
#basis_1 <- splines::bs(df.cases$cases, degree = degree, knots = knot)
model_1 <- lm(
  daily_cases ~ bs(df.cases$cases, degree = 1, knots = 188425), 
  data = df.cases
  )

knot_2 <- 277279
degree_2 <- 2
splineTerm_2 <- "bs(df.cases$cases, degree = 2, knots = 277279)"
#basis_2 <- splines::bs(df.cases$cases, degree = degree, knots = knot)
model_2 <- lm(
  daily_cases ~ bs(df.cases$cases, degree = 2, knots = 277279), 
  data = df.cases
)

# plot with spline
ggplot(data = df.cases, mapping = aes(x = cases, y = daily_cases)) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ splines::bs(x, degree = degree_1, knots = knot_1))
ggplot(data = df.cases, mapping = aes(x = cases, y = daily_cases)) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ splines::bs(x, degree = degree_2, knots = knot_2))

# interpret coefficients
RegBsplineAsPiecePoly(model_1, splineTerm_1)
RegBsplineAsPiecePoly(model_2, splineTerm_2)

## direct for order one
summary(model)$coefficients[2, 1] / (knot - min(df.cases$cases))
(summary(model)$coefficients[3, 1] - summary(model)$coefficients[2,1]) / (max(df.cases$cases - knot))

generateSplineTerm <- function(x_arg, basis, ...){
  knots <- attr(basis, "knots")
  knots_arg <- paste("knots = c(", paste(knots, collapse = ", "), ")")
  degree <- attr(basis, "degree")
  degree_arg <- paste("degree = ", degree)
  term <- paste("bs(", x_arg, ", ", knots_arg, ", ", degree_arg, ")", sep = "")
  return(term)  
}

################################################################################
# 10/7/20 
################################################################################

# Implement truncated power function basis and fit model
## code example 1 from paper
U <- seq(0, 1, length=100)
D <- 2 # degree of series
K <- 3 # number of knots
knots <- (1:K) / (K+1) # creates a series of K equidistant knots
knots <- c(0, knots)
X1 <- outer (U, 1:D, "^") #a s in eq.2
X2 <- outer (U, knots,">") * outer (U, knots, "-")^D # as in eq.3
Bt <- cbind (X1, X2)
matplot(U, Bt, type = "l", lwd=2, col= 1:ncol(Bt))

## knots are blind, equally spaced on full interval
U <- 1:length(df.cases$cases) #seq(df.cases$cases[1], tail(df.cases$cases, n=1), length = length(df.cases$cases))
D <- 3
K <- 10
knots <- (1:K)/(K+1) * 260
X1 <- outer(U, D, "^")
X2 <- outer(U, knots,">") * outer(U, knots, "-")^D
Bt <- cbind (X1, X2)
lm <- lm(daily_cases ~ Bt, data = df.cases)
df.fit <- tibble(
  date = df.cases$date,
  cases = df.cases$cases, 
  daily_cases = df.cases$daily_cases,
  estimated_daily_cases = cbind(1, Bt) %*% t(t(lm$coefficients))
)
df.fit <- df.fit %>% mutate(
  estimated_cases = cumsum(estimated_daily_cases),
  basis_cases = seq(cases[1], tail(cases, n = 1), length = length(cases))
)
ggplot(data = df.fit, mapping = aes(x = 1:260, y = daily_cases)) + 
  geom_point() +
  geom_line(aes(1:260, estimated_daily_cases))

## Loop through dumb knots until perfect fit on one knot, get AIC|BIC
AIC <- c()
BIC <- c()
U <- 1:length(df.cases$cases)
D <- 1
K <- 1:(length(df.cases$cases) - 3)
for(k in K){
  knots <- (1:k)/(k+1) * 260 #(tail(df.cases$cases, n=1) - df.cases$cases[1]) + df.cases$cases[1]
  X1 <- outer(U, D, "^")
  X2 <- outer(U, knots,">") * outer(U, knots, "-")^D
  Bt <- cbind (X1, X2)
  lm <- lm(daily_cases ~ Bt, data = df.cases)
  AIC <- c(AIC, AIC(lm))
  BIC <- c(BIC, BIC(lm))
}
plot(K, AIC)
plot(K, BIC)
### a guess when D = 1, K = 5
knots <-c(55, 77, 140, 175, 225)
X1 <- outer(U, D, "^")
X2 <- outer(U, knots,">") * outer(U, knots, "-")^D
Bt <- cbind (X1, X2)
lm <- lm(daily_cases ~ Bt, data = df.cases)
#### compare to first hundred
AIC(lm)
min(AIC[1:100]) # AIC on guess is worse
BIC(lm)
min(BIC[1:100]) # BIC on guess is best in class for knots up to 100

# Make knots smart, assign to dates, grid search
map_knots <- function(file, knots, depth = 1, ref = knots[[1]][1] - 1, knot_vec = c()){
  for (knot in knots[[depth]]){
    if (knot > ref && depth < length(knots)){
      map_knots(file, knots, depth = depth + 1, ref = knot, knot_vec = c(knot_vec, knot))
    }
    if (knot > ref && depth == length(knots)){
      readr::write_lines(toString(c(knot_vec, knot)), file, append = TRUE)
    }
  }
}


### See if method works... Test on reduced knot region
n <- seq(1, 260, 7)
P <- 1:6
for (p in P){
  knots <- rep(list(n), p)
  filename <- paste("/Users/MisterGrinvalds/Repos/change-point/knots_", toString(p), ".txt", sep = "")
  map_knots(file = filename, knots)
}
### generate models from knot specs
U <- 1:length(df.cases$cases)
D <- 1
knots <- read_csv("/Users/MisterGrinvalds/Repos/change-point/knots_5.txt", col_names = FALSE)
AIC <- c()
BIC <- c()
for (i in 1:nrow(knots)){
  knot <- knots %>% slice(i) %>% unlist()
  X1 <- outer(U, D, "^")
  X2 <- outer(U, knot,">") * outer(U, knot, "-")^D
  Bt <- cbind (X1, X2)
  lm <- lm(daily_cases ~ Bt, data = df.cases)
  AIC <- c(AIC, AIC(lm))
  BIC <- c(BIC, BIC(lm))
}
plot(1:nrow(knots), AIC)
which.min(AIC)
knots[which.min(AIC),]
plot(1:nrow(knots), BIC)
which.min(BIC)
knots[which.min(BIC),]

# Apply to June 18 data
## Results
ggplot(data = df.cases, mapping = aes(x = 1:149, y = daily_cases)) + geom_point()

# Test on smaller region
n <- seq(1, 149, 7)
P <- 1:6
for (p in P){
  knots <- rep(list(n), p)
  filename <- paste("/Users/MisterGrinvalds/Repos/change-point/knots_june18_", toString(p), ".txt", sep = "")
  map_knots(file = filename, knots)
}
## generate models from knot specs
AIC_min <- c()
BIC_min <- c()
for (p in P){
  AIC <- c()
  BIC <- c()
  knots <- read_csv(paste("/Users/MisterGrinvalds/Repos/change-point/knots_june18_", toString(p), ".txt", sep = ''), col_names = FALSE)
  for (i in 1:nrow(knots)){
    U <- 1:length(df.cases$cases)
    D <- 1
    knot <- knots %>% slice(i) %>% unlist()
    X1 <- outer(U, D, "^")
    X2 <- outer(U, knot,">") * outer(U, knot, "-")^D
    Bt <- cbind (X1, X2)
    lm <- lm(daily_cases ~ Bt, data = df.cases)
    AIC <- c(AIC, AIC(lm))
    BIC <- c(BIC, BIC(lm))
    
  }
  AIC_min <- rbind(AIC_min, c(p, which.min(AIC), min(AIC)))
  BIC_min <- rbind(BIC_min, c(p, which.min(BIC), min(BIC)))
  print(qplot(1:nrow(knots), AIC))
  print(qplot(1:nrow(knots), BIC))
}
## Get best fitting model from BIC_min
knots <- read_csv("/Users/MisterGrinvalds/Repos/change-point/knots_june18_4.txt", col_names = FALSE)
U <- 1:length(df.cases$cases)
D <- 1
knot <- knots %>% slice(6387) %>% unlist()
X1 <- outer(U, D, "^")
X2 <- outer(U, knot,">") * outer(U, knot, "-")^D
Bt <- cbind (X1, X2)
lm <- lm(daily_cases ~ Bt, data = df.cases)
df.fit <- tibble(
  date = df.cases$date,
  cases = df.cases$cases, 
  daily_cases = df.cases$daily_cases,
  estimated_daily_cases = cbind(1, Bt) %*% t(t(lm$coefficients))
)
ggplot(data = df.fit, mapping = aes(x = 1:149, y = daily_cases)) + 
  geom_point() +
  geom_line(aes(1:149, estimated_daily_cases))
ggplot(data = df.cases, mapping = aes(x = cases, y = daily_cases)) + 
  geom_point()

################################################################################
# 10/8/20
################################################################################

form_truncated_basis <- function(X, knots, D = 1){
  X1 <- outer(X, D, "^")
  X2 <- outer(X, knot,">") * outer(X, knot, "-")^D
  Bt <- cbind (X1, X2)
  return(Bt)
}


# analysis using `october-07` branch (260 observations)
ggplot(data = df.cases, mapping = aes(x = cases, y = daily_cases)) + geom_point()
## obtain knot map
K <- 30
P <- 1:6
for (p in P){
  knots_list <- list()
  for (knots in split(df.cases$cases, ceiling(seq_along(df.cases$cases)/(length(df.cases$cases)/p)))){
    knots_list <- append(knots_list, list(seq(min(knots), max(knots), length = ceiling(K / p))))
  }
  filename <- paste("/Users/MisterGrinvalds/Repos/change-point/knots_october07_", toString(p), ".txt", sep = "")
  map_knots(file = filename, knots_list)
}
## generate models
P <- 1:6
D <- 1:3
AIC_min <- c()
BIC_min <- c()
for (p in P) {
  for (d in D){
    AIC <- c()
    BIC <- c()
    knots <- read_csv(paste("/Users/MisterGrinvalds/Repos/change-point/knots_october07_", toString(p), ".txt", sep = ''), col_names = FALSE)
    for (i in 1:nrow(knots)){
      knot <- knots %>% slice(i) %>% unlist()
      Bt <- form_truncated_basis(df.cases$cases, knots = knot, D = d)
      lm <- lm(daily_cases ~ Bt, data = df.cases)
      AIC <- c(AIC, AIC(lm))
      BIC <- c(BIC, BIC(lm))
    }
    AIC_min <- rbind(AIC_min, c(p, d, which.min(AIC), min(AIC)))
    BIC_min <- rbind(BIC_min, c(p, d, which.min(BIC), min(BIC)))
    knot <- knots %>% slice(which.min(AIC)) %>% unlist()
    knot[is.na(knot)] <- 0
    Bt <- form_truncated_basis(df.cases$cases, knots = knot, D = d)
    lm <- lm(daily_cases ~ Bt, data = df.cases)
    df.fit <- tibble(
      cases = df.cases$cases, 
      daily_cases = df.cases$daily_cases,
      estimated_daily_cases = cbind(1, Bt) %*% t(t(lm$coefficients))
    )
    title <- ggtitle(paste("Knots: ", toString(p), ", Polynomial Degree: ", toString(d), ", AIC : ", toString(round(AIC(lm))), sep = ""))
    print(ggplot(data = df.fit, mapping = aes(x = cases, y = daily_cases)) + geom_point() + geom_line(aes(x = cases, estimated_daily_cases)) + title)
    print(qplot(1:nrow(knots), AIC) + title)
    print(qplot(1:nrow(knots), BIC) + title)
  }
}

## Heat map on 2D
p = 2
for (d in D){
  AIC <- c()
  BIC <- c()
  knots <- read_csv(paste("/Users/MisterGrinvalds/Repos/change-point/knots_october07_", toString(p), ".txt", sep = ''), col_names = FALSE)
  for (i in 1:nrow(knots)){
    knot <- knots %>% slice(i) %>% unlist()
    Bt <- form_truncated_basis(df.cases$cases, knots = knot, D = d)
    lm <- lm(daily_cases ~ Bt, data = df.cases)
    AIC <- c(AIC, AIC(lm))
    BIC <- c(BIC, BIC(lm))
  }
  df <- cbind(knots, AIC)
  title <- ggtitle(paste("Knots: ", toString(p), ", Polynomial Degree: ", toString(d), ", AIC : ", toString(round(AIC(lm))), sep = ""))
  print(ggplot(data = df, mapping = aes(x = X1, y = X2, z = AIC)) + geom_contour_filled(bins = 14) + title)
}

## Get best fitting models
AIC_min
BIC_min

################################################################################
# 10/23/20
################################################################################
library(tidyverse)

## read in data
file <- file.choose()
df <- read_csv(file)
df.cases <- df %>% 
  mutate(daily_cases = cases  - dplyr::lag(cases, order_by = date, default = first(cases)))

## update knot mapping to include an epsilon
map_knots <- function(file, knots, depth = 1, epsilon = 1, knot_vec = c(), ref = knots[[1]][1]-(epsilon+1)){
  for (knot in knots[[depth]]){
    if (knot > ref + epsilon && depth < length(knots)){
      map_knots(file, knots, depth = depth + 1, epsilon = epsilon, knot_vec = c(knot_vec, knot), ref = knot)
    }
    if (knot > ref + epsilon && depth == length(knots)){
      readr::write_lines(toString(c(knot_vec, knot)), file, append = TRUE)
    }
  }
}

## test knot mapping tool
P <- 1:6
epsilon = 1
for (p in P){
  #knots_list <- rep(list(round(seq(0, max(df.cases$cases), by = 100000))), p)
  knots_list <- rep(list(round(seq(0, max(df.cases$cases), length = 20))), p)
  filename <- paste("/Users/MisterGrinvalds/Repos/change-point/knots_october07_", toString(p), ".txt", sep = "")
  map_knots(file = filename, knots_list, epsilon = epsilon)
}

## run analysis again
form_truncated_basis <- function(X, knots, D = 1){
  X1 <- outer(X, 1:D, "^")
  X2 <- outer(X, knot,">") * outer(X, knot, "-")^D
  Bt <- cbind (X1, X2)
  return(Bt)
}

D <- 1:3
AIC_min <- c()
BIC_min <- c()
for (d in D) {
  for (p in P){
    AIC <- c()
    BIC <- c()
    knots <- read_csv(paste("/Users/MisterGrinvalds/Repos/change-point/knots_october07_", toString(p), ".txt", sep = ''), col_names = FALSE)
    for (i in 1:nrow(knots)){
      knot <- knots %>% slice(i) %>% unlist()
      Bt <- form_truncated_basis(df.cases$cases, knots = knot, D = d)
      lm <- lm(daily_cases ~ Bt, data = df.cases)
      AIC <- c(AIC, AIC(lm))
      BIC <- c(BIC, BIC(lm))
    }
    AIC_min <- rbind(AIC_min, c(p, d, which.min(AIC), min(AIC)))
    BIC_min <- rbind(BIC_min, c(p, d, which.min(BIC), min(BIC)))
    knot <- knots %>% slice(which.min(BIC)) %>% unlist()
    knot[is.na(knot)] <- 0
    Bt <- form_truncated_basis(df.cases$cases, knots = knot, D = d)
    lm <- lm(daily_cases ~ Bt, data = df.cases)
    lm$coefficients[is.na(lm$coefficients)] <- 0
    df.fit <- tibble(
      cases = df.cases$cases, 
      daily_cases = df.cases$daily_cases,
      estimated_daily_cases = cbind(1, Bt) %*% t(t(lm$coefficients))
    )
    title <- ggtitle(paste("Knots: ", toString(p), ", Polynomial Degree: ", toString(d), ", BIC : ", toString(round(AIC(lm))), sep = ""))
    print(
      ggplot(data = df.fit, mapping = aes(x = cases, y = daily_cases)) + 
        geom_point() +
        geom_line(aes(x = cases, estimated_daily_cases), colour = "red") + 
        geom_vline(xintercept = as.vector(knot), colour = "blue") +
        title
    )
  }
}


################################################################################
# 10/26/20
################################################################################

#############
# Load Data #
#############
library(tidyverse)

file <- file.choose()
df <- read_csv(file)
df.cases <- df %>% 
  mutate(daily_cases = cases  - dplyr::lag(cases, order_by = date, default = first(cases)))

##############
# Knot Files #
##############
map_knots <- function(
  file,
  knots,
  buffer_size = 10000,
  epsilon = 0
){
  max_depth <- length(knots)
  buffer <- matrix(nrow = buffer_size, ncol = max_depth)
  i <- 1
  
  fetch_knots <- function(depth, elements, ref){
    for (element in knots[[depth]]){
      if(depth < max_depth & element > ref + epsilon){
        fetch_knots(depth + 1, c(elements, element), element)
      } else if (element > ref + epsilon){
        buffer[i, ] <<- c(elements, element)
        i <<- i + 1
        if (i > buffer_size){
          i <<- 1
          write_csv(data.frame(buffer), file, col_names = FALSE, append = TRUE)
          buffer <<- matrix(nrow = buffer_size, ncol = max_depth)
        }
      }
    }
    if (depth == 1) return(buffer)
  }
  print(buffer)
  fetch_knots(1, c(), knots[[1]][1] - 1)
  write_csv(data.frame(buffer[1:i-1, ]), file, col_names = FALSE, append = TRUE)
}


generate_knots_files <- function(
  knot_range,
  ...,
  buffer_left = 0,
  buffer_right = 0,
  by_length = 20,
  by_resolution = FALSE,
  by_specification = FALSE,
  epsilon = 0,
  extension = ".txt",
  filename_stem = "knot_map_",
  K = 1:1,
  round_knots = TRUE
){
  for (k in K){
    if (by_specification){
      knot_locations <- knot_range
    } else if (by_resolution){
      knot_locations <- seq(min(knot_range) + buffer_left, max(knot_range) - buffer_right, by = by_resolution)
    } else{
      knot_locations <- seq(min(knot_range) + buffer_left, max(knot_range) - buffer_right, length = by_length)
    }
    if (round_knots) knot_locations <- round(knot_locations)
    knot_list <- rep(list(knot_locations), k)
    filename <- paste(filename_stem, toString(k), extension, sep = '')
    print(knot_list)
    print(filename)
    map_knots(file = filename, knots = knot_list, epsilon = epsilon, ...)
  }
}


generate_knots_files(
  knot_range = unique(df.cases$cases),
  buffer_left = 0,
  buffer_right = 0,
  by_specification = TRUE,
  epsilon = 0,
  #filename_stem = "/Users/MisterGrinvalds/Repos/change-point/knots_october07_",
  filename_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_knots_",
  K = 1:3
)

###################
# Generate Models #
###################
form_basis <- function(X, knots, D = 1, derivative = FALSE){
  ifelse(derivative, powers <- 0:(D-1), powers <-1:D)
  B <- outer(X, powers, "^")
  for (j in 1:length(knots)){
    if (derivative){
      #for (d in powers){
      #  B <- cbind(B, outer(X, knots[j], ">") * X^d)
      #}
      if (powers == 0){
        B <- cbind(B, outer(X, knots[j], ">") * 1)
      } else{
        B <- cbind(B, outer(X, knots[j], ">") * 1, outer(X_knots[,j], 1:(D-1), "^"))
      }
    } else{
      X_knots <- outer(X, knots, ">") * outer(X, knots, "-")
      B <- cbind(B, outer(X_knots[,j], powers, "^"))
    }
  }
  return(B)
}
  

GCV <- function(lm, d, knots){
  e <- residuals(lm)
  n <- length(e)
  return(sum(e)^2/((n - d*(length(knots)+1) - 1)^2 / n))
} 


get_results <- function(
  Y,
  basis,
  d,
  knots,
  epsilon = NA
){
  lm <- lm(Y ~ basis)
  if (NA %in% lm$coefficients){
    return(FALSE)
  } else{
    AIC <- AIC(lm)
    BIC <- BIC(lm)
    GCV <- GCV(lm, d = d, knots = knots)
    SSE <- sum(residuals(lm)^2)
    return(c(AIC, BIC, GCV, SSE, d, length(knots), knots))
  }
}


run_modeling <- function(
  X,
  Y,
  ...,
  buffer_size = 10000,
  D = 1:1,
  K = 1:1,
  extension = ".txt",
  filein_stem = "knot_map_",
  fileout_stem = "model_results"
){
  for (k in K){
    buffer <- matrix(nrow = buffer_size, ncol = index_dict[['num_results']] + k)
    i <- 1

    fileout <- paste(fileout_stem, toString(k), extension, sep ='')
    if (!file.exists(fileout)) file.create(fileout)
    connection_in = file(paste(filein_stem, toString(k), extension, sep =''), "r")
    
    while (TRUE){
      line = readLines(connection_in, n = 1)
      if (length(line) == 0) break
      knot <- as.numeric(unlist(strsplit(line, split =",")))
      for (d in D){
        results <- c(FALSE)
        B <- form_basis(X = X, knots = knot, D = d)
        results <- get_results(Y = Y, basis = B, d = d, knots = knot, ...)
        if (!(FALSE %in% results)){
          buffer[i, ] <- results
          i <- i + 1
        }
        if (i > buffer_size){
          write_csv(data.frame(buffer), fileout, col_names = FALSE, append = TRUE)
          buffer <- matrix(nrow = buffer_size, ncol = index_dict[['num_results']] + k)
          i <- 1
        }
      }
    }

    write_csv(data.frame(buffer[1:i-1, ]), fileout, col_names = FALSE, append = TRUE)
    close(connection_in)
  }
}


run_modeling(
  X = df.cases$cases, 
  Y = df.cases$daily_cases, 
  D = 1:3, 
  #K = 1:8, 
  K = 6:10, 
  #filein_stem = "/Users/MisterGrinvalds/Repos/change-point/knots_october07_", 
  #fileout_stem = "/Users/MisterGrinvalds/Repos/change-point/results_october07_"
  filein_stem = "/Users/MisterGrinvalds/Repos/change-point/knots_june18_", 
  fileout_stem = "/Users/MisterGrinvalds/Repos/change-point/results_june18_"
)


############################
# Collect Candidate Models #
############################
index_dict <- list(
  "AIC" = 1,
  "BIC" = 2,
  "GCV" = 3,
  "SSE" = 4,
  "d" = 5,
  "k" = 6,
  "epsilon" = 7,
  "knot_start" = 8,
  "num_criteria" = 4,
  "num_results" = 6
)


get_candidate_models <- function(
  filein_stem,
  extension = ".txt",
  candidates = 1,
  D = 1:1,
  K = 1:1,
  metric = "AIC"
){
  metric = index_dict[[metric]]
  topmodels = rep(list(rep(FALSE, index_dict[["num_criteria"]])), candidates)
  knot_vec = rep(list(rep(FALSE), max(K)), candidates)
  for (k in K){
    connection_in = file(paste(filein_stem, k, extension, sep = ''), "r")
    while (TRUE){
      line = readLines(connection_in, n = 1)
      if (length(line) == 0) break
      line <- as.numeric(unlist(strsplit(line, split =",")))
      if (line[index_dict[["d"]]] %in% D){
        for (i in 1:candidates){
          if (line[metric] < topmodels[[i]][metric] | topmodels[[i]][metric] == FALSE){
            topmodels <- append(topmodels, list(line[1:(index_dict[["knot_start"]]-1)]), after = i-1)
            knot <- rep(NA, max(K))
            knot[1:k] <- line[index_dict[["knot_start"]]:length(line)]
            knots <- append(knots, list(knot), after = i-1)
            break
          }
        }
        topmodels <- topmodels[1:candidates]
        knots <- knots[1:candidates]
      }
    }
    close(connection_in) 
  }
  topmodels <- data.frame(matrix(unlist(topmodels), nrow=length(topmodels), byrow=T))
  colnames(topmodels) = names(index_dict)[1:index_dict[["k"]]]
  knots <- data.frame(knot = matrix(unlist(knots), nrow=length(knots), byrow=T))
  return(list(topmodels = topmodels, knots = knots, metric = metric))
}


get_derivative <- function(X, lm, knots, D){
  coef <- coefficients(lm)[-1] * 1:D
  B <- form_basis(X, knots, D = D, derivative = TRUE)
  return(B %*% t(t(coef)))
}


get_candidate_models(
  filein_stem = "/Users/MisterGrinvalds/Repos/change-point/results_june18_", 
  candidates = 3, 
  D = 1, 
  K = 1:8, 
  metric = "GCV"
)

#######
# Viz #
#######
visualize_candidate_models <- function(
  X,
  Y,
  topmodels_df,
  knot_df,
  metric,
  X_lab =   "total_cases",
  Y_lab = "daily_cases",
  fitted_line = TRUE,
  fitted_line_colour = "red",
  knots_lines = TRUE,
  knots_lines_colour = "blue",
  derivative_line = TRUE,
  derivative_line_colour = "darkgreen"
){
  for (i in 1:nrow(topmodels_df)){
    knots <- knot_df[i,]
    knots <- knots[!is.na(knots)]
    basis <- form_basis(X, knots, topmodels_df[i, index_dict[["d"]]])
    lm <- lm(Y ~ basis)
    lm$coefficients[is.na(lm$coefficients)] <- 0
    estimate <- paste(Y_lab, "_hat", sep = '')
    df.fit <- tibble(
      X_lab = X, 
      Y_lab = Y,
      estimate = cbind(1, basis) %*% t(t(lm$coefficients)),
      derivative = get_derivative(X, lm, knots, topmodels_df[i, index_dict[["d"]]])
    )
    title <- ggtitle(paste(
      "Knots: ", toString(topmodels_df[i, index_dict[["k"]]]), 
      ", Degrees: ", toString(topmodels_df[i, index_dict[["d"]]]), 
      ", Epsilon: ", toString(topmodels_df[i, index_dict[["epsilon"]]]), 
      ", ", metric, ": ", toString(round(topmodels_df[i, index_dict[[metric]]])), 
      sep = ""))
    plot <- ggplot(data = df.fit, mapping = aes(x = X, y = Y)) + geom_point() + title + xlab(X_lab) + ylab(Y_lab)
    if (fitted_line) plot <- plot + geom_line(aes(x = X, estimate), colour = fitted_line_colour)
    if (knots_lines) plot <- plot + geom_vline(xintercept = knots, colour = knots_lines_colour)
    if (derivative_line) plot <- plot + geom_line(aes(x = X, y = derivative), colour = derivative_line_colour)
    print(ggplot(mapping = aes_string(x=X, y=df.fit$derivative))+geom_line())
    print(plot)
  }
}


candidates <- get_candidate_models(
  filein_stem = "/Users/MisterGrinvalds/Repos/change-point/results_june18_", 
  candidates = 1, 
  D = 1:3, 
  K = 1:8, 
  metric = "BIC"
)

visualize_candidate_models(
  df.cases$cases,
  "total_cases",
  df.cases$daily_cases,
  "daily_cases",
  candidates$topmodels,
  candidates$knots,
  names(index_dict)[[candidates$metric]]
)
  
  
################
# Benchmarking #
################
library(lineprof)

n<- 170
k <- 4
Knots <- rep(list(1:n), k)
epsilon <- 0

profile.1 <- lineprof(generate_knots_files(
  knot_range = 1:n,
  buffer_left = 0,
  buffer_right = 0,
  by_specification = TRUE,
  epsilon = 0,
  filename_stem = "/Users/MisterGrinvalds/Repos/change-point/profile1",
  K = k
))

profile.2 <- lineprof(
  for (knot_i in Knots[[1]]){
    for (knot_j in Knots[[2]]){
      for (knot_k in Knots[[3]]){
        if (knot_j > knot_i + epsilon){
          if (knot_k > knot_j + epsilon){
            write_lines(toString(c(knot_i, knot_j, knot_k)), "/Users/MisterGrinvalds/Repos/change-point/profile2.txt", append = TRUE)
          }
        }
      } 
    }
  }
)

index <- 1
buffer = 100000
M <- list()
m <- matrix(nrow = buffer, ncol = 3)
profile.3 <- lineprof(
  {for (knot_i in Knots[[1]]){
    for (knot_j in Knots[[2]]){
      for (knot_k in Knots[[3]]){
        for (knot_l in Knots[[4]]){
          if (knot_j > knot_i + epsilon){
            if (knot_k > knot_j + epsilon){
              if (knot_l > knot_k + epsilon){
                m[index,] <- c(knot_i, knot_j, knot_k)
                index <- index + 1
                if (index > buffer){
                  M <- append(M, list(data.frame(m)))
                  m <- matrix(nrow = buffer, ncol = 3)
                  index <- 1
                }
              }
            }
          }
        }
      }
    }
  }
  M <- append(M, list(data.frame(m)))
  for (df in M){write_csv(df, "/Users/MisterGrinvalds/Repos/change-point/profile3.txt", append = TRUE)}}
)

file <- "/Users/MisterGrinvalds/Repos/change-point/Rcpp.txt"
Rcpp::cppFunction(
  code = '
  void fetch_knots(Rcpp::List& Knots, int depth, int& max_depth, std::vector<int> knot_vec, int ref, int& epsilon){
    List knots = Knots[depth];
    for (int i = 0; i < knots.size(); ++i){
      int knot = knots[i];
      if (depth < max_depth && knot > ref + epsilon){
        std::vector<int> knot_vec_new = knot_vec;
        knot_vec_new.push_back(knot);
        fetch_knots(Knots, depth + 1, max_depth, knot_vec_new, knot, epsilon);
      } else if (depth == max_depth && knot > ref + epsilon){
        std::vector<int> knot_vec_keeper = knot_vec;
        knot_vec_keeper.push_back(knot);
      }
    }
  }
  '
)
profile.4 <- lineprof(fetch_knots(Knots, 0, length(Knots)-1, 0, 0, epsilon))

profile.5 <- lineprof(
  for (knot_i in Knots[[1]]){
    for (knot_j in Knots[[2]]){
      for (knot_k in Knots[[3]]){
        for (knot_l in Knots[[4]]){
          if (knot_j > knot_i + epsilon){
            if (knot_k > knot_j + epsilon){
              if (knot_l > knot_k + epsilon){
                c(knot_i, knot_j, knot_k, knot_l)
              }
            }
          }
        }
      }
    }
  }
)

profile.2.4 <- lineprof(
  for (knot_i in Knots[[1]]){
    for (knot_j in Knots[[2]]){
      for (knot_k in Knots[[3]]){
        for (knot_l in Knots[[3]]){
          if (knot_j > knot_i + epsilon){
            if (knot_k > knot_j + epsilon){
              if (knot_l > knot_k + epsilon){
                write_lines(toString(c(knot_i, knot_j, knot_k, knot_l)), "/Users/MisterGrinvalds/Repos/change-point/profile2.txt", append = TRUE)
                }
              }
            }
          }
        }
      }
    }
)

##############
# 10-28-2020 #
##############


generate_knots_files(
  knot_range = unique(df.cases$cases),
  buffer_left = 0,
  buffer_right = 0,
  by_specification = TRUE,
  epsilon = 250000,
  filename_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_250k_knots_",
  K = 1:8
)

run_modeling(
  X = df.cases$cases, 
  Y = df.cases$daily_cases, 
  D = 1:3, 
  K = 1:8, 
  filein_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_250k_knots_", 
  fileout_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_250k_results_"
)

candidates <- get_candidate_models(
  filein_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_250k_results_", 
  candidates = 1, 
  D = 1, 
  K = 4, 
  metric = "AIC"
)

visualize_candidate_models(
  X = df.cases$cases, 
  Y = df.cases$daily_cases,
  candidates$topmodels,
  candidates$knots,
  names(index_dict)[[candidates$metric]]
)

generate_knots_files(
  knot_range = unique(df.cases$cases),
  buffer_left = 0,
  buffer_right = 0,
  by_specification = TRUE,
  epsilon = 100000,
  filename_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_100k_knots_",
  K = 1:5
)

run_modeling(
  X = df.cases$cases, 
  Y = df.cases$daily_cases, 
  D = 1:3, 
  K = 1:5, 
  filein_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_100k_knots_", 
  fileout_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_100k_results_"
)

generate_knots_files(
  knot_range = unique(df.cases$cases),
  buffer_size = 1000000,
  buffer_left = 0,
  buffer_right = 0,
  by_specification = TRUE,
  epsilon = 0,
  filename_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_0_knots_",
  K = 1:5
)

run_modeling(
  X = df.cases$cases, 
  Y = df.cases$daily_cases,
  buffer_size = 1000000,
  D = 1:3, 
  K = 1:5, 
  filein_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_0_knots_", 
  fileout_stem = "/Users/MisterGrinvalds/Repos/change-point/june18_eps_0_results_"
)

#################
# July 18, 2020 #
#################
batch_dict <- list(
  "june18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/june18_data/june18_",
    "D" = 1:2,
    "K" = 1:8
  ),
  "july18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/july18_data/july18_",
    "D" = 1:2,
    "K" = 1:8
  ),
  "august18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/august18_data/august18_",
    "D" = 1:2,
    "K" = 1:8
  ),
  "september18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/september18_data/september18_",
    "D" = 1:2,
    "K" = 1:8
  ),
  "october18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/october18_data/october18_",
    "D" = 1:2,
    "K" = 1:8
  )
)

library(tidyverse)

for (batch in batch_dict){
  file_nytdata <- paste(batch[["filename_stem"]], "us.csv", sep="")
  file_knots <- paste(batch[["filename_stem"]], "eps_", toString(as.integer(batch_dict[["june18"]][["epsilon"]])), "_knots_", sep="")
  file_results <- paste(batch[["filename_stem"]], "eps_", toString(as.integer(batch_dict[["june18"]][["epsilon"]])), "_results_", sep="")
  
  df <- read_csv(file_nytdata)
  df.cases <- df %>% 
    mutate(daily_cases = cases  - dplyr::lag(cases, order_by = date, default = first(cases)))
  
  generate_knots_files(
    knot_range = unique(df.cases$cases),
    by_specification = TRUE,
    epsilon = batch[["epsilon"]],
    filename_stem = file_knots,
    K = batch[["K"]]
  )
  
  run_modeling(
    X = df.cases$cases, 
    Y = df.cases$daily_cases, 
    D = batch[["D"]],
    K = batch[["K"]], 
    epsilon = batch[["epsilon"]],
    filein_stem = file_knots, 
    fileout_stem = file_results
  )
}
