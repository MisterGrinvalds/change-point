############
# Packages #
############
require(DBI)
require(dbplyr)
require(RSQLite)
require(tidyverse)

#############
# Variables #
#############

cp_dict <- list(
  "AIC" = 1,
  "BIC" = 2,
  "GCV" = 3,
  "SSE" = 4,
  "d" = 5,
  "k" = 6,
  #"epsilon" = 7,
  "knot_start" = 7,
  "num_criteria" = 4,
  "num_results" = 6
)

#############
# Functions #
#############

createKnots <- function(
  db,
  knots,
  buffer_size = 10000,
  epsilon = 0
){
  max_depth <- length(knots)
  buffer <- matrix(nrow = buffer_size, ncol = max_depth)
  i <- 1
  

  
  #      
  fetch_knots <- function(depth, elements, ref){
    for (element in knots[[depth]]){
      if(depth < max_depth & element > ref + epsilon){
        fetch_knots(depth + 1, c(elements, element), element)
      } else if (element > ref + epsilon){
        buffer[i, ] <<- c(elements, element)
        i <<- i + 1
        if (i > buffer_size){
          i <<- 1
          df <- data.frame(buffer)
          colnames(df) <- db[["columns"]]
          dbAppendTable(db[["db"]], db[["table"]], df)
          buffer <<- matrix(nrow = buffer_size, ncol = max_depth)
        }
      }
    }
  }
  
  fetch_knots(1, c(), knots[[1]][1] - epsilon)
  df <- data.frame(buffer[1:i-1, ])
  colnames(df) <- db[["columns"]]
  dbAppendTable(db[["db"]], db[["table"]], df)
}


createKnotsTables <- function(
  database,
  knot_range,
  ...,
  by_length = 20,
  by_resolution = FALSE,
  by_specification = FALSE,
  epsilon = 0,
  K = 1:1,
  round_knots = TRUE
){
  db <- dbConnect(RSQLite::SQLite(), database)
  
  for (k in K){
    print(paste("Generating combinations for", toString(k), "possible knots."))
    
    if (by_specification){
      knot_locations <- knot_range
    } else if (by_resolution){
      knot_locations <- seq(min(knot_range), max(knot_range), by = by_resolution)
    } else{
      knot_locations <- seq(min(knot_range), max(knot_range), length = by_length)
    }
    if (round_knots) knot_locations <- round(knot_locations)
    knot_list <- rep(list(knot_locations), k)
    
    columns <- c()
    for (i in 1:k) columns <- c(columns, paste("knot_", toString(i), sep =''))
    df <- data.frame(rep(list(0.0), k))
    colnames(df) <- columns
    tablename <- paste(toString(k), "_Knots", sep ='')
    dbCreateTable(db, tablename, df)
    
    createKnots(
      db = list(db = db, table = tablename, columns = columns), 
      knots = knot_list, 
      epsilon = epsilon, 
      ...
    )
  }
  
  dbDisconnect(db)
}


createModelsTables <- function(
  database,
  X,
  Y,
  ...,
  buffer_size = 10000,
  D = 1:1,
  K = 1:1,
  extension = ".txt"
){
  db <- dbConnect(RSQLite::SQLite(), database)
  
  for (k in K){
    print(paste("Generating models on", toString(k), "total knots."))
    
    columns <- names(cp_dict)[1:cp_dict[["num_results"]]]
    for (i in 1:k) columns <- c(columns, paste("knot_", toString(i), sep =''))
    table_out <- paste(toString(k), "_Results", sep ='') # change to _Models?
    buffer <- matrix(ncol = cp_dict[['num_results']] + k)
    df <- data.frame(buffer)
    colnames(df) <- columns
    dbCreateTable(db, table_out, df)
    
    table_in <- paste(toString(k), "_Knots", sep ='')
    query <- paste('SELECT COUNT(*) FROM "', table_in, '";', sep ='')
    total_knots <- dbGetQuery(db, query)
    
    rolling_count <- 0
    while(TRUE){
      buffer <- matrix(nrow = buffer_size*max(D), ncol = cp_dict[['num_results']] + k)
      buffer_i <- 1
      query <- paste('SELECT * FROM "', table_in, '" LIMIT ', toString(rolling_count), ', ', toString(buffer_size), ';', sep = '')
      knots <- dbGetQuery(db, query)
      if(nrow(knots) == 0) break
      for (i in 1:nrow(knots)){
        for (d in D){
          B <- formBasis(X = X, knots = as.numeric(knots[i,]), D = d)
          results <- c(getModelResults(Y = Y, basis = B, d = d, knots = as.numeric(knots[i,]), epsilon = epsilon))
          if (!(FALSE %in% results)){
            buffer[buffer_i, ] <- results
            buffer_i <- buffer_i + 1
          }
        }
      }
      df <- data.frame(buffer[1:buffer_i, ])
      colnames(df) <- columns
      dbAppendTable(db, table_out, df)
      buffer <- matrix(nrow = buffer_size*max(D), ncol = cp_dict[['num_results']] + k)
      rolling_count <- rolling_count + buffer_size  
    }
  }
  
  dbDisconnect(db)
}


formBasis <- function(X, knots, D = 1, derivative = FALSE){
  if (derivative){
    if(D ==1){
      B <- cbind(1, outer(X, knots, ">") *1)
    } else{
      powers <- 0:(D-1)
      B <- outer(X, powers, "^")
      for (j in 1:length(knots)){
        X_knots <- c(outer(X, knots[j], ">") * outer(X, knots[j], "-"))
        B <- cbind(B, outer(X, knots[j], ">") * 1, outer(X_knots, 1:(D-1), "^"))
      }
      
    }
  }
  else{
    powers <-1:D
    B <- outer(X, powers, "^")
    for (j in 1:length(knots)){
      X_knots <- c(outer(X, knots[j], ">") * outer(X, knots[j], "-"))
      B <- cbind(B, outer(X_knots, powers, "^"))
    }
  }
  return(B)
}


GCV <- function(lm, d, knots){
  e <- residuals(lm)
  n <- length(e)
  return(sum(e^2)*n/(n - d*(length(knots)+1) - 1)^2)
} 


getModelResults <- function(
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


getCandidateModels <- function(
  database,
  candidates = 1,
  D = 1:1,
  K = 1:1,
  metric = "AIC",
  ascending = TRUE,
  knot_example = FALSE
){
  db <- dbConnect(RSQLite::SQLite(), database)
  
  topmodels <- list()
  knots <- list()
  for (k in K){
    table_in <- paste(k, "_Results", sep = "")
    ifelse (ascending, sort <- "ASC", sort <- "DESC")
    query <- paste(
      'SELECT * FROM "', table_in, 
      '" WHERE ', metric, ' IS NOT NULL AND d BETWEEN ', min(D), ' AND ', max(D), 
      ' ORDER BY ', metric, ' ', sort, 
      ' LIMIT ', candidates, 
      sep = ""
    )
    if (knot_example){
      knot_sort <- c()
      if (k-1 > 0) for (i in 1:(k-1)) knot_sort <- paste(knot_sort, '"knot_', toString(i), '", ', sep ="")
      knot_sort <- paste(knot_sort, '"knot_', toString(k), '"', sep ="")
      query <- paste( 
        'SELECT * FROM "', table_in, 
        '" WHERE d BETWEEN ', min(D), ' AND ', max(D), 
        ' ORDER BY ', knot_sort, ' ', sort, 
        ' LIMIT ', candidates, 
        sep = ""
      )
    }
    results <- dbGetQuery(db, query)
    topmodels <- append(topmodels, list(results[ , 1:cp_dict[["num_results"]]]))
    knots <- append(knots, list(as.data.frame(results[ , cp_dict[["knot_start"]]:dim(results)[2]])))
  }
    
  dbDisconnect(db)
  return(list(topmodels = topmodels, knots = knots, metric = metric))
}


getModelDerivative <- function(X, lm, knots, D){
  coef <- coefficients(lm)[-1] * 1:D
  B <- formBasis(X, knots, D = D, derivative = TRUE)
  return(B %*% t(t(coef)))
}


visualizeCandidateModels <- function(
  X,
  Y,
  candidates,
  output_dict,
  dpi = 300,
  features = list(),
  labels = list()
){
  if (!("derivative_line" %in% names(features)))  features[["derivative_line"]] = list(draw = TRUE, colour = "#FFB000")
  if (!("fitted_line" %in% names(features)))      features[["fitted_line"]] = list(draw = TRUE, colour = "#DC267F")
  if (!("knots_lines" %in% names(features)))      features[["knots_lines"]] = list(draw = TRUE, colour = "#648FFF")
  
  if (!("data_thru" %in% names(labels)))  labels[["data_thru"]] = "no date specified"
  if (!("X_lab" %in% names(labels)))      labels[["X_lab"]] = "total_cases"
  if (!("Y_lab" %in% names(labels)))      labels[["Y_lab"]] = "daily_cases"
  if (!("subtitle" %in% names(labels)))   labels[["subtitle"]] = TRUE
  if (!("title" %in% names(labels)))      labels[["title"]] = TRUE
  
  for (i in 1:length(candidates$topmodels)){
    for (j in 1:nrow(candidates$topmodels[[i]])){
      model <- as.numeric(candidates$topmodels[[i]][j,])
      knots <- as.numeric(candidates$knots[[i]][j,])
      knots <- knots[!is.na(knots)]
      
      metric <- candidates$metric
      subtitle <- c()
      ifelse(
        labels[["subtitle"]] != TRUE, 
        subtitle <- labels[["subtitle"]],
        ifelse(
          labels[["subtitle"]] == TRUE,
          subtitle <- paste("Data Through", labels[["data_thru"]]),
          subtitle <- labels[["subtitle"]]
        )
      )
      title <- c()
      ifelse(
        labels[["title"]] != TRUE, 
        title <- labels[["title"]], 
        ifelse(
          labels[["title"]] == TRUE, 
          title <- paste(
            "Knots: ", toString(model[cp_dict[["k"]]]), 
            ", Degrees: ", toString(model[cp_dict[["d"]]]), 
            #", Epsilon: ", toString(model[cp_dict[["epsilon"]]]), 
            ", ", metric, ": ", toString(round(model[cp_dict[[metric]]])), 
            sep = ""),
          title <- labels[["title"]]
        )
      )
      X_lab <- xlab(labels[["X_lab"]])
      Y_lab <- ylab(labels[["Y_lab"]])
      
      basis <- formBasis(X, knots, model[cp_dict[["d"]]])
      lm <- lm(Y ~ basis)
      lm$coefficients[is.na(lm$coefficients)] <- 0
      estimate <- paste(labels[["Y_lab"]], "_hat", sep = '')
      df.fit <- tibble(
        X_lab = X, 
        Y_lab = Y,
        estimate = cbind(1, basis) %*% t(t(lm$coefficients)),
        derivative = getModelDerivative(X, lm, knots, model[cp_dict[["d"]]])
      )
      
      plot <- ggplot(data = df.fit, mapping = aes(x = X)) + 
        geom_point(aes(y = Y)) + 
        xlab(X_lab) + 
        ylab(Y_lab)
      ifelse(title == FALSE, 0, plot <- plot + labs(title = title))
      ifelse(subtitle == FALSE, 0, plot <- plot + labs(subtitle = subtitle))
      if (features[["fitted_line"]][["draw"]]){
        plot <- plot + 
          geom_line(
            aes(y = estimate), 
            colour = features[["fitted_line"]][["colour"]]
          )
      }
      if (features[["derivative_line"]][["draw"]]){
        plot <- plot + 
          geom_line(
            aes(y = derivative*round((max(X) - min(X))/length(X))+round(max(Y))/2), 
            colour = features[["derivative_line"]][["colour"]]
          ) +
          scale_y_continuous(
            sec.axis = sec_axis(
              ~.-round(max(Y))/2, 
              name = "derivative"
            )
          )
      }
      if (features[["knots_lines"]][["draw"]]){
        plot <- plot + 
          geom_vline(
            xintercept = knots, 
            colour = features[["knots_lines"]][["colour"]]
          )
      }
      
      filename <- paste(output_dict[["filename_stem"]], "_", toString(i), "_", toString(j), ".png", sep = "")
      ggsave(
        filename = filename,
        plot = plot,
        width = output_dict[["width"]],
        height = output_dict[["height"]],
        dpi = dpi
      )
    }
  }
}

