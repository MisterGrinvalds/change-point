w <- 10
h <- 8
delta <- 200/170*w-w

batch_dict = list(
  "june18" = list(
    "data_thru" = "June 18",
    "epsilon" = 250000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/june18_data/june18_",
    "K" = 1:7,
    "width" = w,
    "height" = h
  ),
  "july18" = list(
    "data_thru" = "July 18",
    "epsilon" = 250000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/july18_data/july18_",
    "K" = 1:8,
    "width" = w + delta,
    "height" = h
  ),
  "august18" = list(
    "data_thru" = "August 18",
    "epsilon" = 250000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/august18_data/august18_",
    "K" = 1:8,
    "width" = w + 2*delta,
    "height" = h
  ),
  "september18" = list(
    "data_thru" = "September 18",
    "epsilon" = 250000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/september18_data/september18_",
    "K" = 1:8,
    "width" = w + 3*delta,
    "height" = h
  ),
  "october18" = list(
    "data_thru" = "October 18",
    "epsilon" = 250000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/october18_data/october18_",
    "K" = 1:8,
    "width" = w + 4*delta,
    "height" = h
  )
)

for(batch in batch_dict){
  file_nytdata <- paste(batch[["filename_stem"]], "us.csv", sep="")
  database_name <- paste(batch[["filename_stem"]], ".mysqlite", sep="")
  
  df <- read_csv(file_nytdata)
  df.cases <- df %>% 
    mutate(daily_cases = cases  - dplyr::lag(cases, order_by = date, default = first(cases)))
  
  candidates <- getCandidateModels(
    database = database_name,
    candidates = 1,
    D = 2,
    K = batch[["K"]],
    metric = "BIC"
  )

  visualizeCandidateModels(
    df.cases$cases,
    df.cases$daily_cases,
    candidates,
    labels = list(  
      "data_thru" = batch[["data_thru"]]
    ),
    output_dict =  list(
      "filename_stem" = batch[["filename_stem"]],
      "width" = batch[["width"]],
      "height" = batch[["height"]]
    )
  )
  
}

# free form for generating graphics
candidates <- getCandidateModels(
  database = database_name,
  candidates = 1,
  D = 2,
  K = 1,
  metric = "SSE"
)

visualizeCandidateModels(
  df.cases$cases,
  df.cases$daily_cases,
  candidates,
  labels = list("data_thru" = "June 18th", "title" = "US Covid-19 Data", "subtitle" = FALSE),
  features = list("derivative_line" = list(draw = FALSE, colour = "#FFB000")),
  output_dict =  list(
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/june18_data/june18_",
    "width" = 5,
    "height" = 4
  )
)

plotBasis <- function(X, B){
  plot <- ggplot(mapping = aes(x = X))
  #palette <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
  palette <- c("#648FFF", "#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")
  lines <- c()
  for (j in 1:dim(B)[2]){
    assign(paste("col_", toString(j), sep = ""), B[,j])
    plot <- plot + geom_line(aes_string(y = paste("col_", toString(j), sep = "")), colour = palette[j])
  }
  plot <- plot + xlab("X") + ylab("basis_value")
  print(plot)
}

X <- seq(0, 2, length.out = 100)
B <- formBasis(X, 1, D = 1)
plotBasis(X, cbind(1,B[,1]))
B <- formBasis(X, 1, D = 2)
plotBasis(X, cbind(1,B[,1:2]))
X <- seq(0, 2, length.out = 100)
B <- formBasis(X, 1, D = 1)
plotBasis(X, cbind(1,B))
B <- formBasis(X, 1, D = 2)
plotBasis(X, cbind(1,B))
B <- formBasis(X, c(0.66, 1.33), D = 1)
plotBasis(X, cbind(1,B))
B <- formBasis(X, c(0.66, 1.33), D = 2)
plotBasis(X, cbind(1,B))
B <- formBasis(X, c(0.5, 1, 1.5), D = 1)
plotBasis(X, cbind(1,B))
B <- formBasis(X, c(0.5, 1, 1.5), D = 2)
plotBasis(X, cbind(1,B))
B <- formBasis(X, c(0.4, 0.8, 1.2, 1.6), D = 1)
plotBasis(X, cbind(1,B))
B <- formBasis(X, c(0.4, 0.8, 1.2, 1.6), D = 2)
plotBasis(X, cbind(1,B))


#for fixed n:
n_ <- c()
X <- 1:20

k_<- c()

X <- 1:100
count <- 0
for(i in X){
  count <- count + 1
}
print(count)

X<-1:4
count <- 0
for(i in X){
  for (j in X){
    if(j > i) count <- count + 1
  }
}
print(count)
#n_ <- c(n_, count)
k_ <- c(k_, count)


#k_<- c()
#X <- 1:20
X<-1:15
count <- 0
for(i in X){
  for (j in X){
    for (k in X){
      if(j > i & k > j) count <- count + 1
    }
  }
}
print(count)
#n_ <- c(n_, count)
k_ <- c(k_, count)
X<-1:15


count <- 0
for(i in X){
  for (j in X){
    for (k in X){
      for (l in X){
      if(j > i & k > j & l > k) count <- count + 1
      }
    }
  }
}
print(count)
n_ <- c(n_, count)

X <- 1:10
count <- 0
for(i in X){
  for (j in X){
    for (k in X){
      for (l in X){
        for (m in X){
          if(j > i & k > j & l > k & m > l) count <- count + 1
        }
      }
    }
  }
}
print(count)
n_ <- c(n_, count)

count <- 0
for(i in X){
  for (j in X){
    for (k in X){
      for (l in X){
        for (m in X){
          for (n in X){
          if(j > i & k > j & l > k & m > l & n > m) count <- count + 1
          }
        }
      }
    }
  }
}
print(count)


count <- 0
for(i in X){
  for (j in X){
    for (k in X){
      for (l in X){
        for (m in X){
          for (n in X){
            for (o in X){
              for (p in X){
                for (q in X){
                  if(j > i & k > j & l > k & m > l & n > m & o > n & p > o & q > p) count <- count + 1
                }
              }
            }
          }
        }
      }
    }
  }
}
print(count)

# store
n_10 <- c(45, 120, 210, 252)
n_15 <- c(105, 455, 1365, 3003)
n_20 <- c(190, 1140, 4845, 15504)

plot <- ggplot(data = NULL, aes(x = 2:5))
plot + geom_line(aes(y=n_10)) + geom_line(aes(y=n_15)) + geom_line(aes(y=n_20))  

# k
k_1 <- c(10, 15, 20)
k_2 <- c(45, 105, 190)
k_3 <- c(120, 455, 1140)

plot <- ggplot(data = NULL, aes(x = c(10, 15, 20)))
plot + geom_line(aes(y=k_1)) + geom_line(aes(y=k_2)) + geom_line(aes(y=k_3)) 
n <- 231
d_k0 <- rep(1, n)
fn <- function(vec){
  tmp <- c()
  for (i in 1:length(vec)){
    tmp <- c(tmp, sum(vec[1:i]))
  }
  return(tmp)
}
d_k1 <- fn(d_k0)
d_k2 <- fn(d_k1)
d_k3 <- fn(d_k2)
d_k4 <- fn(d_k3)
d_k5 <- fn(d_k4)
print(cbind(d_k1, d_k2, d_k3, d_k4, d_k5))
