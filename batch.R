batch_dict <- list(
  "june18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/june18_data/june18_",
    "D" = 1:2,
    "K" = 1:6
  ),
  "july18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/july18_data/july18_",
    "D" = 1:2,
    "K" = 1:6
  ),
  "august18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/august18_data/august18_",
    "D" = 1:2,
    "K" = 1:6
  ),
  "september18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/september18_data/september18_",
    "D" = 1:2,
    "K" = 1:6
  ),
  "october18" = list(
    "epsilon" = 200000,
    "filename_stem" = "/Users/MisterGrinvalds/Repos/change-point/data/october18_data/october18_",
    "D" = 1:2,
    "K" = 1:6
  )
)

for (batch in batch_dict){
  file_nytdata <- paste(batch[["filename_stem"]], "us.csv", sep="")
  database_name <- paste(batch[["filename_stem"]], ".mysqlite", sep="")
  
  df <- read_csv(file_nytdata)
  df.cases <- df %>% 
    mutate(daily_cases = cases  - dplyr::lag(cases, order_by = date, default = first(cases)))
  
  createKnotsTables(
    database = database_name,
    knot_range = unique(df.cases$cases),
    by_specification = TRUE,
    epsilon = batch[["epsilon"]],
    K = batch[["K"]]
  )
  
  createModelsTables(
    database = database_name,
    X = df.cases$cases, 
    Y = df.cases$daily_cases, 
    D = batch[["D"]],
    K = batch[["K"]], 
    epsilon = batch[["epsilon"]],
  )
}
