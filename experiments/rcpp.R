n <- 30
eps <- 0


map_knots <- function(
  #file,
  knots,
  buffer_size = 10000,
  epsilon = 0
){
  max_depth <- length(knots)
  out_buffer <- list()
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
          out_buffer <<- append(out_buffer, list(data.frame(buffer)))
          buffer <<- matrix(nrow = buffer_size, ncol = max_depth)
        }
      }
    }
  }

  fetch_knots(1, c(), knots[[1]][1] - 1)
  out_buffer <- append(out_buffer, list(data.frame(buffer)[1:i-1, ]))
  return(out_buffer)
}


profile.5 <- lineprof(map_knots(Knots))

Rcpp::cppFunction(
  code = '
  void fetch_knots(Rcpp::List& Knots, int depth, int& max_depth, std::vector<int> knot_vec, int ref, int& epsilon, const char* file){
    if (!)
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
        for (int i=0; i < knot_vec_keeper.size(); i++){
          Rcout << knot_vec_keeper.at(i) << " " ;
        }
        Rcout << "\\n";
      }
    }
  }
  '
)

#NumericMatrix::Row row = m(m_iter , _ );
#row = knot_vec_keeper;
#m_iter++;
i <- 1
m <- matrix(nrow = 20, ncol = 3)
profile.6 <- lineprof(fetch_knots(Knots, 0, max_depth, 0, 0, epsilon, i, m))

Knots <- rep(list(1:30), 3)
max_depth <- length(Knots) - 1
epsilon = 0
fetch_knots(Knots, 0, max_depth, 0, 0, epsilon, i, m)
