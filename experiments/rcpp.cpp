#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void open_file(const char* file){
  std::ofstream ofs;
  ofs.open(file, std::ios_base::app);
}

// [[Rcpp::export]]
void fetch_knots(Rcpp::List& Knots, int depth, int& max_depth, std::vector<int> knot_vec, int ref, int& epsilon, std::ofstream& ofs){
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
      for (i = 1; i < knot_vec_keeper.size(); ++i){
        
      }
    }
  }
}


// [[Rcpp::export]]
void close_file(const char* file){
  std::ofstream file;
}