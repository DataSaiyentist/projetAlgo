#include <Rcpp.h>
using namespace Rcpp;
#include "FOCuS.h"

// Propriété @gtromano
// https://github.com/gtromano/FOCuS/tree/master/src

// This is the script for exporting the C++ functions to R.

// converting c++ list into R list
std::list<List> convert_output_to_R(const std::list<Quadratic>& c_obj) {
  std::list<List> output;

  //std::cout << "last Q1"<< std::endl;
  for (auto& q:c_obj) {
    // coversion of the data from c++ into R
    // print(q);

    std::list<List> ints(q.ints.size());
    std::transform(q.ints.begin(), q.ints.end(), ints.begin(),
                   [](const auto &i){
                     return List::create(Rcpp::Named("l") = i.l,
                                         Rcpp::Named("u") = i.u);
                   });

    auto l = List::create(Rcpp::Named("a") = q.a,
                          Rcpp::Named("b") = q.b,
                          Rcpp::Named("c") = q.c,
                          Rcpp::Named("ints") = ints,
                          Rcpp::Named("max") = q.max);
    l.attr("class") = "Quadratic";
    output.push_back(l);
  }
  return output;
}


/* ------------------------------------------------------------

 Online version - Rcpp wrapper

 -------------------------------------------------------------- */

// [[Rcpp::export(.FoCUS)]]
List FOCuS (Rcpp::Function dataGen, const double thres, const double& mu0, std::list<double>& grid, const double& K) {
  if (!std::isnan(grid.front())) {
    grid.push_back(INFINITY);
    grid.push_front(-INFINITY);
  }

  long t = 0;
  long cp = -1;

  Quadratic Q0, q1;
  Info info = {Q0, {q1}, 0};

  try {
    // if we don't know the pre-change mean then replace the f with a lambda
    if (std::isnan(mu0)) {
      while(true) {
        t++;
        double y = Rcpp::as<double>(dataGen());
        info = FOCuS_step(std::move(info), y, grid, K);
        if (info.global_max >= thres) {
          cp = t;
          break;
        }
      }
    } else {
      while(true) {
        t++;
        double y = Rcpp::as<double>(dataGen());
        info = FOCuS_step_sim(std::move(info), y - mu0, grid, K);
        if (info.global_max >= thres) {
          cp = t;
          break;
        }
      }
    }

  }
  catch (std::bad_alloc &e) {
    Rcpp::stop("insufficient memory");
  }
  catch (...) {
    auto last_Q1 = convert_output_to_R(info.Q1);
    return List::create(Rcpp::Named("t") = cp,
                        Rcpp::Named("Q1") = last_Q1,
                        Rcpp::Named("warning_message") = "The procedure was interrupted or terminated unexpectedly. The output was successfully returned, however there is a possibility it can be possibly corrupted.");;
  }


  auto last_Q1 = convert_output_to_R(info.Q1);
  return List::create(Rcpp::Named("t") = cp,
                      Rcpp::Named("Q1") = last_Q1);

}
