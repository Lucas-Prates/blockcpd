#include <Rcpp.h>
using namespace Rcpp;

// Auxiliary function
// Computes the minimum distance for the points of A from B.
// Used in Hausdorff Distance
// Note: this function is NOT symmetric!!!
// For example, let A={1,2,3} and B = {2}. Then
// comp_distance(A, B) = 1 but comp_distance(B, A) = 0.
int comp_distance(IntegerVector  A, IntegerVector B){
  int l_A = A.length();
  int l_B = B.length();
  int j = 0;
  int max_dist = 0, dist;
  for(int i=0; i < l_A; i++){
    if(j < l_B - 1){
      while(A[i] > B[j]){
        j++;
        if (j == l_B - 1){
          break;
        }
      }
    }
    if(j == 0){
      dist = abs(A[i]-B[j]);
    }else{
      dist = std::min(abs(A[i]-B[j-1]), abs(A[i]-B[j]));
    }
    max_dist = std::max(max_dist, dist);
  }

  return max_dist;
}

//' @title
//' Rand Index Function for change point detection
//'
//' @description
//' Computes the rand Index (non-adjusted) for the change point sets. A specific
//' equation for change point detection is used to make the computation faster.
//' Proof of correctness of the equation is given in the dissertation.
//'
//' @param cp1 Change point set for model 1 or true change point set.
//' @param cp2 Change point set for model 2 or true change point set.
//' @param m The size of the vector array.
//'
// [[Rcpp::export]]
float compute_rand(IntegerVector cp1,
                         IntegerVector cp2,
                         int const& m){
  cp1.push_front(0);
  cp1.push_back(m);
  cp2.push_front(0);
  cp2.push_back(m);
  int len_cp1 = cp1.length();
  int len_cp2 = cp2.length();
  float M = 0.0; // This is the total sum of M_ij|c1_i - c2_j| from original eq.
  float M_ij;
  for(int i = 0; i < len_cp1 - 1; i++){
    for(int j = 0; j < len_cp2 - 1; j++){
      M_ij = std::max(0, std::min(cp1[i + 1], cp2[j + 1]) - std::max(cp1[i], cp2[j]));
      M += abs(cp1[i]-cp2[j])*M_ij;
    }
  }

  M /= (m*(m - 1)/2);

  return(1 - M);

};


//' @title
//' Hausdorff distance metric
//'
//' @description
//' Computes the Hausdorff distance between change point sets.
//'
//' @param cp1 Change point set for model 1 or true change point set.
//' @param cp2 Change point set for model 2 or true change point set.
//[[Rcpp::export]]
int compute_hausdorff(IntegerVector cp1,
                      IntegerVector cp2){
  int dist = std::max(comp_distance(cp1, cp2), comp_distance(cp2, cp1));
  return dist;

}

//' @title
//' Symmetric difference metric
//'
//' @description
//' Computes the size of the symmetric difference between two change point
//' detection sets
//'
//' @param cp1 Change point set for model 1 or true change point set.
//' @param cp2 Change point set for model 2 or true change point set.
//[[Rcpp::export]]
int compute_symdiff(IntegerVector cp1,
                    IntegerVector cp2){
  IntegerVector cp_inter = intersect(cp1, cp2);

  int sd = cp1.length() + cp2.length() - 2*cp_inter.length();
  return sd;

}

//' @title
//' Jaccard's Index metric
//'
//' @description
//' Computes the Jaccard index between two change point
//' detection sets
//'
//' @param cp1 Change point set for model 1 or true change point set.
//' @param cp2 Change point set for model 2 or true change point set.
//[[Rcpp::export]]
double compute_jaccard(IntegerVector cp1,
                    IntegerVector cp2){
  int size_inter = intersect(cp1, cp2).length();
  int size_union = union_(cp1, cp2).length();
  double jac_ind = (1.0*size_inter)/size_union;
  return jac_ind;

}

//------------------------------------------------------------------------------
