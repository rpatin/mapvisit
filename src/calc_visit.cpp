#include <Rcpp.h>
using namespace Rcpp;



/* from adehabitatHR. Author : C. Calenge
 xi and yi are the center of the circle of radius r
 x1 and y1 are the coordinates of the location inside the circle
 x2 and y2 are the coordinates of the location outside
 */
double interpLoc(double xi, double yi, double x1, double y1,
                 double x2, double y2, double r)
{
  double a, u, v, p, d, g;

  d = hypot(x2-x1, y2-y1);
  a = atan2(y2-y1, x2-x1);
  u = ((xi-x1)*cos(a))+((yi-y1)*sin(a));
  v = ((yi-y1)*cos(a))-((xi-x1)*sin(a));
  g = sqrt(pow(r, 2.0) - pow(v, 2.0));
  p = (u + g)/d;
  return(p);
}



/* Functions modified from adehabitatHR nvisits and HRresidtime */
// [[Rcpp::export]]
List calc_visit_cpp(NumericVector x, NumericVector y, NumericVector t, NumericVector xgrid,
                    NumericVector ygrid, double dist, double maxt)
{
  /* declaring the variables */
  int n = x.size();  /* number of relocations */
  int ngrid = xgrid.size(); /* number of points on the grid */
  IntegerVector deds(n); /* will be used to check whether
   the relocation j
   is within the distance dist
   of the grid point i */
  int sortie = 1;
  IntegerVector nvisit(ngrid); /* the output vector for number of visits */

  double time_enter;
  double time_leave;
  double p;
  NumericVector residtime(ngrid); /* the output vector for the residence time */

  /* Now, calculate the residence time for each relocation */
  for (int i = 0; i < ngrid; ++i) {
    /* No visit at the beginning */
    nvisit(i) = 0;

    /* checks which relocations are within the distance dist from grid point i  */

    for (int j = 0; j < n; ++j) {
      if (hypot(xgrid(i)-x(j), ygrid(i)-y(j))<=dist) {
        deds(j) = 1;
      } else {
        deds(j) = 0;
      }
    }

    /* forward time */
    sortie = 1; /* = 0 when the animal is still
    inside the circle; =1 otherwise */
    time_enter = t(n-1)+maxt;
    time_leave = -maxt-1;
    residtime(i) = 0;
    // Rprintf( "k loop \\ n" );
    /* if this is not the last relocation (otherwise, n forward visit = 0.0) */
    /* for all next relocations */

    for (int k = 0; k<n; ++k) {
      /* if the relocation is outside the circle */
      if (deds[k]==0) {
        /* if this is the first relocation outside */
        if (sortie==0) {
          sortie = 1;
          /* interpolating the time when the animal come out
          of the circle. Same procedure as above
          */
          p = interpLoc(xgrid(i), ygrid(i), x(k-1), y(k-1),
                        x(k), y(k), dist);
          time_leave = t(k-1) + p*(t(k)-t(k-1));
          residtime(i) += time_leave - time_enter;
        }
      } else { /* if the relocation is inside the circle */
          if(sortie>0) {
            sortie = 0;
            if( k == 0){
              time_enter = t(0);
            } else {
              p = interpLoc(xgrid(i), ygrid(i), x(k), y(k),
                            x(k-1), y(k-1), dist);
              time_enter =  t(k) + p*(t(k-1) - t(k));
            }

            if (fabs(time_enter - time_leave) > maxt) {
              nvisit(i)++;
            }
          } else if(k == n-1){
            time_leave = t(k-1);
            residtime(i) += time_leave - time_enter;
          }

      }
    }
  }
  // Rprintf( "end \\ n" );


  return List::create(
    _["nvisit"] = nvisit,
    _["residtime"] = residtime);
}
