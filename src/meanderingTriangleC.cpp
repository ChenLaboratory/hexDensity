#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>
#include <map>
#include <vector>
#include <set>
#include <deque>
#include <cmath> //for std::pow
#include <algorithm> //for std::set_intersection
struct point {
  int x = -1;
  int y = -1;
};

struct point_d {
  double x;
  double y;
};

// for map of point_d
bool operator<(const point_d& l, const point_d& r) {
  // compare x then y
  return (l.x<r.x || (l.x==r.x && l.y<r.y));
}

bool operator==(const point_d& l, const point_d& r) {
  // compare x then y
  return (l.x==r.x && (l.y==r.y));
}

struct triangle {
  point v[3];
};

struct edge {
  point_d e[2];
};
// for set of edge
bool operator<(const edge& l, const edge& r) {
  // compare first edge then second edge
  return (l.e[0]<r.e[0] || ((!(r.e[0]<l.e[0])) && l.e[1]<r.e[1]));
}

// for set intersection with &&
template <class T, class CMP = std::less<T>, class ALLOC = std::allocator<T> >
std::set<T, CMP, ALLOC> operator && (
    const std::set<T, CMP, ALLOC> &s1, const std::set<T, CMP, ALLOC> &s2)
{
  std::set<T, CMP, ALLOC> s;
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::inserter(s, s.begin()));
  return s;
}

//format isoline output as an R list with x,y, and id.
SEXP format_output(double* x, double* y, int* id, int n) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP,3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP,3));
  SET_STRING_ELT(names,0,Rf_mkChar("x"));
  SET_STRING_ELT(names,1,Rf_mkChar("y"));
  SET_STRING_ELT(names,2,Rf_mkChar("id"));
  Rf_setAttrib(out, Rf_install("names"), names);
  
  SEXP x_final = SET_VECTOR_ELT(out,0,Rf_allocVector(REALSXP,n));
  double* x_final_p = REAL(x_final);
  SEXP y_final = SET_VECTOR_ELT(out,1,Rf_allocVector(REALSXP,n));
  double* y_final_p = REAL(y_final);
  SEXP id_final = SET_VECTOR_ELT(out,2,Rf_allocVector(INTSXP,n));
  int* id_final_p = INTEGER(id_final);
  
  for (int i=0;i<n;i++) {
    x_final_p[i] = x[i];
    y_final_p[i] = y[i];
    id_final_p[i] = id[i];
  }
  UNPROTECT(2);
  return out;
}

extern "C" {
  SEXP meanderingTrianglesC(SEXP xleft, SEXP xright, SEXP y,SEXP z,SEXP levels) {
    xleft = PROTECT(Rf_coerceVector(xleft,REALSXP));
    xright = PROTECT(Rf_coerceVector(xright,REALSXP));
    y = PROTECT(Rf_coerceVector(y,REALSXP));
    z = PROTECT(Rf_coerceVector(z,REALSXP));
    levels = PROTECT(Rf_coerceVector(levels,REALSXP));
    
    double *xleft_p,*xright_p,*y_p,*levels_p,*z_p;
    xleft_p = REAL(xleft);
    xright_p = REAL(xright);
    y_p = REAL(y);
    z_p = REAL(z);
    int z_nrows=Rf_nrows(z);
    levels_p = REAL(levels);
    // Get triangles. TODO: Could probably do away with this.
    int n = (Rf_length(xleft)-1)*(Rf_length(y)-1)*2;
    std::vector<triangle> triangles(n);
    int i_triangle = 0;
    
    for (int i=0;i<(Rf_length(xleft)-1);i++) {
      for (int j=0;j<(Rf_length(y)-1);j++) {
        if (j%2==0) {

          triangles[i_triangle++] = {{{i,j},{i+1,j},{i,j+1}}};
          triangles[i_triangle++] = {{{i+1,j},{i,j+1},{i+1,j+1}}};
        } else{
          triangles[i_triangle++] = {{{i,j},{i+1,j+1},{i,j+1}}};
          triangles[i_triangle++] = {{{i,j},{i+1,j+1},{i+1,j}}};
        }
      }
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, Rf_length(levels)));
    for (int l=0;l<Rf_length(levels);l++) {
      std::map<point_d,point_d> interpolatedPos;
      std::vector<edge> contour_segments;
      for(triangle t:triangles) {
        // Find contour line in the triangle
        int score = 0;
        for (int i=0;i<3;i++) {
          // No contour if any corners are NA
          if (!R_finite(z_p[t.v[i].y+t.v[i].x*z_nrows])) {
            score = 0;
            break;
          }
          // R matrix is col-major alignment
          score += (z_p[t.v[i].y+t.v[i].x*z_nrows]>=levels_p[l])*std::pow(2,i);
        }
        point minor;
        point major[2];
        
        switch (score) {
          case 0: // [0,0,0] No contour line
          case 7: // [1,1,1] No contour line
            continue;
          case 1: // [1,0,0]
          case 6: // [0,1,1]
            minor = t.v[0];
            major[0] = t.v[1];
            major[1] = t.v[2];
            break;
          case 2: // [0,1,0]
          case 5: // [1,0,1]
            minor = t.v[1];
            major[0] = t.v[0];
            major[1] = t.v[2];
            break;
          case 3: // [1,1,0]
          case 4: // [0,0,1]
            minor = t.v[2];
            major[0] = t.v[0];
            major[1] = t.v[1];
            break;
        }
        
        point_d crossing_point;
        point_d interpolated_point;
        double how_far;
        point_d contour_point[2];
        for (int m=0;m<2;m++) {
          /* temp coordinates based on row & col to avoid floating point error 
          when joining*/
          crossing_point = {(minor.x-major[m].x)*0.5+major[m].x,
                            (minor.y-major[m].y)*0.5+major[m].y};
          how_far = (levels_p[l]-z_p[major[m].y+major[m].x*z_nrows])/
            (z_p[minor.y+minor.x*z_nrows]-z_p[major[m].y+major[m].x*z_nrows]);
          
          // Actual coordinates
          double minor_x = minor.y%2==0 ? xleft_p[minor.x] : xright_p[minor.x];
          double major_x = major[m].y%2==0 ? xleft_p[major[m].x]:xright_p[major[m].x];
          interpolated_point = {
            (minor_x-major_x)*how_far+major_x,
            (y_p[minor.y]-y_p[major[m].y])*how_far+y_p[major[m].y]
          };
          interpolatedPos[crossing_point] = interpolated_point;
          contour_point[m] = crossing_point;
        }
        contour_segments.push_back({contour_point[0],contour_point[1]});
      }
      // No contour line at this level
      if (contour_segments.size()==0) {
        double x_out[1] = {-1};
        double y_out[1] = {-1};
        int id[1] = {-1};
        SET_VECTOR_ELT(out,l,format_output(x_out,y_out,id,0));
        continue;
      }
      // Joining up
      std::set<edge> unused_segments(contour_segments.begin(),contour_segments.end());
      std::map<point_d,std::set<edge>> segments_by_point;
      for(edge segment:contour_segments) {
        segments_by_point[segment.e[0]].insert(segment);
        segments_by_point[segment.e[1]].insert(segment);
      }
      
      std::vector<std::deque<point_d>> contour_lines;
      int n_out = 0;
      while (!unused_segments.empty()) {
        std::deque<point_d> line;
        // Pop a random segment
        line.push_back((*unused_segments.begin()).e[0]);
        line.push_back((*unused_segments.begin()).e[1]);
        unused_segments.erase(unused_segments.begin());
        // Extend line both ways
        while (true) {
          std::set<edge> tail_candidates = segments_by_point[line.back()]&&unused_segments;
          if (tail_candidates.size()>0) {
            edge tail = *tail_candidates.begin();
            line.push_back(tail.e[1] == line.back() ? tail.e[0] : tail.e[1]);
            unused_segments.erase(tail);
          }
          std::set<edge> head_candidates = segments_by_point[line.front()]&&unused_segments;
          if (head_candidates.size()>0) {
            edge head = *head_candidates.begin();
            line.push_front(head.e[1] == line.front() ? head.e[0] : head.e[1]);
            unused_segments.erase(head);
          }
          if (tail_candidates.size()==0 && head_candidates.size()==0) {
            n_out += line.size();
            contour_lines.push_back(line);
            break;
          }
        }
      }
      
      // Collect result to isoline format
      double *x_out = new double[n_out];
      double *y_out = new double[n_out];
      int *id_out = new int[n_out];
      n_out = 0;
      for (std::size_t i=0;i<contour_lines.size();i++) {
        for(point_d p:contour_lines[i]) {
          // Scale to correct xy.
          point_d correct_point = interpolatedPos[p];
          x_out[n_out] = (correct_point.x);
          y_out[n_out] = (correct_point.y);
          
          id_out[n_out] = i+1;
          n_out++;
        }
      }
      SET_VECTOR_ELT(out,l,format_output(x_out,y_out,id_out,n_out));
      delete [] x_out;
      delete [] y_out;
      delete [] id_out;
    }
    
    UNPROTECT(6);
    return out;
  }
}

