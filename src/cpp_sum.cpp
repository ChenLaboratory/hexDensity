#include <R.h>
#include <Rinternals.h>
#include <map>
#include <vector>
#include <algorithm> 
#include <numeric>
#include <set>
#include <deque>
#include <cmath>
#include <iostream>
struct point {
  int x;
  int y;
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
  return (l.e[0]<r.e[0] || (!(r.e[0]<l.e[0])) && l.e[1]<r.e[1]);
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

struct isoline {
  double *x;
  double *y;
  int *id;
};
//format isoline output as an R list with x,y, and id.
SEXP format_output(double* x, double* y, int* id, int n);
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

// Test matrix and output
extern "C" SEXP expSmooth(SEXP y, SEXP ys, SEXP z, SEXP levels);
SEXP expSmooth(SEXP x, SEXP y, SEXP z, SEXP levels) {
  x = PROTECT(coerceVector(x,REALSXP));
  y = PROTECT(coerceVector(y,REALSXP));
  z = PROTECT(coerceVector(z,REALSXP));
  
  int levels_n = Rf_length(levels);
  int z_nrows = Rf_nrows(z);
  // int z_ncols = Rf_ncols(z);
  double *x_p, *y_p, *r_p;
  
  x_p = REAL(x);
  y_p = REAL(y);
  r_p = REAL(z);

  SEXP out = PROTECT(Rf_allocVector(VECSXP,levels_n));
  
  // Test data
  double x_out[3] = {r_p[0+0*z_nrows],r_p[1+0*z_nrows],r_p[2+0*z_nrows]};
  double y_out[3] = {r_p[2+0*z_nrows],r_p[2+1*z_nrows],r_p[2+2*z_nrows]};
  int id[3] = {0,0,1};
  
  for (int i=0;i<levels_n;i++) {

    SET_VECTOR_ELT(out,i,format_output(x_out,y_out,id,3));
  }
  UNPROTECT(4);
  return out;
}

#include <iostream>
extern "C" {
  void expSmooth2() {
    std::map<int,std::set<int>>test_map;
  }
}

extern "C" {
  SEXP meanderingTrianglesC(SEXP x, SEXP y,SEXP z,SEXP levels) {
    // Not sure if actually need to coerce all these to REALSXP
    x = PROTECT(coerceVector(x,REALSXP));
    y = PROTECT(coerceVector(y,REALSXP));
    z = PROTECT(coerceVector(z,REALSXP));
    levels = PROTECT(coerceVector(levels,REALSXP));
    
    double *rx,*ry,*rlevels,*rz;
    rx = REAL(x);
    ry = REAL(y);
    rz = REAL(z);
    rlevels = REAL(levels);
    for(int i =0;i<length(z);i++) {
      std::cout<<rz[i]<<std::endl;
    }
    // Get triangles
    int n = (length(x)-1)*(length(y)-1)*2;
    std::vector<triangle> triangles(n);
    int i_triangle = 0;
    
    for (int i=0;i<(length(x)-1);i++) {
      for (int j=0;j<(length(y)-1);j++) {
        if (j%2==0) {

          triangles[i_triangle++] = {{{i,j},{i+1,j},{i,j+1}}};
          triangles[i_triangle++] = {{{i+1,j},{i,j+1},{i+1,j+1}}};
        } else{
          triangles[i_triangle++] = {{{i,j},{i+1,j+1},{i,j+1}}};
          triangles[i_triangle++] = {{{i,j},{i+1,j+1},{i+1,j}}};
        }
      }
    }
    std::cout << "There are " << triangles.size() << " triangles" << std::endl;
    
    // isoline res[length(levels)];
    SEXP out = PROTECT(allocVector(VECSXP, length(levels)));
    int out_i=0;
    
    for (int l=0;l<length(levels);l++) {
      std::map<point_d,point_d> interpolatedPos;
      std::vector<edge> contour_segments;
      std::cout << "level is " << rlevels[l] <<std::endl;
      for(triangle t:triangles) {
        std::cout << "Triangle " <<t.v[0].x<<t.v[0].y<<t.v[1].x<<t.v[1].y<<t.v[2].x<<t.v[2].y<<std::endl;
        // Find contour line in the triangle
        int score = 0;
        for (int i=0;i<3;i++) {
          // Matrix is row-major alignment
          std::cout << "At index " << t.v[i].x+t.v[i].y*Rf_nrows(z) << " z is " << (rz[t.v[i].x+t.v[i].y*Rf_nrows(z)]) << std::endl;
          // score += (rz[t.v[i].y+t.v[i].x*Rf_nrows(z)]>=rlevels[l])*pow(2,i);
          score += (rz[t.v[i].x+t.v[i].y*Rf_nrows(z)]>=rlevels[l])*pow(2,i);
        }
        std::cout << "score is: " << score <<std::endl;
        point minority;
        point majority[2];
        
        switch (score) {
          case 0: // [0,0,0] No contour line
          case 7: // [1,1,1] No contour line
            continue;
          case 1: // [1,0,0]
          case 6: // [0,1,1]
            minority = t.v[0];
            majority[0] = t.v[1];
            majority[1] = t.v[2];
            break;
          case 2: // [0,1,0]
          case 5: // [1,0,1]
            minority = t.v[1];
            majority[0] = t.v[0];
            majority[1] = t.v[2];
            break;
          case 3: // [1,1,0]
          case 4: // [0,0,1]
            minority = t.v[2];
            majority[0] = t.v[0];
            majority[1] = t.v[1];
            break;
        }
        
        point_d crossing_point;
        point_d interpolated_point;
        double how_far;
        point_d contour_point[2];
        for (int m=0;m<2;m++) {
          crossing_point = {(minority.x-majority[m].x)*0.5+majority[m].x,
                            (minority.y-majority[m].y)*0.5+majority[m].y};
          how_far = (rlevels[l]-rz[majority[m].x+majority[m].y*Rf_nrows(z)])/
            (rz[minority.x+minority.y*Rf_nrows(z)]-rz[majority[m].x+majority[m].y*Rf_nrows(z)]);
          interpolated_point = {(minority.x-majority[m].x)*how_far+majority[m].x,
                                (minority.y-majority[m].y)*how_far+majority[m].y};
          interpolatedPos[crossing_point] = interpolated_point;
          contour_point[m] = crossing_point;
        }
        contour_segments.push_back({contour_point[0],contour_point[1]});
      }
      
      if (contour_segments.size()==0) {
        // cpp doesnt allow temporary array
        double x_out[0] = {};
        double y_out[0] = {};
        int id[0] = {};
        SET_VECTOR_ELT(out,out_i,format_output(x_out,y_out,id,0));
        out_i++;
        continue;
      }
      std::cout << "Number of contour segments is: " << contour_segments.size() << std::endl;
      // Joining up
      std::cout << "Start joining up" << std::endl; 
      std::set<edge> unused_segments(contour_segments.begin(),contour_segments.end());
      std::map<point_d,std::set<edge>> segments_by_point;
      for(edge segment:contour_segments) {
        segments_by_point[segment.e[0]].insert(segment);
        segments_by_point[segment.e[1]].insert(segment);
      }
      std::vector<std::deque<point_d>> contour_lines;
      int n_out = 0;
      std::cout << "Number of unused segments is: " << unused_segments.size() << std::endl;
      while (!unused_segments.empty()) {
        std::cout <<  unused_segments.size() << " unused segments left" << std::endl;
        std::deque<point_d> line;
        // Pop a random segment
        line.push_back((*unused_segments.begin()).e[0]);
        line.push_back((*unused_segments.begin()).e[1]);
        unused_segments.erase(unused_segments.begin());
        
        while (true) {
          std::set<edge> tail_candidates = segments_by_point[line.back()]&&unused_segments;
          // std::cout <<  tail_candidates.size() << " tail_candiates" << std::endl;
          if (tail_candidates.size()>0) {
            edge tail = *tail_candidates.begin();
            tail_candidates.erase(tail_candidates.begin()); // TODO: Might not need this
            line.push_back(tail.e[1] == line.back() ? tail.e[0] : tail.e[1]);
            unused_segments.erase(tail);
          }
          std::set<edge> head_candidates = segments_by_point[line.front()]&&unused_segments;
          if (head_candidates.size()>0) {
            edge head = *head_candidates.begin();
            head_candidates.erase(head_candidates.begin()); // TODO: might not need this
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
      std::cout << "There are " << contour_lines.size() << " contour lines" << std::endl; 
      std::cout << "Finish joining up" << std::endl; 
      std::cout << n_out << std::endl; 
      double x_out[n_out];
      double y_out[n_out];
      int id_out[n_out];
      n_out = 0;
      std::cout << "Start convert to R" << std::endl; 
      for (int i=0;i<contour_lines.size();i++) {
        for(point_d p:contour_lines[i]) {
          x_out[n_out] = p.x;
          y_out[n_out] = p.y;
          id_out[n_out] = i+1; //TODO: not sure if need +1 here
          n_out++;
        }
      }
      SET_VECTOR_ELT(out,out_i,format_output(x_out,y_out,id_out,n_out));
      out_i++;
    }
    
    UNPROTECT(5);
    return out;
  }
}

// TODO: delete this
// Making & naming list
// library(inline)
//   named <- cfunction(signature(), '
//                        /* allocate and populate list */
//                        SEXP OS = PROTECT(allocVector(VECSXP, 2));
//                        SET_VECTOR_ELT(OS, 0, allocMatrix(REALSXP, 5, 5));
//                        SET_VECTOR_ELT(OS, 1, allocVector(REALSXP, 5));
// 
//                        /* create names */
//                        SEXP nms = PROTECT(allocVector(STRSXP, 2));
//                        SET_STRING_ELT(nms, 0, mkChar("foo"));
//                        SET_STRING_ELT(nms, 1, mkChar("bar"));
// 
//                        /* assign names to list */
//                        setAttrib(OS, R_NamesSymbol, nms);
// 
//                        /* cleanup and return */
//                        UNPROTECT(2);
//                        return OS;')