#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Vector_2                                             Point_2;
typedef std::vector<Point_2>                                    Bezier_curve;

double norm(Point_2 p){

    /*
    Retur: norm of a vector
    */

    return CGAL::to_double(CGAL::sqrt(p.x()*p.x()+p.y()*p.y()));
}

Point_2 normal(Point_2 p){
  return Point_2(p.y(), -p.x());
}

// We check the orientation of two vectors.
double cross(Point_2 a, Point_2 b){

    /*
    Return: cross product between two vectors
    if the result is negative: it is a clockwise rotation
    if the result is positive: it is a anti-clockwise rotation
    */

    return CGAL::to_double(a.x()*b.y()-a.y()*b.x());
}

Point_2 meeting_point(Point_2 P1, Point_2 d1, Point_2 P2, Point_2 d2){

    /*
    Return: meeting points of the line passing through P1, with direction d1
            and the line through P2, with direction d2
    */
    
    double detA = -cross(d1, d2);
    double a = cross(P2-P1, d2)/detA;
    return (P1 + a*d1);
    // double a = (-d2.y() * (P2.x()-P1.x()) + d2.x() * (P2.y()-P1.y()))/detA;
    // return Rat_point_2(P1.x()+a*d1.x(), P1.y()+a*d1.y());
}

double get_distance(Bezier_curve P){

    // Return: the distance between the endpoints

    Point_2 p0 = P[0];
    Point_2 pn = P.back();
    return norm(p0-pn);
}

Point_2 point_on_bezier_curve(Bezier_curve P, double t){
  Point_2 N;
  double d;

  int n = P.size()-1;
  
  std::vector<double> w;
  w.push_back(1);
  for (int i=1; i<n+1; i++) w.push_back(w.back()*(n-i+1.0)/i);

  if (t > 0.5){
    double s = (1.0-t)/t;

    N = P[0];
    d = w[0];

    for (int i=1; i<n+1; i++){
      // w = w * (n-i+1.0)/i;
      N = N*s + w[i]*P[i];
      d = d*s + w[i];
    }
  }else{
    double s = t / (1.0-t);

    N = P[n];
    d = w[n];

    for (int i=0; i<n; i++){
      int j = n-1-i;
      // w = w * (n-j+1.0)/j;
      N = N*s + w[j]*P[j];
      d = d*s + w[j]; 
    }
  }
  return N/d;
}

void subdivide_bezier_curve(Bezier_curve P, Bezier_curve &P1, Bezier_curve &P2, double t=0.5){
  int n = P.size()-1;
  Bezier_curve tmp1;
  Bezier_curve tmp2 = P;
  P1.clear();
  P2.clear();

  for (int i=0; i<n+1; i++){
    tmp1.push_back(P[i]);
    P1.push_back(point_on_bezier_curve(tmp1, t));
    P2.push_back(point_on_bezier_curve(tmp2, t));
    tmp2.erase(tmp2.begin());
  }
}

Bezier_curve degree_elevate(Bezier_curve P){
  int n = P.size()-1;

  Bezier_curve tmp;

  tmp.push_back(P[0]);
  for (int i=0; i<n; i++){
    auto pi = (P[i+1] + P[i])/2.0;
    tmp.push_back(pi);
  }
  tmp.push_back(P.back());
  return tmp;
}

Point_2 circumcenter(Bezier_curve Q, double &t){
  int n;
  double d;
  double di;
  Bezier_curve P = Q;
elevate:
  n = P.size()-1;
  d = 0;
  for(int i=0; i<n; i++){
    di = norm(P[i+1]-P[i]);
    if (di > d) d = di;
  }
  
  
  std::cout << "cout " << d << std::endl;
  if (d > 5){
    P = degree_elevate(P);
    goto elevate;
  }

  Point_2 I = (P.back()+P[0])/2;
  Point_2 bisector = normal(P.back()-P[0]);

  for (int i=0; i<n; i++){
    auto si = P[i+1]-P[i];

    double detA = -cross(si, bisector);
    double a = cross(I-P[i], bisector)/detA;
    if ( 0 < a && a < 1){
      t = i/(n+1.0);
      return point_on_bezier_curve(Q, t);
    }
  }
}

void sample_bezier_curve(Bezier_curve P, double start, double end, int n_sample, std::vector<Point_2> &samples){
  samples.clear();
  double dt = (end-start)/n_sample;

  double t = start;
  while(t < end){
    Point_2 p = point_on_bezier_curve(P, t);

    samples.push_back(p);

    t = t + dt;
  }
}