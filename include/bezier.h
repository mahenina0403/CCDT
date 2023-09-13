#ifndef BEZIER_H
#define BEZIER_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Vector_2                                             Point_2;
typedef std::vector<Point_2>                                    Bezier_curve;

double norm(Point_2 p);
Point_2 normal(Point_2 p);
double cross(Point_2 a, Point_2 b);
Point_2 meeting_point(Point_2 P1, Point_2 d1, Point_2 P2, Point_2 d2);
double get_distance(Bezier_curve P);
Point_2 point_on_bezier_curve(Bezier_curve P, double t);
void subdivide_bezier_curve(Bezier_curve P, Bezier_curve &P1, Bezier_curve &P2, double t);
Point_2 circumcenter(Bezier_curve Q, double &t);
void sample_bezier_curve(Bezier_curve P, double start, double end, int n_sample, std::vector<Point_2> &samples);

#endif