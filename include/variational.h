#include "refinement.h"

double orientation(Point p1, Point p2, Point p3);
bool is_inside(std::vector<Point> polygon, Point p);
bool do_intersect (std::vector<Point> polygon1, std::vector<Point> polygon2);
Point_2 guard_point(Bezier_curve P);
std::vector<Point> guard_envelope(Bezier_curve P);
std::vector<Bezier_curve> make_all_guardable(std::vector<Bezier_curve> curves);
std::vector<Bezier_curve> make_no_guard_intersection(std::vector<Bezier_curve> curves, int count=0);
std::vector<Bezier_curve> target_length(std::vector<Bezier_curve> curves);