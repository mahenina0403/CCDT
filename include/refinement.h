#include <CGAL/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/enum.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/Real_timer.h>

#include <QApplication>
#include <QPixmap>
#include <QPainter>
#include <QPen>
#include <QPoint>
#include <QLabel>
#include <QPolygon>
#include <QBrush>

#include "bezier.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel                                 K;
typedef CGAL::Exact_predicates_tag                                                          Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag>                  CDT;
typedef CDT::Point Point;
typedef CDT::Edge                                                                           Edge;
typedef CDT::Vertex_handle                                                                  Vertex_handle;
typedef CDT::Face_handle                                                                    Face_handle;
typedef CDT::Triangle                                                                       Triangle;
typedef CGAL::Vector_2 < K >                                                                Vector_2;
typedef CGAL::Polygon_2< K >                                              			        Polygon_2;
typedef CGAL::Segment_2 < K >                                                               Segment_2;
#define BOUND CGAL::sqrt(2.0)

Point to_point(Point_2 p);
Point_2 to_point_2(Point p);
double cos_of_angle(Point_2 a, Point_2 b);
double norm(Edge e);
std::vector<Point_2> curve_axis(Bezier_curve P);
bool is_guardable(Bezier_curve P, Point_2* splus, Point_2* sminus);
bool is_valid(Bezier_curve Bezier, Point_2 *splus, Point_2 *sminus);
std::vector<Bezier_curve> bisect(std::vector<Bezier_curve> curves, int i, Point &p, double t=0.5);
std::vector<Bezier_curve> add_boundary(std::vector<Bezier_curve> curves);
std::vector<Bezier_curve> make_all_guard_valid(std::vector<Bezier_curve> curves);
std::vector<Bezier_curve> make_all_valid(std::vector<Bezier_curve> curves);
std::vector<Point_2> get_corners(std::vector<Bezier_curve> curves);
bool is_curved(Edge e, std::vector<Bezier_curve> curves);
CDT triangulation(std::vector<Bezier_curve> curves, bool variational=false);
CDT refine_curve(CDT T, std::vector < Bezier_curve > *curves, std::vector<std::pair<Point_2, double>> corners);
CDT refine_encroached(CDT T, std::vector < Bezier_curve > *curves, std::vector<std::pair<Point_2, double>> corners);
bool is_encroached(Point p, std::vector < Bezier_curve > curves, int &index);
bool is_valid_triangle(CDT T, Face_handle f);
bool comparator(std::pair < Face_handle, double> f1, std::pair < Face_handle, double > f2);
std::vector < std::pair <Face_handle, double>>
insert(std::vector < std::pair <Face_handle, double>> M, std::pair<Face_handle, double> tri);
void
protecting_balls(std::vector < Bezier_curve > *curves, std::vector<Point_2> C , std::vector<std::pair<Point_2, double>> *corners);
CDT
refine_bad_triangle(CDT T, std::vector < Bezier_curve > *curves, std::vector<std::pair<Point_2, double>> corners);
CDT split_constrained_triangles(CDT T, std::vector<std::pair<Point_2, double>> corners);
void
draw(std::vector<Bezier_curve> curves, CDT T, int frame);
int
draw(std::vector<Bezier_curve> curves, int argc, char *argv[]);
bool is_a_corner(Point_2 p, std::vector<Point_2> corners);
CDT refine_edge_length(CDT T, std::vector < Bezier_curve > curves);