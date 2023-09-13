#include "variational.h"

double orientation(Point p1, Point p2, Point p3){
    Triangle T = Triangle(p1,p2,p3);
    return T.orientation();
}

bool is_inside(std::vector<Point> polygon, Point p){

    auto T1 = Triangle(polygon[0], polygon[1], polygon[2]);
    auto T2 = Triangle(polygon[3], polygon[1], polygon[2]);

    if (T1.has_on_bounded_side(p))
        return true;

    if (T2.has_on_bounded_side(p))
        return true;
    
    return false;

    // int n = polygon.size()-1;

    // Point_2 p_pt2 = to_point_2(p);
    // for(auto c : polygon){
    //     Point_2 c_pt2 = to_point_2(c);
    //     if (norm(c_pt2 - p_pt2) == 0){
    //         return false;
    //     }
    // }

    // auto o1 = orientation(polygon[0], polygon[1], p);
    // if (o1==CGAL::COLLINEAR){
    //     return true;
    // }
    // for (int i=1; i<n+1; i++){
    //     int I = (i+1)%n;
    //     double o2  = orientation(polygon[i], polygon[I], p);
    //     if (o2==CGAL::COLLINEAR){
    //         return true;
    //     }
    //     if (o2 != o1){
    //         return false;
    //     }
    // }
    // return true;
}

bool do_intersect (std::vector<Point> polygon1, std::vector<Point> polygon2){
    // for (Point p : polygon2){
    //     if (is_inside(polygon1, p)){
    //         return true;
    //     }
    // }
    // return false;

    std::vector<Point> poly;

    int n1 = polygon1.size();
    int n2 = polygon2.size();
    for (int i=0; i<n1; i++){
        int I = (i+1)%n1;
        for (int j=0; j<n2; j++){
            int J = (j+1)%n2;
            
            if (norm(polygon1[i]-polygon2[j])==0 ||
                norm(polygon1[i]-polygon2[J])==0 ||
                norm(polygon1[I]-polygon2[j])==0 ||
                norm(polygon1[I]-polygon2[J])==0)
                continue;

            auto s1 = Segment_2(polygon1[i], polygon1[I]);
            auto s2 = Segment_2(polygon2[j], polygon2[J]);

            if (CGAL::do_intersect(s1, s2)){
                return true;
            }
        }
    }
    return false;
}

Point_2 guard_point(Bezier_curve P){
    
    /*
    Return: a guard point on the left of the curve
    */

    Point_2 splus;
    Point_2 sminus;

    if (!is_guardable(P, &splus, &sminus)){
        // std::cout << "The curve " << P << " is not guardable." << std::endl;
        Point_2 p2;
        return p2;
    }
    
    Point_2 P0 = P.front();
    Point_2 Pn = P.back();
    
    Point_2 X = meeting_point(P0, splus, Pn, sminus);

    auto direction = Pn - P0;
    direction  = normal(direction);
    if (cross(Pn-P0, direction) > 0)
        direction = - direction;

    direction = direction/norm(direction);
    // std::cout << "Direction: " << norm(direction) << std::endl;
    // return (X + 50 * direction);
    // std::cout << std::endl << "X: " << X.x() << ", " << X.y() << std::endl;
    // std::cout << "d: " << 10*direction.x() << ", " << 10*direction.y() << std::endl;
    return X + direction;
}

std::vector<Point> guard_envelope(Bezier_curve P){

    Point_2 Ol = guard_point(P);
    int n = P.size()-1;
    
    Point_2 P0 = P.front();
    Point_2 Pn = P.back();


    Point_2 Or;
    Or = P0+(Pn-Ol);
    
    std::vector<Point> guard;
    guard.push_back(to_point(P0));
    guard.push_back(to_point(Ol));
    guard.push_back(to_point(Pn));
    guard.push_back(to_point(Or));

    return guard;
}

std::vector<Bezier_curve> make_all_guardable(std::vector<Bezier_curve> curves)
{
    unsigned int            n_curves = curves.size();
    Bezier_curve            Btmp;
    Point_2                 splus;
    Point_2                 sminus;
    int                     k;

    std::cout << "Making all curves guardable... ";

    k = 0;
    while(k < n_curves){
        Btmp = curves[k];
        if (!is_guardable(Btmp, &splus, &sminus)){
            Point I;
            curves = bisect(curves, k, I);
            n_curves = curves.size();
        }
        else{
            k++;
        }
    }
    std::cout << "[OK]" << std::endl;

    return curves;
}

std::vector<Bezier_curve> make_no_guard_intersection(std::vector<Bezier_curve> curves, int count)
{
    int n_curves = curves.size();
    Bezier_curve bezier1;
    Bezier_curve bezier2;
    if (count>-1){
        std::cout << "Making no guard intersection... " << curves.size() << std::endl;
    }
    else if (count > 3){
        return curves;
    }
    std::vector<Point_2> corners = get_corners(curves);

    int i = 0;
    int j;
    while(i<n_curves && i < 30){
        j = 0;
        while(j<n_curves && j < 30){
            // std::cout << i << " " << j << " " << n_curves << std::endl;
            if (i==j){
                j++;
                continue;
            }

            bezier1 = curves[i];
            bezier2 = curves[j];

            int possible_corner = 0;
            
            auto guard1 = guard_envelope(bezier1);
            auto guard2 = guard_envelope(bezier2);

            if (!do_intersect(guard1, guard2)){
                // std::cout << "do not intersect at " << i << " " << j << std::endl;
                j++;
                continue;
            }
            // std::cout << "do intersect ";
            // Triangle T1 = Triangle(guard1[0], guard1[1], guard1[2]);
            // Triangle T2 = Triangle(guard2[0], guard2[1], guard2[2]);
            Triangle T1 = Triangle(guard1[0], guard1[1], guard1[2]);
            Triangle T2 = Triangle(guard2[0], guard2[1], guard2[2]);

            int I;
            if (abs(T1.area()) > abs(T2.area())){
                I = i;
            }
            else if(T1.area()==0 && T2.area()==0){
                j++;
                continue;
            }
            else{
                I = j;
            }
            // std::cout << "and cut " << I << std::endl;
            Point tmp;
            curves = bisect(curves, I, tmp);
            n_curves = curves.size();
        }
        i++;
    }
    std::cout << "[OK]" << curves.size() << std::endl;
    return curves;
}

std::vector<Bezier_curve> target_length(std::vector<Bezier_curve> curves)
{
    unsigned int            n_curves = curves.size();
    Bezier_curve            Btmp;
    Point_2                 splus;
    Point_2                 sminus;
    int                     k;

    std::cout << "Making all curves approx same length... ";

    double L = 1000;
    for (auto bezier : curves){
        auto d = get_distance(bezier);
        if (d < L){
            L = d;
        }
    }

    k = 0;
    while(k < n_curves){
        Btmp = curves[k];
        if (get_distance(Btmp) > 3*L/2){
            Point I;
            curves = bisect(curves, k, I);
            n_curves = curves.size();
        }
        else{
            k++;
        }
    }
    std::cout << "[OK]" << std::endl;

    return curves;
}