
#include "variational.h"

Point to_point(Point_2 p){
    return Point(p.x(), p.y());
}

Point_2 to_point_2(Point p){
    return Point_2(p.x(), p.y());
}

double cos_of_angle(Point_2 a, Point_2 b){

    /*
    Return: cross product between two vectors
    if the result is negative: it is a clockwise rotation
    if the result is positive: it is a anti-clockwise rotation
    */
    double n_a = norm(a);
    double n_b = norm(b);

    a = a/n_a;
    b = b/n_b;

    return CGAL::to_double(a.x()*b.x()+a.y()*b.y());
}

double norm(Edge e){

    /*
    Return: length of an edge
    */
    auto p1 = e.first->vertex(e.first->cw(e.second))->point();      // get one endpoint
    auto p2 = e.first->vertex(e.first->ccw(e.second))->point();     // get the other endpoint
    
    return norm(to_point_2(p1) - to_point_2(p2));                                                 // return its norm
}

std::vector<Point_2> curve_axis(Bezier_curve P){
    /*
    Returs: d, splus, sminus
    
    If a curve is guardable, this function will return
         - the axis d of the curves
         - the steepest control vectors splus, and splus
         (anti-clockwise and clockwise with respect to d, respectively) 
    if at some point, we can not do the calculation anymore
    this mean the curve is not guardable and we return zero vectors.
    */

    unsigned int                n = P.size()-1;
    std::vector<Point_2>    dss;

    // initialise splus as the control vector s0 = P(1) - P(0)
    Point_2 p = P[1];
    Point_2 q = P[0];
    Point_2 splus = p-q;

    // initialise sminus as the control vector s1 = P(2) - P(1)
    p = P[2];
    q = P[1];
    Point_2 sminus = p-q;

    Point_2 si;
    Point_2 zero(0,0);

    // if the rotation splus to sminus is anti-clockwise
    // we swap the values of the two vectors

    if (cross(sminus, splus) < 0){
        Point_2 tmp = sminus;
        sminus = splus;
        splus = tmp;
    }

    for (int i=2; i<n; i++){

        // now take the i-th control vector si = P(i+1) - P(i)
        p = P[i+1];
        q = P[i];
        si = p-q;


        if (cross(si, splus) < 0){          // if the rotation si to splus is clockwise
            if (cross(si, sminus) < 0){     // if the rotation si to sminus is also clockwise
                splus = si;                 
            }
            else{                           // the curve is not guardable and return zeros
                dss.push_back(zero);
                dss.push_back(zero);
                dss.push_back(zero);
                return dss;
            }
        }
        else{                               // if the rotation si to splus is anti-clockwise
            if (cross(si, sminus) < 0){     // if the rotation si to sminus is clockwise
                ;                           // we do nothing as it is just an usual vector inside the two extremes
            }
            else{                           // if the rotation si to sminus is anti-clockwise
                sminus = si;
            }
        }
    }

    double nsp = norm(splus);
    splus = splus / nsp;
    double nsm = norm(sminus);
    sminus = sminus / nsm;

    // the axis d is take as the bisector of splus and sminus

    auto d = splus + sminus;
    dss.push_back(d);
    dss.push_back(splus);
    dss.push_back(sminus);
    return dss;
}

bool is_guardable(Bezier_curve P, Point_2* splus, Point_2* sminus){
    /*
    Return: true or false depending on the guardability of the curve

    if the axis d is the zero vector: it is not guardable
    */
    // std::cout << P.size() << std::endl;
    if (P.size()==2) return true;

    auto caxis = curve_axis(P);
    auto d = caxis[0];
    *splus = caxis[1];
    *sminus = caxis[2];
    if (norm(d) == 0){
        return false;
    }
    return true;
}

bool is_valid(Bezier_curve Bezier, Point_2 *splus, Point_2 *sminus){
    if (Bezier.size()==2) return true;

    if (!is_guardable(Bezier, splus, sminus)) return false;

    // auto dot = splus->x()*sminus->x() + splus->y() * sminus->y();

    // if (dot > 0.97) return true;
    Point_2 p0 = Bezier[0];
    Point_2 pn = Bezier.back();

    Point_2 p0n = pn-p0;

    // the limit angle is 0.361367124 rad
    // and we want the angle of the curve to be less than half of this
    // cos(0.361367124/2) = 0.98372108513
    // this is done so that when we have curve which are valid but have two curved edges,
    // we can just insert the centroid of the triangle to divide into 3,
    // such that each triangle has at most one curved edge, but still has the geometric map

    // if (cos_of_angle(*sminus, p0n) > 0.984 && cos_of_angle(p0n, *splus) > 0.984)return true;

    // no this previous is false, we just need the cos to be less than this so that the curve would be automatically guarded

    if (cos_of_angle(*sminus, p0n) > 0.93542 && cos_of_angle(p0n, *splus) > 0.93542)return true;
    return false;
}

// Rat_point_2 guard_point(Bezier_curve_2 P){
    
//     /*
//     Return: a guard point on the left of the curve
//     */

//     Rat_point_2 splus;
//     Rat_point_2 sminus;

//     if (!is_guardable(P, &splus, &sminus)){
//         std::cout << "The curve " << P << " is not guardable." << std::endl;
//         Rat_point_2 p2;
//         return p2;
//     }
    
//     Rat_point_2 P0 = P.control_point(0);
//     Rat_point_2 Pn = P.control_point(P.number_of_control_points() - 1);
    
//     // X is optimal guard but it does not induce the regular geometric map

//     Rat_point_2 X = meeting_point(P0, splus, Pn, sminus);
//     Rat_point_2 d = Rat_point_2(splus.x()+sminus.x(), splus.y()+sminus.y());
//     Rat_point_2 n = Rat_point_2(-d.y(), d.x());
//     n = Rat_point_2(n.x()/norm(n), n.y()/norm(n));
//     if (cross(d,n)<0) n = Rat_point_2(-n.x(), -n.y());

//     // we shift X following the normal to the axis d
    
//     double shift = 0.01;
//     return Rat_point_2(X.x()+shift*n.x(), X.y()+shift*n.y());
// }

// std::vector<Rat_point_2> guard_envelope(Bezier_curve_2 P){

//     Rat_point_2 Ol = guard_point(P);
//     auto n = P.number_of_control_points()-1;
    
//     auto P0 = P.control_point(0);
//     auto Pn = P.control_point(n);


//     Rat_point_2 Or;
//     Or = Rat_point_2(P0.x()-Ol.x()+Pn.x(), P0.y()-Ol.y()+Pn.y());

//     std::vector<Rat_point_2> guard;
//     guard.push_back(P0);
//     guard.push_back(Ol);
//     guard.push_back(Pn);
//     guard.push_back(Or);

//     return guard;
// }

std::vector<Bezier_curve> bisect(std::vector<Bezier_curve> curves, int i, Point &p, double t){

    Bezier_curve P1;
    Bezier_curve P2;

    Bezier_curve Bezier = curves[i];
    curves.erase(curves.begin()+i);

    subdivide_bezier_curve(Bezier, P1, P2, t=t);

    curves.insert(curves.begin()+i, P1);
    curves.push_back(P2);

    p = to_point(P2[0]);
    return curves;
}

std::vector<Bezier_curve> add_boundary(std::vector<Bezier_curve> curves){

    std::vector<Bezier_curve> boundary;

    Bezier_curve B;
    B.clear();
    B.push_back(Point_2(5,5));
    B.push_back(Point_2(5, 745));
    boundary.push_back(B);

    B.clear();
    B.push_back(Point_2(995, 5));
    B.push_back(Point_2(995, 745));
    boundary.push_back(B);
    
    B.clear();
    B.push_back(Point_2(5,745));
    B.push_back(Point_2(995, 745));
    boundary.push_back(B);

    B.clear();
    B.push_back(Point_2(5, 5));
    B.push_back(Point_2(995,5));
    boundary.push_back(B);

    // unsigned int            n_curves = boundary.size();

    // double L = 1000;
    // for (auto bezier : curves){
    //     auto d = get_distance(bezier);
    //     if (d < L){
    //         L = d;
    //     }
    // }

    // int k = 0;
    // while(k < n_curves){
    //     B = boundary[k];
    //     if (get_distance(B) > 3*L/2){
    //         Point I;
    //         boundary = bisect(boundary, k, I);
    //         n_curves = boundary.size();
    //     }
    //     else{
    //         k++;
    //     }
    // }

    for (auto bezier : boundary){
        curves.push_back(bezier);
    }

    return curves;
}

std::vector<Bezier_curve> make_all_valid(std::vector<Bezier_curve> curves)
{
    unsigned int            n_curves = curves.size();
    Bezier_curve            Btmp;
    Point_2                 splus;
    Point_2                 sminus;
    unsigned int            k;

    std::cout << "Making all curves valid... ";

    bool all_valid = false;
    while(!all_valid){
        all_valid = true;
        n_curves = curves.size();
        for (k=0; k< n_curves; k++){
            Btmp = curves[k];
            if (!is_valid(Btmp, &splus, &sminus)){
                all_valid = false;
                Point I;
                curves = bisect(curves, k, I);
                break;
            }
        }
    }
    std::cout << "[OK]" << std::endl;
    // std::cout << "Number of curves: " << curves.size() << std::endl;
    return curves;
}

std::vector<Bezier_curve> make_all_guard_valid(std::vector<Bezier_curve> curves)
{
    unsigned int            n_curves = curves.size();
    Bezier_curve            Btmp;
    Point_2                 splus;
    Point_2                 sminus;
    int                     k;
    bool is_valid;
    std::cout << "Making all curve guards contained... " << std::endl;
    std::cout << "n_curves in: " << curves.size() << std::endl;
    k = 0;
    while(k < n_curves){
        Btmp = curves[k];
        auto G = guard_envelope(Btmp);
        is_valid = true;

        // destroy triangle that has a vertex outside the boundary
        for (auto p : G){
            if (p.x() < 6 || p.x() > 994)
                is_valid = false;
            if (p.y() < 6 || p.y() > 744)
                is_valid = false;
        }

        // generate perfect sized triangle

        // double d = norm(Btmp.front()-Btmp.back());
        // double d1 = norm(Btmp.front()-to_point_2(G[1]));
        // double d2 = norm(Btmp.back()-to_point_2(G[1]));
        // if ( d1 > d || d2 > d){
        //     is_valid = false;
        // }

        // if ( d1 / d2 > 2 || d2 / d1 > 2){
        //     is_valid = false;
        // }
        
        if (!is_valid){
            Point I;
            curves = bisect(curves, k, I);
            n_curves = curves.size();
        }
        else{
            k++;
        }
    }
    std::cout << "[OK]" << std::endl;
    std::cout << "n_curves out: " << curves.size() << std::endl;
    return curves;
}

std::vector<Point_2> get_corners(std::vector<Bezier_curve> curves){
    int                         n_curves = curves.size();
    std::vector<Point_2>    corners;

    for(int i=0; i < n_curves; i++){
        for (int j=i+1; j < n_curves; j++){
            Bezier_curve bez1 = curves[i];
            Bezier_curve bez2 = curves[j];
            int n1 = bez1.size()-1;
            int n2 = bez2.size()-1;

            Point_2 corner;
            Vector_2 s1;
            Vector_2 s2;
            bool can_be_a_corner = false;

            Point_2 p1 = bez1[0];
            Point_2 q1 = bez1.back();

            Point_2 p2 = bez2[0];
            Point_2 q2 = bez2.back();

            if (p1 == p2){
                corner = p1;
                can_be_a_corner = true;

                s1 = bez1[1]-corner;
                s2 = bez2[1]-corner;
            }

            if (p1 == q2){
                corner = p1;
                can_be_a_corner = true;

                s1 = bez1[1]-corner;
                s2 = bez2[n2-1]-corner;
            }
            if (q1 == p2){
                corner = q1;
                can_be_a_corner = true;

                s1 = bez1[n1-1]-corner;
                s2 = bez2[1]-corner;
            }

            if (q1 == q2){
                corner = q1;
                can_be_a_corner = true;

                s1 = bez1[n1-1]-corner;
                s2 = bez2[n2-1]-corner;
            }

            if (can_be_a_corner){
                if (CGAL::angle(s1, s2) == CGAL::ACUTE && cos_of_angle(s2,s1) > 0.5)
                {
                    std::cout << "Corner" << std::endl;
                    bool is_new = true;
                    for (auto c : corners){
                        if (c==corner) is_new = false;
                    }
                    if (is_new) corners.push_back(corner);
                }
            }
        }
    }
    return corners;
}

bool is_a_corner(Point_2 p, std::vector<std::pair<Point_2, double>> corners){
    for (auto c : corners){
        if (norm(c.first-p)<0.05)
            return true;
    }
    return false;
}

bool is_a_corner(Point_2 p, std::vector<Point_2> corners){
    for (auto c : corners){
        if (norm(c-p)<0.05)
            return true;
    }
    return false;
}

bool is_curved(Edge e, std::vector<Bezier_curve> curves){
    auto p1 = e.first->vertex(e.first->cw(e.second))->point();
    auto p2 = e.first->vertex(e.first->ccw(e.second))->point();

    auto P1 = to_point_2(p1);
    auto P2 = to_point_2(p2);
    for (auto Bezier : curves){
        if (norm(P1-Bezier.front()) == 0 && norm(P2-Bezier.back()) == 0)
            return true;
        if (norm(P2-Bezier.front()) == 0 && norm(P1-Bezier.back()) == 0)
            return true;
    }
    return false;
}

CDT triangulation(std::vector<Bezier_curve> curves, bool variational){
    CDT             T;
    Vertex_handle   v1;
    Vertex_handle   v2;
    
    std::cout << "Starting triangulation... ";

    // Point_2 p0 = Point_2(400,400);
    // Point_2 p1 = Point_2(600,400);
    // Point_2 d0 = Point_2(1,1);
    // Point_2 d1 = Point_2(1,-1);

    // Point_2 x = meeting_point(p0,d0,p1,d1);

    // T.insert(to_point(p0));
    // T.insert(to_point(p1));
    // T.insert(to_point(x));
    // T.insert(to_point(p0+(p1-x)));
    for (auto Bezier : curves){
        
        v1 = T.insert(to_point(Bezier.front()));
        v2 = T.insert(to_point(Bezier.back()));
        
        T.insert_constraint(v1, v2);

        if (variational && Bezier.size()>2){
            auto G = guard_envelope(Bezier);
            for (int i=0; i<4; i++){
                auto I = (i+1)%4;
                v1 = T.insert(G[i]);
                v2 = T.insert(G[I]);
        
                T.insert_constraint(v1, v2);
            }
            // T.insert(G.begin(),G.end());
        }
    }
    std::cout << "[OK]" << std::endl;
    return T;
}

CDT refine_edge_length(CDT T, std::vector < Bezier_curve > curves){
    double L = 1000;
    for (auto bezier : curves){
        auto d = get_distance(bezier);
        if (d < L){
            L = d;
        }
    }

    double max_length = 0;
    Edge longest;
    Point mid;
    for (Edge e : T.finite_edges()){
        if (!T.is_constrained(e)){
            auto p1 = e.first->vertex(e.first->cw(e.second))->point();
            auto p2 = e.first->vertex(e.first->ccw(e.second))->point();

            auto P1 = to_point_2(p1);
            auto P2 = to_point_2(p2);

            if (norm(P1 - P2) > max_length){
                max_length = norm(P1 - P2);
                mid = to_point((P1 + P2)/2);

                longest = e;
            }
        }
    }

    while(max_length > 2*L){
        Vertex_handle v1 = longest.first->vertex(longest.first->cw(longest.second));
        Vertex_handle v2 = longest.first->vertex(longest.first->ccw(longest.second));

        Face_handle f;
        int index;

        T.is_edge(v1, v2, f, index);

        Point c = T.circumcenter(f);
        T.insert(c);

        // T.insert_in_edge(mid, f, index);

        max_length = 0;
        for (Edge e : T.finite_edges()){
            if (!T.is_constrained(e)){
                auto p1 = e.first->vertex(e.first->cw(e.second))->point();
                auto p2 = e.first->vertex(e.first->ccw(e.second))->point();

                auto P1 = to_point_2(p1);
                auto P2 = to_point_2(p2);

                if (norm(P1 - P2) > max_length){
                    max_length = norm(P1 - P2);
                    mid = to_point((P1 + P2)/2);

                    longest = e;
                }
            }
        }
    }
    return T;
}

CDT
refine_curve(CDT T, std::vector < Bezier_curve > *curves, std::vector<std::pair<Point_2, double>> corners){
    // std::cout << "refine_curve" << std::endl;
    // std::cout << "corner size " << corners.size() << std::endl;
    int n_curves = curves->size();
        
    int i = 0;
    while(i < n_curves){

        bool at_i_is_okay = false;
        auto bezier = (*curves)[i];
        
        auto p0 = bezier[0];
        auto pn = bezier.back();
        Vertex_handle vh0 = T.insert(to_point(p0));
        Vertex_handle vhn = T.insert(to_point(pn));
        auto c = CGAL::midpoint(to_point(p0), to_point(pn));

        for(auto ci : corners){
            if (norm(to_point_2(c)-ci.first) < ci.second){
                // std::cout << "here " << i<< std::endl;
                i++;
                at_i_is_okay = true;
                break;
            }
        }

        if (at_i_is_okay) continue;

        at_i_is_okay = true;

        Face_handle fh;
        int index;
        std::vector < Point > points;

        // get the third corner of the faces

        T.is_edge(vh0, vhn, fh, index);
        int i1 = fh->index(vh0);
        int i2 = fh->index(vhn);
        points.push_back(fh->vertex(3-i1-i2)->point());

        T.is_edge(vhn, vh0, fh, index);
        i1 = fh->index(vh0);
        i2 = fh->index(vhn);
        points.push_back(fh->vertex(3-i1-i2)->point());

        
        for (auto p : points){
            if (CGAL::has_larger_distance_to_point(c, to_point(p0), p)){
                at_i_is_okay = false;

                Point I;
                // *curves = bisect(*curves, i, I);
                
                Bezier_curve B = (*curves)[i];
                Point p0 = to_point(B[0]);
                Point pn = to_point(B.back());
                vh0 = T.insert(p0);
                vhn = T.insert(pn);
                
                Face_handle fh;
                int edge_index;
                T.is_edge(vh0, vhn, fh, edge_index);
                T.remove_constraint(fh, edge_index);

                *curves = bisect(*curves, i, I);

                Bezier_curve B1 = (*curves)[i];
                p0 = to_point(B1[0]);
                pn = to_point(B1.back());
                T.insert_constraint(p0, pn);
                Bezier_curve B2 = (*curves).back();
                p0 = to_point(B2[0]);
                pn = to_point(B2.back());
                T.insert_constraint(p0, pn);

                n_curves++;
                break;
            }
        }
        if (at_i_is_okay) i++;
    }
    // std::cout << " [ok]" << std::endl;
    return T;
}

CDT
refine_encroached(CDT T, std::vector < Bezier_curve > *curves, std::vector<std::pair<Point_2, double>> corners) {
    std::cout << "refine encroached...";
    int prev_n;
    int new_n;

    int frame = 0;
constraint_refine:

    // draw(*curves, T, frame);
    // frame++;

    // std::cout << (curves)-> size() << std::endl;
    prev_n = (curves)->size();
    T = refine_curve(T, curves, corners);
    new_n = (curves)->size();
    if (prev_n != new_n) goto constraint_refine;

    std::cout << "[OK]" << std::endl;
    return T;
}

bool is_encroached(Point p, std::vector < Bezier_curve > curves, int &index){

    for(int i=0; i< curves.size(); i++){
        auto bezier = curves[i];
        auto p0 = bezier[0];
        auto pn = bezier.back();

        auto I = CGAL::midpoint(to_point(p0), to_point(pn));
        if (CGAL::has_larger_distance_to_point(I, to_point(p0), p)){
            index = i;
            return true;
        }
    }

    return false;
}

bool is_valid_triangle(CDT T, Face_handle f){
    Point_2 p, q, r;
        
    p = to_point_2(f->vertex(0)->point());
    q = to_point_2(f->vertex(1)->point());
    r = to_point_2(f->vertex(2)->point());

    Point_2 pq = q-p;
    Point_2 qr = r-q;
    Point_2 rp = p-r;

    // double dpq = norm(pq);
    // double dqr = norm(qr);
    // double drp = norm(rp);

    // if (2*dpq < dqr ||  2*dqr < dpq) return false;
    // if (2*dpq < drp ||  2*drp < dpq) return false;
    // if (2*drp < dqr ||  2*dqr < drp) return false;

    
    // limit angle is arcsin(1/2sqrt(2)) = 0.361367124
    // cos(0.361367124) = 0.93541434666

    if (cos_of_angle(pq, -rp) > 0.936) return false;
    if (cos_of_angle(qr, -pq) > 0.936) return false;
    if (cos_of_angle(rp, -qr) > 0.936) return false;

    return true;
}

bool comparator(std::pair < Face_handle, double> f1, std::pair < Face_handle, double > f2){
    return f1.second < f2.second;
}

std::vector < std::pair <Face_handle, double>>
insert(std::vector < std::pair <Face_handle, double>> M, std::pair<Face_handle, double> tri){

    int start = 0;
    int end = M.size()-1;

research:
    int mid = int((start+end)/2);

    if (start==mid){
        M.insert(M.begin()+mid, tri);
        return M;
    }
    else if (comparator(M[mid], tri)) start = mid;
    else end = mid;

    goto research;
}
bool is_endpoint_farthest(Point_2 corner, Bezier_curve bezier, double *dist){
    *dist = 1000;
    double di;
    std::cout << "is_endpoint_farthest ";
    if (norm(corner-bezier.back())==0){
        reverse(bezier.begin(), bezier.end());
    }

    for(int i=0; i<bezier.size(); i++){
        di = norm(corner - bezier[i]);

        if (di > 0 && di < *dist){
            // std::cout << "FALSE" << std::endl;
            // return false;
            *dist = di;
        }
        // *dist = di;
    }
    std::cout << "TRUE" << std::endl;
    return true;
}

void
protecting_balls(std::vector < Bezier_curve > *curves, std::vector<Point_2> C, std::vector<std::pair<Point_2, double>> *corners){
    std::cout << "Generate protecting balls... ";
    int n_corner = C.size();

    for (Point_2 c : C) (*corners).push_back(std::pair<Point_2, double>(c, 1000));

    std::pair<Point_2, double> ci;
    double ri;

    int n_curves = curves->size();
    for (int i=0; i<n_corner; i++){
        ci = (*corners)[i];
        ri = ci.second;

        auto center = ci.first;
        double temp_ri;
        for(int j = 0; j<n_curves; j++){
            auto bezier = (*curves)[j];
            Point p;
            while(!is_endpoint_farthest(center, bezier, &temp_ri)){
                *curves = bisect(*curves, j, p);
                bezier = (*curves)[j];
            }
            if (temp_ri < ri) ri = temp_ri;
        }
        ((*corners)[i]).second = ri-1;
    }

    for(auto c : *corners){
        std::cout << "Radius = " << c.second << std::endl;
    }
    std::cout << "[OK]" << std::endl;
}

CDT
refine_bad_triangle(CDT T, std::vector < Bezier_curve > *curves, std::vector<std::pair<Point_2, double>> corners) {
    std::cout << "Refine bad triangles... ";

    std::vector<std::pair < Face_handle, double >> bad_faces;
    for (Face_handle f: T.finite_face_handles()){
        auto tri = T.triangle(f);
        if(!is_valid_triangle(T,f)){
            bad_faces.push_back(std::pair<Face_handle, double>(f, CGAL::to_double(tri.area())));
        }
    }
    std::sort(bad_faces.begin(), bad_faces.end(), comparator);
    int iter = 0;

    // auto corners = get_corners(*curves);
    // std::cout << std::endl;
    int frame = 1000;
restart:
    
    // draw(*curves, T, frame);
    // frame++;
    
    CGAL::Real_timer t;
    CGAL::Real_timer t_;
    // t_.reset();
    // t_.start();

    if (iter > 5000) return T;
    iter++;

    bool all_is_good = true;
    
    // std::cout << "\rRefine bad triangles..." << iter / 50.0<< "%";
    
    // auto bad = bad_faces[0];
    // bad_faces.erase(bad_faces.begin());

    auto bad = bad_faces.back();
    bad_faces.pop_back();

    Vertex_handle vh1 = bad.first->vertex(0);
    Vertex_handle vh2 = bad.first->vertex(1);
    Vertex_handle vh3 = bad.first->vertex(2);

    Face_handle f;
    auto is_face = T.is_face(vh1, vh2, vh3, f);

    if (is_face){
        Point c = T.circumcenter(f);
        
        for (auto ci : corners){
            if (norm(ci.first-to_point_2(c)) < ci.second){
                if (bad_faces.size()==0){
                    // std::cout << "\rRefine bad triangles... [OK]" << std::endl;
                    return T;
                }
                else
                    goto restart;
            }
        }
        int index;
        
        if (is_encroached(c, *curves, index)){
            Bezier_curve B = (*curves)[index];
            Point p0 = to_point(B[0]);
            Point pn = to_point(B.back());
            Vertex_handle vh0 = T.insert(p0);
            Vertex_handle vhn = T.insert(pn);
            
            Face_handle fh;
            int edge_index;
            T.is_edge(vh0, vhn, fh, edge_index);
            T.remove_constraint(fh, edge_index);

            *curves = bisect(*curves, index, c);

            Bezier_curve B1 = (*curves)[index];
            p0 = to_point(B1[0]);
            pn = to_point(B1.back());
            T.insert_constraint(p0, pn);
            Bezier_curve B2 = (*curves).back();
            p0 = to_point(B2[0]);
            pn = to_point(B2.back());
            T.insert_constraint(p0, pn);

            c = p0;
            
        }
        
        if (c.x() < 50 || c.x() > 950 || c.y() < 50 || c.y() > 950) ;
        else{
            Vertex_handle vhc = T.insert(c);
            auto f = T.incident_faces(vhc);
            auto new_f = f;
            do{
                if(T.is_infinite(new_f)){
                    new_f++;
                    continue;
                }
                Triangle tri = T.triangle(new_f);      
                if(!is_valid_triangle(T,new_f)){
                    bad_faces = insert(bad_faces, std::pair<Face_handle, double>(new_f, CGAL::to_double(tri.area())));
                }
                new_f++;
            }while(new_f != f);
        }
        
    }

    // t_.stop();
    // std::cout << "time t_ " << t_.time() << std::endl;

    if (bad_faces.size()==0){
        // std::cout << "\rRefine bad triangles...[OK]" << std::endl;
        std::cout << "[OK]" << std::endl;
        return T;
    }
    else
        goto restart;

}

CDT split_constrained_triangles(CDT T, std::vector<std::pair<Point_2, double>> corners){

    std::vector<Point> centroids;
    for (Face_handle f : T.finite_face_handles()){
        if (!is_valid_triangle(T,f)) continue;        

        Triangle tri = T.triangle(f);
        bool stop = false;
        bool two_constrained_edges = false;
        for (int i=0; i<3; i++){
            Vertex_handle vh = f->vertex(i);
            if (!T.are_there_incident_constraints(vh)){
                stop = true;
                break;
            }

            int n_constrained = 0;
            auto end = T.incident_edges(vh, f);
            auto eit = end;
            do{
                // if there is an endpoint of the edge that does not belong to the face
                // then we don't need that edge, as we want to stay in the face

                Point p = eit->first->vertex(eit->first->cw(eit->second))->point();      // get one endpoint
                if(!tri.has_on_boundary(p)){
                    eit++;
                    continue;
                }

                p = eit->first->vertex(eit->first->ccw(eit->second))->point();
                if(!tri.has_on_boundary(p)){
                    eit++;
                    continue;
                }
                
                if (T.is_constrained(*eit)) n_constrained++;
                eit++;

            }while(eit!=end);

            bool no_corner = true;
            for(int j=0; j<3; j++){
                Point_2 p = to_point_2(f->vertex(j)->point());
                if (is_a_corner(p, corners)){
                    no_corner = false;
                    break;
                }
            }
            if (n_constrained == 2 && no_corner){
                two_constrained_edges = true;
                break;
            }
        }
        if (stop) continue;
        
        if (two_constrained_edges){
            Triangle tri = T.triangle(f);

            centroids.push_back(CGAL::centroid(tri));
        }
    }

    T.insert(centroids.begin(), centroids.end());
    return T;
}

void
draw(std::vector<Bezier_curve> curves, CDT T, int frame = -1){

    unsigned int                            n_sample = 10;
    unsigned int                            n_curves;
    // unsigned int               n_points;
    std::vector<Point_2>  C;

    // QApplication app(argc, argv);
    QPixmap pixmap(1000, 750);
    pixmap.fill( Qt::white );

    QPainter painter( &pixmap );
    painter.setRenderHint(QPainter::Antialiasing);
    //painter.translate(500, 500); // center figure

    QPen penLines;  // creates a default pen
    penLines.setWidth(1);
    penLines.setBrush(Qt::black);

    QPen penBezier;  // creates a default pen
    penBezier.setWidth(2);
    penBezier.setBrush(Qt::red);

    QPen penPoints;  // creates a default pen
    penPoints.setWidth(1);
    penPoints.setBrush(Qt::black);
    penPoints.setCapStyle(Qt::RoundCap);

    QVector<QPoint> polyPoints;
    QPolygon guard_env;

    std::cout << "Rendering... ";
    for (Edge e : T.finite_edges()){

        if (!T.is_constrained(e)){
            polyPoints.clear();
            auto p1 = e.first->vertex(e.first->cw(e.second))->point();
            auto p2 = e.first->vertex(e.first->ccw(e.second))->point();
            polyPoints.push_back(QPoint(CGAL::to_double(p1.x()), CGAL::to_double(p1.y())));
            polyPoints.push_back(QPoint(CGAL::to_double(p2.x()), CGAL::to_double(p2.y())));
            painter.setPen( penLines );
            painter.drawPolyline(polyPoints);
        }
    }

    for(auto Bezier : curves){
        C.clear();
        polyPoints.clear();
        auto p1 = Bezier[0];
        auto p2 = Bezier.back();
        if (Bezier.size()==2){ 
            polyPoints.push_back(QPoint(CGAL::to_double(p1.x()), CGAL::to_double(p1.y())));
            polyPoints.push_back(QPoint(CGAL::to_double(p2.x()), CGAL::to_double(p2.y())));
        }
        else{
            int sample = int(norm(p1-p2));
            sample_bezier_curve(Bezier, 0,1,sample,C);
            for (auto& p : C) polyPoints << QPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
        }
        
        painter.setPen( penBezier );
        painter.drawPolyline(polyPoints);
    }

    std::cout << "[OK]" << std::endl;
  
    if (frame < 0)
        pixmap.save("render.png");
    else{
        std::stringstream a;
        a << frame;
        a << ".png";
        const std::string tmp = a.str();
        pixmap.save(tmp.c_str());
    }
    
    // QLabel myLabel;
    // myLabel.resize(1000,1000);
    // myLabel.setPixmap(pixmap);
    // myLabel.setWindowTitle("Bezier guarding viewer");
    // myLabel.show();
    // return app.exec();
}

int
draw(std::vector<Bezier_curve> curves, int argc, char *argv[]){

    unsigned int                            n_sample = 2000;
    unsigned int                            n_curves;
    // unsigned int               n_points;
    std::vector<Point_2>  C;

    QApplication app(argc, argv);
    QPixmap pixmap(1000, 1000);
    pixmap.fill( Qt::white );

    QPainter painter( &pixmap );
    //painter.translate(500, 500); // center figure

    QPen penLines;  // creates a default pen
    penLines.setWidth(1);
    penLines.setBrush(Qt::black);

    QPen penBezier;  // creates a default pen
    penBezier.setWidth(1);
    penBezier.setBrush(Qt::red);

    QPen penPoints;  // creates a default pen
    penPoints.setWidth(5);
    penPoints.setBrush(Qt::black);
    penPoints.setCapStyle(Qt::RoundCap);

    QVector<QPoint> polyPoints;
    QPolygon guard_env;

    std::cout << "Rendering... ";
    for(auto Bezier : curves){
        C.clear();
        sample_bezier_curve(Bezier, 0,1,n_sample,C);

        polyPoints.clear();
        for (auto& p : C) polyPoints << QPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
        
        painter.setPen( penBezier );
        painter.drawPolyline(polyPoints);

        // polyPoints.clear();
        // auto p = Bezier.control_point(0);
        // auto q = Bezier.control_point(Bezier.number_of_control_points()-1);
        
        // polyPoints.push_back(QPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y())));
        // polyPoints.push_back(QPoint(CGAL::to_double(q.x()), CGAL::to_double(q.y())));
        // painter.setPen( penPoints );
        // painter.drawPoints(polyPoints);
    }

    auto corners =  get_corners(curves);
    polyPoints.clear();
    for(auto c : corners){
        polyPoints.push_back(QPoint(CGAL::to_double(c.x()), CGAL::to_double(c.y())));
    }
    painter.setPen( penPoints );
    painter.drawPoints(polyPoints);

    std::cout << "[OK]" << std::endl;
  
    std::cout << "[OK]" << std::endl;
  
    std::stringstream a;
    a << "output/";
    auto input = argv[1];
    for(int i=0; i<std::strlen(input); i++){
        if (std::isdigit(input[i])) a << input[i];
    }
    a << "-map.png";
    const std::string tmp = a.str();
    pixmap.save(tmp.c_str());
    
    QLabel myLabel;
    myLabel.resize(1000,1000);
    // pixmap = pixmap.scaled(500, 500);
    myLabel.setPixmap(pixmap);
    
    myLabel.setWindowTitle("Bezier guarding viewer");
    myLabel.show();
    return app.exec();
}


// int main (int argc, char *argv[])
// {

//     if (argc < 2){
//         std::cout << "Error: no input file is given." << std::endl;
//         std::cout << "Usage: guarding filename." << std::endl;
//         return 0;
//     }
//     const char   *filename = argv[1];

//     std::ifstream   in_file (filename);

//     if (! in_file.is_open()) {
//         std::cerr << "Error: failed to open " << filename << std::endl;
//         return 0;
//     }

//     unsigned int                            n_curves;
//     std::vector<Bezier_curve>               curves;
//     Bezier_curve                            B;
//     CDT                                     T;

//     std::cout << "Loading input... ";
//     in_file >> n_curves;
//     for (int k = 0; k < n_curves; k++) {
//         B.clear();
        
//         int n;
//         in_file >> n;

//         for (int k = 0; k < n; k++){
//             int x;
//             int y;
//             in_file >> x;
//             in_file >> y;

//             B.push_back(Point_2(x,y));
//         }

//         curves.push_back (B);
//     }
//     std::cout << "[OK]" << std::endl;
    

//     curves = add_bundary(curves);
//     curves = make_all_valid(curves);
    
//     std::vector<std::pair<Point_2, double>> corners;
//     protecting_balls(curves, &corners);
    
//     T = triangulation(curves);
//     T = refine_encroached(T, &curves, corners);
//     T = refine_bad_triangle(T, &curves, corners);

//     // T = split_constrained_triangles(T);
//     std::cout << "Number of curves: " << curves.size() << std::endl;
//     std::cout << "Number of vertices: " << T.number_of_vertices() << std::endl;
//     std::cout << "Number of triangles: " << T.number_of_faces() << std::endl;
    
//     draw(curves, T, argc, argv);

//     // draw(curves, argc, argv);

//     return 1;
// }