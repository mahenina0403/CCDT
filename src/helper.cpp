/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "helper.h"

#include <QPainter>
#include <QPaintEvent>
#include <QWidget>

//! [0]
Helper::Helper()
{
    background = QBrush(Qt::white);

    circlePen = QPen(Qt::black);
    circlePen.setWidth(7);
    circlePen.setCapStyle(Qt::RoundCap);

    circleSelectedPen = QPen(QColor(50, 102, 168));
    circleSelectedPen.setWidth(10);
    circleSelectedPen.setCapStyle(Qt::RoundCap);

    curvePen = QPen(Qt::black);
    curvePen.setWidth(1);
    curvePen.setCapStyle(Qt::RoundCap);

    textPen = QPen(Qt::black);
    textFont.setPixelSize(50);
}

void Helper::render(QPainter *painter)
{
    std::stringstream ss;

    painter->fillRect(QRectF(0, 0, 1000, 750), background);

    if (withBackground){
        painter->drawPixmap(0,0,bg);
    }

    
    if (showGrid){
        int x = 25;
        while(x<1000){
            if (x==500){
                curvePen.setColor(Qt::green);
            }
            else
                curvePen.setColor(QColor(150,150,150));
            painter->setPen(curvePen);
            painter->drawLine(x, 25, x, 725);
            x = x + 25;
        }
        x = 25;
        while(x<750){
            if (x==375){
                curvePen.setColor(Qt::red);
            }
            else
                curvePen.setColor(QColor(150,150,150));
            painter->setPen(curvePen);
            painter->drawLine(25, x, 975, x);
            x = x + 25;
        }
    }

    if (Curves.empty()){
        // painter->restore();
        return;
    }

    Bezier_curve C;
    double t = 0;
    // double dt = 1 / (double) M;     // equisdistant sample points in [0,1]
    double dt;
    double N;
    for (int i=0; i<Curves.size();i++){
        auto bezier = Curves[i];
        N = (double) norm(bezier.front()-bezier.back());
        dt = 1 / N;
        if (bezier.empty() || bezier.size()<2){
            break;
        }
        C.clear();
        if (generated && straightedge){
            C.push_back(bezier.front());
            C.push_back(bezier.back());
        }
        else if (bezier.size()==2){
            C.push_back(bezier.front());
            C.push_back(bezier.back());
        }
        else{
            t = 0;
            for (int j=0; j<=N; j++, t+=dt){
                C.push_back(point_on_bezier_curve(bezier,t));
            }
        }

        if (generated && (bezier.front().x()!=5 && bezier.front().x()!=995 && bezier.front().y()!=5 && bezier.front().y()!=745))
            curvePen.setColor(Qt::red);
        else if (SelectedCurve == i)
            curvePen.setColor(QColor(50, 102, 168));
        else
            curvePen.setColor(Qt::black);

        painter->setPen(curvePen);
        for (int j=1; j<C.size(); j++) painter->drawLine(C[j-1].x(), C[j-1].y(), C[j].x(), C[j].y());

    }
    
	if (!generated && showPoints){
        for (int i=0; i<Curves.size();i++){
            auto bezier = Curves[i];
            if (bezier.empty()){
                break;
            }
            if (SelectedCurve==i){
                for(int j=0; j<bezier.size(); j++){
                    if (SelectedPoint == j) painter->setPen(circleSelectedPen);
                    else painter->setPen(circlePen);

                    painter->drawPoint(bezier[j].x(), bezier[j].y());
                }
                curvePen.setColor(Qt::black);
                painter->setPen(curvePen);
                
                if(bezier.size()>2)
                    painter->drawLine(bezier[0].x(), bezier[0].y(), bezier[1].x(), bezier[1].y());
                if (bezier.size()==4)
                    painter->drawLine(bezier[2].x(), bezier[2].y(), bezier[3].x(), bezier[3].y());
            }
        }
	}

    if(generated && showMesh){
        QVector<QPoint> polyPoints;
        for (Edge e : T.finite_edges()){
            if (!is_curved(e, Curves)){
            // if (!T.is_constrained(e)){
                polyPoints.clear();
                auto p1 = e.first->vertex(e.first->cw(e.second))->point();
                auto p2 = e.first->vertex(e.first->ccw(e.second))->point();
                polyPoints.push_back(QPoint(CGAL::to_double(p1.x()), CGAL::to_double(p1.y())));
                polyPoints.push_back(QPoint(CGAL::to_double(p2.x()), CGAL::to_double(p2.y())));

                curvePen.setColor(Qt::black);
                painter->setPen( curvePen );
                painter->drawPolyline(polyPoints);
        }
    }
    }
    // painter->restore();
}

void Helper::render_as_png()
{
    
    QPixmap pixmap(1000, 750);
    pixmap.fill( Qt::white );

    QPainter painter( &pixmap );
    // painter.fillRect(QRectF(0, 0, 1000, 750), background);

    curvePen.setColor(QColor(150,150,150));
    painter.setPen(curvePen);
    
    if (Curves.empty()){
        // painter->restore();
        return;
    }

    Bezier_curve C;
    double t = 0;
    double dt = 1 / (double) M;     // equisdistant sample points in [0,1]
    double N;
    for (int i=0; i<Curves.size();i++){
        auto bezier = Curves[i];
        N = (double) norm(bezier.front()-bezier.back());
        dt = 1 / N;

        if (bezier.empty() || bezier.size()<2){
            break;
        }
        C.clear();
        if (generated && straightedge){
            C.push_back(bezier.front());
            C.push_back(bezier.back());
        }
        else if (bezier.size()==2){
            C.push_back(bezier.front());
            C.push_back(bezier.back());
        }
        else{
            t = 0;
            for (int j=0; j<=N; j++, t+=dt){
                C.push_back(point_on_bezier_curve(bezier,t));
            }
        }

        if (generated && (bezier.front().x()!=5 && bezier.front().x()!=995 && bezier.front().y()!=5 && bezier.front().y()!=745))
            curvePen.setColor(Qt::red);
        else if (SelectedCurve == i)
            curvePen.setColor(QColor(50, 102, 168));
        else
            curvePen.setColor(Qt::black);

        painter.setPen(curvePen);
        for (int j=1; j<C.size(); j++) painter.drawLine(C[j-1].x(), C[j-1].y(), C[j].x(), C[j].y());

    }
    
    if(generated && showMesh){
        QVector<QPoint> polyPoints;
        for (Edge e : T.finite_edges()){

            if (!T.is_constrained(e)){
                polyPoints.clear();
                auto p1 = e.first->vertex(e.first->cw(e.second))->point();
                auto p2 = e.first->vertex(e.first->ccw(e.second))->point();
                polyPoints.push_back(QPoint(CGAL::to_double(p1.x()), CGAL::to_double(p1.y())));
                polyPoints.push_back(QPoint(CGAL::to_double(p2.x()), CGAL::to_double(p2.y())));

                curvePen.setColor(Qt::black);
                painter.setPen( curvePen );
                painter.drawPolyline(polyPoints);
            }
        }
    }

    pixmap.save("render.png");
    
    // painter->restore();
}


void Helper::clickLocation(QPainter *painter, QPoint point)
{
    Point_2 tmp = Point_2(point.x(), point.y());
    Point_2 tmp2 = tmp;

    this->selectPoint(painter, point);

    if (SelectedCurve > -1 && SelectedPoint > -1 && editPoint){
        if (snap){
            int X = int(point.x()/25);
            int Y = int(point.y()/25);
            double dmin = norm(tmp-Point_2(25*(X),25*(Y)));
            tmp2 = Point_2(25*(X),25*(Y));
            for (int i=0; i<2; i++){
                for (int j=0; j<2; j++){
                    double d = norm(tmp-Point_2(25*(X+i),25*(Y+j)));
                    if (d < dmin){
                        tmp2 = Point_2(25*(X+i),25*(Y+j));
                        dmin = d;
                    }
                }
            }
        }

        Curves[SelectedCurve][SelectedPoint] = tmp2;
    }

}

void Helper::insertPoint(QPainter *painter, QPoint point)
{
    Point_2 tmp = Point_2(point.x(), point.y());
    Point_2 tmp2 = tmp;

    this->selectPoint(painter, point);
    if (SelectedCurve == -1 && SelectedPoint == -1 && newPoint){
        if (snap){
            int X = int(point.x()/25);
            int Y = int(point.y()/25);
            double dmin = norm(tmp-Point_2(25*(X),25*(Y)));
            tmp2 = Point_2(25*(X),25*(Y));
            for (int i=0; i<2; i++){
                for (int j=0; j<2; j++){
                    double d = norm(tmp-Point_2(25*(X+i),25*(Y+j)));
                    if (d < dmin){
                        tmp2 = Point_2(25*(X+i),25*(Y+j));
                        dmin = d;
                    }
                }
            }
        }

        Bezier_curve bezier;
        if (Curves.empty()){
            Curves.push_back(bezier);
        }

        bezier = Curves.back();
        if (bezier.size()==4){
            Bezier_curve newBezier;
            Curves.push_back(newBezier);
            SelectedPoint = 0;
        }
        else{
            SelectedPoint = bezier.size()-1;
        }
        SelectedCurve = Curves.size()-1;
        Curves[SelectedCurve].push_back(tmp2);
    }
}

void Helper::selectPoint(QPainter *painter, QPoint point)
{
    Point_2 tmp = Point_2(point.x(), point.y());
    SelectedCurve = -1;
    SelectedPoint = -1;

    for(int i=0; i<Curves.size(); i++){
        auto bezier = Curves[i];
        for(int j=0; j<bezier.size(); j++){
            if (norm(tmp-bezier[j]) < 20){
                SelectedCurve = i;
                SelectedPoint = j;
            }
        }
    }
}