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

#ifndef HELPER_H
#define HELPER_H

#include <QBrush>
#include <QFont>
#include <QPen>
#include <QWidget>

//#include <gtk/gtk.h>
//#include <gdk/gdkkeysyms.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

//#include <bits/stdc++.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include "variational.h"

//using namespace std;

//#define pi 3.14159265

#ifndef M
#define M 1000
#endif
//! [0]
class Helper
{
public:
    Helper();

public:
    std::vector<Bezier_curve> Curves;
	
    CDT T;
    std::vector<Point_2> data;
    int SelectedCurve = -1;
    int SelectedPoint = -1;    

    bool showGrid = false;
    bool snap = false;
	bool showPoints = true;
	
    bool newPoint = false;
    bool editPoint = false;
    bool straightedge = false;
    bool generated = false;
    bool showMesh = true;

    bool withBackground = false;

    QPixmap bg;
    
//    void paint(QPainter *painter, QPaintEvent *event, int elapsed);

    void clickLocation(QPainter *painter, QPoint point);
	void selectPoint(QPainter *painter, QPoint point);
	void insertPoint(QPainter *painter, QPoint point);
    void render(QPainter *painter);
    void render_as_png();

private:
    QBrush background;
    QBrush circleBrush;
    QFont textFont;
    QPen circlePen;
    QPen circleSelectedPen;
    QPen curvePen;
    QPen textPen;
};
//! [0]

#endif
