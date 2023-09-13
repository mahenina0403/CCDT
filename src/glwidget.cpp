#include "glwidget.h"

#include <QPainter>
#include <QTimer>
#include <QEvent>
#include <QKeyEvent>

//! [0]
GLWidget::GLWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    elapsed = 0;
//    setFixedSize(1080, 720);
    setAutoFillBackground(false);
    setFocusPolicy(Qt::StrongFocus);
}
//! [0]

//! [1]
void GLWidget::animate()
{
    elapsed = (elapsed + qobject_cast<QTimer*>(sender())->interval()) % 1000;
    update();
}

void GLWidget::gridToggle(){
    helper.showGrid = !helper.showGrid;
    repaint();
}


void GLWidget::showPointsToggle(){
    helper.showPoints = !helper.showPoints;
    repaint();
}

void GLWidget::render(){
    helper.render_as_png();
}

void GLWidget::snapToggle(){
    helper.snap = !helper.snap;
    repaint();
}

void GLWidget::straightedgeToggle(){
    helper.straightedge = !helper.straightedge;
    repaint();
}

void GLWidget::showMeshToggle(){
    helper.showMesh = !helper.showMesh;
    repaint();
}
void GLWidget::drawmode(){
    helper.newPoint = true;
    helper.editPoint = false;
}

void GLWidget::editmode(){
    helper.newPoint = false;
    helper.editPoint = true;
}

void GLWidget::showBgToggle()
{
    helper.withBackground = !helper.withBackground;
    repaint();
}

void GLWidget::cleanData(){
    int n_curves = helper.Curves.size();
    int i = 0;
    while(i<n_curves){
        if(helper.Curves[i].size()<2){
            helper.Curves.erase(helper.Curves.begin()+i);
            n_curves--;
        }
        else
            i++;
    }
}

void GLWidget::generate(){
    if((helper.Curves.back()).empty()){
        helper.Curves.pop_back();
    }

    this->cleanData();

    std::vector<std::pair<Point_2, double>> balls;
    std::vector<Point_2> corners = get_corners(helper.Curves);
    std::cout << corners.size() << std::endl;
    
    // helper.Curves = make_all_guardable(helper.Curves);
    // helper.Curves = make_all_guard_valid(helper.Curves);
    // helper.Curves = make_no_guard_intersection(helper.Curves);
    
    helper.Curves = make_all_valid(helper.Curves);
    std::cout << helper.Curves.size() << std::endl;
    protecting_balls(&helper.Curves, corners, &balls);
    std::cout << helper.Curves.size() << std::endl;
    
    helper.Curves = target_length(helper.Curves);
    helper.Curves = add_boundary(helper.Curves);
    
    
    helper.T = triangulation(helper.Curves, false);

    // helper.T = refine_edge_length(helper.T, helper.Curves);

    helper.T = refine_encroached(helper.T, &helper.Curves, balls);
    helper.T = refine_bad_triangle(helper.T, &helper.Curves, balls);
    // helper.T = split_constrained_triangles(helper.T, balls);

    helper.generated = true;
    helper.newPoint = false;
    helper.editPoint = false;

    helper.showGrid = false;
    update();
}

//! [1]

//! [2]
void GLWidget::paintEvent(QPaintEvent *event)
{
    painter.begin(this);
    painter.setRenderHint(QPainter::Antialiasing);
    helper.render(&painter);
    painter.end();

}
//! [2]


void GLWidget::mousePressEvent(QMouseEvent *event) {
    m_mouse_click = event->pos();
    if (event->button() == Qt::LeftButton)
    {
        painter.begin(this);
        painter.setRenderHint(QPainter::Antialiasing);
        helper.insertPoint(&painter, m_mouse_click);
        helper.render(&painter);
        painter.end();
    }
	else if (event->button() == Qt::RightButton)
    {
        painter.begin(this);
        painter.setRenderHint(QPainter::Antialiasing);
        helper.selectPoint(&painter, m_mouse_click);
        helper.render(&painter);
        painter.end();
    }
    update();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event) {
    m_mouse_move = event->pos();
    if (event->button() == Qt::NoButton)
    {
        if (m_mouse_move != m_mouse_click){
            painter.begin(this);
            painter.setRenderHint(QPainter::Antialiasing);
            helper.clickLocation(&painter, m_mouse_move);
			helper.render(&painter);
            painter.end();
        }
    }

    m_mouse_click = m_mouse_move;
    update();
}

void GLWidget::keyPressEvent(QKeyEvent *event){
    
    if(event->key() == Qt::Key_Return){
        Bezier_curve newBezier;
        helper.Curves.push_back(newBezier);
        helper.SelectedPoint = 0;
        helper.SelectedCurve = helper.Curves.size()-1;

        repaint();
    }
    else if(event->key() == Qt::Key_Backspace || event->key() == Qt::Key_Delete){
        helper.Curves.erase(helper.Curves.begin()+helper.SelectedCurve);
        helper.SelectedCurve = -1;
        helper.SelectedPoint = -1;
    }
    update();
}