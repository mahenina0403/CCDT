#ifndef MAINWINDOW_H
#define MAINWINDOW_H



#include <QMainWindow>
#include <QTimer>
#include <QFileDialog>

#include "../ui_mainwindow.h"
#include "glwidget.h"

#include <filesystem>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
//    QTimer      timer;
public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

public slots:
//    void update();
    void on_actionNew_curve_triggered();
    void on_actionSave_triggered();
    void on_actionRender_triggered();
    void on_actionInsert_image_triggered();
    void on_actionLoad_triggered();
    
private:
    Ui::MainWindow *ui;

};


#endif // MAINWINDOW_H
