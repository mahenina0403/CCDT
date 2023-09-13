#include "mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
//    QTimer *timer = new QTimer(this);
//    connect(timer, &QTimer::timeout, this, &MainWindow::update);
    //    connect(timer, &QTimer::timeout, openGL, &GLWidget::animate);
    //    timer->start(50);
    connect(ui->showGrid, &QCheckBox::stateChanged, ui->openGL, &GLWidget::gridToggle);
    connect(ui->snapOnGrid, &QCheckBox::stateChanged, ui->openGL, &GLWidget::snapToggle);
	connect(ui->showPoints, &QCheckBox::stateChanged, ui->openGL, &GLWidget::showPointsToggle);
	connect(ui->straightEdges, &QCheckBox::stateChanged, ui->openGL, &GLWidget::straightedgeToggle);
    connect(ui->showMesh, &QCheckBox::stateChanged, ui->openGL, &GLWidget::showMeshToggle);
    connect(ui->showBg, &QCheckBox::stateChanged, ui->openGL, &GLWidget::showBgToggle);

    connect(ui->drawMode, &QAbstractButton::clicked, ui->openGL, &GLWidget::drawmode);
    connect(ui->editMode, &QAbstractButton::clicked, ui->openGL, &GLWidget::editmode);

    connect(ui->GenerateButton, &QAbstractButton::clicked, ui->openGL, &GLWidget::generate);
    

//    timer->start(50);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionNew_curve_triggered(){

    ui->openGL->helper.Curves.clear();
    ui->openGL->helper.T.clear();

    ui->openGL->helper.SelectedCurve = -1;
    ui->openGL->helper.SelectedPoint = -1;

	ui->openGL->helper.newPoint = false;
	ui->openGL->helper.editPoint = false;
	ui->openGL->helper.generated = false;

    ui->openGL->repaint();
}

void MainWindow::on_actionSave_triggered(){

    std::ofstream tmpfile;
    tmpfile.open ("save.dat");
    
    ui->openGL->cleanData();

    auto C = ui->openGL->helper.Curves;

    tmpfile << C.size() << std::endl;
    for(auto bezier : C){
        if (bezier.size()>1){
            tmpfile << bezier.size() << " ";
            for (int i=0; i<bezier.size(); i++){
                tmpfile << bezier[i].x() << " ";
                tmpfile << bezier[i].y() << " ";
            }
            tmpfile << std::endl;
        }
    }
    tmpfile.close();
}

void MainWindow::on_actionRender_triggered(){
    
    ui->openGL->render();
}

void MainWindow::on_actionLoad_triggered(){
    auto current_folder = QString::fromUtf8((std::filesystem::current_path().string()).c_str());

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                current_folder,
                                                tr("text (*.dat *.txt)"));

    std::ifstream   in_file (fileName.toUtf8().constData());
    if (in_file.is_open()){
        std::cout << "File opened successfully from " << current_folder.toUtf8().constData() << std::endl;
        
        ui->openGL->helper.Curves.clear();

        int     n_curves;
        in_file >> n_curves;

        Bezier_curve    B;
        for (int k = 0; k < n_curves; k++) {
            B.clear();
        
            int n;
            in_file >> n;

            for (int k = 0; k < n; k++){
                int x;
                int y;
                in_file >> x;
                in_file >> y;

                B.push_back(Point_2(x,y));
            }

        ui->openGL->helper.Curves.push_back(B);
        }

    }
    ui->openGL->cleanData();
    
    ui->openGL->repaint();
}

void MainWindow::on_actionInsert_image_triggered(){

    auto current_folder = QString::fromUtf8((std::filesystem::current_path().string()).c_str());

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                current_folder,
                                                tr("Images (*.png *.jpg)"));

    std::ifstream   in_file (fileName.toUtf8().constData());
    if (in_file.is_open()){
        std::cout << "File opened successfully from " << current_folder.toUtf8().constData() << std::endl;
        ui->showBg->setChecked(true);
        ui->openGL->helper.withBackground = true;
        ui->openGL->helper.bg = QPixmap(fileName);
    }
    
    ui->openGL->repaint();
}