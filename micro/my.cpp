#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

using namespace std;

// Функция источников
double source_temp (double x, unsigned int step){
    double A = 300.0;
    double B = 300.0;
    
    if ((step % 50) < 30)
        A = 400.0;
    if ((step % 50) > 20)
        B = 400.0;
    return A*x + B*(1.0-x);
};

// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    // Температура
    double t;
    // Скорость
    double vx;
    double vy;
    double vz;
    // Ускорение
    double ax;
    double ay;
    double az;
    // Температуропроводность
    double alpha;
    // Давление
    double p;
    // Плотность (масса)
    double rho;

    // Флаг
    bool flag = 0;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), t(1.0), vx(0.0), vy(0.0), vz(0.0), ax(0.0), ay(0.0), az(0.0), alpha(1.0), p(0.0), rho(1.0)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double t, double vx, double vy, double vz, double ax, double ay, double az, double alpha, double p, double rho)
            : x(x), y(y), z(z), t(t), vx(vx), vy(vy), vz(vz), ax(ax), ay(ay), az(az), alpha(alpha), p(p), rho(rho)
    {
    }

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move_a(double tau, double h, unsigned int i, unsigned int j, unsigned int k, vector<vector<vector<CalcNode>>> &points) {
        // Давление смещает точки
        // Соседние точки
        double px_before = points[i-1][j][k].p;
        double px_after = points[i+1][j][k].p;
        double py_before = points[i][j-1][k].p;
        double py_after = points[i][j+1][k].p;
        double pz_before = points[i][j][k-1].p;
        double pz_after = points[i][j][k+1].p;
        // Закон смешения a(rho) = grad p
        ax = (px_before-px_after)/h /rho;
        ay = (py_before-py_after)/h /rho;
        az = (pz_before-pz_after)/h /rho;
    }
    void move_a_edge(double tau, double h, unsigned int i, unsigned int j, unsigned int k, vector<vector<vector<CalcNode>>> &points) {
        // Соседние точки
        double px_before, px_after, py_before, py_after, pz_before, pz_after = 0;
        if (i == 0) {px_before = 0;;}
        else {px_before = points[i-1][j][k].p;}
        if (i == points.size()-1) {px_after = 0;}
        else {px_after = points[i+1][j][k].p;}
        if (j == 0) {py_before = 0;}
        else {py_before = points[i][j-1][k].p;}
        if (j == points[i].size()-1) {py_after = 0;}
        else {py_after = points[i][j+1][k].p;}
        if (k == 0) {pz_before = 0;}
        else {pz_before = points[i][j][k-1].p;}
        if (k == points[i][j].size()-1) {pz_after = 0;}
        else {pz_after = points[i][j][k+1].p;}
        // Закон смешения a(rho) = grad p
        ax = (px_before-px_after)/h /rho;
        ay = (py_before-py_after)/h /rho;
        az = (pz_before-pz_after)/h /rho;
    }
    void move_v(double tau) {
        vx += ax*tau;
        vy += ay*tau;
        vz += az*tau;
    }
    void move_x(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }
    void move_t(double tau, double h, unsigned int i, unsigned int j, unsigned int k, vector<vector<vector<CalcNode>>> &points) {
        //Динамика определяется ур. Теплопроводности: D(T,t) - a*Laplas(T)=0
        // Соседние точки
        double tx_before = points[i-1][j][k].t;
        double tx_after = points[i+1][j][k].t;
        double ty_before = points[i][j-1][k].t;
        double ty_after = points[i][j+1][k].t;
        double tz_before = points[i][j][k-1].t;
        double tz_after = points[i][j][k+1].t;
        //Лаплас
        double Dx = ((tx_after - t)/h - (t - tx_before)/h)/h;
        double Dy = ((ty_after - t)/h - (t - ty_before)/h)/h;
        double Dz = ((tz_after - t)/h - (t - tz_before)/h)/h;
        //Итог
        t += alpha*(Dx+Dy+Dz)*tau;
    }
    void move_t_edge(double l, double step) {
        t = source_temp(l, step);
    }
    void move_alpha(double tau) {
        //Изменение под действием температуры
        alpha += 0.0;
    }
    void move_p(double tau, double h, unsigned int i, unsigned int j, unsigned int k, vector<vector<vector<CalcNode>>> &points) {
        // Сделаем след в закон Гука: tr(epsilon) = 3p/E - 6p(nu)/E
        // Итого p = (E/(1-2nu)/6h * (dx+dy+dz))
        double x_before = points[i-1][j][k].x;
        double x_after = points[i+1][j][k].x;
        double y_before = points[i][j-1][k].y;
        double y_after = points[i][j+1][k].y;
        double z_before = points[i][j][k-1].z;
        double z_after = points[i][j][k+1].z;
        p = 1 *(6*h - x_after + x_before - y_after + y_before - z_after + z_before) /6 /h;
    }
    void move_p_tp() {
        if (t > 370 & flag == 0){
            p += 0;
            flag = true;
        };
    }
    void move_p_edge() {
        p = 0;
    }

};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<vector<vector<CalcNode>>> points;
    vector<vector<vector<CalcNode>>> points_old;

public:
    // Конструктор сетки size x size точек с шагом h по пространству
    CalcMesh(unsigned int size, double h, double alpha0, double rho) {
        points.resize(size);
        points_old.resize(size);
        for(unsigned int i = 0; i < size; i++) {
            points[i].resize(size);
            points_old[i].resize(size);
            for(unsigned int j = 0; j < size; j++) {
                points[i][j].resize(size);
                points_old[i][j].resize(size);
                for(unsigned int k = 0; k < size; k++) {
                    // Начальные координаты зададим равномерно в XYZ
                    double pointX = i * h + h/2;
                    double pointY = j * h + h/2;
                    double pointZ = k * h + h/2;
                    // Начальная температура ~300К
                    double t = 300.0;
                    // Неоднородная теплопроводность
                    double alpha = (pow(pointX, 0.9) + pow(pointY, 0.9) + pow(pointZ, 0.9)) * alpha0;
                    // Задавать скорости в данном примере поленимся, будут нули
                    points[i][j][k] = CalcNode(pointX, pointY, pointZ, t, 0., 0., 0., 0., 0., 0., alpha, 0., rho);
                    points_old[i][j][k] = CalcNode(pointX, pointY, pointZ, t, 0., 0., 0., 0., 0., 0., alpha, 0., rho);
                }
            }
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau, double h, unsigned int step) {
        // Термодинамика
        for(unsigned int i = 0; i < points.size(); i++) {
            for(unsigned int j = 0; j < points[i].size(); j++) {
                for(unsigned int k = 0; k < points[i][j].size(); k++) {
                    if (i==0 | j==0 | k==0 | i==(points.size()-1) | j==(points[i].size()-1) | k==(points[i][j].size()-1)) {
                        // Края
                        double l = double(k)/(points[i][j].size()-1);
                        points[i][j][k].move_t_edge(l, step);
                    }
                    else {
                        points[i][j][k].move_t(tau, h, i, j, k, points_old);
                    }
                    points[i][j][k].move_alpha(tau);
                    points[i][j][k].move_p_tp();
                }
            }
        }
        // Давление
        for(unsigned int i = 0; i < points.size(); i++) {
            for(unsigned int j = 0; j < points[i].size(); j++) {
                for(unsigned int k = 0; k < points[i][j].size(); k++) {
                    if (i==0 | j==0 | k==0 | i==(points.size()-1) | j==(points[i].size()-1) | k==(points[i][j].size()-1)) {
                        // Края                       
                        points[i][j][k].move_p_edge();
                    }
                    else {
                        points[i][j][k].move_p(tau, h, i, j, k, points_old);
                    }
                }
            }
        }
        // Мехническое движение
        for(unsigned int i = 0; i < points.size(); i++) {
            for(unsigned int j = 0; j < points[i].size(); j++) {
                for(unsigned int k = 0; k < points[i][j].size(); k++) {
                    // Динамика
                    if (i==0 | j==0 | k==0 | i==(points.size()-1) | j==(points[i].size()-1) | k==(points[i][j].size()-1)) {
                        // Края                       
                        points[i][j][k].move_a_edge(tau, h, i, j, k, points_old);
                    }
                    else {
                        points[i][j][k].move_a(tau, h, i, j, k, points_old);
                    }
                    // Киниматика
                    points[i][j][k].move_v(tau);
                    points[i][j][k].move_x(tau);
                }
            }
        }
        // Перезапись
        points_old = points;

    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярные поля на точках сетки
        auto tem = vtkSmartPointer<vtkDoubleArray>::New();
        tem->SetName("tem");
        auto pres = vtkSmartPointer<vtkDoubleArray>::New();
        pres->SetName("pressure");

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        // Обходим все точки нашей расчётной сетки
        unsigned int number = (unsigned int)points.size();
        for(unsigned int i = 0; i < number; i++) {
            for(unsigned int j = 0; j < number; j++) {
                for(unsigned int k = 0; k < number; k++) {
                    // Вставляем новую точку в сетку VTK-снапшота
                    dumpPoints->InsertNextPoint(points[i][j][k].x, points[i][j][k].y, points[i][j][k].z);

                    // Добавляем значение векторного поля в этой точке
                    double _vel[3] = {points[i][j][k].vx, points[i][j][k].vy, points[i][j][k].vz};
                    vel->InsertNextTuple(_vel);

                    // И значение скалярных полей тоже
                    tem->InsertNextValue(points[i][j][k].t);
                    pres->InsertNextValue(points[i][j][k].p);
                }
            }
        }

        // Задаём размеры VTK-сетки (в точках, по трём осям)
        structuredGrid->SetDimensions(number, number, number);
        // Грузим точки в сетку
        structuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        structuredGrid->GetPointData()->AddArray(vel);
        structuredGrid->GetPointData()->AddArray(tem);
        structuredGrid->GetPointData()->AddArray(pres);

        // Создаём снапшот в файле с заданным именем
        string fileName = "my" + std::to_string(snap_number) + ".vts";
        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(structuredGrid);
        writer->Write();
    }
};

int main()
{
    // Размер расчётной сетки, точек на сторону
    unsigned int size = 30;
    // Шаг точек по пространству
    double h = 0.01;
    // Шаг по времени
    double tau = 0.01;
    // Температуропроводность, плотность
    double alpha = 0.003;
    double rho = 1;

    // Создаём сетку заданного размера
    CalcMesh mesh(size, h, alpha, rho);

    // Пишем её начальное состояние в VTK
    mesh.snapshot(0);

    // Дальше можно сделать какие-нибудь шаги по времени аналогично 2D-примеру.
    // на каждом шаге считаем новое состояние и пишем его в VTK
    for(unsigned int step = 1; step < 1000; step++) {
        mesh.doTimeStep(tau, h, step);
        mesh.snapshot(step);
    }

    return 0;
}
