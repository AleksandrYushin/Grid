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

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), t(0.0), vx(0.0), vy(0.0), vz(0.0), ax(0.0), ay(0.0), az(0.0), alpha(0.0)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double t, double vx, double vy, double vz, double ax, double ay, double az, double alpha)
            : x(x), y(y), z(z), t(t), vx(vx), vy(vy), vz(vz), ax(ax), ay(ay), az(az), alpha(alpha)
    {
    }
    //get
    double get_t() { return t;}

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move_a(double tau) {
        ax = 0;
        ay = 0;
        az = 0;
    }
    void move_v(double tau) {
        vx += ax*tau;
        vy += ay*tau;
        vz += az*tau;
    }
    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }
    void move_t(double tau, double h, unsigned int i, unsigned int j, unsigned int k, vector<vector<vector<CalcNode>>> &points) {
        //Динамика определяется ур. Теплопроводности: D(T,t) - a*Laplas(T)=0
        // Соседние точки
        double tx_before = points[i-1][j][k].get_t();
        double tx_after = points[i+1][j][k].get_t();
        double ty_before = points[i][j-1][k].get_t();
        double ty_after = points[i][j+1][k].get_t();
        double tz_before = points[i][j][k-1].get_t();
        double tz_after = points[i][j][k+1].get_t();
        //Лаплас
        double Dx = ((tx_after - t)/h - (t - tx_before)/h)/h;
        double Dy = ((ty_after - t)/h - (t - ty_before)/h)/h;
        double Dz = ((tz_after - t)/h - (t - tz_before)/h)/h;
        //Итог
        t += alpha*(Dx+Dy+Dz)*tau;
    }
    void move_t_edge(double t_edge) {
        t = t_edge;
    }

};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<vector<vector<CalcNode>>> points;

public:
    // Конструктор сетки size x size точек с шагом h по пространству
    CalcMesh(unsigned int size, double h, double alpha0) {
        points.resize(size);
        for(unsigned int i = 0; i < size; i++) {
            points[i].resize(size);
            for(unsigned int j = 0; j < size; j++) {
                points[i][j].resize(size);
                for(unsigned int k = 0; k < size; k++) {
                    // Начальные координаты зададим равномерно в XYZ
                    double pointX = i * h + h/2;
                    double pointY = j * h + h/2;
                    double pointZ = k * h + h/2;
                    // Начальная температура ~300К
                    double t = 300.0;
                    // Неоднородная теплопроводность
                    double alpha = (pow(pointX, 0.5) + pow(pointY, 0.5) + pow(pointZ, 0.5)) * alpha0;
                    // Задавать скорости в данном примере поленимся, будут нули
                    points[i][j][k] = CalcNode(pointX, pointY, pointZ, t, 0., 0., 0., 0., 0., 0., alpha);
                }
            }
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau, double h, unsigned int step) {
        // По сути метод просто двигает все точки
        for(unsigned int i = 0; i < points.size(); i++) {
            for(unsigned int j = 0; j < points[i].size(); j++) {
                for(unsigned int k = 0; k < points[i][j].size(); k++) {
                    // Края
                    if (i==0 | j==0 | k==0 | i==(points.size()-1) | j==(points[i].size()-1) | k==(points[i][j].size()-1)) {
                        double x = double(k)/(points[i][j].size()-1);                        
                        points[i][j][k].move_t_edge(source_temp(x, step));
                    }
                    else {
                        points[i][j][k].move_t(tau, h, i, j, k, points);
                    }
                    points[i][j][k].move_a(tau);
                    points[i][j][k].move_v(tau);
                    points[i][j][k].move(tau);
                }
            }
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле на точках сетки
        auto t = vtkSmartPointer<vtkDoubleArray>::New();
        t->SetName("t");

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

                    // И значение скалярного поля тоже
                    t->InsertNextValue(points[i][j][k].t);
                }
            }
        }

        // Задаём размеры VTK-сетки (в точках, по трём осям)
        structuredGrid->SetDimensions(number, number, number);
        // Грузим точки в сетку
        structuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        structuredGrid->GetPointData()->AddArray(vel);
        structuredGrid->GetPointData()->AddArray(t);

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
    // Температуропроводность
    double alpha = 0.003;

    // Создаём сетку заданного размера
    CalcMesh mesh(size, h, alpha);

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
