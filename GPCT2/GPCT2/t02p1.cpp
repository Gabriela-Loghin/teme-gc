
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "glut.h"
#include <vector>
#include <cmath>
#include <iostream>


// dimensiunea ferestrei in pixeli
#define dim 300

unsigned char prevKey;

// concoida lui Nicomede (concoida dreptei)
// $x = a + b \cdot cos(t), y = a \cdot tg(t) + b \cdot sin(t)$. sau
// $x = a - b \cdot cos(t), y = a \cdot tg(t) - b \cdot sin(t)$. unde
// $t \in (-\pi / 2, \pi / 2)$
void Display1() {
    double xmax, ymax, xmin, ymin;
    double a = 1, b = 2;
    double pi = 4 * atan(1.0);
    double ratia = 0.05;
    double t;

    // calculul valorilor maxime/minime ptr. x si y
    // aceste valori vor fi folosite ulterior la scalare
    xmax = a - b - 1;
    xmin = a + b + 1;
    ymax = ymin = 0;
    for (t = -pi / 2 + ratia; t < pi / 2; t += ratia) {
        double x1, y1, x2, y2;
        x1 = a + b * cos(t);
        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        x2 = a - b * cos(t);
        xmax = (xmax < x2) ? x2 : xmax;
        xmin = (xmin > x2) ? x2 : xmin;

        y1 = a * tan(t) + b * sin(t);
        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;

        y2 = a * tan(t) - b * sin(t);
        ymax = (ymax < y2) ? y2 : ymax;
        ymin = (ymin > y2) ? y2 : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    // afisarea punctelor propriu-zise precedata de scalare
    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (t = -pi / 2 + ratia; t < pi / 2; t += ratia) {
        double x1, y1, x2, y2;
        x1 = (a + b * cos(t)) / xmax;
        x2 = (a - b * cos(t)) / xmax;
        y1 = (a * tan(t) + b * sin(t)) / ymax;
        y2 = (a * tan(t) - b * sin(t)) / ymax;

        glVertex2f(x1, y1);
    }
    glEnd();

    glBegin(GL_LINE_STRIP);
    for (t = -pi / 2 + ratia; t < pi / 2; t += ratia) {
        double x1, y1, x2, y2;
        x1 = (a + b * cos(t)) / xmax;
        x2 = (a - b * cos(t)) / xmax;
        y1 = (a * tan(t) + b * sin(t)) / ymax;
        y2 = (a * tan(t) - b * sin(t)) / ymax;

        glVertex2f(x2, y2);
    }
    glEnd();
}

// graficul functiei 
// $f(x) = \bar sin(x) \bar \cdot e^{-sin(x)}, x \in \langle 0, 8 \cdot \pi \rangle$, 
void Display2() {
    double pi = 4 * atan(1.0);
    double xmax = 8 * pi;
    double ymax = exp(1.1);
    double ratia = 0.05;

    // afisarea punctelor propriu-zise precedata de scalare
    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (double x = 0; x < xmax; x += ratia) {
        double x1, y1;
        x1 = x / xmax;
        y1 = (fabs(sin(x)) * exp(-sin(x))) / ymax;

        glVertex2f(x1, y1);
    }
    glEnd();
}

class Punct {
public:
    double x, y;

    Punct(double x, double y) {
        this->x = x;
        this->y = y;
    }
};
//functia f(x)= d(x)/x , 0<x<=100;
//            =1       , x=0;
double Distanta(double x) {
    int i;
    double d1;
    double d2;
    double min;
    for (i = 0; i <= 100; i++) {
        if ((i <= x) and ((i + 1) >= x) and (x - i <= 1) and (x - i >= 0) and ((i + 1) - x) <= 1 and ((i + 1) - x) >= 0) {
            d1 = x - i;
            d2 = (i + 1) - x;
        }
    }
    if (d1 < d2) {
        min = d1;
    }
    else {
        min = d2;
    }

    return min;
}
std::vector<double> get_maximum(double x) {

    double y;
    double ymax = 1;
    double xmax = 100;
    double ymin = 0;
    double xmin = 0;
    std::vector<double> maximum;

    if (x == 0) {
        y = 1;
    }
    else {
        y = Distanta(x) / x;
    }
    xmax = (xmax < x) ? x : xmax;

    ymax = (ymax < y) ? y : ymax;
    ymin = (ymin > y) ? y : ymin;

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    maximum.push_back(xmax);
    maximum.push_back(ymax);

    return maximum;
}

std::vector<Punct> get_set_vertex() {
    double x, xv, yv;
    double ymax = 10;
    double ymin = 0;
    double xmin = 0;
    double xmax = 10;
    double ratia = 0.2;
    double f;
    std::vector<double> maximum;
    std::vector<Punct> puncte;


    for (x = ratia; x <= 75; x = x + ratia) {
            f = Distanta(x) / x;
            maximum = get_maximum(x);
            xv = x / maximum[0];
            yv = f / maximum[1];
            Punct new_point = Punct(xv, yv);
            puncte.push_back(new_point);
    }
    return puncte;
}

void Display3() {

    std::vector<Punct> puncte = get_set_vertex();
    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (Punct p : puncte) {
        glVertex2f(p.x, p.y);
        //glVertex2f(p.x, 0);
    }
    glEnd();
}

//melcul lui Pascal-concoida cercului :
//  x=2*(a*cost+b)*cost,y=2*(a*cost+b)*sint,t->(-pi,pi)

std::vector<double> get_min_max(double a, double b, double ratia, double t) {
    double pi = 4 * atan(1.0);
    double xmax;
    double ymax;
    double ymin, xmin;
    xmax = 2 * a + 2 * b;
    ymax = pi * (a + b);
    ymin = -b * pi;
    xmin = 2 * b;
    std::vector<double> max;

    double x, y;
    x = 2 * (a * cos(t) + b) * cos(t);
    y = 2 * (a * cos(t) + b) * sin(t);

    xmax = (xmax < x) ? x : xmax;
    xmin = (xmin > x) ? x : xmin;

    ymax = (ymax < y) ? y : ymax;
    ymin = (ymin > y) ? y : ymin;

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    max.push_back(xmax);
    max.push_back(ymax);

    return max;
}

void Display4() {

    double xmax, ymax, xmin, ymin;
    double a = 0.2, b = 0.2;
    double pi = 4 * atan(1.0);
    double ratia = 0.05;
    double t;
    std::vector<Punct> puncte;

    for (t = -pi + ratia; t < pi; t += ratia) {
        double x, y;
        std::vector<double> max = get_min_max(0.3, 0.2, 0.05, t);
        x = 2 * (a * cos(t) + b) * cos(t);
        x = x / (max[0] * 2);
        y = 2 * (a * cos(t) + b) * sin(t);
        y = y / max[1];
        Punct punct = Punct(x, y);
        puncte.push_back(punct);

    }
    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_LOOP);
    for (Punct punct : puncte) {
        glVertex2f(punct.x, punct.y);
    }
    glEnd();

}

void Display5() {

    double xmax, ymax, xmin, ymin;
    double a = 0.2;
    double pi = 4 * atan(1.0);
    double ratia = 0.005;
    double t;
    std::vector<Punct> puncte;

    for (t = -pi /2; t < pi/2; t += ratia) {
        if (abs(t) != pi / 6)
        {
            double x, y;
            x = a / (4 * pow(cos(t), 2) - 3);
            y = (a * tan(t)) / (4 * pow(cos(t), 2) - 3);
            Punct punct = Punct(x, y);
            puncte.push_back(punct);
        }
    }
    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_LOOP);
    for (Punct punct : puncte) {
        //glVertex2f(punct.x, punct.y);
        glVertex2f(0, punct.y);
        glVertex2f(punct.x,0);
    }
    glEnd();

}

void Init(void) {

    glClearColor(1.0, 1.0, 1.0, 1.0);

    glLineWidth(1);

    glPointSize(4);

    glPolygonMode(GL_FRONT, GL_LINE);
}

void Display(void) {
    glClear(GL_COLOR_BUFFER_BIT);

    switch (prevKey) {
    case '1':
        Display1();
        break;
    case '2':
        Display2();
        break;
    case '3':
        Display3();
        break;
    case '4':
        Display4();
        break;
    case '5':
        Display5();
        break;
        break;
    default:
        break;
    }

    glFlush();
}

void Reshape(int w, int h) {
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
}

void KeyboardFunc(unsigned char key, int x, int y) {
    prevKey = key;
    if (key == 27) // escape
        exit(0);
    glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) {
}

int main(int argc, char** argv) {

    glutInit(&argc, argv);

    glutInitWindowSize(dim, dim);

    glutInitWindowPosition(100, 100);

    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);

    glutCreateWindow(argv[0]);

    Init();

    glutReshapeFunc(Reshape);

    glutKeyboardFunc(KeyboardFunc);

    glutMouseFunc(MouseFunc);

    glutDisplayFunc(Display);

    glutMainLoop();

    return 0;
}
