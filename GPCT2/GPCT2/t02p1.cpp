
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


/// EXERCITIUL 1 ///

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
    double ratia = 0.05;
    double f;
    std::vector<double> maximum;
    std::vector<Punct> puncte;

    for (x = 0 + ratia; x <= 75; x = x + ratia) {

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
    }
    glEnd();
}

/// EXERCITIUL 2 - melcul lui Pascal-concoida cercului :
//  x=2*(a*cost+b)*cost,y=2*(a*cost+b)*sint,t->(-pi,pi)

std::vector<double> get_min_max(double a, double b, double t) {
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
    double a = 0.3, b = 0.2;
    double pi = 4 * atan(1.0);
    double ratia = 0.05;
    double t;
    std::vector<Punct> puncte;

    for (t = -pi + ratia; t < pi; t += ratia) {
        double x, y;
        std::vector<double> max = get_min_max(0.3, 0.2, t);
        x = 2 * (a * cos(t) + b) * cos(t);
        x = x / (max[0] * 2);
        y = 2 * (a * cos(t) + b) * sin(t);
        y = y / (max[1] / 2);
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

/// EXERCITIUL 3 - TRISECTOAREA LUI LONGCHAMPS 


std::vector<double> get_max_trisectoare() {
    double pi = atan(1) * 4;
    double ratia = 0.025;
    double x, y;
    double xmax, ymax, xmin, ymin;
    std::vector<double> maximum;

    xmin = xmax = (0.2 / (4 * cos(-pi / 2 + ratia) * cos(-pi / 2 + ratia) - 3));
    ymin = ymax = (0.2 * tan(-pi / 2 + ratia)) / (4 * cos(-pi / 2 + ratia) * cos(-pi / 2 + ratia) - 3);
    for (double t = -pi / 2 + ratia; t < pi / 2; t = t + ratia) {
        x = (0.2 / (4 * cos(t) * cos(t) - 3));
        y = (0.2 * tan(t)) / (4 * cos(t) * cos(t) - 3);

        xmax = (xmax < x) ? x : xmax;
        xmin = (xmin > x) ? x : xmin;
        ymax = (ymax < y) ? y : ymax;
        ymin = (ymin > y) ? y : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);


    maximum.push_back(xmax);
    maximum.push_back(ymax);

    return maximum;
}
std::vector < std::pair<double, double>> vf_triunghi(double ratia) {

    std::vector < std::pair<double, double>> varf_triunghi;
    double xmax, ymax;
    double x, y;
    double t;
    double pi = 4 * atan(1);
    std::vector < double> maxx = get_max_trisectoare();

    xmax = maxx[0] / 3.5;
    ymax = maxx[1] / 3.5;

    for (t = -pi / 2 + ratia; t < -pi / 6; t += ratia) {
        if ((t != pi / 6) and (t != -pi / 6)) {
            x = (0.2 / (4 * cos(t) * cos(t) - 3) / xmax) * 2;
            y = (0.2 * tan(t) / (4 * cos(t) * cos(t) - 3)) / ymax;

            for (double r = t; r < t + ratia && r < -pi / 6 - ratia; r += ratia / 9) {
                varf_triunghi.emplace_back((0.2 / (4 * cos(r) * cos(r) - 3) / xmax) * 2,
                    (0.2 * tan(r) / (4 * cos(r) * cos(r) - 3)) / ymax);
            }
       
        }
    }
    return varf_triunghi;
}

void get_margine() {
    double t;
    double const pi = 4 * atan(1);
    double ratia = 0.025;

    double xmax, ymax;
    double x1, y1, x2, y2, x3, y3;

    int inaltime=((-pi / 6 - ratia + pi / 2) / ratia);
    double latime=(-pi / 2 + ratia * (inaltime + 1));

    std::vector < std::pair<double, double>> triunghi = vf_triunghi(ratia);
    std::vector < double> maxx = get_max_trisectoare();

    xmax = maxx[0] / 3.5;
    ymax = maxx[1] / 3.5;

    glColor3f(0.1, 0.1, 1);
    glBegin(GL_LINE_STRIP);

    x1 = (0.2 / (4 * cos(latime) * cos(latime) - 3) / xmax) * 2;
    y1 = (0.2 * tan(latime) / (4 * cos(latime) * cos(latime) - 3)) / ymax;

    x2 = (0.2 / (4 * cos((-pi / 2 + ratia) * cos(-pi / 2 + ratia) - 3) / xmax)) * 2;
    y2 = (0.2 * tan(-pi / 2 + ratia) / (4 * cos(-pi / 2 + ratia) * cos(-pi / 2 + ratia) - 3)) / ymax;

    glVertex2f(x1, y1);
    glVertex2f(x1, y2);
    glVertex2f(x2, y2);

    for (t = -pi / 2 + ratia; t < -pi / 6; t += ratia) {
        if ((t != pi / 6) and (t != -pi / 6)) {
            x3 = (0.2 / (4 * cos(t) * cos(t) - 3) / xmax) * 2;
            y3 = (0.2 * tan(t) / (4 * cos(t) * cos(t) - 3)) / ymax;
            glVertex2f(x3, y3);
        }
    }
    glEnd();
}

void  get_interior_image() {
    double xmax, ymax;
    double x1, y1, x2, y2, x_aux, y_aux,x_vf,y_vf;
    double ratia = 0.025;
    double pi = 4 * atan(1);
    int inaltime = ((-pi / 6 - ratia + pi / 2) / ratia);
    double latime= (-pi / 2 + ratia * (inaltime + 1));
    std::vector < std::pair<double, double>> triunghi = vf_triunghi(ratia);

    std::vector < double> maxx = get_max_trisectoare();
    xmax = maxx[0] / 3.5;
    ymax = maxx[1] / 3.5;

    glColor3f(1, 0.1, 0.1); // rosu
    glPolygonMode(GL_FRONT, GL_FILL);
    glBegin(GL_TRIANGLES);

    x_aux = (0.2 / (4 * cos(latime) * cos(latime) - 3) / xmax) * 2;
    y_aux = (0.2 * tan(-pi / 2 + ratia) / (4 * cos(-pi / 2 + ratia) * cos(-pi / 2 + ratia) - 3)) / ymax;
    x_vf = x_aux;
    y_vf = y_aux;

    for (int i = 0; i < triunghi.size() - 1; i += 3) {
        x1 = triunghi[i].first;
        y1 = triunghi[i].second;

        x2 = triunghi[i + 1].first;
        y2 = triunghi[i + 1].second;

        glVertex2f(x1, y1);
        glVertex2f(x_vf, y_vf);
        glVertex2f(x2, y2);
    }
    glEnd();
}
void Display5() {
  
    get_interior_image();
    get_margine();
  
}


std::vector<double> get_max_cicloida() {
    double pi = atan(1) * 4;
    double ratia = 0.05;
    double x, y;
    double xmax, ymax, xmin, ymin;
    std::vector<double> maximum;
    double a = 0.1;
    double b = 0.2;


    xmin = -1000;
    xmax = 1000;
    ymin = -0.1;
    ymax = 0.1;

    for (double t = -pi / 2 + ratia; t < pi / 2; t = t + ratia) {
        x = a * t - b * sin(t);
        y = a - b * cos(t);

        xmax = (xmax < x) ? x : xmax;
        xmin = (xmin > x) ? x : xmin;
        ymax = (ymax < y) ? y : ymax;
        ymin = (ymin > y) ? y : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);


    maximum.push_back(xmax);
    maximum.push_back(ymax);

    return maximum;
}

std::vector<Punct> get_puncts_cicloida() {
    double ratia = 0.05;
    double x, y, t;
    double a = 0.1;
    double b = 0.2;

    std::vector<Punct> puncte;
    std::vector<double> max = get_max_cicloida();

    double x_max = max[0];
    double y_max = max[1];
    for (t = -100+ratia; t <= 100; t = t + ratia)
    {
        x = 0.1 * t - 0.2 * sin(t);
        y = 0.1 - 0.2 * cos(t);
        Punct punct = Punct(x, y);
        puncte.push_back(punct);

    }
    return puncte;
}

void Display6() {
 
    double const pi = 4 * atan(1);
    double ratia = 0.05;
    double t;
    double const min_v = -100;
    double const max_v = +100;

    std::vector<Punct> vf = get_puncts_cicloida();
   

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (Punct p :vf) {
        glVertex2f(p.x, p.y);
    }
    glEnd();

}

std::vector<double> get_max_epicicloida() {
    double pi = atan(1) * 4;
    double r = 0.3;
    double R = 0.1;
    double ratia = 0.05;
    double x, y;
    double xmax, ymax, xmin, ymin;
    std::vector<double> maximum;
    double a = 0.1;
    double b = 0.2;


    xmin = -10;
    xmax = 10;
    ymin = -0.1;
    ymax = 0.1;

    for (double t = 0+ ratia; t < pi *2; t = t + ratia) {
        x = (R + r) * cos((r * t) / R) - r * cos(t + (r * t) / R);
        y = (R + r) * sin((r * t) / R) - r * sin(t + (r * t) / R);

        xmax = (xmax < x) ? x : xmax;
        xmin = (xmin > x) ? x : xmin;
        ymax = (ymax < y) ? y : ymax;
        ymin = (ymin > y) ? y : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);


    maximum.push_back(xmax);
    maximum.push_back(ymax);

    return maximum;
}

std::vector<Punct>get_puncte_epicicloida() {
    double ratia = 0.05;
    double pi = 4 * atan(1);
    double x, y, t;
    double r = 0.3;
    double R = 0.1;

    std::vector<Punct> puncte;
    std::vector<double> max = get_max_epicicloida();

    double x_max = max[0];
    double y_max = max[1];
 
    for (double t = 0 + ratia; t <2* pi ; t = t + ratia) {
        x = (R + r) * cos((r * t) / R) - r * cos(t + (r * t) / R);
        y = (R + r) * sin((r * t) / R) - r * sin(t + (r * t) / R);
        Punct punct = Punct(x, y);
        puncte.push_back(punct);

    }
    return puncte;
}
void Display7() {
    std::vector<Punct> vf = get_puncte_epicicloida();


    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (Punct p : vf) {
        glVertex2f(p.x, p.y);
    }
    glEnd();

}

std::vector<double> get_max_hipocicloida() {
    double pi = atan(1) * 4;
    double r = 0.3;
    double R = 0.1;
    double ratia = 0.05;
    double x, y;
    double xmax, ymax, xmin, ymin;
    std::vector<double> maximum;
   
    xmin = -10;
    xmax = 10;
    ymin = -0.1;
    ymax = 0.1;

    for (double t = 0 + ratia; t < pi * 2; t = t + ratia) {
        x = (R - r) * cos((r * t) / R) - r * cos(t - (r * t) / R);
        y = (R - r) * sin((r * t) / R) - r * sin(t - (r * t) / R);

        xmax = (xmax < x) ? x : xmax;
        xmin = (xmin > x) ? x : xmin;
        ymax = (ymax < y) ? y : ymax;
        ymin = (ymin > y) ? y : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);


    maximum.push_back(xmax);
    maximum.push_back(ymax);

    return maximum;
}


std::vector<Punct>get_puncte_hipocicloida() {
    double ratia = 0.05;
    double pi = 4 * atan(1);
    double x, y, t;
    double r = 0.3;
    double R = 0.1;

    std::vector<Punct> puncte;
    std::vector<double> max = get_max_hipocicloida();

    double x_max = max[0];
    double y_max = max[1];

    for (double t = 0 + ratia; t < 2 * pi; t = t + ratia) {
        x = (R - r) * cos((r * t) / R) - r * cos(t - (r * t) / R);
        y = (R - r) * sin((r * t) / R) - r * sin(t - (r * t) / R);
        Punct punct = Punct(x, y);
        puncte.push_back(punct);

    }
    return puncte;
}
void Display8() {
    std::vector<Punct> vf = get_puncte_hipocicloida();


    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (Punct p : vf) {
        glVertex2f(p.x, p.y);
    }
    glEnd();

}

std::vector<double> get_max_Bernoulli() {
    double pi = atan(1) * 4;
    double a = 0.4;
    double ratia = pi/500;
    double x, y;
    double r1,r2;
    double xmax, ymax, xmin, ymin;
    std::vector<double> maximum;


    xmin = -0.1;
    xmax = 0.75;
    ymin = -0.1;
    ymax = 0.75;

    for (double t = -pi/4 + ratia; t < pi/4 ; t = t + ratia) {
        r1 = a * sqrt(2 * cos(2*t));
        r2 = -a * sqrt(2 * cos(2 * t));
        x = r1 * cos(t); 
        y = r1 * sin(t);

        xmax = (xmax < x) ? x : xmax;
        xmin = (xmin > x) ? x : xmin;
        ymax = (ymax < y) ? y : ymax;
        ymin = (ymin > y) ? y : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);


    maximum.push_back(xmax);
    maximum.push_back(ymax);

    return maximum;
}

std::vector<Punct>get_puncte_Bernoulli1() {
    
    double pi = 4 * atan(1);
    double x1, y1,x2,y2 ,t;
    double a=0.4;
    double r1, r2;
    double ratia = pi / 500;

    std::vector<Punct> puncte1;
    std::vector<Punct> puncte2;
    std::vector<double> max = get_max_Bernoulli();

    double x_max = 0.75;
    double y_max = 0.75;

    for (double t = -pi/4 + ratia; t < pi/4; t = t + ratia) {
        r1 = a * sqrt(2 * cos(2 * t));
        r2 = -a * sqrt(2 * cos(2 * t));
        x1 = r1 * cos(t)/x_max;
        y1 = r1 * sin(t)/y_max;
       
        Punct punct1 = Punct(x1, y1);
        puncte1.push_back(punct1);

    }
    return puncte1;
}
std::vector<Punct> get_puncte_Bernoulli2() {
    
    double pi = 4 * atan(1);
    double x1, y1, x2, y2, t;
    double a = 0.4;
    double r1, r2;
    double ratia = pi/500;
    std::vector<Punct> puncte2;
    std::vector<double> max = get_max_Bernoulli();

    double x_max = 0.75;
    double y_max = 0.75;

    for (double t = pi / 4-ratia ; t >-pi / 4; t = t - ratia) {
        r1 = a * sqrt(2 * cos(2 * t));
        r2 = -a * sqrt(2 * cos(2 * t));
        x1 = r1 * cos(t)/x_max;
        y1 = r1 * sin(t)/y_max;
        x2 = -x1;
        y2 = -y1;

        Punct punct2 = Punct(x2, y2);
        puncte2.push_back(punct2);

    }
    return  puncte2;
}
void Display9() {
    std::vector<Punct> vf1 = get_puncte_Bernoulli1();
    std::vector<Punct> vf2 = get_puncte_Bernoulli2();

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (Punct p : vf2) {
        glVertex2f(p.x, p.y);
    }
    for (Punct p : vf1) {
        glVertex2f(p.x, p.y);
    }
    glEnd();

}

std::vector<double> get_max_spirala() {
    double pi = atan(1) * 4;
    double a = 0.02;
    double ratia = 0.05;
    double x, y;
    double r ;
    double xmax, ymax, xmin, ymin;
    std::vector<double> maximum;


    xmin = -25;
    xmax = 50;
    ymin = -25;
    ymax = 50;

    for (double t = 0 + ratia; t < 100; t = t + ratia) {
        r = a * exp(1 + t);
        x = r * cos(t);
        y = r * sin(t);

        xmax = (xmax < x) ? x : xmax;
        xmin = (xmin > x) ? x : xmin;
        ymax = (ymax < y) ? y : ymax;
        ymin = (ymin > y) ? y : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);


    maximum.push_back(xmax);
    maximum.push_back(ymax);

    return maximum;
}

std::vector<Punct>get_puncte_spirala() {

    double pi = 4 * atan(1);
    double x, y, t;
    double a = 0.4;
    double r;
    double ratia = 0.05;

    std::vector<Punct> puncte;
    std::vector<double> max = get_max_spirala();

    double x_max = 50;
    double y_max = 50;

    for (double t = 0 + ratia; t < 100; t = t + ratia) {
        r = a * exp(t+1);

        x = r * sin(t) / y_max;
        y = - r * cos(t) / x_max;
       
        Punct punct = Punct(x, y);
        puncte.push_back(punct);

    }
    return puncte;
}

void Display0() {

    std::vector<Punct> vf = get_puncte_spirala();
    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (Punct p : vf) {
        glVertex2f(p.x, p.y);
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
    case '6':
        Display6();
        break;
    case '7':
        Display7();
        break;
    case '8':
        Display8();
        break;
    case '9':
        Display9();
        break;
    case '0':
        Display0();
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
