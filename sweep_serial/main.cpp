#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "timer.h"

using namespace std;

double function(double x)
{
    return cos(x) / x;
}

double derFunction(double x)
{
    return (x * sin(x) + cos(x)) / (x * x);
}

void serialSweep(int n, double *a, double *c, double *b, double *f, double *p)
{
    double m;
    for (int i = 1; i < n; i++) {
        m = a[i] / c[i - 1];
        c[i] = c[i] - m * b[i - 1];
        f[i] = f[i] - m * f[i - 1];
    }

    p[n - 1] = f[n - 1] / c[n - 1];

    for (int i = n - 2; i >= 0; i--)
        p[i]=( f[i] - b[i] * p[i + 1] ) / c[i];
}

int main(int argc, char *argv[])
{
    int gridSize = 100;

    if (argc < 4) {
        cout << "Input Error: too few arguments." << endl;
        cout << "Usage: programm inputfile outputfile timefile" << endl;
        return -1;
    }
    ifstream fileIn(argv[1]);
    string input;
    getline(fileIn, input);
    istringstream iss(input);
    iss >> gridSize;
    if (gridSize < 2) {
        cout << "Input Error: grid size is wrong." << endl;
        return -1;
    }

    string fileOutName = argv[2];
    string fileTimeName = argv[3];

    double leftPoint = 1.;
    double rightPoint = 10.;

    double *subDiagonal = new double[gridSize - 1];
    double *supDiagonal = new double[gridSize - 1];
    double *mainDiagonal = new double[gridSize];
    double *f = new double[gridSize];

    double *p = new double[gridSize];

    double h = (rightPoint - leftPoint) / (gridSize - 1);
    for (int i = 0; i < gridSize - 1; i ++) {
        if (i == 0)
        {
            mainDiagonal[i] = 3./8 * h;//(h / 2);
            subDiagonal[i] = 1./8 * h;//(h / 2);
        }
        else if (i == 1)
        {
            subDiagonal[i] = 1./8 * h;
            supDiagonal[i - 1] = 1./8 * h;//(h / 2);
            mainDiagonal[i] = 3./8 * 2 * h;//(h + (h / 2));
        }
        else
        {
            subDiagonal[i] = 1./8 * h;
            supDiagonal[i - 1] = 1./8 * h;
            mainDiagonal[i] = 3./8 * (h + h);
        }
        if (i == 0)
        {
            f[i] = (function(leftPoint + (h / 2)) - function(leftPoint));
        }
        else
        {
            f[i] = (function(leftPoint + (i + 1) * h - (h / 2)) - function(leftPoint + i * h - (h / 2)));
        }
    }
    supDiagonal[gridSize - 2] = 1./8 * (h / 2);
    mainDiagonal[gridSize - 1] = 3./8 * (h / 2);
    f[gridSize - 1] = (function(rightPoint) - function(rightPoint - (h / 2)));


    Timer timer;
    timer.start();

    serialSweep(gridSize, subDiagonal, mainDiagonal, supDiagonal, f, p);

    timer.stop();

    double *realFunc = new double[4 * gridSize + 1];
    double *splineFunc = new double[4 * gridSize + 1];

    double *realFuncDer = new double[4 * gridSize + 1];
    double *splineFuncDer = new double[4 * gridSize + 1];

    for (int i = 0; i < 4 * (gridSize - 1) + 1; ++i) {
        double point = leftPoint + (h * i) / 4;
        realFunc[i] = function(point);
        realFuncDer[i] = derFunction(point);
    }

    for (int i = 0; i < (gridSize - 1); ++i) {
        double x = leftPoint + i * h;
        double piplus1 = leftPoint + (i + 1) * h - (h / 2);
        double currH = h;
        splineFuncDer[i * 4] = ((p[i + 1] + p[i]) / 2) + ((p[i + 1] + p[i]) * (x - piplus1) / currH);
        splineFunc[i * 4] = function(piplus1) + ((p[i + 1] + p[i]) * (x - piplus1) / 2) + ((p[i + 1] - p[i]) * (x - piplus1) *  (x - piplus1) / (2 * currH));
        x += h / 4;
        splineFunc[i * 4 + 1] = function(piplus1) + ((p[i + 1] + p[i]) * (x - piplus1) / 2) + ((p[i + 1] - p[i]) * (x - piplus1) *  (x - piplus1) / (2 * currH));
        splineFuncDer[i * 4 + 1] = ((p[i + 1] + p[i]) / 2) + ((p[i + 1] + p[i]) * (x - piplus1) / currH);
        x += h / 4;
        splineFunc[i * 4 + 2] = function(piplus1) + ((p[i + 1] + p[i]) * (x - piplus1) / 2) + ((p[i + 1] - p[i]) * (x - piplus1) *  (x - piplus1) / (2 * currH));
        splineFuncDer[i * 4 + 2] = ((p[i + 1] + p[i]) / 2) + ((p[i + 1] + p[i]) * (x - piplus1) / currH);
        x += h / 4;
        splineFunc[i * 4 + 3] = function(piplus1) + ((p[i + 1] + p[i]) * (x - piplus1) / 2) + ((p[i + 1] - p[i]) * (x - piplus1) *  (x - piplus1) / (2 * currH));
        splineFuncDer[i * 4 + 3] = ((p[i + 1] + p[i]) / 2) + ((p[i + 1] + p[i]) * (x - piplus1) / currH);
    }
    double x = rightPoint;
    splineFunc[(gridSize - 1) * 4] = function(x);
    splineFuncDer[(gridSize - 1) * 4] = p[gridSize];

    double maxError = 0.0;
    for (int i = 0; i < 4 * (gridSize - 1) + 1; ++i) {
        cout << realFunc[i] << " " << splineFunc[i] << endl;
        double currError = abs(realFunc[i] - splineFunc[i]);
        if (currError > maxError) {
            maxError = currError;
        }
    }
    //cout << maxError <<endl;

    double maxErrorD = 0.0;
    for (int i = 0; i < 4 * (gridSize - 1) + 1; ++i) {
        //cout << realFuncDer[i] << " " << splineFuncDer[i] << endl;
        double currErrorDer = abs(realFuncDer[i] - splineFuncDer[i]);
        if (currErrorDer > maxErrorD) {
            maxErrorD = currErrorDer;
        }
    }
    //cout << maxErrorD<<endl;

    ofstream fileOut(fileOutName.c_str());
    fileOut << maxError << endl << maxErrorD;
    ofstream fileTime(fileTimeName.c_str());
    fileTime << timer.getElapsed();

    delete []subDiagonal;
    delete []supDiagonal;
    delete []mainDiagonal;
    delete []f;
    delete []p;
    delete []realFunc;
    delete []realFuncDer;
    delete []splineFunc;
    delete []splineFuncDer;

    return 0;
}

