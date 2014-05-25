#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "timer.h"

using namespace std;

double function(double x)
{
    //return sin(x) / x;
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
    if (argc < 4) {
        cout << "Input Error: too few arguments." << endl;
        cout << "Usage: programm inputfile outputfile timefile" << endl;
        return -1;
    }
    ifstream fileIn(argv[1]);
    string input;
    getline(fileIn, input);
    istringstream iss(input);

    int gridSize;

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

    double *a = new double[gridSize - 1];
    double *b = new double[gridSize - 1];
    double *c = new double[gridSize - 1];

    double *p = new double[gridSize];

    double h = (rightPoint - leftPoint) / (gridSize - 1);
    for (int i = 0; i < gridSize - 1; i ++) {
        subDiagonal[i] = 1./8;
        supDiagonal[i] = 1./8;
        mainDiagonal[i] = 3./4;
        f[i] = (function(leftPoint + (i + 1) * h) - function(leftPoint + i * h)) / h;
    }
    mainDiagonal[gridSize - 1] = 3./4;
    f[gridSize - 1] = (function(rightPoint + h) - function(rightPoint)) / h;

    Timer timer;
    timer.start();

    serialSweep(gridSize, subDiagonal, mainDiagonal, supDiagonal, f, p);

    timer.stop();

    for (int j = 1; j < gridSize; ++j) {
        double xj = leftPoint + (j - 1) * h;
        a[j - 1] = (p[j] - p[j - 1]) / (2 * h);
        b[j - 1] = (xj / h) * (p[j - 1] - p[j]) + (p[j] + p[j - 1]) / 2;
        c[j - 1] = function(xj) + ((p[j] - p[j - 1]) * xj * xj) / (2 * h) - (p[j] + p[j - 1]) / 2 * xj;
    }

    double *realFunc = new double[4 * gridSize + 1];
    double *splineFunc = new double[4 * gridSize + 1];

    double *realFuncDer = new double[4 * gridSize + 1];
    double *splineFuncDer = new double[4 * gridSize + 1];

    double hDer = (double)h / 3.0;
    for (int i = 0; i < 4 * (gridSize - 1) + 1; ++i) {
        double point = leftPoint + (h * i) / 4;
        realFunc[i] = function(point);
        realFuncDer[i] = derFunction(point);
    }

    for (int i = 0; i < (gridSize - 1); ++i) {
        double x = leftPoint + i * h;
        splineFunc[i * 4] = a[i] * x * x + b[i] * x + c[i];
        splineFuncDer[i * 4] = a[i] * x + b[i];
        x += h / 4;
        splineFunc[i * 4 + 1] = a[i] * x * x + b[i] * x + c[i];
        splineFuncDer[i * 4 + 1] = a[i] * x + b[i];
        x += h / 4;
        splineFunc[i * 4 + 2] = a[i] * x * x + b[i] * x + c[i];
        splineFuncDer[i * 4 + 2] = a[i] * x + b[i];
        x += h / 4;
        splineFunc[i * 4 + 3] = a[i] * x * x + b[i] * x + c[i];
        splineFuncDer[i * 4 + 3] = a[i] * x + b[i];
    }
    double x = rightPoint;
    splineFunc[(gridSize - 1) * 4] = a[gridSize - 2] * x * x + b[gridSize - 2] * x + c[gridSize - 2];
    splineFuncDer[(gridSize - 1) * 4] = (a[gridSize - 2] * (x + hDer) * (x + hDer) + b[gridSize - 2] * (x + hDer) - (a[gridSize - 2] * (x - hDer) * (x - hDer) + b[gridSize - 2] * (x - hDer) )) / hDer;

    double maxError = 0.0;
    for (int i = 0; i < 4 * (gridSize - 1) + 1; ++i) {
        //cout << realFunc[i] << " " << splineFunc[i] << endl;
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

    //cout<<timer.getElapsed()<<endl;

    ofstream fileOut(fileOutName.c_str());
    fileOut << maxError << endl << maxErrorD;
    ofstream fileTime(fileTimeName.c_str());
    fileTime << timer.getElapsed();

    return 0;
}

