#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "timer.h"
#include <cilk/cilk.h>

#define DEBUG 0
#define THREAD_COUNT 4

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
        m = a[i - 1] / c[i - 1];
        c[i] = c[i] - m * b[i - 1];
        f[i] = f[i] - m * f[i - 1];
    }

    p[n - 1] = f[n - 1] / c[n - 1];

    for (int i = n - 2; i >= 0; i--){
        p[i] = ( f[i] - b[i] * p[i + 1] ) / c[i];
    }
}

double getSplinePointY(int knotCount, double* knots,
                       double* nodesX, double* nodesY,
                       double* coeffs,
                       double x)
{
    int i;
    double c;
    for (i = 1; i < knotCount; ++i){
        if (x <= knots[i]){
            break;
        }
    }

    if (i < knotCount){
        c = x - nodesX[i];
        return nodesY[i] +
                0.5 * (coeffs[i] + coeffs[i - 1]) * c +
                (coeffs[i] - coeffs[i - 1]) * c * c / (2 * (knots[i] - knots[i - 1]));
    }

    return 0.0;
}

double getSplineDerPointY(int knotCount, double* knots,
                       double* nodesX, double* nodesY,
                       double* coeffs,
                       double x)
{
    int i;
    double c;
    for (i = 1; i < knotCount; ++i){
        if (x <= knots[i]){
            break;
        }
    }

    if (i < knotCount){
        c = x - nodesX[i];
        return 0.5 * (coeffs[i] + coeffs[i - 1]) +
                (coeffs[i] - coeffs[i - 1]) * c / (2 * (knots[i] - knots[i - 1]));
    }

    return 0.0;
}

void doForward(double *subDiagonal, double *mainDiagonal, double *supDiagonal, double *b, int blockRowCount, int blockIndex)
{
    int blockStartRow = blockIndex * blockRowCount;
    int blockEndRow   = blockStartRow + blockRowCount - 1;

    for(int blockRowIndex = blockStartRow; blockRowIndex < blockEndRow; blockRowIndex++){
        double factor = subDiagonal[blockRowIndex] / mainDiagonal[blockRowIndex];

        if(blockIndex > 0){
            subDiagonal[blockRowIndex] = -subDiagonal[blockRowIndex - 1] * factor;// Non-diagonal elements
        }else{
            subDiagonal[blockRowIndex] = 0.0;
        }

        mainDiagonal[blockRowIndex + 1] -= supDiagonal[blockRowIndex] * factor;
        b[blockRowIndex + 1] -= b[blockRowIndex] * factor;
    }
}

void doBackward(double *subDiagonal, double *mainDiagonal, double *supDiagonal, double *b, int blockRowCount, int blockIndex)
{
    int blockStartRow = blockIndex * blockRowCount;
    int blockEndRow   = blockStartRow + blockRowCount - 1;
        blockStartRow = max(blockStartRow, 1);

    for (int blockRowIndex = blockEndRow - 1; blockRowIndex >= blockStartRow; blockRowIndex--){
        double factor = supDiagonal[blockRowIndex - 1] / mainDiagonal[blockRowIndex];

        if(blockIndex > 0){
            if(blockRowIndex > blockStartRow){
                subDiagonal[blockRowIndex - 2] -= subDiagonal[blockRowIndex - 1] * factor;
            }else{
                mainDiagonal[blockRowIndex - 1] -= subDiagonal[blockRowIndex - 1] * factor;
            }
        }
        supDiagonal[blockRowIndex - 1] = -supDiagonal[blockRowIndex] * factor;
        b[blockRowIndex - 1] -= b[blockRowIndex] * factor;
    }
}

void doAuxilaryMatrix(double *subDiagonal, double *mainDiagonal, double *supDiagonal, double *b, int blockRowCount, int blockIndex, double **auxiliaryMatrix, int blockCount)
{
    int blockStartRow = blockIndex * blockRowCount;
    int blockEndRow   = blockStartRow + blockRowCount - 1;

    if(blockIndex > 0)
    {
        auxiliaryMatrix[0][blockIndex - 1] = subDiagonal[blockEndRow - 1];
    }

    if(blockIndex < blockCount - 1)
    {
        auxiliaryMatrix[2][blockIndex] = supDiagonal[blockEndRow];
    }

    auxiliaryMatrix[1][blockIndex] = mainDiagonal[blockEndRow];
    auxiliaryMatrix[3][blockIndex] = b[blockEndRow];
}

void doX(double *auxiliaryX, double *x, int blockRowCount, int blockIndex)
{
    int blockStartRow = blockIndex * blockRowCount;
    int blockEndRow   = blockStartRow + blockRowCount - 1;

    x[blockEndRow] = auxiliaryX[blockIndex];
}

void doCalculate(double *subDiagonal, double *mainDiagonal, double *supDiagonal, double *b, double *x, int blockRowCount, int blockIndex)
{
    int blockStartRow = blockIndex * blockRowCount;
    int blockEndRow   = blockStartRow + blockRowCount - 1;

    for(int blockRowIndex = blockStartRow; blockRowIndex < blockEndRow; blockRowIndex++){
        x[blockRowIndex] = b[blockRowIndex];

        if(blockIndex > 0){
            x[blockRowIndex] -= subDiagonal[blockRowIndex - 1] * x[blockStartRow - 1];
        }

        x[blockRowIndex] -= supDiagonal[blockRowIndex] * x[blockEndRow];
        x[blockRowIndex] /= mainDiagonal[blockRowIndex];
    }
}

void cilkSweep(int equationCount,double* subDiagonal, double* mainDiagonal, double* supDiagonal,
              double* b, int blockCount, double* x)
{
    int blockRowCount = equationCount / blockCount;
    double* auxiliaryMatrix[5];
    auxiliaryMatrix[0] = new double [blockCount - 1];
    auxiliaryMatrix[1] = new double [blockCount];
    auxiliaryMatrix[2] = new double [blockCount - 1];
    auxiliaryMatrix[3] = new double [blockCount];
    auxiliaryMatrix[4] = new double [blockCount];

//    omp_set_num_threads(blockCount);

//#pragma omp parallel for shared(subDiagonal, mainDiagonal, supDiagonal, b, blockRowCount)
    cilk_for (int i = 0; i < blockCount; ++i)
    {
        doForward(subDiagonal, mainDiagonal, supDiagonal, b, blockRowCount, i);
    }
//#pragma omp parallel for shared(subDiagonal, mainDiagonal, supDiagonal, b, blockRowCount)
    cilk_for (int i = 0; i < blockCount; ++i)
    {
        doBackward(subDiagonal, mainDiagonal, supDiagonal, b, blockRowCount, i);
    }
//#pragma omp parallel for shared(subDiagonal, mainDiagonal, supDiagonal, b, blockRowCount, auxiliaryMatrix)
    cilk_for (int i = 0; i < blockCount; ++i)
    {
        doAuxilaryMatrix(subDiagonal, mainDiagonal, supDiagonal, b, blockRowCount, i, auxiliaryMatrix, blockCount);
    }
    serialSweep(blockCount,auxiliaryMatrix[0], auxiliaryMatrix[1], auxiliaryMatrix[2], auxiliaryMatrix[3], auxiliaryMatrix[4]);
//#pragma omp parallel for shared(x, blockRowCount, auxiliaryMatrix)
    cilk_for (int i = 0; i < blockCount; ++i)
    {
        doX(auxiliaryMatrix[4], x, blockRowCount, i);
    }
//#pragma omp parallel for shared(subDiagonal, mainDiagonal, supDiagonal, b, blockRowCount, x)
    cilk_for (int i = 0; i < blockCount; ++i)
    {
        doCalculate(subDiagonal, mainDiagonal, supDiagonal, b, x, blockRowCount, i);
    }

    delete []auxiliaryMatrix[0];
    delete []auxiliaryMatrix[1];
    delete []auxiliaryMatrix[2];
    delete []auxiliaryMatrix[3];
    delete []auxiliaryMatrix[4];
}

int main(int argc, char *argv[])
{
#if DEBUG
    int gridSize = 1000;
#else
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
#endif

    double leftPoint = 1.;
    double rightPoint = 10.;

    double *subDiagonal = new double[gridSize - 1];
    double *supDiagonal = new double[gridSize - 1];
    double *mainDiagonal = new double[gridSize];
    double *f = new double[gridSize];

    double *p = new double[gridSize];

    int i = 0;
    double h = (rightPoint - leftPoint) / (gridSize - 1);

    double *knots = new double[gridSize];
    knots[0] = leftPoint;
    knots[gridSize - 1] = rightPoint;

    double *nodesX = new double[gridSize + 1];
    double *nodesY = new double[gridSize + 1];
    nodesX[0] = knots[0];
    nodesY[0] = function(nodesX[0]);

    for (i = 1; i < gridSize - 1; i++)
    {
        knots[i] = i * h + leftPoint;

        nodesX[i] = (knots[i] + knots[i - 1]) / 2;
        nodesY[i] = function(nodesX[i]);
    }

    i = gridSize - 1;
    nodesX[i] = (knots[i] + knots[i - 1]) / 2;
    nodesY[i] = function(nodesX[i]);

    nodesX[gridSize] = knots[gridSize - 1];
    nodesY[gridSize] = function(nodesX[gridSize]);


    for (i = 1; i < gridSize - 1; i++)
    {
        subDiagonal[i - 1] = knots[i] - knots[i - 1];
        supDiagonal[i - 1] = knots[i] - knots[i - 1];
        mainDiagonal[i]    = 3 * (knots[i + 1] - knots[i - 1]);
        f[i - 1]           = 8 * (nodesY[i] - nodesY[i - 1]);
    }

    mainDiagonal[0] = 3 * (knots[1] - knots[0]);
    mainDiagonal[gridSize - 1] = 3 * (knots[gridSize - 1] - knots[gridSize - 2]);
    subDiagonal[gridSize - 2]  = knots[gridSize - 1] - knots[gridSize - 2];
    supDiagonal[gridSize - 2]  = knots[gridSize - 1] - knots[gridSize - 2];
    f[gridSize - 2]            = 8 * (nodesY[gridSize - 1] - nodesY[gridSize - 2]);
    f[gridSize - 1]            = 8 * (nodesY[gridSize] - nodesY[gridSize - 1]);

    Timer timer;
    timer.start();

//    serialSweep(gridSize, subDiagonal, mainDiagonal, supDiagonal, f, p);
    cilkSweep(gridSize, subDiagonal, mainDiagonal, supDiagonal, f, THREAD_COUNT, p);

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

    for (int i = 0; i < 4 * (gridSize - 1) + 1; ++i) {
        double point = leftPoint + (h * i) / 4;
        splineFunc[i] = getSplinePointY(gridSize, knots, nodesX, nodesY, p, point);
        splineFuncDer[i] = getSplineDerPointY(gridSize, knots, nodesX, nodesY, p, point);
    }

    double maxError = 0.0;
    for (int i = 0; i < 4 * (gridSize - 1) + 1; ++i) {
        double currError = abs(realFunc[i] - splineFunc[i]);
        if (currError > maxError) {
            maxError = currError;
        }
    }
    cout << maxError <<endl;

    double maxErrorD = 0.0;
    for (int i = 0; i < 4 * (gridSize - 1) + 1; ++i) {
        //cout << realFuncDer[i] << " " << splineFuncDer[i] << endl;
        double currErrorDer = abs(realFuncDer[i] - splineFuncDer[i]);
        if (currErrorDer > maxErrorD) {
            maxErrorD = currErrorDer;
        }
    }
    //cout << maxErrorD<<endl;
#if DEBUG

#else
    ofstream fileOut(fileOutName.c_str());
    fileOut << maxError << endl << maxErrorD;
    ofstream fileTime(fileTimeName.c_str());
    fileTime << timer.getElapsed();
#endif

    delete []subDiagonal;
    delete []supDiagonal;
    delete []mainDiagonal;
    delete []f;
    delete []p;
    delete []realFunc;
    delete []realFuncDer;
    delete []splineFunc;
    delete []splineFuncDer;
    delete []nodesX;
    delete []nodesY;
    delete []knots;

    return 0;
}

