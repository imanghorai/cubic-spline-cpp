#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void sort(const int n, vector<double>& x, vector<double>& y) {
    for (int i = 0; i < n; i++) {
        int minIndex = i;
        for (int j = i + 1; j < n; j++) {
            if (x[j] < x[minIndex]) {
                minIndex = j;
            }
        }
        swap(x[minIndex], x[i]);
        swap(y[minIndex], y[i]);
    }
}

MatrixXd CubicSpline(int n, const vector<double>& x, const vector<double>& y) {
    MatrixXd A(4*(n-1), 4*(n-1));
    A.setZero();
    VectorXd cons(4*(n-1));
    MatrixXd coeffs(n-1, 4); 
    int row = 0;

    for (int i = 0; i < n - 1; i++) {
        A(row, 4 * i) = 1;  
        cons(row) = y[i];
        row++;

        A(row, 4 * i) = 1;  
        A(row, 4 * i + 1) = x[i + 1] - x[i];  
        A(row, 4 * i + 2) = pow(x[i + 1] - x[i], 2);  
        A(row, 4 * i + 3) = pow(x[i + 1] - x[i], 3); 
        cons(row) = y[i + 1];
        row++;
    }

    for (int i = 0; i < n - 2; i++) {
        A(row, 4 * i + 1) = 1;  // b[i]
        A(row, 4 * i + 2) = 2 * (x[i + 1] - x[i]);  // c[i]
        A(row, 4 * i + 3) = 3 * pow(x[i + 1] - x[i], 2);  // d[i]

        A(row, 4 * (i + 1) + 1) = -1;  // b[i+1]
        cons(row) = 0;
        row++;
    }

    for (int i = 0; i < n - 2; i++) {
        A(row, 4 * i + 2) = 2;  // c[i]
        A(row, 4 * i + 3) = 6 * (x[i + 1] - x[i]);  // d[i]

        A(row, 4 * (i + 1) + 2) = -2;  // c[i+1]
        cons(row) = 0;
        row++;
    }

    A(row, 2) = 2;
    cons(row) = 0;
    row++;

    A(row, 4 * (n - 2) + 2) = 2;
    A(row, 4 * (n - 2) + 3) = 6 * (x[n - 1] - x[n - 2]);
    cons(row) = 0;

    VectorXd solution = A.colPivHouseholderQr().solve(cons);

    for (int i = 0; i < n - 1; i++) {
        coeffs(i, 0) = solution[4 * i];
        coeffs(i, 1) = solution[4 * i + 1];  
        coeffs(i, 2) = solution[4 * i + 2];
        coeffs(i, 3) = solution[4 * i + 3];
    }

    return coeffs;
}

double interpolate(const double xval, const int n, const vector<double>& x, const MatrixXd& coeffs) {
    double yval = 0;
    if (xval < x[0] || xval > x[n-1]) {
        cout << "X value out of range" << endl;
        return 0;
    } else {
        for (int i = 0; i < n; i++) {
            if (xval == x[i]) {
                return coeffs(i, 0);
            }
        }
        for (int i = 0; i < n - 1; i++) {
            if (xval > x[i] && xval < x[i + 1]) {
                double dx = xval - x[i];
                yval = coeffs(i, 0) + coeffs(i, 1) * dx + coeffs(i, 2) * pow(dx, 2) + coeffs(i, 3) * pow(dx, 3);
                return yval;
            }
        }
    }
    return 0;
}

int main() {
    cout << "Enter the number of points: ";
    int n;
    cin >> n;

    vector<double> x(n), y(n);
    cout << "Enter the X coordinates: ";
    for (int i = 0; i < n; i++) cin >> x[i];

    cout << "Enter the Y coordinates: ";
    for (int i = 0; i < n; i++) cin >> y[i];

    sort(n, x, y);

    MatrixXd coeffs = CubicSpline(n, x, y);

    cout << "i-->";
    for (int i = 0; i < n - 1; i++) cout << "\t\t[" << i << "]";
    cout << "\na[i]\t\t";
    for (int i = 0; i < n - 1; i++) cout << coeffs(i, 0) << "\t\t";
    cout << "\nb[i]\t\t";
    for (int i = 0; i < n - 1; i++) cout << coeffs(i, 1) << "\t\t";
    cout << "\nc[i]\t\t";
    for (int i = 0; i < n - 1; i++) cout << coeffs(i, 2) << "\t\t";
    cout << "\nd[i]\t\t";
    for (int i = 0; i < n - 1; i++) cout << coeffs(i, 3) << "\t\t";
    cout << "\n";

    cout << "Enter the value of x which you want to interpolate: ";
    double xval;
    cin >> xval;
    double yval = interpolate(xval, n, x, coeffs);
    cout << "The interpolated value is: " << yval << endl;

    return 0;
}
