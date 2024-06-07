#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
double func(double x, double a, double b, double c, double d) { //integrand
    return a*exp(-d*pow(x-c, 2))+5.14;
}
//left rule (identical to the Newton-Cotes rule with zero interpolation)
double LeftRectangleMethod(double xn, double xk, double a, double b, double c, double d, int n) {
    double h = (xk - xn) / n;
    double sum = 0.; 
    for (int i = 0; i < n; i++)
        sum += func((xn + i*h), a, b, c, d);
    return h * sum;
}
//right rule
double RightRectangleMethod(double xn, double xk, double a, double b, double c, double d, int n) {
    double h = (xk - xn) / n;
    double sum = 0.; 
    for (int i = 0; i < n; i++)
        sum += func((xn + i*h + h), a, b, c, d);
    return h * sum;
}
//midpoint rule
double MiddleRectangleMethod(double xn, double xk, double a, double b, double c, double d, int n) {
    double h = (xk - xn) / n;
    double sum = 0.; 
    for (int i = 0; i < n; i++)
        sum += func((xn + i*h + h/2), a, b, c, d);
    return h * sum;
}
//trapezoid rule (identical to the Newton-Cotes rule with the first degree of interpolation)
double TrapezoidalMethod(double xn, double xk, double a, double b, double c, double d, int n) {
    double h = (xk - xn) / n;
    double sum = func(xn, a, b, c, d) + func((xn + h), a, b, c, d);
    for (int i = 1; i < n; i++)
        sum += func((xn + i*h), a, b, c, d) + func((xn + i*h + h), a, b, c, d);
    return sum*h/2;
}
//Newton-Cotes rule (Simpson's rule is implemented due to the second degree of interpolation)
double KotesMethod(double xn, double xk, double a, double b, double c, double d, int n, int m){
    int koef[9][9] = {  {1, 0, 0, 0, 0, 0, 0, 0, 0},
					 	{1, 1, 0, 0, 0, 0, 0, 0, 0},
						{1, 4, 1, 0, 0, 0, 0, 0, 0},
						{1, 3, 3, 1, 0, 0, 0, 0, 0},
						{7, 32, 12, 32, 7, 0, 0, 0, 0},
						{19, 75, 50, 50, 75, 19, 0, 0, 0},
						{41, 216, 27, 272, 27, 216, 41, 0, 0},
						{751, 3577, 1323, 2989, 2989, 1323, 3577, 751, 0},
						{989, 5888, -928, 10496, -4540, 10496, -928, 5888, 989}};
	double sk[9] = {1, 2, 6, 8, 90, 288, 840, 17280, 28350};
	if ((m < 0) || (m > 9)) std::cerr << "WARNING\nIncorrect interpolation degree value in Kotes's method\nPlease double check the value\n\n";
	double h = (xk - xn)/n;
	double sum = 0.;
	for (int i = 0; i < n; i++){
       double x = xn + i*h;
        for (int j = 0; j <= m; j++)
		    if (m >= 1) sum += koef[m][j]*func((x + j*h/m), a, b, c, d);
            else sum += koef[m][j]*func((x + j*h), a, b, c, d); 
	}
	return h*sum/sk[m];
}
int main() {
    double xn = -2.4, xk = 0.0, a = 2.14, b = 5.14, c = -0.14, d = 0.34; //Integration limits and coefficients of the integrand
    const int n = 50; const int m = 6; //n - number of splits; m - degree of interpolation
    double abs1, abs2, abs3, abs4, abs5;
    double meth1 = LeftRectangleMethod(xn, xk, a, b, c, d, n);
    double meth2 = RightRectangleMethod(xn, xk, a, b, c, d, n);
    double meth3 = MiddleRectangleMethod(xn, xk, a, b, c, d, n);
    double meth4 = TrapezoidalMethod(xn, xk, a, b, c, d, n);
    double meth5 = KotesMethod(xn, xk, a, b, c, d, n, 2);
    double meth6 = KotesMethod(xn, xk, a, b, c, d, n, m);
    abs1 = fabs(meth6 - meth1); //calculation of the absolute error of the method, taking into account that the Newton-Cotes rule is the most accurate
    abs2 = fabs(meth6 - meth2);
    abs3 = fabs(meth6 - meth3);
    abs4 = fabs(meth6 - meth4);
    abs5 = fabs(meth6 - meth5);
    ///*output to the terminal
    std::cout << "Results\n"
    << "Left rectangle method result: " << std::setprecision(5) << meth1 << std::endl
    << "Right rectangle method: " << std::setprecision(5) << meth2 << std::endl
    << "Middle rectangle method: " << std::setprecision(5) << meth3 << std::endl
    << "Trapezoidal method result: " << std::setprecision(5) << meth4 << std::endl
    << "Simpson's method: " << std::setprecision(5) << meth5 << std::endl
    << "Kotes's method result: " << std::setprecision(5) << meth6 << std::endl
    << "\nAbsolute error\n"
    << "Left rectangle method: " << std::scientific << std::setprecision(3) << abs1 << std::endl
    << "Right rectangle method: " << std::scientific <<  std::setprecision(3) << abs2 << std::endl
    << "Middle rectangle method: " << std::scientific <<  std::setprecision(3) << abs3 << std::endl
    << "Trapezoidal method: " << std::scientific <<  std::setprecision(3) << abs4 << std::endl
    << "Simpson's method: " << std::scientific <<  std::setprecision(3) << abs5 << std::endl
    ;//*/
    //
    /*output to HTML
    std::ofstream outFile("D:\\Temp\\test1.txt"); //create a file in a folder
    if (!outFile) {
        std::cout << "Error" << std::endl;
        exit(1);
    }
    //tabel of results
    if ((m < 0) || (m > 9)) {outFile << "Incorrect interpolation degree value in Kotes's method"
    << "</br> Please double check the value";}
    else {
    outFile << "<style> table, th, td {border: 1px solid black;} </style>";
    outFile << "<table>\n"
        << "<thead>\n"
        << "<tr>\n"
        << "<th colspan=\"2\"> <b>Результаты расчета определенного интеграла</b> </th>\n"
        << "</tr>\n"
        << "<tr>\n"
        << "<th> <b>Метод</b> </th>\n"
        << "<th> <b>Результат</b> </th>\n"
        << "</tr>\n"
        << "</thead>\n"
        << "<tbody>\n";
        outFile << "<tr>\n";
        outFile << "<td>" << "Левых прямогульников" << "</td>\n";
        outFile << "<td>" << std::setprecision(6) << meth1 << "</td>\n";
        outFile << "</tr>\n";
        outFile << "<td>" << "Правых прямогульников" << "</td>\n";
        outFile << "<td>" << std::setprecision(6) << meth2 << "</td>\n";
        outFile << "</tr>\n";
        outFile << "<td>" << "Средних прямогульников" << "</td>\n";
        outFile << "<td>" << std::setprecision(6) << meth3 << "</td>\n";
        outFile << "</tr>\n";
        outFile << "<td>" << "Трапеций" << "</td>\n";
        outFile << "<td>" << std::setprecision(6) << meth4 << "</td>\n";
        outFile << "</tr>\n";
        outFile << "<td>" << "Симпсона" << "</td>\n";
        outFile << "<td>" << std::setprecision(6) << meth5 << "</td>\n";
        outFile << "</tr>\n";
        outFile << "<td>" << "Котеса 6 порядка" << "</td>\n";
        outFile << "<td>" << std::setprecision(6) << meth6 << "</td>\n";
        outFile << "</tr>\n";
    outFile << "</tbody>\n" << "</table>\n" << "</body>\n"<< "</br>\n" << "</html>";
    //tabel of errors
    outFile << "<style> table, th, td {border: 1px solid black;} </style>";
    outFile << "<table>\n"
        << "<thead>\n"
        << "<tr>\n"
        << "<th colspan=\"2\"> <b>Абсолютная ошибка метода</b> </th>\n"
        << "</tr>\n"
        << "<tr>\n"
        << "<th> <b>Метод</b> </th>\n"
        << "<th> <b>Абсолютная ошибка</b> </th>\n"
        << "</tr>\n"
        << "</thead>\n"
        << "<tbody>\n";
        outFile << "<tr>\n";
        outFile << "<td>" << "Левых прямогульников" << "</td>\n";
        outFile << "<td>" << std::scientific << std::setprecision(3) << abs1 << "</td>\n";
        outFile << "</tr>\n";
        outFile << "<td>" << "Парвых прямогульников" << "</td>\n";
        outFile << "<td>" << std::scientific << std::setprecision(3) << abs2 << "</td>\n";
        outFile << "</tr>\n";
        outFile << "<td>" << "Средних прямогульников" << "</td>\n";
        outFile << "<td>" << std::scientific << std::setprecision(3) << abs3 << "</td>\n";
        outFile << "</tr>\n";
        outFile << "<td>" << "Трапеций" << "</td>\n";
        outFile << "<td>" << std::scientific << std::setprecision(3) << abs4 << "</td>\n";
        outFile << "</tr>\n";
        outFile << "<td>" << "Симпсона" << "</td>\n";
        outFile << "<td>" << std::scientific << std::setprecision(3) << abs5 << "</td>\n";
        outFile << "</tr>\n";
    outFile << "</tbody>\n" << "</table>\n" << "</body>\n" << "</html>";}//
    //
    */
    return 0;
} 
