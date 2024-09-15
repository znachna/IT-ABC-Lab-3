#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>
#include <nlohmann/json.hpp>

using namespace std;

using json = nlohmann::json;

struct Point {
    double x;
    double y;

    Point(double serialNumber, double a, double b, double (*Func)(double)) {
        x = a + (serialNumber-1) * (b - a) / 20;
        y = Func(x);
    }
};

struct PrecDerInfo {

    double d1t;
    double d2t;

    PrecDerInfo(double (*FirstDer)(double), double (*SecondDer)(double), double x) {
        d1t = FirstDer(x);
        d2t = SecondDer(x);
    }
};

struct UnprecDerInfo {

    double xValue;

    double hParam;

    double firstDer;
    double secondDer;

    double firstDeviation;
    double secondDeviation;
};

double GetPrecFirstDer(double x) {
    return 3 * pow(x, 2) - 10 * x;
}

double GetPrecSecondDer(double x) {
    return 6 * x - 10;
}

double GetPrecisionIntergral(double a, double b) {
    double Fa = pow(a, 4) / 4 - 5 * pow(a, 3) / 3;
    double Fb = pow(b, 4) / 4 - 5 * pow(b, 3) / 3;
    return Fb - Fa;
}

double GetFunctionResult(double x) {
    return pow(x, 3) - 5 * pow(x, 2);
}

double GetFirstDerivative(double x, double (*Func)(double), double h) {
    double y3 = Func(x + h);
    double y1 = Func(x - h);

    return (y3 - y1) / (2 * h);
}

double GetFirstDerForFirstElem(double x, double (*Func)(double), double h) {
    double y3 = Func(x + h);
    double y2 = Func(x);
    double y1 = Func(x - h);

    double result = (-1 * (3 * y1 - 4 * y2 + y3) ) / (2 * h);

    return result;
}

double GetFirstDerForLastElem(double x, double (*Func)(double), double h) {
    double y3 = Func(x + h);
    double y2 = Func(x);
    double y1 = Func(x - h);

    double result = (y1 - 4 * y2 + 3*y3) / (2 * h);

    return result;
}


double GetSecondDerivative(double x, double (*Func)(double), double h) {
    double y1 = Func(x-h);
    double y2 = Func(x);
    double y3 = Func(x + h);

    return (y1 - 2 * y2 + y3) / pow(h, 2);
}

double GetIntegral(double a, double b, double precision, double (*Func)(double), double& m) {

    m = 2; 

    double h = (b - a) / 2; 
    double h2 = h / 2.0;

    double st = Func(a + h) + (Func(a) + Func(b)) / 2.0; 

    double ss = Func(b - h2) + Func(a + h2);

    double counter = 0;

    while (((h * fabs(ss - st)) > precision)) {

        counter++;

        st = st + ss;
        ss = 0;
        h = h2;
        h2 = h2 / 2.0;
        m = 2 * m;

        double x = a + h2;

        do {
            ss = ss + Func(x);
            x = x + h;
        } while (x < b);
    }

    double result = h * (st + 2 * ss) / 3.0;
    return result;
}


int main()
{
    json j;

    const int VARIANTS_FOR_TESTING = 3;
    const int POINTS_NUMBER = 21;
    const vector <double> H_PARAMS{ 0.2, 0.1, 0.05 };

    double a, b;
    cout << "Print variables A and B: " << endl;
    cin >> a >> b;

    vector <Point> points;
    vector <PrecDerInfo> precDerevationInfo;

    vector <UnprecDerInfo> unprecDerivations;

    for (int i = 0; i < POINTS_NUMBER; i++) {
        points.push_back(Point(i, a, b, GetFunctionResult));

        PrecDerInfo precInfo(GetPrecFirstDer, GetPrecSecondDer, points.back().x);
        precDerevationInfo.push_back(precInfo);

        cout << "X_value: " << points.back().x << endl;


        for (auto H_PARAM : H_PARAMS) {
            UnprecDerInfo unprecDerv;

            unprecDerv.xValue = points.back().x;
            unprecDerv.hParam = H_PARAM;
            
            cout << "H_PARAM: " << H_PARAM << endl;

            switch (i) {
                case 0: 
                    unprecDerv.firstDer = GetFirstDerForFirstElem(points.back().x, GetFunctionResult, H_PARAM);

                    break;
                case POINTS_NUMBER - 1:
                    unprecDerv.firstDer = GetFirstDerForLastElem(points.back().x, GetFunctionResult, H_PARAM);

                    break;
                default:
                    unprecDerv.secondDer = GetSecondDerivative(points.back().x, GetFunctionResult, H_PARAM);
                    unprecDerv.secondDeviation = precInfo.d2t - unprecDerv.secondDer;
                    cout << "Second der deviation: " << unprecDerv.secondDeviation << endl;
                    unprecDerv.firstDer = GetFirstDerivative(points.back().x, GetFunctionResult, H_PARAM);

                    break;
            }
            unprecDerv.firstDeviation = precInfo.d1t - unprecDerv.firstDer;
            unprecDerivations.push_back(unprecDerv);

            cout << "First der deviation: " << unprecDerv.firstDeviation << endl;
        }
        cout << endl;
    }

    vector <double> unprecIntergrals;
    const vector <double> PRECISION_VALUES{ 0.1, 0.01, 0.001 };
    vector <vector <double> > intergralSteps;
    vector <double> integralDevitation;

    double precIntergral = GetPrecisionIntergral(a, b);

    for (auto const PRECISION_VALUE : PRECISION_VALUES) {
        double m;
        unprecIntergrals.push_back(GetIntegral(a, b, PRECISION_VALUE, GetFunctionResult, m));
        intergralSteps.push_back(vector<double>{PRECISION_VALUE, m});
        integralDevitation.push_back(precIntergral - unprecIntergrals.back());

        cout << endl << "PRECISION_VALUE: " << PRECISION_VALUE << endl;
        cout << "Integral Devitation: " << integralDevitation.back() << endl;
        cout << "Steps: " << m << endl;

        j["Intergral deviation"].push_back({ {"precision", PRECISION_VALUE}, {"steps", m} });
    }

    for (const auto H_PARAM : H_PARAMS) {

        for (const auto& unprecDerv : unprecDerivations) {
            if (unprecDerv.hParam == H_PARAM) {
                j["First derivative"][to_string(H_PARAM)].push_back({ {"x",unprecDerv.xValue},{"deviation",unprecDerv.firstDeviation} });
                j["Second derivative"][to_string(H_PARAM)].push_back({ {"x",unprecDerv.xValue},{"deviation",unprecDerv.secondDeviation} });
            }
        }
    }

    ofstream file("plots.json");
    file << j.dump(1);
    file.close();

    system("python plots.py");
}