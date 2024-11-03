#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

struct Point {
    double x;
    double y;
};

struct PrecDerInfo {
    double d1t;
    double d2t;
};

struct UnprecDerInfo {
    double xValue;
    double hParam;
    double firstDer;
    double secondDer;
    double firstDeviation;
    double secondDeviation;
};

class Config {
public:
    int POINTS_NUMBER;
    vector<double> H_PARAMS;
    vector<double> PRECISION_VALUES;
    string PLOTS_JSON_FILE;
    double a;
    double b;

    Config(const string& configFilePath) {
        ifstream configFile(configFilePath.c_str());
        if (!configFile) {
            cerr << "Не удалось открыть конфигурационный файл: " << configFilePath << endl;
            exit(EXIT_FAILURE);
        }

        json configJson;
        configFile >> configJson;
        configFile.close();

        POINTS_NUMBER = configJson.value("POINTS_NUMBER", 21);
        H_PARAMS = configJson["H_PARAMS"].get<vector<double>>();
        PRECISION_VALUES = configJson["PRECISION_VALUES"].get<vector<double>>();
        PLOTS_JSON_FILE = configJson.value("PLOTS_JSON_FILE", "plots.json");
        a = configJson.value("A", -2.0);
        b = configJson.value("B", 5.0);
    }
};

class Calculator {
public:
    Config config;
    vector<Point> points;
    vector<PrecDerInfo> precDerivationInfo;
    vector<UnprecDerInfo> unprecDerivations;
    json resultJson;

    Calculator(const Config& cfg) : config(cfg){
    }

    double GetFunctionResult(double x) {
        return pow(x, 3) - 5.0 * pow(x, 2);
    }

    double GetPrecFirstDer(double x) {
        return 3.0 * x * x - 10.0 * x;
    }

    double GetPrecSecondDer(double x) {
        return 6.0 * x - 10.0;
    }

    double GetPrecisionIntegral(double a, double b) {
        double Fa = pow(a, 4) / 4.0 - (5.0 * pow(a, 3)) / 3.0;
        double Fb = pow(b, 4) / 4.0 - (5.0 * pow(b, 3)) / 3.0;
        return Fb - Fa;
    }

    double GetFirstDerivative(double x, double h) {
        return (GetFunctionResult(x + h) - GetFunctionResult(x - h)) / (2.0 * h);
    }

    double GetFirstDerForFirstElem(double x, double h) {
        return (-3.0 * GetFunctionResult(x - h) + 4.0 * GetFunctionResult(x) - GetFunctionResult(x + h)) / (2.0 * h);
    }

    double GetFirstDerForLastElem(double x, double h) {
        return (3.0 * GetFunctionResult(x + h) - 4.0 * GetFunctionResult(x) + GetFunctionResult(x - h)) / (2.0 * h);
    }

    double GetSecondDerivative(double x, double h) {
        return (GetFunctionResult(x - h) - 2.0 * GetFunctionResult(x) + GetFunctionResult(x + h)) / (h * h);
    }

    double GetIntegral(double a, double b, double precision, double& m) {
        double h = (b - a) / 2.0;
        double h2 = h / 2.0;
        double st = GetFunctionResult(a + h) + (GetFunctionResult(a) + GetFunctionResult(b)) / 2.0;
        double ss = GetFunctionResult(a + h2) + GetFunctionResult(b - h2);
        m = 2.0;

        while (h * fabs(ss - st) > precision) {
            st += ss;
            ss = 0.0;
            h /= 2.0;
            h2 /= 2.0;
            m *= 2.0;
            for (double x = a + h2; x < b; x += h) {
                ss += GetFunctionResult(x);
            }
        }

        return h * (st + 2.0 * ss) / 3.0;
    }

    string formatDouble(double d) {
        string s = to_string(d);
        size_t dot_pos = s.find('.');
        if (dot_pos != string::npos) {
            size_t last_non_zero = s.find_last_not_of('0');
            if (last_non_zero != string::npos && last_non_zero > dot_pos) {
                s.erase(last_non_zero + 1);
            }
            if (s.back() == '.') {
                s.pop_back();
            }
        }
        return s;
    }

    void computeDerivatives() {
        for (int i = 0; i < config.POINTS_NUMBER; i++) {
            Point p;
            p.x = config.a + i * (config.b - config.a) / 20.0;
            p.y = GetFunctionResult(p.x);
            points.push_back(p);

            PrecDerInfo precInfo;
            precInfo.d1t = GetPrecFirstDer(p.x);
            precInfo.d2t = GetPrecSecondDer(p.x);
            precDerivationInfo.push_back(precInfo);

            cout << "Point " << i << ": x = " << p.x << endl;

            for (size_t j = 0; j < config.H_PARAMS.size(); j++) {
                double h = config.H_PARAMS[j];
                UnprecDerInfo unprecDerv;
                unprecDerv.xValue = p.x;
                unprecDerv.hParam = h;

                if (i == 0) {
                    unprecDerv.firstDer = GetFirstDerForFirstElem(p.x, h);
                    unprecDerv.secondDer = NAN;
                    unprecDerv.secondDeviation = NAN;
                }
                else if (i == config.POINTS_NUMBER - 1) {
                    unprecDerv.firstDer = GetFirstDerForLastElem(p.x, h);
                    unprecDerv.secondDer = NAN;
                    unprecDerv.secondDeviation = NAN;
                }
                else {
                    unprecDerv.firstDer = GetFirstDerivative(p.x, h);
                    unprecDerv.secondDer = GetSecondDerivative(p.x, h);
                    unprecDerv.secondDeviation = precInfo.d2t - unprecDerv.secondDer;
                }

                unprecDerv.firstDeviation = precInfo.d1t - unprecDerv.firstDer;
                unprecDerivations.push_back(unprecDerv);

                cout << "  H_PARAM = " << h << endl;
                cout << "    First Derivative Deviation: " << unprecDerv.firstDeviation << endl;
                if (!isnan(unprecDerv.secondDeviation)) {
                    cout << "    Second Derivative Deviation: " << unprecDerv.secondDeviation << endl;
                }
            }
        }
    }

    void computeIntegrals() {
        double precIntegral = GetPrecisionIntegral(config.a, config.b);
        cout << "\nExact Integral Value: " << precIntegral << endl;

        for (size_t i = 0; i < config.PRECISION_VALUES.size(); i++) {
            double precision = config.PRECISION_VALUES[i];
            double m = 0.0;
            double integralResult = GetIntegral(config.a, config.b, precision, m);
            json integral_entry;
            integral_entry["precision"] = precision;
            integral_entry["steps"] = m;
            resultJson["Integral"].push_back(integral_entry);
            resultJson["IntegralDeviation"].push_back(precIntegral - integralResult);

            cout << "Precision: " << precision
                << ", Integral Deviation: " << (precIntegral - integralResult)
                << ", Steps: " << m << endl;
        }
    }

    void prepareJSON() {
        for (size_t i = 0; i < config.H_PARAMS.size(); i++) {
            double h = config.H_PARAMS[i];
            string h_str = formatDouble(h);

            vector<json> firstDerEntries;
            vector<json> secondDerEntries;

            for (size_t j = 0; j < unprecDerivations.size(); j++) {
                if (unprecDerivations[j].hParam == h) {
                    json first_der_entry;
                    first_der_entry["x"] = unprecDerivations[j].xValue;
                    first_der_entry["deviation"] = unprecDerivations[j].firstDeviation;
                    firstDerEntries.push_back(first_der_entry);

                    if (!isnan(unprecDerivations[j].secondDeviation)) {
                        json second_der_entry;
                        second_der_entry["x"] = unprecDerivations[j].xValue;
                        second_der_entry["deviation"] = unprecDerivations[j].secondDeviation;
                        secondDerEntries.push_back(second_der_entry);
                    }
                }
            }

            resultJson["First derivative"][h_str] = firstDerEntries;
            resultJson["Second derivative"][h_str] = secondDerEntries;
        }

        resultJson["H_PARAMS"] = config.H_PARAMS;
        resultJson["PRECISION_VALUES"] = config.PRECISION_VALUES;
    }

    void saveResults() {
        ofstream file(config.PLOTS_JSON_FILE.c_str());
        if (file) {
            file << resultJson.dump(4);
            cout << "\nData saved to " << config.PLOTS_JSON_FILE << endl;
        }
        else {
            cerr << "Не удалось открыть файл " << config.PLOTS_JSON_FILE << " для записи." << endl;
            exit(EXIT_FAILURE);
        }
    }

    void run() {
        computeDerivatives();
        computeIntegrals();
        prepareJSON();
        saveResults();
        system("python plots.py");
    }
};

int main() {
    Config config("config.json");
    cout << "Using A = " << config.a << ", B = " << config.b << endl;

    Calculator calculator(config);
    calculator.run();

    return 0;
}
