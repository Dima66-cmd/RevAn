//
// Created by Dima on 08.10.2021.
//

#include <iostream>
#include <cmath>
#include <vector>

struct Point {
    double x;
    double y;

    Point(double x, double y) : x(x), y(y) {}

    Point operator+(const Point& p) const {
        return {x + p.x, y + p.y};
    }

    Point operator-(const Point& p) const {
        return {x - p.x, y - p.y};
    }

    Point operator*(const double& k) const {
        return {x * k, y * k};
    }
};

class CassiniLine {
public:
    CassiniLine(const Point& f1, const Point& f2, const double& a) : f1_(f1), f2_(f2), a_(a) {}

    void SetF1(const Point& f1) {
        f1_=f1;
    }

    void SetF2(const Point& f2) {
        f2_ = f2;
    }

    void SetA(const double& a) {
        a_ = a;
    }

    Point GetF1() const {
        return f1_;
    }

    Point GetF2() const {
        return f2_;
    }

    double GetA() const {
        return a_;
    }

    double GetDistancePolar(const double& angle) const {
        double phi_1 = PolarPoint(f1_).phi;
        double phi_2 = PolarPoint(f2_).phi;

        // 4 * \rho^4 * (1 - cos(\phi - \phi_1) * (1 - cos(\phi - \phi_2)) = a^4
        return pow((pow(a_, 4) / (4 * (1 - cos(angle - phi_1)) * (1 - cos(angle - phi_2)))), 0.25);
    }

    double GetCurvature(const double& angle) const {
        // https://ru.wikipedia.org/wiki/Овал_Кассини#Свойства

        double rho = GetDistancePolar(angle);
        double c = GetC();

        return (2 * pow(a_, 2) * pow(rho, 3)) / (pow(c, 4) - pow(a_, 4) - 3 * pow(rho, 3));
    }

    std::vector<Point> GetInflectionPoints() const {
        // https://ru.wikipedia.org/wiki/Овал_Кассини#Свойства

        std::vector<Point> res;
        double c = GetC();

        if (!(c < a_ && a_ < c * sqrt(2))) {
            return res;
        }

        //  Shift from (-c, 0)(c, 0) to F1F2
        Point center = (f2_ + f1_) * 0.5;
        double rho = pow((pow(a_, 4) - pow(c, 4)) / 3, 0.25) - Norm(center);
        double x = -sqrt((pow(a_, 4) / pow(c, 4) - 1) / 3);

        // cos(2 * \phi) = x
        std::vector<double> phis;
        phis.emplace_back(acos(x) / 2);
        phis.emplace_back(2 * M_PI - acos(x) / 2);
        phis.emplace_back(acos(x) / 2 + M_PI_2);
        phis.emplace_back(2 * M_PI - acos(x) / 2 - M_PI_2);


        //  Shift from (-c, 0)(c, 0) to F1F2
        Point delta = f2_ - f1_;
        for (auto& phi: phis) {
            phi -= PolarPoint(delta).phi;
            res.emplace_back(PolarPoint(rho, phi).ToPoint());
        }

        return res;
    }

    std::string GetType() const {
        // https://ru.wikipedia.org/wiki/Овал_Кассини#Особенности_формы

        double c = GetC();
        if (a_ == 0 && c != 0) {
            return "Две точки";
        }

        if (0 < a_ && a_ < c) {
            return "Два овала";
        }

        if (a_ == c) {
            return "Лемниската";
        }

        if (c < a_ && a_ < c * sqrt(2)) {
            return "Лемниската с перегибами";
        }

        if (a_ >= c * sqrt(2)) {
            return "Овал";
        }

        return "Окружность";
    }

    char* GetEquation() const {
        char* res = new char[255];
        snprintf(res,
                 255,
                 "((x - %f)^2 + (y - %f)^2) * ((x - %f)^2 + (y - %f)^2) = %f ^ 4",
                 f1_.x,
                 f1_.y,
                 f2_.x,
                 f2_.y,
                 a_);
        return res;
    }

private:
    struct PolarPoint {
        PolarPoint(double rho, double phi) : rho(rho), phi(phi) {}

        explicit PolarPoint(const Point& p) {
            rho = Norm(p);


            // Edge cases
            if (p.x == 0 && p.y == 0) {
                phi = 0;
            } else if (p.y == 0) {
                phi = p.x > 0 ? 0 : M_PI;
            } else if (p.x == 0) {
                phi = p.y > 0 ? M_PI_2 : 3 * M_PI_2;
            } else {
                // Main case
                double abs_phi = atan(abs(p.y) / abs(p.x));
                if (p.x < 0) {
                    abs_phi = M_PI - abs_phi;
                }

                if (p.y < 0) {
                    abs_phi = 2 * M_PI - abs_phi;
                }

                phi = abs_phi;
            }
        }

        Point ToPoint() const {
            return {rho * cos(phi), rho * sin(phi)};
        }

        double rho;

        // [0, 2 * pi]
        double phi;
    };

    static double Distance(const Point& p1, const Point& p2) {
        return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
    }

    static double Norm(const Point& p) {
        return Distance(p, {0, 0});
    }

    double GetC() const {
        return Distance(f1_, f2_) / 2;
    }

    Point f1_;
    Point f2_;
    double a_;
};

void Print(const std::string& s) {
    std::cout << s << std::endl;
}

void Test(Point p1, Point p2, double a, double angle) {
    CassiniLine l(p1, p2, a);
    auto s1 = l.GetEquation();
    auto s2 = l.GetType();

    for (const auto& p: l.GetInflectionPoints()) {
    }

    auto d1 = l.GetDistancePolar(angle);
    auto d2 = l.GetCurvature(angle);
}

int test() {
    const int RANGE = 3;
    for (int x1 = -RANGE; x1 <= RANGE; ++x1) {
        for (int y1 = -RANGE; y1 <= RANGE; ++y1) {
            for (int x2 = -RANGE; x2 <= RANGE; ++x2) {
                for (int y2 = -RANGE; y2 <= RANGE; ++y2) {
                    for (int a = -RANGE; a <= RANGE; ++a) {
                        for (int angle = 0; angle <= 360; ++angle) {
                            Test(
                                    {(double) x1 / RANGE, (double) y1 / RANGE},
                                    {(double) x2 / RANGE, (double) y2 / RANGE},
                                    (double) a / RANGE,
                                    (double) angle / 180 * M_PI
                            );
                        }
                    }
                }
            }
        }
    }
    return 0;
}
int main()
{

}