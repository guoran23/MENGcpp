#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

class Equilibrium {
private:
    double pressure;
    double temperature;

public:
    Equilibrium(double p, double t);
    double getPressure() const;
    double getTemperature() const;
    void adjustEquilibrium(double dPressure, double dTemperature);
    void displayEquilibrium() const;
};

#endif // EQUILIBRIUM_H

