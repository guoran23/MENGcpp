#include "Equilibrium.h"
#include <iostream>

Equilibrium::Equilibrium(double p, double t) : pressure(p), temperature(t) {}

double Equilibrium::getPressure() const {
    return pressure;
}

double Equilibrium::getTemperature() const {
    return temperature;
}

void Equilibrium::adjustEquilibrium(double dPressure, double dTemperature) {
    pressure += dPressure;
    temperature += dTemperature;
}

void Equilibrium::displayEquilibrium() const {
    std::cout << "Pressure: " << pressure << " Pa, Temperature: " << temperature << " K\n";
}

