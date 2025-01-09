#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

class Equilibrium
{
private:
    double pressure;
    double temperature;
    double Btor;
    double aminor;
    double qbar_axis;
    double qbar_edge;

public:
    Equilibrium(double p, double t);
    double getPressure() const;
    double getTemperature() const;
    void adjustEquilibrium(double dPressure, double dTemperature);
    void displayEquilibrium() const;
    double getBR()
    {
    }
    double get_qbar(double rho)
    {
        double q2 = (qbar_edge - qbar_axis) / aminor / aminor;
        double qbar = qbar_axis + q2 * rho * rho;
        return qbar;
    }
};

#endif // EQUILIBRIUM_H
