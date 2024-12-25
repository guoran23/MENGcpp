#ifndef PARTICLE_H
#define PARTICLE_H

#include "Field.h"
#include <vector>

class Particle {
private:
    double mass;
    double charge;
    std::vector<double> position;
    std::vector<double> velocity;

public:
    Particle(double m, double q, std::vector<double> pos, std::vector<double> vel);
    void move(const Field& field, double dt);
    void displayState() const;
};

#endif // PARTICLE_H

