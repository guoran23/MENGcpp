#include "Particle.h"
#include <iostream>

Particle::Particle(double m, double q, std::vector<double> pos, std::vector<double> vel)
    : mass(m), charge(q), position(pos), velocity(vel) {}

void Particle::move(const Field& field, double dt) {
    double force = field.calculateForce(charge);
    double acceleration = force / mass;

    // Update velocity and position in 1D for simplicity
    velocity[0] += acceleration * dt;
    position[0] += velocity[0] * dt;
}

void Particle::displayState() const {
    std::cout << "Position: (" << position[0] << ", " << position[1] << ", " << position[2] << ")"
              << "\nVelocity: (" << velocity[0] << ", " << velocity[1] << ", " << velocity[2] << ")\n";
}

