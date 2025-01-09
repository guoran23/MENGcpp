#include <iostream>
#include "Field.h"
#include "Particle.h"
#include "Equilibrium.h"

int main() {
    // Create a magnetic field
    Field field(10.0);

    // Create a particle
    Particle particle(1.0, 1.0, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});

    // Create an equilibrium state
    Equilibrium eq(101325.0, 300.0);

    // Simulate particle motion
    std::cout << "Initial State:\n";
    particle.displayState();
    eq.displayEquilibrium();
    // ---------------------------------------------------------
    eq.readInput("input.txt");
    eq.printState();
    std::cout<<"qbar(0)="<<eq.get_qbar(0.0)<<std::endl;

    double timeStep = 0.1;
    for (int i = 0; i < 10; ++i) {
        particle.move(field, timeStep);
    }

    std::cout << "\nFinal State:\n";
    particle.displayState();
    eq.displayEquilibrium();

    return 0;
}

