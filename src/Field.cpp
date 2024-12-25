#include "Field.h"

Field::Field(double strength) : fieldStrength(strength) {}

double Field::getStrength() const {
    return fieldStrength;
}

void Field::setStrength(double strength) {
    fieldStrength = strength;
}

double Field::calculateForce(double charge) const {
    return charge * fieldStrength;
}

