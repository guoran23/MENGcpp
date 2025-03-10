/*
 * Author: Guo Meng
 * Email: guo.meng@ipp.mpg.de
 * Created Date: 2024-12-10
 * Last Modified: 2025-03-10
 * License: MIT License
 *
 * Description:
 * This file is part of MENG++ project.
 */

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

