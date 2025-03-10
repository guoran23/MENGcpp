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

#ifndef FIELD_H
#define FIELD_H

class Field {
private:
    double fieldStrength;

public:
    Field(double strength = 1.0);
    double getStrength() const;
    void setStrength(double strength);
    double calculateForce(double charge) const;
};

#endif // FIELD_H

