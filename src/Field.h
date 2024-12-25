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

