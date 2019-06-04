#ifndef _POINT_H
#define _POINT_H

// Class point used to define the Schrodinger Potential.
// The Schrodinger Potential is a vector of Points (x, y)
class Point{
    private:
        double x, y;
    public:
        // Default class constructor
        Point();
        // Class constructor
        Point(const double xx, const double yy);
        // Setter of x
        void setX(const double xx);
        // Setter of y
        void setY(const double yy);
        // Getter of x
        double getX();
        // Getter of y
        double getY();
};

;
#endif