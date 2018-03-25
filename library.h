#ifndef MATHLIB_LIBRARY_H
#define MATHLIB_LIBRARY_H

#include <cmath>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace Math
{
    using std::sin;
    using std::cos;
    using std::tan;
    using std::asin;
    using std::acos;
    using std::atan;
    using std::atan2;
    using std::sinh;
    using std::cosh;
    using std::tanh;
    using std::asinh;
    using std::acosh;
    using std::atanh;
    using std::exp;
    using std::logb;
    using std::log2;
    using std::log10;
    using std::pow;
    using std::sqrt;
    using std::abs;
    using std::fabs;
    
    struct Point;
    struct Vector;
    struct Line;
    struct Circle;

    struct Point
    {
        double x, y;
        inline explicit Point(double = 0, double = 0);
        Point& operator +=(const Vector&);
    };

    struct Vector
    {
        double x, y;
        inline explicit Vector(double = 0, double = 0);
        Vector& operator +=(const Vector&);
        Vector& operator *=(double);
        Vector& operator /=(double);
    };
    
    struct Line
    {
        Point P;
        Vector v;
        inline explicit Line(const Point& = Point(), const Vector& = Vector());
        const Point GetPoint(double) const;
    };
    
    struct Circle
    {
        Point C;
        double r;
        inline explicit Circle(const Point& = Point(), double = 0);
        const Point GetPoint(double) const;
    };
    
    typedef std::vector<Point> Polygon;

    inline const Vector operator +(Vector, const Vector&);
    inline const Point operator +(Point, const Vector&);
    const Point operator +(const Vector&, const Point&);
    const Vector operator -(const Point&, const Point&);
    inline const Vector operator *(Vector, double);
    inline const Vector operator /(Vector, double);
    
    const bool operator <(const Point&, const Point&);
    const bool operator <(const Vector&, const Vector&);
    const bool operator <(const Line&, const Line&);
    const bool operator ==(const Point&, const Point&);
    const bool operator ==(const Vector&, const Vector&);
    const bool operator ==(const Line&, const Line&);
    const bool operator ==(const Circle&, const Circle&);
    inline const bool operator <=(const Point&, const Point&);
    inline const bool operator <=(const Vector&, const Vector&);
    inline const bool operator <=(const Line&, const Line&);
    inline const bool operator >(const Point&, const Point&);
    inline const bool operator >(const Vector&, const Vector&);
    inline const bool operator >(const Line&, const Line&);
    inline const bool operator >=(const Point&, const Point&);
    inline const bool operator >=(const Vector&, const Vector&);
    inline const bool operator >=(const Line&, const Line&);
    inline const bool operator !=(const Point&, const Point&);
    inline const bool operator !=(const Vector&, const Vector&);
    inline const bool operator !=(const Line&, const Line&);
    inline const bool operator !=(const Circle&, const Circle&);
    
    const double eps = 1e-10;
    const double E = exp(1);
    const double PI = acos(-1);
    
    const int dcmp(double);
    
    inline const double Dot(const Vector&, const Vector&);
    inline const double Length(const Vector&);
    inline const double Angle(const Vector&, const Vector&);
    inline const double Angle(const Vector&);
    inline const double Cross(const Vector&, const Vector&);
    inline const double Area2(const Point&, const Point&, const Point&);
    inline const Vector Rotate(const Vector&, double);
    inline const Vector Normal(const Vector&); // Please ensure Vector is not zero vector
    
    const int Intersection(const Line&, const Line&, std::vector<Point>&);
    const int Intersection(const Line&, const Point&, const Point&, std::vector<Point>&);
    inline const int Intersection(const Point&, const Point&, const Line&, std::vector<Point>&);
    const int Intersection(const Point&, const Point&, const Point&, const Point&, std::vector<Point>&);
    const int Intersection(const Circle&, const Point&, const Point&, std::vector<Point>&);
    inline const int Intersection(const Point&, const Point&, const Circle&, std::vector<Point>&);
    const int Intersection(const Line&, const Circle&, std::vector<Point>&);
    inline const int Intersection(const Circle&, const Line&, std::vector<Point>&);
    const int Intersection(const Circle&, const Circle&, std::vector<Point>&);
    const int Intersection(const Line&, const Polygon&, std::vector<Point>&);
    inline const int Intersection(const Polygon&, const Line&, std::vector<Point>&);
    const int Intersection(const Polygon&, const Point&, const Point&, std::vector<Point>&);
    inline const int Intersection(const Point&, const Point&, const Polygon&, std::vector<Point>&);
    const int Intersection(const Polygon&, const Circle&, std::vector<Point>&);
    inline const int Intersection(const Circle&, const Polygon&, std::vector<Point>&);
    const int Intersection(const Polygon&, const Polygon&, std::vector<Point>&);
    const int Intersection(std::vector<Line>, Polygon&);
    
    const bool ProperIntersection(const Line&, const Line&);
    const bool ProperIntersection(const Point&, const Point&, const Line&);
    inline const bool ProperIntersection(const Line&, const Point&, const Point&);
    const bool ProperIntersection(const Point&, const Point&, const Point&, const Point&);
    
    inline const double Distance(const Point&, const Point&);
    const double Distance(const Point&, const Line&);
    inline const double Distance(const Line&, const Point&);
    const double Distance(const Point&, const Point&, const Point&);
    
    const Point Project(const Point&, const Line&);
    inline const Point Project(const Line&, const Point&);
    
    inline const bool On(const Point&, const Line&);
    inline const bool On(const Line&, const Point&);
    const bool On(const Point&, const Point&, const Point&);
    inline const bool On(const Point&, const Circle&);
    inline const bool On(const Circle&, const Point&);
    const bool On(const Point&, const Polygon&);
    inline const bool On(const Polygon&, const Point&);

    inline const bool In(const Point&, const Line&);
    inline const bool In(const Line&, const Point&);
    inline const bool In(const Point&, const Circle&);
    inline const bool In(const Circle&, const Point&);
    const bool In(const Point&, const Polygon&);
    inline const bool In(const Polygon&, const Point&);
    
    const double Area(const Polygon&);
    inline const double Area(const Circle&);
    
    const int Tangents(const Point&, const Circle&, std::vector<Line>&);
    inline const int Tangents(const Circle&, const Point&, std::vector<Line>&);
    const int Tangents(Circle, Circle, std::vector<Line>&);

    const Circle CircumscribedCircle(const Point&, const Point&, const Point&);
    const Circle InscribedCircle(const Point&, const Point&, const Point&);

    const int ConvexHull(std::vector<Point>, Polygon&);

    const Polygon Cut(const Polygon&, const Line&);
    inline const Polygon Cut(const Line&, const Polygon&);
}

#endif