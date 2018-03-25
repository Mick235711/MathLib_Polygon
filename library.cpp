#include "library.h"

namespace Math
{
    inline Point::Point(double _x, double _y)
        : x(_x), y(_y)
    {}

    inline Vector::Vector(double _x, double _y)
        : x(_x), y(_y)
    {}

    inline Line::Line(const Point& _P, const Vector & _v)
        : P(_P), v(_v)
    {}

    inline Circle::Circle(const Point& _C, double _r)
        : C(_C), r(_r)
    {}

    Point& Point::operator +=(const Vector& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    Vector& Vector::operator +=(const Vector& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    Vector& Vector::operator *=(double rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    Vector& Vector::operator /=(double rhs)
    {
        x /= rhs;
        y /= rhs;
        return *this;
    }

    const Point operator +(const Vector& lhs, const Point& rhs)
    {
        return Point(lhs.x + rhs.x, lhs.y + rhs.y);
    }

    const Vector operator -(const Point& lhs, const Point& rhs)
    {
        return Vector(lhs.x - rhs.x, lhs.y + rhs.y);
    }
    
    const bool operator <(const Point& lhs, const Point& rhs)
    {
        return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y);
    }
    
    const bool operator ==(const Point& lhs, const Point& rhs)
    {
        return dcmp(lhs.x - rhs.x) == 0 && dcmp(lhs.y - rhs.y) == 0;
    }

    const bool operator <(const Vector& lhs, const Vector& rhs)
    {
        return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y);
    }

    const bool operator ==(const Vector& lhs, const Vector& rhs)
    {
        return dcmp(lhs.x - rhs.x) == 0 && dcmp(lhs.y - rhs.y) == 0;
    }

    const bool operator <(const Line& lhs, const Line& rhs)
    {
        return dcmp(Angle(lhs.v) - Angle(rhs.v)) < 0;
    }

    const bool operator ==(const Line& lhs, const Line& rhs)
    {
        return dcmp(Cross(lhs.v, rhs.v)) == 0 && On(lhs.P, rhs);
    }

    const bool operator ==(const Circle& lhs, const Circle& rhs)
    {
        return lhs.C == rhs.C && dcmp(lhs.r - rhs.r) == 0;
    }

    inline const Vector operator +(Vector lhs, const Vector& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    inline const Point operator +(Point lhs, const Vector& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    inline const Vector operator *(Vector lhs, double rhs)
    {
        lhs *= rhs;
        return lhs;
    }

    inline const Vector operator /(Vector lhs, double rhs)
    {
        lhs /= rhs;
        return lhs;
    }

    inline const bool operator <=(const Point& lhs, const Point& rhs)
    {
        return !(rhs < lhs);
    }

    inline const bool operator >(const Point& lhs, const Point& rhs)
    {
        return rhs < lhs;
    }

    inline const bool operator >=(const Point& lhs, const Point& rhs)
    {
        return !(lhs < rhs);
    }

    inline const bool operator !=(const Point& lhs, const Point& rhs)
    {
        return !(lhs == rhs);
    }

    inline const bool operator <=(const Vector& lhs, const Vector& rhs)
    {
        return !(rhs < lhs);
    }

    inline const bool operator >(const Vector& lhs, const Vector& rhs)
    {
        return rhs < lhs;
    }

    inline const bool operator >=(const Vector& lhs, const Vector& rhs)
    {
        return !(lhs < rhs);
    }

    inline const bool operator !=(const Vector& lhs, const Vector& rhs)
    {
        return !(lhs == rhs);
    }

    inline const bool operator <=(const Line& lhs, const Line& rhs)
    {
        return !(rhs < lhs);
    }

    inline const bool operator >(const Line& lhs, const Line& rhs)
    {
        return rhs < lhs;
    }

    inline const bool operator >=(const Line& lhs, const Line& rhs)
    {
        return !(lhs < rhs);
    }

    inline const bool operator !=(const Line& lhs, const Line& rhs)
    {
        return !(lhs == rhs);
    }

    inline const bool operator !=(const Circle& lhs, const Circle& rhs)
    {
        return !(lhs == rhs);
    }

    const Point Circle::GetPoint(double ang) const
    {
        return Point(C.x + r * cos(ang), C.y + r * sin(ang));
    }

    const Point Line::GetPoint(double t) const
    {
        return P + v * t;
    }
    
    const int dcmp(double d)
    {
        if (fabs(d) < eps) return 0;
        else return d < 0 ? -1 : 1;
    }

    inline const double Dot(const Vector& lhs, const Vector& rhs)
    {
        return lhs.x * rhs.x + lhs.y * rhs.y;
    }

    inline const double Length(const Vector& v)
    {
        return sqrt(Dot(v, v));
    }

    inline const double Angle(const Vector& lhs, const Vector& rhs)
    {
        return acos(Dot(lhs, rhs) / Length(lhs) / Length(rhs));
    }

    inline const double Angle(const Vector& v)
    {
        return atan2(v.y, v.x);
    }

    inline const double Cross(const Vector& lhs, const Vector& rhs)
    {
        return lhs.x * rhs.y - lhs.y * rhs.x;
    }

    inline const double Area2(const Point& A, const Point& B, const Point& C)
    {
        return Cross(B - A, C - A);
    }

    inline const Vector Rotate(const Vector& v, double rad)
    {
        return Vector(v.x * cos(rad) - v.y * sin(rad), v.x * sin(rad) + v.y * (cos(rad)));
    }

    inline const Vector Normal(const Vector& v)
    {
        double L = Length(v);
        return Vector(-v.y / L, v.x / L);
    }

    const int Intersection(const Line& lhs, const Line& rhs, std::vector<Point>& sol)
    {
        if (dcmp(Cross(lhs.v, rhs.v)) == 0)
        {
            if (On(lhs.P, rhs))
                return -1;
            return 0;
        }
        Vector u = lhs.P - rhs.P;
        double t = Cross(rhs.v, u) / Cross(lhs.v, rhs.v);
        sol.push_back(lhs.GetPoint(t));
        return 1;
    }

    const int Intersection(const Line& l, const Point& A, const Point& B, std::vector<Point>& sol)
    {
        std::vector<Point> temp;
        int res = Intersection(l, Line(A, B - A), sol);
        if (res <= 0) return res;
        if (ProperIntersection(l, A, B))
        {
            sol.push_back(temp[0]);
            return 1;
        }
        else return 0;
    }

    inline const int Intersection(const Point& A, const Point& B, const Line& l, std::vector<Point>& sol)
    {
        return Intersection(l, A, B, sol);
    }

    const int Intersection(const Point& A1, const Point& A2, const Point& B1, const Point& B2, std::vector<Point>& sol)
    {
        std::vector<Point> temp;
        int res = Intersection(Line(A1, A2 - A1), Line(B1, B2 - B1), sol);
        if (res <= 0) return res;
        if (ProperIntersection(A1, A2, B1, B2))
        {
            sol.push_back(temp[0]);
            return 1;
        }
        else return 0;
    }

    const int Intersection(const Circle& C, const Point& A, const Point& B, std::vector<Point>& sol)
    {
        std::vector<Point> temp;
        int res = Intersection(C, Line(A, B - A), sol);
        if (res <= 0) return res;
        int cnt = 0;
        for (std::size_t i = 0; i < temp.size(); i++)
            if (On(temp[i], A, B))
            {
                sol.push_back(temp[i]);
                cnt++;
            }
        return cnt;
    }

    inline const int Intersection(const Point& A, const Point& B, const Circle& C, std::vector<Point>& sol)
    {
        return Intersection(C, A, B, sol);
    }

    const int Intersection(const Line& lhs, const Circle& rhs, std::vector<Point>& sol)
    {
        double a = lhs.v.x, b = lhs.P.x - rhs.C.x, c = lhs.v.y, d = lhs.P.y - rhs.C.y;
        double e = a * a + c * c, f = 2 * (a * b + c * d), g = b * b + d * d - rhs.r * rhs.r;
        double delta = f * f - 4 * e * g;
        if (dcmp(delta) < 0) return 0;
        else if (dcmp(delta) == 0)
        {
            sol.push_back(lhs.GetPoint(-f / (2 * e)));
            return 1;
        }
        else
        {
            sol.push_back(lhs.GetPoint((-f - sqrt(delta)) / (2 * e)));
            sol.push_back(lhs.GetPoint((-f + sqrt(delta)) / (2 * e)));
            return 2;
        }
    }

    inline const int Intersection(const Circle& lhs, const Line& rhs, std::vector<Point>& sol)
    {
        return Intersection(rhs, lhs, sol);
    }

    const int Intersection(const Circle& lhs, const Circle& rhs, std::vector<Point>& sol)
    {
        double d = Distance(lhs.C, rhs.C);
        if (dcmp(d) == 0)
        {
            if (dcmp(lhs.r - rhs.r) == 0) return -1;
            else return 0;
        }
        if (dcmp(lhs.r + rhs.r - d) < 0) return 0;
        if (dcmp(fabs(lhs.r - rhs.r) - d) > 0) return 0;
        double a = Angle(rhs.C - lhs.C);
        double da = acos((lhs.r * lhs.r + d * d - rhs.r * rhs.r) / (2 * lhs.r * d));
        Point P1 = lhs.GetPoint(a - da), P2 = lhs.GetPoint(a + da);
        sol.push_back(P1);
        if (P1 == P2) return 1;
        sol.push_back(P2);
        return 2;
    }

    const int Intersection(const Line& lhs, const Polygon& rhs, std::vector<Point>& sol)
    {
        int cnt = 0;
        std::vector<Point> temp;
        for (std::size_t i = 0; i < rhs.size(); i++)
        {
            int res = Intersection(lhs, rhs[i], rhs[(i + 1) % rhs.size()], temp);
            if (res == -1) return -1;
            cnt += res;
        }
        for (std::size_t i = 0; i < temp.size(); i++)
            sol.push_back(temp[i]);
        return cnt;
    }

    inline const int Intersection(const Polygon& lhs, const Line& rhs, std::vector<Point>& sol)
    {
        return Intersection(rhs, lhs, sol);
    }

    const int Intersection(const Polygon& p, const Point& A, const Point& B, std::vector<Point>& sol)
    {
        std::vector<Point> temp;
        int res = Intersection(p, Line(A, B - A), sol);
        if (res <= 0) return res;
        int cnt = 0;
        for (std::size_t i = 0; i < temp.size(); i++)
            if (On(temp[i], A, B))
            {
                sol.push_back(temp[i]);
                cnt++;
            }
        return cnt;
    }

    inline const int Intersection(const Point& A, const Point& B, const Polygon& p, std::vector<Point>& sol)
    {
        return Intersection(p, A, B, sol);
    }

    const int Intersection(const Circle& lhs, const Polygon& rhs, std::vector<Point>& sol)
    {
        int cnt = 0;
        std::vector<Point> temp;
        for (std::size_t i = 0; i < rhs.size(); i++)
        {
            int res = Intersection(lhs, rhs[i], rhs[(i + 1) % rhs.size()], temp);
            if (res == -1) return -1;
            cnt += res;
        }
        for (std::size_t i = 0; i < temp.size(); i++)
            sol.push_back(temp[i]);
        return cnt;
    }

    inline const int Intersection(const Polygon& lhs, const Circle& rhs, std::vector<Point>& sol)
    {
        return Intersection(rhs, lhs, sol);
    }

    const int Intersection(const Polygon& lhs, const Polygon& rhs, std::vector<Point>& sol)
    {
        int cnt = 0;
        std::vector<Point> temp;
        for (std::size_t i = 0; i < rhs.size(); i++)
        {
            int res = Intersection(lhs, rhs[i], rhs[(i + 1) % rhs.size()], temp);
            if (res == -1) return -1;
            cnt += res;
        }
        for (std::size_t i = 0; i < temp.size(); i++)
            sol.push_back(temp[i]);
        return cnt;
    }

    const int Intersection(std::vector<Line> L, Polygon& p)
    {
        std::sort(L.begin(), L.end());
        int first, last;
        std::size_t n = L.size();
        Point* P = new Point[n];
        Line* l = new Line[n];
        l[first = last = 0] = L[0];
        for (std::size_t i = 1; i < n; i++)
        {
            while (first < last && !In(L[i], P[last - 1])) last--;
            while (first < last && !In(L[i], P[first])) first++;
            l[++last] = L[i];
            if (dcmp(Cross(l[last].v, l[last - 1].v)) == 0)
            {
                last--;
                if (In(l[last], L[i].P)) l[last] = L[i];
            }
            if (first < last)
            {
                std::vector<Point> temp;
                Intersection(l[last - 1], l[last], temp);
                P[last - 1] = temp[0];
            }
        }
        while (first < last && !In(l[first], P[last - 1])) last--;
        if (last - first <= 1) return 0;
        std::vector<Point> temp;
        Intersection(l[last], l[first], temp);
        P[last] = temp[0];
        p.clear();
        p.resize(n);
        int m = 0;
        for (int i = first; i <= last; i++) p[m++] = P[i];
        p.resize(m);
        return m;
    }

    const bool ProperIntersection(const Line& lhs, const Line& rhs)
    {
        return Cross(lhs.v, rhs.v) != 0;
    }

    const bool ProperIntersection(const Point& P1, const Point& P2, const Line& l)
    {
        double c1 = Cross(l.v, P1 - l.P), c2 = Cross(l.v, P2 - l.P);
        return dcmp(c1) * dcmp(c2) < 0;
    }

    inline const bool ProperIntersection(const Line& l, const Point& P1, const Point& P2)
    {
        return ProperIntersection(P1, P2, l);
    }

    const bool ProperIntersection(const Point& P1, const Point& P2, const Point& Q1, const Point& Q2)
    {
        double c1 = Cross(P2 - P1, Q1 - P1), c2 = Cross(P2 - P1, Q2 - P1),
               c3 = Cross(Q2 - Q1, P1 - Q1), c4 = Cross(Q2 - Q1, P2 - Q1);
        return dcmp(c1) * dcmp(c2) < 0 && dcmp(c3) * dcmp(c4) < 0;
    }

    inline const double Distance(const Point& A, const Point& B)
    {
        return Length(A - B);
    }

    const double Distance(const Point& P, const Line& l)
    {
        return fabs(Cross(l.v, P - l.P)) / Length(l.v);
    }

    inline const double Distance(const Line& l, const Point& P)
    {
        return Distance(P, l);
    }

    const double Distance(const Point& P, const Point& A, const Point& B)
    {
        if (A == B) return Distance(P, A);
        Vector v1 = B - A, v2 = P - A, v3 = P - B;
        if (dcmp(Dot(v1, v2)) < 0) return Length(v2);
        else if (dcmp(Dot(v1, v3)) > 0) return Length(v3);
        else return fabs(Cross(v1, v2)) / Length(v1);
    }

    const Point Project(const Point& P, const Line& l)
    {
        return l.GetPoint(Dot(l.v, P - l.P) / Dot(l.v, l.v));
    }

    inline const Point Project(const Line& l, const Point& P)
    {
        return Project(P, l);
    }

    inline const bool On(const Point& P, const Line& l)
    {
        return dcmp(Cross(P - l.P, l.v)) == 0;
    }

    inline const bool On(const Line& l, const Point& P)
    {
        return On(P, l);
    }

    const bool On(const Point& P, const Point& A, const Point& B)
    {
        return dcmp(Cross(A - P, B - P)) == 0 && dcmp(Dot(A - P, B - P)) < 0;
    }

    inline const bool On(const Point& P, const Circle& C)
    {
        return dcmp(Distance(P, C.C) - C.r) == 0;
    }

    inline const bool On(const Circle& C, const Point& P)
    {
        return On(P, C);
    }

    const bool On(const Point& P, const Polygon& p)
    {
        for (std::size_t i = 0; i < p.size(); i++)
            if (On(P, p[i], p[(i + 1) % p.size()]))
                return true;
        return false;
    }

    inline const bool On(const Polygon& p, const Point& P)
    {
        return On(P, p);
    }

    inline const bool In(const Point& P, const Line& l)
    {
        return Cross(l.v, P - l.P) > 0;
    }

    inline const bool In(const Line& l, const Point& P)
    {
        return In(P, l);
    }

    inline const bool In(const Point& P, const Circle& C)
    {
        return dcmp(Distance(P, C.C) - C.r) < 0;
    }

    inline const bool In(const Circle& C, const Point& P)
    {
        return In(P, C);
    }

    const bool In(const Point& P, const Polygon& p)
    {
        int wn = 0;
        std::size_t n = p.size();
        for (size_t i = 0; i < n; i++)
        {
            if (On(P, p[i], p[(i + 1) % n])) return false;
            int k = dcmp(Cross(p[(i + 1) % n] - p[i], P - p[i]));
            int d1 = dcmp(p[i].y - P.y);
            int d2 = dcmp(p[(i + 1) % n].y - P.y);
            if (k > 0 && d1 <= 0 && d2 > 0) wn++;
            if (k < 0 && d2 <= 0 && d1 > 0) wn--;
        }
        return wn != 0;
    }

    inline const bool In(const Polygon& p, const Point& P)
    {
        return In(P, p);
    }

    const double Area(const Polygon& p)
    {
        double area = 0;
        for (std::size_t i = 1; i < p.size() - 1; i++)
            area += Cross(p[i] - p[0], p[i + 1] - p[0]);
        return area / 2;
    }

    inline const double Area(const Circle& C)
    {
        return PI * C.r * C.r;
    }

    const int Tangents(const Point& P, const Circle& C, std::vector<Line>& sol)
    {
        Vector u = C.C - P;
        double dist = Length(u);
        if (dcmp(dist - C.r) < 0) return 0;
        else if (dcmp(dist - C.r) == 0)
        {
            sol.push_back(Line(P, Rotate(u, PI / 2)));
            return 1;
        }
        else
        {
            double ang = asin(C.r / dist);
            sol.push_back(Line(P, Rotate(u, -ang)));
            sol.push_back(Line(P, Rotate(u, ang)));
            return 2;
        }
    }

    inline const int Tangents(const Circle& C, const Point& P, std::vector<Line>& sol)
    {
        return Tangents(P, C, sol);
    }

    const int Tangents(Circle C1, Circle C2, std::vector<Line>& sol)
    {
        int cnt = 0;
        std::vector<Point> A, B;
        if (C1.r < C2.r)
        {
            using std::swap;
            swap(C1, C2);
            swap(A, B);
        }
        double d2 = Distance(C1.C, C2.C);
        double rdiff = C1.r - C2.r;
        double rsum = C1.r + C2.r;
        if (dcmp(d2 - rdiff * rdiff) < 0) return 0;
        double base = Angle(C2.C - C1.C);
        if (dcmp(d2) == 0 && dcmp(C1.r - C2.r) == 0) return -1;
        if (dcmp(d2 - rdiff * rdiff) == 0)
        {
            A.push_back(C1.GetPoint(base));
            B.push_back(C2.GetPoint(base));
            cnt++;
        }
        else
        {
            double ang = acos(rdiff / sqrt(d2));
            A.push_back(C1.GetPoint(base + ang));
            B.push_back(C2.GetPoint(base + ang));
            cnt++;
            A.push_back(C1.GetPoint(base - ang));
            B.push_back(C2.GetPoint(base - ang));
            cnt++;
            if (dcmp(d2 - rsum * rsum) == 0)
            {
                A.push_back(C1.GetPoint(base));
                B.push_back(C2.GetPoint(base));
                cnt++;
            }
            else if (dcmp(d2 - rsum * rsum) > 0)
            {
                double ang2 = acos(rsum / sqrt(d2));
                A.push_back(C1.GetPoint(base + ang2));
                B.push_back(C2.GetPoint(PI + base + ang2));
                cnt++;
                A.push_back(C1.GetPoint(base - ang2));
                B.push_back(C2.GetPoint(PI + base - ang2));
                cnt++;
            }
        }
        for (std::size_t i = 0; i < A.size(); i++)
            sol.push_back(Line(A[i], B[i] - A[i]));
        return cnt;
    }

    const Circle CircumscribedCircle(const Point& A, const Point& B, const Point& C)
    {
        Vector b = B - A, c = C - A;
        double D = 2 * Cross(b, c);
        double cx = (c.y * (b.x * b.x + b.y * b.y) - b.y * (c.x * c.x + c.y * c.y)) / D + A.x;
        double cy = (b.x * (c.x * c.x + c.y * c.y) - c.x * (b.x * b.x + b.y * b.y)) / D + A.y;
        return Circle(Point(cx, cy), Distance(A, Point(cx, cy)));
    }

    const Circle InscribedCircle(const Point& A, const Point& B, const Point& C)
    {
        double a = Distance(B, C), b = Distance(A, C), c = Distance(A, B);
        Point P = Point() + ((A - Point()) * a + (B - Point()) * b + (C - Point()) * c) / (a + b + c);
        return Circle(P, Distance(P, Line(A, B - A)));
    }

    const int ConvexHull(std::vector<Point> P, Polygon& p)
    {
        std::sort(P.begin(), P.end());
        int m = 0;
        p.clear();
        p.resize(P.size());
        for (std::size_t i = 0; i < P.size(); i++)
        {
            while (m > 1 && dcmp(Cross(p[m - 1] - p[m - 2], P[i] - p[m - 2])) < 0) m--;
            p[m++] = P[i];
        }
        int k = m;
        for (int i = (int)P.size() - 2; i >= 0; i--)
        {
            while (m > k && dcmp(Cross(p[m - 1] - p[m - 2], P[i] - p[m - 2])) < 0) m--;
            p[m++] = P[i];
        }
        if (P.size() > 1) m--;
        p.resize(m);
        return m;
    }

    const Polygon Cut(const Polygon& p, const Line& l)
    {
        Polygon newp;
        std::size_t n = p.size();
        for (std::size_t i = 0; i < n; i++)
        {
            Point C = p[i], D = p[(i + 1) % n];
            if (dcmp(Cross(l.v, C - l.P)) >= 0) newp.push_back(C);
            if (dcmp(Cross(l.v, C - D)) != 0)
            {
                std::vector<Point> temp;
                Intersection(l, Line(C, D - C), temp);
                if (On(temp[0], C, D)) newp.push_back(temp[0]);
            }
        }
        return newp;
    }

    inline const Polygon Cut(const Line& l, const Polygon& p)
    {
        return Cut(p, l);
    }
}
