//usage : for points use pt<int> or pt<double> or pt<long long> or whatever type
//        for lines use line<int> or line<double> or line<long long> or whatever
//        not properly optimized
//        take from cpgeo book

namespace geometry{


    #ifndef M_PI
        #define M_PI acosl(-1)
    #endif
 
    template<typename T>
    struct pt {
        T x,y;
        pt operator+(pt p) {return {x+p.x, y+p.y};}
        pt operator-(pt p) {return {x-p.x, y-p.y};}
        pt operator*(T d) {return {x*d, y*d};}
        pt operator/(T d) {return {x/d, y/d};} // only for floating-point
    };
 

    template<typename T>
    bool operator==(pt<T> a, pt<T> b) {return a.x == b.x && a.y == b.y;}
    template<typename T>
    bool operator!=(pt<T> a, pt<T> b) {return !(a == b);}
 
    template<typename T>
    T sq(pt<T> p) {return p.x*p.x + p.y*p.y;}

    template<typename T>
    double abs(pt<T> p) {return sqrt(sq(p));}
 
 
    template<typename T>
    ostream& operator<<(ostream& os, pt<T> p) {
        #ifdef DBG_LOCAL
            return os << "("<< p.x << "," << p.y << ")";
        #else
            return os<<p.x<<" "<<p.y;
        #endif
    }
 
    template<typename T>
    istream& operator>>(istream& os, pt<T>& p) {
        return os>>p.x>>p.y;
    }
 
    template <typename T> 
    int sgn(T x) {
        return (T(0) < x) - (x < T(0));
    }
 
    template <typename T> 
    pt<T> translate(pt<T> v, pt<T> p) {return p+v;}
 
    template <typename T> 
    pt<T> scale(pt<T> c, double factor, pt<T> p) {
        return c + (p-c)*factor;
    }
 
    template <typename T> 
    pt<T> rot(pt<T> p, double a) {
        return {p.x*cos(a) - p.y*sin(a), p.x*sin(a) + p.y*cos(a)};
    }
 
    template <typename T> 
    pt<T> perp(pt<T> p) {return {-p.y, p.x};}
 
    template <typename T> 
    T dot(pt<T> v, pt<T> w) {return v.x*w.x + v.y*w.y;}
 
    template <typename T> 
    T cross(pt<T> v, pt<T> w) {return v.x*w.y - v.y*w.x;}
 
    template <typename T> 
    bool isPerp(pt<T> v, pt<T> w) {return dot(v,w) == 0;}
 
    template <typename T> 
    pt<T> linearTransfo(pt<T> p, pt<T> q, pt<T> r, pt<T> fp, pt<T> fq) {
        pt<T> pq = q-p, num{cross(pq, fq-fp), dot(pq, fq-fp)};
        return fp + pt<T>{cross(r-p, num), dot(r-p, num)} / sq(pq);
    }
 
    template <typename T> 
    double angle(pt<T> v, pt<T> w) {
        return acos(clamp(dot(v,w) / abs(v) / abs(w), (double)-1.0, (double)1.0));
    }
 
    template <typename T> 
    T orient(pt<T> a, pt<T> b, pt<T> c) {return cross(b-a,c-a);}
 
    template <typename T> 
    bool inAngle(pt<T> a, pt<T> b, pt<T> c, pt<T> p) {
        assert(orient(a,b,c) != 0);
        if (orient(a,b,c) < 0) swap(b,c);
        return orient(a,b,p) >= 0 && orient(a,c,p) <= 0;
    }
 
 
    template <typename T> 
    double orientedAngle(pt<T> a, pt<T> b, pt<T> c) {
        if (orient(a,b,c) >= 0)
        return angle(b-a, c-a);
        else
        return 2*M_PI - angle(b-a, c-a);
    }
 
 
    template <typename T> 
    bool isConvex(vector<pt<T>> p) {
        bool hasPos=false, hasNeg=false;
        for (int i=0, n=p.size(); i<n; i++) {
        int o = orient(p[i], p[(i+1)%n], p[(i+2)%n]);
        if (o > 0) hasPos = true;
        if (o < 0) hasNeg = true;
        }
        return !(hasPos && hasNeg);
    }
    template <typename T> 
    bool half(pt<T> p, pt<T> v = {0,1}) { // true if in blue half
        return cross(v,p) < 0 || (cross(v,p) == 0 && dot(v,p) <
0);
    }
 
    template <typename T> 
    void polarSort(vector<pt<T>> &v) {
        sort(v.begin(), v.end(), [](pt<T> s, pt<T> t) {
            return  make_tuple(half(s), 0, sq(s)) <
                   make_tuple(half(t), cross(s,t), sq(t));
        });
    }
 
    template <typename T> 
    void polarSortAround(pt<T> o, vector<pt<T>> &v) {
        sort(v.begin(), v.end(), [o](pt<T> s, pt<T> t) {
            return  make_tuple(half(s-o), 0, sq(s-o)) <
              make_tuple(half(t-o), cross(s-o, t-o), sq(t-o));
        });
    }
 
 
    template <typename T> 
    struct line {
        pt<T> v; T c;
        // From direction vector v and offset c
        line(pt<T> _v, T _c) : v(_v), c(_c) {}
        // From equation ax+by=c
        line(T _a, T _b, T _c) : v({_b,-_a}), c(_c) {}
        // From points P and Q
        line(pt<T> p, pt<T> q) : v(q-p), c(cross(v,p)) {}
        // Will be defined later:
        // - these work with T = int
        T side(pt<T> p) {return cross(v,p)-c;}
        double dist(pt<T> p) {return std::abs(side(p)) / abs(v);}
        double sqDist(pt<T> p) {return side(p)*side(p) / (double)sq(v);}
        line<T> perpthrough(pt<T> p) {return {p, p + perp(v)};}
        bool cmpProj(pt<T> p, pt<T> q) {
            return dot(v,p) < dot(v,q);
        }
        line<T> translate(pt<T> t) {return {v, c + cross(v,t)};}
        // - these require T = double
        line<T> shiftLeft(double dist) {return {v, c + dist*abs(v)};}
 
        pt<T> proj(pt<T> p) {return p - perp(v)*side(p)/sq(v);}
        pt<T> refl(pt<T> p) {return p - perp(v)*2*side(p)/sq(v);}
    };
 
    template <typename T> 
    ostream& operator<<(ostream& os, line<T> p) {
            return os<<-p.v.y<<" "<<p.v.x<<" "<<p.c<<endl;
    }
 
    template <typename T> 
    bool inter(line<T> l1, line<T> l2, pt<T> &out) {
        T d = cross(l1.v, l2.v);
        if (d == 0) return false;
        out = (l2.v*l1.c - l1.v*l2.c) / d; // requires floating-point coordinates
        return true;
    }
    template <typename T> 
    line<T> bisector(line<T> l1, line<T> l2, bool interior) {
        assert(cross(l1.v, l2.v) != 0); // l1 and l2 cannot be parallel!
        double sign = interior ? 1 : -1;
        return {l2.v/abs(l2.v) + l1.v/abs(l1.v) * sign,
                l2.c/abs(l2.v) + l1.c/abs(l1.v) * sign};
    }
    template <typename T> 
    bool inDisk(pt<T> a, pt<T> b, pt<T> p) {
        return dot(a-p, b-p) <= 0;
    }
    template <typename T> 
    bool onSegment(pt<T> a, pt<T> b, pt<T> p) {
        return orient(a,b,p) == 0 && inDisk(a,b,p);
    }
    template <typename T> 
    bool properInter(pt<T> a, pt<T> b, pt<T> c, pt<T> d, pt<T> &out) {
        double oa = orient(c,d,a),
        ob = orient(c,d,b),
        oc = orient(a,b,c),
        od = orient(a,b,d);
        // Proper intersection exists iff opposite signs
        if (oa*ob < 0 && oc*od < 0) {
            out = (a*ob - b*oa) / (ob-oa);
            return true;
        }
        return false;
    }
 
    template <typename T> 
    struct cmpX {
        bool operator()(const pt<T>& a, const pt<T>& b) const {
            return make_pair(a.x, a.y) < make_pair(b.x, b.y);
        }
    };
    template <typename T> 
    set<pt<T>,cmpX<T>> inters(pt<T> a, pt<T> b, pt<T> c, pt<T> d) {
        pt<T> out;
        if (properInter(a,b,c,d,out)) return {out};
        set<pt<T>,cmpX<T>> s;
        if (onSegment(c,d,a)) s.insert(a);
        if (onSegment(c,d,b)) s.insert(b);
        if (onSegment(a,b,c)) s.insert(c);
        if (onSegment(a,b,d)) s.insert(d);
        return s;
    }
 
    template <typename T> 
    double segPoint(pt<T> a, pt<T> b, pt<T> p) {
        if (a != b) {
            line l(a,b);
            if (l.cmpProj(a,p) && l.cmpProj(p,b)) // if closest to projection
            return l.dist(p);
            // output distance to line
        }
        return min(abs(p-a), abs(p-b)); // otherwise distance to A or B
    }
    template <typename T> 
    double segSeg(pt<T> a, pt<T> b, pt<T> c, pt<T> d) {
        pt<T> dummy;
        if (properInter(a,b,c,d,dummy))
            return 0;
        return min({segPoint(a,b,c), segPoint(a,b,d),
                    segPoint(c,d,a), segPoint(c,d,b)});
    }
 
 
    template <typename T> 
    double areaTriangle(pt<T> a, pt<T> b, pt<T> c) {
        return std::abs(cross(b-a, c-a)) / 2.0;
    }
    template <typename T> 
    double areaPolygon(const vector<pt<T>>& p) {
        double area = 0.0;
        for (int i = 0, n = p.size(); i < n; i++) {
            area += cross(p[i], p[(i+1)%n]); // wrap back to 0 if i == n-1
        }
        return std::abs(area) / 2.0;
    }
    // true if P at least as high as A (blue part)
    template <typename T> 
    bool above(pt<T> a, pt<T> p) {
        return p.y >= a.y;
    }
    // check if [PQ] crosses ray from A
    template <typename T> 
    bool crossesRay(pt<T> a, pt<T> p, pt<T> q) {
        return (above(a,q) - above(a,p)) * orient(a,p,q) > 0;
    }
    // if strict, returns false when A is on the boundary
    template <typename T> 
    bool inPolygon(vector<pt<T>> p, pt<T> a, bool strict = true) {
        int numCrossings = 0;
        for (int i = 0, n = p.size(); i < n; i++) {
            if (onSegment(p[i], p[(i+1)%n], a))
                return !strict;
            numCrossings += crossesRay(a, p[i], p[(i+1)%n]);
        }
        return numCrossings & 1; // inside if odd number of crossings
    }
    // amplitude travelled around point A, from P to Q
    template <typename T> 
    double angleTravelled(pt<T> a, pt<T> p, pt<T> q) {
        double ampli = angle(p-a, q-a);
        if (orient(a,p,q) > 0) return ampli;
        else return -ampli;
    }
 
    template <typename T> 
    int windingNumber(vector<pt<T>> p, pt<T> a) {
        double ampli = 0;
        for (int i = 0, n = p.size(); i < n; i++)
            ampli += angleTravelled(a, p[i], p[(i+1)%n]);
        return round(ampli / (2*M_PI));
    }
 
 
    //Circle
 
 
 
    template <typename T> 
    pt<T> circumCenter(pt<T> a, pt<T> b, pt<T> c) {
        b = b-a, c = c-a; // consider coordinates relative to A
        assert(cross(b,c) != 0); // no circumcircle if A,B,C aligned
        return a + perp(b*sq(c) - c*sq(b))/cross(b,c)/2;
    }
 
 
    template <typename T> 
    int circleLine(pt<T> o, double r, line<T> l, pair<pt<T>,pt<T>> &out) {
        double h2 = r*r - l.sqDist(o);
        if (h2 >= 0) { // the line touches the circle
            pt<T> p = l.proj(o); // point P
            pt<T> h = l.v*sqrt(h2)/abs(l.v); // vector parallel to l, of length h
            out = {p-h, p+h};
        }
        return 1 + sgn(h2);
    }
 
    template <typename T> 
    int circleCircle(pt<T> o1, double r1, pt<T> o2, double r2, pair<pt<T>,pt<T>> &out) {
        pt<T> d=o2-o1; double d2=sq(d);
        if (d2 == 0) {assert(r1 != r2); return 0;} // concentric circles
        double pd = (d2 + r1*r1 - r2*r2)/2; // = |O_1P| * d
        double h2 = r1*r1 - pd*pd/d2; // = hˆ2
        if (h2 >= 0) {
            pt<T> p = o1 + d*pd/d2, h = perp(d)*sqrt(h2/d2);
            out = {p-h, p+h};
        }
        return 1 + sgn(h2);
    }
 
}
 
using namespace geometry;
using pti = pt<int>;

