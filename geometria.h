#ifndef GEOMETRIA_H
#define GEOMETRIA_H

#include <iostream>
#include <cmath>
#include <utility>
#include <vector>

using std::min;
using std::max;

class Point {
private:
    double x, y;
public:
    Point() : x(0), y(0) {}
    Point(double _x, double _y) : x(_x), y(_y) {}

    double get_x() const{
        return this->x;
    }

    double get_y() const{
        return this->y;
    }

    void show() const{ 
        std::cerr << "(" << x << ", " << y << ")";
    }
};

using Polygon = std::vector<Point>;

double distance(const Point &a, const Point &b) {
    double dx = a.get_x() - b.get_x();
    double dy = a.get_y() - b.get_y();
    return std::sqrt(dx*dx + dy*dy);
}

int orientation(const Point &a, const Point &b, const Point &c) {
    //vetorial_product(AB, BC) = ABx.BCy - ABy.BCx
    std::pair<double, double> ab = {b.get_x() - a.get_x(), b.get_y() - a.get_y()};
    std::pair<double, double> bc = {c.get_x() - b.get_x(), c.get_y() - b.get_y()};
    double v = (ab.first * bc.second) - (ab.second * bc.first);

    if (std::abs(v) < 1e-9) return 0; // collinear
    return (v < 0) ? 1 : 2;           // 1: cw, 2: ccw
}

bool on_segment(const Point &a, const Point &b, const Point &c) { //a, b, c collinear points
    return (min(a.get_x(), c.get_x()) <= b.get_x() && b.get_x() <= max(a.get_x(), c.get_x())) 
            && (min(a.get_y(), c.get_y()) <= b.get_y() && b.get_y() <= max(a.get_y(), c.get_y()));
}

bool is_segments_intersection(const Point &p1, const Point &q1, const Point &p2, const Point &q2) {
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    if (o1 != o2 && o3 != o4) return true;
    //colinnear cases
    if (o1 == 0 && on_segment(p1, p2, q1)) return true;
    if (o2 == 0 && on_segment(p1, q2, q1)) return true;
    if (o3 == 0 && on_segment(p2, p1, q2)) return true;
    if (o4 == 0 && on_segment(p2, q1, q2)) return true;
    return false;
}

bool is_intersection_segment_polygon(const Point &p1, const Point &p2, const Polygon &polygon) { // Retorna true se o segmento (p1–p2) intercepta qualquer aresta do polígono
    int n = polygon.size();
    
    for (int i = 0; i < n; ++i) {
        const Point &a = polygon[i];
        const Point &b = polygon[(i+1)%n];
        if(is_segments_intersection(p1, p2, a, b))
            return true;
    }
    return false;
}

/*
using Poligono = std::vector<Point>;

// Calcula a distância euclidiana entre dois Points
double distancia(Point p1, Point p2) {
    return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
}
*/
#endif // GEOMETRIA_H