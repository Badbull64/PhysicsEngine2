#pragma once
#ifndef QTRSTRUCT_H
#define QTRSTRUCT_H
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GLFW/glfw3.h>
#include<iostream>
#include <cmath>
#include <vector>
#include "particleC.h"
#include "Vector2D.h"

struct Point {
    Vector2D pos;
};

struct QTRect {
    //potential issue btw
    Vector2D pos;
    Vector2D dimensions; 
    Vector2D vertices[4];

    QTRect(const Vector2D& p = { -1.0f, -1.0f }, const Vector2D& s = { 2.0f, 2.0f }) : pos(p), dimensions(s)
    {
    };


    bool operator==(const QTRect& other) const
    {
        return (pos == other.pos) && (dimensions == other.dimensions);
    };

    bool operator!=(const QTRect& other) const
    {
        return (pos != other.pos) || (dimensions != other.dimensions);
    };

    bool containsp(const Point& point) const {
        return (point.pos.x >= pos.x - dimensions.x && point.pos.x <= pos.x + dimensions.x &&
            point.pos.y >= pos.y - dimensions.y && point.pos.y <= pos.y + dimensions.y);
    };
    template<typename T>
    bool contains(int id,std::vector<T> objvec) const {
        QTRect rangei = objvec[id].getBoundingBox();
        return (rangei.pos.x >= pos.x) && (rangei.pos.x + rangei.dimensions.x < pos.x + dimensions.x) &&
            (rangei.pos.y >= pos.y) && (rangei.pos.y + rangei.dimensions.y < pos.y + dimensions.y);
    };
    bool containr(const QTRect& rangei) const {
        return (rangei.pos.x >= pos.x) && (rangei.pos.x + rangei.dimensions.x < pos.x + dimensions.x) &&
            (rangei.pos.y >= pos.y) && (rangei.pos.y + rangei.dimensions.y < pos.y + dimensions.y);
    };

    bool intersects(const QTRect& range) const {
        return !(range.pos.x - range.dimensions.x > pos.x + dimensions.x ||
            range.pos.x + range.dimensions.x < pos.x - dimensions.x ||
            range.pos.y - range.dimensions.y > pos.y + dimensions.y ||
            range.pos.y + range.dimensions.y < pos.y - dimensions.y);
    };
    bool overlaps(const QTRect& range) const {
        return (pos.x < range.pos.x + range.dimensions.x && pos.x + dimensions.x >= range.pos.x && pos.y < range.pos.y + range.dimensions.y && pos.y + dimensions.y >= range.pos.y);
    }

};

struct QTCircle {
    Vector2D pos;
    float radius;
    float rsqrd = radius * radius;

    bool contains(const Point& point) const {
        float d = (point.pos.x - pos.x) * (point.pos.x - pos.x) + (point.pos.y - pos.y) * (point.pos.y - pos.y);
        return (d <= rsqrd);
    };

    bool intersects(const QTRect& range) const {
        double xDist = std::abs(range.pos.x - pos.x);
        double yDist = std::abs(range.pos.y - pos.y);

        double w = range.dimensions.x;
        double h = range.dimensions.y;

        double edges = (xDist - w) * (xDist - w) + (yDist - h) * (yDist - h);

        // no intersection
        if (xDist > (radius + w) || yDist > (radius + h))
            return false;

        // intersection within the circle
        if (xDist <= w || yDist <= h)
            return true;

        // intersection on the edge of the circle
        return (edges <= radius * radius);
    };
};


#endif //QTRSTRUCT_H