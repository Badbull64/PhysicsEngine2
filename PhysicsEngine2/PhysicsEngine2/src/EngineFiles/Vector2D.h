#pragma once
#ifndef VECTOR2D_H
#define VECTOR2D_H

struct Vector2D {
    float x;
    float y;

    Vector2D();
    Vector2D(float xVal, float yVal);
    Vector2D& operator+=(const Vector2D& other);
    Vector2D operator+(const Vector2D& other) const;
    Vector2D& operator-=(const Vector2D& other);
    Vector2D operator-(const Vector2D& other) const;
    Vector2D& operator*=(float scalar);
    Vector2D operator*(float scalar) const;
    Vector2D& operator/=(float scalar);
    Vector2D operator/(float scalar) const;
    bool operator!=(const Vector2D& other) const;
    bool operator ==(const Vector2D& other) const;
    float VecHorizontalAngle();
    void normalize();
    friend Vector2D operator*(float scalar, const Vector2D& vector);
    float operator*(const Vector2D& other) const;
    float operator%(const Vector2D& other) const;
};

#endif