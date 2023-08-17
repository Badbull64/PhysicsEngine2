#include "Vector2D.h"
#include <cmath>

Vector2D::Vector2D() : x(0.0f), y(0.0f) {}

Vector2D::Vector2D(float xVal, float yVal) : x(xVal), y(yVal) {}

Vector2D& Vector2D::operator+=(const Vector2D& other) {
    x += other.x;
    y += other.y;
    return *this;
}

Vector2D Vector2D::operator+(const Vector2D& other) const {
    Vector2D result(*this);
    result += other;
    return result;
}

bool Vector2D::operator!=(const Vector2D& other) const
{
    return (x != other.x) || (y != other.y);
}

bool Vector2D::operator==(const Vector2D& other) const
{
    return (x == other.x) || (y == other.y);
}

Vector2D& Vector2D::operator-=(const Vector2D& other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

Vector2D Vector2D::operator-(const Vector2D& other) const {
    Vector2D result(*this);
    result -= other;
    return result;
}

Vector2D& Vector2D::operator*=(float scalar) {
    x *= scalar;
    y *= scalar;
    return *this;
}

Vector2D Vector2D::operator*(float scalar) const {
    Vector2D result(*this);
    result *= scalar;
    return result;
}

Vector2D& Vector2D::operator/=(float scalar) {
    x /= scalar;
    y /= scalar;
    return *this;
}

Vector2D Vector2D::operator/(float scalar) const {
    Vector2D result(*this);
    result /= scalar;
    return result;
}

float Vector2D::VecHorizontalAngle() {
    return atan2(y, x);
}
void Vector2D::normalize() {
    float n = std::sqrtf(x * x + y * y);
    x /= n;
    y /= n;
}

Vector2D operator*(float scalar, const Vector2D& vector) {
    Vector2D result(vector);
    result *= scalar;
    return result;
} 

float Vector2D::operator*(const Vector2D& other) const {
    return x * other.y - y * other.x;
}

float Vector2D::operator%(const Vector2D& other) const {
    return x * other.x + y * other.y;
}