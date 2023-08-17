#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include "Vector2D.h";
using namespace std;
//use Verlet integration 
struct Object {
    Vector2D position;
    Vector2D position_last;
    Vector2D acceleration;
    float        radius = 10.0f;
    Object() = default;
    Object(Vector2D position_, float radius_)
        : position{ position_ }
        , position_last{ position_ }
        , acceleration{ 0.0f, 0.0f }
        , radius{ radius_ }
    {}

    void update(float dt)
    {
        // Compute how much we moved
        const Vector2D displacement = position - position_last;
        // Update position
        position_last = position;
        position = position + displacement + acceleration * (dt * dt);
        // Reset acceleration
        acceleration = {};
    }

    void accelerate(Vector2D a)
    {
        acceleration += a;
    }

    void setVelocity(Vector2D v, float dt)
    {
        position_last = position - (v * dt);
    }

    void addVelocity(Vector2D v, float dt)
    {
        position_last -= v * dt;
    }

    [[nodiscard]]
    Vector2D getVelocity(float dt) const
    {
        return (position - position_last) / dt;
    }
};

class Solver {
    
};