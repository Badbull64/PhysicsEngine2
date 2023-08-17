#pragma once
#ifndef TIMESTEP_H
#define TIMESTEP_H
#include <GLFW/glfw3.h>

class TimeStep {
public:
    TimeStep() : deltaTime(0.0f), lastFrame(0.0f) {}

    // Call this function to start the time stepping
    void start() {
        lastFrame = static_cast<float>(glfwGetTime());
    }

    // Call this function at the end of each frame to update the time step
    void update() {
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
    }

    // Get the time elapsed between the current and last frame
    float dt() const {
        return deltaTime;
    }

private:
    float deltaTime; // Time between current frame and last frame
    float lastFrame; // Time of the last frame
};

#endif // TIMESTEP_H
