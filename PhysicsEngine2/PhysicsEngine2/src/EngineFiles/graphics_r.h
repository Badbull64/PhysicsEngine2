#pragma once
#ifndef GRAPHICS_R_H
#define GRAPHICS_R_H
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>
const float PIs = 3.1415926;
void bCircle_gr(float xPoss, float yPoss, float rad, float x, float y, float z) {
    const int steps = 42;
    const float angle = 2.0f * PIs / steps;
    glBegin(GL_TRIANGLE_FAN);
    glColor3f(x, y, z);
    glVertex2f(xPoss, yPoss);
    for (int i = 0; i <= steps; i++) {
        float nx = xPoss + rad * sin(angle * i);
        float ny = yPoss + rad * cos(angle * i);
        glVertex2f(nx, ny);
    }
    glEnd();
}
void DrawLines_gr(float xPos, float yPos, float nxPos, float nyPos, float width, float x, float y, float z) {
    //float x1 = (xPos / SCREEN_WIDTH) * 2.0 - 1.0;
    //float y1 = 1.0 - (yPos / SCREEN_HEIGHT) * 2.0;
    //float x2 = (nxPos / SCREEN_WIDTH) * 2.0 - 1.0;
    //float y2 = 1.0 - (nyPos / SCREEN_HEIGHT) * 2.0;
    glEnable(GL_LINE_SMOOTH);
    glColor3f(x, y, z);
    glLineWidth(width);
    glBegin(GL_LINE_LOOP);
    //glVertex2f(x1, y1);
    //glVertex2f(x2, y2);        
    glVertex2f(xPos, yPos);
    glVertex2f(nxPos, nyPos);
    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glFlush();
}
void DrawRectangle_gr(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4, float r, float g, float b, float alpha) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_QUADS);
    glColor4f(r, g, b, alpha);
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
    glVertex2f(x3, y3);
    glVertex2f(x4, y4);
    glEnd();
    glDisable(GL_BLEND);
    glFlush();
}

#endif // GRAPHICS_R_H
