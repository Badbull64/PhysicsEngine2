#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <cmath>
#include <algorithm>
#include <math.h>
#include <random>
#include <atomic>
#include <map>
#include <stdlib.h>
#include <list>
#include <cstdint>
#include <array>
#include "EngineFiles/Id_Math.h"
#include "EngineFiles/Vector2D.h"
#include "EngineFiles/Timestep.h"
#include "EngineFiles/particleC.h"
#include "EngineFiles/qutree.h"
//#include <queue>
//elapsed time correct the spelling
//Do the waveFront project NEXT
//hey
//create a general object type;
using namespace std;
#define SCREEN_HEIGHT 850
#define SCREEN_WIDTH 850
#define SC_DIST 10000
GLFWwindow* window;
const float PI = 3.1415926;
bool check = true;
double felapsed;
double xPos, yPos;
float deltatime;
float elapsed = 0.0f;
bool m_mouseButtonPressed[2] = { false, false };
bool m_mouseButtonReleased[2] = { false, false };

void sleepMilliseconds(int milliseconds) {
    std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
}

struct wBall {
    float radius, E_store;
    float px, py;
    float vx, vy;
    float ax, ay;
    float ox, oy;
    float prev_px;
    float prev_py;
    float friction;
    int mass, Id;
    int Kx = 0.5f * mass * (vx) * (vx);
    int Ky = 0.5f * mass * (vy) * (vy);
    float fSimTimeRemaining;
};

struct wLineSegment {
    float sx, sy;
    //float sxV, syV;
    float ex, ey;
    //float exV, eyV;
    float radius;
};

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if (action == GLFW_PRESS)
        {
            m_mouseButtonPressed[0] = true;
        }
        else if (action == GLFW_RELEASE)
        {
            m_mouseButtonPressed[0] = false;
            m_mouseButtonReleased[0] = true;
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
        if (action == GLFW_PRESS)
        {
            m_mouseButtonPressed[1] = true;
        }
        else if (action == GLFW_RELEASE)
        {
            m_mouseButtonPressed[1] = false;
            m_mouseButtonReleased[1] = true;
        }
    }
}

void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    xPos = (xpos / width) * 2.0 - 1.0;
    yPos = 1.0 - (ypos / height) * 2.0;
}

void windowResizeCallback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

class BallPhysics {
public:
    std::vector<wBall> vectorBalls;
    std::vector<wLineSegment> vectorLines;
    wLineSegment* LSelected = nullptr;
    wBall* pSelected = nullptr;
    wBall* pTemp = nullptr;
    bool wSelctedLineStart = false;
    bool CSelectedLine = false;
    const float fStable = 0.005f;
    void bCircle(float xPos, float yPos, float rad, float x, float y, float z) {
        const int steps = 42;
        const float angle = 2.0f * PI / steps;
        glBegin(GL_TRIANGLE_FAN);
        glColor3f(x, y, z);
        glVertex2f(xPos, yPos);
        for (int i = 0; i <= steps; i++) {
            float nx = xPos + rad * sin(angle * i);
            float ny = yPos + rad * cos(angle * i);
            glVertex2f(nx, ny);
        }
        glEnd();
    }

    void DrawLines(float xPos, float yPos, float nxPos, float nyPos, float width) {
        glEnable(GL_LINE_SMOOTH);
        glColor3f(0, 0, 1.0f);
        glLineWidth(width);
        glBegin(GL_LINE_LOOP);
        glVertex2f(xPos, yPos);
        glVertex2f(nxPos, nyPos);
        glEnd();
        glDisable(GL_LINE_SMOOTH);
        glFlush();
    }

    void AddBalls(float x, float y, float r) {
        wBall b{};
        b.px = x; b.py = y;
        b.vx = 0; b.vy = 0;
        b.ax = 0; b.ay = 0;
        b.ox = 0; b.oy = 0;
        b.radius = r;
        //take mass to be in kg
        b.mass = 10.0f;
        b.friction = 0;
        b.fSimTimeRemaining = 0.0f;
        b.Id = vectorBalls.size();
        b.prev_px = x + 0.003;
        b.prev_py = y + 0.03;
        vectorBalls.emplace_back(b);
    }
    void RotateLines(float rotationSpeed)
    {
        static float angle = 0.0f;
        static float previousTime = 0.0f;

        // Calculate the elapsed time since the last rotation update

        float currentTime = glfwGetTime() * 0.001;
        float elapsedTime = currentTime - previousTime;

        // Calculate the change in angle based on the rotation speed and elapsed time
        float deltaAngle = rotationSpeed * elapsedTime;

        for (auto& line : vectorLines) {
            float centerX = 0.5f * (line.sx + line.ex);
            float centerY = 0.5f * (line.sy + line.ey);

            float dxStart = line.sx - centerX;
            float dyStart = line.sy - centerY;
            float dxEnd = line.ex - centerX;
            float dyEnd = line.ey - centerY;

            float radiusStart = std::sqrt(dxStart * dxStart + dyStart * dyStart);
            float radiusEnd = std::sqrt(dxEnd * dxEnd + dyEnd * dyEnd);

            float newAngleStart = std::atan2(dyStart, dxStart) + deltaAngle;
            float newAngleEnd = std::atan2(dyEnd, dxEnd) + deltaAngle;

            line.sx = centerX + radiusStart * std::cos(newAngleStart);
            line.sy = centerY + radiusStart * std::sin(newAngleStart);
            line.ex = centerX + radiusEnd * std::cos(newAngleEnd);
            line.ey = centerY + radiusEnd * std::sin(newAngleEnd);
        }

        previousTime = currentTime;
    }

    void BallCollisions() {
        int nSimulationUpdates = 1;
        float fSimElapsedTime = felapsed / (float)nSimulationUpdates;
        int nMaxSimulationSteps = 1;
        for (int i = 0; i < nSimulationUpdates; i++) {
            for (auto& ball : vectorBalls) {
                ball.fSimTimeRemaining = fSimElapsedTime;
            }
            for (int i = 0; i < nMaxSimulationSteps; i++)
            {
                for (auto& ball : vectorBalls)
                {
                    // Add Drag to emulate rolling friction
                    if (ball.fSimTimeRemaining > 0.0f) {
                        /* float dt = ball.fSimTimeRemaining;
                        float vx = (ball.px - ball.prev_px)*0.999;
                        float vy = (ball.py - ball.prev_py) * 0.999;

                        ball.ax = 0;
                        ball.ay = -0.0981;

                        ball.prev_px = ball.px;
                        ball.prev_py = ball.py;
                        ball.px += vx + ball.ax * dt * dt;
                        ball.py += vy + ball.ay * dt * dt;
                        // Wrap the balls around the screen
                        if (ball.px < -1.0f + ball.radius) {
                            ball.px = -1.0f + ball.radius;
                            ball.prev_px = (ball.px + (vx)* 0.899);
                        }
                        if (ball.px > 1.0f- ball.radius) {
                            ball.px = 1.0f- ball.radius;
                            ball.prev_px = (ball.px + (vx) * 0.899);
                        }
                        if (ball.py < -1.0f + ball.radius) {
                            ball.py = -1.0f + ball.radius;
                            ball.prev_py = (ball.py + (vy) * 0.899);
                        }
                        if (ball.py > 1.0f- ball.radius) {
                            ball.py = 1.0f - ball.radius;
                            ball.prev_py = (ball.py + (vy) * 0.899);
                        }*/
                        ball.ox = ball.px;								// Store original position this epoch
                        ball.oy = ball.py;

                        ball.ax = -ball.vx * 0.8f;						// Apply drag and gravity
                        ball.ay = -ball.vy * 0.8f + 100.0f;

                        ball.vx += ball.ax * ball.fSimTimeRemaining;	// Update Velocity
                        ball.vy += ball.ay * ball.fSimTimeRemaining;

                        ball.px += ball.vx * ball.fSimTimeRemaining;	// Update position
                        ball.py += ball.vy * ball.fSimTimeRemaining;

                        // Crudely wrap balls to screen - note this cause issues when collisions occur on screen boundaries
                        if (ball.px < -1 + ball.radius) ball.vx *= -1;
                        if (ball.px >= 1 - ball.radius) ball.vx *= -1;
                        if (ball.py < -1 + ball.radius) ball.vy *= -1;
                        if (ball.py >= 1 - ball.radius) ball.vy *= -1;

                        // Stop ball when velocity is neglible
                        if (fabs(ball.vx * ball.vx + ball.vy * ball.vy) < fStable)
                        {
                            ball.vx = 0;
                            ball.vy = 0;
                        }
                    }
                }

                std::vector<pair<wBall*, wBall*>> vecCollidingPairs;
                vector<wBall*> vecFakeBalls;
                // Static collisions, i.e. overlap
                for (auto& ball : vectorBalls)
                {
                    for (auto& edge : vectorLines) {
                        float fLineX1 = edge.ex - edge.sx;
                        float fLineY1 = edge.ey - edge.sy;
                        float fLineX2 = ball.px - edge.sx;
                        float fLineY2 = ball.py - edge.sy;
                        float fEdgeLength = fLineX1 * fLineX1 + fLineY1 * fLineY1;
                        //dot product
                        float dotProduct = fLineX1 * fLineX2 + fLineY1 * fLineY2;
                        float t = (dotProduct >= 0 && dotProduct <= fEdgeLength) ? dotProduct / fEdgeLength : 0;

                        float fClosestPointX = edge.sx + t * fLineX1;
                        float fClosestPointY = edge.sy + t * fLineY1;
                        float fDist = sqrtf(((ball.px - fClosestPointX) * (ball.px - fClosestPointX)) + ((ball.py - fClosestPointY) * (ball.py - fClosestPointY)));
                        if (fDist <= (ball.radius + edge.radius)) {
                            wBall* fakeBall = new wBall();
                            fakeBall->radius = edge.radius;
                            fakeBall->mass = ball.mass * 1.0f;
                            fakeBall->px = fClosestPointX;
                            fakeBall->py = fClosestPointY;
                            fakeBall->vx = -ball.vx;
                            fakeBall->vy = -ball.vy;
                            vecFakeBalls.push_back(fakeBall);
                            vecCollidingPairs.push_back({ &ball,fakeBall });
                            float fOverlap = 1.0f * (fDist - ball.radius - fakeBall->radius);
                            ball.px -= fOverlap * (ball.px - fakeBall->px) / fDist;
                            ball.py -= fOverlap * (ball.py - fakeBall->py) / fDist;

                        }

                    }
                    for (auto& target : vectorBalls)
                    {
                        if (ball.Id != target.Id)
                        {
                            if (((ball.px - target.px) * (ball.px - target.px)) + ((ball.py - target.py) * (ball.py - target.py)) <= (ball.radius + target.radius) * (ball.radius + target.radius))
                            {
                                // Collision has occured
                                vecCollidingPairs.push_back({ &ball, &target });

                                // Distance between ball centers
                                float fDistance = sqrtf((ball.px - target.px) * (ball.px - target.px) + (ball.py - target.py) * (ball.py - target.py));

                                // Calculate displacement required
                                float fOverlap = 0.5f * (fDistance - ball.radius - target.radius);

                                // Displace Current Ball away from collision
                                ball.px -= fOverlap * (ball.px - target.px) / fDistance;
                                ball.py -= fOverlap * (ball.py - target.py) / fDistance;
                                // Displace Target Ball away from collision
                                target.px += fOverlap * (ball.px - target.px) / fDistance;
                                target.py += fOverlap * (ball.py - target.py) / fDistance;
                            }
                        }
                    }
                    float fIntendedSpeed = sqrtf(ball.vx * ball.vx + ball.vy * ball.vy);
                    float fIntendedDist = fIntendedSpeed * ball.fSimTimeRemaining;
                    float fActualDist = sqrtf((ball.px - ball.ox) * (ball.px - ball.ox) + (ball.py - ball.oy) * (ball.py - ball.oy));
                    float fActualTime = fActualDist / fIntendedSpeed;

                    ball.fSimTimeRemaining = ball.fSimTimeRemaining - fActualTime;
                }

                // Now work out dynamic collisions
                for (auto c : vecCollidingPairs)
                {
                    wBall* b1 = c.first;
                    wBall* b2 = c.second;

                    // Distance between balls
                    float fDistance = sqrtf((b1->px - b2->px) * (b1->px - b2->px) + (b1->py - b2->py) * (b1->py - b2->py));

                    // Normal
                    float nx = (b2->px - b1->px) / fDistance;
                    float ny = (b2->py - b1->py) / fDistance;

                    // Tangent
                    float tx = -ny;
                    float ty = nx;

                    // Dot Product Tangent
                    float dpTan1 = b1->vx * tx + b1->vy * ty;
                    float dpTan2 = b2->vx * tx + b2->vy * ty;

                    // Dot Product Normal
                    float dpNorm1 = b1->vx * nx + b1->vy * ny;
                    float dpNorm2 = b2->vx * nx + b2->vy * ny;

                    // Conservation of momentum in 1D
                    float m1 = (dpNorm1 * (b1->mass - b2->mass) + 2.0f * b2->mass * dpNorm2) / (b1->mass + b2->mass);
                    float m2 = (dpNorm2 * (b2->mass - b1->mass) + 2.0f * b1->mass * dpNorm1) / (b1->mass + b2->mass);

                    // Update ball velocities
                    b1->vx = tx * dpTan1 + nx * m1;
                    b1->vy = ty * dpTan1 + ny * m1;
                    b2->vx = tx * dpTan2 + nx * m2;
                    b2->vy = ty * dpTan2 + ny * m2;
                }

                for (auto& b : vecFakeBalls) {
                    delete b;
                    vecFakeBalls.clear();

                    vecCollidingPairs.clear();
                }
            }
        }

    }

    void AddLine(float startx, float starty, float endx, float endy, float radius) {
        wLineSegment line;
        line.ex = startx;
        line.ey = starty;
        line.sx = endx;
        line.sy = endy;
        line.radius = radius;
        vectorLines.emplace_back(line);
    }

    void Conditions() {
        //can change
        float fRad = 0.02f;
        float fLineRad = 0.10f;
        //for (int i = 0; i < 1000; i++) {
        //    AddBalls(rand() % SCREEN_WIDTH-10, rand() % SCREEN_HEIGHT-10, fRad);
        //}
        AddBalls(0.2f, 0.2f, fRad);
        //AddBalls(0.3f,0.5f, fRad);
        //rand() % 40 + 2
        //vectorLines.push_back({-0.4f, -0.4f, 0.3f, 0.3f, fLineRad});
        //AddLine(-0.7f,0.3f,-0.7f,-0.3f,0.007f);
        //vectorLines.push_back({ 500.0f, 400.0f, 600.0f, 700.0f, (float)(fLineRad * 0.4) });

    }

    //0.5 * (line.sx + line.ex
    //loop func
    void loopPH() {
        //glfwSetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS, GLFW_TRUE);
        //bll.AddBalls(centerX + 200, centerY, fRad);
        for (auto& ball : vectorBalls) {
            bCircle(ball.px, ball.py, ball.radius, 0.02f, 0, 1.0f);
        }
        for (auto& line : vectorLines) {
            /* //line.mx = 0.5f * (line.ex + line.sx);
             //line.my = 0.5f * (line.ey + line.sy);
             bCircle(line.sx,line.sy,line.radius,0,0.1f,0.1f);
             bCircle(line.ex,line.ey,line.radius,0, 0.1f, 0.1f);
             float nx = -(line.ey - line.sy);
             float ny = (line.ex - line.sx);
             float d = sqrt(nx*nx + ny * ny);
             nx /= d;
             ny /= d;
             DrawLines((line.sx+nx*line.radius),(line.sy + ny * line.radius),(line.ex + nx * line.radius),(line.ey + ny * line.radius),8);
             DrawLines((line.sx-nx*line.radius),(line.sy - ny * line.radius),(line.ex - nx * line.radius),(line.ey - ny * line.radius),8); */
            bCircle(line.sx, line.sy, line.radius * 3, 0, 1.0f, 0);
            bCircle(line.ex, line.ey, line.radius * 3, 0, 1.0f, 0);
            DrawLines(line.ex, line.ey, line.sx, line.sy, 5.0f);
        }
        //draw the velcoity line here so if slelected != nullptr
        //end of velcority if statement
        bool currentLeftButtonState = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS);
        bool currentRightState = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS);
        //comment for the original || (currentRightState && !prevRightState)
        LSelected = nullptr;
        if (m_mouseButtonPressed[0] || m_mouseButtonPressed[1]) {
            pSelected = nullptr;
            for (auto& ball : vectorBalls) {
                if (((ball.px - xPos) * (ball.px - xPos)) + ((ball.py - yPos) * (ball.py - yPos)) < (ball.radius * ball.radius)) {
                    pSelected = &ball;
                    pTemp = pSelected;
                    break;
                }
            }
        }

        if (currentLeftButtonState) {
            if (pSelected != nullptr && LSelected == nullptr) {
                pSelected->px = xPos;
                pSelected->py = yPos;
                pSelected->prev_px = xPos;
                pSelected->prev_py = yPos;

            }
        }
        //check the nullptr state
        if (m_mouseButtonReleased[0]) {
            pSelected = nullptr;
            LSelected = nullptr;
            //Lcentre = nullptr;
            m_mouseButtonReleased[0] = false;
        }
        //DRAW QUE LINE
        if (m_mouseButtonReleased[1]) {
            if (pTemp != nullptr) {
                pTemp->vx = 5.0f * ((pTemp->px) - xPos);
                pTemp->vy = 5.0f * ((pTemp->py) - yPos);
            }
            pSelected = nullptr;
            LSelected == nullptr;
            m_mouseButtonReleased[1] = false;
            pTemp = nullptr;
        }
        if (m_mouseButtonPressed[0] || m_mouseButtonPressed[1]) {
            for (auto& line : vectorLines) {
                if (((line.sx - xPos) * (line.sx - xPos)) + ((line.sy - yPos) * (line.sy - yPos)) < (line.radius * line.radius)) {
                    LSelected = &line;
                    wSelctedLineStart = true;
                    break;
                }
                if (((line.ex - xPos) * (line.ex - xPos)) + ((line.ey - yPos) * (line.ey - yPos)) < (line.radius * line.radius)) {
                    LSelected = &line;
                    wSelctedLineStart = false;
                    break;
                }
            }
        }

        if (currentLeftButtonState) {
            if (LSelected != nullptr) {
                if (wSelctedLineStart) {
                    LSelected->sx = xPos;
                    LSelected->sy = yPos;
                }
                else {
                    LSelected->ex = xPos;
                    LSelected->ey = yPos;
                }
            }
        }
        if (currentRightState) {
            if (LSelected != nullptr) {
                if (wSelctedLineStart) {
                    int deltaX = xPos - LSelected->sx;
                    int deltaY = yPos - LSelected->sy;
                    LSelected->sx = xPos;
                    LSelected->sy = yPos;
                    LSelected->ex += deltaX;
                    LSelected->ey += deltaY;
                }
                else {
                    int deltaX = xPos - LSelected->ex;
                    int deltaY = yPos - LSelected->ey;
                    LSelected->ex = xPos;
                    LSelected->ey = yPos;
                    LSelected->sx += deltaX;
                    LSelected->sy += deltaY;
                }
            }
        }
        //distance between ball centre and object centre
        //velocity slider
        //if (pSelected != nullptr) cout << "current Ball kinetic Energy,  {Kx: " << pSelected->Kx <<" , Ky:" << pSelected->Kx << "}" << '\n';
        BallCollisions();
        RotateLines(1000);
    }



};
struct Ray {
    Vector2D PivotPos;
    Vector2D DragPos;
    Vector2D intersectionPos;
    Vector2D reflectionDir;
    float radius;
};

struct mirror {
    Vector2D edge1;
    Vector2D edge2;
    bool isplane;
    float radius = 20.0f;
    float id;
    float m_linewidth;
    //refractive stuff
    bool isrefractive;
    float refractiveIndex;
    bool isside1;
    bool isside2;
    bool isside3;
    bool isside4;
    bool ismirror;
};

struct plane {
    Vector2D edge1;
    Vector2D edge2;
    Vector2D edge3;
    Vector2D edge4;
    float linewidthw;
};

struct obj_r {
    Vector2D pos;
    float radius;
};

//for this function remeber to do memeory management

//later today start work on the 
class WaveFront {
private:

public:
    vector<mirror> mVec;
    vector<Ray> Rvec;
    vector<plane> pVec;
    Ray* dragSelected = nullptr;
    mirror* mSelected = nullptr;
    plane* planeselected = nullptr;
    vector<Vector2D> intersectionPoints;
    vector<float> refractiveStates;

    //outside snells law define as air
    float AirIdx = 1.0f;
    bool dintsct(const Vector2D& p3, const Vector2D& temp, const Vector2D& p1, const Vector2D& p2, Vector2D& intersect) {
        float x1 = p1.x;
        float y1 = p1.y;
        float x2 = p2.x;
        float y2 = p2.y;
        float x3 = p3.x;
        float y3 = p3.y;
        float x4 = p3.x + temp.x;
        float y4 = p3.y + temp.y;

        float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (den == 0) {
            return false;
        }

        float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
        float u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / den;
        if (t > 0.01 && t < 1 && u > 0.01) {
            intersect.x = (x1 + t * (x2 - x1));
            intersect.y = (y1 + t * (y2 - y1));
            return true;
        }
        else {
            return false;
        }
        //return false;
    }

    Vector2D NormalVectorTransform(const Vector2D& edge, const Vector2D& LightPos) {
        Vector2D normalVec;
        normalVec.x = -edge.y;
        normalVec.y = edge.x;
        float dP = LightPos.x * normalVec.x + LightPos.y * normalVec.y;
        if (dP > 0) {
            normalVec.x = -normalVec.x;
            normalVec.y = -normalVec.y;
        }

        return normalVec;
    }

    Vector2D reflectVector(const Vector2D& vector, const Vector2D& wallNormal) {
        double dotProduct = vector.x * wallNormal.x + vector.y * wallNormal.y;

        Vector2D reflectedVector;
        reflectedVector.x = vector.x - (2 * dotProduct * wallNormal.x) + 0.001f;
        reflectedVector.y = vector.y - (2 * dotProduct * wallNormal.y) + 0.001f;

        return reflectedVector;
    }

    Vector2D refractVector(const Vector2D& vectorI, const Vector2D& normalVec, float Rid) {
        //rid is the ration of n1 > n2 where rid = n1/n2 so rid is should allways be less than 1
        Vector2D k = -1 * vectorI;
        float dot = k.x * normalVec.x + k.y * normalVec.y;
        Vector2D refractedVector;
        float sqrtval = 1 - ((Rid) * (Rid)) * (1 - ((dot) * (dot)));
        float w;
        if (sqrtval < 0) {
            float dot2 = vectorI.x * normalVec.x + vectorI.y * normalVec.y;
            refractedVector.x = vectorI.x - (2 * dot2 * normalVec.x) + 0.0000099f;
            refractedVector.y = vectorI.y - (2 * dot2 * normalVec.y) + 0.0000099f;
        }
        else {
            w = (Rid * dot - std::sqrtf(sqrtval));
            refractedVector = w * normalVec - Rid * k;
        }


        return refractedVector;
    }


    void bCircle(float xPoss, float yPoss, float rad, float x, float y, float z) {
        const int steps = 42;
        const float angle = 2.0f * PI / steps;
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
    //change the xpos and ypos tings
    void DrawLines(float xPos, float yPos, float nxPos, float nyPos, float width, float x, float y, float z) {
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

    void DrawRectangle(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4, float r, float g, float b, float alpha) {
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

    int windingNumber(Vector2D p, Vector2D a, Vector2D b) {
        if (a.y <= p.y) {
            if (b.y > p.y && (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x) > 0)
                return 1;
        }
        else {
            if (b.y <= p.y && (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x) < 0)
                return -1;
        }
        return 0;
    }

    // function to check if a point P lies inside a quadrilateral 
    bool m_pointcheck(Vector2D A, Vector2D B, Vector2D C, Vector2D D, Vector2D P) {
        int wn = 0;
        wn += windingNumber(P, A, B);
        wn += windingNumber(P, B, C);
        wn += windingNumber(P, C, D);
        wn += windingNumber(P, D, A);
        return wn != 0;
    }

    void AddRay(float x, float y, float dragradius) {
        Ray r;
        r.DragPos.x = x;
        r.DragPos.y = y;
        r.PivotPos.x = x + 0.09f;
        r.PivotPos.y = y + 0.09f;
        r.intersectionPos.x = 0;
        r.intersectionPos.y = 0;
        r.reflectionDir.x = 0;
        r.reflectionDir.y = 0;
        r.radius = dragradius;
        Rvec.emplace_back(r);
    }

    void AddMirror(Vector2D edgestart, Vector2D edgeEnd, float mirrorId, float linews, bool isrefractives) {
        mirror m;
        m.edge1 = edgestart;
        m.edge2 = edgeEnd;
        m.id = mirrorId;
        m.isrefractive = isrefractives;
        m.m_linewidth = linews;
        m.isplane = false;
        m.ismirror = true;
        mVec.emplace_back(m);
    }

    void AddRefractionSurface(Vector2D edgestart, Vector2D edgeEnd, float refractiveID, bool isrefractives, float linew, bool isplanes, bool s1, bool s2, bool s3, bool s4) {
        mirror r;
        r.edge1 = edgestart;
        r.edge2 = edgeEnd;
        r.refractiveIndex = refractiveID;
        r.isrefractive = isrefractives;
        r.m_linewidth = linew;
        r.isplane = isplanes;
        r.isside1 = s1;
        r.isside2 = s2;
        r.isside3 = s3;
        r.isside4 = s4;
        r.ismirror = false;
        mVec.emplace_back(r);
    }

    void AddRefractionPlane(Vector2D edge1, Vector2D edge2, Vector2D edge3, Vector2D edge4, float lw, float refractiveindex) {
        plane p1;
        p1.edge1 = edge1;
        p1.edge2 = edge2;
        p1.edge3 = edge3;
        p1.edge4 = edge4;
        p1.linewidthw = lw;
        pVec.emplace_back(p1);
        AddRefractionSurface(p1.edge1, p1.edge2, refractiveindex, true, 3.0f, true, true, false, false, false);
        AddRefractionSurface(p1.edge2, p1.edge3, refractiveindex, true, 3.0f, true, false, true, false, false);
        AddRefractionSurface(p1.edge3, p1.edge4, refractiveindex, true, 3.0f, true, false, false, true, false);
        AddRefractionSurface(p1.edge1, p1.edge4, refractiveindex, true, 3.0f, true, false, true, false, true);
    }

    void Conditions() {

        AddRay(-0.3f, -0.3f, 0.035f);

        //limit it too 1 refraction plane
        AddRefractionPlane(Vector2D(0.4f, 0.4f), Vector2D(0.6f, 0.4f), Vector2D(0.6f, 0.55f), Vector2D(0.4f, 0.55f), 0.05f, 1.33f);
        //AddRefractionSurface(Vector2D(0.3f, 0.3f), Vector2D(0.7f, 0.1f),1.33f,true,0.05f,false,false,false,false,false);
        //AddRefractionSurface(Vector2D(0.2f, 0.2f), Vector2D(0.2f, 0.4f),1.33f,true,0.05f,false,false,false,false,false);
        AddMirror(Vector2D(-0.7f, 0.8f), Vector2D(-0.7f, -0.1f), 0, 0.02f, false);
        AddMirror(Vector2D(-0.3f, 0.8f), Vector2D(-0.3f, -0.1f), 0, 0.02f, false);
        AddMirror(Vector2D(-0.4f, 0.8f), Vector2D(-0.4f, -0.1f), 0, 0.02f, false);
        //AddMirror(Vector2D(-0.5f, 0.8f), Vector2D(-0.5f, -0.1f), 0, 0.02f, false);
        //AddMirror(Vector2D(-0.1f, 0.8f), Vector2D(-0.1f, -0.1f), 0, 0.02f, false);
    }

    void initInteReflec(Vector2D& direction, Vector2D& intersect, bool& istrue, bool& ismirrorcheck) {
        Vector2D intersectedo;
        float dist = std::numeric_limits<float>::max();
        float ClosestWall = std::numeric_limits<float>::max();
        for (size_t i = 0; i < Rvec.size(); i++) {
            Vector2D edge1m;
            Vector2D edge2m;
            Vector2D cInt;
            bool intersectionfoundm = false;
            bool isrefractivem = false;
            float ridxm = 0;
            Vector2D dir = Rvec[i].PivotPos - Rvec[i].DragPos;
            float nD = std::sqrtf(dir.x * dir.x + dir.y * dir.y);
            dir /= nD;
            bool isaplane = false;
            for (size_t j = 0; j < mVec.size(); j++) {
                mVec[j].id = 0;
                if (dintsct(Rvec[i].PivotPos, dir, mVec[j].edge1, mVec[j].edge2, intersectedo)) {
                    dist = (intersectedo.x - Rvec[i].PivotPos.x) * (intersectedo.x - Rvec[i].PivotPos.x) + (intersectedo.y - Rvec[i].PivotPos.y) * (intersectedo.y - Rvec[i].PivotPos.y);
                    if (dist < ClosestWall) {
                        ClosestWall = dist;
                        edge1m = mVec[j].edge1;
                        edge2m = mVec[j].edge2;
                        isaplane = mVec[j].isplane;
                        intersectionfoundm = true;
                        cInt = intersectedo;
                        isrefractivem = mVec[j].isrefractive;
                        ridxm = mVec[j].refractiveIndex;
                        mVec[j].id = 1;
                        ismirrorcheck = mVec[j].ismirror;
                    }
                }
            }
            if (intersectionfoundm) {
                Vector2D wallnormal = NormalVectorTransform(edge2m - edge1m, dir);
                float nW = std::sqrtf(wallnormal.x * wallnormal.x + wallnormal.y * wallnormal.y);
                istrue = true;
                wallnormal /= nW;
                if (isrefractivem) {
                    //implement the point polygon program
                    float nr = 0;
                    if (isaplane) {
                        for (auto& planes : pVec) {
                            if (m_pointcheck(planes.edge1, planes.edge2, planes.edge3, planes.edge4, Rvec[i].PivotPos)) {
                                nr = ridxm;
                            }
                            else {
                                nr = 1 / ridxm;
                            }
                        }
                    }
                    else {
                        nr = 1 / ridxm;
                    }
                    direction = refractVector(dir, wallnormal, nr);
                }
                else {
                    direction = reflectVector(dir, wallnormal);

                }
                float dti = std::sqrtf(direction.x * direction.x + direction.y * direction.y);
                direction /= dti;
                intersect = cInt;

            }

        }
    }

    void RayRecursion(int& reflectedraycount, Vector2D& fdir, int nRay) {
        if (nRay <= 0) {
            cout << "input Valid Ray number" << endl;
            return;
        }

        Vector2D rDir, intersection, nextIntersections, tempIntersection;
        bool ifIntersects = false;
        bool ismirrortrue = false;
        initInteReflec(rDir, intersection, ifIntersects, ismirrortrue);
        const Vector2D rintersection = intersection;
        bool ist = false;
        if (ifIntersects) {
            intersectionPoints.emplace_back(rintersection);
            reflectedraycount = 1;
            for (int i = 0; i < nRay; i++) {
                bool currentmirrorstate = false;
                float nq = 0;
                Vector2D edge1store;
                Vector2D edge2store;
                bool isrefractiveM_r = false;
                float ridx = 0;
                bool istrueIntersection = false;
                float ClosestMr = std::numeric_limits<float>::max();;
                bool isplanecheck = false;
                for (size_t j = 0; j < mVec.size(); j++) {
                    if (dintsct(intersection, rDir, mVec[j].edge2, mVec[j].edge1, nextIntersections)) {
                        float distrr = (intersection.x - nextIntersections.x) * (intersection.x - nextIntersections.x) + (intersection.y - nextIntersections.y) * (intersection.y - nextIntersections.y);
                        if (distrr < ClosestMr) {
                            ClosestMr = distrr;
                            isrefractiveM_r = mVec[j].isrefractive;
                            ridx = mVec[j].refractiveIndex;
                            isplanecheck = mVec[j].isplane;
                            currentmirrorstate = mVec[j].ismirror;
                            tempIntersection = nextIntersections;
                            edge1store = mVec[j].edge1;
                            edge2store = mVec[j].edge2;
                        }
                        istrueIntersection = true;
                    }
                }
                if (istrueIntersection) {
                    intersection = tempIntersection;
                    Vector2D q;
                    Vector2D wsallnormals = NormalVectorTransform(edge2store - edge1store, rDir);
                    float nsWss = std::sqrtf(wsallnormals.x * wsallnormals.x + wsallnormals.y * wsallnormals.y);
                    wsallnormals /= nsWss;
                    if (isrefractiveM_r) {
                        if (isplanecheck && !ismirrortrue) {
                            nq = ridx;
                            q = refractVector(rDir, wsallnormals, nq);
                        }
                        else {
                            nq = 1 / ridx;
                            q = refractVector(rDir, wsallnormals, nq);
                        }

                    }
                    else {
                        //nq = 1 / ridx;
                        q = reflectVector(rDir, wsallnormals);
                    }
                    float Rn = std::sqrtf(q.x * q.x + q.y * q.y);
                    q /= Rn;
                    rDir = q;
                    intersectionPoints.emplace_back(tempIntersection);

                }
                ismirrortrue = currentmirrorstate;
                if (i == nRay - 1 && istrueIntersection == false) {
                    fdir = rDir;
                }
            }
        }
    }

    void loopWF() {
        Vector2D fdir;
        bool edge1selcted = false;
        int reflectedRayCount = 0;
        bool edge1s = false;
        bool edge2s = false;
        bool edge3s = false;
        bool edge4s = false;
        //final arg is the number of rays you want
        RayRecursion(reflectedRayCount, fdir, 100);
        //draw rays
        for (auto& r : Rvec) {
            bCircle(r.DragPos.x, r.DragPos.y, r.radius, 0.1f, 0.5f, 1.0f);
            bCircle(r.PivotPos.x, r.PivotPos.y, 0.01f, 0.7f, 0.2f, 0.2f);
            //yellow colour 1.02f,1.02f,0
            DrawLines(r.DragPos.x, r.DragPos.y, r.PivotPos.x, r.PivotPos.y, 0.009f, 1.02f, 1.02f, 0);
            if (intersectionPoints.size() > 1) {
                DrawLines(r.PivotPos.x, r.PivotPos.y, intersectionPoints[0].x, intersectionPoints[0].y, 0.009f, 1.02f, 1.02f, 0);

            }
            else if (intersectionPoints.size() == 1) {
                //this case kinda ignored for the plane
                Vector2D dirse = r.PivotPos - r.DragPos;
                float norms = std::sqrtf(dirse.x * dirse.x + dirse.y * dirse.y);
                dirse /= norms;
                Vector2D wallnorm;
                Vector2D dirr;
                for (auto& m : mVec) {
                    if (m.id == 1) {
                        wallnorm = NormalVectorTransform(m.edge2 - m.edge1, dirse);

                        float wm = std::sqrtf(wallnorm.x * wallnorm.x + wallnorm.y * wallnorm.y);
                        wallnorm /= wm;
                        //DrawLines(intersectionPoints[0].x, intersectionPoints[0].y, intersectionPoints[0].x + wallnorm.x * 100.0f, intersectionPoints[0].y + wallnorm.y * 100.0f, 1.0f, 1.0f, 0, 0);
                        float dotpdt = wallnorm.x * dirse.x + wallnorm.y * dirse.y;
                        if (m.isrefractive) {
                            float dq = 1 / m.refractiveIndex;
                            dirr = refractVector(dirse, wallnorm, dq);
                        }
                        else {
                            dirr = reflectVector(dirse, wallnorm);
                        }
                    }
                }
                float dqp = std::sqrtf(dirr.x * dirr.x + dirr.y * dirr.y);
                dirr /= dqp;
                //yellow colour 1.02f, 1.02f, 0
                DrawLines(r.PivotPos.x, r.PivotPos.y, intersectionPoints[0].x, intersectionPoints[0].y, 0.009f, 1.02f, 1.02f, 0);
                DrawLines(intersectionPoints[0].x, intersectionPoints[0].y, intersectionPoints[0].x + dirr.x * SC_DIST, intersectionPoints[0].y + dirr.y * SC_DIST, 0.009f, 1.02f, 1.02f, 0);
            }
            else {
                Vector2D direction = r.PivotPos - r.DragPos;
                float norm = std::sqrtf(direction.x * direction.x + direction.y * direction.y);
                direction /= norm;;
                Vector2D endPos = r.PivotPos + direction * SC_DIST;
                //uncomment for the refractive stuff

                //test case line yellow colour 1.02f, 1.02f, 0
                DrawLines(r.PivotPos.x, r.PivotPos.y, endPos.x, endPos.y, 0.009f, 1.02f, 1.02f, 0);
            }

        }

        if (intersectionPoints.size() > 1) {
            for (size_t i = 0; i < intersectionPoints.size() - 1; i++) {
                DrawLines(intersectionPoints[i].x, intersectionPoints[i].y, intersectionPoints[i + 1].x, intersectionPoints[i + 1].y, 0.009f, 1.02f, 1.02f, 0);

            }
            if (fdir.x != 0.0f && fdir.y != 0.0f) {
                DrawLines(intersectionPoints[intersectionPoints.size() - 1].x, intersectionPoints[intersectionPoints.size() - 1].y, intersectionPoints[intersectionPoints.size() - 1].x + fdir.x * SC_DIST, intersectionPoints[intersectionPoints.size() - 1].y + fdir.y * SC_DIST, 0.009f, 1.02f, 1.02f, 0);
            }
        }

        for (auto& m : mVec)
        {
            if (m.isrefractive && m.isplane == false) {
                //idea for later change the colour corresponding to refractive index to the higher the index the darker the colour is 
                DrawLines(m.edge1.x, m.edge1.y, m.edge2.x, m.edge2.y, m.m_linewidth, 0.12f, 0.77f, 0.64f);
                //DrawRectangle(m.edge1.x,m.edge1.y,m.edge2.x,m.edge2.y,2.0f, 0.3685f, 0.67f, 0.5193f);

                //DrawRectangle(m.edge1.x, m.edge1.y, m.edge2.x, m.edge2.y, m.edge3.x, m.edge3.y, m.edge4.x, m.edge4.y, 0.3685f, 0.67f, 0.5193f);
            }
            else if (!m.isrefractive) {
                DrawLines(m.edge1.x, m.edge1.y, m.edge2.x, m.edge2.y, m.m_linewidth, 0.658f, 0.6916f, 0.7f);
            }
        }

        for (auto& pl : pVec) {
            bCircle(pl.edge1.x, pl.edge1.y, 0.01f, 0.1f, 0.2f, 0.7f);
            bCircle(pl.edge2.x, pl.edge2.y, 0.01f, 0.1f, 0.2f, 0.7f);
            bCircle(pl.edge3.x, pl.edge3.y, 0.01f, 0.1f, 0.2f, 0.7f);
            bCircle(pl.edge4.x, pl.edge4.y, 0.01f, 0.1f, 0.2f, 0.7f);
            DrawRectangle(pl.edge1.x, pl.edge1.y, pl.edge2.x, pl.edge2.y, pl.edge3.x, pl.edge3.y, pl.edge4.x, pl.edge4.y, 0.12f, 0.77f, 0.64f, 0.5f);
            for (auto& m2 : mVec) {
                if (m2.isplane && m2.isrefractive) {
                    if (m2.isside1) {
                        m2.edge1 = pl.edge1;
                        m2.edge2 = pl.edge2;
                    }
                    if (m2.isside2) {
                        m2.edge1 = pl.edge2;
                        m2.edge2 = pl.edge3;
                    }
                    if (m2.isside3) {
                        m2.edge1 = pl.edge3;
                        m2.edge2 = pl.edge4;
                    }
                    if (m2.isside4) {
                        m2.edge1 = pl.edge1;
                        m2.edge2 = pl.edge4;
                    }
                }
            }
        }

        if (m_mouseButtonPressed[0] || m_mouseButtonPressed[1]) {
            dragSelected = nullptr;
            mSelected = nullptr;
            planeselected = nullptr;
            for (auto& r : Rvec) {
                if (((r.DragPos.x - xPos) * (r.DragPos.x - xPos)) + ((r.DragPos.y - yPos) * (r.DragPos.y - yPos)) < (r.radius * r.radius)) {
                    dragSelected = &r;
                    break;
                }
                if (((r.PivotPos.x - xPos) * (r.PivotPos.x - xPos)) + ((r.PivotPos.y - yPos) * (r.PivotPos.y - yPos)) < (r.radius * r.radius)) {
                    dragSelected = &r;
                    break;
                }
            }
            for (auto& m : mVec) {
                if (!m.isplane) {
                    if (((m.edge1.x - xPos) * (m.edge1.x - xPos)) + ((m.edge1.y - yPos) * (m.edge1.y - yPos)) < (m.m_linewidth * m.m_linewidth)) {
                        mSelected = &m;
                        edge1selcted = true;
                        break;
                    }
                    if (((m.edge2.x - xPos) * (m.edge2.x - xPos)) + ((m.edge2.y - yPos) * (m.edge2.y - yPos)) < (m.m_linewidth * m.m_linewidth)) {
                        mSelected = &m;
                        edge1selcted = false;
                        break;
                    }
                }
            }

            for (auto& s : pVec) {
                //xpos and ypos is mouse position
                if (((s.edge1.x - xPos) * (s.edge1.x - xPos)) + ((s.edge1.y - yPos) * (s.edge1.y - yPos)) < (s.linewidthw * s.linewidthw)) {
                    planeselected = &s;
                    edge1s = true;
                    break;
                }
                if (((s.edge2.x - xPos) * (s.edge2.x - xPos)) + ((s.edge2.y - yPos) * (s.edge2.y - yPos)) < (s.linewidthw * s.linewidthw)) {
                    planeselected = &s;
                    edge2s = true;
                    break;
                }
                if (((s.edge3.x - xPos) * (s.edge3.x - xPos)) + ((s.edge3.y - yPos) * (s.edge3.y - yPos)) < (s.linewidthw * s.linewidthw)) {
                    planeselected = &s;
                    edge3s = true;
                    break;
                }
                if (((s.edge4.x - xPos) * (s.edge4.x - xPos)) + ((s.edge4.y - yPos) * (s.edge4.y - yPos)) < (s.linewidthw * s.linewidthw)) {
                    planeselected = &s;
                    edge4s = true;
                    break;
                }
            }
        }
        if (m_mouseButtonPressed[1]) {
            if (dragSelected != nullptr && mSelected == nullptr && planeselected == nullptr) {
                dragSelected->PivotPos.x = xPos;
                dragSelected->PivotPos.y = yPos;
            }
        }
        if (m_mouseButtonPressed[0]) {
            if (dragSelected != nullptr && mSelected == nullptr && planeselected == nullptr) {
                dragSelected->DragPos.x = xPos;
                dragSelected->DragPos.y = yPos;
            }
            if (mSelected != nullptr && dragSelected == nullptr && planeselected == nullptr) {
                if (edge1selcted) {
                    mSelected->edge1.x = xPos;
                    mSelected->edge1.y = yPos;
                }
                else {
                    mSelected->edge2.x = xPos;
                    mSelected->edge2.y = yPos;
                }
            }
            if (planeselected != nullptr && dragSelected == nullptr && mSelected == nullptr) {
                if (edge1s) {
                    planeselected->edge1.x = xPos;
                    planeselected->edge1.y = yPos;
                }
                if (edge2s) {
                    planeselected->edge2.x = xPos;
                    planeselected->edge2.y = yPos;
                }
                if (edge3s) {
                    planeselected->edge3.x = xPos;
                    planeselected->edge3.y = yPos;
                }
                if (edge4s) {
                    planeselected->edge4.x = xPos;
                    planeselected->edge4.y = yPos;
                }
            }
        }

        if (m_mouseButtonReleased[0] || m_mouseButtonReleased[1]) {
            dragSelected = nullptr;
            mSelected = nullptr;
            planeselected = nullptr;
        }
        cout << "number of rays: " << intersectionPoints.size() << '\n';
        intersectionPoints.clear();
    }
};

//quad tree tester
struct GridObj {
    Vector2D pos = Vector2D(-1,1);
    float h;
    int id;
};

class WaveGen {
public:
    //create the grid structure
    //first find a way to display the grid
    // the second value is the height of the grid at that point store the items as a long list
    //let the grid cells have a width and height of the screen or in this case the ndc coords set an inital w and h to 0.1
   
    //the pair stores the height and the cell type
    std::vector<GridObj> m_Grid;
    
    //so if i have a normalized width and height of 2 and i want to have 100 cells
    
    //allways initailize 100 cells world height is grid 1 starts at -1,1 and then 0.1 along is grid 2, and then 0.1 along is grid 3 and 
    void gridInit(int n) {
        Vector2D pos(-1, 1);
        for (int i = 0; i < 21; i++) {
            for (int j = 0; j < 21; j++) {
                pos.x += 0.1f;
                m_Grid[i].pos = pos;
                m_Grid.emplace_back({pos,100,1});
            }
            pos.x = -1;
            pos.y -= 0.1f;
        }
    }


    //change to draw the grid cells
    void DrawLine(Vector2D& c, Vector2D& q, float width, float x, float y, float z) {
        glEnable(GL_LINE_SMOOTH);
        glColor3f(x, y, z);
        glLineWidth(width);
        glBegin(GL_LINE_LOOP);
        glVertex2f(c.x, c.y);
        glVertex2f(q.x, q.y);
        glEnd();
        glDisable(GL_LINE_SMOOTH);
        glFlush();
    }
    void conditions() {

    }

    void loop() {

    }

};

class IdealGas {
public:
    int checks = 0;
    int size;
    int wc = 0;
    std::vector<particle> objects;
    std::vector<IdealBox> boxVec;
    QuadTreeContainer<particle*> grid;
    const int targetFPS = 100;
    const double targetFrameTime = 1.0 / targetFPS;
    TimeStep ts;
    bool use_heatmap = true;
    bool wasHKeyPressed = false;
    float maxT = -1.0f * FLT_MAX, minT = FLT_MAX;
    void DrawLine(Vector2D& c, Vector2D& q, float width, float x, float y, float z) {
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
        glVertex2f(c.x, c.y);
        glVertex2f(q.x, q.y);
        glEnd();
        glDisable(GL_LINE_SMOOTH);
        glFlush();
    }

    float gaussianKernel(double r, double h) {
        if (r >= 0 && r <= h) {
            double q = r / h;
            return (315.0 / (64.0 * PI * pow(h, 9))) * pow(h * h - r * r, 3);
        }
        else {
            return 0.0f;
        }
    }

    void CircleFunc(float xPoss, float yPoss, float rad, Vec3 c) {
        const int steps = 13;
        const float angle = 2.0f * PI / steps;
        glBegin(GL_TRIANGLE_FAN);
        glColor3f(c.x, c.y, c.z);
        glVertex2f(xPoss, yPoss);
        for (int i = 0; i <= steps; i++) {
            float nx = xPoss + rad * sin(angle * i);
            float ny = yPoss + rad * cos(angle * i);
            glVertex2f(nx, ny);
        }
        glEnd();
    }

    void processInput() {
        // Check for user input here
        if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS && !wasHKeyPressed) {
            use_heatmap = !use_heatmap;
            wasHKeyPressed = true;
        }
        if (glfwGetKey(window, GLFW_KEY_H) == GLFW_RELEASE) {
            wasHKeyPressed = false;
        }
    }

    Vec3 HeatGradient(float normalized_density, float min_density, float max_density) {
        float t = ((normalized_density - min_density) / (max_density - min_density));
        float r = t*t*1.2f;
        float b = 1.0f-t * t * 1.2f;
        return Vec3(r, 0.0f, b);
    }

    void normalizeVector(Vector2D& vector) {
        float nm = std::sqrtf(vector.x * vector.x + vector.y * vector.y);
        vector.x /= nm;
        vector.y /= nm;
    }
    void drawQTRect(const QTRect& rect) {
        glBegin(GL_LINE_LOOP);
        glLineWidth(0.1f);
        glColor4f(0.8f,0.7f,0.65f,0.5f);
        for (int i = 0; i < 4; ++i) {
            glVertex2f(rect.vertices[i].x, rect.vertices[i].y);
        }
        glEnd();
    }


    Vec3 GradientTemp(float val,float max,float min) {
        float t = (val - min) / (max-min);
        float r = t;
        float b = 1.0f-t;
        return Vec3(r,0,b);
    }

    void conditions() {
        ts.start();
        for (int i = 0; i < 1; i++) {
            std::random_device rd;
            std::mt19937 gen(rd());
            float lower_bound = 0.07f;
            float upper_bound = 0.08f;
            std::uniform_real_distribution<float> dist(lower_bound, upper_bound);
            float random_value = dist(gen);
            double lower_boundd = -0.45;
            double upper_boundd = 0.45;
            std::uniform_real_distribution<float> distr(lower_boundd, upper_boundd);
            float random_x = static_cast<float>(distr(gen));
            float random_y = static_cast<float>(distr(gen));
            particle pt;
            pt.radius = random_value;
            pt.mass = pt.radius * 1000.0f;
            pt.pos.x = random_x;
            pt.pos.y = random_y;
            pt.idx = i;
            pt.velocity = Vector2D(random_x+1, random_y-1);
            objects.emplace_back(pt);
        }
        IdealBox b(Vector2D(-0.4f,0.4f),Vector2D(0.4f,0.4f),Vector2D(0.4f,-0.4f),Vector2D(-0.4f,-0.4f),0.04f);
        b.rVel = Vector2D(0,0);
        b.initLength = b.getLength();
        boxVec.emplace_back(b);
        size = objects.size();
        cout << "Press H to toggle Kinetic Energy Heatmap!" << '\n';
    }
    float currTemp;
    float prevTemp;
    void loop() {
        //construct a search area  
        processInput();
        ts.update();
        grid.clear();
        std::map<int, bool> flag;
        double averageKE = 0;
        
        if (use_heatmap) {
            for (auto& box : boxVec) {
                Vector2D e1 = Vector2D(box.p1.x+0.01,box.p1.y);
                Vector2D e2 = Vector2D(box.p2.x-0.01,box.p2.y);
                DrawLine(e1, e2, 0.6f, 1.0f, 0.5f, 1.0f);
                //drawQTRect(box.r1);
                drawQTRect(box.r1);
                drawQTRect(box.r2);
                drawQTRect(box.r3);
                drawQTRect(box.r4);
                box.render();
            }
        }


        for (auto& objv : objects) {
            if (use_heatmap) {
                Vec3 vCol(1.0f,0.5f,0.25f);

                //cout << objv.velocity.x << "," << objv.velocity.y << endl;
                //CircleFunc(objv.pos.x, objv.pos.y, objv.radius, vCol);
                if (objv.highlight) {
                    CircleFunc(objv.pos.x, objv.pos.y, objv.radius, Vec3(1.0f,0,0));
                }
                else {
                    CircleFunc(objv.pos.x, objv.pos.y, objv.radius, Vec3(0,1.0f,0));
                }
                objv.setHighlight(false);
            }
            else {
                CircleFunc(objv.pos.x, objv.pos.y, objv.radius, Vec3(0.77f, 0.80f, 0.3f));
            }
            if (size < 1000 && !use_heatmap) {
                Vector2D velocity = objv.velocity;
                float nm = std::sqrtf(objv.velocity.x * objv.velocity.x + objv.velocity.y * objv.velocity.y);
                velocity /= nm;
                Vector2D q = objv.pos + velocity * objv.radius * 1.7f;
                DrawLine(objv.pos, q, 0.11f, 0.84f, 0.87f, 0.78f);
            }
            grid.insert(&objv, objv.getBoundingBox());
        }

       maxT = -1.0f * FLT_MAX;
       minT = FLT_MAX;
        if (use_heatmap) {
            for (int i = 0; i < size; i++) {
                if (flag[i] == false) {
                    for (const auto& obj : grid.search(objects[i].getSearchRadius())) {
                        if (flag[i] == false) {
                            if (objects[i] != *obj->item && objects[i].checkCollision(*obj->item)) {
                                objects[i].handleCollisions(obj->item);
                                flag[obj->item->idx] = true;
                            }
                        }
                    }
                }
                maxT = max(objects[i].gettempContribution(),maxT);
                minT = min(objects[i].gettempContribution(),minT);
            }

            for (int i = 0; i < size; i++) {
                for (auto& box : boxVec) {
                    if (box.detectCollision_(objects[i],box.p1,box.p2)) {
                        objects[i].setHighlight(true);
                        box.HandleCollision(objects[i], box.p1, box.p2);
                    }
                }
            }
        } else {
            for (int i = 0; i < size; i++) {
                if (flag[i] == false) {
                    for (const auto& obj : grid.search(objects[i].getSearchRadius())) {
                        if (objects[i] != *obj->item && objects[i].checkCollision(*obj->item)) {
                            objects[i].handleCollisions(obj->item);
                            flag[obj->item->idx] = true;
                        }
                    }
                }
            }
        }


        for (auto& obj : objects) obj.move(ts.dt());
        for (auto& b : boxVec)  {
            if (m_mouseButtonPressed[0]) {
                b.p1.x = xPos;
                b.p1.y = yPos;
            }
            b.move(ts.dt());
            b.applyStickConstraints(b.p1, b.p2);
            b.ShapeOverlap_SAT_STATIC(b.r1, b.r3);
            b.ShapeOverlap_SAT_STATIC(b.r1, b.r4);
            b.ShapeOverlap_SAT_STATIC(b.r1,b.r2);

        }
        double frameTime = glfwGetTime();
        double elapsedTime = frameTime - ts.dt();
        double remainingTime = targetFrameTime - elapsedTime;
        if (remainingTime > 0) {

            glfwWaitEventsTimeout(remainingTime);
        }
        wc++;
    }
};
//includes pendulum sim
/*class SpringSystem {
public:

};
//fluid simulation stuff
class FluidSystem {
public:
}; */



int main() {
    //newBallPhysics bf;
    //BallPhysics bll;
    //WaveFront wF;
    IdealGas g;
    /* Initialize the library */
    if (!glfwInit())
        return -1;
    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Physics Engine", NULL, NULL);
    glfwSetWindowPos(window, 300, 80);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    /* Make the window's context current */
    glfwMakeContextCurrent(window);
    if (glewInit() != GLEW_OK) {
        glfwTerminate();
        return -1;
    }


    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    //glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
    glEnable(GL_DEPTH_TEST);
    // Set the projection matrix
    glfwSetFramebufferSizeCallback(window, windowResizeCallback);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0); // Set the orthographic projection using NDC

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();


    // Enable anti-aliasing
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);


    //initialise the conditions of the program
    //wF.Conditions();
    //bll.Conditions();
    //bf.conditions();

    g.conditions();
    //qTree.resize(bd);

    //could remove this line
   
    // Insert all objects into the quadtree at the beginning of the simulation
    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


        
        //cout << "object count: " << count << endl;
        /* camera setup */

        // 
        //All Physics Functions here

        //loop func
       // std::chrono::steady_clock::time_point currentTime = std::chrono::steady_clock::now();
        //std::chrono::duration<double> elapsed = currentTime - previousTime;
        //felapsed = glfwGetTime() * 0.001;



        //.bf.loopbp();
        g.loop();
        //console view
        /*end = std::chrono::high_resolution_clock::now();

        // Calculate the elapsed time in seconds
        std::chrono::duration<float> duration = end - start;
        deltatime = duration.count();
        felapsed = duration.count()*0.001;
        elapsed += deltatime;*/
        //wF.loopWF();
        //bll.loopPH();
        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
    }



    glfwTerminate();
    return 0;
}
