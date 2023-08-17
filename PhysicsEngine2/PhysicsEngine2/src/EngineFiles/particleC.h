#pragma once
#ifndef PARTICLEC_H
#include <iostream>
#include <cmath>
#include <immintrin.h> 
#include <vector>
#include "QTRstruct.h"
#include "particleC.h"
#include "Vector2D.h"
const double PI_V = 3.14159;

bool isBetween(double start, double middle, double end) {
	// Check if a value lies between two other values
	return start <= middle && middle <= end;
}

float GaussianKernel(float r, float h) {
	float q = r / h;
	return (1.0f / (PI_V * pow(h, 2))) * exp(-q * q);
}

float length(const Vector2D& v) {
	return std::sqrtf(v.x * v.x + v.y * v.y);
}
Vector2D normal(const Vector2D& l1, Vector2D& l2) {
	Vector2D dir = l2-l1;
	return Vector2D(-dir.y, dir.x);
}

float dot(const Vector2D& a, const Vector2D& b) {
	return a.x * b.x + a.y * b.y;
}

struct Rect {
	Vector2D e1;
	Vector2D e3;
	Vector2D e2;
	Vector2D e4;
};

struct Vec3 {
	float x;
	float y;
	float z;
	Vec3(float xi,float yi, float zi) :x(xi), y(yi), z(zi) {};
};
struct Vec4 {
	float x;
	float y;
	float z;
	float w;
	Vec4(float xi, float yi, float zi,float wi) :x(xi), y(yi), z(zi),w(wi) {};
};

struct particle {
	Vector2D pos;
	float radius;
	bool highlight = false;
	Vector2D velocity;
	Vector2D acceleration;
	//implement later
	//float sml = 3 * radius;
	float ela = 0.895f;
	float mass;
	float idx;
	float density;
	float p = 1.0f; 
	float temp;
	float getKineticE() const {
		float v = velocity.x * velocity.x + velocity.y * velocity.y;
		float actualmass = mass*1e-21f;
		return 0.5 *actualmass * v;
	};

	void conductHeat(particle*&other, float dt,float thermalConductivity) {
		float distance = sqrtf((pos.x - other->pos.x) * (pos.x - other->pos.x) + (pos.y - other->pos.y) * (pos.y - other->pos.y));; // Implement this function

	};

	float gettempContribution() const {
		float boltzman = 1.38064e-23;
		return ((2.0f * getKineticE()) / (3.0f * boltzman));
	};

	Vec3 KE_col() const {
		float x = 1 - getKineticE();
		return Vec3(1.0f, x, x);
	};
	QTRect getBoundingBox() const {
		Vector2D position(pos.x, pos.y);
		Vector2D size(2 * radius, 2 * radius);
		QTRect bs;
		bs.pos = position;
		bs.dimensions = size;
		return bs;
	};
	QTRect getSearchRadius() const {
		Vector2D position(pos.x, pos.y);
		Vector2D size(3 * radius,3 * radius);
		QTRect bs;
		bs.pos = position;
		bs.dimensions = size;
		return bs;
	};

	void move(float dt) {
		
		//acceleration.x = -velocity.x * 0.99f;
		//acceleration.y = -velocity.y * 0.99f;
		//init = 0.033
		//velocity.x += acceleration.x * dt;
		//velocity.y += acceleration.y * dt;

		pos.x += velocity.x * dt;
		pos.y += velocity.y * dt;
		//test bool overlap theory tomorrow
		if (pos.x <= -p + radius) {
			pos.x = -p + radius;
			velocity.x *= -1;

		}
		if (pos.x >= p - radius) {
			pos.x = p - radius;
			velocity.x *= -1;
		}
		if (pos.y <= -p + radius) {
			pos.y = -p + radius;
			velocity.y *= -1;
		}
		if (pos.y >= p - radius) {
			pos.y = p - radius;
			velocity.y *= -1;
		}
	
		/*if (fabs(velocity.x * velocity.x + velocity.y * velocity.y) < fStable)
		{
			velocity.x = 0;
			velocity.y = 0;
		}*/
		
	};
	bool checkCollision(particle& other) {
		const float distft = (pos.x - other.pos.x) * (pos.x - other.pos.x) + (pos.y - other.pos.y) * (pos.y - other.pos.y);
		const float radiusSumSquared = (radius + other.radius) * (radius + other.radius);
		return (distft <= radiusSumSquared);
	};

	void handleCollisions(particle*& other) {
		float Distance = sqrtf((pos.x - other->pos.x) * (pos.x - other->pos.x) + (pos.y - other->pos.y) * (pos.y - other->pos.y));
		float fDistance = 1 / Distance;
		float fOverlap = 0.5f * (Distance - radius - other->radius);
		float nx = (pos.x - other->pos.x) * fDistance;
		float ny = (pos.y - other->pos.y) * fDistance;
		float tx = -ny;
		float ty = nx;
		float dpTan1 = other->velocity.x * tx + other->velocity.y * ty;
		float dpTan2 = velocity.x * tx + velocity.y * ty;
		float dpNorm1 = other->velocity.x * nx + other->velocity.y * ny;
		float dpNorm2 = velocity.x * nx + velocity.y * ny;
		float m1 = 1.00f * (dpNorm1 * (other->mass - mass) + 2.0f * mass * dpNorm2) / (other->mass + mass);
		float m2 = 1.00f * (dpNorm2 * (mass - other->mass) + 2.0f * other->mass * dpNorm1) / (mass + other->mass);
		//add an efficiency slider thing
		other->velocity.x = (tx * dpTan1 + nx * m1);
		other->velocity.y = (ty * dpTan1 + ny * m1);
		velocity.x = (tx * dpTan2 + nx * m2);
		velocity.y = (ty * dpTan2 + ny * m2);
		pos.x -= fOverlap * (pos.x - other->pos.x) * fDistance;
		pos.y -= fOverlap * (pos.y - other->pos.y) * fDistance;
		other->pos.x += fOverlap * (pos.x - other->pos.x) * fDistance;
		other->pos.y += fOverlap * (pos.y - other->pos.y) * fDistance;
	 

	};
	//probaly a smarter check but oh well
	bool operator!=(const particle& other) const {
		return (pos.x != other.pos.x) || (pos.y != other.pos.y);
	};
	float dist(float x1, float y1, float x2, float y2) {
		float deltaX = x2 - x1;
		float deltaY = y2 - y1;
		float distance = std::sqrt(deltaX * deltaX + deltaY * deltaY);
		return distance;
	};

	bool intersects(const particle& other) const {
		double dx = other.pos.x - pos.x;
		double dy = other.pos.y - pos.y;
		double distanceSq = dx * dx + dy * dy;
		double radiusSum = radius + other.radius;

		return distanceSq < radiusSum * radiusSum;
	}

	void setHighlight(bool val) {
		highlight = val;
	};

};

struct IdealBox {
	Vector2D p1, p2, p3, p4, v1, v2, v3, v4;
	Vector2D l1Centre;
	Vector2D rVel = Vector2D(0.00001f,0.00001f);
	Vector2D v1a;
	Vector2D v2a;
	Vector2D cPoint;
	Vector2D op1;
	Vector2D op2;
	Vector2D globalnormal;
	//the 4 rectangles
	QTRect r1;
	QTRect r2;
	QTRect r3;
	QTRect r4;
	//define the other 4 as rectangles define one as the wierd normal ting

	float p = 1.0f;
	float initLength;
	bool isOpen = true;
	float Pressure;
	float thickness;
	//adjust to see if it reacts with large balls
	float mass = 550.0f;
	//line segements p1-p2 and p2-p3 and p3-p2 and p1-p4
	IdealBox(Vector2D edge1,Vector2D edge2,Vector2D edge3,Vector2D edge4,float f) : p1(edge1), p2(edge2),p3(edge3),p4(edge4),thickness(f) {
		op1 = p1;
		op2 = p2;
		p1.x += 0.08f;
		p2.x -= 0.08f;
		p1.y += 0.05f;
		p2.y += 0.05f;
		r1.pos = 0.5f * (p1 + p2);
		r2.pos = 0.5f * (p3 + op2);
		r3.pos = 0.5f * (p3 + p4);
		r4.pos = 0.5f * (op1 + p4);
		Vector2D edge_ = p2 - p1;
		globalnormal.x = -edge_.y;
		globalnormal.y = edge_.x;
		float nw = std::sqrtf(globalnormal.x*globalnormal.x + globalnormal.y*globalnormal.y);
		globalnormal /= nw;
		globalnormal *= 0.008f;
		//the p1 edge
		r1.vertices[0] = p1 - globalnormal;
		r1.vertices[1] = p2 - globalnormal;
		r1.vertices[2] = p2 + globalnormal;
		r1.vertices[3] = p1 + globalnormal;

		r2.vertices[0] = Vector2D(op2.x - thickness, op2.y);
		r2.vertices[1] = Vector2D(p3.x - thickness, p3.y);
		r2.vertices[2] = Vector2D(p3.x + thickness, p3.y);
		r2.vertices[3] = Vector2D(op2.x + thickness, op2.y);

		r3.vertices[0] = Vector2D(p3.x, p3.y - thickness);
		r3.vertices[1] = Vector2D(p4.x, p4.y - thickness);
		r3.vertices[2] = Vector2D(p4.x, p4.y + thickness);
		r3.vertices[3] = Vector2D(p3.x, p3.y + thickness);

		r4.vertices[0] = Vector2D(p4.x - thickness, p4.y);
		r4.vertices[1] = Vector2D(op1.x - thickness, op1.y);
		r4.vertices[2] = Vector2D(op1.x + thickness, op1.y);
		r4.vertices[3] = Vector2D(p4.x + thickness, p4.y);
		
	};
	

	float getLength() {
		return std::sqrtf((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
	};

	void applyStickConstraints(Vector2D s1, Vector2D s2) {
		double dx = s1.x - s2.x;
		double dy = s1.y - s2.y;
		double distance = std::sqrt(dx * dx + dy * dy);
		double difference = initLength - distance;
		double percent = difference / distance / 2;
		double offsetX = dx * percent;
		double offsetY = dy * percent;
		p2 -= Vector2D(offsetX, offsetY);
		p1 += Vector2D(offsetX, offsetY);
	};

	void move(float dt) {

		v1a.x = -v1.x * 0.9f;
		v1a.y = -v1.y * 0.9f-0.5f;
		// -0.689f is the grav constant
		v2a.x= -v2.x * 0.9f;
		v2a.y = -v2.y * 0.9f-0.5f;

		v1 += v1a * dt;
		v2 += v2a * dt;

		p1 += v1 * dt;
		p2 += v2 * dt;

		if (p1.x <= -p) {
			p1.x = -p;
			v1.x *= -0.95f;

		}
		if (p1.x >= p) {
			p1.x = p;
			v1.x *= -0.95f;
		}
		if (p1.y <= -p) {
			p1.y = -p;
			v1.y *= -0.95f;
		}
		if (p1.y >= p) {
			p1.y = p;
			v1.y *= -0.95f;
		}
		//p2bit
		if (p2.x <= -p) {
			p2.x = -p;
			v2.x *= -0.95f;

		}
		if (p2.x >= p) {
			p2.x = p;
			v2.x *= -0.95f;
		}
		if (p2.y <= -p) {
			p2.y = -p;
			v2.y *= -0.95f;
		}
		if (p2.y >= p) {
			p2.y = p;
			v2.y *= -0.95f;
		}
		//but the rest here...
		//put the constraints here
		Vector2D edge_m1 = p2 - p1;
		globalnormal.x = -edge_m1.y;
		globalnormal.y = edge_m1.x;
		float nl = std::sqrtf(globalnormal.x * globalnormal.x + globalnormal.y * globalnormal.y);
		globalnormal /= nl;
		globalnormal *= 0.008;
		r1.vertices[0] = p1 - globalnormal;
		r1.vertices[1] = p2 - globalnormal;
		r1.vertices[2] = p2 + globalnormal;
		r1.vertices[3] = p1 + globalnormal;
		r1.pos = 0.5f * (p1 + p2);
	};
	
	bool detectCollision_(particle& target, Vector2D& e1, Vector2D& e2) {
		//stuff here 
		Vector2D lineDir = Vector2D(e2.x - e1.x, e2.y - e1.y);
		Vector2D toTarget = Vector2D(target.pos.x - e1.x, target.pos.y - e1.y);
		double lineLength = std::sqrt(lineDir.x * lineDir.x + lineDir.y * lineDir.y);
		Vector2D lineDirNormalized = Vector2D(lineDir.x / lineLength, lineDir.y / lineLength);
		double projection = toTarget.x * lineDirNormalized.x + toTarget.y * lineDirNormalized.y;
		projection = std::max(0.0, std::min(projection, lineLength));
		Vector2D closestPoint = Vector2D(e1.x + lineDirNormalized.x * projection, e1.y + lineDirNormalized.y * projection);
		cPoint = closestPoint;
		double distance = std::sqrt((closestPoint.x - target.pos.x) * (closestPoint.x - target.pos.x) +(closestPoint.y - target.pos.y) * (closestPoint.y - target.pos.y));
		if (distance <= target.radius) {
			return true;
		}
		return false;
	};
	//handle the collisions
	
	void HandleCollision(particle& target, Vector2D& m_edge1, Vector2D& m_edge2) {
		//handle the scaling take edge1 to be A and edge2 to be B
		//Vector2D centre = 0.5f*(m_edge2 + m_edge1); //let this be x2 and target centre x1
		float a_m = std::sqrtf((m_edge1.x-cPoint.x)*(m_edge1.x-cPoint.x)+(m_edge1.y-cPoint.y)*(m_edge1.y-cPoint.y));
		float b_m = std::sqrtf((m_edge2.x - cPoint.x) * (m_edge2.x - cPoint.x) + (m_edge2.y - cPoint.y) * (m_edge2.y - cPoint.y));
		float mu_a = b_m/initLength;
		float mu_b = a_m/initLength;
		//distance between the centres
		//retard there is no distance between centre
		float dist = std::sqrtf((cPoint.x-target.pos.x)*(cPoint.x-target.pos.x) + (cPoint.y-target.pos.y)*(cPoint.y-target.pos.y));
		float fdist = 1 / dist;
		float overlapDistance = target.radius - dist;
		Vector2D n = target.pos-cPoint;
		n *= fdist;
		Vector2D t;
		t.x = -n.y;
		t.y = n.x;
		//normal and edge
		//Vector2D t = m_edge2 - m_edge1;
		//float nm = std::sqrtf(t.x*t.x + t.y*t.y);
		//t /= nm;
		//Vector2D n = Vector2D(-1.0f * t.y,t.x);
		float dpTan1 = target.velocity.x * t.x + target.velocity.y * t.y;
		float dpTan2 = rVel.x * t.x + rVel.y * t.y;
		float dpNorm1 = target.velocity.x * n.x + target.velocity.y * n.y;
		float dpNorm2 = rVel.x* n.x + rVel.y * n.y;
		float m1 = 1.00f * (dpNorm1 * (target.mass - mass) + 2.0f * mass * dpNorm2) / (target.mass + mass);
		float m2 = 1.00f * (dpNorm2 * (mass - target.mass) + 2.0f * target.mass * dpNorm1) / (mass + target.mass);
		target.velocity.x = (t.x * dpTan1 + n.x * m1);
		target.velocity.y = (t.y * dpTan1 + n.y * m1);
		rVel.x = (t.x * dpTan2 + n.x * m2);
		rVel.y = (t.y * dpTan2 + n.y * m2);
		//the edge1 vel
		v1 = rVel * mu_a;
		//the edge2 vel
		v2 = rVel * mu_b;
		 // Calculate the overlap distance
		target.pos += n * overlapDistance;
		m_edge1 -= n * (overlapDistance * mu_a); // Adjust edge1 based on mu_a
		m_edge2 -= n * (overlapDistance * mu_b);
	};
	//render the box set edge thickness
	bool doIntersect_l(const Vector2D& wp1, const Vector2D& q1,const Vector2D& wp2, const Vector2D& q2,Vector2D& intersection) {
		Vector2D v1 = q1 - wp1;
		Vector2D v2 = q2 - wp2;

		double cross1 = v1 * (wp2 - wp1);
		double cross2 = v1 * (q2 - wp1);
		double cross3 = v2 * (wp1 - wp2);
		double cross4 = v2 * (q1 - wp2);

		if ((cross1 * cross2 <= 0) && (cross3 * cross4 <= 0)) {
			// Calculate the intersection point
			double t = cross1 / (cross1 - cross2);
			intersection = wp2 + v2 * t;
			return true;
		}

		return false;
	};

	void resolveCol_w(Vector2D& r1, Vector2D& r2, Vector2D r3, Vector2D r4,int sideID) {
		Vector2D intersectionPoint;
		if (doIntersect_l(r1,r2,r3,r4,intersectionPoint)) {
			//define normal
			std::cout << "wall intersection happening!" << std::endl;
			//bit of a silly method but i dont see the point in doing SAT or some other collision response stuff cuz i dont need it 
			/*if (sideID == 2) {
				if (r1.x >= p3.x) {
					r1.x = p3.x;
					v1.x *= -0.95f;
				}
				if (r2.x >= p3.x) {
					r2.x = p3.x;
					v2.x *= -0.95f;
				}
			}
			if (sideID == 3) {

			}
			if (sideID == 4) {

			}*/
			float a_m = std::sqrtf((r1.x - intersectionPoint.x) * (r1.x - intersectionPoint.x) + (r1.y - intersectionPoint.y) * (r1.y - intersectionPoint.y));
			float b_m = std::sqrtf((r2.x - intersectionPoint.x) * (r2.x - intersectionPoint.x) + (r2.y - intersectionPoint.y) * (r2.y - intersectionPoint.y));
			float lengthO = std::sqrtf((r4.x-r3.x)*(r4.x-r3.x)+(r4.y-r3.y)*(r4.y-r3.y));
			float mu_a = b_m / lengthO;
			float mu_b = a_m / lengthO;
			Vector2D q = r4 - r3;
			Vector2D normal;
			normal.x = -q.y;
			normal.y = q.x;
			normal /= lengthO;
			float dotProduct = rVel.x * normal.x + rVel.y * normal.y;
			Vector2D prevVel = rVel;
			rVel.x = prevVel.x - (2 * dotProduct * normal.x);
			rVel.y = prevVel.y - (2 * dotProduct * normal.y);
			v1 = rVel * mu_a;
			v2 = rVel * mu_b;


		}
	};
	//second attempt at an intersection function
	bool findRectangleCollision(const QTRect& rect1, const QTRect& rect2,float& overlap, Vector2D& collisionPoint, Vector2D& collisionNormal) {
		float minOverlap = INFINITY;
		bool collision = true; // Assume collision by default

		for (int i = 0; i < 4; ++i) {
			Vector2D edge = rect1.vertices[(i + 1) % 4] - rect1.vertices[i];
			Vector2D normal(-edge.y, edge.x);  // Perpendicular normal (outward)

			// Normalize the normal
			float len = std::sqrt(normal.x * normal.x + normal.y * normal.y);
			normal.x /= len;
			normal.y /= len;

			// Project vertices onto the normal
			float min1 = INFINITY, max1 = -INFINITY;
			float min2 = INFINITY, max2 = -INFINITY;

			for (int j = 0; j < 4; ++j) {
				float projection = rect1.vertices[j].x * normal.x + rect1.vertices[j].y * normal.y;
				min1 = std::min(min1, projection);
				max1 = std::max(max1, projection);
			}

			for (int j = 0; j < 4; ++j) {
				float projection = rect2.vertices[j].x * normal.x + rect2.vertices[j].y * normal.y;
				min2 = std::min(min2, projection);
				max2 = std::max(max2, projection);
			}

			// Check for separation
			if (max1 < min2 || max2 < min1) {
				collision = false;
				break;
			}

			// Calculate overlap
			float o = std::min(max1, max2) - std::max(min1, min2);
			if (o < minOverlap) {
				minOverlap = o;
				collisionNormal = normal;
			}
		}

		if (collision) {
			overlap = minOverlap;

			// Calculate point of intersection
			collisionPoint = rect1.pos + collisionNormal * (overlap * 0.5f);

			return true;
		}

		return false;
	};

	void handleRectCol(QTRect& rect1, QTRect& rect2) {
		float overlap;
		Vector2D collisionNorm;
		Vector2D colPoint;
		if (findRectangleCollision(rect1, rect2, overlap, colPoint, collisionNorm)) {
			float moveDistance = overlap * 0.5f;
			std::cout << "wall intersection happening!" << std::endl;
			// Move rect1 away from rect2 along collisionNormal
	
			/*only use this section if i need dynamic rectangle collsion*/
			// Move rect2 away from rect1 along collisionNormal
			//rect2.pos.x += collisionNorm.x * moveDistance;
			//rect2.pos.y += collisionNorm.y * moveDistance;

			
			//float a_m = std::sqrtf((p1.x - colPoint.x) * (p1.x - colPoint.x) + (p1.y - colPoint.y) * (p1.y - colPoint.y));
			//float b_m = std::sqrtf((p2.x - colPoint.x) * (p2.x - colPoint.x) + (p2.y - colPoint.y) * (p2.y - colPoint.y));
			//float mu_a = b_m / initLength;
			//float mu_b = a_m / initLength;
			//float dotProduct = rVel.x * collisionNorm.x + rVel.y * collisionNorm.y;
			//rVel.x = rVel.x - (2 * dotProduct * collisionNorm.x);
			//rVel.y = rVel.y - (2 * dotProduct * collisionNorm.y);
			//v1 = rVel * mu_a * 0.2f;
			//v2 = rVel * mu_b * 0.2f;
			 
			p1.x -= collisionNorm.x * moveDistance;
			p1.y -= collisionNorm.y * moveDistance;
			p2.x -= collisionNorm.x * moveDistance;
			p2.y -= collisionNorm.y * moveDistance;
		}
	};

	void CircleFunc2(float xPoss, float yPoss, float rad, Vec3 c) {
		const int steps = 13;
		const float angle = 2.0f * 3.14159 / steps;
		glBegin(GL_TRIANGLE_FAN);
		glColor3f(c.x, c.y, c.z);
		glVertex2f(xPoss, yPoss);
		for (int i = 0; i <= steps; i++) {
			float nx = xPoss + rad * sin(angle * i);
			float ny = yPoss + rad * cos(angle * i);
			glVertex2f(nx, ny);
		}
		glEnd();
	};
	//might still use but too tedious
	bool ShapeOverlap_SAT_STATIC(QTRect& pwr1, QTRect& pwr2) {
		Vector2D collisionPoint;
		QTRect* rect1 = &pwr1;
		QTRect* rect2 = &pwr2;
		Vector2D axisProj;
		float overlap = INFINITY;

		for (int shape = 0; shape < 2; shape++) {
			if (shape == 1) {
				std::swap(rect1, rect2); // Swap rectangles for second loop iteration
			}

			for (int a = 0; a < 4; a++) {
				int b = (a + 1) % 4;
				Vector2D edge = rect1->vertices[b] - rect1->vertices[a];
				axisProj = Vector2D(-edge.y, edge.x);
				float nkl = std::sqrtf(axisProj.x * axisProj.x + axisProj.y * axisProj.y);
				axisProj /= nkl;
				// Calculate min and max 1D projections for rect1 along axisProj
				float min_r1 = INFINITY, max_r1 = -INFINITY;
				for (int p = 0; p < 4; p++) {
					float q = (rect1->vertices[p].x * axisProj.x + rect1->vertices[p].y * axisProj.y);
					min_r1 = std::min(min_r1, q);
					max_r1 = std::max(max_r1, q);
				}

				// Calculate min and max 1D projections for rect2 along axisProj
				float min_r2 = INFINITY, max_r2 = -INFINITY;
				for (int p = 0; p < 4; p++) {
					float q = (rect2->vertices[p].x * axisProj.x + rect2->vertices[p].y * axisProj.y);
					min_r2 = std::min(min_r2, q);
					max_r2 = std::max(max_r2, q);
				}

				// Calculate actual overlap along projected axis, and store the minimum
				overlap = std::min(std::min(max_r1, max_r2) - std::max(min_r1, min_r2), overlap);
				if (overlap > 0.0f) {
					collisionPoint = (rect1->pos + rect2->pos) * 0.5f + axisProj * (overlap * 0.5f);
				}
				// Check for separation along the current axis
				if (!(max_r2 >= min_r1 && max_r1 >= min_r2)) {
					return false; // No collision
				}
			}
		}
		// If we got here, the objects have collided. Displace r1 by overlap along the vector between the two object centers.
		
		//THE DISTANCE TOO CENTRE METHOD DOESNT WORK RETRy
		
		CircleFunc2(collisionPoint.x,collisionPoint.y,0.03f,Vec3(0.5f,0.25f,0.9f));


		//Vector2D d = Vector2D(rect1->pos.x - rect2->pos.x, rect1->pos.y - rect2->pos.y);
		//Vector2D d = Vector2D(r1.pos.x - pwr2.pos.x,r1.pos.y - pwr2.pos.y);
		Vector2D d = Vector2D(rect1->pos.x - rect2->pos.x, rect1->pos.y - rect2->pos.y);
		//d = Vector2D(-d.y,d.x);
		//std::cout << "the object centres so r2 first or wall: " << r2.pos.x << "," << r2.pos.y << " now R1: " << r1.pos.x << "," << r1.pos.y << std::endl;
		float s = sqrtf(d.x * d.x + d.y * d.y);
		d /= s;
		//overlap += 0.0005;
		//std::cout << 1 << std::endl;
		//std::cout << rect2->pos.x << "," << rect2->pos.y << std::endl;
		//distance between the centres


		//FOR D1 and D2 use the perpendicular distance between two lines to figure out which one is in contact or which one to make a response
		float d1 = std::sqrtf((p1.x- pwr2.pos.x)*(p1.x - pwr2.pos.x) + (p1.y - pwr2.pos.y)*(p1.y - pwr2.pos.y));
		float d2 = std::sqrtf((p2.x- pwr2.pos.x)*(p2.x - pwr2.pos.x) + (p2.y - pwr2.pos.y)*(p2.y - pwr2.pos.y));
		//float d1 = std::sqrtf((p1.x- collisionPoint.x)*(p1.x - collisionPoint.x) + (p1.y - collisionPoint.y)*(p1.y - collisionPoint.y));
		//float d2 = std::sqrtf((p2.x- collisionPoint.x)*(p2.x - collisionPoint.x) + (p2.y - collisionPoint.y)*(p2.y - collisionPoint.y));
		/*if (d1 < d2) {
			p1 -= overlap * d;
			if (pwr2.pos == r2.pos) {
				//std::cout << "r2" << std::endl;
				v1.x *= -0.95f;
			}
			if (pwr2.pos == r3.pos) {
				//std::cout << "r3" << std::endl;
				v1.y *= -0.95f;
			}
			if (pwr2.pos == r4.pos) {
				//std::cout << "r4" << std::endl;
				v1.x *= -0.95f;
			}

		}
		else {
			p2 -= overlap * d;
			if (pwr2.pos == r2.pos) {
				//std::cout << "r2" << std::endl;
				v2.x *= -0.95f;
			}
			if (pwr2.pos == r3.pos) {
				//std::cout << "r3" << std::endl;
				v2.y *= -0.95f;
			}
			if (pwr2.pos == r4.pos) {
				//std::cout << "r4" << std::endl;
				v2.x *= -0.95f;
			}
		}*/
		std::cout <<d.x << "," << d.y << " now with overlap: " << overlap * d.x << "," << overlap*d.y << std::endl;
		p1 -= overlap * d;
		p2 -= overlap * d;
		return true; // Collision detected
	};

	/*
	bool findCollisionPoint(const Vector2D& vertexPos, const Vector2D& vertexVelocity,
		const Vector2D& wallStart, const Vector2D& wallEnd, Vector2D& collisionPoint) {
		// Calculate the direction vector of the line segment
		Vector2D wallDirection = wallEnd - wallStart;

		// Calculate the direction vector of the vertex's trajectory
		Vector2D vertexDirection = vertexVelocity.normalize();

		// Calculate the dot product between the two direction vectors
		float dotProduct = vertexDirection.x * wallDirection.y - vertexDirection.y * wallDirection.x;

		// Check if the dot product is close to zero, indicating parallel lines
		if (std::abs(dotProduct) < 0.0001f) {
			return false; // No collision
		}

		// Calculate the vector from wall start to vertex position
		Vector2D wallToVertex = vertexPos - wallStart;

		// Calculate the parameter t for the collision point
		float t = (wallToVertex.x * wallDirection.y - wallToVertex.y * wallDirection.x) / dotProduct;

		// Check if t is within the range [0, 1], indicating that the collision point lies on the line segment
		if (t >= 0.0f && t <= 1.0f) {
			collisionPoint = wallStart + t * wallDirection;
			return true; // Collision detected
		}

		return false; // No collision
	}

	void handleWallCollisions(const QTRect& rect1,const Vector2D& wallStart, const Vector2D& wallEnd, float restitution) {
		for (int i = 0; i < 4; ++i) {
			Vector2D vertexPos(rect1.vertices[i].x, rect1.vertices[i].y);
			Vector2D collisionPoint(0.0f, 0.0f);

			if (findCollisionPoint(vertexPos, vertexVelocity, wallStart, wallEnd, collisionPoint)) {
				// Calculate collision normal
				Vector2D collisionNormal = (collisionPoint - vertexPos).normalize();

				// Calculate relative velocity
				float relativeVelocity = vertexVelocity.dot(collisionNormal);

				// Calculate impulse
				float impulse = -(1.0f + restitution) * relativeVelocity;
				Vector2D impulseVector = collisionNormal * impulse;

				// Apply impulse to velocity
				vertices[i].vx += impulseVector.x;
				vertices[i].vy += impulseVector.y;
			}
		}
	}*/

	void render() {

		if (isOpen) {
			glColor3f(1.0f, 1.0f, 1.0f);  // Set color to white
			glLineWidth(thickness);
			glBegin(GL_LINES);

			glColor3f(1.0f,0.5f,1.0f);
			glVertex2f(p3.x, p3.y);  // Top-left vertex
			glVertex2f(op2.x, op2.y);   // Top-right vertex
			glEnd();
			glBegin(GL_LINES);
			glVertex2f(p3.x, p3.y);  // Top-left vertex
			glVertex2f(p4.x, p4.y);   // Top-right vertex
			glEnd();
			glBegin(GL_LINES);
			glVertex2f(op1.x, op1.y);  // Top-left vertex
			glVertex2f(p4.x, p4.y);   // Top-right vertex
			glEnd();
			glLineWidth(1.0f);
		}
		else {
			glLineWidth(thickness);
			glBegin(GL_LINE_LOOP);
			glVertex2f(p1.x, p1.y);
			glVertex2f(p2.x, p2.y);
			glVertex2f(p3.x, p3.y);
			glVertex2f(p4.x, p4.y);
			glEnd();
			glLineWidth(1.0);
		}
	}
};



#endif // PARTICLEC_H

