/*
 * Vector3D.cpp
 *
 *  Created on: 03.06.2010
 *      Author: ilya
 */

#include "Vector3D.h"

Vector3D::Vector3D() {
	// TODO Auto-generated constructor stub
	x = y = z = 0;
}

void Vector3D::setZ(double z) {
	this->z = z;
}

void Vector3D::setY(double y) {
	this->y = y;
}

double Vector3D::getZ() const {
	return z;
}

double Vector3D::getX() const {
	return x;
}

double Vector3D::getY() const {
	return y;
}

Vector3D::Vector3D(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

void Vector3D::setX(double x) {
	this->x = x;
}

Vector3D operator +(Vector3D& a, Vector3D& b) {
	return Vector3D(a.getX() + b.getX(), a.getY() + b.getY(), a.getZ()
			+ b.getZ());
}
Vector3D operator -(Vector3D& a, Vector3D& b) {
	return Vector3D(a.getX() - b.getX(), a.getY() - b.getY(), a.getZ()
			- b.getZ());
}
Vector3D operator*(Vector3D& a, double k) {
	return Vector3D(a.getX() * k, a.getY() * k, a.getZ() * k);
}
Vector3D operator*(double k, Vector3D& a) {
	return Vector3D(a.getX() * k, a.getY() * k, a.getZ() * k);
}

double Vector3D::abs(void) const {
	return sqrt(x * x + y * y + z * z);
}
Vector3D scalarMul(Vector3D& a, Vector3D& b) {
	return Vector3D(a.getX() * b.getX(), a.getY() * b.getY(), a.getZ()
			* b.getZ());
}
Vector3D vectorMul(Vector3D& a, Vector3D& b) {
	return Vector3D(a.getY() * b.getZ() - a.getZ() * b.getY(), a.getZ()
			* b.getX() - a.getX() * b.getZ(), a.getX() * b.getY() - a.getY()
			* b.getX());
}
Vector3D::~Vector3D() {
	// TODO Auto-generated destructor stub
}
