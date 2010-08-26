/*
 * Vector3D.cpp
 *
 *  Created on: 03.06.2010
 *  Updated on: 26.08.2010
 *      Author: Ilya Kryukov
 */

#include "Vector3D.h"

Vector3D::Vector3D() 
{
	x = y = z = 0;
	w = 1;
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
	w = 1;
}

void Vector3D::setX(double x) {
	this->x = x;
}

Vector3D operator +(Vector3D& a, Vector3D& b) 
{
	return Vector3D(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector3D operator -(Vector3D& a, Vector3D& b) 
{
	return Vector3D(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector3D operator*(Vector3D& a, double k) 
{
	return Vector3D(a.x * k, a.y * k, a.z * k);
}

Vector3D operator*(double k, Vector3D& a) 
{
	return Vector3D(a.x * k, a.y * k, a.z * k);
}

Vector3D operator+=(Vector3D& a, Vector3D& b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}

Vector3D operator-=(Vector3D& a, Vector3D& b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

double Vector3D::abs(void) const 
{
	return sqrt(x * x + y * y + z * z);
}

Vector3D scalarMul(Vector3D& a, Vector3D& b) 
{
	return Vector3D(a.x * b.x, a.y * b.y, a.z * b.z);
}

Vector3D vectorMul(Vector3D& a, Vector3D& b) 
{
	return Vector3D(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

Vector3D::~Vector3D() {
}
