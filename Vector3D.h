/*
 * Vector3D.h
 *
 *  Created on: 03.06.2010
 *  Updated on: 26.08.2010
 *      Author: Ilya Kryukov
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_
#include <cmath>

// Vector or point in 3d space
// aligned by 16 byte w - padding element
class Vector3D {
public:
	double x;
	double y;
	double z;
	double w; // padding

	//constructors
	Vector3D();
	Vector3D(double x, double y, double z);

	//getters-setters
	double getX() const;
	double getY() const;
	double getZ() const;
	void setX(double x);
	void setY(double y);
	void setZ(double z);

	double abs(void) const;

	friend Vector3D operator+(Vector3D& a, Vector3D& b);
	friend Vector3D operator-(Vector3D& a, Vector3D& b);
	friend Vector3D operator*(Vector3D& a, double k);
	friend Vector3D operator*(double k, Vector3D& a);

	friend Vector3D operator+=(Vector3D& a, Vector3D& b);
	friend Vector3D operator-=(Vector3D& a, Vector3D& b);
	//destructor
	virtual ~Vector3D();
};

//Math functions for vector
Vector3D scalarMul(Vector3D& a, Vector3D& b);
Vector3D vectorMul(Vector3D& a, Vector3D& b);

#endif /* VECTOR3D_H_ */
