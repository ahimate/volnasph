/*
 * Vector3D.h
 *
 *  Created on: 03.06.2010
 *      Author: ilya
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_
#include <cmath>

//Объект представляющий вектор или точку в пространстве
class Vector3D {
private:
	double x;
	double y;
	double z;
public:
	Vector3D();
	Vector3D(double x, double y, double z);
	double getX() const;
	double getY() const;
	double getZ() const;
	void setX(double x);
	void setY(double y);
	void setZ(double z);
	friend Vector3D operator+(Vector3D& a, Vector3D& b);
	friend Vector3D operator-(Vector3D& a, Vector3D& b);
	friend Vector3D operator*(Vector3D& a, double k);
	friend Vector3D operator*(double k, Vector3D& a);
	friend Vector3D scalarMul(Vector3D& a, Vector3D& b);
	friend Vector3D vectorMul(Vector3D& a, Vector3D& b);
	double abs(void) const;
	virtual ~Vector3D();
};

#endif /* VECTOR3D_H_ */
