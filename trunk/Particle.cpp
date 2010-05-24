#include "Particle.h"
#include <cmath>
#include <iostream>
#include <limits>
using namespace std;
double Pi=3.1415926358;
double koefWroo=1.566681480;////315/(64*Pi);
double koefWpres=14.32394496;////45/(Pi
//////////����
double Particle::fluM=0.06; // �����
double Particle::fluNu=3.5; //���� ��������
double Particle::fluROO=998.29; //��������� ���������
double Particle::fluK=3; //���� � ����� ���������
double Particle::wallM=0.04; // �����

Particle::Particle(double cpX, double cpY, double cpZ)
{
	cX=cpX;
	cY=cpY;
	cZ=cpZ;
	cU=0;
	cV=0;
	cW=0;
	cRo=0;////////////////////////////////////////
	cP=0;//////////////////////////////////////////
	type=0;
}
Particle::Particle(double x,double y,double z,double u,double v,double w)
{
	cX=x;
	cY=y;
	cZ=z;
	cU=u;
	cV=v;
	cW=w;
	cRo=fluROO;
	cP=0;
	type=1;
}


double Particle::parDistanse(Particle& a,Particle& b)
{
	double tmp1=pow(b.cX-a.cX,2);
	double tmp2=pow(b.cY-a.cY,2);
	double tmp3=pow(b.cZ-a.cZ,2);
	double tmp4=tmp1+tmp2+tmp3;
	return pow(tmp4,0.5);
};

double Particle::fluFpress(Particle& a, Particle& b, double re,double ba)
{
	double modR=Particle::parDistanse(a,b);
	if(modR>re)
	{
		return 0;
	}
	else
	{
		double tmp=koefWpres/(pow(re,6));
		double Wpress=tmp*pow((re-modR),2)*ba/modR;
		double tmp10=a.cP/(pow(a.cRo,2));
		double tmp20=b.cP/(pow(b.cRo,2));
		return b.fluM*Wpress*(tmp10+tmp20);
	}
};
double Particle::fluFvis(Particle& a, Particle& b, double re,double ba)
{
	double modR=Particle::parDistanse(a,b);
	if(modR>re)
	{
		return 0;
	}
	else
	{
		double tmp=koefWpres/(pow(re,6));
		double Wvis=tmp*(re-modR);
		return b.fluM*Wvis*ba/b.cRo;
	}
};
double Particle::fluRoWall(double riw,double wallmass,double re)///////////////////��������
{
	double tmp=315/(64*Pi*pow(re,9));
	double tmp2=pow(riw,2);
	double W=tmp*(pow(abs(re*re-tmp2),3));//������
	return wallmass*W;
};
double Particle::fluRo(Particle& a,Particle& b, double re)
{
	double modR=a.parDistanse(a,b);
	if(modR>re)
	{
		return 0;
	}
	else
	{
		double tmp=koefWroo/pow(re,9);
		double tmp3=pow(modR,2);
		double W=tmp*(pow((re*re-tmp3),3));
		return b.fluM*W;
	}
};


/*
double static fluFpress(Fluid& a, Fluid& b, double re,double ba); //������� �������� ba=b.cX(y,z)-a.cX(y,z, a=i,b=j
double static fluFvis(Fluid& a, Fluid& b, double re,double ba); //������� ��� �������� ba=b.cU(v,w)-a.cU(v,w), a=i,b=j
double static fluRo(Fluid& a,Fluid& b, double re); //������� ��� ������� ��������� a-�������,b-�� ������

static double fluM; // �����
static double fluNu; //���� ��������
static double fluROO; //��������� ���������
static double fluPO; //��������� ��������
static double fluK; //���� � ����� ���������*/


//get

/*double Particle::getX()
{
return cX;
}; 
double Particle::getY()
{
return cY;
}; 
double Particle::getZ()
{
return cZ;
}; 

double Fluid::getU()
{
return cU;
}; 
double Fluid::getV()
{
return cV;
}; 
double Fluid::getW()
{
return cW;
}; 
double Fluid::getRo()
{
return cRo;
};
double Fluid::getP()
{
return fluPO+fluK*(cRo-fluROO);
};
//set
void Fluid::setU(double Spe)
{
cU=Spe;
}; 
void Fluid::setV(double Spe)
{
cV=Spe;
}; 
void Fluid::setW(double Spe)
{
cW=Spe;
}; 
void Fluid::setRo(double Den)
{
cRo=Den;
}; 
void Fluid::setP(double press)
{
cP=press;
}; 
void Particle::setX(double Spe)
{
cX=Spe;
}; 
void Particle::setY(double Spe)
{
cY=Spe;
}; 
void Particle::setZ(double Spe)
{
cZ=Spe;
}; 
*/