#include <vector>
#include "Particle.h"
#include <fstream> 
using namespace std;
class Bucket
{
	public:
	vector<Particle> part;
	bool status;//1-пустой
	double bX; // х-координата начала €чейки 
	double bY; // y-координата начала €чейки
	double bZ; // z-координата начала €чейки
// constructor
	Bucket(){bX=0;bY=0;bZ=0;}
	Bucket(double buX,double buY,double buZ);
	
	void BucketInFile(ofstream &out);
	void BucketVelInFile(ofstream &out);
	void BucketRoInFile(ofstream &out);
	void BucketPInFile(ofstream &out);
	
	void FluInFile(ofstream &out);
	void FluVelInFile(ofstream &out);
	void FluRoInFile(ofstream &out);
	void FluPInFile(ofstream &out);

	double bucRiw(Particle& a);
	double bucDen(Particle& a,double re);
	void bucForse(Particle& a,double re,double res[6]);
	void bucPosition(Particle& a,double timestep);
};


	