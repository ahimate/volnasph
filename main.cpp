#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream> 
#include <map>
#include <stdio.h>
#include <stdlib.h>

#include "Particle.h"
#include "Bucket.h"
#include "Vector3D.h"
#include "constants.h"

using namespace std;

double grdXmin, grdXmax, grdYmin, grdYmax, grdZmin, grdZmax;
int parXcount, parYcount, parZcount; //a number of particles under x,y,z directions
double fluXStep, fluYStep, fluZStep; //steps under x,y,z directions
double wallXStep, wallYStep, wallZStep; //steps under x,y,z directions
int wallXcount, wallYcount;
double wallXmin, wallXmax, wallYmin, wallYmax, wallZmin, wallZmax;
char fluidSurfaceFileName[80], wallSurfaceFileName[80], outDataName[80];

int gridXstep, gridYstep;
double** wall; //height of wall/slope particles
double** flu;
//int filename_count = 5;
double timeStep; // = 0.1;
int steptimeout;
int stepCount; // = 100;

//Bucket
double bucketRadius; // radius of bucket
map<int, Bucket> mapa;
map<int, Bucket> mapa_new;
int particlesCount;// a number of particles at all
int kolvoParInRadX;// = 40; //a number of particles in support radius

double sigma; // surface tension coefficient
double wallRoo = 1000;//initial density of wall particles

int getHash(double x, double y, double z, ofstream &out)//c нуля,а количество с 1
{
	//	if ((x < wallXmin) || (x > (wallXmax + bucketRadius)) || (y < wallYmin)
	//			|| (y > (wallYmax + bucketRadius))) {
	//		cout << "Error getHash(): point is outside the area!" << x << ", " << y
	//				<< ", " << z << " )" << endl;
	//		out << "Error getHash(): point is outside the area!" << x << ", " << y
	//				<< ", " << z << " )" << endl;
	//		return 0;
	//	} else {
	int xx = (int) ((x - wallXmin) / (2 * bucketRadius));
	int yy = (int) ((y - wallYmin) / (2 * bucketRadius));
	int zz = (int) ((z - wallZmin) / (2 * bucketRadius));
	if ((xx < 9999) && (yy < 9999) && (zz < 9999)) {
		return xx * 100000000 + yy * 10000 + zz;
	} else {
		cout << "Error getHash(): bucket counter is too big" << endl;
		out << "Error getHash(): bucket counter is too big" << endl;
		return 0;
	}
	//	}
}
void reverseHash(int key, double array[3], ofstream &out) { //returns bucket coordinates! not particles coordinates
	if (key < 0) {
		cout << "Error reverseHash: point is outside the area! " << key << endl;
		out << "Error reverseHash: point is outside the area! " << key << endl;
	} else {
		array[0] = (int) (key / 100000000);
		array[1] = (int) ((key - array[0] * 100000000) / 1000);
		array[2] = (int) (key - array[0] * 100000000 - array[1] * 10000);
	}
}

void outputParticles(map<int, Bucket>& mapama, const char* outname) {
	ofstream out(outname);
	out << "x" << "\t";
	out << "y" << "\t";
	out << "z" << "\t";
	out << "U" << "\t";
	out << "V" << "\t";
	out << "W" << "\t";
	out << "Ro" << "\t";
	out << "P" << "\t";
	out << endl;
	for (map<int, Bucket>::iterator it = mapama.begin(); it != mapama.end(); it++) {
		(*it).second.BucketParticlesInFile(out);
	}
	out.close();
	cout << outname << " created" << endl;
}
void outputFluidParticles(map<int, Bucket>& mapama, const char* outname) {
	ofstream out(outname);
	out << "x" << "\t";
	out << "y" << "\t";
	out << "z" << "\t";
	out << "U" << "\t";
	out << "V" << "\t";
	out << "W" << "\t";
	out << "Ro" << "\t";
	out << "P" << "\t";
	out << endl;
	for (map<int, Bucket>::iterator it = mapama.begin(); it != mapama.end(); it++) {
		if ((*it).second.status == THERE_IS_FLUID) {
			(*it).second.FluidParticlesInFile(out);
		}
	}
	out.close();
	cout << outname << " created" << endl;
}

void inputData(char* fluidSurfaceFileName, char* wallSurfaceFileName,
		const char* outnm) // input initial parameters from file
{
	cout << endl << "Input data... " << fluidSurfaceFileName << " , "
			<< wallSurfaceFileName << endl;
	ofstream out(outnm, ios::trunc);
	ifstream fluidInput(fluidSurfaceFileName);
	if (fluidInput.fail()) {
		cout << "Error reading " << fluidSurfaceFileName << endl;
		out << "Error reading " << fluidSurfaceFileName << endl;
	}
	ifstream wallInput(wallSurfaceFileName);
	if (wallInput.fail()) {
		cout << "Error reading " << wallSurfaceFileName << endl;
		out << "Error reading " << wallSurfaceFileName << endl;
	}
	string tmpString;
	wallInput >> tmpString >> wallXcount >> wallYcount >> wallXmin >> wallXmax
			>> wallYmin >> wallYmax >> wallZmin >> wallZmax;
	fluidInput >> tmpString >> parXcount >> parYcount >> grdXmin >> grdXmax
			>> grdYmin >> grdYmax >> grdZmin >> grdZmax;

	cout << "wallXmin" << "\t" << "wallXmax" << "\t" << "wallYmin" << "\t"
			<< "wallYmax" << "\t" << "wallZmin" << "\t" << "wallZmax" << endl;
	cout << wallXmin << "\t" << wallXmax << "\t" << wallYmin << "\t"
			<< wallYmax << "\t" << wallZmin << "\t" << wallZmax << endl;
	out << "wallXmin" << " " << "wallXmax" << " " << "wallYmin" << " "
			<< "wallYmax" << " " << "wallZmin" << " " << "wallZmax" << endl;
	out << wallXmin << "\t" << wallXmax << "\t" << wallYmin << "\t" << wallYmax
			<< "\t" << wallZmin << "\t" << wallZmax << endl;

	fluXStep = (grdXmax - grdXmin) / (parXcount - 1);
	fluYStep = (grdYmax - grdYmin) / (parYcount - 1);
	wallXStep = (wallXmax - wallXmin) / (wallXcount - 1);
	wallYStep = (wallYmax - wallYmin) / (wallYcount - 1);

	wall = new double*[wallXcount];
	for (int i = 0; i < wallXcount; i++)
		wall[i] = new double[wallYcount];
	flu = new double*[parXcount];
	for (int i = 0; i < parXcount; i++)
		flu[i] = new double[parYcount];

	double elementsCount = 0; //needed to count the volume

	for (int j = 0; j < wallYcount; j++) {
		for (int i = 0; i < wallXcount; i++) {
			wallInput >> wall[i][j];
			elementsCount += 1;
		}
	}
	double Vwall = elementsCount * wallXStep * wallYStep;
	elementsCount = 0;
	for (int j = 0; j < parYcount; j++) {
		for (int i = 0; i < parXcount; i++) {
			fluidInput >> flu[i][j];
			elementsCount += flu[i][j];
		}
	}

	double Vfluid = elementsCount * fluXStep * fluYStep;
	if (Particle::fluNu > 10) {
		Particle::fluM = 0.04;
	} else {
		Particle::fluM = 0.02;
	}

	particlesCount = (int) ((Particle::fluROO * Vfluid) / Particle::fluM);
	bucketRadius = pow((3 * Vfluid * kolvoParInRadX)
			/ (4 * pi * particlesCount), 1.0 / 3.0);
	Vwall = Vwall * 2 * bucketRadius;
	Particle::wallM = (Particle::wallROO * Vwall) / particlesCount;
	double parStep = pow(Vfluid / particlesCount, 0.333333333);

	out << "kolvoParN flu   " << particlesCount << endl;
	out << "parStep   " << parStep << endl;
	out << "bucRadius   " << bucketRadius << endl;
	out << "fluM " << Particle::fluM << endl;
	out << "wallM  " << Particle::wallM << endl;

	cout << "kolvoParN flu   " << particlesCount << endl;
	cout << "parStep   " << parStep << endl;
	cout << "bucRadius   " << bucketRadius << endl;
	cout << "fluM " << Particle::fluM << endl;
	cout << "wallM  " << Particle::wallM << endl;

	if (bucketRadius <= parStep) {
		cout << "Error kolvoParN and input data" << endl;
		out << "Error kolvoParN and input data" << endl;
	}
	//breaking slope into particles
	for (double tX = wallXmin; tX <= wallXmax; tX += parStep) {
		for (double tY = wallYmin; tY <= wallYmax; tY += parStep) {
			double tZ = wall[(int) ((tX - wallXmin) / wallXStep)][(int) ((tY
					- wallYmin) / wallYStep)];
			elementsCount = tZ - bucketRadius;

			double tmpP;
			if ((tX < grdXmin) || (tX > grdXmax) || (tY < grdYmin) || (tY
					> grdYmax)) {
				tmpP = 0;
			} else {
				tmpP = Particle::fluROO * g * flu[(int) ((tX - grdXmin)
						/ fluXStep)][(int) ((tX - grdYmin) / fluYStep)];
			}
			for (double k = (tZ - parStep); k >= elementsCount; k -= parStep) {
				Particle a(tX, tY, k);
				a.cP = tmpP;

				if (mapa.find(getHash(tX, tY, k, out)) == mapa.end()) {
					Bucket b(int((tX - wallXmin) / (2 * bucketRadius)) * (2
							* bucketRadius) + wallXmin, int((tY - wallYmin)
							/ (2 * bucketRadius)) * (2 * bucketRadius)
							+ wallYmin,
							int((k - wallZmin) / (2 * bucketRadius)) * (2
									* bucketRadius) + wallZmin);
					mapa[getHash(tX, tY, k, out)] = b;
				}
				mapa[getHash(tX, tY, k, out)].particlesVector.push_back(a);
				//mapa[getHash(tX, tY, k)].status = THERE_NO_FLUID; //???
			}
			if ((tY < (wallYmin + bucketRadius)) || (tY > (wallYmax
					- bucketRadius))) {
				for (double k = tZ; k <= (tZ + 2 * wallZmax + 2 * grdZmax); k
						+= parStep) {
					Particle a(tX, tY, k);
					a.cP = 0;
					if (mapa.find(getHash(tX, tY, k, out)) == mapa.end()) {
						Bucket b(int((tX - wallXmin) / (2 * bucketRadius)) * (2
								* bucketRadius) + wallXmin, int((tY - wallYmin)
								/ (2 * bucketRadius)) * (2 * bucketRadius)
								+ wallYmin, int((k - wallZmin) / (2
								* bucketRadius)) * (2 * bucketRadius)
								+ wallZmin);
						mapa[getHash(tX, tY, k, out)] = b;
					}
					mapa[getHash(tX, tY, k, out)].particlesVector.push_back(a);
					//mapa[getHash(tX, tY, k)].status = THERE_NO_FLUID; //???
				}
			}

			if (((tX < (wallXmin + bucketRadius)) || (tX > (wallXmax
					- bucketRadius))) && (tY > (wallYmin + bucketRadius))
					&& (tY < (wallYmax - bucketRadius))) {
				for (double k = tZ; k <= (tZ + 2 * wallZmax + 2 * grdZmax); k
						+= parStep) {
					Particle a(tX, tY, k);
					a.cP = 0;
					if (mapa.find(getHash(tX, tY, k, out)) == mapa.end()) {
						Bucket b(int((tX - wallXmin) / (2 * bucketRadius)) * (2
								* bucketRadius) + wallXmin, int((tY - wallYmin)
								/ (2 * bucketRadius)) * (2 * bucketRadius)
								+ wallYmin, int((k - wallZmin) / (2
								* bucketRadius)) * (2 * bucketRadius)
								+ wallZmin);
						mapa[getHash(tX, tY, k, out)] = b;
					}
					mapa[getHash(tX, tY, k, out)].particlesVector.push_back(a);
					//mapa[getHash(tX, tY, k)].status = THERE_NO_FLUID; //???
				}
			}
		}
	}
	//breaking fluid into particles
	int tmp5 = 0;
	for (double tX = grdXmin; tX <= grdXmax; tX += parStep) {
		for (double tY = grdYmin; tY <= grdYmax; tY += parStep) {
			double tmp1 = wall[(int) ((tX - wallXmin) / wallXStep)][(int) ((tY
					- wallYmin) / wallYStep)];
			double tmp4 = flu[(int) ((tX - grdXmin) / fluXStep)][(int) ((tY
					- grdYmin) / fluYStep)];
			elementsCount = tmp1 + tmp4;
			for (double k = tmp1; k <= elementsCount; k += parStep) {
				if ((tY >= (wallYmin + bucketRadius)) && (tY <= (wallYmax
						- bucketRadius)) && (tX >= (wallXmin + bucketRadius))
						&& (tX <= (wallXmax - bucketRadius))) {
					Particle a(tX, tY, k, 0, 0, 0);
					tmp5++;
					a.cP = Particle::fluROO * g * (elementsCount - k);
					a.cU5 = a.cU;
					a.cV5 = a.cV;
					a.cW5 = a.cW;
					if (mapa.find(getHash(tX, tY, k, out)) == mapa.end()) {
						Bucket b(int((tX - wallXmin) / (2 * bucketRadius)) * (2
								* bucketRadius) + wallXmin, int((tY - wallYmin)
								/ (2 * bucketRadius)) * (2 * bucketRadius)
								+ wallYmin, int((k - wallZmin) / (2
								* bucketRadius)) * (2 * bucketRadius)
								+ wallZmin);
						mapa[getHash(tX, tY, k, out)] = b;
					}
					mapa[getHash(tX, tY, k, out)].particlesVector.push_back(a);
					mapa[getHash(tX, tY, k, out)].status = THERE_IS_FLUID;
				}
			}
		}
	}

	out << "created number of fluid particles " << tmp5 << endl;
	Particle::fluM = (Particle::fluROO * Vfluid) / tmp5;
	out << "real fluid particle mass " << Particle::fluM << endl;
	cout << "Finish input data " << endl;
	fluidInput.close();
	wallInput.close();
	out.close();
}
void inputParams(char* input) {
	ifstream inp(input);
	if (inp.fail()) {
		cout << "Error reading " << input << "'." << endl;
	}

	inp >> timeStep >> stepCount >> steptimeout;
	inp >> Particle::fluNu >> sigma >> Particle::fluK;
	inp >> Particle::fluROO >> Particle::wallROO;
	inp >> kolvoParInRadX;
	inp >> fluidSurfaceFileName >> wallSurfaceFileName >> outDataName;
	inputData(fluidSurfaceFileName, wallSurfaceFileName, outDataName);
}

//бежит по мапа_нью и в пустой мапа добавляет ячейки и частицы
void updateBucket(int time, ofstream &out) {
	mapa_new.clear();
	mapa_new = mapa;
	mapa.clear();
	int count = 0;
	int tmpcount = 0;
	for (map<int, Bucket>::iterator it = mapa_new.begin(); it != mapa_new.end(); it++) {
		int part_size = (*it).second.particlesVector.size();
		for (int i = 0; i < part_size; i++) {
			double tmpx = (*it).second.particlesVector[i].x;
			double tmpy = (*it).second.particlesVector[i].y;
			double tmpz = (*it).second.particlesVector[i].z;

			if ((*it).second.particlesVector[i].type == PARTICLE_TYPE_FLUID) {
				//			if ((tmpx < wallXmin) || (tmpx > wallXmax) || (tmpy < wallYmin)
				//					|| (tmpy > wallYmax) || (tmpz < (wallZmin - bucketRadius))
				//					|| (tmpz > (wallZmax + grdZmax + 2 * bucketRadius))) {
				if ((tmpx < wallXmin) || (tmpx > wallXmax) || (tmpy < wallYmin)
						|| (tmpy > wallYmax) || (tmpz < (wallZmin
						- bucketRadius))) {
					//точку вне области не считаем

					count++;

				} else {
					if (mapa.find(getHash(tmpx, tmpy, tmpz, out)) == mapa.end()) {
						Bucket a(int((tmpx - grdXmin) / (2 * bucketRadius))
								* (2 * bucketRadius) + grdXmin, int((tmpy
								- grdYmin) / (2 * bucketRadius)) * (2
								* bucketRadius) + grdYmin, int((tmpz - grdZmin)
								/ (2 * bucketRadius)) * (2 * bucketRadius)
								+ grdZmin);
						mapa[getHash(tmpx, tmpy, tmpz, out)] = a;
					}
					mapa[getHash(tmpx, tmpy, tmpz, out)].particlesVector.push_back(
							(*it).second.particlesVector[i]);
					if ((*it).second.particlesVector[i].type
							== PARTICLE_TYPE_FLUID) {
						tmpcount++;
						mapa[getHash(tmpx, tmpy, tmpz, out)].status
								= THERE_IS_FLUID;
					}
				}
			} else {
				if (mapa.find(getHash(tmpx, tmpy, tmpz, out)) == mapa.end()) {
					Bucket a(int((tmpx - grdXmin) / (2 * bucketRadius)) * (2
							* bucketRadius) + grdXmin, int((tmpy - grdYmin)
							/ (2 * bucketRadius)) * (2 * bucketRadius)
							+ grdYmin, int((tmpz - grdZmin)
							/ (2 * bucketRadius)) * (2 * bucketRadius)
							+ grdZmin);
					mapa[getHash(tmpx, tmpy, tmpz, out)] = a;
				}
				mapa[getHash(tmpx, tmpy, tmpz, out)].particlesVector.push_back(
						(*it).second.particlesVector[i]);
				if ((*it).second.particlesVector[i].type == PARTICLE_TYPE_FLUID) {
					tmpcount++;
					mapa[getHash(tmpx, tmpy, tmpz, out)].status
							= THERE_IS_FLUID;
				}
			}
		}
	}
	cout << "T: " << time << "  particles: " << tmpcount << " free: " << count
			<< endl;
	out << "T: " << time << "  particles: " << tmpcount << " free: " << count
			<< endl;
	mapa_new.clear();
}
void BucketStep(ofstream &out) {
	//first calculating density and pressure
	int countZeroDensity = 0;
	for (map<int, Bucket>::iterator it = mapa.begin(); it != mapa.end(); it++) {
		int part_size = (*it).second.particlesVector.size();
		for (int i = 0; i < part_size; i++) {
			if ((*it).second.particlesVector[i].type == PARTICLE_TYPE_FLUID) {

				double tmp1 = 0;

				for (double x1 = -1; x1 < 2; x1++) {
					for (double x2 = -1; x2 < 2; x2++) {
						for (double x3 = -1; x3 < 2; x3++) {
							int forHash = (*it).first + getHash(x1 * 2
									* bucketRadius + wallXmin, x2 * 2
									* bucketRadius + wallYmin, x3 * 2
									* bucketRadius + wallZmin, out);
							if (mapa.find(forHash) == mapa.end()) {
								//do nothing for the bucket does not exist
							} else {
								//calculate density
								tmp1 += mapa[forHash].bucDen(
										(*it).second.particlesVector[i],
										bucketRadius);
							}
						}
					}
				}

				(*it).second.particlesVector[i].cRo = tmp1;

				if (tmp1 <= epsilon) {
					//if ((*it).second.particlesVector[i].type == PARTICLE_TYPE_FLUID) {
					countZeroDensity++;
					//}
					//Particle tmpPart = (*it).second.particlesVector[i];
					//cout << "Ro==0: " << (*it).second.particlesVector[i].cRo << " "
					//		<< (*it).second.particlesVector[i].type << " "
					//		<< bucketRadius << " " << tmpPart.x << " " << tmpPart.y
					//		<< " " << tmpPart.z << " " << " " << endl;
					//	 out << "Ro==0: " << (*it).second.particlesVector[i].cRo << " "
					//	 << (*it).second.particlesVector[i].type << " " << (*it).first << endl;
				}
			}
			//need to calculate density for wall particles or not

			//calculating pressure for particles
			if ((*it).second.particlesVector[i].type == PARTICLE_TYPE_FLUID) {
				(*it).second.particlesVector[i].cP
						= (*it).second.particlesVector[i].fluK
								* ((*it).second.particlesVector[i].cRo
										- (*it).second.particlesVector[i].fluROO);
				//			} else {
				//				(*it).second.particlesVector[i].cP
				//						= (*it).second.particlesVector[i].fluK
				//								* ((*it).second.particlesVector[i].cRo
				//										- (*it).second.particlesVector[i].wallROO);
			}
			//need to calculate new pressure for wall particles
			//now it is consider pressure of wall particles does not change
			//	}
		}
	}
	cout << "Zero density particles number = " << countZeroDensity << endl;

	cout << "FINISHING density and press" << endl;
	out << "FINISHING density and press" << endl;

	//	mapa_new.clear(); //first delete all elements of mapa_new
	//	mapa_new = mapa; // then use it for new values calculation

	for (map<int, Bucket>::iterator it = mapa.begin(); it != mapa.end(); it++) {
		int part_size = (*it).second.particlesVector.size();
		for (int i = 0; i < part_size; i++) {
			if ((*it).second.particlesVector[i].type == PARTICLE_TYPE_FLUID) {
				double forse[6] = { 0, 0, 0, 0, 0, 0 };
				double nsurf[3] = { 0, 0, 0 };
				double fsurf = 0;
				for (int x1 = -1; x1 < 2; x1++) {
					for (int x2 = -1; x2 < 2; x2++) {
						for (int x3 = -1; x3 < 2; x3++) {
							int forHash = (*it).first + getHash(x1 * 2
									* bucketRadius + wallXmin, x2 * 2
									* bucketRadius + wallYmin, x3 * 2
									* bucketRadius + wallZmin, out);
							if (mapa.find(forHash) == mapa.end()) {
							} else {
								double f[6] = { 0, 0, 0, 0, 0, 0 };
								double fs[3] = { 0, 0, 0 };
								mapa[forHash].bucForse(
										(*it).second.particlesVector[i],
										bucketRadius, f);
								mapa[forHash].bucFnorm(
										(*it).second.particlesVector[i],
										bucketRadius, fs);

								for (int k = 0; k < 6; k++) {
									forse[k] += f[k];
								}
								for (int k = 0; k < 3; k++) {
									nsurf[k] += fs[k];
								}
							}
						}

					}
				}
				double ni2 = pow(nsurf[0], 2) + pow(nsurf[1], 2) + pow(
						nsurf[2], 2);
				double surL2 = sqrt((Particle::fluROO) / kolvoParInRadX);
				if (ni2 >= surL2) {
					for (int x1 = -1; x1 < 2; x1++) {
						for (int x2 = -1; x2 < 2; x2++) {
							for (int x3 = -1; x3 < 2; x3++) {
								int forHash = (*it).first + getHash(x1 * 2
										* bucketRadius + wallXmin, x2 * 2
										* bucketRadius + wallYmin, x3 * 2
										* bucketRadius + wallZmin, out);
								if (mapa.find(forHash) != mapa.end()) {
									//double fs = 0;
									fsurf += mapa[forHash].bucFsurf(
											(*it).second.particlesVector[i],
											bucketRadius);
								}
							}
						}
					}
					fsurf /= sqrt(ni2);
				}
				double acceleration[3] = { 0, 0, 0 };
				double gravity[3] = { 0, 0, g };
				for (int k = 0; k < 3; k++) {
					acceleration[k] = forse[k]
							+ ((*it).second.particlesVector[i].fluNu * forse[k
									+ 3] - sigma * nsurf[k] * fsurf
									- gravity[k])
									/ (*it).second.particlesVector[i].cRo;
				}
				//save velocity/2 on (t-0.5dt)
				(*it).second.particlesVector[i].cU = 0.5
						* (*it).second.particlesVector[i].cU5;
				(*it).second.particlesVector[i].cV = 0.5
						* (*it).second.particlesVector[i].cV5;
				(*it).second.particlesVector[i].cW = 0.5
						* (*it).second.particlesVector[i].cW5;

				//calculate new velocity on (t+0.5dt)
				(*it).second.particlesVector[i].cU5 += timeStep
						* acceleration[0];
				(*it).second.particlesVector[i].cV5 += timeStep
						* acceleration[1];
				(*it).second.particlesVector[i].cW5 += timeStep
						* acceleration[2];

				//calculate new coordinates on (t+dt)
				(*it).second.particlesVector[i].x += timeStep
						* (*it).second.particlesVector[i].cU5;
				(*it).second.particlesVector[i].y += timeStep
						* (*it).second.particlesVector[i].cV5;
				(*it).second.particlesVector[i].z += timeStep
						* (*it).second.particlesVector[i].cW5;

				//calculate new velocity on (t+dt)
				(*it).second.particlesVector[i].cU += 0.5
						* (*it).second.particlesVector[i].cU5;
				(*it).second.particlesVector[i].cV += 0.5
						* (*it).second.particlesVector[i].cV5;
				(*it).second.particlesVector[i].cW += 0.5
						* (*it).second.particlesVector[i].cW5;
			}
		}
	}
	//	mapa.clear();
	//	mapa = mapa_new;
	cout << "FINISH forse and velosity" << endl;
}

int main(int argc, char* argv[]) {

	if (argc != 2) {
		cout << "Error argv" << endl;
		cout << "Usage: volnasph inputParametersFileName" << endl;
		return 0;
	}
	inputParams(argv[1]);
	ofstream out(outDataName, ios::app);
	outputFluidParticles(mapa, "Flu_000.txt");
	outputParticles(mapa, "All_000.txt");
	for (int t = 1; t <= stepCount; t++) {
		cout << "t   " << t << endl;
		out << "t   " << t << endl;
		BucketStep(out);
		if ((t % steptimeout) == 0) {
			char outname[30];
			sprintf(outname, "\%05d_flu.txt", t);
			outputFluidParticles(mapa, outname);
			sprintf(outname, "\%05d_all.txt", t);
			outputParticles(mapa, outname);
		}
		updateBucket(t, out);
		//		char tmpch;
		//		cin >> tmpch;
	}
	cout << "Finished" << endl;
	out << "Finished" << endl;
	out.close();
	return 0;
}

