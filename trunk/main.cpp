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

int getHash(double x, double y, double z)//c нуля,а количество с 1
{
	if ((x < wallXmin) || (x > (wallXmax + bucketRadius)) || (y < wallYmin)
			|| (y > (wallYmax + bucketRadius))) {
		cout << "Error Hash: point is outside the area!" << x << ", " << y
				<< ", " << z << " )" << endl;
		return 0;
	} else {
		int xx = (int) ((x - wallXmin) / (2 * bucketRadius));
		int yy = (int) ((y - wallYmin) / (2 * bucketRadius));
		int zz = (int) ((z - wallZmin) / (2 * bucketRadius));
		if ((xx < 9999) && (yy < 9999) && (zz < 9999)) {
			return xx * 100000000 + yy * 10000 + zz;
		} else {
			cout << "Error Hash: bucket counter is too big" << endl;
			return 0;
		}
	}
}

void outputData(map<int, Bucket>& mapama, const char* outname) {
	ofstream out(outname);
	for (map<int, Bucket>::iterator it = mapama.begin(); it != mapama.end(); it++) {
		(*it).second.BucketInFile(out);
	}
	out.close();
	cout << outname << " created" << endl;
}

void outputParticles(map<int, Bucket>& mapama, const char* outname) {
	ofstream out(outname);
	for (map<int, Bucket>::iterator it = mapama.begin(); it != mapama.end(); it++) {
		(*it).second.BucketParticlesInFile(out);
	}
	out.close();
	cout << outname << " created" << endl;
}

void outputFluidData(map<int, Bucket>& mapama, const char* outname) {
	ofstream out(outname);
	for (map<int, Bucket>::iterator it = mapama.begin(); it != mapama.end(); it++) {
		if ((*it).second.status == THERE_IS_FLUID) {
			(*it).second.FluInFile(out);
		}
	}
	out.close();
	cout << outname << " created" << endl;
}
void outputFluidParticles(map<int, Bucket>& mapama, const char* outname) {
	ofstream out(outname);
	for (map<int, Bucket>::iterator it = mapama.begin(); it != mapama.end(); it++) {
		if ((*it).second.status == THERE_IS_FLUID) {
			(*it).second.FluidParticlesInFile(out);
		}
	}
	out.close();
	cout << outname << " created" << endl;
}

void inputData(char* fluidSurfaceFileName,char* wallSurfaceFileName, const char* outnm) // input initial parameters from file
{
	cout << endl << "Input data... " << fluidSurfaceFileName << " , "
			<< wallSurfaceFileName << endl;
	ofstream out(outnm);
	ifstream fluidInput(fluidSurfaceFileName);
	if (fluidInput.fail()) {
		cout << "Error reading " << fluidSurfaceFileName << endl;
	}
	ifstream wallInput(wallSurfaceFileName);
	if (wallInput.fail()) {
		cout << "Error reading " << wallSurfaceFileName << endl;
	}
	string tmpString;
	wallInput >> tmpString >> wallXcount >> wallYcount >> wallXmin >> wallXmax
			>> wallYmin >> wallYmax >> wallZmin >> wallZmax;
	fluidInput >> tmpString >> parXcount >> parYcount >> grdXmin >> grdXmax
			>> grdYmin >> grdYmax >> grdZmin >> grdZmax;

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

	for (int j = 0; j < wallYcount; j++) {
		for (int i = 0; i < wallXcount; i++) {
			wallInput >> wall[i][j];
		}
	}
	double elementsCount = 0; //needed to count the volume
	for (int j = 0; j < parYcount; j++) {
		for (int i = 0; i < parXcount; i++) {
			fluidInput >> flu[i][j];
			elementsCount += flu[i][j];
		}
	}

	double Vall = elementsCount * fluXStep * fluYStep;
	if (Particle::fluNu > 10) {
		Particle::fluM = 0.04;
	} else {
		Particle::fluM = 0.02;
	}

	particlesCount = (int) ((Particle::fluROO * Vall) / Particle::fluM);
	Particle::wallM = (Particle::wallROO * Vall) / particlesCount;
	double Vpar = Vall / particlesCount;
	double parStep = pow(Vpar, 0.333333333);

	elementsCount = (3 * Vall * kolvoParInRadX) / (4 * pi * particlesCount);
	bucketRadius = pow(elementsCount, 0.33333333333);
	out << "kolvoParN flu   " << particlesCount << endl;
	out << "parStep   " << parStep << endl;
	out << "bucRadius   " << bucketRadius << endl;
	out << "fluM " << Particle::fluM << endl;
	out << "wallM  " << Particle::wallM << endl;

	if (bucketRadius <= parStep) {
		cout << "Error kolvoParN and input data" << endl;
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
				tmpP = Particle::fluROO * 9.82 * flu[(int) ((tX - grdXmin)
						/ fluXStep)][(int) ((tX - grdYmin) / fluYStep)];
			}
			for (double k = (tZ - parStep); k >= elementsCount; k -= parStep) {
				Particle a(tX, tY, k);
				a.cP = tmpP;

				if (mapa.find(getHash(tX, tY, k)) == mapa.end()) {
					Bucket b(int((tX - wallXmin) / (2 * bucketRadius)) * (2
							* bucketRadius) + wallXmin, int((tY - wallYmin)
							/ (2 * bucketRadius)) * (2 * bucketRadius)
							+ wallYmin,
							int((k - wallZmin) / (2 * bucketRadius)) * (2
									* bucketRadius) + wallZmin);
					mapa[getHash(tX, tY, k)] = b;
				}
				mapa[getHash(tX, tY, k)].part.push_back(a);
				mapa[getHash(tX, tY, k)].status = THERE_NO_FLUID; //???
			}
			if ((tY < (wallYmin + bucketRadius)) || (tY > (wallYmax
					- bucketRadius))) {
				for (double k = tZ; k <= (tZ + 2 * wallZmax + 2 * grdZmax); k
						+= parStep) {
					Particle a(tX, tY, k);
					a.cP = 0;
					if (mapa.find(getHash(tX, tY, k)) == mapa.end()) {
						Bucket b(int((tX - wallXmin) / (2 * bucketRadius)) * (2
								* bucketRadius) + wallXmin, int((tY - wallYmin)
								/ (2 * bucketRadius)) * (2 * bucketRadius)
								+ wallYmin, int((k - wallZmin) / (2
								* bucketRadius)) * (2 * bucketRadius)
								+ wallZmin);
						mapa[getHash(tX, tY, k)] = b;
					}
					mapa[getHash(tX, tY, k)].part.push_back(a);
					mapa[getHash(tX, tY, k)].status = THERE_NO_FLUID; //???
				}
			}

			if (((tX < (wallXmin + bucketRadius)) || (tX > (wallXmax
					- bucketRadius))) && (tY > (wallYmin + bucketRadius))
					&& (tY < (wallYmax - bucketRadius))) {
				for (double k = tZ; k <= (tZ + 2 * wallZmax + 2 * grdZmax); k
						+= parStep) {
					Particle a(tX, tY, k);
					a.cP = 0;
					if (mapa.find(getHash(tX, tY, k)) == mapa.end()) {
						Bucket b(int((tX - wallXmin) / (2 * bucketRadius)) * (2
								* bucketRadius) + wallXmin, int((tY - wallYmin)
								/ (2 * bucketRadius)) * (2 * bucketRadius)
								+ wallYmin, int((k - wallZmin) / (2
								* bucketRadius)) * (2 * bucketRadius)
								+ wallZmin);
						mapa[getHash(tX, tY, k)] = b;
					}
					mapa[getHash(tX, tY, k)].part.push_back(a);
					mapa[getHash(tX, tY, k)].status = THERE_NO_FLUID; //???
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
					a.cP = Particle::fluROO * 9.82 * (elementsCount - k);
					a.cU5 = a.cU;
					a.cV5 = a.cV;
					a.cW5 = a.cW;
					if (mapa.find(getHash(tX, tY, k)) == mapa.end()) {
						Bucket b(int((tX - wallXmin) / (2 * bucketRadius)) * (2
								* bucketRadius) + wallXmin, int((tY - wallYmin)
								/ (2 * bucketRadius)) * (2 * bucketRadius)
								+ wallYmin, int((k - wallZmin) / (2
								* bucketRadius)) * (2 * bucketRadius)
								+ wallZmin);
						mapa[getHash(tX, tY, k)] = b;
					}
					mapa[getHash(tX, tY, k)].part.push_back(a);
					mapa[getHash(tX, tY, k)].status = THERE_IS_FLUID;
				}
			}
		}
	}

	out << "created number of fluid particles " << tmp5 << endl;
	Particle::fluM = (Particle::fluROO * Vall) / tmp5;
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
//	inp >> epsilon;
	inp >> Particle::fluNu >> sigma >> Particle::fluK;
	inp >> Particle::fluROO >> Particle::wallROO;
	inp >> kolvoParInRadX;
	char fluidSurfaceFileName[80], wallSurfaceFileName[80], outDataName[80];
	inp >> fluidSurfaceFileName >> wallSurfaceFileName >> outDataName;
	inputData(fluidSurfaceFileName, wallSurfaceFileName, outDataName);
}

//бежит по мапа_нью и в пустой мапа добавляет ячейки и частицы
void updateBucket(int time) {
	mapa.clear();
	int count = 0;
	int tmpcount = 0;
	for (map<int, Bucket>::iterator it = mapa_new.begin(); it != mapa_new.end(); it++) {
		int part_size = (*it).second.part.size();
		for (int i = 0; i < part_size; i++) {
			double tmpx = (*it).second.part[i].x;
			double tmpy = (*it).second.part[i].y;
			double tmpz = (*it).second.part[i].z;

			if ((tmpx < wallXmin) || (tmpx > wallXmax) || (tmpy < wallYmin)
					|| (tmpy > wallYmax) || (tmpz < (wallZmin - bucketRadius))
					|| (tmpz > (wallZmax + grdZmax + 2 * bucketRadius))) {
				//точку вне области не считаем
				if ((*it).second.part[i].type == PARTICLE_TYPE_FLUID) {
					count++;
				}
			} else {
				if (mapa.find(getHash(tmpx, tmpy, tmpz)) == mapa.end()) {
					Bucket a(int((tmpx - grdXmin) / (2 * bucketRadius)) * (2
							* bucketRadius) + grdXmin, int((tmpy - grdYmin)
							/ (2 * bucketRadius)) * (2 * bucketRadius)
							+ grdYmin, int((tmpz - grdZmin)
							/ (2 * bucketRadius)) * (2 * bucketRadius)
							+ grdZmin);
					mapa[getHash(tmpx, tmpy, tmpz)] = a;
				}
				mapa[getHash(tmpx, tmpy, tmpz)].part.push_back(
						(*it).second.part[i]);
				if ((*it).second.part[i].type == PARTICLE_TYPE_FLUID) {
					tmpcount++;
					mapa[getHash(tmpx, tmpy, tmpz)].status = THERE_IS_FLUID;
				}
			}
		}
	}
	cout << "T: " << time << "  part: " << tmpcount << " free: " << count
			<< endl;
	mapa_new.clear();
}

void BucketStep() {
	for (map<int, Bucket>::iterator it = mapa.begin(); it != mapa.end(); it++) {
		int part_size = (*it).second.part.size();
		for (int i = 0; i < part_size; i++) {
			if ((*it).second.part[i].type == 1) {
				double tmp1 = 0;
				for (int x1 = -1; x1 < 2; x1++) {
					for (int x2 = -1; x2 < 2; x2++) {
						for (int x3 = -1; x3 < 2; x3++) {
							if (mapa.find((*it).first + x1 * 100000000 + x2
									* 10000 + x3) == mapa.end()) {

							} else {
								tmp1 += mapa[((*it).first + x1 * 100000000 + x2
										* 10000 + x3)].bucDen(
										(*it).second.part[i], bucketRadius);
							}
						}
					}
				}
				(*it).second.part[i].cRo = tmp1;
				if (tmp1 <= epsilon) {
					cout << "Ro==0: " << (*it).second.part[i].cRo;
					int tmp11;
					cin >> tmp11;
				}
				if ((*it).second.part[i].type == PARTICLE_TYPE_FLUID) {
					(*it).second.part[i].cP = (*it).second.part[i].fluK
							* ((*it).second.part[i].cRo
									- (*it).second.part[i].fluROO);
				}
			}
		}
	}
	cout << "FINISH density and press" << endl;
	mapa_new = mapa;

	for (map<int, Bucket>::iterator it = mapa.begin(); it != mapa.end(); it++) {
		int part_size = (*it).second.part.size();
		for (int i = 0; i < part_size; i++) {
			if ((*it).second.part[i].type == 1) {
				double forse[6] = { 0, 0, 0, 0, 0, 0 };
				double nsurf[3] = { 0, 0, 0 };
				double fsurf = 0;
				for (int x1 = -1; x1 < 2; x1++) {
					for (int x2 = -1; x2 < 2; x2++) {
						for (int x3 = -1; x3 < 2; x3++) {
							if (mapa.find((*it).first + x1 * 100000000 + x2
									* 10000 + x3) == mapa.end()) {
							} else {
								double f[6] = { 0, 0, 0, 0, 0, 0 };
								double fs[3] = { 0, 0, 0 };
								mapa[((*it).first + x1 * 100000000 + x2 * 10000
										+ x3)].bucForse((*it).second.part[i],
										bucketRadius, f);
								mapa[((*it).first + x1 * 100000000 + x2 * 10000
										+ x3)].bucFnorm((*it).second.part[i],
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

				double surL2 = (Particle::fluROO) / kolvoParInRadX;
				if (ni2 >= surL2) {
					for (int x1 = -1; x1 < 2; x1++) {
						for (int x2 = -1; x2 < 2; x2++) {
							for (int x3 = -1; x3 < 2; x3++) {
								if (mapa.find((*it).first + x1 * 100000000 + x2
										* 10000 + x3) == mapa.end()) {
								} else {
									double fs = 0;
									fsurf += mapa[((*it).first + x1 * 100000000
											+ x2 * 10000 + x3)].bucFsurf(
											(*it).second.part[i], bucketRadius);
								}
							}
						}
					}
				}
				double acceleration[3] = { 0, 0, 0 };
				for (int k = 0; k < 3; k++) {
					double tmp8 = 0;
					if (fsurf == 0) {
						tmp8 = (*it).second.part[i].cRo * forse[k]
								+ (*it).second.part[i].fluNu * forse[k + 3];
					} else {
						tmp8 = (*it).second.part[i].cRo * forse[k]
								+ (*it).second.part[i].fluNu * forse[k + 3]
								- sigma * nsurf[k] * fsurf / (sqrt(ni2));
					}
					if (k == 2) {
						tmp8 += (*it).second.part[i].cRo * (-9.82);
					}
					acceleration[k] = tmp8 / ((*it).second.part[i].cRo);
				}
				mapa_new[(*it).first].part[i].cU5 = (*it).second.part[i].cU5
						+ timeStep * acceleration[0];
				mapa_new[(*it).first].part[i].cV5 = (*it).second.part[i].cV5
						+ timeStep * acceleration[1];
				mapa_new[(*it).first].part[i].cW5 = (*it).second.part[i].cW5
						+ timeStep * acceleration[2];

				mapa_new[(*it).first].part[i].x = (*it).second.part[i].x
						+ timeStep * mapa_new[(*it).first].part[i].cU5;
				mapa_new[(*it).first].part[i].y = (*it).second.part[i].y
						+ timeStep * mapa_new[(*it).first].part[i].cV5;
				mapa_new[(*it).first].part[i].y = (*it).second.part[i].z
						+ timeStep * mapa_new[(*it).first].part[i].cW5;

				mapa_new[(*it).first].part[i].cU = 0.5
						* ((*it).second.part[i].cU5
								+ mapa_new[(*it).first].part[i].cU5);
				mapa_new[(*it).first].part[i].cV = 0.5
						* ((*it).second.part[i].cV5
								+ mapa_new[(*it).first].part[i].cV5);
				mapa_new[(*it).first].part[i].cW = 0.5
						* ((*it).second.part[i].cW5
								+ mapa_new[(*it).first].part[i].cW5);
			}
		}
	}
	mapa = mapa_new;
	cout << "FINISH forse and velosity" << endl;
}

int main(int argc, char* argv[]) {

	if (argc != 2) {
		cout << "Error argv" << endl;
		cout
				<< "Usage: volnasph inputParametersFileName"
				<< endl;
		return 0;
	}
	inputParams(argv[1]);

	outputFluidData(mapa, "Flu_000.txt");

	for(int t=1;t<=stepCount;t++)
	{
		cout<<"t   "<<t<<endl;
		BucketStep();
		if((t%steptimeout)==0)
		{
			char outname[30];
			sprintf(outname,"Movie\\%05d_flu.txt", t);
			outputFluidData(mapa,outname);
		}
		updateBucket(t);
	}
	cout<<"Finished"<<endl;
	return 0;
}

