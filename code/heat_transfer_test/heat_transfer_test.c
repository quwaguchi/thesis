#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double df_coef = 0.05;
double sliderSpeed = 1.0;
double FINISH_TIME;
double iceTemp = -1.0;
double sliderTemp = -1.0;

#define DIM                  2
#define PARTICLE_DISTANCE    0.00025 //0.025
#define DT                   0.001
#define OUTPUT_INTERVAL      1000

#define ARRAY_SIZE           20000
#define KINEMATIC_VISCOSITY  (1.0E-6)
#define FLUID_DENSITY        1000.0 
#define G_X  0.0      
#define G_Y  -9.8
#define G_Z  0.0      
#define RADIUS_FOR_NUMBER_DENSITY  (2.1*PARTICLE_DISTANCE) 
#define RADIUS_FOR_GRADIENT        (2.1*PARTICLE_DISTANCE) 
#define RADIUS_FOR_LAPLACIAN       (3.1*PARTICLE_DISTANCE) 
#define COLLISION_DISTANCE         (0.5*PARTICLE_DISTANCE)
#define THRESHOLD_RATIO_OF_NUMBER_DENSITY  0.97   
#define COEFFICIENT_OF_RESTITUTION 0.2
#define COMPRESSIBILITY (0.45E-9)
#define EPS             (0.01 * PARTICLE_DISTANCE)     
#define ON              1
#define OFF             0
#define RELAXATION_COEFFICIENT_FOR_PRESSURE 0.2
#define GHOST  -1
#define FLUID   0
#define SOLID_AND_FLUID   1
#define SOLID 2
#define AIR 3
// #define WALL    3
// #define DUMMY_WALL  4
#define SLIDER  5
#define GHOST_OR_DUMMY  -1
#define SURFACE_PARTICLE 1    
#define INNER_PARTICLE   0      
#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0 
#define DIRICHLET_BOUNDARY_IS_CONNECTED     1 
#define DIRICHLET_BOUNDARY_IS_CHECKED       2 

void initializeParticlePositionAndVelocity_for2dim( void );
void initializeParticlePositionAndVelocity_for3dim( void );
void calConstantParameter( void );
void calNZeroAndLambda( void );
double weight( double distance, double re );
void mainLoopOfSimulation( void );
void calGravity( void );
void calViscosity( void );
void moveParticle( void );
void collision( void );
void calPressure( void );
void calNumberDensity( void );
void setBoundaryCondition( void );
void setSourceTerm( void );
void setMatrix( void );
void exceptionalProcessingForBoundaryCondition( void );
void checkBoundaryCondition( void );
void increaseDiagonalTerm( void );
void solveSimultaneousEquationsByGaussEliminationMethod( void );
void removeNegativePressure( void );
void setMinimumPressure( void );
void calPressureGradient( void );
void moveParticleUsingPressureGradient( void );
void writeData_inProfFormat( void );
void writeData_inVtuFormat( void );
void writeData_inCsvFormat( void );


void calTemperature(void);
void moveSlider(void);

static double Acceleration[3*ARRAY_SIZE];
static int    ParticleType[ARRAY_SIZE];
static double Position[3*ARRAY_SIZE];
static double Velocity[3*ARRAY_SIZE];
static double Pressure[ARRAY_SIZE];
static double NumberDensity[ARRAY_SIZE];
static int    BoundaryCondition[ARRAY_SIZE];
static double SourceTerm[ARRAY_SIZE];
static int    FlagForCheckingBoundaryCondition[ARRAY_SIZE];
static double CoefficientMatrix[ARRAY_SIZE * ARRAY_SIZE];
static double MinimumPressure[ARRAY_SIZE];

int    FileNumber;
double Time;  
int    NumberOfParticles;
double Re_forNumberDensity,Re2_forNumberDensity; 
double Re_forGradient,     Re2_forGradient; 
double Re_forLaplacian,    Re2_forLaplacian; 
double N0_forNumberDensity;
double N0_forGradient;
double N0_forLaplacian;
double Lambda;
double collisionDistance,collisionDistance2;
double FluidDensity;
double Temperature[ARRAY_SIZE] = {};

int main( int argc, char** argv ) {

  printf("\n*** START PARTICLE-SIMULATION ***\n");
  if( DIM == 2 ){
    initializeParticlePositionAndVelocity_for2dim();
  }else{
    initializeParticlePositionAndVelocity_for3dim();
  }
  calConstantParameter();
  mainLoopOfSimulation();

  printf("*** END ***\n\n");
  return 0;
}


void initializeParticlePositionAndVelocity_for2dim( void ){
  int iX, iY;
  int nX, nY;
  double x, y, z;
  int i = 0;
  int flagOfParticleGeneration;

  nX = (int)(1.5/PARTICLE_DISTANCE);  
  nY = (int)(1.5/PARTICLE_DISTANCE);

  for(iX= 0; iX<nX; iX++){
      x = PARTICLE_DISTANCE * (double)(iX);
      flagOfParticleGeneration = OFF;
      
      if (x==0.0){
         ParticleType[i] = SOLID;
         flagOfParticleGeneration = ON;
         Temperature[i] = -15.0;
      }else if(x>0.0 && x<0.1){
         ParticleType[i] = FLUID;
         flagOfParticleGeneration = ON;
         Temperature[i] = 10.0;
      }

      if( flagOfParticleGeneration == ON){
        Position[i*3]=x; Position[i*3+1]=0.0; Position[i*3+2]=0.0;
        i++;
      }
  }
  NumberOfParticles = i;
}

void calConstantParameter( void ){
  Re_forNumberDensity  = RADIUS_FOR_NUMBER_DENSITY;  
  Re_forGradient       = RADIUS_FOR_GRADIENT;  
  Re_forLaplacian      = RADIUS_FOR_LAPLACIAN;  
  Re2_forNumberDensity = Re_forNumberDensity*Re_forNumberDensity;
  Re2_forGradient      = Re_forGradient*Re_forGradient;
  Re2_forLaplacian     = Re_forLaplacian*Re_forLaplacian;
  calNZeroAndLambda();
  FluidDensity       = FLUID_DENSITY;
  collisionDistance  = COLLISION_DISTANCE; 
  collisionDistance2 = collisionDistance*collisionDistance;
  FileNumber=0;
  Time=0.0;
}

void calNZeroAndLambda( void ){
  int iX, iY, iZ;
  int iZ_start, iZ_end;
  double xj, yj, zj, distance, distance2;
  double xi, yi, zi;

  if( DIM == 2 ){
    iZ_start = 0; iZ_end = 1;
  }else{
    iZ_start = -4; iZ_end = 5;
  }

  N0_forNumberDensity = 0.0;
  N0_forGradient      = 0.0;
  N0_forLaplacian     = 0.0;
  Lambda              = 0.0;
  xi = 0.0;  yi = 0.0;  zi = 0.0;

  for(iX= -4; iX<5; iX++){
    for(iY= -4;iY<5;iY++){
      for(iZ= iZ_start;iZ<iZ_end;iZ++){
        if( ((iX==0.0)&&(iY==0.0)) && (iZ==0.0) )continue;
        xj = PARTICLE_DISTANCE * (double)(iX);
        yj = PARTICLE_DISTANCE * (double)(iY);
        zj = PARTICLE_DISTANCE * (double)(iZ);
        distance2 = (xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi);
        distance = sqrt(distance2);
        N0_forNumberDensity += weight(distance, Re_forNumberDensity);
        N0_forGradient      += weight(distance, Re_forGradient);
        N0_forLaplacian     += weight(distance, Re_forLaplacian);
        Lambda              += distance2 * weight(distance, Re_forLaplacian);
      }
    }
  }
  Lambda = Lambda/N0_forLaplacian;
}

double weight( double distance, double re ){
  double weightIJ;

  if( distance >= re ){
    weightIJ = 0.0;
  }else{
    weightIJ = (re/distance) - 1.0;
  }
  return weightIJ;
}

void mainLoopOfSimulation( void ){
	FINISH_TIME = 1500.0;
  int iTimeStep = 0;

  // writeData_inVtuFormat();
  // writeData_inProfFormat();
  writeData_inCsvFormat();

  while(1){
    // calPressure();
    calTemperature();
    // moveSlider();

    iTimeStep++;
    Time += DT;
    if( (iTimeStep % OUTPUT_INTERVAL) == 0 ){
      printf("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticles: %d\n", iTimeStep, Time, NumberOfParticles);
      // writeData_inVtuFormat();
      // writeData_inProfFormat();
      writeData_inCsvFormat();
      FileNumber++;
    }

    if( Time >= FINISH_TIME ){break;}
  }
}


double rhoForTempCalc(int ParticleType, double Temp){
  // ParticleType can't be SLIDER
  //rho[kg/m^3]：密度
  
  double rho;

  if (ParticleType == AIR){
    rho = 1.27;
  }else{

    if (Temp < -0.01){
      rho = 917;
    }else if(0.01 < Temp){
      rho = 1000;
    }else{
      rho = 917 + (1000-917)*((Temp+0.01)/0.02);
    }

  }

  return rho;
}

double lambdaForTempCalc(int ParticleType, double Temp){
  // ParticleType can't be SLIDER
  //lambda[W/(m*K)]：熱伝導率
  
  double lambda;

  if (ParticleType == AIR){
    lambda = 0.0245;
  }else{
    if (Temp < -0.01){
      lambda = 2.19;
    }else if(0.01 < Temp){
      lambda = 0.576;
    }else{
      lambda = 2.19 + (0.576-2.19)*((Temp+0.01)/0.02);
    }
  }

  return lambda;
}

double cForTempCalc(int ParticleType, double Temp){
  // ParticleType can't be SLIDER
  //c[J/(kg*K)]：比熱
  
  double c;

  if (ParticleType == AIR){
    c = 1006;
  }else{
    if (Temp < -0.01){
      c = 2040;
    }else if(0.01 < Temp){
      c = 4200;
    }else{
      c = 2040 + (4200-2040)*((Temp+0.01)/0.02) + 334000/0.02;
    }
  }

  return c;
}


void calTemperature( void ){
  int i,j;
  double distance, distance2;
  double w;
  double xij, yij, zij;
  double dgd = (2.0*DIM)/(N0_forLaplacian*Lambda);
  double rho, lambda, c;

  for (i=0; i<NumberOfParticles; i++){
   if (i==0){
      Temperature[i] = -15.0;
      continue;
   }

    rho = rhoForTempCalc(ParticleType[i], Temperature[i]);
    lambda = lambdaForTempCalc(ParticleType[i], Temperature[i]);
    c = cForTempCalc(ParticleType[i], Temperature[i]);

    double UpdatedTemp = Temperature[i];

    for (j=0; j<NumberOfParticles; j++){
      if (j==i) continue;

      xij = Position[j*3] - Position[i*3];
      distance2 = xij*xij;
      distance = sqrt(distance2);
      w =  weight(distance, Re_forLaplacian);    

      UpdatedTemp += lambda/(rho*c)*dgd*(Temperature[j]-Temperature[i])*w*DT;
    }

    Temperature[i] = UpdatedTemp;

    if(i == NumberOfParticles-1){
      Temperature[i] = Temperature[i-1];
    }
  
    if (Temperature[i] < -0.01){
      ParticleType[i] = SOLID;
    }else if (0.01 < Temperature[i]){
      ParticleType[i] = FLUID;
    }else{
      ParticleType[i] = SOLID_AND_FLUID;
    }
  }
}


void writeData_inProfFormat( void ){
  int i;
  FILE *fp;
  char fileName[256];
  char dirName[] = "prof";

  sprintf(fileName, "%s/output_%04d.prof", dirName, FileNumber);
  fp = fopen(fileName, "w");
  fprintf(fp,"%lf\n",Time);
  fprintf(fp,"%d\n",NumberOfParticles);
  for(i=0;i<NumberOfParticles;i++) {
    fprintf(fp,"%d %lf %lf %lf %f %f\n"
	    ,ParticleType[i], Position[i*3], Position[i*3+1], Position[i*3+2], Pressure[i], Temperature[i]);
  }
  fclose(fp);
}


void writeData_inVtuFormat( void ){
  int i;
  double absoluteValueOfVelocity;
  FILE *fp;
  char fileName[1024];
  char dirName[] = "vtu";

  sprintf(fileName, "%s/output_%04d.vtu", dirName, FileNumber);
  fp=fopen(fileName,"w");
  fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",NumberOfParticles,NumberOfParticles);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%lf %lf %lf\n",Position[i*3],Position[i*3+1],Position[i*3+2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",ParticleType[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Temperature' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(double)Temperature[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='Int32' Name='offsets' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i+1);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='UInt8' Name='types' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"1\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</UnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void writeData_inCsvFormat( void ){
  int i;
  FILE *fp;

  fp = fopen("output.csv", "a");

  for(i=0;i<NumberOfParticles;i++) {
    if (ParticleType[i]==SOLID && (ParticleType[i+1]==SOLID_AND_FLUID || ParticleType[i+1]==FLUID)){
      fprintf(fp, "%f, %f\n", Time, Position[i*3]);
    }
  }
  fclose(fp);
}