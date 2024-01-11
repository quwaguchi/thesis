#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DIM                  2
#define PARTICLE_DISTANCE    0.0125 //0.025
#define DT                   0.001
#define OUTPUT_INTERVAL      20

#define ARRAY_SIZE           10000
#define FINISH_TIME          20.0 
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
#define WALL    3
#define DUMMY_WALL  4
#define SLIDER  5
#define GHOST_OR_DUMMY  -1
#define SURFACE_PARTICLE 1    
#define INNER_PARTICLE   0      
#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0 
#define DIRICHLET_BOUNDARY_IS_CONNECTED     1 
#define DIRICHLET_BOUNDARY_IS_CHECKED       2 

// parameters for friction heat
double df_coef = 0.05;
double load = 60;
double stableSpeed = 1.0;

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

static double sliderPosition[ARRAY_SIZE];

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
double FlagForUpdate[ARRAY_SIZE] = {};

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

  nX = (int)(2.0/PARTICLE_DISTANCE)+5;  
  nY = (int)(0.5/PARTICLE_DISTANCE)+5;
  for(iX= -4; iX<nX; iX+=3){
    for(iY= -4;iY<nY;iY++){
      x = PARTICLE_DISTANCE * (double)(iX);
      y = PARTICLE_DISTANCE * (double)(iY);
      z = 0.0;
      flagOfParticleGeneration = OFF;
      
      /* solid region */
      if( ((x>=0.0)&&(x<2.00)) &&((y>=0.0)&&(y<0.40)) ){ 
        ParticleType[i]=SOLID;
        flagOfParticleGeneration = ON;
        FlagForUpdate[i] = 1;
        Temperature[i] = -5;
      }

      /* slider region */
      if( ((x>=-4.0*PARTICLE_DISTANCE)&&(x<0.50))&&( (y>=0.4)&&(y<0.5)) ){
        ParticleType[i] = SLIDER;
        Temperature[i] = -3;
        flagOfParticleGeneration = ON;
      }  

      if( flagOfParticleGeneration == ON){
        Position[i*3]=x; Position[i*3+1]=y; Position[i*3+2]=z;
        i++;
      }
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

  for(iX= -4;iX<5;iX++){
    for(iY= -4;iY<5;iY++){
      for(iZ= iZ_start;iZ<iZ_end;iZ++){
        if( ((iX==0)&&(iY==0)) && (iZ==0) )continue;
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

double sliderLeft = -4.0*PARTICLE_DISTANCE;
double sliderRight = 0.5;

void mainLoopOfSimulation( void ){
  int iTimeStep = 0;

  writeData_inVtuFormat();
  writeData_inProfFormat();

  while(1){
    calTemperature();
    moveSlider();

    iTimeStep++;
    Time += DT;
    if( (iTimeStep % OUTPUT_INTERVAL) == 0 ){
      printf("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticles: %d\n", iTimeStep, Time, NumberOfParticles);
      writeData_inVtuFormat();
      writeData_inProfFormat();
    }

    sliderLeft += stableSpeed*DT;
    sliderRight += stableSpeed*DT;

    if( Time >= FINISH_TIME ){break;}
  }
}


void calTemperature( void ){
  int i,j;
  double distance, distance2;
  double w;
  double xij, yij, zij;
  double dgd = (2.0*DIM)/(N0_forLaplacian*Lambda);
  double Ts = -0.01;
  double Tl = 0.01;
  double lambda, rho, c;

  for (i=0; i<NumberOfParticles; i++){
    if (!FlagForUpdate[i]) continue;

    double UpdatedTemp = Temperature[i];

    for (j=0; j<NumberOfParticles; j++){
      if (j==i) continue;

      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      w =  weight(distance, Re_forLaplacian);

      if (UpdatedTemp < Ts){
        lambda = 2.19;
        rho = 917;
        c = 2040;
      }else if(Tl < UpdatedTemp){
        lambda = 0.576;
        rho = 1000;
        c = 4200;
      }else{
        lambda = 2.19 + (0.576-2.19)/(Tl-Ts)*(UpdatedTemp-Ts);
        rho = 917 + (1000-917)/(Tl-Ts)*(UpdatedTemp-Ts);
        c = 2040 + (4200-2040)/(Tl-Ts)*(UpdatedTemp-Ts) + 334000/(Tl-Ts);
      }

      //consider friction heat and change the formula below later
      UpdatedTemp += lambda/(rho*c)*10e3*dgd*(Temperature[j]-Temperature[i])*w*DT;

      if ((sliderLeft<=Position[i*3])&&(Position[i*3]<=sliderRight)&&(0.4-PARTICLE_DISTANCE<=Position[i*3+1])&&(Position[i*3+1]<=0.4+PARTICLE_DISTANCE)){
        UpdatedTemp += df_coef*load*stableSpeed*DT;
      }
    }

    Temperature[i] = UpdatedTemp;

    if (Temperature[i] < Ts){
      ParticleType[i] = SOLID;
    }else if (Tl < Temperature[i]){
      ParticleType[i] = FLUID;
    }else{
      ParticleType[i] = SOLID_AND_FLUID;
    }
  }
}

void moveSlider(void){   
   for (int i=0; i<NumberOfParticles; i++){
     if (ParticleType[i] == SLIDER){
       Position[i*3] += stableSpeed*DT;
     }
   }
}

void writeData_inProfFormat( void ){
  int i;
  FILE *fp;
  char fileName[256];

  sprintf(fileName, "output_%04d.prof",FileNumber);
  fp = fopen(fileName, "w");
  fprintf(fp,"%lf\n",Time);
  fprintf(fp,"%d\n",NumberOfParticles);
  for(i=0;i<NumberOfParticles;i++) {
    fprintf(fp,"%d %lf %lf %lf %f\n"
	    ,ParticleType[i], Position[i*3], Position[i*3+1], Position[i*3+2], Temperature[i]);
  }
  fclose(fp);
  FileNumber++;
}


void writeData_inVtuFormat( void ){
  int i;
  double absoluteValueOfVelocity;
  FILE *fp;
  char fileName[1024];

  sprintf(fileName, "particle_%04d.vtu", FileNumber);
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