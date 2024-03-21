#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DIM                  2
#define PARTICLE_DISTANCE    0.001
#define DT                   0.001
#define OUTPUT_INTERVAL      10

#define ARRAY_SIZE           10000
#define FINISH_TIME          1.0 
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
#define AIR    3
#define SLIDER 4
#define WALL 5
#define DUMMY_WALL 6
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
double stayInBoundary(double x, double min, double max);
void writeData_inProfFormat( void );
void writeData_inVtuFormat( void );

void calTemperature(void);
void moveSlider( void );

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

// Parameters

#define df_coef  0.05
#define sliderSpeed  5.0
#define iceTemp -1.0
#define sliderTemp -1.0

#define x_ice_left  0.0
#define x_ice_right  0.3-PARTICLE_DISTANCE
#define y_ice_bottom  0.0
#define y_ice_top  0.005-PARTICLE_DISTANCE
#define x_slider_left   -0.3
#define x_slider_right   0.0-PARTICLE_DISTANCE
#define y_slider_bottom   0.005
#define y_slider_top   0.007-PARTICLE_DISTANCE


void initializeParticlePositionAndVelocity_for2dim( void ){
  int iX, iY;
  int nX, nY;
  double x, y, z;
  int i = 0;
  int flagOfParticleGeneration;

  nX = (int)(1.0/PARTICLE_DISTANCE);  
  nY = (int)(0.1/PARTICLE_DISTANCE);

  for(iX= 0; iX<nX; iX++){
    for(iY= 0; iY<nY; iY++){
      x = PARTICLE_DISTANCE * (double)(iX) - 0.3;
      y = PARTICLE_DISTANCE * (double)(iY) - 0.005;
      z = 0.0;
      flagOfParticleGeneration = OFF;

      /*
      // dummy wall region
      if ((x>=x_ice_left-4*PARTICLE_DISTANCE-EPS && x<=x_ice_right+4*PARTICLE_DISTANCE+EPS) && (y>=y_ice_bottom-4*PARTICLE_DISTANCE-EPS && y<=y_ice_top-EPS)){
        ParticleType[i]=DUMMY_WALL;
        flagOfParticleGeneration = ON;
        Temperature[i] = iceTemp;
      }

      // wall region 
      if ((x>=x_ice_left-2*PARTICLE_DISTANCE-EPS && x<=x_ice_right+2*PARTICLE_DISTANCE+EPS) && (y>=y_ice_bottom-2*PARTICLE_DISTANCE-EPS && y<=y_ice_top-EPS)){
        ParticleType[i]=DUMMY_WALL;
        flagOfParticleGeneration = ON;
        Temperature[i] = iceTemp;
      }
      
      // solid region 
      if ((x>=x_ice_left-EPS && x<=x_ice_right+EPS) && ((y>=y_ice_bottom-EPS && y<y_ice_top-EPS)||(iX%100==0 && fabs(y-y_ice_top)<EPS))){ 
        ParticleType[i]=SOLID;
        flagOfParticleGeneration = ON;
        Temperature[i] = iceTemp;
      }

      // slider region 
      if( (x>=x_slider_left-EPS && x<=x_slider_right+EPS) && (y>=y_slider_bottom-EPS && y<=y_slider_top+EPS) ){
        ParticleType[i] = SLIDER;
        Temperature[i] = sliderTemp;
        flagOfParticleGeneration = ON;
      } 

      // air region 
      if((x>=x_ice_left-EPS && x<=x_ice_right+EPS) && (y>=y_slider_bottom-EPS && y<=y_slider_top+EPS)){
        ParticleType[i] = AIR;
        Temperature[i] = sliderTemp;
        flagOfParticleGeneration = ON;
      } 
      */

      // dummy wall region
      if ((x>=-0.004 && x<=0.304) && (y>=-0.004 && y<=0.004)){
        ParticleType[i]=DUMMY_WALL;
        flagOfParticleGeneration = ON;
        Temperature[i] = iceTemp;
      }

      // wall region 
      if ((x>=-0.002 && x<=0.302) && (y>=-0.002 && y<=0.004)){
        ParticleType[i]=WALL;
        flagOfParticleGeneration = ON;
        Temperature[i] = iceTemp;
      }
      
      // solid region 
      if ((x>=0.0 && x<=0.3) && ((y>=0.0 && y<0.004)||(iX%50==0 && fabs(y-0.004)<EPS))){ 
        ParticleType[i]=SOLID;
        flagOfParticleGeneration = ON;
        Temperature[i] = iceTemp;
      }

      // slider region 
      if( (x>=-0.3 && x<0.0) && (y>=0.005 && y<=0.007) ){
        ParticleType[i] = SLIDER;
        Temperature[i] = sliderTemp;
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
  xi = -0.004;  yi = -0.004;  zi = 0.0;

  for(iX= -4; iX<5; iX++){
    for(iY= -4; iY<5; iY++){  
      for(iZ= iZ_start;iZ<iZ_end;iZ++){
        if((iX==-4 && iY==-4) && (iZ==0)) continue;
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

double sliderLeft = x_slider_left;
double sliderRight = x_slider_right;

void mainLoopOfSimulation( void ){
  int iTimeStep = 0;

  writeData_inVtuFormat();
  writeData_inProfFormat();

  while(1){
    moveSlider();
    calTemperature();
    calGravity();
    calViscosity();
    moveParticle();
    collision();
    calPressure();
    calPressureGradient();
    moveParticleUsingPressureGradient();

    iTimeStep++;
    Time += DT;
    if( (iTimeStep % OUTPUT_INTERVAL) == 0 ){
      printf("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticles: %d\n", iTimeStep, Time, NumberOfParticles);
      writeData_inVtuFormat();
      writeData_inProfFormat();
    }
    if( Time >= FINISH_TIME ){break;}
  }
}


void calGravity( void ){
  int i;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      Acceleration[i*3  ]=G_X;
      Acceleration[i*3+1]=G_Y;
      Acceleration[i*3+2]=G_Z;
    }else{
      Acceleration[i*3  ]=0.0;
      Acceleration[i*3+1]=0.0;
      Acceleration[i*3+2]=0.0;
    }
  }
}


void calViscosity( void ){
  int i,j;
  double viscosityTerm_x, viscosityTerm_y, viscosityTerm_z;
  double distance, distance2;
  double w;
  double xij, yij, zij;
  double a;

  a = (KINEMATIC_VISCOSITY)*(2.0*DIM)/(N0_forLaplacian*Lambda);

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] != FLUID) continue;
    viscosityTerm_x = 0.0;  viscosityTerm_y = 0.0;  viscosityTerm_z = 0.0;

    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      if(distance<Re_forLaplacian){
        w =  weight(distance, Re_forLaplacian);
        viscosityTerm_x +=(Velocity[j*3  ]-Velocity[i*3  ])*w;
        viscosityTerm_y +=(Velocity[j*3+1]-Velocity[i*3+1])*w;
        viscosityTerm_z +=(Velocity[j*3+2]-Velocity[i*3+2])*w;
      }
    }
    viscosityTerm_x = viscosityTerm_x * a;
    viscosityTerm_y = viscosityTerm_y * a;
    viscosityTerm_z = viscosityTerm_z * a;
    Acceleration[i*3  ] += viscosityTerm_x;
    Acceleration[i*3+1] += viscosityTerm_y;
    Acceleration[i*3+2] += viscosityTerm_z;
  }
}

double stayInBoundary(double x, double min, double max){
  if (x < min){
    return min;
  }else if (x > max){
    return max;
  }else{
    return x;
  }
}

void moveParticle( void ){
  int i;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      Velocity[i*3  ] += Acceleration[i*3  ]*DT; 
      Velocity[i*3+1] += Acceleration[i*3+1]*DT; 
      Velocity[i*3+2] += Acceleration[i*3+2]*DT;

      Position[i*3] += Velocity[i*3]*DT;
      Position[i*3+1] += Velocity[i*3+1]*DT;
      Position[i*3  ] = stayInBoundary(Position[i*3], 0.0, 0.3);
      Position[i*3+1] = stayInBoundary(Position[i*3+1], 0.0, 0.005);
      Position[i*3+2] += Velocity[i*3+2]*DT;
    }
   
    Acceleration[i*3  ]=0.0;
    Acceleration[i*3+1]=0.0;
    Acceleration[i*3+2]=0.0;
  } 
}


void collision( void ){
  int    i,j;
  double xij, yij, zij;
  double distance,distance2;
  double forceDT; /* forceDT is the impulse of collision between particles */
  double mi, mj;
  double velocity_ix, velocity_iy, velocity_iz;
  double e = COEFFICIENT_OF_RESTITUTION;
  static double VelocityAfterCollision[3*ARRAY_SIZE];

  for(i=0;i<3*NumberOfParticles;i++){ 
    VelocityAfterCollision[i] = Velocity[i]; 
  }

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      mi = FluidDensity;
      velocity_ix = Velocity[i*3  ];  
      velocity_iy = Velocity[i*3+1];  
      velocity_iz = Velocity[i*3+2];

  
      for(j=0;j<NumberOfParticles;j++){
        if( (j==i) || (ParticleType[j]==GHOST) ) continue;
        xij = Position[j*3  ] - Position[i*3  ];
        yij = Position[j*3+1] - Position[i*3+1];
        zij = Position[j*3+2] - Position[i*3+2];
        distance2 = (xij*xij) + (yij*yij) + (zij*zij);
        if(distance2<collisionDistance2){
          distance = sqrt(distance2);
          forceDT = (velocity_ix-Velocity[j*3  ])*(xij/distance)
                  +(velocity_iy-Velocity[j*3+1])*(yij/distance)
                  +(velocity_iz-Velocity[j*3+2])*(zij/distance);
          if(forceDT > 0.0){
            mj = FluidDensity;
            forceDT *= (1.0+e)*mi*mj/(mi+mj);
            velocity_ix -= (forceDT/mi)*(xij/distance); 
            velocity_iy -= (forceDT/mi)*(yij/distance); 
            velocity_iz -= (forceDT/mi)*(zij/distance);
            /*
            if(j>i){ fprintf(stderr,"WARNING: Collision occurred between %d and %d particles.\n",i,j); }
            */
          }
        }
      }
      VelocityAfterCollision[i*3  ] = velocity_ix; 
      VelocityAfterCollision[i*3+1] = velocity_iy; 
      VelocityAfterCollision[i*3+2] = velocity_iz;
    }
  }

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      Position[i*3] += Velocity[i*3]*DT;
      Position[i*3+1] += Velocity[i*3+1]*DT;
      // Position[i*3  ] = stayInBoundary(Position[i*3], 0.0, 0.3);
      // Position[i*3+1] = stayInBoundary(Position[i*3+1], 0.0, 0.005);
      Position[i*3+2] += Velocity[i*3+2]*DT;

      Velocity[i*3  ] = VelocityAfterCollision[i*3  ]; 
      Velocity[i*3+1] = VelocityAfterCollision[i*3+1]; 
      Velocity[i*3+2] = VelocityAfterCollision[i*3+2];
    }
  }
}

void calPressure( void ){
  calNumberDensity();
  setBoundaryCondition();
  setSourceTerm();
  setMatrix();
  solveSimultaneousEquationsByGaussEliminationMethod();
  removeNegativePressure();
  setMinimumPressure();
}


void calNumberDensity( void ){
  int    i,j;
  double xij, yij, zij;
  double distance, distance2;
  double w;

  for(i=0;i<NumberOfParticles;i++){
    NumberDensity[i] = 0.0;
    if(ParticleType[i] == GHOST) continue;


    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      w =  weight(distance, Re_forNumberDensity);
      NumberDensity[i] += w;
    }
  }
}


void setBoundaryCondition( void ){
  int i;
  double n0 = N0_forNumberDensity;
  double beta = THRESHOLD_RATIO_OF_NUMBER_DENSITY;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i]==GHOST){
      BoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if( NumberDensity[i] < beta * n0 ){
      BoundaryCondition[i]=SURFACE_PARTICLE;
    }else{
      BoundaryCondition[i]=INNER_PARTICLE;
    }
  }
}


void setSourceTerm( void ){
  int i;
  double n0    = N0_forNumberDensity;
  double gamma = RELAXATION_COEFFICIENT_FOR_PRESSURE;

  for(i=0;i<NumberOfParticles;i++){
    SourceTerm[i]=0.0;
    if(ParticleType[i]==GHOST) continue;
    if(BoundaryCondition[i]==INNER_PARTICLE){
      SourceTerm[i] = gamma * (1.0/(DT*DT))*((NumberDensity[i]-n0)/n0);
    }else if(BoundaryCondition[i]==SURFACE_PARTICLE){
      SourceTerm[i]=0.0;
    }
  }
}


void setMatrix( void ){
  double xij, yij, zij;
  double distance, distance2;
  double coefficientIJ;
  double n0 = N0_forLaplacian;
  int    i,j;
  double a;
  int n = NumberOfParticles;

  for(i=0;i<NumberOfParticles;i++){
    for(j=0;j<NumberOfParticles;j++){
      CoefficientMatrix[i*n+j] = 0.0;
    }
  }

  a = 2.0*DIM/(n0*Lambda);

  for(i=0;i<NumberOfParticles;i++){
    if(BoundaryCondition[i] != INNER_PARTICLE) continue;

    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (BoundaryCondition[j]==GHOST_OR_DUMMY) ) continue;
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij)+(yij*yij)+(zij*zij);
      distance  = sqrt(distance2);
      if(distance>=Re_forLaplacian)continue;
      coefficientIJ = a * weight(distance, Re_forLaplacian)/FluidDensity;
      CoefficientMatrix[i*n+j]  = (-1.0)*coefficientIJ;
      CoefficientMatrix[i*n+i] += coefficientIJ;
    }
    CoefficientMatrix[i*n+i] += (COMPRESSIBILITY)/(DT*DT);
  }
  exceptionalProcessingForBoundaryCondition();
}


void exceptionalProcessingForBoundaryCondition( void ){
  /* If there is no Dirichlet boundary condition on the fluid, 
     increase the diagonal terms of the matrix for an exception. This allows us to solve the matrix without Dirichlet boundary conditions. */
  checkBoundaryCondition();
  increaseDiagonalTerm();
}


void checkBoundaryCondition( void ){
  int i,j,count;
  double xij, yij, zij, distance2;

  for(i=0;i<NumberOfParticles;i++){
    if (BoundaryCondition[i]==GHOST_OR_DUMMY){
      FlagForCheckingBoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if (BoundaryCondition[i]==SURFACE_PARTICLE){
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CONNECTED;
    }else{
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
    }
  }

  do {
    count=0;

    for(i=0;i<NumberOfParticles;i++){
      if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_CONNECTED){
        for(j=0;j<NumberOfParticles;j++){
          if( j==i ) continue;
          if(ParticleType[j]==GHOST) continue;
          if(FlagForCheckingBoundaryCondition[j]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
            xij = Position[j*3  ] - Position[i*3  ];
            yij = Position[j*3+1] - Position[i*3+1];
            zij = Position[j*3+2] - Position[i*3+2];
            distance2 = (xij*xij)+(yij*yij)+(zij*zij);
            if(distance2>=Re2_forLaplacian)continue;
            FlagForCheckingBoundaryCondition[j]=DIRICHLET_BOUNDARY_IS_CONNECTED;
          }
        }
        FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CHECKED;
        count++;
      }
    }
  } while (count!=0); /* This procedure is repeated until the all fluid or wall particles (which have Dirichlet boundary condition in the particle group) are in the state of "DIRICHLET_BOUNDARY_IS_CHECKED".*/

  for(i=0;i<NumberOfParticles;i++){
    if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
      fprintf(stderr,"WARNING: There is no dirichlet boundary condition for %d-th particle.\n",i );
    }
  }
}


void increaseDiagonalTerm( void ){
  int i;
  int n = NumberOfParticles;

  for(i=0;i<n;i++) {
    if(FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED ){
      CoefficientMatrix[i*n+i] = 2.0 * CoefficientMatrix[i*n+i];
    }
  }
}


void solveSimultaneousEquationsByGaussEliminationMethod( void ){
  int    i,j,k;
  double c;
  double sumOfTerms;
  int    n = NumberOfParticles;

  for(i=0; i<n; i++){ 
    Pressure[i] = 0.0; 
  }

  for(i=0; i<n-1; i++){
    if ( BoundaryCondition[i] != INNER_PARTICLE ) continue;

    for(j=i+1; j<n; j++){
      if(BoundaryCondition[j]==GHOST_OR_DUMMY) continue;
      c = CoefficientMatrix[j*n+i]/CoefficientMatrix[i*n+i];
  
      for(k=i+1; k<n; k++){
	      CoefficientMatrix[j*n+k] -= c * CoefficientMatrix[i*n+k];
      }
      SourceTerm[j] -= c*SourceTerm[i];
    }
  }

  for( i=n-1; i>=0; i--){
    if ( BoundaryCondition[i] != INNER_PARTICLE ) continue;
    sumOfTerms = 0.0;

    for( j=i+1; j<n; j++ ){
      if(BoundaryCondition[j]==GHOST_OR_DUMMY) continue;
      sumOfTerms += CoefficientMatrix[i*n+j] * Pressure[j];
    }
    Pressure[i] = (SourceTerm[i] - sumOfTerms)/CoefficientMatrix[i*n+i];
  }
}


void removeNegativePressure( void ){
  int i;

  for(i=0;i<NumberOfParticles;i++) {
    if(Pressure[i]<0.0)Pressure[i]=0.0;
  }
}


void setMinimumPressure( void ){
  double xij, yij, zij, distance2;
  int i,j;

  for(i=0;i<NumberOfParticles;i++) {
    if(ParticleType[i]==GHOST)continue;
    MinimumPressure[i]=Pressure[i];


    for(j=0;j<NumberOfParticles;j++) {
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;

      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij)+(yij*yij)+(zij*zij);

      if(distance2>=Re2_forGradient)continue;

      if( MinimumPressure[i] > Pressure[j] ){
	      MinimumPressure[i] = Pressure[j];
      }
    }
  }
}


void calPressureGradient( void ){
  int    i,j;
  double gradient_x, gradient_y, gradient_z;
  double xij, yij, zij;
  double distance, distance2;
  double w,pij;
  double a;

  a =DIM/N0_forGradient;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] != FLUID) continue;
    gradient_x = 0.0;  gradient_y = 0.0;  gradient_z = 0.0;


    for(j=0;j<NumberOfParticles;j++){
      if( j==i ) continue;
      if( ParticleType[j]==GHOST ) continue;

      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      if(distance<Re_forGradient){
        w =  weight(distance, Re_forGradient);
        pij = (Pressure[j] - MinimumPressure[i])/distance2;
        gradient_x += xij*pij*w;
        gradient_y += yij*pij*w;
        gradient_z += zij*pij*w;
      }
    }
    gradient_x *= a;
    gradient_y *= a;
    gradient_z *= a;
    Acceleration[i*3  ]= (-1.0)*gradient_x/FluidDensity;
    Acceleration[i*3+1]= (-1.0)*gradient_y/FluidDensity;
    Acceleration[i*3+2]= (-1.0)*gradient_z/FluidDensity;
  }
}


void moveParticleUsingPressureGradient( void ){
  int i;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      Velocity[i*3  ] +=Acceleration[i*3  ]*DT;
      Velocity[i*3+1] +=Acceleration[i*3+1]*DT;
      Velocity[i*3+2] +=Acceleration[i*3+2]*DT;

      Position[i*3  ] +=Acceleration[i*3  ]*DT*DT;
      Position[i*3+1] +=Acceleration[i*3+1]*DT*DT;
      Position[i*3+2] +=Acceleration[i*3+2]*DT*DT;

      Position[i*3  ] = stayInBoundary(Position[i*3], 0.0, 0.3);
      Position[i*3+1] = stayInBoundary(Position[i*3+1], 0.0, 0.005);
    }
    Acceleration[i*3  ]=0.0;
    Acceleration[i*3+1]=0.0;
    Acceleration[i*3+2]=0.0;
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
    if (ParticleType[i]==SLIDER){
      Temperature[i] = sliderTemp;
      continue;
    }

    rho = rhoForTempCalc(ParticleType[i], Temperature[i]);
    lambda = lambdaForTempCalc(ParticleType[i], Temperature[i]);
    c = cForTempCalc(ParticleType[i], Temperature[i]);

    double UpdatedTemp = Temperature[i];

    if ((Position[i*3]<=sliderRight) && (Position[i*3]>=sliderLeft) && (Position[i*3+1]>0.005-2*PARTICLE_DISTANCE) && (ParticleType[i]!=WALL) && (ParticleType[i]!=DUMMY_WALL)){
      UpdatedTemp += df_coef*0.5*9.8*2.0*sliderSpeed*DT/(2*rho*c*pow(PARTICLE_DISTANCE, 3)); 
    }  

    for (j=0; j<NumberOfParticles; j++){
      if (j==i || (ParticleType[j]==SLIDER && sliderRight<0)) continue;

      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      w =  weight(distance, Re_forLaplacian);    

      UpdatedTemp += 2.7*lambda/(rho*c)*dgd*(Temperature[j]-Temperature[i])*w*DT;
    }

    Temperature[i] = UpdatedTemp;

    if (ParticleType[i]==AIR || ParticleType[i]==DUMMY_WALL || ParticleType[i]==WALL){
      continue;
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


void moveSlider(void){
	double travelDistance = sliderSpeed*DT;

  for (int i=0; i<NumberOfParticles; i++){
		if (ParticleType[i] == SLIDER){
      Position[i*3] += travelDistance;
    }
  }

  sliderLeft += travelDistance;
  sliderRight += travelDistance;  
}


void writeData_inProfFormat( void ){
  int i;
  FILE *fp;
  char fileName[256];

  sprintf(fileName, "prof/output_%04d.prof",FileNumber);
  fp = fopen(fileName, "w");
  fprintf(fp,"%lf\n",Time);
  fprintf(fp,"%d\n",NumberOfParticles);
  for(i=0;i<NumberOfParticles;i++) {
    fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n"
	    ,ParticleType[i], Position[i*3], Position[i*3+1], Position[i*3+2]
	    ,Velocity[i*3], Velocity[i*3+1], Velocity[i*3+2], Pressure[i], NumberDensity[i], Temperature[i]);
  }
  fclose(fp);
  FileNumber++;
}


void writeData_inVtuFormat( void ){
  int i;
  double absoluteValueOfVelocity;
  FILE *fp;
  char fileName[1024];

  sprintf(fileName, "vtu/particle_%04d.vtu", FileNumber);
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
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    absoluteValueOfVelocity=
      sqrt( Velocity[i*3]*Velocity[i*3] + Velocity[i*3+1]*Velocity[i*3+1] + Velocity[i*3+2]*Velocity[i*3+2] );
    fprintf(fp,"%f\n",(float)absoluteValueOfVelocity);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Pressure' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(float)Pressure[i]);
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
