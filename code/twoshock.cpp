#include "iostream"
#include <vector>
#include <fstream>
#include <math.h>
#include <ctime>
using namespace std ;

void pause(int dur)
{
int temp = time(NULL) + dur;

while(temp > time(NULL));
}

// Function to convert premetive(W) to conserved(U)
void W2U(vector<float> & W,vector<float> & U)
{
	U[0] = W[0];
	U[1] = W[0]*W[1];
	U[2] = W[2]/(W[4]-1) + 0.5*W[0]*W[1]*W[1];
	U[3] = W[3];
	U[4] = W[4];

}

// Function to convert conserved(U) to premetive(W) 
void U2W(vector<float> & U,vector<float> & W)
{
	W[0] = U[0];
	W[1] = U[1]/U[0];
	W[2] = (U[4]-1)*(U[2] - 0.5*W[0]*W[1]*W[1]);
	W[3] = U[3];
	W[4] = U[4];
}

// Calculating Euler Flux
void getEulerFlux(vector<float> & W, vector<float> & F)
{	
	float rho = W[0];
	float u = W[1];
	float p = W[2];
	float a2 = W[4]*p/rho;
	
	F[0] = rho*u;
	F[1] = rho*u*u + p;
	F[2] = rho*u*( a2/(W[4]-1.0) + 0.5*u*u );
}

// Use the Roe approximate Riemann solver to calculate Fluxes.
void getRoeFlux(vector<float> & WL,vector<float> & WR, vector<float> & FRoe, float inter)
{
	vector<float> UL(5);
	vector<float> UR(5);

	W2U(WL,UL);
	W2U(WR,UR);

	float gamma = inter*WR[4] + (1.0-inter)*WL[4];

	// Primitive and other variables.
    // Left state
    float rhoL = WL[0];
	float uL = WL[1];
	float pL = WL[2];
	float aL = sqrt(gamma*pL/rhoL);
	float eL = UL[2];
	float HL = (eL+pL)/rhoL;

	// Right state
    float rhoR = WR[0];
	float uR = WR[1];
	float pR = WR[2];
	float aR = sqrt(gamma*pR/rhoR);
	float eR = UR[2];
	float HR = (eR+pR)/rhoR;


	//Roe Averages
	float RT = sqrt(rhoR/rhoL);
	float rhoInt = RT*rhoL;
	float uInt = (uL+RT*uR)/(1+RT);
	float HInt = (HL+RT*HR)/(1+RT);
	float aInt = sqrt((gamma-1.0)*(HInt-0.5*uInt*uInt));

	// Differences in primitive variables	
	float drho = rhoR - rhoL;
	float du = uR - uL;
	float dp = pR - pL;

	// Wave strangth (Characterstic variables)
	float dU[3];
	dU[0] = 0.5*(dp - rhoInt*aInt*du)/(aInt*aInt);
	dU[1] = -( dp/(aInt*aInt) - drho );
	dU[2] = 0.5*(dp + rhoInt*aInt*du)/(aInt*aInt);

	// Absolute values of the wave speeds (Eigenvalues)
	float ws[3];
	ws[0] = fabs(uInt - aInt);
	ws[1] = fabs(uInt);
	ws[2] = fabs(uInt + aInt);

	// Modified wave speeds for nonlinear fields (the so-called entropy fix,
	// which is often implemented to remove non-physical expansion shocks).
    // There are various ways to implement the entropy fix. This is one of them
	
	#if 1
	float Da = max(0.0, 4.0*((uR-aR)-(uL-aL)));
	if (ws[0] < 0.5*Da)
	{
		ws[0] = ws[0]*ws[0]/Da + 0.25*Da ;
	}
	Da = max(0.0, 4.0*((uR+aR)-(uL+aL)));
	if (ws[2]<0.5*Da)
	{
		ws[2] = ws[2]*ws[2]/Da + 0.25*Da;
	}
	#endif

	// Right engenvectors
	float R[3][3];

	R[0][0] = 1.0;
	R[1][0] = uInt- aInt;
	R[2][0] = HInt - uInt*aInt;

	R[0][1] = 1.0;
	R[1][1] = uInt;
	R[2][1] = 0.5*uInt*uInt;

	R[0][2] = 1.0;
	R[1][2] = uInt + aInt;
	R[2][2] = HInt + uInt*aInt;

	// Compute the average Flux
	vector<float> FEulerL(3);
	vector<float> FEulerR(3);
	
	getEulerFlux(WL,FEulerL);
	getEulerFlux(WR,FEulerR);

	FRoe[0] = 0.5*(FEulerL[0]+FEulerR[0]);
	FRoe[1] = 0.5*(FEulerL[1]+FEulerR[1]);
	FRoe[2] = 0.5*(FEulerL[2]+FEulerR[2]);
	
	// Add the matrix dissipation term to complete the Roe Flux
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			FRoe[i] = FRoe[i] - 0.5*ws[j]*dU[j]*R[i][j];
		}
	}

}

// CFL condition to maintain stability
float getTimeStep(float CFL, float dx, int NCells, vector< vector<float> > & U)
{
	float maxSpeed = -1.0;
	vector<float> W(5);
	
	float u;
	float a;
	float lembdaMax;

	for (int i = 1; i < NCells-1; ++i)
	{
		U2W(U[i],W);
		u = W[1];
		a = sqrt(W[4]*W[2]/W[0]);
		lembdaMax = fabs(u)+a;
		maxSpeed = max(maxSpeed,lembdaMax);	
	}
	return CFL*dx/maxSpeed; // CFL condition
}

// Update solution
void update(vector<vector<float> > & U, vector<vector<float> > & Flux, float dt, float dx, int NCells)
{
	vector<float> WL(5);
	vector<float> WR(5);

	for (int i = 0; i < NCells-1; ++i) // Flux using gamma_L
	{
		U2W(U[i],WL);
		U2W(U[i+1],WR);
		getRoeFlux(WL,WR,Flux[i],0.); // FRoe[i]

	}

	#if 1 // Single gas
	for (int i = 1; i < NCells-1; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			U[i][j] = U[i][j] + (dt/dx) * (Flux[i-1][j] - Flux[i][j]);
		}
	}
	#endif

	#if 0 // Multi gas
	// For two gases
	vector<float> Rho(NCells);
	vector<float> RhoU(NCells);
	vector<float> Phi(NCells);

	int interface; // Location of the interface
	for (int i = 0; i < NCells; ++i)
	{
		Rho[i] = U[i][0];
		RhoU[i] = U[i][1];
		Phi[i] = U[i][3];
		
		//searching for interface
		if (U[i][3]*U[i+1][3] < 0)
		{
			interface = i;
		}
	}

	U2W(U[interface],WL);
	U2W(U[interface+1],WR);
	
	vector<float> FL(3); // Left flux for (interface + 1)th cell using gamma_R	
	getRoeFlux(WL,WR,FL,1.); 

	for (int i = 1; i < NCells-1; ++i)
	{
		if (i == interface+1) // use the FL for 
		{
			for (int j = 0; j < 3; ++j)
			{
				U[i][j] = U[i][j] + (dt/dx) * (FL[j] - Flux[i][j]);
			}	
		}
		else
		{
			for (int j = 0; j < 3; ++j)
			{
				U[i][j] = U[i][j] + (dt/dx) * (Flux[i-1][j] - Flux[i][j]);
			}
		}
	}
	
	// Updating Phi 
	for (int i = 1; i < NCells-2; ++i)
	{
		U[i][3] = (Rho[i]*Phi[i]+(dt/dx)*((RhoU[i+1]*Phi[i+1])-(RhoU[i]*Phi[i])))/(U[i][0]);
	}

	Updating Gamma
	for (int i = 0; i < NCells; ++i)
	{
		if (U[i][3] > 0)
		{
			U[i][4] = 1.4; 
		}
		else
		{
			U[i][4] = 1.4;
		}
	}
	#endif

	// Applying Boundary Condition
	for (int j = 0; j < 5; ++j)
	{
		U[0][j] = U[1][j];	
		U[NCells-1][j] = U[NCells-2][j]; 
	}	

}

void TestCase(vector<vector<float> > & U)
{
	std::vector<float> WL0(5);
	std::vector<float> WM0(5);
	std::vector<float> WR0(5);

	// Case 1 : Simple moving 2 Shocks BVP
	// Left initial values
	WL0[0] = 1. ; // rho
	WL0[1] = 0. ; // u
	WL0[2] = 100000.0 ; // p
	WL0[3] = 1. ; // psi
	WL0[4] = 1.4 ; // gamma


	// initial values inside 
	WM0[0] = 0.125 ; 
	WM0[1] = 0. ; 
	WM0[2] = 10000.0 ; 
	WM0[3] = -1. ; 
	WM0[4] = 1.4 ;

	// Right initial values
	WR0[0] = 2. ; 
	WR0[1] = 0. ; 
	WR0[2] = 200000.0 ; 
	WR0[3] = -1. ; 
	WR0[4] = 1.4 ;

	// Case 2 : Multigas 1D shock
	// // Left initial values
	// WL0[0] = 1. ; // rho
	// WL0[1] = 1. ; // u
	// WL0[2] = 1. ; // p
	// WL0[3] = 1. ; // psi
	// WL0[4] = 1.4; // gamma


	// // Right initial values
	// WR0[0] = 0.1 ; 
	// WR0[1] = 1. ; 
	// WR0[2] = 1. ; 
	// WR0[3] = -1. ; 
	// WR0[4] = 1.4 ;

	// Case 3 : More stiff multi phase 
	// Left initial values
	// WL0[0] = 1000. ; // rho
	// WL0[1] = 0. ; // u
	// WL0[2] = 1e9 ; // p
	// WL0[3] = 1. ; // psi
	// WL0[4] = 4.4 ; // gamma


	// // Right initial values
	// WR0[0] = 50. ; 
	// WR0[1] = 0. ; 
	// WR0[2] = 1e5 ; 
	// WR0[3] = -1. ; 
	// WR0[4] = 4.4 ;

	int NCells = U.size();

	for (int i = 0; i < NCells; ++i)
	{
		if(i<NCells/4)
		{
			W2U(WL0,U[i]);
		}
		else if(i>=NCells/4 && i<3*NCells/4)
		{
			W2U(WM0,U[i]);
		}
		else
		{
			W2U(WR0,U[i]);
		}
	}
}

int main()
{
	time_t start, end ;
	time(&start); // noteing the starting time

	//INPUTS = 1.4 ;
	int NCells = 10000;
	float dx = 5000.0/NCells;

	float CFL = 0.8;
	float t = 0.0; 
	float tEnd = 3.8;
	int nSteps = 0; 
	float dt;

	typedef vector<float> Dim1;
	typedef vector<Dim1> Dim2;

	Dim2 U(NCells,Dim1(5)); // To store the conserved variables 
	Dim2 W(NCells,Dim1(5)); // To store the conserved variables 
	Dim2 Flux(NCells,Dim1(3)); // To store the Flux


	TestCase(U); // setting the initial values
	
	// Solver
	int framcout = 0;
	int stepcount = 0;
	while(t<tEnd)
	{
		stepcount++;
		dt = getTimeStep(CFL,dx,NCells,U);
		t = t+dt;
		update(U,Flux,dt,dx,NCells);
		cout<< "t" << t << endl; 
		
		if (stepcount%100 == 0)
		{
			framcout ++;
			cout<< framcout << endl;
			ofstream WriteW ;
			WriteW.open("premetive.csv");
			// WriteW << "x" << "," << "rho" << "," << "u" << ","<< "e" << "," << "psi" << endl ;
			for (int i = 0; i < NCells; ++i)
			{
				U2W(U[i],W[i]);
				WriteW << dx*i << "," << W[i][0] << "," << W[i][1] <<","<< W[i][2] <<","<< W[i][3] <<","<< W[i][4] << endl ;
			}
			pause(1);
		}
		
	
	}
	
	//OUTPUTS IN CSV FORMET
	ofstream WriteU ;
	WriteU.open("conserved.csv");
	// WriteU << "x" << "," << "rho" << "," << "u" << ","<< "e" << "," << "psi" << endl ;
	for (int i = 0; i < NCells; ++i)
	{
		WriteU << dx*i << "," << U[i][0] << "," << U[i][1] <<","<< U[i][2] <<","<< U[i][3] <<","<< U[i][4] << endl ;
	}
	
	time(&end) ;
	double diff = difftime (end,start);
	cout << "Time taken by the solver in secs = " << diff << endl ;	

	return 0;
}
