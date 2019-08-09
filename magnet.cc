#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <math.h>       /* exp */
using namespace std;
using namespace itensor;
int main()
    {

int Nx;

double Jx;

double Jy;

double Jz;

double hhx;

double hhy;

double hhz;

int MSW;

int jx1;
int jy1;


cout << "Number of X unit cells Nx";
cin >> Nx;
cout << "The magnetic coupling Jx";
cin >> Jx;
cout << "The magnetic coupling Jy";
cin >> Jy;
cout << "The magnetic coupling Jz";
cin >> Jz;
cout << "The magnetic field along x direction";
cin >> hhx;
cout << "The magnetic field along y direction";
cin >> hhy;
cout << "The magnetic field along z direction";
cin >> hhz;

ofstream myfile;
myfile.open ("Nx_"+std::to_string(Nx)+"_hhx_"+std::to_string(hhx)+"_hhy_"+std::to_string(hhy)+"_hhz_"+std::to_string(hhz)+"f.dat");

int N = 2*Nx*Nx;  

auto sites = SpinHalf(N);

auto ampoH = AutoMPO(sites);
for(int jx = 0; jx <Nx; ++jx)
  for(int jy = 0; jy <Nx; ++jy)
    {
    int NA=2*(jx+Nx*jy)+1;
    int NBX=2*(jx+Nx*jy)+2;

    if (jx-1<0) {jx1=Nx-1;} else {jx1=jx-1;}
    if (jy+1>=Nx) {jy1=0;} else {jy1=jy+1;}

    int NBY=2*(jx1+Nx*jy)+2;
    int NBZ=2*(jx1+Nx*jy1)+2;

    ampoH +=   Jx,  "Sx",NA,"Sx",NBX;

    ampoH +=   Jy,  "Sy",NA,"Sy",NBY;

    ampoH +=   Jz, "Sz",NA,"Sz",NBZ;

    println(NA,NBX);
    println(NA,NBY);
    println(NA,NBZ);
}

auto H0=MPO(ampoH);



//Set up random initial wavefunction

auto state = InitState(sites);
for(int i = 1; i <= N/2; ++i) 
        {
            state.set(2*i,"Up");
            state.set(2*(i-1)+1,"Dn");
        }

auto psi = MPS(state);

normalize(psi);

auto sweeps = Sweeps(20);
    sweeps.maxm() = 10,20,40,80,100,200,400,600,800,1000;
    sweeps.cutoff() = 1E-6;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    auto energy = dmrg(psi,H0,sweeps,"Quiet");

println("Ground State Energy = ",energy);


myfile << energy << "\t";

myfile << energy/N*4 << "\n";

for(int j = 1; j <= N; ++j)
    {
    //Make site j the MPS "orthogonality center"
    psi.position(j);
    //Measure magnetization
    Real Szj = (psi.A(j)
               * sites.op("Sz",j)
               * dag(prime(psi.A(j),Site))).real();
    println("Sz_",j," = ",Szj);

    auto Sxj = (psi.A(j)
               * sites.op("Sx",j)
               * dag(prime(psi.A(j),Site))).cplx();
    println("Sx",j," = ",Sxj);

    auto Syj = (psi.A(j)
               * sites.op("Sy",j)
               * dag(prime(psi.A(j),Site))).cplx();
    println("Sy",j," = ",Syj);
myfile << j << "\t";
myfile << Sxj << "\t";
myfile << Syj << "\t";
myfile << Szj << "\n";
    }


println("Ground State Energy = ",energy/N*4);

   return 0;
    }


