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

double hh;

int MSW;

int i1;

int j1;

cout << "Number of X unit cells Nx";
cin >> Nx;
cout << "The magnetic coupling Jx";
cin >> Jx;
cout << "The magnetic coupling Jy";
cin >> Jy;
cout << "The magnetic coupling Jz";
cin >> Jz;
cout << "The magnetic field along z direction";
cin >> hh;
cout << "The number of sweep MSW";
cin >> MSW;

// read in the parameters



ofstream myfile;
myfile.open ("Jx_"+std::to_string(Jx)+"_Jy_"+std::to_string(Jy)+"_Jz_"+std::to_string(Jz)+"_Nx_"+std::to_string(Nx)+"_hh_"+std::to_string(hh)+".dat");

ofstream myfile1;
myfile1.open ("Correlation_Jx_"+std::to_string(Jx)+"_Jy_"+std::to_string(Jy)+"_Jz_"+std::to_string(Jz)+"_Nx_"+std::to_string(Nx)+"_hh_"+std::to_string(hh)+".dat");


ofstream myfile2;
myfile2.open ("Magnetization_Jx_"+std::to_string(Jx)+"_Jy_"+std::to_string(Jy)+"_Jz_"+std::to_string(Jz)+"_Nx_"+std::to_string(Nx)+"_hh_"+std::to_string(hh)+".dat");


Jx=Jx*4;
Jy=Jy*4; // The spin operator here gives eigenvalues of +/- 1/2 rather than +/- 1
Jz=Jz*4;
hh=hh*2;

int N=2;

int N_=4*N*Nx;
        
int NDD;

int NDY=4*Nx-1;
int NDZ=8*Nx-2;

auto sites_ = SpinHalf(N_);

auto    H = MPO(sites_);



std::vector<Index> links(N_+1); // This is the vector of bond dimension


    for(int l = 0; l <= N_; ++l) 
        {
        links.at(l) = Index(nameint("hl",l),5+NDY+NDZ);
        }

    Index const& last = (links.at(0));

    for(int n = 1; n <= N_; ++n)
        {

        auto& W = H.Anc(n);
        auto row = dag(links.at(n-1));
        auto col = (n==N_ ? last : links.at(n));

        W = ITensor(dag(sites_(n)),prime(sites_(n)),row,col);

        W += sites_.op("Id",n) * row(1) * col(1); //ending state

        W += sites_.op("Sp",n) * row(2) * col(1)*0.5;	        	//X
        W += sites_.op("Sm",n) * row(2) * col(1)*0.5;	        	//X

        W += sites_.op("Sp",n) * row(3) * col(1)*0.5;			//Y
        W += -sites_.op("Sm",n) * row(3) * col(1)*0.5;			//Y

        W += sites_.op("Sp",n) * row(3+NDY) * col(1)*0.5;		//Y
        W +=  -sites_.op("Sm",n) * row(3+NDY) * col(1)*0.5;		//Y

        W += sites_.op("Sz",n) * row(4+NDY) * col(1);			//Z

        W += sites_.op("Id",n) * row(5+NDY+NDZ) * col(5+NDY+NDZ);	//starting state

        W += 2*hh*sites_.op("Sz",n) * row(5+NDY+NDZ) * col(1);	//Hubbard interaction

        for (int jj1 = 4; jj1 <=NDY+2; ++jj1)		{W += sites_.op("Id",n) * row(jj1) * col(jj1-1); }		//ok
        for (int jj2 = 5+NDY; jj2 <=NDY+NDZ+4; ++jj2)	{W += sites_.op("Id",n) * row(jj2) * col(jj2-1); }

        if (n%2==0 )	//ok
{       W += sites_.op("Sp",n) * row(5+NDY+NDZ) * col(2) * Jx *0.5;		//NN	
        W += sites_.op("Sm",n) * row(5+NDY+NDZ) * col(2) * Jx *0.5;		//NN
}

        if (n%2==1 && n%(4*Nx)!=1)	//ok
{       W += sites_.op("Sp",n) * row(5+NDY+NDZ) * col(3+NDY) * Jy *0.5;		//NN
        W += -sites_.op("Sm",n) * row(5+NDY+NDZ) * col(3+NDY) * Jy *0.5;		//NN
}

        if (n%(4*Nx)==2)	//ok
{       W += sites_.op("Sp",n) * row(5+NDY+NDZ) * col(2+NDY) * Jy *0.5;
        W += -sites_.op("Sm",n) * row(5+NDY+NDZ) * col(2+NDY) * Jy *0.5;    
}


	if (n%2==0 && (n<=4*Nx) && n>2)		{ NDD=N_-2*n+4;   W += Jz*sites_.op("Sz",n) * row(5+NDY+NDZ)* col(3+NDY+NDD);}  //ok

        if (n==2)				{ W += Jz*sites_.op("Sz",n) * row(5+NDY+NDZ)* col(3+NDY+4*Nx); }   //ok

        if (n==1)				{ W += Jz*sites_.op("Sz",n) * row(5+NDY+NDZ)* col(3+NDY+2); }   //ok

	if (n%2==1 && (n>4*Nx+1))		{ NDD=2*(N_-n)+4; W += Jz*sites_.op("Sz",n) *row(5+NDY+NDZ)* col(3+NDY+NDD);}	//ok


}


    auto sweeps = Sweeps(MSW);
    sweeps.maxm() = 40,80,120,200,300,400,600,800,1000,1200;
    sweeps.cutoff() = 2E-7;
    sweeps.niter() = 3,2;
    sweeps.noise() = 1E-7,1E-8,1E-9,0.0;

    auto state = InitState(sites_);
    for(int i = 1; i <= N_; ++i) 
        {
    state.set(i,"Up");
        }
    auto psi = MPS(state);


    auto LH = setElt(links.at(0)(5+NDY+NDZ));
    auto RH = setElt(dag(last)(1));

    H.Anc(0) = LH;
    H.Anc(N_+1) = RH;

auto res = idmrg(psi,H,sweeps,{"OutputLevel",1});

normalize(psi);

double energy=res.energy;



// create a hole



writeToFile(format("sites_file_Jx_%d_Jy_%d_Jz_%d_hh_%d_Nx_%df1",Jx,Jy,Jz,hh,Nx),sites_);

writeToFile(format("psi_Jx_%d_Jy_%d_Jz_%d_hh_%d_Nx_%df1",Jx,Jy,Jz,hh,Nx),psi);     //file name will be psi_100

myfile << energy << "\t";

myfile << energy/N_ << "\n";


psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation

for (int i=1;i<N_;++i)
{
for (int j=i+1;j<=N_;++j)
{

auto op_i = sites_.op("Sz",i);
auto op_j = sites_.op("Sz",j);

//below we will assume j > i

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
psi.position(i); 


//index linking i to i+1:
auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);

auto C = psi.A(i)*op_i*dag(prime(psi.A(i),Site,ir));
for(int k = i+1; k < j; ++k)
    {
    C *= psi.A(k);
    C *= dag(prime(psi.A(k),Link));
    }
C *= psi.A(j);
C *= op_j;
//index linking j to j-1:
auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);
C *= dag(prime(psi.A(j),jl,Site));

auto result = C.real(); //or C.cplx() if expecting complex

myfile1 << i <<"\t";
myfile1 << j <<"\t";
myfile1 << result <<"\n";

}
}

for(int j = 1; j <= N_; ++j)
    {
    //Make site j the MPS "orthogonality center"
    psi.position(j);
    //Measure magnetization
    Real Szj = (psi.A(j)
               * sites_.op("Sz",j)
               * dag(prime(psi.A(j),Site))).real();
    println("Sz_",j," = ",Szj);
myfile2 << j <<"\t";
myfile2 << Szj <<"\n";
    }

for(int j = 1; j <= N_; ++j)
    {
    //Make site j the MPS "orthogonality center"
    psi.position(j);
    //Measure magnetization
    Real Spj = (psi.A(j)
               * sites_.op("Sp",j)
               * dag(prime(psi.A(j),Site))).real();
    println("Sp_",j," = ",Spj);
myfile2 << j <<"\t";
myfile2 << Spj <<"\n";
    }

for(int j = 1; j <= N_; ++j)
    {
    //Make site j the MPS "orthogonality center"
    psi.position(j);
    //Measure magnetization
    Real Smj = (psi.A(j)
               * sites_.op("Sm",j)
               * dag(prime(psi.A(j),Site))).real();
    println("Sm_",j," = ",Smj);
myfile2 << j <<"\t";
myfile2 << Smj <<"\n";
    }

    return 0;
    }
