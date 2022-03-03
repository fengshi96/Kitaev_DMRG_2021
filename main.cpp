#include <iostream>
#include "itensor/all.h"
#include "src/honeycomb.h"
#include "src/entanglement.h"
#include <Eigen/Dense>

using namespace itensor;

/* check out
 * http://itensor.org/docs.cgi?page=formulas/2d_dmrg&vers=cppv3
 * https://www.itensor.org/support/2532/how-itensor-simulate-kitaev-honeycomb-model
 * https://github.com/ITensor/ITensor/blob/v3/itensor/mps/lattice/latticebond.h
 * */

int main(int argc, char* argv[]) {
    if(argc < 2)  {
        printfln("Usage: %s input_file",argv[0]);
        return 0;
    }
    auto input = InputGroup(argv[1],"input");
    int Nx = input.getInt("Nx");
    int Ny = input.getInt("Ny");
    int N = Nx*Ny*2;
    bool xperiodic = input.getInt("IsPeriodicX");
    bool yperiodic = input.getInt("IsPeriodicY");

    double Kx = input.getInt("Kx");
    double Ky = input.getInt("Ky");
    double Kz = input.getInt("Kz");

    double Hx = input.getInt("Hx");
    double Hy = input.getInt("Hy");
    double Hz = input.getInt("Hz");

    // sweep parameters
    auto totalSweeps = input.getInt("totalSweeps");
    auto sw_table = InputGroup(input,"sweep_table");
    auto sweeps = Sweeps(totalSweeps,sw_table);
    println(sweeps);

//    std::string cutting = "0,1,2,3";
//    int cutLabel = 1;
//    for (char const &c : cutting) if(c == ',') cutLabel++;
//    std::cout << "cutLabel=" << cutLabel << std::endl;


    //  Initialize the site degrees of freedom.

    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    // auto lattice = honeycombLattice(Nx,Ny,{"YPeriodic=",yperiodic, "XPeriodic=",xperiodic, "Cutting=",cutting});
    auto lattice = honeycombLattice(Nx,Ny,{"YPeriodic=",yperiodic, "XPeriodic=",xperiodic});
    std::cout << lattice;

    auto ampo = AutoMPO(sites);
    // Kitaev interaction
    for(auto bnd : lattice) {
        if (bnd.type == "Sx")
            ampo +=  Kz, bnd.type, bnd.s1, bnd.type, bnd.s2;
        else if (bnd.type == "Sy")
            ampo += Ky, bnd.type, bnd.s1, bnd.type, bnd.s2;
        else if (bnd.type == "Sz")
            ampo += Kz, bnd.type, bnd.s1, bnd.type, bnd.s2;
    }


    // Add magnetic field
    for (int i = 1; i <=N; ++i) {
        ampo += Hx, "Sx", i;
        ampo += Hy, "Sy", i;
        ampo += Hz, "Sz", i;
    }

    auto H = toMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i)
    {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
    }

    auto psi0 = MPS(state);

    // overlap calculates matrix elements of MPO's with respect to MPS's
    // inner(psi0,H,psi0) = <psi0|H|psi0>
    printfln("Initial energy = %.5f", innerC(psi0,H,psi0) );

    // Begin the DMRG calculation
    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing overlap = %.10f", innerC(psi,H,psi) );
    writeToFile(std::string("sites.dat"),sites);
    writeToFile(std::string("psi.dat"),psi);



//    // Begin entanglment
//    int bond = cutLabel;  // index of which bond to cut
//    psi.position(bond);
//
//    //SVD this wavefunction to get the spectrum of density-matrix eigenvalues
//    auto l = leftLinkIndex(psi,bond);
//    auto s = siteIndex(psi,bond);
//    auto [U,S,V] = svd(psi(bond),{l,s});
//    auto u = commonIndex(U,S);
//
//    //Apply von Neumann formula
//    //to the squares of the singular values
//    Real SvN = 0.;
//    for(auto n : range1(dim(u)))
//    {
//        auto Sn = elt(S,n,n);
//        auto p = sqr(Sn);
//        if(p > 1E-12) SvN += -p*log(p);
//    }
//    printfln("Across bond b=%d, SvN = %.10f",bond,SvN);


    return 0;
}
