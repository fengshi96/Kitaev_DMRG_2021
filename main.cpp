#include <iostream>
#include "itensor/all.h"
#include "src/honeycomb.h"
#include "src/entanglement.h"

using namespace itensor;

/* check out
 * http://itensor.org/docs.cgi?page=formulas/2d_dmrg&vers=cppv3
 * https://www.itensor.org/support/2532/how-itensor-simulate-kitaev-honeycomb-model
 * https://github.com/ITensor/ITensor/blob/v3/itensor/mps/lattice/latticebond.h
 * */

int main(int argc, char* argv[]) {

    double Kz = std::stod(argv[1]);
    double field = std::stod(argv[2]);
    std::cout << "Kz = " << Kz << std::endl;
    std::cout << "Field = " << field << std::endl;

    int Nx = 3;
    int Ny = 3;
    int N = Nx*Ny*2;
    bool yperiodic = true;
    bool xperiodic = true;
    std::string cutting = "0,1,2,3";

    int cutLabel = 1;
    for (char const &c : cutting) if(c == ',') cutLabel++;
    std::cout << "cutLabel=" << cutLabel << std::endl;


    //  Initialize the site degrees of freedom.

    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    //
    // Use the AutoMPO feature to create the
    // next-neighbor Heisenberg model.
    //

    auto lattice = honeycombLattice(Nx,Ny,{"YPeriodic=",yperiodic,
                                                "XPeriodic=",xperiodic,
                                                "Cutting=",cutting});
    std::cout << lattice;

    auto ampo = AutoMPO(sites);
    // Isotropic Kitaev interaction
    for(auto bnd : lattice) {
        if (bnd.type == "Sz")
            ampo +=  Kz, bnd.type, bnd.s1, bnd.type, bnd.s2;
        else
            ampo += bnd.type, bnd.s1, bnd.type, bnd.s2;
    }


    // sx, 1, sx, 2

    // Add isotropic magnetic field
    for (int i = 1; i <=N; ++i) {
        ampo += field, "Sx", i;
        ampo += field, "Sy", i;
        ampo += field, "Sz", i;
    }

    auto H = toMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    // This choice implicitly sets the global Sz quantum number
    // of the wavefunction to zero. Since it is an MPS
    // it will remain in this quantum number sector.
    //
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i)
    {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
    }

    auto psi0 = MPS(state);

    //
    // overlap calculates matrix elements of MPO's with respect to MPS's
    // inner(psi0,H,psi0) = <psi0|H|psi0>
    //
    printfln("Initial energy = %.5f", innerC(psi0,H,psi0) );

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep.
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(12);
    sweeps.mindim() = 100,100,100,100,100,100,100,100,100,100,100,100;
    sweeps.maxdim() = 100,200,200,300,500,500,500,500,500,500,700,800;
    sweeps.cutoff() = 1E-8;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-5,1E-6,1E-7,1E-8;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing overlap = %.10f", innerC(psi,H,psi) );
    writeToFile(std::string("sites.dat"),sites);
    writeToFile(std::string("psi.dat"),psi);



    // Begin entanglment
    int bond = cutLabel;  // index of which bond to cut
    psi.position(bond);

    //SVD this wavefunction to get the spectrum of density-matrix eigenvalues
    auto l = leftLinkIndex(psi,bond);
    auto s = siteIndex(psi,bond);
    auto [U,S,V] = svd(psi(bond),{l,s});
    auto u = commonIndex(U,S);

    //Apply von Neumann formula
    //to the squares of the singular values
    Real SvN = 0.;
    for(auto n : range1(dim(u)))
    {
        auto Sn = elt(S,n,n);
        auto p = sqr(Sn);
        if(p > 1E-12) SvN += -p*log(p);
    }
    printfln("Across bond b=%d, SvN = %.10f",bond,SvN);


    return 0;
}
