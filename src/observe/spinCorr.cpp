//
// Created by shifeng on 3/3/22.
//
//#include <iostream>
//#include "itensor/all.h"
//#include "../honeycomb.h"
//using namespace itensor;

//int main(int argc, char* argv[]) {
//    if(argc < 2)  {
//        printfln("Usage: %s input_file",argv[0]);
//        return 0;
//    }
//    auto input = InputGroup(argv[1],"input");
//    auto orbitals = input.getInt("Orbitals");
//    auto LLX = input.getInt("LLX");
//    auto totalSweeps = input.getInt("totalSweeps");
//    auto N = LLX;
//    auto Norb = LLX*orbitals;
//    std::cout << std::setprecision(12);
//
//    double JExch = input.getReal("JExch");
//    double Lambda = input.getReal("Lambda");
//
//    // -- setting up sweeping paramaters --
//    auto sw_table = InputGroup(input,"table_name");
//    auto sweeps = Sweeps(totalSweeps,sw_table);
//    println(sweeps);
//    Args args;
//    args.add("Verbose",true);
//
//    Eigen::MatrixXi N1neigh_;
//    Eigen::VectorXi indx_,Nc_;
//    Chain1(LLX, orbitals, Nc_, N1neigh_, indx_);
//
//    Eigen::MatrixXd Connx(Norb,Norb), Conny(Norb,Norb), Connz(Norb,Norb);
//    Eigen::MatrixXd Sorx(Norb,Norb), Sory(Norb,Norb), Sorz(Norb,Norb);
//    Connx.setZero();
//
//    // --- Kitaev Model --- // Neighbors for each site
//    auto sites = SpinOne(Norb); //make a chain of N spin 1/2's
//    auto ampo = AutoMPO(sites);
//    std::cout << " Reading from File " << std::endl;
//    readFromFile("../sites.txt",sites);
//    auto psi = readFromFile<MPS>(std::string("../psi.txt"),sites);
//
//    // =========================================================
//    // =================== Observables =========================
//    // =========================================================
//    std::cout << " Starting calculations of Observables " << std::endl;
//    std::ofstream outfile;
//    outfile.open("Observables.dat");
//
//    Cplx SzTot=0, LzTot=0, JzTot=0;
//    println("\nj Lx Ly Lz = ");
//    for(int j=0; j < N; ++j) {
//        //re-gauge psi to get ready to measure at position j
//        //int ja = j*orbitals + 0 + 1;
//        int jb = j*orbitals + 1;
//
//        // ==== upper orbital - orbital index ======
//        psi.position(jb);
//        ITensor ket = psi.A(jb);
//        ITensor bra = dag(prime(ket,Site));
//
//        ITensor Lxjop = sites.op("Sx",jb); //*sites.op("Sx",j);
//        ITensor Lyjop = sites.op("Sy",jb); //*sites.op("Sy",j);
//        ITensor Lzjop = sites.op("Sz",jb); //*sites.op("Sz",j);
//
//        //take an inner product
//        auto Lxj = (bra*Lxjop*ket).cplx();
//        auto Lyj = (bra*Lyjop*ket).cplx();
//        auto Lzj = (bra*Lzjop*ket).cplx();
//
//
//        if((j>10) && (j<=N-10)) {
//            //SzTot += Cplx(szj); //.cplx();
//            LzTot += Cplx(Lzj);
//            //JzTot += Cplx(szj)+Cplx(Lzj);
//        }
//
//        printfln("%d %.12f %.12f %.12f",j,Lxj,Lyj,Lzj);
//        outfile << Lxj << "  " << Lyj << "  " << Lzj << "  \n";
//
//    }
//    outfile << std::endl;
//    outfile.close();

//    return 0;
//}

