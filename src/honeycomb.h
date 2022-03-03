//
// Created by shifeng on 6/12/21.
//

#ifndef ITENSOR_2D_HONEYCOMB_H
#define ITENSOR_2D_HONEYCOMB_H

#include <iostream>
#include <string>
#include "matrix.h"


std::vector<size_t> labelParser(const std::string& cutting = "", const std::string& delimiter = ",") {
    size_t pos = 0;
    std::string token;

    std::vector<size_t> sysIndx;
    sysIndx.reserve(12);

    std::string str2Cut = cutting;
    while ((pos = str2Cut.find(delimiter)) != std::string::npos) {
        token = str2Cut.substr(0, pos);
        sysIndx.emplace_back(std::stoi(token, nullptr, 10));
        str2Cut.erase(0, pos + delimiter.length());
    }
    sysIndx.emplace_back(std::stoi(str2Cut, nullptr, 10));
    for (auto const &l : sysIndx) {std::cout << "sysIndx:" << l << std::endl;}

    return sysIndx;
}


alg::Matrix<int> Honeycomb(size_t Nx, size_t Ny, bool xPBC = true, bool yPBC = true,
                           const std::string& cutting = "") {
    std::vector<size_t> sysIndx;
    if (cutting != ""){
        sysIndx = labelParser(cutting);
    }

    int Number1Neigh_ = 3;
    int Nsite_ = Nx * Ny * 2;

    alg::Matrix<int> Nc_, N1neigh_(Nsite_, Number1Neigh_);
    std::vector<int> indx_, indy_;


    int LLX = Nx;
    int LLY = Ny;
    std::cout << "creating Honeycomb ... " << std::endl;
    double scalex = 2.0, scaley = 4.0 / sqrt(3.0);
    std::vector<double> t1(2), t2(2);
    t1[0] = 1.0 * scalex;
    t1[1] = 0;
    t2[0] = 0.5 * scalex;
    t2[1] = sqrt(3.0) / 2.0 * scaley;

    // Site labeling
    indx_.resize(Nsite_);
    indy_.resize(Nsite_);
    Nc_.resize(LLX * 2 + LLY, LLY * 2);
    Nc_.fill(-1);

    int xv = 0, counter = 0;
    for (int i = 0; i < LLX; i++) {
        if (i != 0) { xv += t1[0]; }
        int y0 = 0;
        int x1, y1;
        x1 = xv + 1.0;
        y1 = 1.0;


        for (int j = 0; j < LLY; j++) {
            int cxa = int(xv + j * t2[0]);
            int cxb = int(x1 + j * t2[0]);
            int cya = int(y0 + j * t2[1]);
            int cyb = int(y1 + j * t2[1]);

            indx_[counter] = cxa;
            indy_[counter] = cya;
            Nc_(cxa, cya) = counter;
            counter++;

            indx_[counter] = cxb;
            indy_[counter] = cyb;
            Nc_(cxb, cyb) = counter;
            counter++;
        }
    }

    if (sysIndx.size() > 0) {

        int newCounter = 0;
        for (auto site : sysIndx) {

            int& x = indx_[newCounter];
            int& y = indy_[newCounter];
            int& ix = indx_[site];
            int& iy = indy_[site];

            // permutate in Nc_
            int tmp = Nc_(ix, iy);
            Nc_(ix, iy) = newCounter;
            Nc_(x, y) = tmp;

            // permutate in indx_ and indy_
            int xtmp = x;
            int ytmp = y;
            x = ix;
            y = iy;
            ix = xtmp;
            iy = ytmp;

            newCounter ++;
        }
    }
//    for (auto i : indx_) std::cout << i << " ";
//    std::cout <<std::endl;
//    for (auto i : indy_) std::cout << i << " ";
//    std::cout << std::endl;


    int xmax = *std::max_element(indx_.begin(), indx_.end());
    int ymax = *std::max_element(indy_.begin(), indy_.end());

    int jx, jy;
    for (int i = 0; i < Nsite_; i++) {    // ith site
        int ix = indx_[i];
        int iy = indy_[i];

        // SxSx (kitaev) - neighbor 0
        jx = ix + 1;
        jy = iy + 1;
        if (jx <= xmax && jy <= ymax && Nc_(jx, jy) != -1) {
            int j = Nc_(jx, jy);
            N1neigh_(i, 2) = j;
            N1neigh_(j, 2) = i;
        }

        // SySy (kitaev) - neighbor 1
        jx = ix + 1;
        jy = iy - 1;
        if(jx<=xmax && jy<=ymax && jy>=0 && Nc_(jx,jy)!=-1)  {
            int j = Nc_(jx,jy);
            N1neigh_(i,1) = j;
            N1neigh_(j,1) = i;
        }

        // SzSz (kitaev) - neighbor 2
        jx = ix;
        jy = iy + 1;
        if(jx<=xmax && jy<=ymax && Nc_(jx,jy)!=-1)  {
            int j = Nc_(jx,jy);
            N1neigh_(i,0) = j;
            N1neigh_(j,0) = i;
        }



        // Apply PBC
        if(yPBC) {
            jx = int(ix - LLY);
            jy = 0;
            if(jx>=0 && iy==ymax && Nc_(jx,jy)!=-1) {
                int j = Nc_(jx,jy);
                N1neigh_(i,0) = j;
                N1neigh_(j,0) = i;
            }
        }

        if(xPBC) {
            jx = int(ix + LLX*2 - 1);
            jy = int(iy + 1);
            if (jx<=xmax && iy<=ymax && iy % 2 == 0 && Nc_(jx,jy)!=-1) {
                int j = Nc_(jx,jy);
                N1neigh_(i,1) = j;
                N1neigh_(j,1) = i;
            }
        }
        //        std::cout << "(ix,iy)=(" << ix << "," << iy << "), "
//                  << "(jx,jy)=(" << jx << "," << jy << "), "
//                  << "j=" << Nc_(jx,jy) << " i=" << i << std::endl;
    }

    // draw grid, start from 1
    for (int i = 0; i < (LLX * 2 + LLY) * LLY * 2; ++i)
        Nc_[i]++;
    Nc_.print();
    std::cout << std::endl;

    return N1neigh_;
} // end Honeycomb



namespace itensor {
    // using LatticeGraph = std::vector<LatticeBond>;
    LatticeGraph inline honeycombLattice(int Nx, int Ny, Args const& args = Args::global()) {

        int N = Nx * Ny * 2;

        auto xperiodic = args.getBool("XPeriodic", true);
        auto yperiodic = args.getBool("YPeriodic", true);
        auto cutting = args.getString("Cutting", "");
        // Periodicity on x/y is meaningless for one dimensional chain
        xperiodic = xperiodic && (Ny > 1);
        yperiodic = yperiodic && (Ny > 1);

        alg::Matrix<int> N1neigh = Honeycomb(Nx, Ny, xperiodic, yperiodic, cutting);

        LatticeGraph latt;
        latt.reserve(int(N * 1.5) + 1);

        // z bond
        for (int i = 0; i < N; ++i) {
            int site1 = i + 1;
            int site2 = N1neigh(i, 0) + 1;  // + 1 since itensor counts from 1
            if (site2 != -1 && site2 < site1) {
                //std::cout << "Push z-bond  (" << site1 << "," << site2 << ")" << std::endl;
                latt.emplace_back(site1, site2, "Sz");
            }
        }

        // y bond
        for (int i = 0; i < N; ++i) {
            int site1 = i + 1;
            int site2 = N1neigh(i, 1) + 1;  // + 1 since itensor counts from 1
            if (site2 != -1 && site2 < site1) {
                latt.emplace_back(site1, site2, "Sy");
            }
        }

        // x bond
        for (int i = 0; i < N; ++i) {
            int site1 = i + 1;
            int site2 = N1neigh(i, 2) + 1;  // + 1 since itensor counts from 1
            if (site2 != -1 && site2 < site1) {
                latt.emplace_back(site1, site2, "Sx");
            }
        }

        for (int i = 0; i < N1neigh.rows() * N1neigh.cols(); ++i)
            N1neigh[i]++;
        N1neigh.print();
        std::cout << std::endl;

        return latt;
    }
}  //namespace itensor










#endif //ITENSOR_2D_HONEYCOMB_H
