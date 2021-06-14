//
// Created by shifeng on 6/13/21.
//
/* http://itensor.org/docs.cgi?page=formulas/entanglement_mps = entanglement of an MPS
 * http://itensor.org/support/580/entanglement-entropy-of-2d-mps?show=580#q580 = advice on 2D entanglement
 * http://itensor.org/docs.cgi?page=classes/bondgate&vers=cppv3 = Bond gate
 * http://itensor.org/docs.cgi?vers=cppv3&page=formulas/tevol_trotter
 * https://iopscience-iop-org.proxy.lib.ohio-state.edu/article/10.1088/1367-2630/12/5/055026 = swap gates
 */

#ifndef ITENSOR_2D_ENTANGLEMENT_H
#define ITENSOR_2D_ENTANGLEMENT_H


namespace itensor {

    std::vector<BondGate> swapGateList(const std::vector<size_t>& sysIndx, const std::vector<LatticeBond>& lattice) {
        int Nsys = sysIndx.size();  // number of sites in the sub-system

        int indxCounter = 0;
        for (size_t indx : sysIndx) {

        }

    }

}  //namespace itensor

#endif //ITENSOR_2D_ENTANGLEMENT_H
