#ifndef MODEL_H
#define MODEL_H
#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <tuple>
#include <map>
#include <malloc.h>
#include "common.h"
#include "Lattice.h"

using namespace std;

class spinModel
{

public:
    spinModel(Lattice *lattice, double theta, double Hmag);

    void build_localSpin();

    void build_interAll();
    void build_interAll(int oneD);

    void build_oneBodyG();

    void build_twoBodyG();

    void build_threeBodyG();

    void build_fourBodyG();

    void build_sixBodyG();

    void build_calcmod();

    void build_modpara();

    void build_namelist();

    void build_trans();

    void twositeInteraction(const int siteIndexa, const int siteIndexb);

    void onesiteInteraction(const int siteIndex);

    void transInteraction(const int siteIndex);

    void singleIonAnisotropy();

    ~spinModel();
    

private:
    Lattice *lattice;

    double theta;

    double h;

    int Nsites;

    int Nbonds;

    int Nplaqs;

    double SpinVal;

    int Sdim;

    double J, K, G1, G2, A; // 哈密顿量矩阵参数

    double Hmag, Htheta, Hpsi; // magnetic field strength, angles, components

    std::vector<std::tuple<int, int, int>> Bonds;

    std::vector<std::tuple<int, int, int, int, int, int>> Plaquettes;

    double Xbond[9];

    double Ybond[9];

    double Zbond[9];

    double bond[9];

    double SIA[9];

    double trans[3];

    my_cplx *interAllMatrix;

    my_cplx *transVec;

    int interactionInterAll = 0;
    int interactionTrans = 0;

    void clearInter();
    void clearTrans();

    void twoOpterm(string termType, double coeff);
    void oneOpterm(string termType, double coeff);

    void setBondInteraction();

    void setBondInteraction(string bondType);

    void setField(); // 磁场

    void interMatrix();

    void transMatrix();

    void ZZterm(const my_cplx coeff);

    void PPterm(const my_cplx coeff);

    void PMterm(const my_cplx coeff);

    void PZterm(const my_cplx coeff);

    void MPterm(const my_cplx coeff);

    void MMterm(const my_cplx coeff);

    void MZterm(const my_cplx coeff);

    void ZPterm(const my_cplx coeff);

    void ZMterm(const my_cplx coeff);

    void Pterm(const my_cplx coeff);
    void Mterm(const my_cplx coeff);
    void Zterm(const my_cplx coeff);

    double CGcoeff(double sigma1, string opa, double sigma2, string opb);
    double CGcoeff(double sigma, string op);

    void addToInterAll(string fileName, int indexa, double sigma1, double sigma2,
                       int indexb, double sigma3, double sigma4, double realpart, double imagpart);

    void addToTrans(string fileName, int const indexa, double const sigma1, int indexb, double sigma2, double realpart, double imagpart);
};

#endif
