#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <tuple>
#include <map>
#include <malloc.h>
#include "common.h"

using namespace std;

class Lattice
{
private:
    int cluster_type;
    int Nsites;

    int Nbonds;
    int Nplaqs;

    std::string cluster_name;

    std::vector<std::tuple<int, int, int>> Bonds;

    std::vector<std::tuple<int, int, int, int, int, int>> Plaquettes;

    double LatticeA = 1.0;
    double *a1 = new double[3]{LatticeA, 0.0, 0.0};
    double *a2 = new double[3]{LatticeA / 2.0, sqrt(3.0) * LatticeA / 2.0, 0.0};
    double *a3 = new double[3]{0.0, 0.0, LatticeA};

    double *b1 = new double[3];
    double *b2 = new double[3];

public:
    Lattice(int cluster_type);
    void BuildCluster(int cluster_type);
    void PrintClusterDetails();

    //*********************** */
    int get_Nsites(void) { return Nsites; }

    int get_Nbonds(void) { return Nbonds; }

    int get_Nplaqs(void) { return Nplaqs; }

    std::array<double, 3> get_a1() const { return {a1[0], a1[1], a1[2]}; }
    std::array<double, 3> get_a2() const { return {a2[0], a2[1], a2[2]}; }
    std::array<double, 3> get_b1() const { return {b1[0], b1[1], b1[2]}; }
    std::array<double, 3> get_b2() const { return {b2[0], b2[1], b2[2]}; }

    std::vector<std::tuple<int, int, int>> get_Bonds(void) { return Bonds; }
    std::vector<std::tuple<int, int, int, int, int, int>> get_Plaquettes(void) { return Plaquettes; }

    void Get_Rspace_point(int cluster_type);
    void Cal_b();
    void Cal_BZ();
    void Map_sites(int old_value, int new_value, const std::vector<std::tuple<int, double, double>> &old_sites, std::vector<std::tuple<int, double, double>> &new_sites);

    void BuildMomLine();
    void Output_Mpoint_to_file(vector<pair<double, double>> &Discrete_Mpoints, const string &Filename);
    void Discrete_DSSF_line(std::pair<double, double> a, std::pair<double, double> b, int n, int per, std::vector<std::pair<double, double>> &Discrete_line);
    void BuildMomSpace();

    //
    ~Lattice();
};
