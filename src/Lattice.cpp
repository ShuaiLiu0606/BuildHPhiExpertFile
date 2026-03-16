#include <Model.h>
#include <common.h>

Lattice::~Lattice()
{

    delete[] a1;
    delete[] a2;
    delete[] a3;

    delete[] b1;
    delete[] b2;
}

Lattice::Lattice(int cluster_type) : cluster_type(cluster_type)
{

    BuildCluster(cluster_type);
    PrintClusterDetails();
    Cal_b();
    Get_Rspace_point(cluster_type);
    Cal_BZ();

    BuildMomLine();
    BuildMomSpace();
}

void Lattice::PrintClusterDetails()
{
    std::cout << std::endl
              << "Cluster: " << cluster_name << std::endl
              << std::endl;
    std::cout << "# of cells/plaquettes:      " << Nplaqs << std::endl;
    std::cout << "# of sites:                 " << Nsites << std::endl;
    std::cout << "# of unique bonds:          " << Nbonds << std::endl;

    std::cout << std::endl;
    std::cout << "Bonds:" << std::endl;
    std::cout << "x:";

    ofstream my_fiel("check_bond.dat");

    for (int b = 0; b < Nbonds; b++)
    {
        const auto &bond = Bonds[b];
        if (std::get<2>(bond) == 0) // Type 0
        {
            std::cout << " (" << std::get<0>(bond) << "," << std::get<1>(bond) << ")";

            my_fiel << "0" << setw(4) << right << std::get<0>(bond) << setw(5) << right << std::get<1>(bond) << endl;
        }
    }
    std::cout << std::endl;

    std::cout << "y:";
    for (int b = 0; b < Nbonds; b++)
    {
        const auto &bond = Bonds[b];
        if (std::get<2>(bond) == 1) // Type 1
        {
            std::cout << " (" << std::get<0>(bond) << "," << std::get<1>(bond) << ")";

            my_fiel << "1" << setw(4) << right << std::get<0>(bond) << setw(5) << right << std::get<1>(bond) << endl;
        }
    }
    std::cout << std::endl;

    std::cout << "z:";
    for (int b = 0; b < Nbonds; b++)
    {
        const auto &bond = Bonds[b];
        if (std::get<2>(bond) == 2) // Type 2
        {
            std::cout << " (" << std::get<0>(bond) << "," << std::get<1>(bond) << ")";

            my_fiel << "2" << setw(4) << right << std::get<0>(bond) << setw(5) << right << std::get<1>(bond) << endl;
        }
    }
    std::cout << std::endl
              << std::endl;
}

void Lattice::BuildCluster(int cluster_type)
{
    switch (cluster_type)
    {
    case 0: // 6-cluster
    {
        cluster_name = "6sites";

        Nsites = 6;

        Bonds.push_back(std::make_tuple(1, 2, 0));
        Bonds.push_back(std::make_tuple(5, 4, 0));
        Bonds.push_back(std::make_tuple(3, 0, 0));

        Bonds.push_back(std::make_tuple(0, 5, 1));
        Bonds.push_back(std::make_tuple(2, 3, 1));
        Bonds.push_back(std::make_tuple(1, 4, 1));

        Bonds.push_back(std::make_tuple(0, 1, 2));
        Bonds.push_back(std::make_tuple(3, 4, 2));
        Bonds.push_back(std::make_tuple(2, 5, 2));

        Nbonds = Bonds.size();

        Plaquettes.push_back(std::make_tuple(0, 1, 2, 3, 4, 5));
        Plaquettes.push_back(std::make_tuple(2, 5, 4, 1, 0, 3));
        Plaquettes.push_back(std::make_tuple(4, 3, 0, 5, 2, 1));

        Nplaqs = Plaquettes.size();
    }
    break;
    case 1: // 8-C3 cluster 2 * 2 * 2
    {
        cluster_name = "2-2-2sites";

        Nsites = 8;

        Bonds.push_back(std::make_tuple(0, 1, 0));
        Bonds.push_back(std::make_tuple(2, 3, 0));
        Bonds.push_back(std::make_tuple(4, 5, 0));
        Bonds.push_back(std::make_tuple(6, 7, 0));

        Bonds.push_back(std::make_tuple(1, 2, 1));
        Bonds.push_back(std::make_tuple(3, 0, 1));
        Bonds.push_back(std::make_tuple(5, 6, 1));
        Bonds.push_back(std::make_tuple(7, 4, 1));

        Bonds.push_back(std::make_tuple(1, 4, 2));
        Bonds.push_back(std::make_tuple(3, 6, 2));
        Bonds.push_back(std::make_tuple(5, 0, 2));
        Bonds.push_back(std::make_tuple(7, 2, 2));

        Nbonds = Bonds.size();

        Plaquettes.push_back(std::make_tuple(2, 7, 6, 5, 0, 1));
        Plaquettes.push_back(std::make_tuple(0, 5, 4, 7, 2, 3));
        Plaquettes.push_back(std::make_tuple(6, 3, 2, 1, 4, 5));
        Plaquettes.push_back(std::make_tuple(4, 1, 0, 3, 6, 7));

        Nplaqs = Plaquettes.size();
    }
    break;
    case 2: // 12-C3 cluster 3 * 2 * 2
    {
        cluster_name = "3-2-2sites";
        Nsites = 12;

        Bonds.push_back(std::make_tuple(0, 1, 0));
        Bonds.push_back(std::make_tuple(2, 3, 0));
        Bonds.push_back(std::make_tuple(4, 5, 0));
        Bonds.push_back(std::make_tuple(6, 7, 0));
        Bonds.push_back(std::make_tuple(8, 9, 0));
        Bonds.push_back(std::make_tuple(10, 11, 0));

        Bonds.push_back(std::make_tuple(1, 2, 1));
        Bonds.push_back(std::make_tuple(3, 4, 1));
        Bonds.push_back(std::make_tuple(5, 0, 1));
        Bonds.push_back(std::make_tuple(7, 8, 1));
        Bonds.push_back(std::make_tuple(9, 10, 1));
        Bonds.push_back(std::make_tuple(11, 6, 1));

        Bonds.push_back(std::make_tuple(1, 6, 2));
        Bonds.push_back(std::make_tuple(3, 8, 2));
        Bonds.push_back(std::make_tuple(5, 10, 2));
        Bonds.push_back(std::make_tuple(7, 0, 2));
        Bonds.push_back(std::make_tuple(9, 2, 2));
        Bonds.push_back(std::make_tuple(11, 4, 2));

        Nbonds = Bonds.size();

        Plaquettes.push_back(std::make_tuple(2, 9, 8, 7, 0, 1));
        Plaquettes.push_back(std::make_tuple(4, 11, 10, 9, 2, 3));
        Plaquettes.push_back(std::make_tuple(0, 7, 6, 11, 4, 5));
        Plaquettes.push_back(std::make_tuple(8, 3, 2, 1, 6, 7));
        Plaquettes.push_back(std::make_tuple(10, 5, 4, 3, 8, 9));
        Plaquettes.push_back(std::make_tuple(6, 1, 0, 5, 10, 11));

        Nplaqs = Plaquettes.size();
    }
    break;
    case 3: // 18-cluster 3 * 3 * 2
    {
        cluster_name = "3-3-2sites";

        Nsites = 18;

        Bonds.push_back(std::make_tuple(0, 1, 0));
        Bonds.push_back(std::make_tuple(2, 3, 0));
        Bonds.push_back(std::make_tuple(4, 5, 0));
        Bonds.push_back(std::make_tuple(6, 7, 0));
        Bonds.push_back(std::make_tuple(8, 9, 0));
        Bonds.push_back(std::make_tuple(10, 11, 0));
        Bonds.push_back(std::make_tuple(12, 13, 0));
        Bonds.push_back(std::make_tuple(14, 15, 0));
        Bonds.push_back(std::make_tuple(16, 17, 0));

        Bonds.push_back(std::make_tuple(1, 2, 1));
        Bonds.push_back(std::make_tuple(3, 4, 1));
        Bonds.push_back(std::make_tuple(5, 0, 1));
        Bonds.push_back(std::make_tuple(7, 8, 1));
        Bonds.push_back(std::make_tuple(9, 10, 1));
        Bonds.push_back(std::make_tuple(11, 6, 1));
        Bonds.push_back(std::make_tuple(13, 14, 1));
        Bonds.push_back(std::make_tuple(15, 16, 1));
        Bonds.push_back(std::make_tuple(17, 12, 1));

        Bonds.push_back(std::make_tuple(1, 6, 2));
        Bonds.push_back(std::make_tuple(3, 8, 2));
        Bonds.push_back(std::make_tuple(5, 10, 2));
        Bonds.push_back(std::make_tuple(7, 12, 2));
        Bonds.push_back(std::make_tuple(9, 14, 2));
        Bonds.push_back(std::make_tuple(11, 16, 2));
        Bonds.push_back(std::make_tuple(13, 0, 2));
        Bonds.push_back(std::make_tuple(15, 2, 2));
        Bonds.push_back(std::make_tuple(17, 4, 2));

        Nbonds = Bonds.size();

        Plaquettes.push_back(std::make_tuple(1, 6, 7, 8, 3, 2));
        Plaquettes.push_back(std::make_tuple(3, 8, 9, 10, 5, 4));
        Plaquettes.push_back(std::make_tuple(5, 10, 11, 6, 1, 0));
        Plaquettes.push_back(std::make_tuple(7, 12, 13, 14, 9, 8));
        Plaquettes.push_back(std::make_tuple(9, 14, 15, 16, 11, 10));
        Plaquettes.push_back(std::make_tuple(11, 16, 17, 12, 7, 6));
        Plaquettes.push_back(std::make_tuple(13, 0, 1, 2, 15, 14));
        Plaquettes.push_back(std::make_tuple(15, 2, 3, 4, 17, 16));

        Nplaqs = Plaquettes.size();
    }
    break;
    case 4: // 24-C6 cluster
    {
        cluster_name = "24-C6";

        Nsites = 24;

        Bonds.push_back(std::make_tuple(0, 1, 0));
        Bonds.push_back(std::make_tuple(2, 3, 0));
        Bonds.push_back(std::make_tuple(4, 5, 0));
        Bonds.push_back(std::make_tuple(6, 7, 0));
        Bonds.push_back(std::make_tuple(8, 9, 0));
        Bonds.push_back(std::make_tuple(10, 11, 0));
        Bonds.push_back(std::make_tuple(14, 13, 0));
        Bonds.push_back(std::make_tuple(16, 15, 0));
        Bonds.push_back(std::make_tuple(18, 17, 0));
        Bonds.push_back(std::make_tuple(20, 19, 0));
        Bonds.push_back(std::make_tuple(22, 21, 0));
        Bonds.push_back(std::make_tuple(23, 12, 0));

        Bonds.push_back(std::make_tuple(1, 2, 1));
        Bonds.push_back(std::make_tuple(3, 4, 1));
        Bonds.push_back(std::make_tuple(5, 6, 1));
        Bonds.push_back(std::make_tuple(7, 8, 1));
        Bonds.push_back(std::make_tuple(9, 10, 1));
        Bonds.push_back(std::make_tuple(11, 0, 1));
        Bonds.push_back(std::make_tuple(13, 12, 1));
        Bonds.push_back(std::make_tuple(15, 14, 1));
        Bonds.push_back(std::make_tuple(17, 16, 1));
        Bonds.push_back(std::make_tuple(19, 18, 1));
        Bonds.push_back(std::make_tuple(21, 20, 1));
        Bonds.push_back(std::make_tuple(23, 22, 1));

        Bonds.push_back(std::make_tuple(0, 17, 2));
        Bonds.push_back(std::make_tuple(2, 15, 2));
        Bonds.push_back(std::make_tuple(4, 13, 2));
        Bonds.push_back(std::make_tuple(6, 23, 2));
        Bonds.push_back(std::make_tuple(8, 21, 2));
        Bonds.push_back(std::make_tuple(10, 19, 2));
        Bonds.push_back(std::make_tuple(12, 11, 2));
        Bonds.push_back(std::make_tuple(14, 9, 2));
        Bonds.push_back(std::make_tuple(16, 7, 2));
        Bonds.push_back(std::make_tuple(18, 5, 2));
        Bonds.push_back(std::make_tuple(20, 3, 2));
        Bonds.push_back(std::make_tuple(22, 1, 2));

        Nbonds = Bonds.size();

        Plaquettes.push_back(std::make_tuple(22, 1, 0, 11, 12, 23));
        Plaquettes.push_back(std::make_tuple(20, 3, 2, 1, 22, 21));
        Plaquettes.push_back(std::make_tuple(18, 5, 4, 3, 20, 19));
        Plaquettes.push_back(std::make_tuple(6, 23, 12, 13, 4, 5));
        Plaquettes.push_back(std::make_tuple(8, 21, 22, 23, 6, 7));
        Plaquettes.push_back(std::make_tuple(10, 19, 20, 21, 8, 9));
        Plaquettes.push_back(std::make_tuple(12, 11, 10, 9, 14, 13));
        Plaquettes.push_back(std::make_tuple(4, 13, 14, 15, 2, 3));
        Plaquettes.push_back(std::make_tuple(14, 9, 8, 7, 16, 15));
        Plaquettes.push_back(std::make_tuple(2, 15, 16, 17, 0, 1));
        Plaquettes.push_back(std::make_tuple(16, 7, 6, 5, 18, 17));
        Plaquettes.push_back(std::make_tuple(0, 17, 18, 19, 10, 11));

        Nplaqs = Plaquettes.size();
    }
    break;
    case 5: //
    {
        cluster_name = "4-3-2sites";

        Nsites = 24;

        Bonds.push_back(std::make_tuple(0, 1, 0));
        Bonds.push_back(std::make_tuple(2, 3, 0));
        Bonds.push_back(std::make_tuple(4, 5, 0));
        Bonds.push_back(std::make_tuple(6, 7, 0));
        Bonds.push_back(std::make_tuple(8, 9, 0));
        Bonds.push_back(std::make_tuple(10, 11, 0));
        Bonds.push_back(std::make_tuple(12, 13, 0));
        Bonds.push_back(std::make_tuple(14, 15, 0));
        Bonds.push_back(std::make_tuple(16, 17, 0));
        Bonds.push_back(std::make_tuple(18, 19, 0));
        Bonds.push_back(std::make_tuple(20, 21, 0));
        Bonds.push_back(std::make_tuple(22, 23, 0));

        Bonds.push_back(std::make_tuple(1, 2, 1));
        Bonds.push_back(std::make_tuple(3, 4, 1));
        Bonds.push_back(std::make_tuple(5, 6, 1));
        Bonds.push_back(std::make_tuple(7, 0, 1));
        Bonds.push_back(std::make_tuple(9, 10, 1));
        Bonds.push_back(std::make_tuple(11, 12, 1));
        Bonds.push_back(std::make_tuple(13, 14, 1));
        Bonds.push_back(std::make_tuple(15, 8, 1));
        Bonds.push_back(std::make_tuple(17, 18, 1));
        Bonds.push_back(std::make_tuple(19, 20, 1));
        Bonds.push_back(std::make_tuple(21, 22, 1));
        Bonds.push_back(std::make_tuple(23, 16, 1));

        Bonds.push_back(std::make_tuple(1, 8, 2));
        Bonds.push_back(std::make_tuple(3, 10, 2));
        Bonds.push_back(std::make_tuple(5, 12, 2));
        Bonds.push_back(std::make_tuple(7, 14, 2));
        Bonds.push_back(std::make_tuple(9, 16, 2));
        Bonds.push_back(std::make_tuple(11, 18, 2));
        Bonds.push_back(std::make_tuple(13, 20, 2));
        Bonds.push_back(std::make_tuple(15, 22, 2));
        Bonds.push_back(std::make_tuple(17, 0, 2));
        Bonds.push_back(std::make_tuple(19, 2, 2));
        Bonds.push_back(std::make_tuple(21, 4, 2));
        Bonds.push_back(std::make_tuple(23, 6, 2));

        Nbonds = Bonds.size();

        Plaquettes.push_back(std::make_tuple(2, 19, 18, 17, 0, 1));
        Plaquettes.push_back(std::make_tuple(4, 21, 20, 19, 2, 3));
        Plaquettes.push_back(std::make_tuple(6, 23, 22, 21, 4, 5));
        Plaquettes.push_back(std::make_tuple(0, 17, 16, 23, 6, 7));
        Plaquettes.push_back(std::make_tuple(10, 3, 2, 1, 8, 9));
        Plaquettes.push_back(std::make_tuple(12, 5, 4, 3, 10, 11));
        Plaquettes.push_back(std::make_tuple(14, 7, 6, 5, 12, 13));
        Plaquettes.push_back(std::make_tuple(8, 1, 0, 7, 14, 15));
        Plaquettes.push_back(std::make_tuple(18, 11, 10, 9, 16, 17));
        Plaquettes.push_back(std::make_tuple(20, 13, 12, 11, 18, 19));
        Plaquettes.push_back(std::make_tuple(22, 15, 14, 13, 20, 21));
        Plaquettes.push_back(std::make_tuple(16, 9, 8, 15, 22, 23));

        Nplaqs = Plaquettes.size();
    }
    break;
    }
}

void Lattice::Get_Rspace_point(int cluster_type)
{

    int XX = 4; // 列数
    int YY = 4; // 行数

    std::vector<std::tuple<int, double, double>> old_sites_location;
    // 子格B相对于子格A的偏移量
    double *delta = new double[2]{a1[0] / 3.0 + a2[0] / 3.0, a1[1] / 3.0 + a2[1] / 3.0};

    // 长方形晶格，格点对高度求商得到横坐标(列数)，高度取余得到纵坐标(行数数)//
    for (int i = 0; i < XX; i++)
    {
        for (int j = 0; j < (2 * YY); j++)
        {
            int sites = i * (2 * YY) + j;

            int N1, N2;
            double *temp_sites_location = new double[2];
            if (sites < (XX * YY))
            {
                N1 = sites / (YY);
                N2 = sites % (YY);

                temp_sites_location[0] = N1 * a1[0] + N2 * a2[0];
                temp_sites_location[1] = N1 * a1[1] + N2 * a2[1];
            }
            else
            {
                int new_index = sites - (XX * YY);
                N1 = new_index / (YY);
                N2 = new_index % (YY);

                temp_sites_location[0] = N1 * a1[0] + N2 * a2[0] + delta[0];
                temp_sites_location[1] = N1 * a1[1] + N2 * a2[1] + delta[1];
            }
            old_sites_location.push_back(std::make_tuple(sites, temp_sites_location[0], temp_sites_location[1]));

            delete[] temp_sites_location;
        }
    }

    delete[] delta;

    std::vector<std::tuple<int, double, double>> new_sites_location;
    for (int i = 0; i < Nsites; i++)
    {
        new_sites_location.emplace_back(i, 0.0, 0.0);
    }

    std::vector<int> new_sites;
    for (int i = 0; i < Nsites; i++)
    {
        new_sites.push_back(i);
    }

    // 12sites 2*3
    // vector<int> old_sites = {0, 16, 4, 20, 1, 17, 5, 21, 2, 18, 6, 22};

    std::vector<int> old_sites(Nsites);

    if (cluster_type == 0)
    {
        // 6sites
        old_sites = {16, 1, 17, 5, 20, 4};
    }
    else if (cluster_type == 1)
    {
        // 8sites 2*2
        old_sites = {0, 16, 4, 20, 1, 17, 5, 21};
    }
    else if (cluster_type == 2)
    {
        // 12sites 3*2
        old_sites = {0, 16, 4, 20, 8, 24, 1, 17, 5, 21, 9, 25};
    }
    else if (cluster_type == 3)
    {
        // 18sites 3*3
        old_sites = {0, 16, 4, 20, 8, 24, 1, 17, 5, 21, 9, 25, 2, 18, 6, 22, 10, 26};
    }
    else if (cluster_type == 4)
    {
        // 24 sites C6
        old_sites = {28, 12, 24, 8, 20, 14, 26, 10, 22, 6, 18, 2, 17, 5, 21, 9, 25, 13, 29, 3, 19, 7, 23, 11};
    }
    else if (cluster_type == 5)
    {
        // 4*3
        old_sites = {0, 16, 4, 20, 8, 24, 12, 28, 1, 17, 5, 21, 9, 25, 13, 29, 2, 18, 6, 22, 10, 26, 14, 30};
    }

    for (int j = 0; j < old_sites.size(); j++)
    {
        int old_value0 = old_sites[j];
        int new_value0 = new_sites[j];
        Map_sites(old_value0, new_value0, old_sites_location, new_sites_location);
    }

    std::ofstream my_file1("old_Rspace_location.dat");
    for (auto elem : old_sites_location)
    {
        my_file1 << std::setw(4) << std::fixed << std::right << std::get<0>(elem) << "    " << std::setw(20) << std::fixed << std::right << std::setprecision(16) << std::get<1>(elem) << "    " << std::get<2>(elem) << std::endl;
    }

    std::ofstream my_file2("new_Rspace_location.dat");

    for (auto elem : new_sites_location)
    {

        my_file2 << std::setw(4) << std::fixed << std::right << std::get<0>(elem) << "    " << std::setw(20) << std::fixed << std::right << std::setprecision(16) << std::get<1>(elem) << "    " << std::get<2>(elem) << std::endl;
    }

}

void Lattice::Map_sites(int old_value, int new_value, const std::vector<std::tuple<int, double, double>> &old_sites, std::vector<std::tuple<int, double, double>> &new_sites)
{
    for (auto &tupleA : old_sites)
    {
        int intValueA = std::get<0>(tupleA);
        double doubleValueA1 = std::get<1>(tupleA);
        double doubleValueA2 = std::get<2>(tupleA);

        if (intValueA == old_value)
        {
            for (auto &tupleB : new_sites)
            {
                int intValueB = std::get<0>(tupleB);
                if (intValueB == new_value)
                {
                    std::get<1>(tupleB) = doubleValueA1;
                    std::get<2>(tupleB) = doubleValueA2;
                }
            }
        }
    }
}

void Lattice::Cal_b()
{
    double *vects = new double[9];
    std::fill(vects, vects + 9, 0.0);

    double *a2a3 = new double[3];
    double *a3a1 = new double[3];
    double *a1a2 = new double[3];

    Vec_Cross(a2, a3, a2a3);
    Vec_Cross(a3, a1, a3a1);
    Vec_Cross(a1, a2, a1a2);

    double area = Vec_Dot(3, a1, 1, a2a3, 1);

    if (fabs(area) < 1e-10)
    {
        std::cerr << "ERROR! Area is zero" << std::endl;
        delete[] a2a3;
        delete[] a3a1;
        delete[] a1a2;
        return;
    }

    double temp = 2 * PI / area;

    Vec_ax(3, temp, a2a3, vects);
    Vec_ax(3, temp, a3a1, vects + 3);
    Vec_ax(3, temp, a1a2, vects + 6);

    delete[] a2a3;
    delete[] a3a1;
    delete[] a1a2;

    for (int i = 0; i < 3; i++)
    {
        std::cout << "b" << i + 1 << ": {";
        for (int j = 0; j < 3; j++)
        {
            std::cout << std::fixed << std::setprecision(8) << std::setw(11) << std::right
                      << vects[i * 3 + j];
            if (j < 2)
            {
                std::cout << "  ";
            }
        }
        std::cout << " }" << std::endl;
    }
    std::cout << std::endl;

    b1[0] = vects[0];
    b1[1] = vects[1];

    b2[0] = vects[3];
    b2[1] = vects[4];

    delete[] vects;
}

void Lattice::Cal_BZ()
{

    // 枚举 n1, n2 的所有组合
    std::vector<std::pair<std::vector<double>, double>> points;
    for (int n1 = -2; n1 <= 2; n1++)
    {
        for (int n2 = -2; n2 <= 2; n2++)
        {
            if (n1 == 0 && n2 == 0)
                continue; // 跳过原点

            // 计算格点 G = n1*b1 + n2*b2
            std::vector<double> G = {n1 * b1[0] + n2 * b2[0], n1 * b1[1] + n2 * b2[1]};

            // 计算模长平方
            double dist2 = G[0] * G[0] + G[1] * G[1];
            points.emplace_back(G, dist2);
            // 输出调试信息
            // cout << "G: (" << G[0] << ", " << G[1] << ") -> dist2: " << dist2 << "\n";
        }
    }

    // 按距离平方排序
    sort(points.begin(), points.end(), [](const std::pair<std::vector<double>, double> &a, const std::pair<std::vector<double>, double> &b)
         { return a.second < b.second; });

    // 输出排序后的格点和距离
    // cout << "Sorted points by distance squared:\n";
    // for (const auto &[vec, dist2] : points)
    // cout << "G: (" << vec[0] << ", " << vec[1] << ") -> dist2: " << dist2 << "\n";

    // 查找最小和第二小距离平方值
    double min_dist2 = std::numeric_limits<double>::max();
    double second_dist2 = std::numeric_limits<double>::max();

    for (const auto &point : points)
    {
        const auto &vec = point.first;    // 解构 pair 的第一个元素
        const auto &dist2 = point.second; // 解构 pair 的第二个元素

        if (dist2 < min_dist2)
        {
            second_dist2 = min_dist2;
            min_dist2 = dist2;
        }
        else if (dist2 < second_dist2 && dist2 != min_dist2)
        {
            second_dist2 = dist2;
        }
    }

    // 提取所有最小距离平方值的格点
    std::vector<std::vector<double>> min_distance_points;
    for (const auto &point : points)
    {
        const auto &vec = point.first;    // 手动解构 pair 的第一个元素
        const auto &dist2 = point.second; // 手动解构 pair 的第二个元素

        if (abs(dist2 - min_dist2) < 1e-8)
        {
            //  cout << "min_dist2 = " << min_dist2 << endl;
            min_distance_points.push_back(vec);
        }
    }

    // 提取所有第二小距离平方值的格点
    std::vector<std::vector<double>> second_distance_points;
    for (const auto &point : points)
    {
        const auto &vec = point.first;    // 手动解构 pair 的第一个元素
        const auto &dist2 = point.second; // 手动解构 pair 的第二个元素

        if (std::abs(dist2 - second_dist2) < 1e-8)
        {
            // cout << "second_dist2 = " << second_dist2 << endl;
            second_distance_points.push_back(vec);
        }
    }

    // 将最小距离平方值的格点存入第一个矩阵中
    std::vector<std::vector<double>> nearest_neighbors;
    for (const auto &point : min_distance_points)
    {
        nearest_neighbors.push_back(point);
    }

    // 将第二小距离平方值的格点存入第二个矩阵中
    std::vector<std::vector<double>> second_neighbors;
    for (const auto &point : second_distance_points)
    {
        second_neighbors.push_back(point);
    }

    //*********************first*************************

    // F_计算每个坐标的角度，并存储在一个向量中
    std::vector<std::pair<double, std::vector<double>>> F_anglesAndCoords;
    for (const auto &coord : nearest_neighbors)
    {
        double angle = atan2(coord[1], coord[0]);
        if (angle < 0)
        {
            angle += 2 * M_PI; // 调整角度到 [0, 2π] 范围
        }
        F_anglesAndCoords.push_back({angle, coord});
    }
    // F_根据角度对坐标进行排序
    sort(F_anglesAndCoords.begin(), F_anglesAndCoords.end(),
         [](const std::pair<double, std::vector<double>> &a, const std::pair<double, std::vector<double>> &b)
         {
             return a.first < b.first;
         });
    // F_提取排序后的坐标
    std::vector<std::vector<double>> F_sorted_neighbors;
    for (const auto &pair : F_anglesAndCoords)
    {
        F_sorted_neighbors.push_back(pair.second);
    }

    std::vector<std::vector<double>> tmp_F_BZ(7, std::vector<double>(2, 0.0));
    tmp_F_BZ[0] = {0, 0}; // 起点
    for (int i = 1; i <= 6; i++)
    {
        tmp_F_BZ[i] = F_sorted_neighbors[i - 1];
    }

    // 计算垂直平分线的中点
    std::vector<std::vector<double>> F_midpoints(6, std::vector<double>(2));
    for (int i = 1; i <= 6; i++)
    {
        F_midpoints[i - 1][0] = tmp_F_BZ[i][0] / 2.0;
        F_midpoints[i - 1][1] = tmp_F_BZ[i][1] / 2.0;
        cout << std::fixed << std::setprecision(8)
             << "FBZ points" << i << ": (" << std::setw(11) << std::right << F_midpoints[i - 1][0] << std::setw(13) << F_midpoints[i - 1][1] << ")" << endl;
    }
    cout << endl;

    // 计算垂直平分线的斜率 m = y/x;则 垂直平分线的斜率为：-x/y
    std::vector<double> F_perp_slopes(6);
    for (int i = 1; i <= 6; i++)
    {
        if ((tmp_F_BZ[i][1]) != 0) // y是否为0
        {
            F_perp_slopes[i - 1] = -tmp_F_BZ[i][0] / tmp_F_BZ[i][1];
        }
        else
        {
            F_perp_slopes[i - 1] = 1e100; // 平分线垂直于 x 轴，斜率无穷大
        }
    }

    ofstream F_BZ("FBZ_location.dat");

    std::vector<std::pair<double, double>> F_intersection_points;
    // 计算相邻垂直平分线的交点
    for (int i = 0; i < 6; i++)
    {
        int j = (i + 1) % 6; // 只考虑相邻的垂直平分线
        if (F_perp_slopes[i] > 1e80)
        {                                                                                          // 如果第 i 条是竖直线
            double x_inter = F_midpoints[i][0];                                                    // 竖直线的 x 坐标
            double y_inter = F_perp_slopes[j] * (x_inter - F_midpoints[j][0]) + F_midpoints[j][1]; // 计算对应的 y 坐标
            F_intersection_points.push_back({x_inter, y_inter});

            F_BZ << fixed << setw(14) << setprecision(9) << x_inter << "    " << setw(14) << y_inter;
            F_BZ << endl;
        }
        else if (F_perp_slopes[j] > 1e80)
        {                                                                                          // 如果第 j 条是竖直线
            double x_inter = F_midpoints[j][0];                                                    // 竖直线的 x 坐标
            double y_inter = F_perp_slopes[i] * (x_inter - F_midpoints[i][0]) + F_midpoints[i][1]; // 计算对应的 y 坐标
            F_intersection_points.push_back({x_inter, y_inter});

            F_BZ << fixed << setw(14) << setprecision(9) << x_inter << "    " << setw(14) << y_inter;
            F_BZ << endl;
        }
        else if (abs(F_perp_slopes[i] - F_perp_slopes[j]) > 1e-6)
        { // 一般情况
            double x_inter = (F_perp_slopes[i] * F_midpoints[i][0] - F_perp_slopes[j] * F_midpoints[j][0] + F_midpoints[j][1] - F_midpoints[i][1]) / (F_perp_slopes[i] - F_perp_slopes[j]);
            double y_inter = F_perp_slopes[i] * (x_inter - F_midpoints[i][0]) + F_midpoints[i][1];
            F_intersection_points.push_back({x_inter, y_inter});

            F_BZ << fixed << setw(14) << setprecision(9) << x_inter << "    " << setw(14) << y_inter;
            F_BZ << endl;
        }
    }
    // 将第一个交点重新写入文件，形成循环
    if (!F_intersection_points.empty())
    {
        F_BZ << fixed << setw(14) << setprecision(9) << F_intersection_points[0].first << "    " << setw(14) << F_intersection_points[0].second;
        F_BZ << endl;
    }

    // ********************Second******************************

    std::vector<std::pair<double, std::vector<double>>> S_anglesAndCoords;
    for (const auto &coord : second_neighbors)
    {
        double angle = atan2(coord[1], coord[0]);
        if (angle < 0)
        {
            angle += 2 * M_PI; // 调整角度到 [0, 2π] 范围
        }
        S_anglesAndCoords.push_back({angle, coord});
    }
    sort(S_anglesAndCoords.begin(), S_anglesAndCoords.end(),
         [](const std::pair<double, std::vector<double>> &a, const std::pair<double, std::vector<double>> &b)
         {
             return a.first < b.first;
         });

    std::vector<std::vector<double>> S_sorted_neighbors;
    for (const auto &pair : S_anglesAndCoords)
    {
        S_sorted_neighbors.push_back(pair.second);
    }

    std::vector<std::vector<double>> tmp_S_BZ(7, std::vector<double>(2, 0.0));
    tmp_S_BZ[0] = {0, 0}; // 起点
    for (int i = 1; i <= 6; i++)
    {
        tmp_S_BZ[i] = S_sorted_neighbors[i - 1];
    }

    std::vector<std::vector<double>> S_midpoints(6, std::vector<double>(2));
    for (int i = 1; i <= 6; i++)
    {
        S_midpoints[i - 1][0] = tmp_S_BZ[i][0] / 2.0;
        S_midpoints[i - 1][1] = tmp_S_BZ[i][1] / 2.0;
        cout << std::fixed << std::setprecision(8)
             << "SBZ points" << i << ": (" << std::setw(11) << std::right << S_midpoints[i - 1][0] << std::setw(13) << std::right << S_midpoints[i - 1][1] << ")" << endl;
    }
    cout << endl;

    std::vector<double> S_perp_slopes(6);
    for (int i = 1; i <= 6; i++)
    {
        if ((tmp_S_BZ[i][1]) != 0)
        {
            S_perp_slopes[i - 1] = -tmp_S_BZ[i][0] / tmp_S_BZ[i][1];
        }
        else
        {
            S_perp_slopes[i - 1] = 1e100; // 平分线垂直于 x 轴
        }
    }

    ofstream S_BZ("SBZ_location.dat");

    std::vector<std::pair<double, double>> S_intersection_points;
    // 计算相邻垂直平分线的交点  m = y/x;则 垂直平分线的斜率为：-x/y
    for (int i = 0; i < 6; i++)
    {
        int j = (i + 1) % 6; // 只考虑相邻的垂直平分线
        if (S_perp_slopes[i] > 1e80)
        {                                                                                          // 如果第 i 条是竖直线
            double x_inter = S_midpoints[i][0];                                                    // 竖直线的 x 坐标
            double y_inter = S_perp_slopes[j] * (x_inter - S_midpoints[j][0]) + S_midpoints[j][1]; // 计算对应的 y 坐标
            S_intersection_points.push_back({x_inter, y_inter});
            S_BZ << fixed << setw(14) << setprecision(9) << x_inter << "    " << setw(14) << y_inter;
            S_BZ << endl;
        }
        else if (S_perp_slopes[j] > 1e80)
        {                                                                                          // 如果第 j 条是竖直线
            double x_inter = S_midpoints[j][0];                                                    // 竖直线的 x 坐标
            double y_inter = S_perp_slopes[i] * (x_inter - S_midpoints[i][0]) + S_midpoints[i][1]; // 计算对应的 y 坐标
            S_intersection_points.push_back({x_inter, y_inter});
            S_BZ << fixed << setw(14) << setprecision(9) << x_inter << "    " << setw(14) << y_inter;
            S_BZ << endl;
        }
        else if (abs(S_perp_slopes[i] - S_perp_slopes[j]) > 1e-6)
        { // 一般情况
            double x_inter = (S_perp_slopes[i] * S_midpoints[i][0] - S_perp_slopes[j] * S_midpoints[j][0] + S_midpoints[j][1] - S_midpoints[i][1]) / (S_perp_slopes[i] - S_perp_slopes[j]);
            double y_inter = S_perp_slopes[i] * (x_inter - S_midpoints[i][0]) + S_midpoints[i][1];
            S_intersection_points.push_back({x_inter, y_inter});
            S_BZ << fixed << setw(14) << setprecision(9) << x_inter << "    " << setw(14) << y_inter;
            S_BZ << endl;
        }
    }
    // 将第一个交点重新写入文件，形成循环
    if (!S_intersection_points.empty())
    {
        S_BZ << fixed << setw(14) << setprecision(9) << S_intersection_points[0].first << "    " << setw(14) << S_intersection_points[0].second;
        S_BZ << endl;
    }
}

void Lattice::BuildMomLine()
{
    pair<double, double> b1p = {b1[0], b1[1]};
    pair<double, double> b2p = {0, b2[1]};

    pair<double, double> _b1p = {-b1[0], -b1[1]};
    pair<double, double> _b2p = {-0, -b2[1]};
    const vector<pair<double, double>> mom_line;

    pair<double, double> gamma = {0.0, 0.0};
    pair<double, double> gammap1 = {b2[0] + b1[0], b2[1] + b1[1]};
    pair<double, double> gammap2 = {-b2[0], -b2[1]};
    pair<double, double> xp = {(-b2[0] + b1[0]) / 2.0, (-b2[1] + b1[1]) / 2.0};
    pair<double, double> x = {b1[0], 0};

    vector<pair<double, double>> Discrete_line1;
    Discrete_DSSF_line(gamma, gammap1, 2, 0, Discrete_line1);
    Output_Mpoint_to_file(Discrete_line1, "DSSF_Mline_location.dat");

    vector<pair<double, double>> Discrete_line2;
    Discrete_DSSF_line(gammap1, gammap2, 6, 0, Discrete_line2);
    Output_Mpoint_to_file(Discrete_line2, "DSSF_Mline_location.dat");

    vector<pair<double, double>> Discrete_line3;
    Discrete_DSSF_line(gammap2, xp, 1, 0, Discrete_line3);
    Output_Mpoint_to_file(Discrete_line3, "DSSF_Mline_location.dat");

    vector<pair<double, double>> Discrete_line4;
    Discrete_DSSF_line(xp, gamma, 3, 0, Discrete_line4);
    Output_Mpoint_to_file(Discrete_line4, "DSSF_Mline_location.dat");

    vector<pair<double, double>> Discrete_line5;
    Discrete_DSSF_line(gamma, x, 3, 0, Discrete_line5);
    Output_Mpoint_to_file(Discrete_line5, "DSSF_Mline_location.dat");

    vector<pair<double, double>> Discrete_line6;
    Discrete_DSSF_line(x, xp, 3, 0, Discrete_line6);
    Discrete_line6.push_back(xp);
    Output_Mpoint_to_file(Discrete_line6, "DSSF_Mline_location.dat");
}

void Lattice::Discrete_DSSF_line(std::pair<double, double> a, std::pair<double, double> b, int n, int per, std::vector<std::pair<double, double>> &Discrete_line)
{
    // 计算每段的x和y的增量
    double deltaX = (b.first - a.first) / n;
    double deltaY = (b.second - a.second) / n;

    if (per != 0 && per != 1)
    {
        cout << "per must be 0 or 1" << endl;
        return;
    }

    int num = n + per;
    for (int i = 0; i < num; i++)
    {
        // 不包含 b 点
        std::pair<double, double> segmentPoint;
        segmentPoint.first = a.first + i * deltaX;
        segmentPoint.second = a.second + i * deltaY;
        Discrete_line.push_back(segmentPoint);
    }
}

void Lattice::BuildMomSpace()
{
    pair<double, double> gamma = {0.0, 0.0};
    pair<double, double> gammap1 = {b2[0] + b1[0], b2[1] + b1[1]};
    pair<double, double> gammap2 = {-b2[0], -b2[1]};
    pair<double, double> xp = {(-b2[0] + b1[0]) / 2.0, (-b2[1] + b1[1]) / 2.0};
    pair<double, double> x = {b1[0], 0};

    pair<double, double> b1p = {b1[0], b1[1]};
    pair<double, double> b2p = {0, b2[1]};

    pair<double, double> _b1p = {-b1[0], -b1[1]};
    pair<double, double> _b2p = {-0, -b2[1]};

    ofstream MspaceF("DSSF_Mspace_location.dat", ios::app);
    MspaceF << fixed << setw(12) << setprecision(8) << right << b2p.first << setw(16) << right << b2p.second << endl;

    pair<double, double> site1 = {(_b1p.first + b2p.first) / 2, (_b1p.second + b2p.second) / 2};
    pair<double, double> site2 = {(gammap1.first + b2p.first) / 2, (gammap1.second + b2p.second) / 2};
    vector<pair<double, double>> Discrete_space1;
    Discrete_DSSF_line(site1, site2, 3, 1, Discrete_space1);
    Output_Mpoint_to_file(Discrete_space1, "DSSF_Mspace_location.dat");

    pair<double, double> site3 = _b1p;
    pair<double, double> site4 = gammap1;
    vector<pair<double, double>> Discrete_space2;
    Discrete_DSSF_line(site3, site4, 6, 1, Discrete_space2);
    Output_Mpoint_to_file(Discrete_space2, "DSSF_Mspace_location.dat");

    pair<double, double> site5 = {-5.23598776, 1.81379936};
    pair<double, double> site6 = {5.23598776, 1.81379936};
    vector<pair<double, double>> Discrete_space3;
    Discrete_DSSF_line(site5, site6, 5, 1, Discrete_space3);
    Output_Mpoint_to_file(Discrete_space3, "DSSF_Mspace_location.dat");

    pair<double, double> site7 = {_b1p.first, 0};
    pair<double, double> site8 = {b1p.first, 0};
    vector<pair<double, double>> Discrete_space4;
    Discrete_DSSF_line(site7, site8, 6, 1, Discrete_space4);
    Output_Mpoint_to_file(Discrete_space4, "DSSF_Mspace_location.dat");

    pair<double, double> site9 = {site5.first, -site5.second};
    pair<double, double> site10 = {site6.first, -site6.second};
    vector<pair<double, double>> Discrete_space5;
    Discrete_DSSF_line(site9, site10, 5, 1, Discrete_space5);
    Output_Mpoint_to_file(Discrete_space5, "DSSF_Mspace_location.dat");

    pair<double, double> site11 = {site3.first, -site3.second};
    pair<double, double> site12 = {site4.first, -site4.second};
    vector<pair<double, double>> Discrete_space6;
    Discrete_DSSF_line(site11, site12, 6, 1, Discrete_space6);
    Output_Mpoint_to_file(Discrete_space6, "DSSF_Mspace_location.dat");

    pair<double, double> site13 = {site1.first, -site1.second};
    pair<double, double> site14 = {site2.first, -site2.second};
    vector<pair<double, double>> Discrete_space7;
    Discrete_DSSF_line(site13, site14, 3, 1, Discrete_space7);
    Output_Mpoint_to_file(Discrete_space7, "DSSF_Mspace_location.dat");

    MspaceF << fixed << setw(12) << setprecision(8) << right << -b2p.first << setw(16) << right << -b2p.second << endl;
}

void Lattice::Output_Mpoint_to_file(vector<pair<double, double>> &Discrete_Mpoints, const string &Filename)
{
    ofstream OutFile(Filename, ios::app);
    if (!OutFile)
    {
        cout << "Failed to open:" << Filename << endl;
        return;
    }

    for (const auto &point : Discrete_Mpoints)
    {
        OutFile << fixed << setw(12) << setprecision(8) << right << point.first << setw(16) << right << point.second << endl;
    }
}