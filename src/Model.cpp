#include <Model.h>
#include <common.h>

using namespace std;

spinModel::spinModel(Lattice *lattice_, double theta_, double Hmag_) : lattice(lattice_), theta(theta_), Hmag(Hmag_)
{

    Nsites = lattice->get_Nsites();

    Nbonds = lattice->get_Nbonds();

    Nplaqs = lattice->get_Nplaqs();

    Plaquettes = lattice->get_Plaquettes();

    Bonds = lattice->get_Bonds();

    //****************************** */

    SpinVal = 0.5;
    Sdim = 2 * SpinVal + 1;

    interAllMatrix = new my_cplx[Sdim * Sdim * Sdim * Sdim];
    transVec = new my_cplx[Sdim * Sdim];

    cout << "-------Spin Parameters--------" << endl;
    cout << "Spin: " << setprecision(1) << right << SpinVal << endl;
    cout << "Sdim: " << right << Sdim << endl;
    cout << endl;

    J = 0.0;
    A = 0.0;
    G2 = 0.0;

    K = cos(theta * PI);
    G1 = -sin(theta * PI);

    Htheta = 0.0; // polar angle of field from c* direction
    Hpsi = 0.0;   // azimuthal angle of field from a -> c direction

    cout << "-------Ham Parameters--------" << endl;
    cout << fixed << setprecision(6)
         << " J:" << setw(10) << right << J << endl
         << " K:" << setw(10) << right << K << endl
         << "G1:" << setw(10) << right << G1 << endl
         << "G2:" << setw(10) << right << G2 << endl
         << " A:" << setw(10) << right << A << endl
         << " h:" << setw(10) << right << Hmag << endl;

    cout << endl;
}

// 每个键方向的各种相互作用的类型 SxSx; SxSy; SxSz; SySx...
string intTerms[9] = {"XX", "XY", "XZ",
                      "YX", "YY", "YZ",
                      "ZX", "ZY", "ZZ"};

string intTrans[3] = {"X", "Y", "Z"};

void spinModel::build_calcmod()
{
    ofstream calcmodF;

    calcmodF.open("calcmod.def"); // ios::app表示叠加

    if (!calcmodF)
    {
        cout << "failed to open calcmod.def" << endl;
    }

    calcmodF << "#CalcType = 0:Lanczos 1:TPQCalc 2:FullDiag 3:CG 4:Time-evolution 5:cTPQ" << endl;
    calcmodF << "#CalcModel = 0:Hubbard 1:Spin 2:Kondo 3:HubbardGC 4:SpinGC 5:KondoGC" << endl;
    calcmodF << "#Restart = 0:None 1:Restart_out 2:Restart_in 3:Restart_in + Restart_out" << endl;
    calcmodF << "#calcSpec = 0:None 1:Normal 2:No H*Phi 3:Save 4: Restart 5:Restart&Save" << endl;
    calcmodF << "# if CalcType = 1 or 5, InitialVecType shoule be -1." << endl;
    calcmodF << "CalcType         3" << endl;
    calcmodF << "CalcModel        4" << endl;
    calcmodF << "ReStart          0" << endl;
    calcmodF << "CalcSpec         0" << endl;
    calcmodF << "CalcEigenVec     0" << endl;
    calcmodF << "InitialVecType   0" << endl;
    calcmodF << "InputEigenVec    0" << endl;
    calcmodF << "OutputEigenVec   1" << endl;
    calcmodF << "InputHam         0" << endl;
    calcmodF << "OutputHam        0" << endl;
    calcmodF << "OutputExVec      0" << endl;

    calcmodF.close();
}

void spinModel::build_modpara()
{

    ofstream modparaF;

    modparaF.open("modpara.def");

    if (!modparaF)
    {
        cout << "failed to open modpara.def" << endl;
    }

    modparaF << "--------------------" << endl;
    modparaF << "Model_Parameters   0" << endl;
    modparaF << "--------------------" << endl;
    modparaF << "HPhi_Cal_Parameters" << endl;
    modparaF << "--------------------" << endl;
    modparaF << "CDataFileHead    zvo" << endl;
    modparaF << "CParaFileHead    zqp" << endl;
    modparaF << "--------------------" << endl;
    modparaF << "Nsite            " << Nsites << endl;
    modparaF << "Lanczos_max      10000" << endl;
    modparaF << "initial_iv       -1" << endl;
    modparaF << "exct             5" << endl;
    modparaF << "LanczosEps       16" << endl;
    modparaF << "LanczosTarget    5" << endl;
    modparaF << "LargeValue       1" << endl;
    modparaF << "NumAve           3" << endl;
    modparaF << "ExpecInterval    5" << endl;
    modparaF << "ExpandCoef       12" << endl;
    modparaF << "NOmega           200" << endl;
    modparaF << "OmegaMax         3.900000000000000e+01     0.000000000000000e+00" << endl;
    modparaF << "OmegaMin         -3.900000000000000e+01    0.000000000000000e+00" << endl;
    modparaF << "OmegaOrg         0.000000000000000e+00     0.000000000000000e+00" << endl;
    // modparaF << "PreCG            1  " << endl;

    modparaF.close();
}

void spinModel::build_namelist()
{
    ofstream namelistF;

    namelistF.open("namelist.def"); // ios::app表示叠加

    if (!namelistF)
    {
        cout << "failed to open namelist.def" << endl;
    }

    namelistF << "CalcMod         calcmod.def" << endl;
    namelistF << "ModPara         modpara.def" << endl;
    namelistF << "LocSpin         locspn.def" << endl;
    namelistF << "Trans           trans.def" << endl;
    namelistF << "InterAll        interAll.def" << endl;
    namelistF << "OneBodyG        greenone.def" << endl;
    namelistF << "TwoBodyG        greentwo.def" << endl;
    namelistF << "ThreeBodyG      greenthree.def" << endl;
    namelistF << "fourBodyG       greenfour.def" << endl;
    namelistF << "SixBodyG        greensix.def" << endl;
    namelistF << "SpectrumVec     zvo_eigenvec_0" << endl;

    namelistF.close();
}

//==================================================================
//   Type of  interaction  matrix between site i and j
//
//                Jxx  Jxy  Jxz
//                Jyx  Jyy  Jyz
//                Jzx  Jzy  Jzz
//
//==================================================================

void spinModel::interMatrix()
{

    /*     interaction matrix between two sites       */

    Xbond[0] = J + K, Xbond[1] = G2, Xbond[2] = G2;
    Xbond[3] = G2, Xbond[4] = J, Xbond[5] = G1;
    Xbond[6] = G2, Xbond[7] = G1, Xbond[8] = J;

    Ybond[0] = J, Ybond[1] = G2, Ybond[2] = G1;
    Ybond[3] = G2, Ybond[4] = J + K, Ybond[5] = G2;
    Ybond[6] = G1, Ybond[7] = G2, Ybond[8] = J;

    Zbond[0] = J, Zbond[1] = G1, Zbond[2] = G2;
    Zbond[3] = G1, Zbond[4] = J, Zbond[5] = G2;
    Zbond[6] = G2, Zbond[7] = G2, Zbond[8] = J + K;

    SIA[0] = A, SIA[1] = A, SIA[2] = A;
    SIA[3] = A, SIA[4] = A, SIA[5] = A;
    SIA[6] = A, SIA[7] = A, SIA[8] = A;

    /* Heisenberg Model*/
    bond[0] = J, bond[1] = 0.0, bond[2] = 0.0;
    bond[3] = 0.0, bond[4] = J, bond[5] = 0.0;
    bond[6] = 0.0, bond[7] = 0.0, bond[8] = J;

    cout << "-------Model Parameters--------" << endl;
    cout << endl;
    cout << "X-bond" << endl;
    for (auto i = 0; i < 3; i++)
    {
        for (auto j = 0; j < 3; j++)
        {
            cout << setw(10) << fixed << setprecision(6) << Xbond[i * 3 + j] << '\t';
        }

        cout << endl;
    }

    cout << endl;
    cout << "Y-bond" << endl;
    for (auto i = 0; i < 3; i++)
    {
        for (auto j = 0; j < 3; j++)
        {
            cout << setw(10) << fixed << setprecision(6) << Ybond[i * 3 + j] << '\t';
        }

        cout << endl;
    }

    cout << endl;
    cout << "Z-bond" << endl;
    for (auto i = 0; i < 3; i++)
    {
        for (auto j = 0; j < 3; j++)
        {
            cout << setw(10) << fixed << setprecision(6) << Zbond[i * 3 + j] << '\t';
        }

        cout << endl;
    }
    cout << endl;

    // Xbond[i]*= 0.5似乎是因为Interall要求是厄米共轭，所以每两个点之间的相互作用算了两次
    for (auto i = 0; i < 9; i++)
    {
        Xbond[i] *= 0.5;
        Ybond[i] *= 0.5;
        Zbond[i] *= 0.5;
    }
}

void spinModel::transMatrix()
{

    vector<double> ha = {1.0, 1.0, -2.0};
    Normalize(ha.size(), ha.data()); // crystallographic a direction of honeycomb

    vector<double> hb = {-1.0, 1.0, 0.0};
    Normalize(hb.size(), hb.data()); // crystallographic b direction of honeycomb

    vector<double> hc = {1.0, 1.0, 1.0};
    Normalize(hc.size(), hc.data()); // crystallographic c* direction of honeycomb

    // 假定 Htheta 和 Hpsi 是度数，转换为弧度
    trans[0] = Hmag * (sin(PI * Htheta / 180.0) * (cos(PI * Hpsi / 180.0) * ha[0] + sin(PI * Hpsi / 180.0) * hb[0]) + cos(PI * Htheta / 180.0) * hc[0]);
    trans[1] = Hmag * (sin(PI * Htheta / 180.0) * (cos(PI * Hpsi / 180.0) * ha[1] + sin(PI * Hpsi / 180.0) * hb[1]) + cos(PI * Htheta / 180.0) * hc[1]);
    trans[2] = Hmag * (sin(PI * Htheta / 180.0) * (cos(PI * Hpsi / 180.0) * ha[2] + sin(PI * Hpsi / 180.0) * hb[2]) + cos(PI * Htheta / 180.0) * hc[2]);

    cout << "-------Magnetic Parameters--------" << endl;
    cout << "      H = " << setw(9) << Hmag << endl;
    cout << "  theta = " << setw(9) << Htheta / 180.0 << " PI" << endl;
    cout << "    psi = " << setw(9) << Hpsi / 180.0 << " PI" << endl;
    cout << "   hvec =  (" << trans[0] << ", " << trans[1] << ", " << trans[2] << ") xyz" << endl;
    cout << "   hvec =  (" << sin(PI * Htheta / 180.0) * cos(PI * Hpsi / 180.0) << ", " << sin(PI * Htheta / 180.0) * sin(PI * Hpsi / 180.0) << ", " << cos(PI * Htheta / 180.0) << ") abc*" << endl;
    cout << endl;
}

void spinModel::build_interAll()
{
    // 先给各个格点的之间的相互作用赋值
    interMatrix();

    for (int b = 0; b < Nbonds; b++)
    {
        const auto &bond = Bonds[b];
        if (std::get<2>(bond) == 0) // Type 0 X
        {
            setBondInteraction("X");
            int xbond_x = std::get<0>(bond);
            int xbond_y = std::get<1>(bond);

            twositeInteraction(xbond_x, xbond_y);
            twositeInteraction(xbond_y, xbond_x);

            clearInter();
        }
        else if (std::get<2>(bond) == 1) // Type 1 Y
        {
            setBondInteraction("Y");
            int ybond_x = std::get<0>(bond);
            int ybond_y = std::get<1>(bond);

            twositeInteraction(ybond_x, ybond_y);
            twositeInteraction(ybond_y, ybond_x);

            clearInter();
        }
        else if (std::get<2>(bond) == 2) // Type 2 Z
        {
            setBondInteraction("Z");
            int zbond_x = std::get<0>(bond);
            int zbond_y = std::get<1>(bond);

            twositeInteraction(zbond_x, zbond_y);
            twositeInteraction(zbond_y, zbond_x);

            clearInter();
        }
        else
        {
            std::cerr << "ERROR!!! Lattice" << std::endl;
        }
    }

    // 单离子各向异性
    singleIonAnisotropy();

    /* add the number of two spin interaction term to interAll.def file  */

    fstream of("interAll.def", ios::in | ios::out | ios::binary);

    if (!of.is_open())
    {
        cout << "Failed to open the file." << endl;
    }

    // 移动文件指针到第58个字节位置
    of.seekp(58, ios::beg);

    string strNum = to_string(interactionInterAll);

    // 写入字符串到文件
    of.write(strNum.c_str(), strNum.size());

    of.close();
}

void spinModel::build_interAll(int oneD)
{
    // 一维系统
    //  先给各个格点的之间的相互作用赋值
    interMatrix();

    for (int i = 0; i < Nsites; i++)
    {
        int j = (i + 1) % Nsites;

        // 读取每个格点的相互作用类型以及相互作用强度，然后计算各个格点的相互作用(twoOpterm)
        setBondInteraction();

        // 生成一系列 interAll.def 文件所需要的数字
        twositeInteraction(i, j); // 这里不计算 twositeInteraction(i, j); 是正确的，为什么？

        // 避免发生相互作用结果发生累加，所以没作用一次，就清除一次。
        clearInter();
    }

    /* add the number of two spin interaction term to interAll.def file  */

    fstream of("interAll.def", ios::in | ios::out | ios::binary);

    if (!of.is_open())
    {
        cout << "Failed to open the file." << endl;
    }

    // 移动文件指针到第58个字节位置
    of.seekp(58, ios::beg);

    string strNum = to_string(interactionInterAll);

    // 写入字符串到文件
    of.write(strNum.c_str(), strNum.size());

    of.close();
}

void spinModel::build_trans()
{
    // 设置磁场参数：磁场大小，方向，角度...

    transMatrix();

    for (int i = 0; i < Nsites; i++)
    {
        // 将自旋算符转化为升降算符的系数
        setField();
        // 遍历 interactionTrans，然后往 trans.def 文件里面填入值
        transInteraction(i);

        clearTrans();
    }

    /* add the number of two spin interaction term to interAll.def file  */

    fstream of("trans.def", ios::in | ios::out | ios::binary);

    if (!of.is_open())
    {
        cout << "Failed to open the trans.def." << endl;
    }

    // 移动文件指针到第58个字节位置
    of.seekp(64, ios::beg);

    string strNum = to_string(interactionTrans);

    // 写入字符串到文件
    of.write(strNum.c_str(), strNum.size());

    of.close();
}
/*    localSpin file  */
void spinModel::build_localSpin()
{
    ofstream localSpin("locspn.def");

    if (!localSpin)
    {
        cerr << "Failed to open locspn.def" << endl;
    }

    localSpin << "=================================================" << endl;
    localSpin << "NlocalSpin" << "         " << Nsites << endl;
    localSpin << "=================================================" << endl;
    localSpin << "================" << "i_0IteElc_halfLocSpn" << "================" << endl;
    localSpin << "=================================================" << endl;

    for (auto i = 0; i < Nsites; i++)
    {
        localSpin << left << setw(2) << i << setw(8) << right << int(SpinVal * 2) << endl;
    }

    localSpin.close();
}

void spinModel::build_oneBodyG()
{
    ofstream oneBodyG;

    oneBodyG.open("greenone.def");

    if (!oneBodyG)
    {
        cout << "failed to open oneBody.def" << endl;
    }

    oneBodyG << "================================================" << endl;
    oneBodyG << "NCisAjs" << "             " << Nsites * Sdim * Sdim << endl;
    oneBodyG << "================================================" << endl;
    oneBodyG << "========" << "One Body Green Function" << "==========" << endl;
    oneBodyG << "================================================" << endl;

    for (auto i = 0; i < Nsites; i++)
    {
        for (auto s1 = 0; s1 < Sdim; s1++)
        {
            for (auto s2 = 0; s2 < Sdim; s2++)
            {
                oneBodyG << '\t' << i << '\t' << s1
                         << '\t' << i << '\t' << s2 << endl;
            }
        }
    }

    oneBodyG.close();
}

void spinModel::build_twoBodyG()
{
    ofstream twoBodyG("greentwo.def");

    if (!twoBodyG)
    {
        cout << "failed to open twoBody.def" << endl;
    }

    twoBodyG << "===============================================" << endl;
    twoBodyG << "NCisAjsCktAltDC" << "               " << pow(Nsites, 2) * pow(Sdim, 4) << endl;
    twoBodyG << "===============================================" << endl;
    twoBodyG << "========" << "Two Body Green Function" << "==========" << endl;
    twoBodyG << "===============================================" << endl;

    for (auto i = 0; i < Nsites; i++)
    {
        for (auto s1 = 0; s1 < Sdim * Sdim; s1++)
        {
            auto sigma1 = s1 / Sdim;
            auto sigma2 = s1 % Sdim;

            for (auto j = 0; j < Nsites; j++)
            {
                for (auto s2 = 0; s2 < Sdim * Sdim; s2++)
                {
                    auto sigma3 = s2 / Sdim;
                    auto sigma4 = s2 % Sdim;

                    twoBodyG << '\t' << i << '\t' << sigma1 << '\t' << i << '\t' << sigma2
                             << '\t' << j << '\t' << sigma3 << '\t' << j << '\t' << sigma4 << endl;
                }
            }
        }
    }
    twoBodyG.close();
}

void spinModel::build_threeBodyG()
{
    ofstream threeBodyG("greenthree.def");

    if (!threeBodyG)
    {
        cout << "failed to open greenthree.def" << endl;
        return;
    }

    threeBodyG << "===============================================" << endl;
    threeBodyG << "NCisAjsCktAltDCEmuFnu" << "               "
               << fixed << setprecision(0) << Nsites * (Nsites - 1) * (Nsites - 2) / 6 * pow(Sdim, 6) << endl;
    threeBodyG << "===============================================" << endl;
    threeBodyG << "==========" << "Three Body Green Function" << "============" << endl;
    threeBodyG << "===============================================" << endl;

    // Nc3 = N(N-1)(N-2)/3! terms, 2024 for 24 sites
    // 确保 i < j < k

    for (auto i = 0; i < Nsites; i++)
    {
        for (auto j = i + 1; j < Nsites; j++)
        {
            for (auto k = j + 1; k < Nsites; k++)
            {
                for (auto s1 = 0; s1 < Sdim * Sdim; s1++)
                {
                    auto sigma1 = s1 / Sdim;
                    auto sigma2 = s1 % Sdim;
                    for (auto s2 = 0; s2 < Sdim * Sdim; s2++)
                    {
                        auto sigma3 = s2 / Sdim;
                        auto sigma4 = s2 % Sdim;

                        for (auto s3 = 0; s3 < Sdim * Sdim; s3++)
                        {

                            auto sigma5 = s3 / Sdim;
                            auto sigma6 = s3 % Sdim;

                            threeBodyG << '\t' << i << '\t' << sigma1 << '\t' << i << '\t' << sigma2
                                       << '\t' << j << '\t' << sigma3 << '\t' << j << '\t' << sigma4
                                       << '\t' << k << '\t' << sigma5 << '\t' << k << '\t' << sigma6 << endl;
                        }
                    }
                }
            }
        }
    }
    threeBodyG.close();
}

void spinModel::build_fourBodyG()
{
    // 读取键信息（bondType bond1 bond2）
    my_mat BondIndex;
    Mat_Read("check_bond.dat", BondIndex);
    int numBonds = BondIndex.size();

    // 打开输出文件
    ofstream fourBodyG("greenfour.def");
    if (!fourBodyG)
    {
        cerr << "Error: Failed to open greenfour.def" << endl;
        return;
    }

    // 文件头（包含精确条目数）
    fourBodyG << "===============================================" << endl;
    fourBodyG << "NCisAjsCktAltDCEmuFnuGxvHyw" << "       "
              << fixed << setprecision(0) << numBonds * numBonds * pow(Sdim, 8) << endl;
    fourBodyG << "===============================================" << endl;
    fourBodyG << "=========== Four Body Green Function ==========" << endl;
    fourBodyG << "===============================================" << endl;

    // 遍历所有键对
    for (auto b1 = 0; b1 < numBonds; b1++)
    {
        uint32_t i = BondIndex[b1][1]; // 键b1的起点格点
        uint32_t j = BondIndex[b1][2]; // 键b1的终点格点

        for (auto b2 = 0; b2 < numBonds; b2++)
        {
            uint32_t k = BondIndex[b2][1]; // 键b2的起点格点
            uint32_t l = BondIndex[b2][2]; // 键b2的终点格点

            for (auto s1 = 0; s1 < Sdim * Sdim; s1++)
            {
                auto sigma1 = s1 / Sdim;
                auto sigma2 = s1 % Sdim;

                for (auto s2 = 0; s2 < Sdim * Sdim; s2++)
                {
                    auto sigma3 = s2 / Sdim;
                    auto sigma4 = s2 % Sdim;

                    for (auto s3 = 0; s3 < Sdim * Sdim; s3++)
                    {
                        auto sigma5 = s3 / Sdim;
                        auto sigma6 = s3 % Sdim;
                        for (auto s4 = 0; s4 < Sdim * Sdim; s4++)
                        {
                            auto sigma7 = s4 / Sdim;
                            auto sigma8 = s4 % Sdim;

                            // 写入完整索引：i σ1 σ2, j σ3 σ4, k σ5 σ6, l σ7 σ8
                            fourBodyG << '\t' << i << '\t' << sigma1 << '\t' << i << '\t' << sigma2
                                      << '\t' << j << '\t' << sigma3 << '\t' << j << '\t' << sigma4
                                      << '\t' << k << '\t' << sigma5 << '\t' << k << '\t' << sigma6
                                      << '\t' << l << '\t' << sigma7 << '\t' << l << '\t' << sigma8 << endl;
                        }
                    }
                }
            }
        }
    }

    fourBodyG.close();
}

void spinModel::build_sixBodyG()
{

    ofstream plaqFile("check_plaq.dat");
    if (!plaqFile)
    {
        cout << "failed to open check_plaq.dat" << endl;
        return;
    }

    ofstream sixBodyG("greensix.def");
    if (!sixBodyG)
    {
        cerr << "Error: Failed to open greensix.def" << endl;
        return;
    }

    sixBodyG << "================================================" << endl;
    sixBodyG << "NCisAjsCktAltDCEmuFnuGxvHywIzvJawKbxLcy" << "    "
             << fixed << setprecision(0) << Nsites / 2 * pow(Sdim, 12) << endl; // 6格点 × 2自旋/格点 → 12自旋分量
    sixBodyG << "=================================================" << endl;
    sixBodyG << "=========== Six Body Green Function ============" << endl;
    sixBodyG << "=================================================" << endl;

    // 遍历所有plaquettes
    for (int p = 0; p < Plaquettes.size(); p++)
    {

        // 提取6个格点索引
        int i = std::get<0>(Plaquettes[p]);
        int j = std::get<1>(Plaquettes[p]);
        int k = std::get<2>(Plaquettes[p]);
        int l = std::get<3>(Plaquettes[p]);
        int m = std::get<4>(Plaquettes[p]);
        int n = std::get<5>(Plaquettes[p]);

        plaqFile << i << "  " << j << "  " << k << "  " << l << "  " << m << "  " << n << endl;

        for (auto s1 = 0; s1 < Sdim * Sdim; s1++)
        {
            auto sigma1 = s1 / Sdim;
            auto sigma2 = s1 % Sdim;

            for (auto s2 = 0; s2 < Sdim * Sdim; s2++)
            {
                auto sigma3 = s2 / Sdim;
                auto sigma4 = s2 % Sdim;

                for (auto s3 = 0; s3 < Sdim * Sdim; s3++)
                {
                    auto sigma5 = s3 / Sdim;
                    auto sigma6 = s3 % Sdim;

                    for (auto s4 = 0; s4 < Sdim * Sdim; s4++)
                    {
                        auto sigma7 = s4 / Sdim;
                        auto sigma8 = s4 % Sdim;

                        for (auto s5 = 0; s5 < Sdim * Sdim; s5++)
                        {
                            auto sigma9 = s5 / Sdim;
                            auto sigma10 = s5 % Sdim;

                            for (auto s6 = 0; s6 < Sdim * Sdim; s6++)
                            {
                                auto sigma11 = s6 / Sdim;
                                auto sigma12 = s6 % Sdim;

                                // 写入索引：6个格点 × 2自旋分量
                                sixBodyG << '\t' << i << '\t' << sigma1 << '\t' << i << '\t' << sigma2
                                         << '\t' << j << '\t' << sigma3 << '\t' << j << '\t' << sigma4
                                         << '\t' << k << '\t' << sigma5 << '\t' << k << '\t' << sigma6
                                         << '\t' << l << '\t' << sigma7 << '\t' << l << '\t' << sigma8
                                         << '\t' << m << '\t' << sigma9 << '\t' << m << '\t' << sigma10
                                         << '\t' << n << '\t' << sigma11 << '\t' << n << '\t' << sigma12 << endl;
                            }
                        }
                    }
                }
            }
        }
    }
    plaqFile.close();
    sixBodyG.close();
}

void spinModel::setBondInteraction()
{
    // chain model

    // 两个嵌套的循环遍历 3x3 的矩阵
    // 读取每个格点的相互作用强度，然后计算各个格点的相互作用
    for (auto i = 0; i < 3; i++)
    {
        for (auto j = 0; j < 3; j++)
        {
            twoOpterm(intTerms[i * 3 + j], bond[i * 3 + j]);
        }
    }

    // cout << "Bond Interactions" << endl;
}

void spinModel::setField()
{
    // 磁场
    for (auto i = 0; i < 3; i++)
    {
        oneOpterm(intTrans[i], trans[i]);
    }
}

void spinModel::singleIonAnisotropy()
{
    // 单离子各向异性
    clearInter();

    for (auto j = 0; j < 3; j++)
    {
        for (auto k = 0; k < 3; k++)
        {
            twoOpterm(intTerms[j * 3 + k], SIA[j * 3 + k]);
        }
    }

    for (auto i = 0; i < Nsites; i++)
    {

        onesiteInteraction(i);
    }
}

void spinModel::setBondInteraction(string bondType)
{
    // 二维系统，需要先确定键的作用类型

    assert(bondType == "X" || bondType == "Y" || bondType == "Z");

    for (auto i = 0; i < 3; i++)
    {
        // 先选择键的方向(因为不同键方向的相同耦合作用的耦合强度可能是不同的)，再选择相互作用类型，再选择键方向相互作用类型的耦合强度
        for (auto j = 0; j < 3; j++)
        {
            if (bondType == "X")
            {
                twoOpterm(intTerms[i * 3 + j], Xbond[i * 3 + j]);
            }
            else if (bondType == "Y")
            {
                twoOpterm(intTerms[i * 3 + j], Ybond[i * 3 + j]);
            }
            else
            {
                twoOpterm(intTerms[i * 3 + j], Zbond[i * 3 + j]);
            }
        }
    }
}

void spinModel::twoOpterm(string termType, double coeff)
{

    /********************************
            Sp = Sx + i Sy
            Sm = Sx - i Sy
            Sx = 1/2 (Sp + Sm)
            Sy = 1/2i (Sp - Sm)
    ********************************/
    // 算符的系数 * 耦合强度
    if (termType == "XX")
    {
        PPterm(0.25 * coeff);
        PMterm(0.25 * coeff);
        MPterm(0.25 * coeff);
        MMterm(0.25 * coeff);
    }
    else if (termType == "XY")
    {
        PPterm(-0.25 * cplxi * coeff);
        PMterm(0.25 * cplxi * coeff);
        MPterm(-0.25 * cplxi * coeff);
        MMterm(0.25 * cplxi * coeff);
    }
    else if (termType == "XZ")
    {
        PZterm(0.5 * coeff);
        MZterm(0.5 * coeff);
    }
    else if (termType == "YX")
    {
        PPterm(-0.25 * cplxi * coeff);
        PMterm(-0.25 * cplxi * coeff);
        MPterm(0.25 * cplxi * coeff);
        MMterm(0.25 * cplxi * coeff);
    }
    else if (termType == "YY")
    {
        PPterm(-0.25 * coeff);
        PMterm(0.25 * coeff);
        MPterm(0.25 * coeff);
        MMterm(-0.25 * coeff);
    }
    else if (termType == "YZ")
    {
        PZterm(-0.5 * cplxi * coeff);
        MZterm(0.5 * cplxi * coeff);
    }
    else if (termType == "ZX")
    {
        ZPterm(0.5 * coeff);
        ZMterm(0.5 * coeff);
    }
    else if (termType == "ZY")
    {
        ZPterm(-0.5 * cplxi * coeff);
        ZMterm(0.5 * cplxi * coeff);
    }
    else if (termType == "ZZ")
    {
        ZZterm(coeff);
    }
    else if (termType == "PP")
    {
        PPterm(coeff);
    }
    else if (termType == "PM")
    {
        PMterm(coeff);
    }
    else if (termType == "PZ")
    {
        PZterm(coeff);
    }
    else if (termType == "MP")
    {
        MPterm(coeff);
    }
    else if (termType == "MM")
    {
        MMterm(coeff);
    }
    else if (termType == "MZ")
    {
        MZterm(coeff);
    }
    else if (termType == "ZP")
    {
        ZPterm(coeff);
    }
    else if (termType == "ZM")
    {
        ZMterm(coeff);
    }
    else
    {
        cout << "Interaction term " << termType << " is not defined!";
    }
}

void spinModel::oneOpterm(string termType, double coeff)
{
    /********************************
           Sp = Sx + i Sy
           Sm = Sx - i Sy
           Sx = 1/2 (Sp + Sm)
           Sy = 1/2i (Sp - Sm)
   ********************************/

    // 算符的系数 * 耦合强度
    if (termType == "X")
    {
        Pterm(0.5 * coeff);
        Mterm(0.5 * coeff);
    }
    else if (termType == "Y")
    {
        Pterm(-0.5 * cplxi * coeff);
        Mterm(0.5 * cplxi * coeff);
    }
    else if (termType == "Z")
    {
        Zterm(coeff);
    }
}

void spinModel::twositeInteraction(const int siteIndexa, const int siteIndexb)
{
    for (auto i = 0; i < Sdim * Sdim; i++)
    {
        auto sigma1 = i / Sdim;
        auto sigma2 = i % Sdim;
        for (auto j = 0; j < Sdim * Sdim; j++)
        {
            auto sigma3 = j / Sdim;
            auto sigma4 = j % Sdim;

            auto interReal = interAllMatrix[i * Sdim * Sdim + j].real();
            auto interImag = interAllMatrix[i * Sdim * Sdim + j].imag();

            // cout <<"Real" << '\t' << interReal <<'\t'<< "Imag" << '\t' << interImag << endl;

            if (abs(interReal) > ACCURATECUT || abs(interImag) > ACCURATECUT)
            {
                interactionInterAll += 1;

                addToInterAll("interAll.def", siteIndexa, sigma1, sigma2, siteIndexb, sigma3, sigma4, interReal, interImag);
            }
        }
    }
}

void spinModel::onesiteInteraction(const int siteIndex)
{
    for (auto i = 0; i < Sdim * Sdim; i++)
    {
        auto sigma1 = i / Sdim;
        auto sigma2 = i % Sdim;
        for (auto j = 0; j < Sdim * Sdim; j++)
        {
            auto sigma3 = j / Sdim;
            auto sigma4 = j % Sdim;

            auto interReal = interAllMatrix[i * Sdim * Sdim + j].real();
            auto interImag = interAllMatrix[i * Sdim * Sdim + j].imag();

            if (abs(interReal) > ACCURATECUT || abs(interImag) > ACCURATECUT)
            {
                interactionInterAll += 1;

                addToInterAll("interAll.def", siteIndex, sigma1, sigma2, siteIndex, sigma3, sigma4, interReal, interImag);
            }
        }
    }
}

void spinModel::transInteraction(const int siteIndex)
{
    for (auto i = 0; i < Sdim * Sdim; i++)
    {
        auto sigma1 = i / Sdim; // 0 0 1 1
        auto sigma2 = i % Sdim; // 0 1 0 1

        auto interReal = transVec[i].real();
        auto interImag = transVec[i].imag();

        // if (abs(interReal) > ACCURATECUT || abs(interImag) > ACCURATECUT)
        //{
        interactionTrans += 1;

        addToTrans("trans.def", siteIndex, sigma1, siteIndex, sigma2, interReal, interImag);
        //}
    }
}

double spinModel::CGcoeff(double sigma1, string opa, double sigma2, string opb)
{
    assert(opa == "P" || opa == "M" || opa == "Z");
    assert(opb == "P" || opb == "M" || opb == "Z");

    double coeffa, coeffb;

    if (opa == "P" || opa == "M")
    {
        coeffa = sqrt(SpinVal * SpinVal + SpinVal - sigma1 * sigma1 - sigma1);
    }
    else
    {
        coeffa = -sigma1;
    }

    if (opb == "P" || opb == "M")
    {
        coeffb = sqrt(SpinVal * SpinVal + SpinVal - sigma2 * sigma2 - sigma2);
    }
    else
    {
        coeffb = -sigma2;
    }

    return coeffa * coeffb;
}

double spinModel::CGcoeff(double sigma, string op)
{
    assert(op == "P" || op == "M" || op == "Z");

    double coeff;

    if (op == "P" || op == "M")
    {
        coeff = sqrt(SpinVal * SpinVal + SpinVal - sigma * sigma - sigma);
    }
    else
    {
        coeff = -sigma;
    }

    return coeff;
}

void spinModel::Pterm(const my_cplx coeff)
{
    double sigma;

    for (auto i = 0; i < Sdim - 1; i++)
    {
        sigma = i - (Sdim - 1) * 0.5;
        // 如果是自旋0.5，那么sigma1 = -0.5

        transVec[i * Sdim + i + 1] += CGcoeff(sigma, "P") * coeff; // index = 1; CGcoeff= 1.0
    }
}

void spinModel::Mterm(const my_cplx coeff)
{
    double sigma;

    for (auto i = 0; i < Sdim - 1; i++)
    {
        sigma = i - (Sdim - 1) * 0.5;
        // 如果是自旋0.5，那么sigma1 = -0.5

        transVec[((i + 1) * Sdim + i)] += CGcoeff(sigma, "M") * coeff; // index = 2; CGcoeff= 1.0
    }
}

void spinModel::Zterm(const my_cplx coeff)
{
    double sigma;

    for (auto i = 0; i < Sdim; i++)
    {
        sigma = i - (Sdim - 1) * 0.5;
        // 如果是自旋0.5，那么sigma1 = -0.5

        transVec[i * Sdim + i] += CGcoeff(sigma, "Z") * coeff;

        // index = 0; CGcoeff= 0.5
        // index = 3; CGcoeff= -0.5
    }
}

void spinModel::PPterm(const my_cplx coeff) // here coeff= j * 0.25
{
    double sigma1, sigma2;

    for (auto i = 0; i < Sdim - 1; i++)
    {
        sigma1 = i - (Sdim - 1) * 0.5;
        // 如果是自旋0.5，那么sigma1 = -0.5
        // 如果是自旋1，那么sigma1 = -1, 0
        for (auto j = 0; j < Sdim - 1; j++)
        {
            sigma2 = j - (Sdim - 1) * 0.5;
            // 如果是自旋0.5，那么sigma2 = -0.5
            // 如果是自旋1，那么sigma2 = -1, 0

            interAllMatrix[(i * Sdim + i + 1) * Sdim * Sdim + j * Sdim + j + 1] += CGcoeff(sigma1, "P", sigma2, "P") * coeff;
            //  5
        }
    }
}

void spinModel::PMterm(const my_cplx coeff)
{
    double sigma1, sigma2;

    for (auto i = 0; i < Sdim - 1; i++)
    {
        sigma1 = i - (Sdim - 1) * 0.5;
        for (auto j = 0; j < Sdim - 1; j++)
        {

            sigma2 = j - (Sdim - 1) * 0.5;
            // CGcoeff = 1
            interAllMatrix[(i * Sdim + i + 1) * Sdim * Sdim + (j + 1) * Sdim + j] += CGcoeff(sigma1, "P", sigma2, "M") * coeff;
            // 6
        }
    }
}

void spinModel::PZterm(const my_cplx coeff)
{
    double sigma1, sigma2;
    for (auto i = 0; i < Sdim - 1; i++)
    {
        sigma1 = i - (Sdim - 1) * 0.5;
        for (auto j = 0; j < Sdim; j++) // 0，1
        {
            sigma2 = j - (Sdim - 1) * 0.5;

            interAllMatrix[(i * Sdim + i + 1) * Sdim * Sdim + j * Sdim + j] += CGcoeff(sigma1, "P", sigma2, "Z") * coeff;
        }
    }
}

void spinModel::MPterm(const my_cplx coeff)
{
    double sigma1, sigma2;

    for (auto i = 0; i < Sdim - 1; i++)
    {
        sigma1 = i - (Sdim - 1) * 0.5;

        for (auto j = 0; j < Sdim - 1; j++)
        {
            sigma2 = j - (Sdim - 1) * 0.5;

            interAllMatrix[((i + 1) * Sdim + i) * Sdim * Sdim + j * Sdim + j + 1] += CGcoeff(sigma1, "M", sigma2, "P") * coeff;
        }
    }
}

void spinModel::MMterm(const my_cplx coeff)
{
    double sigma1, sigma2;

    for (auto i = 0; i < Sdim - 1; i++)
    {
        sigma1 = i - (Sdim - 1) * 0.5;

        for (auto j = 0; j < Sdim - 1; j++)
        {
            sigma2 = j - (Sdim - 1) * 0.5;
            // 5
            // CG:1
            interAllMatrix[((i + 1) * Sdim + i) * Sdim * Sdim + (j + 1) * Sdim + j] += CGcoeff(sigma1, "M", sigma2, "M") * coeff;
            // 10
        }
    }
}

void spinModel::MZterm(const my_cplx coeff)
{
    double sigma1, sigma2;

    for (auto i = 0; i < Sdim - 1; i++)
    {
        sigma1 = i - (Sdim - 1) * 0.5;

        for (auto j = 0; j < Sdim; j++)
        {
            sigma2 = j - (Sdim - 1) * 0.5;

            interAllMatrix[((i + 1) * Sdim + i) * Sdim * Sdim + j * Sdim + j] += CGcoeff(sigma1, "M", sigma2, "Z") * coeff;
            // 8,11
        }
    }
}
void spinModel::ZPterm(const my_cplx coeff)
{
    double sigma1, sigma2;

    for (auto i = 0; i < Sdim; i++)
    {
        sigma1 = i - (Sdim - 1) * 0.5;
        for (auto j = 0; j < Sdim - 1; j++)
        {
            sigma2 = j - (Sdim - 1) * 0.5;
            // 2
            interAllMatrix[(i * Sdim + i) * Sdim * Sdim + j * Sdim + j + 1] += CGcoeff(sigma1, "Z", sigma2, "P") * coeff;
            // 1,13
        }
    }
}

void spinModel::ZMterm(const my_cplx coeff)
{
    double sigma1, sigma2;

    for (auto i = 0; i < Sdim; i++)
    {
        sigma1 = i - (Sdim - 1) * 0.5;
        for (auto j = 0; j < Sdim - 1; j++)
        {
            sigma2 = j - (Sdim - 1) * 0.5;

            interAllMatrix[(i * Sdim + i) * Sdim * Sdim + (j + 1) * Sdim + j] += CGcoeff(sigma1, "Z", sigma2, "M") * coeff;
            // 2,14
            /* cout << "interAllMatrix: " << (i * Sdim + i) * Sdim * Sdim + (j + 1) * Sdim + j << endl;
            cout << "CG: " << CGcoeff(sigma1, "Z", sigma2, "M") << endl; */
        }
    }
}

void spinModel::ZZterm(const my_cplx coeff)
{
    double sigma1, sigma2;

    for (auto i = 0; i < Sdim; i++)
    {
        sigma1 = i - (Sdim - 1) * 0.5;

        for (auto j = 0; j < Sdim; j++)
        {
            sigma2 = j - (Sdim - 1) * 0.5;

            interAllMatrix[(i * Sdim + i) * Sdim * Sdim + j * Sdim + j] += CGcoeff(sigma1, "Z", sigma2, "Z") * coeff;
            // 0,3,12,15
        }
    }
}

void spinModel::addToInterAll(const string fileName, int indexa, double sigma1, double sigma2,
                              int indexb, double sigma3, double sigma4, double realpart, double imagpart)
{
    bool fileExists = std::filesystem::exists(fileName);

    ofstream outfile(fileName, ios::app);

    if (!outfile.is_open())
    {
        cerr << "Error: Unable to open file " << fileName << endl;
        return;
    }

    if (!fileExists)
    {
        // 如果文件不存在，写入表头内容
        outfile << "===========================================" << endl
                << "NInterAll" << "                 " << endl
                << "===========================================" << endl
                << "================zInterAll================" << endl
                << "===========================================" << endl;
    }

    // 写入数据行
    outfile << indexa << '\t' << sigma1 << '\t' << indexa << '\t' << sigma2 << '\t'
            << indexb << '\t' << sigma3 << '\t' << indexb << '\t' << sigma4 << '\t';
    outfile << fixed << setprecision(16) << setw(22) << right << realpart << setw(22) << right << imagpart << endl;
}

void spinModel::addToTrans(const string fileName, int indexa, double sigma1, int indexb, double sigma2, double realpart, double imagpart)
{
    bool fileExists = std::filesystem::exists(fileName);

    ofstream outfile(fileName, ios::app);
    if (!outfile.is_open())
    {
        // 错误处理：无法打开文件
        cerr << "Error: Unable to open file " << fileName << endl;
        return;
    }

    if (!fileExists)
    {

        // 如果文件不存在，写入表头内容
        outfile << "===========================================" << endl
                << "NTransfer" << "                 " << endl
                << "===========================================" << endl
                << "================i_j_s_tijs=================" << endl
                << "===========================================" << endl;
    }

    // 写入数据行
    outfile << indexa << '\t' << sigma1 << '\t' << indexb << '\t' << sigma2 << '\t'
            << fixed << setprecision(16) << setw(22) << right << realpart << setw(22) << right << imagpart << endl;
}

/* void spinModel::addToTrans(string fileName, int const indexa, double const sigma1, int indexb, double sigma2, double realpart, double imagpart)
{
    ofstream outfile;

    struct stat buffer;

    if (stat(fileName.c_str(), &buffer) != 0)
    {
        outfile.open(fileName, ios::app);

        if (!outfile)
        {
            return;
        }

        outfile << "===========================================" << endl;
        outfile << "NTransfer" << '\t' << "                " << endl;
        outfile << "===========================================" << endl;
        outfile << "================" << "i_j_s_tijs" << "=================" << endl;
        outfile << "===========================================" << endl;
        outfile << indexa << '\t' << sigma1 << '\t' << indexb << '\t' << sigma2 << '\t';

        outfile << setprecision(16) << fixed
                << setw(22) << right << realpart << setw(22) << right << imagpart << endl;

        outfile.close();
    }
    else
    {
        outfile.open(fileName, ios::app);

        if (!outfile)
        {
            return;
        }

        outfile << indexa << '\t' << sigma1 << '\t' << indexb << '\t' << sigma2 << '\t';
        outfile << setprecision(16) << fixed
                << setw(22) << right << realpart << setw(22) << right << imagpart << endl;

        outfile.close();
    }
} */

void spinModel::clearInter()
{
    for (auto i = 0; i < Sdim * Sdim; i++)
    {
        for (auto j = 0; j < Sdim * Sdim; j++)
        {
            interAllMatrix[i * Sdim * Sdim + j] = complex<double>{0, 0};
        }
    }
}

void spinModel::clearTrans()
{
    // transVec 向量的大小是 Sdim+1
    for (auto i = 0; i < Sdim * Sdim; i++)
    {
        transVec[i] = complex<double>{0, 0};
    }
}

/* void spinModel::build_sixBodyG(vector<vector<int>> plaqs)
{
    int num_Wp = plaqs.size();
    int num = 64 * num_Wp;
    // W_p = sigma_i^x * sigma_j^y * sigma_k^z * sigma_l^x * sigma_m^x * sigma_n^z
    ofstream sixBodyG("greensix.def");
    sixBodyG << "=============================================" << endl
             << "SixbodyG" << "       " << num << endl
             << "=============================================" << endl
             << "====== i s1 j s2 k s3 m s4 n s5 l s6 ========" << endl
             << "=============================================" << endl;

    for (int n = 0; n < num_Wp; n++)
    {
        int site0 = plaqs[n][0];
        int site1 = plaqs[n][1];
        int site2 = plaqs[n][2];
        int site3 = plaqs[n][3];
        int site4 = plaqs[n][4];
        int site5 = plaqs[n][5];

        for (int ii = 0; ii < 4; ii++)
        {
            for (int jj = 0; jj < 4; jj++)
            {

                sixBodyG << site0 << " " << 0 << "   " << site0 << " " << 1 << "   "
                         << site1 << " " << 0 << "   " << site1 << " " << 1 << "   ";

                if (jj == 0 || jj == 1)
                {
                    sixBodyG << site2 << " " << 0 << "   " << site2 << " " << 0 << "   ";
                    if (ii == 0)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        // if (inde1 % 2 == 0)
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        // else if (inde1 % 2 == 1)
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 1)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 2)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 3)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                }
                else if (jj == 2 || jj == 3)
                {
                    sixBodyG << site2 << " " << 1 << "   " << site2 << " " << 1 << "   ";
                    if (ii == 0)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 1)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 2)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 3)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                }
                sixBodyG << endl;
            }
        }

        for (int ii = 0; ii < 4; ii++)
        {
            for (int jj = 0; jj < 4; jj++)
            {

                sixBodyG << site0 << " " << 0 << "   " << site0 << " " << 1 << "   "
                         << site1 << " " << 1 << "   " << site1 << " " << 0 << "   ";

                if (jj == 0 || jj == 1)
                {
                    sixBodyG << site2 << " " << 0 << "   " << site2 << " " << 0 << "   ";
                    if (ii == 0)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 1)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 2)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 3)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                }
                else if (jj == 2 || jj == 3)
                {
                    sixBodyG << site2 << " " << 1 << "   " << site2 << " " << 1 << "   ";
                    if (ii == 0)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 1)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 2)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 3)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                }
                sixBodyG << endl;
            }
        }

        for (int ii = 0; ii < 4; ii++)
        {
            for (int jj = 0; jj < 4; jj++)
            {

                sixBodyG << site0 << " " << 1 << "   " << site0 << " " << 0 << "   "
                         << site1 << " " << 0 << "   " << site1 << " " << 1 << "   ";

                if (jj == 0 || jj == 1)
                {
                    sixBodyG << site2 << " " << 0 << "   " << site2 << " " << 0 << "   ";
                    if (ii == 0)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 1)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 2)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 3)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                }
                else if (jj == 2 || jj == 3)
                {
                    sixBodyG << site2 << " " << 1 << "   " << site2 << " " << 1 << "   ";
                    if (ii == 0)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 1)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 2)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 3)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                }
                sixBodyG << endl;
            }
        }

        for (int ii = 0; ii < 4; ii++)
        {
            for (int jj = 0; jj < 4; jj++)
            {

                sixBodyG << site0 << " " << 1 << "   " << site0 << " " << 0 << "   "
                         << site1 << " " << 1 << "   " << site1 << " " << 0 << "   ";

                if (jj == 0 || jj == 1)
                {
                    sixBodyG << site2 << " " << 0 << "   " << site2 << " " << 0 << "   ";
                    if (ii == 0)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 1)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 2)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 3)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 0)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 1)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                }
                else if (jj == 2 || jj == 3)
                {
                    sixBodyG << site2 << " " << 1 << "   " << site2 << " " << 1 << "   ";
                    if (ii == 0)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 1)
                    {
                        sixBodyG << site3 << " " << 0 << "   " << site3 << " " << 1 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 2)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 0 << "   " << site4 << " " << 1 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                    else if (ii == 3)
                    {
                        sixBodyG << site3 << " " << 1 << "   " << site3 << " " << 0 << "   "
                                 << site4 << " " << 1 << "   " << site4 << " " << 0 << "   ";
                        if (jj == 2)
                        {
                            sixBodyG << site5 << " " << 0 << "   " << site5 << " " << 0 << "   ";
                        }
                        else if (jj == 3)
                        {
                            sixBodyG << site5 << " " << 1 << "   " << site5 << " " << 1 << "   ";
                        }
                    }
                }
                sixBodyG << endl;
            }
        }
    }
} */

spinModel::~spinModel()
{
    delete[] interAllMatrix;
    delete[] transVec;
    // 注意：不要删除 lattice，因为它可能在其他地方使用
}