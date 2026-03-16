#include <common.h>
#include <Model.h>
#include <Lattice.h>

int main(int argc, char **argv)
{

    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " <theta> <field>" << endl;
        return 1;
    }

    // 获取命令行参数
    double theta_ = atof(argv[1]); // 参数 1：theta
    double Field = atof(argv[2]);  // 参数 2：field

    cout << "theta: " << fixed << setprecision(4) << theta_ << endl;
    cout << "Field: " << fixed << setprecision(4) << Field << endl;

    //*************************************
    Lattice lattice(1); // 0:6;  1:2*2*2; 2:3*2*2; 3:3*3*2; 4:24_C6; 5:4*3*2
    //************************************
    spinModel Model(&lattice, theta_, Field);

    Model.build_localSpin();

    Model.build_calcmod();

    Model.build_modpara();

    Model.build_namelist();

    Model.build_interAll();

    Model.build_trans();

    Model.build_oneBodyG();

    Model.build_twoBodyG();

    Model.build_threeBodyG();

    Model.build_fourBodyG();

    Model.build_sixBodyG();

    return 0;
}
