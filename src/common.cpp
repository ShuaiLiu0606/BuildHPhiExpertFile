#include <common.h>

//****************** 向量归一 *********************
void Normalize(const int &dim, double *wf)
{
    if (wf == nullptr || dim <= 0)
    {
        return;
    }

    double norm = cblas_dnrm2(dim, wf, 1);
    if (norm > 0.0)
    {
        cblas_dscal(dim, 1.0 / norm, wf, 1);
    }
}

void Normalize(const int &dim, complex<double> *wf)
{
    if (wf == nullptr || dim <= 0)
    {
        return;
    }

    double norm = cblas_dznrm2(dim, wf, 1);
    if (norm > 0.0)
    {
        cblas_zdscal(dim, 1.0 / norm, wf, 1);
    }
}

double Vec_Dot(const int &dim, double *x, const int incx, double *y, const int incy)
{
    return cblas_ddot(dim, x, incx, y, incy);
}

my_cplx Vec_Dot(const int &dim, my_cplx *x, const int incx, my_cplx *y, const int incy)
{
    my_cplx aux;
    cblas_zdotc_sub(dim, x, incx, y, incy, &aux);
    return aux;
}

void Vec_Cross(const double *vec1, const double *vec2, double *result)
{
    fill(result, result + 3, 0.0);
    if (vec1 == nullptr || vec2 == nullptr || result == nullptr)
    {
        cerr << "Error: Null pointer passed to Vec_cross." << std::endl;
        return;
    }

    // 计算叉积
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

void Vec_Cross(const my_cplx *vec1, const my_cplx *vec2, my_cplx *result)
{
    fill(result, result + 3, complex<double>(0.0, 0.0));
    if (vec1 == nullptr || vec2 == nullptr || result == nullptr)
    {
        cerr << "Error: Null pointer passed to Vec_cross." << endl;
        return;
    }

    // 计算叉积
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

// y = a*x
void Vec_ax(const int &dim, const my_cplx &a, my_cplx *x, my_cplx *y)
{
    for (int i = 0; i < dim; i++)
    {
        y[i] = a * x[i];
    }
}

void Vec_ax(const int &dim, const double &a, double *x, double *y)
{
    for (int i = 0; i < dim; i++)
    {
        y[i] = a * x[i];
    }
}

// mat
void Mat_Read(const string &filename, vector<vector<double>> &matrix)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed open to file: " << filename << endl;
        return;
    }

    matrix.clear(); // 关键：清空旧数据，防止内存累积

    string line;
    while (getline(file, line))
    {
        vector<double> row;
        istringstream iss(line);
        double value;
        while (iss >> value)
        {
            row.push_back(value);
        }
        matrix.push_back(row);
    }

    file.close();
}

void Mat_Read(const string &filename, vector<vector<complex<double>>> &matrix)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed open to file: " << filename << endl;
        return;
    }

    matrix.clear(); // 关键：清空旧数据，防止内存累积

    string line;
    while (getline(file, line))
    {
        vector<complex<double>> row;
        istringstream iss(line);
        string value_str;

        // 从行中读取每个复数值
        while (iss >> value_str)
        {
            // 去除括号
            value_str.erase(remove(value_str.begin(), value_str.end(), '('), value_str.end());
            value_str.erase(remove(value_str.begin(), value_str.end(), ')'), value_str.end());

            // 找到逗号的位置
            size_t comma_pos = value_str.find(',');

            if (comma_pos != string::npos)
            {
                try
                {
                    // 提取实部和虚部字符串
                    string real_str = value_str.substr(0, comma_pos);
                    string imag_str = value_str.substr(comma_pos + 1);

                    // 转换为 double 并创建复数
                    double real = stod(real_str);
                    double imag = stod(imag_str);

                    complex<double> value(real, imag);
                    row.push_back(value);
                }
                catch (const std::invalid_argument &e)
                {
                    cerr << "Invalid number format in: " << value_str << endl;
                }
                catch (const std::out_of_range &e)
                {
                    cerr << "Number out of range in: " << value_str << endl;
                }
            }
        }

        matrix.push_back(row);
    }
    file.close();
}



