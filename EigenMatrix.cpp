#include <iostream>
#include <ctime>
using namespace std;
#include <Eigen/Core>
#include <Eigen/Dense>
#define MATRIX_SIZE 50

int main(int argc, char** argv)
{
    Eigen::Matrix<float, 2, 3> matrix_23;
    Eigen::Vector3d v_3d;
    Eigen::Matrix3d matrix_33 = Eigen::Matrix3d::Zero();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_dynamic;
    Eigen::MatrixXd matrix_x;
    matrix_23 << 1,2,3,4,5,6;
    cout << "matrix_23:\n" << matrix_23 <<"\n" <<endl;
    for(int i=0; i<2; i++)
        for(int j=0; j<3; j++)
            cout << matrix_23(i,j) << endl;
    v_3d << 3,2,1;
    Eigen::Matrix<double, 2, 1> result = matrix_23.cast<double>() * v_3d;
    cout << "matrix multi vector = \n" << result << endl;
    matrix_33 = Eigen::Matrix3d::Random();
    cout << "matrix_33:\n" << matrix_33 << endl << endl;
    cout << "transpose of matrix_33:\n" << matrix_33.transpose() << endl;
    cout << "sum of matrix_33:\n" << matrix_33.sum() << endl;
    cout << "ji:\n" << matrix_33.trace() << endl;
    cout<<"Hello SLAM !"<<endl;
    return 0;
}
