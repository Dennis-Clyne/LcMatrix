#include <vector>

namespace LcMatrix {
        class Matrix {
                std::vector<std::vector<double>> matrix;
                int row, col;
                void setRC();

                public:
                Matrix();
                Matrix(std::vector<std::vector<double>>);
                Matrix(int, int, double x = 0.0);
                double get(int, int);
                void set(int, int, double);
                std::vector<std::vector<double>> getMatrix();
                int getRow();
                int getCol();
                double max();
                Matrix max(bool);
                Matrix argMax(bool = true);
                double sum();
                Matrix sum(bool);
                Matrix dot(Matrix);
                Matrix t();
                Matrix operator * (Matrix);
                Matrix operator * (double);
                Matrix operator + (Matrix);
                Matrix operator + (double);
                Matrix operator - (Matrix);
                Matrix operator - (double);
                Matrix operator - ();
                Matrix operator / (Matrix);
                Matrix operator / (double);
                Matrix pow(double);
                Matrix log();
                Matrix exp();
                void print();
        }; // class

        double max(Matrix);
        Matrix max(Matrix, bool);
        Matrix argMax(Matrix, bool = true);
        double sum(Matrix);
        Matrix sum(Matrix, bool);
        Matrix dot(Matrix, Matrix);
        Matrix t(Matrix);
        Matrix pow(Matrix, double);
        Matrix log(Matrix);
        Matrix exp(Matrix);
        void print(Matrix);
} // namespace
