#include <vector>

// コンパイル時にlcmatrix.cppを一緒にコンパイルしなくてもいいようにする.
#define INCLUDE_ONLY_HEADER

#ifdef INCLUDE_ONLY_HEADER
        #ifndef LCMATRIX
        #define LCMATRIX
        #endif
#endif

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
                Matrix abs();
                double ave();
                Matrix dot(const Matrix&) const;
                Matrix t();
                Matrix operator * (const Matrix&) const;
                Matrix operator * (const double) const;
                Matrix operator + (const Matrix&) const;
                Matrix operator + (const double) const;
                Matrix operator - (const Matrix&) const;
                Matrix operator - (const double) const;
                Matrix operator - () const;
                Matrix operator / (const Matrix&) const;
                Matrix operator / (const double) const;
                Matrix operator <= (const Matrix&) const;
                Matrix operator <= (const double) const;
                Matrix pow(double) const;
                Matrix log() const;
                Matrix exp() const;
                void print(int n = 5) const;
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
        void print(Matrix, int n = 5);
} // namespace

#ifdef INCLUDE_ONLY_HEADER
#include "lcmatrix.cpp"
#endif
