#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#ifndef LCMATRIX
#define LCMATRIX
#include "lcmatrix.hpp"
#endif

/*
 * Matrix class
 */
namespace LcMatrix {
        /*
         * matrixを正規分布で初期化する
         */
        void Matrix::initRand() {
                std::random_device rnd;
                std::mt19937_64 mt(rnd());
                // 平均0.0, 分散値1.0の正規分布
                std::normal_distribution<> norm(0.0, 1.0);

                matrix.clear();
                for (int i = 0; i < size; i++) {
                        matrix.push_back(norm(mt));
                }
        }

        /*
         * 宣言だけするとき
         */
        Matrix::Matrix() {
        }

        /*
         * 行列受け取ってmatrixを作る
         */
        Matrix::Matrix(std::vector<std::vector<double>> x) {
                row = x.size();
                col = x[0].size();
                size = row * col;
                matrix.reserve(size);

                for (auto i : x) {
                        for (double j : i) {
                                matrix.push_back(j);
                        }
                }
        }

        /*
         * i行j列のmatrixを要素分だけメモリ確保して作る
         */
        Matrix::Matrix(int i, int j) {
                row = i;
                col = j;
                size = row * col;
                matrix.reserve(size);
        }

        /*
         * i行j列のmatrixをxで初期化して作る
         */
        Matrix::Matrix(int i, int j, double x) {
                matrix = std::vector<double>(i * j, x);
                row = i;
                col = j;
                size = row * col;
        }

        /*
         * matrixのi行j列目の要素を返す
         * idnexは0から
         */
        double Matrix::get(int i, int j) {
                return matrix[col * i + j];
        }

        /*
         * matrixのi行j列目にxを代入する
         */
        void Matrix::set(int i, int j, double x) {
                matrix[col * i + j] = x;
        }

        /*
         * matrix自体を返す
         */
        std::vector<double> Matrix::getMatrix() {
                return matrix;
        }

        /*
         * matrixの行数を返す
         */
        int Matrix::getRow() {
                return row;
        }

        /*
         * matrixの列数を返す
         */
        int Matrix::getCol() {
                return col;
        }

        /*
         * max関数
         * matrixで最も大きな値を返す
         */
        double Matrix::max() {
                return *std::max_element(std::begin(matrix), std::end(matrix));
        }

        /*
         * max関数(bool)
         * 行または列内で最大の要素をMatrixオブジェクトとして返す.
         * true : 行
         * false : 列
         */
        Matrix Matrix::max(bool x) {
                auto ite = std::begin(matrix);
                auto end = std::end(matrix);

                // 行方向
                if (x) {
                        Matrix result(row, 1);
                        for (; ite != end; ite += col) {
                                result.matrix.push_back(*std::max_element(ite, ite + col));
                                //result.set(i, 0, *std::max_element(matrix[i].begin(), matrix[i].end()));
                        }

                        return result;
                } 
                // 列方向
                else {
                        Matrix result(1, col, 0);
                        auto r_ite = std::begin(result.matrix);
                        auto r_end = std::end(result.matrix);
                        auto r_ite_st = r_ite; // r_iteをとっとく

                        for (; ite != end; ite++, r_ite++) {
                                if (r_ite == r_end - 1)
                                        r_ite = r_ite_st;

                                if (*r_ite < *ite)
                                        *r_ite = *ite;
                        }

                        return result;
                }
        }

        /*
         * argMax関数
         * 行または列内で最大の要素のインデックスを行列として返す.
         * インデックスは0から始まる.
         * true : 行
         * false : 列
         * デフォルトは行.
         */
        Matrix Matrix::argMax(bool x) {
                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                // 行方向
                if (x) {
                        Matrix result(row, 1);
                        for (; ite != end; ite += col) {
                                result.matrix.push_back((int)std::distance(ite, std::max_element(ite, ite + col)));
                                //result.set(i, 0, (int)std::distance(matrix[i].begin(), std::max_element(matrix[i].begin(), matrix[i].end())));
                        }

                        return result;
                } 
                // 列方向
                else {
                        Matrix result(1, col, 0);
                        auto r_ite = std::begin(result.matrix);
                        auto r_end = std::end(result.matrix);
                        auto r_ite_st = r_ite; // r_iteをとっとく

                        Matrix max_mat(1, col, 0); // 最大値とっとく用
                        auto m_ite = std::begin(max_mat.matrix);
                        auto m_end = std::end(max_mat.matrix);
                        auto m_ite_st = m_ite; // m_iteをとっとく

                        int count = 0;
                        for (; ite != end; ite++, m_ite++, r_ite++) {
                                if (m_ite == m_end) {
                                        m_ite = m_ite_st;
                                        r_ite = r_ite_st;
                                        count++;
                                }

                                if (*m_ite < *ite) {
                                        *m_ite = *ite;
                                        *r_ite = count;
                                }
                        }

                        return result;
                }
        }

        /*
         * sum関数
         * 行列の全要素の和
         */
        double Matrix::sum() {
                double sum = std::accumulate(std::begin(matrix), std::end(matrix), 0.0);
                return sum;
        }

        /*
         * 行列の列または行方向の和
         * true : 行方向
         * false : 列方向
         */
        Matrix Matrix::sum(bool x) {
                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                // 行方向
                if (x) {
                        Matrix result(row, 1);
                        for (; ite != end; ite += col) {
                                result.matrix.push_back(std::accumulate(ite, ite + col, 0.0));
                        }

                        return result;
                } 
                // 列方向
                else {
                        Matrix result(1, col, 0);
                        auto r_ite = std::begin(result.matrix);
                        auto r_end = std::end(result.matrix);
                        auto r_ite_st = r_ite; // r_iteをとっとく

                        for (; ite != end; ite++, r_ite++) {
                                if (r_ite == r_end)
                                        r_ite = r_ite_st;

                                *r_ite += *ite;
                        }

                        return result;
                }
        }

        /*
         * matrixに含まれるxの数を返す関数
         */
        int Matrix::count(double x) {
                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                return std::count(ite, end, x);
        }

        /*
         * 行列の要素を絶対値にして返す
         */
        Matrix Matrix::abs() {
                Matrix result(row, col);
                for (double i : matrix) {
                        if (i < 0)
                                result.matrix.push_back(-i);
                        else
                                result.matrix.push_back(i);
                }

                return result;
        }

        /*
         * 行列の全要素の平均を返す
         */
        double Matrix::ave() {
                double sum = this->sum();
                return sum / size;
        }

        /*
         * 行列の内積
         */
        Matrix Matrix::dot(const Matrix &x) const {
                int xcol = x.col;
                Matrix result(row, xcol);

                double ans = 0;
                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < xcol; j++) {
                                for (int z = 0; z < col; z++) {
                                        ans += matrix[col * i + z] * x.matrix[xcol * z + j];
                                        //result.matrix[i][j] += matrix[i][z] * x.matrix[z][j];
                                }
                                result.matrix.push_back(ans);
                                ans = 0;
                        }
                }

                return result;
        }

        /*
         * 行列の転置
         */
        Matrix Matrix::t() {
                Matrix result(col, row);
                auto ite = std::begin(matrix);
                auto end = ite + col; // 1行目の終わりまで
                for (; ite !=end; ite++) {
                        for (int i = 0; i < row; i++) {
                                result.matrix.push_back(*(ite + col * i));
                                //result.matrix[j][i] = matrix[i][j];
                        }
                }
                return result;
        }

        /*
         * 行列の要素どうしの掛け算 (行列 * 行列)
         */
        Matrix Matrix::operator * (const Matrix &pre_x) const {
                Matrix result(row, col);

                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                auto pre_x_ite = std::begin(pre_x.matrix);

                // 行と列が同数のとき
                if (row == pre_x.row && col == pre_x.col) {
                        for (; ite != end; ite++, pre_x_ite++) {
                                result.matrix.push_back(*ite * *pre_x_ite);
                        }
                        return result;
                }

                auto pre_x_end = std::end(pre_x.matrix);

                // numpyのブロードキャストみたいなやつ
                // 行数同じで列が1つ
                auto pre_x_ite_str = pre_x_ite;
                if (row == pre_x.row && pre_x.col == 1) {
                        for (; ite != end; ite++, pre_x_ite++) {
                                if (pre_x_ite == pre_x_end - 1)
                                        pre_x_ite = pre_x_ite_str;
                                result.matrix.push_back(*ite * *pre_x_ite);
                        }
                }
                // 列数同じで行が1つ
                if (col == pre_x.col && pre_x.row == 1) {
                        int count = 0;
                        for (; ite != end; ite++, count++) {
                                if (count == row - 1) {
                                        pre_x_ite++;
                                        count = 0;
                                }
                                result.matrix.push_back(*ite * *pre_x_ite);
                        }
                }

                return result;
        }

        /*
         * 行列の各要素とdoubleの掛け算 (行列 * double)
         */
        Matrix Matrix::operator * (const double x) const {
                Matrix result(row, col);

                for (double i : matrix) {
                        result.matrix.push_back(i * x);
                }

                return result;
        }

        /*
         * 行列の各要素をx倍する.
         */
        void Matrix::operator *= (const double x) {
                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                for (; ite != end; ite++) {
                        *ite *= x;
                }
        }

        /*
         * 行列の和 (行列 + 行列)
         */
        Matrix Matrix::operator + (const Matrix &pre_x) const {
                Matrix result(row, col);

                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                auto pre_x_ite = std::begin(pre_x.matrix);

                // 行と列が同数のとき
                if (row == pre_x.row && col == pre_x.col) {
                        for (; ite != end; ite++, pre_x_ite++) {
                                result.matrix.push_back(*ite + *pre_x_ite);
                        }
                        return result;
                }

                auto pre_x_end = std::end(pre_x.matrix);

                // numpyのブロードキャストみたいなやつ
                // 行数同じで列が1つ
                auto pre_x_ite_str = pre_x_ite;
                if (row == pre_x.row && pre_x.col == 1) {
                        for (; ite != end; ite++, pre_x_ite++) {
                                if (pre_x_ite == pre_x_end - 1)
                                        pre_x_ite = pre_x_ite_str;
                                result.matrix.push_back(*ite + *pre_x_ite);
                        }
                }
                // 列数同じで行が1つ
                if (col == pre_x.col && pre_x.row == 1) {
                        int count = 0;
                        for (; ite != end; ite++, count++) {
                                if (count == row - 1) {
                                        pre_x_ite++;
                                        count = 0;
                                }
                                result.matrix.push_back(*ite + *pre_x_ite);
                        }
                }

                return result;
        }

        /*
         * 行列の和 (行列 + double)
         */
        Matrix Matrix::operator + (const double x) const {
                Matrix result(row, col);

                for (double i : matrix) {
                        result.matrix.push_back(i + x);
                }

                return result;
        }

        /*
         * 行列の差 (行列 - 行列)
         */
        Matrix Matrix::operator - (const Matrix &pre_x) const {
                Matrix result(row, col);

                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                auto pre_x_ite = std::begin(pre_x.matrix);

                // 行と列が同数のとき
                if (row == pre_x.row && col == pre_x.col) {
                        for (; ite != end; ite++, pre_x_ite++) {
                                result.matrix.push_back(*ite - *pre_x_ite);
                        }
                        return result;
                }

                auto pre_x_end = std::end(pre_x.matrix);

                // numpyのブロードキャストみたいなやつ
                // 行数同じで列が1つ
                auto pre_x_ite_str = pre_x_ite;
                if (row == pre_x.row && pre_x.col == 1) {
                        for (; ite != end; ite++, pre_x_ite++) {
                                if (pre_x_ite == pre_x_end - 1)
                                        pre_x_ite = pre_x_ite_str;
                                result.matrix.push_back(*ite - *pre_x_ite);
                        }
                }
                // 列数同じで行が1つ
                if (col == pre_x.col && pre_x.row == 1) {
                        int count = 0;
                        for (; ite != end; ite++, count++) {
                                if (count == row - 1) {
                                        pre_x_ite++;
                                        count = 0;
                                }
                                result.matrix.push_back(*ite - *pre_x_ite);
                        }
                }

                return result;
        }

        /*
         * 行列の差 (行列 - double)
         */
        Matrix Matrix::operator - (const double x) const {
                Matrix result(row, col);

                for (double i : matrix) {
                        result.matrix.push_back(i - x);
                }

                return result;
        }

        /*
         * 行列の各要素の符号反転
         */
        Matrix Matrix::operator - () const {
                Matrix result(row, col);

                for (double i : matrix) {
                        result.matrix.push_back(-i);
                }

                return result;
        }

        /*
         * 行列の各要素どうしの割り算 (行列 / 行列)
         */
        Matrix Matrix::operator / (const Matrix &pre_x) const {
                Matrix result(row, col);

                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                auto pre_x_ite = std::begin(pre_x.matrix);

                // 行と列が同数のとき
                if (row == pre_x.row && col == pre_x.col) {
                        for (; ite != end; ite++, pre_x_ite++) {
                                result.matrix.push_back(*ite / *pre_x_ite);
                        }
                        return result;
                }

                auto pre_x_end = std::end(pre_x.matrix);

                // numpyのブロードキャストみたいなやつ
                // 行数同じで列が1つ
                auto pre_x_ite_str = pre_x_ite;
                if (row == pre_x.row && pre_x.col == 1) {
                        for (; ite != end; ite++, pre_x_ite++) {
                                if (pre_x_ite == pre_x_end - 1)
                                        pre_x_ite = pre_x_ite_str;
                                result.matrix.push_back(*ite / *pre_x_ite);
                        }
                }
                // 列数同じで行が1つ
                if (col == pre_x.col && pre_x.row == 1) {
                        int count = 0;
                        for (; ite != end; ite++, count++) {
                                if (count == row - 1) {
                                        pre_x_ite++;
                                        count = 0;
                                }
                                result.matrix.push_back(*ite / *pre_x_ite);
                        }
                }

                return result;
        }

        /*
         * 行列の割り算 (行列 / double)
         */
        Matrix Matrix::operator / (const double x) const {
                Matrix result(row, col);

                for (double i : matrix) {
                        result.matrix.push_back(i / x);
                }

                return result;
        }

        /*
         * matrixとpre_xの各要素を比較して,
         * 等しければ1, そうでなければ0の行列を返す関数
         */
        Matrix Matrix::operator == (const Matrix &pre_x) const {
                Matrix result(row, col);

                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                auto pre_x_ite = std::begin(pre_x.matrix);

                // 行と列が同数のとき
                for (; ite != end; ite++, pre_x_ite++) {
                        if (*ite == *pre_x_ite)
                                result.matrix.push_back(1);
                        else 
                                result.matrix.push_back(0);
                }

                return result;
        }

        /*
         * matrixとxの各要素を比較して,
         * 等しければ1, そうでなければ0の行列を返す関数
         */
        Matrix Matrix::operator == (const double x) const {
                Matrix result(row, col);

                // 行と列が同数のとき
                for (double i : matrix) {
                        if (i == x)
                                result.matrix.push_back(1);
                        else 
                                result.matrix.push_back(0);
                }

                return result;
        }


        /*
         * matrixとpre_xの各要素を比較して,
         * x以下なら1, それ以外は0の行列を返す関数
         */
        Matrix Matrix::operator <= (const Matrix &pre_x) const {
                Matrix result(row, col);

                auto ite = std::begin(matrix);
                auto end = std::end(matrix);
                auto pre_x_ite = std::begin(pre_x.matrix);

                // 行と列が同数のとき
                for (; ite != end; ite++, pre_x_ite++) {
                        if (*ite <= *pre_x_ite)
                                result.matrix.push_back(1);
                        else 
                                result.matrix.push_back(0);
                }

                return result;
        }

        /*
         * matrixとxを比較して,
         * x以下なら1, それ以外は0の行列を返す関数
         */
        Matrix Matrix::operator <= (const double x) const {
                Matrix result(row, col);

                // 行と列が同数のとき
                for (double i : matrix) {
                        if (i <= x)
                                result.matrix.push_back(1);
                        else 
                                result.matrix.push_back(0);
                }

                return result;
        }

        /*
         * 行列の各要素のn乗
         */
        Matrix Matrix::pow(double n) const {
                Matrix result(row, col);

                for (double i : matrix) {
                        result.matrix.push_back(std::pow(i, n));
                }

                return result;
        }

        /*
         * 行列の各要素のlog() (自然対数)
         */
        Matrix Matrix::log() const {
                Matrix result(row, col);

                for (double i : matrix) {
                        result.matrix.push_back(std::log(i));
                }

                return result;
        }

        /*
         * 行列の各要素のexp()
         */
        Matrix Matrix::exp() const {
                Matrix result(row, col);

                for (double i : matrix) {
                        result.matrix.push_back(std::exp(i));
                }

                return result;
        }

        /*
         * 行列の表示
         * 小数点以下n桁まで表示(デフォルトはn = 5).
         */
        void Matrix::print(int n) const {
                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                //printf("%.15lf, ", j);
                                std::cout << std::fixed << std::setprecision(n) << matrix[col * i + j] << ", ";
                        }
                        std::cout << std::endl;
                }
                std::cout << std::endl;
        }
} // namespace Matrix class

namespace LcMatrix {
        /*
         * max関数
         */
        double max(Matrix x) {
                return x.max();
        }

        /*
         * max関数(bool)
         * 行または列内で最大の要素を行列として返す.
         * true : 行
         * false : 列
         * デフォルトは行.
         */
        Matrix max(Matrix x, bool i) {
                return x.max(i);
        }

        /*
         * argMax関数
         * 行または列内で最大の要素のインデックスを行列として返す.
         * インデックスは0から始まる.
         * true : 行
         * false : 列
         * デフォルトは行.
         */
        Matrix argMax(Matrix x, bool i) {
                return x.argMax(i);
        }

        /*
         * sum関数
         */
        double sum(Matrix x) {
                return x.sum();
        }

        /*
         * 行列の列または行方向の和
         * true : 行方向
         * false : 列方向
         */
        Matrix sum(Matrix x, bool i) {
                return x.sum(i);
        }

        /*
         * 行列の内積
         */
        Matrix dot(Matrix x, Matrix y) {
                return x.dot(y);
        }

        /*
         * 行列の転置
         */
        Matrix t(Matrix x) {
                return x.t();
        }

        /*
         * 行列の各要素のn乗
         */
        Matrix pow(Matrix x, double n) {
                return x.pow(n);
        }

        /*
         * 行列の各要素のlog() (自然対数)
         */
        Matrix log(Matrix x) {
                return x.log();
        }

        /*
         * 行列の各要素のexp()
         */
        Matrix exp(Matrix x) {
                return x.exp();
        }

        /*
         * print関数
         */
        void print(Matrix x, int n) {
                return x.print(n);
        }
} // namespace
