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

namespace LcMatrix {
        /*
        std::vector<std::vector<double>> Matrix::matrix;
        double Matrix::row, Matrix::col;
        */

        void Matrix::setRC() {
                row = matrix.size();
                col = matrix[0].size();
        }

        /*
         * 宣言だけするとき
         */
        Matrix::Matrix() {
        }

        /*
         * 2次元vectorの行列のmatrixを作る
         */
        Matrix::Matrix(std::vector<std::vector<double>> x) {
                matrix = x;
                setRC();
        }

        /*
         * i行j列のmatrixをxで初期化して作る
         * xはデフォルト引数でその値は0.0
         */
        Matrix::Matrix(int i, int j, double x) {
                matrix = std::vector<std::vector<double>>(i, std::vector<double>(j, x));
                setRC();
        }

        /*
         * matrix[i][j]を返す
         */
        double Matrix::get(int i, int j) {
                return matrix[i][j];
        }

        void Matrix::set(int i, int j, double x) {
                matrix[i][j] = x;
        }

        /*
         * matrix自体を返す
         */
        std::vector<std::vector<double>> Matrix::getMatrix() {
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
         * 行列内で最も大きな値を返す
         */
        double Matrix::max() {
                std::vector<double> max_vec;
                for (auto i : matrix) {
                        max_vec.push_back(*std::max_element(i.begin(), i.end()));
                }

                return *std::max_element(max_vec.begin(), max_vec.end());
        }

        /*
         * max関数(bool)
         * 行または列内で最大の要素を行列として返す.
         * true : 行
         * false : 列
         */
        Matrix Matrix::max(bool x) {
                // 行方向
                if (x) {
                        Matrix result(row, 1);
                        for (int i = 0; i < row; i++) {
                                result.set(i, 0, *std::max_element(matrix[i].begin(), matrix[i].end()));
                        }
                        return result;
                } 
                // 列方向
                else {
                        Matrix result(1, col);
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        if (result.get(0, j) < matrix[i][j])
                                                result.set(0, j, matrix[i][j]);
                                }
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
                // 行方向
                if (x) {
                        Matrix result(row, 1);
                        for (int i = 0; i < row; i++) {
                                result.set(i, 0, (int)std::distance(matrix[i].begin(), std::max_element(matrix[i].begin(), matrix[i].end())));
                        }
                        return result;
                } 
                // 列方向
                else {
                        Matrix result(1, col);
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        if (result.get(0, j) < matrix[i][j])
                                                result.set(0, j, i);
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
                double sum = 0;
                for (auto i : matrix) {
                        sum += std::accumulate(i.begin(), i.end(), 0.0);
                }
                return sum;
        }

        /*
         * 行列の列または行方向の和
         * true : 行方向
         * false : 列方向
         */
        Matrix Matrix::sum(bool x) {
                // 行方向
                if (x) {
                        Matrix result(row, 1);
                        for (int i = 0; i < row; i++) {
                                result.set(i, 0, std::accumulate(matrix[i].begin(), matrix[i].end(), 0.0));
                        }
                        return result;
                } 
                // 列方向
                else {
                        Matrix result(1, col);
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        result.set(0, j, result.get(0, j) + matrix[i][j]);
                                }
                        }
                        return result;
                }
        }

        /*
         * 行列の内積
         */
        Matrix Matrix::dot(Matrix x) {
                Matrix result(row, x.getCol());

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < x.getCol(); j++) {
                                for (int z = 0; z < col; z++) {
                                        result.matrix[i][j] += (matrix[i][z] * x.matrix[z][j]);
                                }
                        }
                }

                return result;
        }

        /*
         * 行列の転置
         */
        Matrix Matrix::t() {
                Matrix result(col, row);
                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[j][i] = matrix[i][j];
                        }
                }
                return result;
        }

        /*
         * 行列の要素どうしの掛け算 (行列 * 行列)
         */
        Matrix Matrix::operator * (Matrix pre_x) {
                Matrix result(row, col);

                // 行と列が同数のとき
                if (row == pre_x.getRow() && col == pre_x.getCol()) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        result.matrix[i][j] = matrix[i][j] * pre_x.matrix[i][j];
                                }
                        }
                        return result;
                }

                Matrix x(row, col);

                // numpyのブロードキャストみたいなやつ
                // 行数同じで列が1つ
                if (row == pre_x.getRow() && pre_x.getCol() == 1) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        x.set(i, j, pre_x.get(i, 0));
                                }
                        }
                }
                // 列数同じで行が1つ
                if (col == pre_x.getCol() && pre_x.getRow() == 1) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        x.set(i, j, pre_x.get(0, j));
                                }
                        }
                }

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = matrix[i][j] * x.matrix[i][j];
                        }
                }

                return result;
        }

        /*
         * 行列の各要素とdoubleの掛け算 (行列 * double)
         */
        Matrix Matrix::operator * (double x) {
                Matrix result(row, col);

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = matrix[i][j] * x;
                        }
                }

                return result;
        }


        /*
         * 行列の和 (行列 + 行列)
         */
        Matrix Matrix::operator + (Matrix pre_x) {
                Matrix result(row, col);

                // 行と列が同数のとき
                if (row == pre_x.getRow() && col == pre_x.getCol()) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        result.matrix[i][j] = matrix[i][j] + pre_x.matrix[i][j];
                                }
                        }
                        return result;
                }

                Matrix x(row, col);

                // numpyのブロードキャストみたいなやつ
                // 行数同じで列が1つ
                if (row == pre_x.getRow() && pre_x.getCol() == 1) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        x.set(i, j, pre_x.get(i, 0));
                                }
                        }
                }
                // 列数同じで行が1つ
                if (col == pre_x.getCol() && pre_x.getRow() == 1) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        x.set(i, j, pre_x.get(0, j));
                                }
                        }
                }

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = matrix[i][j] + x.matrix[i][j];
                        }
                }
                return result;
        }

        /*
         * 行列の和 (行列 + double)
         */
        Matrix Matrix::operator + (double x) {
                Matrix result(row, col);

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = matrix[i][j] + x;
                        }
                }
                return result;
        }

        /*
         * 行列の差 (行列 - 行列)
         */
        Matrix Matrix::operator - (Matrix pre_x) {
                Matrix result(row, col);

                // 行と列が同数のとき
                if (row == pre_x.getRow() && col == pre_x.getCol()) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        result.matrix[i][j] = matrix[i][j] - pre_x.matrix[i][j];
                                }
                        }
                        return result;
                }

                Matrix x(row, col);

                // numpyのブロードキャストみたいなやつ
                // 行数同じで列が1つ
                if (row == pre_x.getRow() && pre_x.getCol() == 1) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        x.set(i, j, pre_x.get(i, 0));
                                }
                        }
                }
                // 列数同じで行が1つ
                if (col == pre_x.getCol() && pre_x.getRow() == 1) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        x.set(i, j, pre_x.get(0, j));
                                }
                        }
                }

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = matrix[i][j] - x.matrix[i][j];
                        }
                }
                return result;
        }

        /*
         * 行列の差 (行列 - double)
         */
        Matrix Matrix::operator - (double x) {
                Matrix result(row, col);

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = matrix[i][j] - x;
                        }
                }
                return result;
        }

        /*
         * 行列の各要素の符号反転
         */
        Matrix Matrix::operator - () {
                Matrix result(row, col);

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = -matrix[i][j];
                        }
                }
                return result;
        }

        /*
         * 行列の各要素どうしの割り算 (行列 / 行列)
         */
        Matrix Matrix::operator / (Matrix pre_x) {
                Matrix result(row, col);

                // 行と列が同数のとき
                if (row == pre_x.getRow() && col == pre_x.getCol()) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        result.matrix[i][j] = matrix[i][j] / pre_x.matrix[i][j];
                                }
                        }
                        return result;
                }

                Matrix x(row, col);

                // numpyのブロードキャストみたいなやつ
                // 行数同じで列が1つ
                if (row == pre_x.getRow() && pre_x.getCol() == 1) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        x.set(i, j, pre_x.get(i, 0));
                                }
                        }
                }
                // 列数同じで行が1つ
                if (col == pre_x.getCol() && pre_x.getRow() == 1) {
                        for (int i = 0; i < row; i++) {
                                for (int j = 0; j < col; j++) {
                                        x.set(i, j, pre_x.get(0, j));
                                }
                        }
                }

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = matrix[i][j] / x.matrix[i][j];
                        }
                }
                return result;
        }

        /*
         * 行列の割り算 (行列 / double)
         */
        Matrix Matrix::operator / (double x) {
                Matrix result(row, col);

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = matrix[i][j] / x;
                        }
                }
                return result;
        }

        /*
         * 行列の各要素のn乗
         */
        Matrix Matrix::pow(double n) {
                Matrix result(row, col);

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = std::pow(matrix[i][j], n);
                        }
                }
                return result;
        }

        /*
         * 行列の各要素のlog() (自然対数)
         */
        Matrix Matrix::log() {
                Matrix result(row, col);

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = std::log(matrix[i][j]);
                        }
                }
                return result;
        }

        /*
         * 行列の各要素のexp()
         */
        Matrix Matrix::exp() {
                Matrix result(row, col);

                for (int i = 0; i < row; i++) {
                        for (int j = 0; j < col; j++) {
                                result.matrix[i][j] = std::exp(matrix[i][j]);
                        }
                }
                return result;
        }

        /*
         * 行列の表示
         * 小数点以下n桁まで表示(デフォルトはn = 5).
         */
        void Matrix::print(int n) {
                for (auto i : matrix) {
                        for (double j : i) {
                                //printf("%.15lf, ", j);
                                std::cout << std::fixed << std::setprecision(n) << j << ", ";
                        }
                        std::cout << std::endl;
                }
                std::cout << std::endl;
        }
} // namespace

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
