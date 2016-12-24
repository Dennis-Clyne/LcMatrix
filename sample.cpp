// g++ -std=c++11 sample.cpp

#include "lcmatrix.hpp"

void sample01() {
        // (1)
        printf("(1-1)\n");
        LcMatrix::Matrix zero(2, 3);    // 2行3列の0行列
        zero.print();                   // 出力

        using namespace LcMatrix;

        // (2)
        printf("(1-2)\n");
        Matrix m(2, 2, 2.5);   // 2行2列で全要素が2.5の行列
        m.print(2);            // 小数点以下2桁まで出力

        int row = m.getRow();  // mの行数を取得
        int col = m.getCol();  // mの列数を取得
        
        // (3)
        printf("(1-3)\n");
        std::cout << "row, col : ";
        std::cout << row << ", " << col << std::endl << std::endl;

        // (4)
        printf("(1-4)\n");
        Matrix x({ {1, 2, 3}, 
                   {4, 5, 6} });        // 任意の値で行列を作る
        x.print(0); 
        
        // (5)
        printf("(1-5)\n");
        double value = x.get(1, 1);     // xの2行2列目の値を取得
        std::cout << value << std::endl << std::endl;

        // (6)
        printf("(1-6)\n");
        x.set(1, 1, 7.5);        // xの2行2列目を7.5に変更
        x.print(2);
}

void sample02() {
        using namespace LcMatrix;

        Matrix m({ {1, 2, 3},
                   {4, 5, 6} });

        m.print(0);
        (-m).print(0);          // 符号を反転して出力
        m.t().print(0);         // mの転置行列を出力

        (m - 2).print(0);       // mの全要素から2を引いて出力

        // mと同じ次元数で要素がすべて2の行列から, mを引いて出力
        (Matrix(m.getRow(), m.getCol(), 2) - m).print(0);

        (m + 2).print(0);       // mの全要素から2を引いて出力

        // mと同じ次元数で要素がすべて2の行列に, mを足して出力
        (Matrix(m.getRow(), m.getCol(), 2) + m).print(0);

        pow(m, 2).print(0);     // mの全要素を2乗して出力
        exp(m).print();         // mの全要素のexpを出力
        log(m).print();         // mの全要素の自然対数を出力
}

void sample03() {
        using namespace LcMatrix;

        Matrix m1({ {1, 2},
                    {3, 4} });

        Matrix m2({ {5, 6},
                    {7, 8} });

        (m1 + m2).print(0);     // m1とm2の各要素を足して出力
        (m1 - m2).print(0);     // m1とm2の各要素を引いて出力
        (m1 * m2).print(0);     // m1とm2の各要素を掛けて出力
        (m1 / m2).print();      // m1とm2の各要素を割って出力

        dot(m1, m2).print(0);   // m1とm2の内積をとって出力
}

int main()
{
        std::cout << "\n--- sample01 ---\n" << std::endl;
        sample01();

        std::cout << "\n--- sample02 ---\n" << std::endl;
        sample02();

        std::cout << "\n--- sample03 ---\n" << std::endl;
        sample03();
}
