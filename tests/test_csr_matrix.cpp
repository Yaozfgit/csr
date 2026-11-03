#include <gtest/gtest.h>
#include <csr/csr_matrix.hpp>
#include <csr/matrix.hpp>
#include <vector>
#include <cmath>

//------1.空矩阵：全零时SpMV结果必须是全灵----
TEST(CSRMatrixTest, EmptyMatrix) {
    CSRMatrix m(3);
    auto y = m.multiply(std::vector<double>{1.0, 2.0, 3.0});
    ASSERT_EQ(y.size(), 3u);
    EXPECT_DOUBLE_EQ(y[0], 0.0);
    EXPECT_DOUBLE_EQ(y[1], 0.0);
    EXPECT_DOUBLE_EQ(y[2], 0.0);
    EXPECT_EQ(m.nnz(), 0u);    
}

//-----2. 教科书4x4：手工验证三个数组的精确内容----
//直接检查内部表示是否就是CSR定义本身——比只看结果更严格
TEST(CSRMatrixTest, DataLayout) {
    CSRMatrix m(4);
    m.insert(0,0,10.0); m.insert(0,3,12.0);
    m.insert(1,1,20.0);
    m.insert(2,2,30.0);
    m.insert(3,0,14.0); m.insert(3,3,15.0);

    EXPECT_EQ(m.nnz(), 6u);
    //row_ptr = [0,2,3,4,6],包含隐式的最后一行行尾
    std::vector<std::size_t> expect_rp = {0,2,3,4,6};
    std::vector<std::size_t> expect_ci = {0,3,1,2,0,3};
    std::vector<double> expect_v = {10,12,20,30,14,15};
    for (int i = 0; i<4; ++i)
        EXPECT_EQ(m.row_ptr_at(i), expect_rp[i])<<"row_ptr["<<i<<"]";
    for (int k = 0; k< 4; ++k) {
        EXPECT_EQ(m.col_idx_at(k), expect_ci[k])<<"col_idx["<<k<<"]";
        EXPECT_DOUBLE_EQ(m.values_at(k), expect_v[k])<<"values["<<k<<"]";
    }

}

//------3.SpMV vs 手算-----
TEST(CSRMatrixTest, SpMVHand) {
    CSRMatrix m(4);
    m.insert(0,0,10);m.insert(0,3,12.0);
    m.insert(1,1,20.0);
    m.insert(2,2,30.0);
    m.insert(3,0,14.0); m.insert(3,3,15.0);
    auto y = m.multiply(std::vector<double>{1.0, 2.0, 3.0, 4.0});
    EXPECT_DOUBLE_EQ(y[0], 10*1+12*4); //58
    EXPECT_DOUBLE_EQ(y[1], 20*2); //40
    EXPECT_DOUBLE_EQ(y[2], 30*3); //90
    EXPECT_DOUBLE_EQ(y[3], 14*1+15*4); //74
}

//---4.空行：第2行全零，不能错位——CSR bug高发区---
TEST(CSRMatrixTest, EmptyRows) {
    CSRMatrix m(3);
    m.insert(0,1,5.0);  //第0行只有一个非零元
    //第一行整行时零
    m.insert(2,0,7.0);   //这是第二行
    auto y = m.multiply(std::vector<double>{1.0,1.0,1.0});
    EXPECT_DOUBLE_EQ(y[0],5.0);
    EXPECT_DOUBLE_EQ(y[1],0.0);  //空行必须精确为0
    EXPECT_DOUBLE_EQ(y[2],7.0);
}

//----5. 1x1 边界（前面做过）-----
TEST(CSRMatrixTest, OneByOne) {
    CSRMatrix m(1);
    m.insert(0,0,3.0);
    auto y = m.multiply(std::vector<double>{2.0});
    EXPECT_DOUBLE_EQ(y[0],6.0);
}

//----6. 对拍：CSR与dense结果必须逐元素一致（重要）-----
// 构造32x32的三对角矩阵，同一向量分别用两种实现算 SpMV
TEST(CSRMatrixTest, AgainstDense) {
    const std::size_t n =32;
    CSRMatrix s(n);
    Matrix d(n);
    std::vector<double> x(n);
    for (std::size_t i = 0; i<n; ++i) {
        x[i] = 1.0 +i*0.1;
        if (i>0) {
            s.insert(i,i-1,1.0); d(i,i-1) =1.0;
        }
        s.insert(i,i,-2.0); d(i,i)=-2.0;
        if (i+1<n) {
            s.insert(i,i+1,1.0); d(i,i+1) =1.0;
        }
    }
    auto ys = s.multiply(x);
    auto yd = d.multiply(x);
    ASSERT_EQ(ys.size(), yd.size());
    for (std::size_t i = 0; i<n; ++i)
    EXPECT_DOUBLE_EQ(ys[i],yd[i]) <<"row "<<i;
}

//---7. insert 乱序被拒绝（保护“按行序”的约定）----
TEST(CSRMatrixTest, RejectOutOfOrderInsert) {
    CSRMatrix m(3);
    m.insert(1,0,1.0);
    EXPECT_THROW(m.insert(0,2,5.0), std::invalid_argument);//回头：拒绝
}
    