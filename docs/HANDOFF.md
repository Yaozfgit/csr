# 项目交接快照 — csr（1D CFD 求解器，C++/CMake/gtest）

## 架构
- src/: matrix, csr_matrix, poisson, linalg, grid, flux, advect1d, burgers
- include/: 头文件；thomas.h 为新增的 O(n) 三对角求解器（约定：
  a[i]·x[i-1], d[i]·x[i], c[i]·x[i+1]，a[0] 与 c[n-1] 不参与运算，置零属防御性）
- tests/: 每模块一个 gtest；test_thomas.cpp 含 Tiny2x2 + MatchesDenseSolver，
  共用 assemble_b() 组装右端项
- CMakeLists：FetchContent googletest v1.14.0，enable_testing + add_test 全部就位

## 已完成
- Thomas 算法实现与对拍测试全绿（2025-xx-xx）
- 根因复盘：测试中 c[0] 被误置零（应为 c[n-1]）+ 两个用例组装逻辑不一致

## 当前任务：BurgersSolver 集成收尾
1. src/burgers.cpp step()：删除 CSRMatrix M 死代码；
   三对角系数 diag 边界行必须为 1.0+3.0*r（复刻原 CSR 矩阵的边界离散），
   内部 1.0+2.0*r，lower=-r/upper=-r，lower[0]=upper[n-1]=0
2. 全量 ctest 验收：DiffusionConvergence/SelfConvergence 误差应与
   替换前一致（~1e-8 档），总耗时下降
3. git commit

## 约定与坑（重要）
- 诊断方法：用报错数字反推中间量定位（指纹法），测试先断言 b 再断言 x
- ctest 需在 build 目录下运行；改 CMakeLists 后必须重新 cmake configure
- 曾踩坑：改代码后忘记重编导致旧二进制误报；双测试用例修一漏一
