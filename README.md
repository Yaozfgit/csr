 # CSR Poisson Solver

用现代 C++ 从零实现的 CSR 稀疏矩阵 + poisson方程求解器，配套完整 TDD 测试与收敛阶验证。


## 项目背景

P1 结业项目：不依赖 Eigen 等现成线性代数库，手写核心数据结构与算法，
重点训练 CMake / Git / GoogleTest 工程实践与数值方法验证能力。

## 数值方法

求解 1D 泊松方程 -u'' = f， Dirichlet 边界，二阶精度中心差分（离散）：
  (-1/h²)u_{i-1} + (2/h²)u_i + (-1/h²)u_{i+1} = f(x_i)
线性系统用带部分主元的高斯消元求解。

## 项目结构

include/csr/
    csr_matrix.hpp   CSRMatrix —— CSR 格式稀疏矩阵（存储、spMV、at）
    poisson.hpp      build_poisson —— 三对角泊松矩阵组装
    linalg.hpp       Matrix（稠密矩阵）+ gauss_solve（高斯消元）
src/                 对应的 .cpp 实现
tests/               GoogleTest 测试（每个模块一个文件）


## 构建与运行

cmake -S . -B build
cmake --build build
ctest --test-dir build          # 或逐个跑 build/test_*

## 验证结果：二阶收敛

误差 max|u_h - sin(pi x)|，h 每次减半：

n=   16  err=2.839e-03  order=1.825
n=   32  err=7.547e-04  order=1.911
n=   64  err=1.946e-04  order=1.955
n=  128  err=4.942e-05  order=1.978
n=  256  err=1.245e-05  order=1.989
n=  512  err=3.125e-06  order=1.994
n= 1024  err=7.828e-07  order=1.997
n= 2048  err=1.959e-07  order=1.999

实测阶数单调逼近 2.00，单调递增，与中心差分理论精度一致。

## 学到的工程实践

- TDD 红绿循环：每个功能先写失败测试
- 分支管理：一个任务一个分支，merge 后删除
- 测试设计：分支覆盖选抽检点；每条断言须贡献新信息
- 数值验证：用收敛阶（数学性质）而非逐点比对验收代码
- 坐标只出自 interior_points()，刻度只出自 g.h()，禁手写 (i+1)/(n+1)


2026.11.09
# CSR — 一维双曲守恒律求解器（学习工程）

## 这是什么
一维有限体积求解库：线性对流方程（Rusanov 通量）、泊松方程（P3 中）。
核心类：Grid1D（均匀网格）、LinearAdvection（推进器）、rusanov_flux（数值通量）。

## 怎么构建与验收
cmake -S . -B build
cmake --build build
cd build && ctest && cd ..

## 数学性质如何被验证
- 通量三条性质（守恒形式/镜像对称/符号齐次）：test_flux
- 解的三条性质（质量守恒/整周期回位/CFL 防护）：test_advection
- 守恒性容差 1e-12（机器精度量级）；回位测试同时给出数值耗散证据
演示：./build/advect_demo 200 0.9

数值实验记录（Rusanov 基线 vs MUSCL）：见 docs/muscl_convergence.md