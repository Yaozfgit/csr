#pragma once
namespace csr {
    //界面通量：uL/uR是界面左右两侧的状态
    double rusanov_flux(double uL, double uR, double a);
}