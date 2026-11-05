#pragma once
#include <vector>
#include <cstddef>

class CSRMatrix {
public:
    CSRMatrix(std::size_t n);
    void insert(std::size_t i, std::size_t j, double v);
    std::vector<double> multiply(const std::vector<double>& x) const;    
    std::size_t nnz() const;

    std::size_t row_ptr_at(std::size_t k) const;
    std::size_t col_idx_at(std::size_t k) const;
    double values_at(std::size_t k) const;
    std::size_t row_end(std::size_t i) const;

    double at(std::size_t i, std::size_t j) const;    
         
private:
    std::size_t n_;
    std::size_t current_row_ = 0;
    std::vector<std::size_t> row_ptr_ = {0};
    std::vector<std::size_t> col_idx_;
    std::vector<double> values_;
};