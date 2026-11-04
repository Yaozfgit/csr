#include "csr/csr_matrix.hpp"

CSRMatrix::CSRMatrix(std::size_t n) : n_(n), current_row_(0), row_ptr_{0} {}
std::size_t CSRMatrix::row_end(std::size_t i) const {
        return (i+1< row_ptr_.size()) ? row_ptr_[i+1] : values_.size();
    }
std::vector<double> CSRMatrix::multiply(const std::vector<double>& x) const {
    if (x.size() !=n_) throw std::invalid_argument("size mismatch");
    std::vector<double> y(n_,0.0);     
    for (std::size_t i=0; i < n_; ++i) {
        double sum = 0.0;        
        std::size_t begin = row_ptr_[i];
        std::size_t end = row_end (i);
        for (std::size_t k = begin; k<end; ++k)
        sum +=values_[k] * x[col_idx_[k]];
    y[i] = sum;
    }
    return y; 
}

std::size_t CSRMatrix::nnz() const { return values_.size(); }

std::size_t CSRMatrix::row_ptr_at(std::size_t k) const {
    if (k>=row_ptr_.size()) throw std::out_of_range("row_ptr");
    return row_ptr_[k];
}

std::size_t CSRMatrix::col_idx_at(std::size_t k) const {
    if (k>=col_idx_.size()) throw std::out_of_range("col_idx");
    return col_idx_[k];
}

double CSRMatrix::values_at(std::size_t k) const {
    if (k>=values_.size()) throw std::out_of_range("values");
    return values_[k];
}

double CSRMatrix::at(std::size_t i, std::size_t j) const {
    if (i >= n_) throw std::out_of_range ("row index");
    for (std::size_t k=row_ptr_[i]; k<row_end(i); ++k) {
        if (col_idx_[k]==j) return values_[k];
    }
    return 0.0;
}

void CSRMatrix::insert(std::size_t i, std::size_t j, double v) {
    if (i < current_row_ || j >= n_)
       throw std::invalid_argument("insert must be in row order");
    while (current_row_ < i) {
        row_ptr_.push_back(values_.size());
        ++current_row_;
    }
    col_idx_.push_back(j);
    values_.push_back(v);
}