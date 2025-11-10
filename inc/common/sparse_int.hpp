#ifndef SPARSE_INT_HPP
#define SPARSE_INT_HPP

#include <set>
#include <string>

class SparseInt {
    private:
        long long msb;
        long long lsb;
        long long offset;
        std::set<unsigned long long> indices;
        void insert_index(unsigned long long idx);

    public:
        SparseInt();
        SparseInt(unsigned long long n);
        SparseInt operator+(SparseInt other);
        SparseInt operator*(SparseInt other);
        void operator+=(SparseInt other);
        SparseInt operator<<(long long shift);
        SparseInt operator>>(long long shift);
        std::string to_string();
};

#endif