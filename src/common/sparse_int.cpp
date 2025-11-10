#include "common/sparse_int.hpp"

#include <cassert>
#include <cmath>
#include <ranges>

SparseInt::SparseInt() {
    msb = -1;
    lsb = 0;
    offset = 0;
}

SparseInt::SparseInt(unsigned long long n) {
    msb = -1;
    lsb = 0;
    offset = 0;

    if (n == 0) {
        return;
    }

    long long idx = 0;
    while (n > 0) {
        if (n % 2) {
            indices.insert(idx);
            if (msb == -1) {
                lsb = idx;
            }
            msb = idx;
        }
        n >>= 1;
        idx++;
    }
    return;
}

void SparseInt::insert_index(unsigned long long idx) {
    if (indices.find(idx) == indices.end()) {
        indices.insert(idx);
        if (idx > msb) {
            msb = idx;
        }
        if (idx < lsb) {
            lsb = idx;
        }
    } else {
        indices.erase(idx);
        if (idx == lsb) {
            lsb++;
        }
        insert_index(idx + 1);
    }
}

SparseInt SparseInt::operator+(SparseInt other) {
    if (msb == -1) {
        return other;
    }
    if (other.msb == -1) {
        return *this;
    }
    SparseInt result;

    result.offset = 0;
    for (auto idx : indices) {
        result.indices.insert(idx + offset);
    }
    result.lsb = lsb + offset;
    result.msb = msb + offset;

    for (auto idx : other.indices) {
        result.insert_index(idx + other.offset);
    }
    return result;
}

SparseInt SparseInt::operator*(SparseInt other) {
    if (msb == -1) {
        return *this;
    }
    if (other.msb == -1) {
        return other;
    }
    SparseInt result;

    for (auto & b: indices) {
        result += other << (b + offset);
    }
    return result;
}

void SparseInt::operator+=(SparseInt other) {
    if (msb == -1) {
        *this = other;
        return;
    }
    if (other.msb == -1) {
        return;
    }
    *this = *this + other;
    return;
}

SparseInt SparseInt::operator<<(long long shift) {
    if (msb == -1) {
        return *this;
    }
    SparseInt result(*this);
    result.offset += shift;
    assert(result.offset + lsb >= 0);
    return result;
}

SparseInt SparseInt::operator>>(long long shift) {
    if (msb == -1) {
        return *this;
    }
    SparseInt result(*this);
    result.offset -= shift;
    assert(result.offset + lsb >= 0);
    return result;
}

std::string SparseInt::to_string() {
    std::string result = "[";
    for (auto& idx : std::views::reverse(indices)) {
        result += std::to_string(idx + offset) + (idx == lsb ? "" : " ");
    }
    result += "]";
    return result;
}