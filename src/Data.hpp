#pragma once

#include <cstddef>
#include <vector>

union kTableField {
    long long int i64;
    double fp64;

    kTableField() {
        i64 = 0;
        fp64 = 0.l;
    }
    kTableField(const long long& i) {
        i64 = i;
    }
    kTableField(const double& d) {
        fp64 = d;
    }
};

class Table {
    std::vector<std::vector<kTableField>> _data;

public:
    Table(const size_t& m, const size_t& n) {
        _data.resize(m);
        for (auto& r : _data) {
            r.resize(n);
        }
    }

    kTableField at(const size_t& i, const size_t& j) {
        return _data[i][j];
    }

    kTableField& operator()(const size_t& i, const size_t& j) { 
        return _data[i][j];
    }

    std::vector<kTableField>& operator[](const size_t& i) {
        return _data[i];
    }
};