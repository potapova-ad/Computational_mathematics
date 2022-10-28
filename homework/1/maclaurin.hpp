#include <matplot/matplot.h>
#include <functional>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>

using function_type = std::function<double(const double)>;

int fact(int n)
{
    int res = 1;
    for (int i = 2; i <= n; ++i)
        res *= i;
    return res;
}

int C(int n, int k)
{
    return fact(n) / fact(k) / fact(n - k);
}


class derivative
{
public:
    derivative(const function_type &f) : _f(f){};
    virtual double value(const int n, const double x, const double dh) = 0;
    int size() { return _val.size(); }
    virtual ~derivative() {}

protected:
    const function_type &_f;
    std::vector<std::pair<double, std::map<int, double>>> _val;
    int __find_it_key(const double dh) 
    {
    for (int i = 0; i < _val.size(); ++i)
        if (std::abs(_val[i].first - dh) <= 1e-6)
        {
            return i;
        }
    return -1;
}

    
};

class central_der final : public derivative
{
public:
    central_der(const function_type &f) : derivative(f) {}
    double value(const int n, const double x, const double dh) override
    {
    int key = __find_it_key(dh);
    if (key != -1 && _val[key].second.size() > n)
        return _val[key].second[n];
    if (key == -1)
    {
        key = _val.size();
        _val.push_back({dh, {}});
    }

    if (_val[key].second.size() == 0)
    {
        _val[key].second[0] = _f(x);
        return _val[key].second[0];
    }

    double t_sum = 0;
    for (int i = 0; i <= n; ++i)
    {
        t_sum += double(C(n, i)) * ((i % 2) ? -1 : 1) * _f(x + (n - double(2 * i)) * dh);
    }

    _val[key].second[n] = t_sum / std::pow(2 * dh, n);
    return _val[key].second[n];
}
};


class back_der final : public derivative
{
public:
    back_der(const function_type &f) : derivative(f) {}
    double value(const int n, const double x, const double dh) override
    {
    int key = __find_it_key(dh);
    if (key != -1 && _val[key].second.size() > n)
        return _val[key].second[n];
    if (key == -1)
    {
        key = _val.size();
        _val.push_back({dh, {}});
    }

    if (_val[key].second.size() == 0)
    {
        _val[key].second[0] = _f(x);
        return _val[key].second[0];
    }

    double t_sum = 0;
    for (int i = 0; i <= n; ++i)
    {
        t_sum += double(C(n, i)) * ((i % 2) ? -1 : 1) * _f(x - (double(i)) * dh);
    }
    _val[key].second[n] = t_sum / std::pow(dh, n);
    return _val[key].second[n];
}
};

class forward_der final : public derivative
{
public:
    forward_der(const function_type &f) : derivative(f) {}
    double value(const int n, const double x, const double dh) override
    {
    int key = __find_it_key(dh);
    if (key != -1 && _val[key].second.size() > n)
        return _val[key].second[n];
    if (key == -1)
    {
        key = _val.size();
        _val.push_back({dh, {}});
    }

    if (_val[key].second.size() == 0)
    {
        _val[key].second[0] = _f(x);
        return _val[key].second[0];
    }

    double t_sum = 0;
    for (int i = 0; i <= n; ++i)
    {
        t_sum += double(C(n, i)) * ((i % 2) ? -1 : 1) * _f(x + (n - double(i)) * dh);
    }
    _val[key].second[n] = t_sum / std::pow(dh, n);
    return _val[key].second[n];
}
};

class Maclaurin final
{
private:
    derivative *_der = nullptr;
    std::vector<std::pair<double, std::vector<double>>> _val = {};

    int __find_key(const double x)
    {
    for (int i = 0; i < _val.size(); ++i)
        if (std::abs(_val[i].first - x) <= 1e-6)
        {
            return i;
        }
    return -1;
}

public:

    Maclaurin(const function_type &f);
    int size() { return _der->size(); }

    double value(const int n, const double x, const double dt)
{
    int key = __find_key(x);
    if (key != -1 && _val[key].second.size() > n)
        return _val[key].second[n];
    if (key == -1)
    {
        key = _val.size();
        _val.push_back({x, {}});
    }

    if (_val[key].second.size() == 0)
        _val[key].second.push_back(_der->value(0, 0, dt));
    for (int k = _val[key].second.size(); k <= n; ++k)
    {
        _val[key].second.push_back(_val[key].second[k - 1] + std::pow(x, k) * _der->value(k, 0, dt) / fact(k));
    }
    return _val[key].second[n];
}
};







