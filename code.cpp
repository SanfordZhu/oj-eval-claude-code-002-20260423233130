#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

// Do not use "using namespace std;"

typedef unsigned int uint;
typedef long long ll;

const int BASE = 1000000000;
const int BASE_DIGITS = 9;

namespace sjtu {

struct Complex {
    long double x, y;
    Complex(long double x = 0, long double y = 0) : x(x), y(y) {}
    Complex operator+(const Complex &other) const {
        return Complex(x + other.x, y + other.y);
    }
    Complex operator-(const Complex &other) const {
        return Complex(x - other.x, y - other.y);
    }
    Complex operator*(const Complex &other) const {
        return Complex(x * other.x - y * other.y, x * other.y + y * other.x);
    }
};

void fft(std::vector<Complex> &a, bool invert) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        long double ang = 2 * 3.141592653589793L / len * (invert ? -1 : 1);
        Complex wlen(cosl(ang), sinl(ang));
        for (int i = 0; i < n; i += len) {
            Complex w(1);
            for (int j = 0; j < len / 2; j++) {
                Complex u = a[i + j];
                Complex v = a[i + j + len/2] * w;
                a[i + j] = u + v;
                a[i + j + len/2] = u - v;
                w = w * wlen;
            }
        }
    }

    if (invert) {
        for (int i = 0; i < n; i++) {
            a[i].x /= n;
            a[i].y /= n;
        }
    }
}

std::vector<ll> multiply(const std::vector<ll> &a, const std::vector<ll> &b) {
    std::vector<Complex> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < (int)(a.size() + b.size() - 1))
        n <<= 1;
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++)
        fa[i] = fa[i] * fb[i];
    fft(fa, true);

    std::vector<ll> result(n);
    ll carry = 0;
    for (int i = 0; i < n; i++) {
        ll val = (ll)(fa[i].x + 0.5);
        val += carry;
        result[i] = val % BASE;
        carry = val / BASE;
    }
    while (result.size() && result.back() == 0)
        result.pop_back();
    return result;
}

std::vector<ll> multiply_naive(const std::vector<ll> &a, const std::vector<ll> &b) {
    std::vector<ll> result(a.size() + b.size(), 0);
    for (size_t i = 0; i < a.size(); i++) {
        ll carry = 0;
        for (size_t j = 0; j < b.size() || carry; j++) {
            ll cur = result[i + j] + a[i] * (j < b.size() ? b[j] : 0) + carry;
            result[i + j] = cur % BASE;
            carry = cur / BASE;
        }
    }
    while (result.size() && result.back() == 0)
        result.pop_back();
    return result;
}

struct int2048_private {
    std::vector<ll> data;
    bool negative;

    int2048_private() : negative(false) {
        data.push_back(0);
    }
    int2048_private(ll n) {
        if (n < 0) {
            negative = true;
            n = -n;
        } else {
            negative = false;
        }
        if (n == 0) {
            data.push_back(0);
        } else {
            while (n > 0) {
                data.push_back(n % BASE);
                n /= BASE;
            }
        }
    }
    int2048_private(const std::string &s) {
        read(s);
    }
    int2048_private(const int2048_private &other) = default;

    void trim() {
        while (data.size() > 1 && data.back() == 0)
            data.pop_back();
        if (data.empty())
            data.push_back(0);
        if (is_zero())
            negative = false;
    }

    bool is_zero() const {
        return data.size() == 1 && data[0] == 0;
    }

    int compare_abs(const int2048_private &other) const {
        if (data.size() != other.data.size())
            return data.size() > other.data.size() ? 1 : -1;
        for (int i = data.size() - 1; i >= 0; i--) {
            if (data[i] != other.data[i])
                return data[i] > other.data[i] ? 1 : -1;
        }
        return 0;
    }

    bool abs_less(const int2048_private &other) const {
        return compare_abs(other) < 0;
    }

    bool abs_greater(const int2048_private &other) const {
        return compare_abs(other) > 0;
    }

    void add_abs(const int2048_private &other) {
        ll carry = 0;
        for (size_t i = 0; i < std::max(data.size(), other.data.size()) || carry; i++) {
            if (i < data.size())
                carry += data[i];
            if (i < other.data.size())
                carry += other.data[i];
            if (i < data.size())
                data[i] = carry % BASE;
            else
                data.push_back(carry % BASE);
            carry /= BASE;
        }
        trim();
    }

    void sub_abs(const int2048_private &other) {
        ll carry = 0;
        for (size_t i = 0; i < other.data.size() || carry; i++) {
            data[i] -= carry;
            if (i < other.data.size())
                data[i] -= other.data[i];
            if (data[i] < 0) {
                data[i] += BASE;
                carry = 1;
            } else {
                carry = 0;
            }
        }
        trim();
    }

    void mul_abs(const int2048_private &other) {
        if (is_zero() || other.is_zero()) {
            data.clear();
            data.push_back(0);
            negative = false;
            return;
        }
        if (data.size() > 100 || other.data.size() > 100) {
            data = multiply(data, other.data);
        } else {
            data = multiply_naive(data, other.data);
        }
        trim();
    }

    ll mul_single(ll x) {
        ll carry = 0;
        for (size_t i = 0; i < data.size(); i++) {
            ll product = data[i] * x + carry;
            data[i] = product % BASE;
            carry = product / BASE;
        }
        if (carry > 0)
            data.push_back(carry);
        trim();
        return carry;
    }

    void read(const std::string &s) {
        data.clear();
        if (s.empty()) {
            negative = false;
            data.push_back(0);
            return;
        }
        int start = 0;
        if (s[0] == '-') {
            negative = true;
            start = 1;
        } else if (s[0] == '+') {
            negative = false;
            start = 1;
        } else {
            negative = false;
            start = 0;
        }
        for (int i = s.size() - 1; i >= start; i -= BASE_DIGITS) {
            int begin = std::max(start, i - BASE_DIGITS + 1);
            std::string num = s.substr(begin, i - begin + 1);
            data.push_back(stoll(num));
        }
        trim();
    }

    std::string to_string() const {
        if (is_zero())
            return "0";
        std::string s;
        if (negative)
            s += '-';
        bool first = true;
        for (int i = data.size() - 1; i >= 0; i--) {
            if (first) {
                char buf[20];
                snprintf(buf, sizeof(buf), "%lld", data[i]);
                s += buf;
                first = false;
            } else {
                char buf[20];
                snprintf(buf, sizeof(buf), "%09lld", data[i]);
                s += buf;
            }
        }
        return s;
    }
};

static int2048_private div_approx(const int2048_private &a, const int2048_private &b) {
    if (a.is_zero())
        return int2048_private(0);
    int2048_private res;
    res.negative = a.negative != b.negative;
    int n = a.data.size();
    int m = b.data.size();
    if (n < m) {
        res.data.push_back(0);
        return res;
    }
    ll highest_a = a.data.back();
    ll highest_b = b.data.back();
    ll approx = highest_a / (highest_b + 1);
    if (approx == 0) approx = 1;
    res.data.push_back(approx);
    return res;
}

static int2048_private divide(int2048_private a, const int2048_private &b);

static int2048_private divide(int2048_private a, const int2048_private &b) {
    if (a.is_zero())
        return int2048_private(0);

    bool a_neg = a.negative;
    bool b_neg = b.negative;
    int2048_private a_abs = a;
    int2048_private b_abs = b;
    a_abs.negative = false;
    b_abs.negative = false;
    int2048_private one(1);

    // Long division: build quotient one digit at a time from MSB to LSB
    int2048_private quotient;
    quotient.data.resize(a_abs.data.size());
    int2048_private current;

    for (int i = a_abs.data.size() - 1; i >= 0; i--) {
        current.data.insert(current.data.begin(), a_abs.data[i]);
        current.trim();
        ll low = 0, high = BASE - 1;
        ll best = 0;
        // Binary search for the largest digit q where b_abs * q <= current
        while (low <= high) {
            ll mid = low + (high - low) / 2;
            int2048_private tmp = b_abs;
            tmp.mul_single(mid);
            if (tmp.compare_abs(current) <= 0) {
                best = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        quotient.data[i] = best;
        int2048_private product = b_abs;
        product.mul_single(best);
        current.sub_abs(product);
    }

    quotient.trim();
    int2048_private q_abs = quotient;

    // Check if exact division
    int2048_private product = q_abs;
    product.mul_abs(b_abs);
    bool has_remainder = (product.compare_abs(a_abs) != 0);

    int2048_private result;
    if (a_neg == b_neg) {
        result = q_abs;
        result.negative = false;
    } else {
        if (has_remainder) {
            // floor(-x) = -ceil(x) = -(q_abs + 1)
            q_abs.add_abs(one);
        }
        result = q_abs;
        result.negative = true;
    }
    result.trim();
    return result;
}

static int2048_private mod(const int2048_private &a, const int2048_private &b) {
    int2048_private div_result = divide(a, b);
    int2048_private product = div_result;
    product.mul_abs(b);
    product.negative = (div_result.negative != b.negative);

    // x % y = x - (x/y) * y
    int2048_private result;
    if (a.negative == product.negative) {
        // same sign: a - product → subtract absolute values
        if (a.compare_abs(product) >= 0) {
            result = a;
            result.sub_abs(product);
        } else {
            result = product;
            result.sub_abs(a);
            if (!result.is_zero())
                result.negative = !a.negative;
        }
    } else {
        // different signs: a - product = a + (-product) → add absolute values
        result = a;
        result.add_abs(product);
    }
    result.trim();
    return result;
}

class int2048 {
  int2048_private data;
public:
  // Constructors
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &);
  // Output the stored big integer, no need for newline
  void print();

  // Add a big integer
  int2048 &add(const int2048 &);
  // Return the sum of two big integers
  friend int2048 add(int2048, const int2048 &);

  // Subtract a big integer
  int2048 &minus(const int2048 &);
  // Return the difference of two big integers
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};

// Implementation

int2048::int2048() : data() {}

int2048::int2048(long long n) : data(n) {}

int2048::int2048(const std::string &s) : data(s) {}

int2048::int2048(const int2048 &other) : data(other.data) {}

// Integer1 methods

void int2048::read(const std::string &s) {
    data.read(s);
}

void int2048::print() {
    if (data.is_zero()) {
        std::cout << "0";
        return;
    }
    if (data.negative)
        std::cout << "-";
    bool first = true;
    for (int i = data.data.size() - 1; i >= 0; i--) {
        if (first) {
            printf("%lld", data.data[i]);
            first = false;
        } else {
            printf("%09lld", data.data[i]);
        }
    }
}

int2048 &int2048::add(const int2048 &other) {
    if (data.negative == other.data.negative) {
        data.add_abs(other.data);
    } else {
        if (data.compare_abs(other.data) >= 0) {
            data.sub_abs(other.data);
        } else {
            int2048_private tmp = other.data;
            tmp.sub_abs(data);
            data = tmp;
            // Keep the sign from other (tmp already has the correct sign from other)
        }
    }
    data.trim();
    return *this;
}

int2048 add(int2048 a, const int2048 &b) {
    return a.add(b);
}

int2048 &int2048::minus(const int2048 &other) {
    if (data.negative != other.data.negative) {
        data.add_abs(other.data);
    } else {
        if (data.compare_abs(other.data) >= 0) {
            data.sub_abs(other.data);
        } else {
            int2048_private tmp = other.data;
            tmp.sub_abs(data);
            data = tmp;
            data.negative = !data.negative;
        }
    }
    data.trim();
    return *this;
}

int2048 minus(int2048 a, const int2048 &b) {
    return a.minus(b);
}

// Integer2 methods

int2048 int2048::operator+() const {
    return int2048(*this);
}

int2048 int2048::operator-() const {
    int2048 result(*this);
    if (!result.data.is_zero())
        result.data.negative = !result.data.negative;
    return result;
}

int2048 &int2048::operator=(const int2048 &other) {
    if (this != &other) {
        data = other.data;
    }
    return *this;
}

int2048 &int2048::operator+=(const int2048 &other) {
    return add(other);
}

int2048 operator+(int2048 a, const int2048 &b) {
    return a.add(b);
}

int2048 &int2048::operator-=(const int2048 &other) {
    return minus(other);
}

int2048 operator-(int2048 a, const int2048 &b) {
    return a.minus(b);
}

int2048 &int2048::operator*=(const int2048 &other) {
    data.mul_abs(other.data);
    data.negative = data.negative != other.data.negative;
    data.trim();
    return *this;
}

int2048 operator*(int2048 a, const int2048 &b) {
    a *= b;
    return a;
}

int2048 &int2048::operator/=(const int2048 &other) {
    data = divide(data, other.data);
    data.trim();
    return *this;
}

int2048 operator/(int2048 a, const int2048 &b) {
    a /= b;
    return a;
}

int2048 &int2048::operator%=(const int2048 &other) {
    data = mod(data, other.data);
    data.trim();
    return *this;
}

int2048 operator%(int2048 a, const int2048 &b) {
    a %= b;
    return a;
}

std::istream &operator>>(std::istream &is, int2048 &num) {
    std::string s;
    is >> s;
    num.read(s);
    return is;
}

std::ostream &operator<<(std::ostream &os, const int2048 &num) {
    if (num.data.is_zero()) {
        os << "0";
        return os;
    }
    if (num.data.negative)
        os << "-";
    bool first = true;
    for (int i = num.data.data.size() - 1; i >= 0; i--) {
        if (first) {
            os << num.data.data[i];
            first = false;
        } else {
            char buf[20];
            snprintf(buf, sizeof(buf), "%09lld", num.data.data[i]);
            os << buf;
        }
    }
    return os;
}

bool operator==(const int2048 &a, const int2048 &b) {
    if (a.data.negative != b.data.negative)
        return false;
    return a.data.compare_abs(b.data) == 0;
}

bool operator!=(const int2048 &a, const int2048 &b) {
    return !(a == b);
}

bool operator<(const int2048 &a, const int2048 &b) {
    if (a.data.negative != b.data.negative)
        return a.data.negative;
    if (a.data.negative)
        return a.data.compare_abs(b.data) > 0;
    else
        return a.data.compare_abs(b.data) < 0;
}

bool operator>(const int2048 &a, const int2048 &b) {
    return b < a;
}

bool operator<=(const int2048 &a, const int2048 &b) {
    return !(a > b);
}

bool operator>=(const int2048 &a, const int2048 &b) {
    return !(a < b);
}

} // namespace sjtu

#endif
