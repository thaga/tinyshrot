#include "SH.h"

#define TINYSHROT_IMPLEMENTATION
#include "../tinyshrot.h"

#include <cmath>
#include <cassert>
#include <algorithm>


namespace {

constexpr float π = 3.14159265f;

// 階乗
constexpr int factorial(int n) { return n >= 2 ? n * factorial(n - 1) : 1; }

// 二重階乗 (ただし n >= -1 のときしか正しくない) 
constexpr int double_factorial(int n) { return n >= 2 ? n * double_factorial(n - 2) : 1; }

// (-1)^n
constexpr int pow_m1(int n) { return n % 2 == 0 ? 1 : -1; }


// 参照資料 http://t-pot.com/program/88_SH/
// 使ってる文字などほぼそのまま

 // ルジャンドル陪関数 P_l^m
//-----------------------------------------------
// l=0, m=0 → float P<0,0>(float);
// l=1, m=0 → float P<1,0>(float);
// ...

// ルジャンドル陪関数 P_l^m (l == m の場合)
template <unsigned int l, unsigned int m>
constexpr std::enable_if_t<(l == m), float>
P(float x) { return pow_m1(m) * double_factorial(2 * m - 1) * std::pow(1 - x * x, m / 2.0f); }

// ルジャンドル陪関数 P_l^m (l == m+1 の場合)
template <unsigned int l, unsigned int m>
constexpr std::enable_if_t<(l == m + 1), float>
P(float x) { return (2 * m + 1) * x * P<m, m>(x); }

// ルジャンドル陪関数 P_l^m (l > m+1 の場合)
template <unsigned int l, unsigned int m>
constexpr std::enable_if_t<(l > m + 1), float>
P(float x) { return ((2 * l - 1) * x * P<l - 1, m>(x) - (l + m - 1) * P<l - 2, m>(x)) / (l - m); }


 // 球面調和関数 Y_l^m
//--------------------------------------------
// l=0, m= 0 → float Y<0,0>(float);
// l=1, m=-1 → float Y<1,-1>(float);
// l=1, m= 0 → float Y<1,0>(float);
// ...

// 球面調和関数に出てくる係数 K_l^m
template <unsigned int l, int m>
const float K = std::sqrt((2 * l + 1) / (4 * π) * factorial(l - std::abs(m)) / factorial(l + std::abs(m)));

// 球面調和関数 Y_l^m (m > 0 の場合)
template <unsigned int l, int m>
std::enable_if_t<(m > 0), float>
Y(float θ, float φ) { return std::sqrt(2.0f) * K<l, m> * std::cos(m * φ) * P<l, m>(std::cos(θ)); }

// 球面調和関数 Y_l^m (m < 0 の場合)
template <unsigned int l, int m>
std::enable_if_t<(m < 0), float>
Y(float θ, float φ) { return std::sqrt(2.0f) * K<l, m> * std::sin(-m * φ) * P<l, -m>(std::cos(θ)); }

// 球面調和関数 Y_l^m (m == 0 の場合)
template <unsigned int l, int m>
std::enable_if_t<(m == 0), float>
Y(float θ, float φ) { return K<l, 0> * P<l, 0>(std::cos(θ)); }


 // 球面調和関数を取得する(実行時にl,mを指定する版)
//---------------------------------------------------------

// 球面調和関数をstaticに作ってテーブルに置いておく (とりあえず l = 5 まで) 
static const std::vector<std::function<float(float, float)>> Ylm_table = {
	Y<0, 0>,
	Y<1, -1>, Y<1, 0>, Y<1, 1>,
	Y<2, -2>, Y<2, -1>, Y<2, 0>, Y<2, 1>, Y<2, 2>,
	Y<3, -3>, Y<3, -2>, Y<3, -1>, Y<3, 0>, Y<3, 1>, Y<3, 2>, Y<3, 3>,
	Y<4, -4>, Y<4, -3>, Y<4, -2>, Y<4, -1>, Y<4, 0>, Y<4, 1>, Y<4, 2>, Y<4, 3>, Y<4, 4>,
	Y<5, -5>, Y<5, -4>, Y<5, -3>, Y<5, -2>, Y<5, -1>, Y<5, 0>, Y<5, 1>, Y<5, 2>, Y<5, 3>, Y<5, 4>, Y<5, 5>
};

static const unsigned int MAX_L = 5;

// 球面調和関数の取得
// Y(0, 0) は Y<0, 0>を返す
// Y(1, -1) は Y<1, -1>を返す
// ...
const std::function<float(float, float)> &Y(unsigned int l, int m) {
	assert(l <= MAX_L);
	assert(-static_cast<int>(l) <= m && m <= static_cast<int>(l));
	return Ylm_table[l * l + l + m];
}

} // end of anonymous namespace



namespace sh {

 // SHSクラスメンバ関数
//=========================

SHS::SHS(unsigned int max_l) : _max_l(max_l), _c((max_l + 1) * (max_l + 1)) {
	assert(max_l <= MAX_L);
}

SHS::SHS(const std::vector<float> &c) {
	// c.size() == (_max_l + 1) * (_max_l + 1) の関係を満たすべき
	_max_l = static_cast<unsigned int>(std::sqrt(c.size()) - 1);
	assert(c.size() == (_max_l + 1) * (_max_l + 1));
	_c = c;
}

unsigned int SHS::max_l() const {
	return _max_l;
}

float SHS::coeff(unsigned int l, int m) const {
	assert(l <= _max_l);
	assert(-static_cast<int>(l) <= static_cast<int>(m) && m <= static_cast<int>(l));
	return _c[l * l + l + m];
}

unsigned int SHS::num_coeffs() const {
	return _c.size();
}

const float *SHS::coeffs() const {
	return _c.data();
}

void SHS::set_coeff(unsigned int l, int m, float c) {
	assert(l <= _max_l);
	assert(-static_cast<int>(l) <= static_cast<int>(m) && m <= static_cast<int>(l));
	_c[l * l + l + m] = c;
}

void SHS::set_coeffs(const float *v, int n) {
	_c.assign(v, v + n);
	_c.resize((_max_l + 1) * (_max_l + 1));
}

void SHS::set_coeffs(const std::vector<float> &c) {
	_c = c;
	_c.resize((_max_l + 1) * (_max_l + 1));
}



void SHS::spherical_harmonics_transform(const std::function<float(float, float)> &f, unsigned int div_theta, unsigned int div_phi) {
	for (unsigned int l = 0; l <= _max_l; ++l) {
		for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m) {
			const float c = spherical_integration(l, m, div_theta, div_phi, f);
			this->set_coeff(l, m, c);
		}
	}
}

float SHS::spherical_integration(
	unsigned int l, int m,
	unsigned int dt, unsigned int dp,
	const std::function<float(float, float)> &f
) {
	const std::function<float(float, float)> &Ylm = Y(l, m);
	const unsigned int NUM_DIV_θ = dt;
	const unsigned int NUM_DIV_φ = dp;
	const float dθ = π / NUM_DIV_θ;
	const float dφ = 2 * π / NUM_DIV_φ;

	float val = 0; // 積分結果

	// θ方向の積分
	for (unsigned int t = 0; t < NUM_DIV_θ; ++t) {
		const float θ = dθ * (t + 0.5f);

		float val_t = 0; // θにおける積分結果

		// φ方向の積分
		for (unsigned int p = 0; p < NUM_DIV_φ; ++p) {
			const float φ = dφ * (p + 0.5f);
			val_t += f(θ, φ) * Ylm(θ, φ);
		}

		val += val_t * std::sin(θ);
	}

	val *= dθ * dφ;

	return val;
}

float SHS::eval(float θ, float φ) const {
	float val = 0;
	for (unsigned int l = 0; l <= _max_l; ++l) {
		for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m) {
			const float c_l_m = this->coeff(l, m);
			val += c_l_m * Y(l, m)(θ, φ);
		}
	}
	return val;
}

float SHS::eval(float x, float y, float z) const {
	const float r = std::sqrt(x*x + y*y);
	if (r == 0.0f) {
		return eval(z > 0 ? 0 : π, 0);
	} else {
		const float θ = std::atan2(r, z);
		const float φ = std::atan2(y, x);
		return eval(θ, φ);
	}
}




 // SHSオブジェクトの比較
//-----------------------------------------------------------

bool operator==(const SHS &shs1, const SHS &shs2) {
	// 係数の数は同じでないといけない
	const int nc1 = shs1.num_coeffs();
	const int nc2 = shs2.num_coeffs();
	if (nc1 != nc2) return false;

	// 係数が同じならOK
	for (int i = 0; i < nc1; ++i) {
		const float c0 = shs1.coeffs()[i];
		const float c1 = shs2.coeffs()[i];
		if (c0 != c1) return false;
	}
	return true;
}

bool operator!=(const SHS &shs1, const SHS &shs2) {
	return !(shs1 == shs2);
}

bool near_equal(const SHS &shs1, const SHS &shs2, float error_rate) {
	// 係数の数は同じでないといけない
	const int nc1 = shs1.num_coeffs();
	const int nc2 = shs2.num_coeffs();
	if (nc1 != nc2) return false;

	// 係数が大体同じならOK
	for (int i = 0; i < nc1; ++i) {
		const float c0 = shs1.coeffs()[i];
		const float c1 = shs2.coeffs()[i];

		if (c0 == 0) if (std::abs(c1) > error_rate) return false; else continue;
		if (c1 == 0) if (std::abs(c0) > error_rate) return false; else continue;

		const float e = std::abs(c0 - c1) / std::max(c0, c1);
		if (e > error_rate) return false;
	}
	return true;
}




 // 球面調和関数の回転
//--------------------------------------------------------------

// z軸周りにγ、y軸周りにβ、z軸周りにαの順に回転させる
SHS rotate(const SHS &shs, float α, float β, float γ) {
	const int L = shs.max_l();

	// 回転行列を求める
	std::vector<std::vector<double>> matrices(L+1);
	std::vector<double *> mat_of_l(L+1);
	for (int l = 0; l <= L; ++l) {
		matrices[l].resize((2*l+1)*(2*l+1));
		mat_of_l[l] = matrices[l].data();
	}
	tinysh_rotation(mat_of_l.data(), L, α, -β, γ); // (tinysh_rotationはy軸回転の方向が逆なのでマイナスを付けておく)(もしかして左手系？)

	// 新しい球面調和級数を返す
	SHS shs2(L);

	// 引数shsの中の係数に、tinysh_rotationで作った回転行列を掛けて、shs2に格納する
	for (int l = 0; l <= L; ++l) {
		const std::vector<double> &mat_l = matrices[l];
		const int vlen = 2 * l + 1;
		const float *v = shs.coeffs() + l * l;
		for (int m = -l; m <= l; ++m) {
			float c = 0;
			for (int i = 0; i < vlen; ++i) {
				c += static_cast<float>(mat_l[vlen * (l + m) + i]) * v[i];
			}
			shs2.set_coeff(l, m, c);
		}
	}

	return std::move(shs2);
}

// x軸周りにangleだけ回転させる
SHS rotate_x(const SHS &shs, float angle) {
	return std::move(rotate(shs, -π/2, angle, π/2));
}

// y軸周りにangleだけ回転させる
SHS rotate_y(const SHS &shs, float angle) {
	return std::move(rotate(shs, 0, angle, 0));
}

// z軸周りにangleだけ回転させる
SHS rotate_z(const SHS &shs, float angle) {

	SHS shs2(shs.max_l());

	const unsigned int ml = shs.max_l();
	for (unsigned int l = 0; l <= ml; ++l) {

		const float c0 = shs.coeff(l, 0);
		shs2.set_coeff(l, 0, c0);

		for (int m = 1; m <= static_cast<int>(l); ++m) {
			const float c = std::cos(m * angle);
			const float s = std::sin(m * angle);

			const float x = shs.coeff(l, m);
			const float y = shs.coeff(l, -m);

			const float x2 = c * x - s * y;
			const float y2 = s * x + c * y;

			shs2.set_coeff(l, m, x2);
			shs2.set_coeff(l, -m, y2);
		}
	}

	return std::move(shs2);
}




 // テスト
//===============================================

void test() {
	// ベースになる関数類のテスト
	assert(factorial(0) == 1);
	assert(factorial(1) == 1);
	assert(factorial(2) == 2);
	assert(factorial(3) == 6);
	assert(factorial(4) == 24);
	assert(factorial(5) == 120);
	assert(factorial(6) == 720);
	assert(factorial(7) == 5040);

	assert(double_factorial(-1) == 1);
	assert(double_factorial(0) == 1);
	assert(double_factorial(1) == 1);
	assert(double_factorial(2) == 2);
	assert(double_factorial(3) == 3);
	assert(double_factorial(4) == 8);
	assert(double_factorial(5) == 15);
	assert(double_factorial(6) == 48);
	assert(double_factorial(7) == 105);
	assert(double_factorial(8) == 384);

	assert((P<0, 0>(-1) == 1));
	assert((P<0, 0>(0) == 1));
	assert((P<0, 0>(1) == 1));

	assert((P<1, 0>(-1) == -1));
	assert((P<1, 0>(0) == 0));
	assert((P<1, 0>(1) == 1));

	assert((P<1, 1>(-1) == 0));
	assert((P<1, 1>(0) == -1));
	assert((P<1, 1>(1) == 0));

	assert((P<2, 0>(-1) == 1));
	assert((P<2, 0>(0) == -0.5));
	assert((P<2, 0>(1) == 1));

	assert((P<2, 1>(-1) == 0));
	assert((P<2, 1>(0) == 0));
	assert((P<2, 1>(1) == 0));

	assert((P<2, 2>(-1) == 0));
	assert((P<2, 2>(0) == 3));
	assert((P<2, 2>(1) == 0));


	// 球面調和関数のテスト

	// 適当な係数のSHを作る(3次)
	SHS sh(3);
	sh.set_coeffs({
		16,
		1, 2, 3,
		4, 5, 6, 7, 8,
		9, 10, 11, 12, 13, 14, 15
	});

	assert(sh == sh);

#if 1
	// 上ので作ったSH関数を評価したものをSH変換して、元のSH関数と同じになるかテスト
	SHS sh2(sh.max_l());
	sh2.spherical_harmonics_transform([&](float t, float p){return sh.eval(t, p);});

	// 係数が大体同じならOK
	assert(near_equal(sh, sh2, 0.01f));
#endif

#if 1
	// z軸回転のテスト
	// rotate_z関数で係数を回転させたものと、元のSH関数の座標系を回転させて評価したものをSH変換して、
	// 大体同じになることを確かめる
	const float zRotAngle = 2;
	const SHS sh3 = rotate_z(sh, zRotAngle);
	//const SHS sh3 = rotate(sh, zRotAngle, 0, 0);
	SHS sh4(sh.max_l());
	sh4.spherical_harmonics_transform([&](float t, float p){return sh.eval(t, p - zRotAngle);});
	assert(near_equal(sh3, sh4));
#endif

#if 1
	// y軸回転のテスト
	// rotate_y関数で係数を回転させたものと、元のSH関数の座標系を回転させて評価したものをSH変換して、
	// 大体同じになることを確かめる
	const float yRotAngle = 2;
	const SHS sh5 = rotate_y(sh, yRotAngle);
	//const SHS sh5 = rotate(sh, 0, yRotAngle, 0);
	SHS sh6(sh.max_l());
	sh6.spherical_harmonics_transform([&](float t, float p) {
		const float x = std::cos(p) * std::sin(t);
		const float y = std::sin(p) * std::sin(t);
		const float z = std::cos(t);

		const float c = std::cos(-yRotAngle);
		const float s = std::sin(-yRotAngle);

		const float rz = z * c - x * s;
		const float rx = z * s + x * c;
		return sh.eval(rx, y, rz);
	});
	assert(near_equal(sh5, sh6));
#endif

#if 1
	// z軸周りにγ、y軸周りにβ、z軸周りにαの順に回転させるテスト
	// rotate関数で係数を回転させたものと、元のSH関数の座標系を回転させて評価したものをSH変換して、
	// 大体同じになることを確かめる
	const float α = 1;
	const float β = 2;
	const float γ = 3;
	const SHS sh7 = rotate(sh, α, β, γ);
	SHS sh8(sh.max_l());
	sh8.spherical_harmonics_transform([&](float t, float p) {
		const float x0 = std::cos(p) * std::sin(t);
		const float y0 = std::sin(p) * std::sin(t);
		const float z0 = std::cos(t);

		const float c0 = std::cos(-α);
		const float s0 = std::sin(-α);
		const float x1 = x0 * c0 - y0 * s0;
		const float y1 = x0 * s0 + y0 * c0;
		const float z1 = z0;

		const float c1 = std::cos(-β);
		const float s1 = std::sin(-β);
		const float z2 = z1 * c1 - x1 * s1;
		const float x2 = z1 * s1 + x1 * c1;
		const float y2 = y1;

		const float c2 = std::cos(-γ);
		const float s2 = std::sin(-γ);
		const float x3 = x2 * c2 - y2 * s2;
		const float y3 = x2 * s2 + y2 * c2;
		const float z3 = z2;

		return sh.eval(x3, y3, z3);
	});
	assert(near_equal(sh7, sh8));
#endif

#if 1
	// x軸周りにxRotAngle回転させるテスト
	// rotate_x関数で係数を回転させたものと、元のSH関数の座標系を回転させて評価したものをSH変換して、
	// 大体同じになることを確かめる
	const float xRotAngle = 1;
	const SHS sh9 = rotate_x(sh, xRotAngle);
	SHS sh10(sh.max_l());
	sh10.spherical_harmonics_transform([&](float t, float p) {
		const float x = std::cos(p) * std::sin(t);
		const float y = std::sin(p) * std::sin(t);
		const float z = std::cos(t);

		const float c = std::cos(-xRotAngle);
		const float s = std::sin(-xRotAngle);

		const float ry = y * c - z * s;
		const float rz = y * s + z * c;
		return sh.eval(x, ry, rz);
	});
	assert(near_equal(sh9, sh10));
#endif
}


} // end of namespace sh
