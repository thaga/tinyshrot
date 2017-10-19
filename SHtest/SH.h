#pragma once

#include <vector>
#include <functional>

// 参考資料
// http://t-pot.com/program/88_SH/

// θやφなど座標系の取り方も、上記ページの下方の図の通り
// http://t-pot.com/program/88_SH/image107.png


namespace sh {

// 球面調和級数クラス
// Spherical Harmonics Series
class SHS {
public:
	// 次数だけを指定するコンストラクタ
	SHS(unsigned int max_l = 4);
	// コピーコンストラクタ
	SHS(const SHS &s2) : _max_l(s2._max_l), _c(s2._c) {}
	// ムーブコンストラクタ
	SHS(const SHS &&s2) : _max_l(s2._max_l), _c(std::move(s2._c)) {}
	// std::vectorで係数を指定するコンストラクタ(係数は1, 4, 9, 16, ...個である必要がある)
	SHS(const std::vector<float> &c);
	// 代入
	SHS &operator=(const SHS &s2) { _max_l = s2._max_l; _c = s2._c; return *this; }

	// 次数を取得
	unsigned int max_l() const;
	// 係数を取得(l, mを指定)
	float coeff(unsigned int l, int m) const;

	// 係数の数を取得
	unsigned int num_coeffs() const;
	// 係数の配列を取得
	const float *coeffs() const;

	// 係数を設定する(l, mを指定)
	void set_coeff(unsigned int l, int m, float c);
	// 係数を配列で設定する
	void set_coeffs(const float *v, int n);
	// 係数をstd::vectorで設定する
	void set_coeffs(const std::vector<float> &c);

	// 球面調和関数への変換
	// f: (θ,φ)を受け取って値を返す関数
	// div_theta: 積分時θ方向の分割数
	// div_phi: 積分時φ方向の分割数
	void spherical_harmonics_transform(const std::function<float(float, float)> &f, unsigned int div_theta = 128, unsigned int div_phi = 256);

	// 球面調和関数の評価
	float eval(float θ, float φ) const;
	float eval(float x, float y, float z) const;

private:
	float spherical_integration(unsigned int l, int m, unsigned int dt, unsigned int dp, const std::function<float(float, float)> &f);

	unsigned int _max_l;
	std::vector<float> _c;
};

// 等しい、等しくない、だいたい等しい
bool operator==(const SHS &shs1, const SHS &shs2);
bool operator!=(const SHS &shs1, const SHS &shs2);
bool near_equal(const SHS &shs1, const SHS &shs2, float error_rate = 0.01f);

// z軸周りにγ、y軸周りにβ、z軸周りにαの順に回転させる
SHS rotate(const SHS &shs, float α, float β, float γ);

// x軸周りにangleだけ回転させる
SHS rotate_x(const SHS &shs, float angle);

// y軸周りにangleだけ回転させる
SHS rotate_y(const SHS &shs, float angle);

// z軸周りにangleだけ回転させる
SHS rotate_z(const SHS &shs, float angle);


// テスト
void test();

} // end of namespace sh
