# 11.25.1
QR
// homework2.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include<math.h>
#include <iomanip>
#include<stdlib.h>
using namespace std;
const double epsilon =1.0e-12;
const int n = 10;
double lambda[n];//存实数特征值//
double lambdaf[n][2];//存复数特征值//
double A[n][n];
double Ai[n][n];//存储未经拟上三角化的矩阵A//
double At[n][n];//A的转置矩阵//
double Aqr[n][n];//执行QR分解后的矩阵A//
double Aqrt[n][n];//Aqr的转置矩阵//
double Q[n][n], R[n][n];
double Mk[n][n];
double Mkt[n][n];//Mk的转置矩阵//
void Hessenberg();//矩阵的拟上三角化//
void QRdecomposeA();//对经过拟上三角化后的矩阵A进行QR分解//
void QRdecompose(); //带双步位移的QR分解求矩阵的特征值//
void Transposition(int m);//矩阵的转置//
void TranspositionMk(int m);//Mk矩阵的转置//
void Gauss(double lambda);//列主元素Gauss消去法求解特征向量//
int Judge(int r,int n);//判断矩阵A某列的某些值是否全为零//
int main()
{
	//给矩阵A赋初值//
	for (int i = 0; i <= n-1; i++) {
		for (int j = 0; j <= n-1; j++) {
			if (i == j) {
				A[i][j] = 1.52*cos((i + 1) + 1.2*(j + 1));
				Ai[i][j] = 1.52*cos((i + 1) + 1.2*(j + 1));
			}

			else {
				A[i][j] = sin(0.5*(i + 1) + 0.2*(j + 1));
				Ai[i][j] = sin(0.5*(i + 1) + 0.2*(j + 1));
			}
		}
	}
	//对矩阵A进行拟上三角化//
	Hessenberg();
	//输出拟上三角化后的矩阵A//
	cout << setiosflags(ios::scientific) << setprecision(12);
	cout << endl << "拟上三角化后的矩阵为：" << endl;
	for (int i = 0; i <= n-1; i++) {
		for (int j = 0; j <= n-1; j++) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	//对经过拟上三角化后的矩阵A进行QR分解//
	QRdecomposeA();
	cout << endl << "拟上三角矩阵经QR分解后的正交矩阵Q为：" << endl;
	for (int i = 0; i <= n - 1; i++) {
		for (int j = 0; j <= n - 1; j++) {
			cout << Q[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl << "拟上三角矩阵经QR分解后的上三角矩阵R为：" << endl;
	for (int i = 0; i <= n - 1; i++) {
		for (int j = 0; j <= n - 1; j++) {
			cout << R[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl << "矩阵RQ为：" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double c = 0;
			for (int k = 0; k < n; k++) {
				c += R[i][k] * Q[k][j];
			}
			cout << c << " ";
		}
		cout << endl;
	}
	//带双步位移的QR分解求特征值//
	QRdecompose();
	cout << endl << "特征值如下：";
	for (int i = 0; i < n; i++) {
		cout << endl;
		if (lambda[i] != 0) {
			cout << lambda[i] << " ";
		}
		else {
			cout << lambdaf[i][0] << " ,"<< lambdaf[i][1]<<" *i";
		}
	}
	//求实数特征值对应的特征向量//
	for (int i = 0; i < n; i++) {
		cout << endl;
		if (lambda[i] != 0) {
			cout << endl<< "特征值:" << lambda[i] << " 的特征向量为：";
			Gauss(lambda[i]);
		}
	}
    return 0;
}

void Hessenberg() {
	for (int r = 1; r <= n-2; r++) {
		//判断A(i,r)是否全为零
		if (Judge(r,n) == 1) {//当不全为零时，执行以下计算
			double dr=0, cr=0, hr=0;
			double ur[n];
			double pr[n], qr[n], tr=0, wr[n];
			for (int i = r + 1; i <= n; i++) {
				dr += A[i - 1][r - 1] * A[i - 1][r - 1];
			}
			dr = sqrt(dr);
			if (A[r][r - 1] > epsilon) cr = -dr;
			else  cr = dr;
			hr = cr*cr - cr*A[r][r - 1];
			for (int i = 1; i <= n; i++) {
				if (i <= r) ur[i - 1] = 0;
				else if (i == r + 1) ur[i - 1] = A[i - 1][r - 1] - cr;
				else if (i > r + 1) ur[i - 1] = A[i - 1][r - 1];
			}
			Transposition(n);//得到转置后的矩阵At
			for (int i = 1; i <= n; i++) {
				pr[i-1] = 0;
				for (int j = 1; j <= n; j++) {
					pr[i - 1] += At[i - 1][j - 1] * ur[j - 1];
				}
				pr[i - 1] = pr[i - 1] / hr;
			}
			for (int i = 1; i <= n; i++) {
				qr[i - 1] = 0;
				for (int j = 1; j <= n; j++) {
					qr[i - 1] += A[i - 1][j - 1] * ur[j - 1];
				}
				qr[i - 1] = qr[i - 1] / hr;
			}
			for (int i = 1; i <= n; i++) {
				tr += pr[i - 1] * ur[i - 1];
			}
			tr = tr / hr;
			for (int i = 1; i <= n; i++) {
				wr[i - 1] = qr[i - 1] - tr*ur[i - 1];
			}
			for (int i = 1; i <= n; i++) {
				for (int j = 1; j <= n; j++) {
					A[i - 1][j - 1] = A[i - 1][j - 1] - wr[i - 1] * ur[j - 1] - ur[i - 1] * pr[j - 1];
				}
			}
		}
	}
}
int Judge(int r,int n) {
	for (int i = r + 2; i <= n; i++) {
		if (abs(A[i - 1][r - 1]) >= epsilon) return 1;//返回1则说明有数不为零
	}
	return 0;//全为零，则返回零
}
void Transposition(int m) {
	for (int i = 0; i < m; i++) {//将矩阵A目前的值赋矩阵At
		for (int j = 0; j < m; j++) {
			At[i][j] = A[i][j];
		}
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < i; j++) { //j<i保证了只进行一半的元素的转换
			swap(At[i][j], At[j][i]);
		}
	}
}
void TranspositionMk(int m) {
	for (int i = 0; i < m; i++) {//将矩阵Mk目前的值赋矩阵Mkt
		for (int j = 0; j < m; j++) {
			Mkt[i][j] = Mk[i][j];
		}
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < i; j++) { //j<i保证了只进行一半的元素的转换
			swap(Mkt[i][j], Mkt[j][i]);
		}
	}
}
void TranspositionAqr() {
	for (int i = 0; i < n; i++) {//将矩阵Aqr目前的值赋矩阵Aqrt
		for (int j = 0; j < n; j++) {
			Aqrt[i][j] = Aqr[i][j];
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) { //j<i保证了只进行一半的元素的转换
			swap(Aqrt[i][j], Aqrt[j][i]);
		}
	}
}
void QRdecompose() {
	int k = 1;
	const int L = 500;
	for (int m = n; m > 0;) {
		if (abs(A[m - 1][m - 2]) <= epsilon) {
			lambda[m - 1] = A[m - 1][m - 1];
			m = m - 1;
		}
		else {
			double s1, s2;
			double a, b, c;
			a = 1; b = A[m - 2][m - 2] + A[m - 1][m - 1];
			c = A[m - 2][m - 2] * A[m - 1][m - 1] - A[m - 2][m - 1] * A[m - 1][m - 2];
			if (abs(A[m - 2][m - 3]) <= epsilon) {
				if ((b*b - 4 * a*c) > 0) {
					s1 = (b + sqrt(b*b - 4 * a*c)) / (2 * a);
					s2 = (b - sqrt(b*b - 4 * a*c)) / (2 * a);
					lambda[m - 1] = s1;
					lambda[m - 2] = s2;
					m = m - 2;
				}
				else if((b*b - 4 * a*c) <= 0){
					lambda[m - 1] = 0;
					lambda[m - 2] = 0;
					lambdaf[m - 1][0] = b / (2 * a);
					lambdaf[m - 1][1] = sqrt(-(b*b - 4 * a*c)) / (2 * a);
					lambdaf[m - 2][0] = b / (2 * a);
					lambdaf[m - 2][1] = -sqrt(-(b*b - 4 * a*c)) / (2 * a);
					m = m - 2;
				}
			}
			else if (m == 2) {
				if ((b*b - 4 * a*c) > 0) {
					s1 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
					s2 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
					lambda[m - 1] = s1;
					lambda[m - 2] = s2;
				}
				else {
					lambdaf[m - 1][0] = -b / (2 * a);
					lambdaf[m - 1][1] = sqrt(-(b*b - 4 * a*c));
					lambdaf[m - 2][0] = -b / (2 * a);
					lambdaf[m - 2][1] = -sqrt(-(b*b - 4 * a*c));
				}
				break;
			}
			else if (k == L) {
				cout << "未能得到A的全部特征值" << endl;
				break;
			}
			else {
				double s, t;
				s = A[m - 2][m - 2] + A[m - 1][m - 1];
				t = A[m - 2][m - 2] * A[m - 1][m - 1] - A[m - 2][m - 1] * A[m - 1][m - 2];
				double Ak2[n][n], I[n][n];
				for (int i = 1; i <= m; i++) {
					for (int j = 1; j <= m; j++) {
						Ak2[i - 1][j - 1] = 0;
						for (int k = 0; k < m; k++) {
								Ak2[i - 1][j - 1] += A[i - 1][k] * A[k][j - 1];
						}
					}
				}
				for (int i = 1; i <= m; i++) {
					for (int j = 1; j <= m; j++) {
						if (i == j) I[i - 1][j - 1] = 1;
						else I[i - 1][j - 1] = 0;
					}
				}
				for (int i = 1; i <= m; i++) {
					for (int j = 1; j <= m; j++) {
						Mk[i - 1][j - 1] = Ak2[i - 1][j - 1] - s*A[i - 1][j - 1] + t*I[i - 1][j - 1];
					}
				}
				//对矩阵Mk做QR分解//
				for (int r = 1; r <= m - 1; r++) {
					if (Judge(r - 1, n) == 1) {//当不全为零时，执行以下计算
						double dr = 0, cr = 0, hr = 0;
						double ur[n];
						double vr[n],pr[n], qr[n], tr = 0, wr[n];
						for (int i = r; i <= m; i++) {
							dr += Mk[i - 1][r - 1] * Mk[i - 1][r - 1];
						}
						dr = sqrt(dr);
						if (Mk[r - 1][r - 1] > epsilon) cr = -dr;
						else  cr = dr;
						hr = cr*cr - cr*Mk[r - 1][r - 1];
						for (int i = 1; i <= m; i++) {
							if (i <= r - 1) ur[i - 1] = 0;
							else if (i == r) ur[i - 1] = Mk[i - 1][r - 1] - cr;
							else if (i > r) ur[i - 1] = Mk[i - 1][r - 1];
						}
						TranspositionMk(m);//得到转置后的矩阵Mkt
						for (int i = 1; i <= m; i++) {
							vr[i - 1] = 0;
							for (int j = 1; j <= m; j++) {
								vr[i - 1] += Mkt[i - 1][j - 1] * ur[j - 1];
							}
							vr[i - 1] = vr[i - 1] / hr;
						}
						for (int i = 1; i <= m; i++) {
							for (int j = 1; j <= m; j++) {
								Mk[i - 1][j - 1] = Mk[i - 1][j - 1] - ur[i - 1] * vr[j - 1];
							}
						}
						Transposition(m);//得到转置后的矩阵At
						for (int i = 1; i <= m; i++) {
							pr[i - 1] = 0;
							for (int j = 1; j <= m; j++) {
								pr[i - 1] += At[i - 1][j - 1] * ur[j - 1];
							}
							pr[i - 1] = pr[i - 1] / hr;
						}
						for (int i = 1; i <= m; i++) {
							qr[i - 1] = 0;
							for (int j = 1; j <= m; j++) {
								qr[i - 1] += A[i - 1][j - 1] * ur[j - 1];
							}
							qr[i - 1] = qr[i - 1] / hr;
						}
						for (int i = 1; i <= m; i++) {
							tr += pr[i - 1] * ur[i - 1];
						}
						tr = tr / hr;
						for (int i = 1; i <= m; i++) {
							wr[i - 1] = qr[i - 1] - tr*ur[i - 1];
						}
						for (int i = 1; i <= m; i++) {
							for (int j = 1; j <= m; j++) {
								A[i - 1][j - 1] = A[i - 1][j - 1] - wr[i - 1] * ur[j - 1] - ur[i - 1] * pr[j - 1];
							}
						}
					}
				}
				k = k + 1;
			}
		}
		if (m == 1) {
			lambda[m - 1] = A[m - 1][m - 1];
			m = 0;
		}
	}
}
void QRdecomposeA() {
	for (int i = 0; i < n; i++) {//将矩阵A目前的值赋矩阵Aqr
		for (int j = 0; j < n; j++) {
			Aqr[i][j] = A[i][j];
		}
	}
	for (int i = 1; i <= n; i++) {//给矩阵Q赋初值
		for (int j = 1; j <= n; j++) {
			if (i == j) Q[i - 1][j - 1] = 1;
			else Q[i - 1][j - 1] = 0;
		}
	}
	for (int r = 1; r <= (n - 1); r++) {
		int judge=0;
		for (int i = r + 1; i <= n; i++) {
			if (abs(Aqr[i - 1][r - 1]) >= epsilon) judge = 1;
		}
		if (judge == 1) {
			double dr = 0, cr = 0, hr = 0;
			double ur[n];
			double pr[n], qr[n], tr = 0, wr[n];
			for (int i = r; i <= n; i++) {
				dr += Aqr[i - 1][r - 1] * Aqr[i - 1][r - 1];
			}
			dr = sqrt(dr);
			if (Aqr[r - 1][r - 1] > epsilon) cr = -dr;
			else  cr = dr;
			hr = cr*cr - cr*Aqr[r - 1][r - 1];
			for (int i = 1; i <= n; i++) {
				if (i < r) ur[i - 1] = 0;
				else if (i == r) ur[i - 1] = Aqr[i - 1][r - 1] - cr;
				else if (i > r) ur[i - 1] = Aqr[i - 1][r - 1];
			}
			for (int i = 1; i <= n; i++) {
				wr[i - 1] = 0;
				for (int j = 1; j <= n; j++) {
					wr[i - 1] += Q[i - 1][j - 1] * ur[j - 1];
				}
			}
			for (int i = 1; i <= n; i++) {
				for (int j = 1; j <= n; j++) {
					Q[i - 1][j - 1] = Q[i - 1][j - 1] - wr[i - 1] * ur[j - 1] / hr;
				}
			}
			TranspositionAqr();
			for (int i = 1; i <= n; i++) {
				pr[i - 1] = 0;
				for (int j = 1; j <= n; j++) {
					pr[i - 1] += Aqrt[i - 1][j - 1] * ur[j - 1];
				}
				pr[i - 1] = pr[i - 1] / hr;
			}
			for (int i = 1; i <= n; i++) {
				for (int j = 1; j <= n; j++) {
					Aqr[i - 1][j - 1] = Aqr[i - 1][j - 1] - ur[i - 1] * pr[j - 1];
				}
			}
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			R[i][j] = Aqr[i][j];
		}
	}
}
void Gauss(double lambda) {
	double Aig[n][n], I[n][n];
	double x[n], b[n];
	double m[n][n];
	for (int i = 0; i < n; i++) {//单位矩阵I
		for (int j = 0; j < n; j++) {
			if (i == j) I[i][j] = 1;
			else I[i][j] = 0;
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Aig[i][j] = Ai[i][j] - lambda*I[i][j];
		}
	}
	for (int k = 1; k <= (n - 1); k++) {//消元过程
		int ik=0;
		double Max = 0;
		for (int i = k; i <= n; i++) {
			if (abs(Aig[i - 1][k - 1]) > Max) {
				Max = abs(Aig[i - 1][k - 1]);
				ik = i;
			}
		}
		for (int j = k; j <= n; j++) {
			swap(Aig[k - 1][j - 1], Aig[ik - 1][j - 1]);
		}
		for (int i = k + 1; i <= n; i++) {
			m[i - 1][k - 1] = Aig[i - 1][k - 1] / Aig[k - 1][k - 1];
			for (int j = k + 1; j <= n; j++) {
				Aig[i - 1][j - 1] = Aig[i - 1][j - 1] - m[i - 1][k - 1] * Aig[k - 1][j - 1];
			}
		}
	}
	x[n - 1] = 1.0;
	for (int k = n - 1; k >= 1; k--) {
		double xx=0;
		for (int j = k + 1; j <= n; j++) {
			xx += Aig[k - 1][j - 1] * x[j - 1];
		}
		x[k - 1] = -xx / Aig[k - 1][k - 1];
	}
	for (int i = 0; i < n; i++) {
		cout << endl;
		cout << x[i];
	}
}
