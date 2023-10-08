---
layout:     post
title:      "Black-Scholes-Merton期权定价公式推导"
subtitle:   "风险中性测度方法"
date:       2021-06-11
author:     "Alex"
catalog: true
tags:
  - 金融数学
  - 期权定价
  - C++
  - Python
---

十几年前上大学期间，因为机缘巧合结识了期权，对期权的定价无疑充满了金融学和数学的集体智慧，初学时囿于知识所限，只知定价公式，无法从头开始进行推导，后来学习了《Stochastic Calculus for Finance》，对其大致过程有了初步了解，书中对期权定价公式导出的基本思路是基于无风险套利原则在货币市场和股票市场进行复制对冲未来行权风险，涉及具体方法主要分两种方法（i）导出BSM偏微分方程，然后直接求解析解（p153-159），比较可惜的是，导出偏微分方程后，未给出方程求解过程而是一笔带过直接给出了解析解，让人意犹未尽，（ii）借助风险中性测度进行推导（p210-220），相比于求解偏微分方程所涉及的技巧，该种方法相对简单直观。

本文主要针对方法（ii）进行总结介绍，其核心思想是借助*Girsanov Theorem*进行测度变换，求得风险中性测度下投资组合折现过程为鞅，然后求条件期望得到定价公式。后面C++代码抄自[C++ for Financial Mathematics](https://nms.kcl.ac.uk/john.armstrong/cppbook/cpp-website.html)，`boost`配置可参考网文。


## 1. 风险中性测度下的股价过程

### 1.1 概率测度$\mathbb{P}$下的股价过程

假设 $W_t, 0 \leq t \leq T $是概率空间$(\Omega, \mathcal{F}, \mathbb{P})$上的布朗运动，$\mathcal{F}_t, 0\leq t \leq T$是该布朗运动的域流，即$W_t$是$\mathcal{F}_t$可测函数，说人话就是在时刻$t$我们获得了信息集$\mathcal{F}_t$，该信息集可以决定$W_t$的值。

考虑股票市场一只股票的股价$S_t$，假设其满足几何布朗运动：
$$
d S_t = \alpha_t S_t d t+\sigma_t S_t d W_t, \quad 0 \leq t \leq T.  \tag{1.1}
$$

其中平均回报率$\alpha_t$和股价波动率$\sigma_t$为适应性过程即$\mathcal{F}_t$可测，同时$S_t$等价于：

$$
{S_t} = S_0 \exp \left[ \int^{t}_{0} \sigma_s dW_s + \int^{t}_{0} (\alpha_s - \frac{\sigma_s^2}{2}) ds \right] \tag{1.2}
$$

公式$(1.1)$与公式$(1.2)$等价性推导：

从公式$(1.1)$得到公式$(1.2)$，我们令
$$
f(t, x) = \log x \tag{1.3}\\
\frac{\partial f(t, x)}{\partial t} = 0 \\
\frac{\partial f(t, x)}{\partial x} = \frac{1}{x} \\
\frac{\partial ^2 f(t, x)}{\partial x^2} = -\frac{1}{x^2}
$$
根据伊藤公式：
$$
\begin{aligned}
d f(t, S_t) &= \frac{\partial f(t, S_t)}{\partial t} dt + \frac{\partial f(t, S_t)}{\partial S_t} dS_t + \frac{1}{2} \frac{\partial ^2 f(t, x)}{\partial {S_t}^2} dS_t dS_t \\
&= 0 dt + \frac{1}{S_t} (\alpha_t S_t d t+\sigma_t S_t d W_t) - \frac{1}{2 S_t^2} \sigma_t^2 S_t^2 dt \\
&= \alpha_t dt + \sigma dW_t - \frac{\sigma_t^2}{2} dt \\
&= \sigma_t dW_t + (\alpha_t - \frac{\sigma_t^2}{2}) dt
\end{aligned} \tag{1.4}
$$

>其中$dW_t dt = 0, dtdt=0,dW_t dW_t = dt$

对公式$\text{1.4}$两边求积分
$$
\int^{t}_{0} d \log f(s, S_s) = \int^{t}_{0} \sigma_s dW_s + \int^{t}_{0} (\alpha_s - \frac{\sigma_s^2}{2}) dt \\
\log S_t - \log S_0 = \int^{t}_{0} \sigma_s dW_s + \int^{t}_{0} (\alpha_s - \frac{\sigma_s^2}{2}) dt \\
{S_t} = S_0 \exp \left[ \int^{t}_{0} \sigma_s dW_s + \int^{t}_{0} (\alpha_s - \frac{\sigma_s^2}{2}) ds \right] \tag{1.5}
$$
以上我们推导了从公式$(1.1)$得到公式$(1.2)$，反过来公式$(1.2)$得到公式$(1.1)$，可以令$f(t,x) = e^x$，利用伊藤公式对其求微分得到。$\square$

假设我们有适应性利率过程$R_t$，定义其折现过程
$$
D_t = e^{-\int_{0}^{t} R_s d s} \tag{1.6}
$$

则其微分为

$$
d D_t = -R_t D_t d t \tag{1.7}
$$

$D_t S_t$及其微分$d (D_t S_t)$分别为：

$$
D_t S_t = S_0 \exp \left[ \int_{0}^{t} \sigma_s d W_s+\int_{0}^{t} ( \alpha_s - R_s - \frac{1}{2} \sigma^{2}_s ) d s\right] \tag{1.8}
$$

$$
\begin{aligned}
d(D_t S_t) &=(\alpha_t - R_t) D_t S_t d t+\sigma_t D_t S_t d W_t \\
&= \sigma_t D_t S_t(\Theta_t d t+d W_t)
\end{aligned} \tag{1.9}
$$

其中定义风险的市场价格$\Theta_t = \frac{\alpha_t - R_t}{\sigma_t} $。

### 1.2 *Girsanov Theorem*

假设 $W_t, 0 \leq t \leq T $是概率空间$(\Omega, \mathcal{F}, \mathbb{P})$上的布朗运动，$\mathcal{F}_t, 0\leq t \leq T$是该布朗运动的域流，$\Theta_t, 0 < t < T$为$\mathcal{F}_t$适应性过程。定义

$$
Z_t =\exp \left( -\int_{0}^{t} \Theta_s d W_s - \frac{1}{2} \int_{0}^{t} \Theta^{2}_s d s\right) \tag{1.10}
$$

$$
\widetilde{W}_t =W_t + \int_{0}^{t} \Theta_s d s \tag{1.11}
$$

并假设

$$
\mathbb{E} \int_{0}^{T} \Theta^{2}_s Z^{2}_s d s < \infty \tag{1.12}
$$

让$Z = Z_T$。那么$\mathbb{E}Z = 1$且在由

$$
\widetilde{\mathbb{P}}(A) = \int_{A}Z(\omega) d \mathbb{P}(\omega) \quad A \in \mathcal{F}_t. \tag{1.13}
$$

定义的概率测度$\widetilde{\mathbb{P}}$下，$\widetilde{W}_t, 0 < t < T$是布朗运动。$\Box$

### 1.3 概率测度$\widetilde{\mathbb{P}}$下的股价过程

根据*Girsanov Theorem*，在概率测度$\widetilde{\mathbb{P}}$下，$d \widetilde W_t = \Theta_t d t + d W_t$，因此公式$\text(1.9)$也可以写为

$$
d(D_t S_t) = \sigma_t D_t S_t d \widetilde W_t \tag{1.14}
$$

两边同时积分

$$
D_t S_t = S_0 + \int_{0}^{t} \sigma_s D_s S_s d \widetilde{W}_s \tag{1.15}
$$

由于在概率测度$\widetilde{\mathbb{P}}$下，$\int_{0}^{t} \sigma_s D_s S_s d \widetilde{W}_s$是伊藤过程，因此是一个鞅。

**因此我们称*Girsanov Theorem*下的概率测度$\widetilde{\mathbb{P}}$为风险中性测度(risk-neutral measure)**。

将$d \widetilde W_t = \Theta_t d t+d W_t$带入公式$(1.1)$，可以得到在概率测度$\widetilde{\mathbb{P}}$下，公式$(1.1)$和$(1.2)$分别可以改写成公式$(1.16)$和公式$(1.17)$的形式

$$
d S_t = R_t S_t d t+\sigma_t S_t d \widetilde {W}_t \tag{1.16}
$$

$$
S_t = S_0 \exp \left[ \int_{0}^{t} \sigma_s d \widetilde{W}_s + \int_{0}^{t} ( R_s - \frac{1}{2} \sigma^{2}_s ) d s\right] \tag{1.17}
$$

## 2. 风险中性测度下投资组合的价值过程

假设投资者初始资本为$X_0$，在任意时间$t, 0 < t < T $持有$\Delta_t$份股票，同时以利率$R_t$投资或借贷于货币市场，以维持自融资状态，则投资组合价值的微分为

$$
\begin{aligned}
d X_t &= \Delta_t d S_t + R_t ( X_t - \Delta_t S_t ) d t \\
\\
&= \Delta_t (\alpha_t S_t d t+\sigma_t S_t d W_t) + R_t (X_t - \Delta_t S_t ) d t \\
\\
&= R_t X_t d t+\Delta_t (\alpha_t - R_t) S_t d t + \Delta_t \sigma_t S_t d W_t \\
\\
&= R_t X_t d t + \Delta_t \sigma_t S_t (\Theta_t d t+d W_t)
\end{aligned} \tag{2.1}
$$

根据$\text{Ito}$乘法法则
$$
d(X_t Y_t) = X_t dY_t + Y_t dX_t + dX_t dY_t \tag{2.2}
$$


由公式$(1.7),(1.9)$和$(2.1)$可得
$$
\begin{aligned}
d(D_t X_t) &= \Delta_t \sigma_t D_t S_t (\Theta_t d t+d W_t) \\
\\
&=\Delta_t d(D_t S_t) \\
\\
&= \Delta_t \sigma_t D_t S_t d \widetilde W_t
\end{aligned} \tag{2.3}
$$

由此投资者有两种选择：

1. 以利率$R_t$投资于货币市场；

2. 在风险中性测度$\widetilde{\mathbb{P}}$下投资于平均回报率为$R_t$的股票。

但由于在风险中性测度$\widetilde{\mathbb{P}}$下，投资组合的折现价值$D_t X_t$是鞅，因此不管投资者如何选择，其投资组合的平均回报率均为为$R_t$。

## 3. 风险中性测度下的期权定价

我们令$\mathcal{F}_T$可测的随机变量$V_T$表示在时刻$T$衍生证券空头的潜在支付(*payoff*) $(S_T-K)^+$，投资者为了对冲看涨期权空头即未来所面临的潜在支出$V_T$，那么其持有的投资组合$X_t$需要使以下等式几乎必然成立(*almost surely*)

$$
X_T = V_T \tag{3.1}
$$

由$D_t X_t$在测度$\widetilde{\mathbb{P}}$是鞅的事实我们有

$$
D_t X_t = \widetilde{\mathbb{E}}(D_T X_T \vert \mathcal{F}_t) = \widetilde{\mathbb{E}} (D_T V_T \vert  \mathcal{F}_t) \tag{3.2}
$$

$X_t$表示在时刻$t$为完全对冲衍生证券支付$V_T$所持有的投资组合价值，我们将其称之为衍生证券在时刻 $t$ 的价格并用$V_t$表示，那么公式$(3.2)$可以写成

$$
D_t V_t = \widetilde{\mathbb{E}}(D_T V_T \vert \mathcal{F}_t), \quad 0 \leq t \leq T \tag{3.3}
$$

由于$D_t$是$\mathcal{F}_t$可测的，因此我们可以将其移到公式右侧，得到

$$
V_t = \widetilde{\mathbb{E}} \left(e^{-\int_{t}^{T} R_s d s} V_T \vert \mathcal{F}_t \right), \quad 0 \leq t \leq T \tag{3.4}
$$

我们将公式$(3.3)$和$(3.4)$称为连续时间下风险中性定价公式(*risk-neutral pricing formula*)。

## 4. 推导$\text{Black-Scholes-Merton}$公式

为简单起见，我们假设$\sigma _t$和$R_t$分别为常数 $\sigma$ 和 $r$，则公式$(3.4)$简化为

$$
\widetilde{\mathbb{E}}\left[e^{-r(T-t)} (S_T - K)^{+} \vert \mathcal{F}_t \right] \tag{4.1}
$$

公式$(4.1)$仅依赖于时刻$t$和股价$S_t$，由于几何布朗运动是马尔可夫过程，因此存在$c(t, S_t)$满足

$$
c(t, S_t) = \widetilde{\mathbb{E}}\left[e^{-r(T-t)}(S_T-K)^{+} \mid \mathcal{F}_t\right] \tag{4.2}
$$

公式$(1.17)$简化为

$$
S_t = S_0 \exp \left[ \sigma \widetilde{W}_t+\left(r-\frac{1}{2} \sigma^{2}\right) t\right] \tag{4.3}
$$

则$S_T$等于

$$
\begin{aligned}
S_T &= S_t \exp \left[ \sigma(\widetilde{W}_T - \widetilde{W}_t) + \left(r-\frac{1}{2} \sigma^{2}\right) \tau \right] \\
\\
&= S_t \exp \left[ -\sigma \sqrt{\tau} Y+\left(r-\frac{1}{2} \sigma^{2}\right) \tau \right]
\end{aligned} \tag{4.4}
$$

其中$\tau = T - t$，$Y$是标准正态随机变量

$$
Y = - \frac{\widetilde{W}_T - \widetilde{W}_t} {\sqrt{T-t}} \tag{4.5}
$$

公式$(4.2)$可以写成如下形式

$$
\begin{aligned}
c(t, x) &= \widetilde{\mathbb{E}}\left [ e^{-r \tau}\left(x \exp \left [ -\sigma \sqrt{\tau} Y+\left(r-\frac{1}{2} \sigma^{2}\right) \tau \right] -K \right)^{+} \right ] \\
\\
&= \frac{e^{-r \tau}}{\sqrt{2 \pi}} \int_{-\infty}^{\infty} \left\{ x \exp \left[ -\sigma \sqrt{\tau} y+\left(r-\frac{\sigma^{2}}{2} \right) \tau \right] - K \right\}^{+} e^{-\frac{1}{2} y^{2}} d y
\end{aligned} \tag{4.6}
$$

其中被积函数

$$
\left\{ x \exp \left[ -\sigma \sqrt{\tau} y+\left(r-\frac{\sigma^{2}}{2} \right) \tau \right] - K \right\}^{+}\tag{4.7}
$$

为正，当且仅当

$$
y < d_{-}(\tau, x)=\frac{1}{\sigma \sqrt{\tau}}\left[\log \frac{x}{K}+\left(r-\frac{\sigma^{2}}{2} \right) \tau\right] \tag{4.8}
$$

由此

$$
\begin{aligned}
c(t, x) &= \frac{e^{-r \tau}}{\sqrt{2 \pi}} \int_{-\infty}^{d_{-}(\tau, x)} \left(x \exp \left\{-\sigma \sqrt{\tau} y+\left(r-\frac{\sigma^{2}}{2} \right) \tau\right\}-K\right) e^{-\frac{1}{2} y^{2}} d y \\
\\
&=\frac{1}{\sqrt{2 \pi}} \int_{-\infty}^{d_{-}(\tau, x)} x \exp \left\{-\frac{y^{2}}{2}-\sigma \sqrt{\tau} y-\frac{\sigma^{2} \tau}{2}\right\} d y \\
\\
& \quad \ - \frac{1}{\sqrt{2 \pi}} \int_{-\infty}^{d_{-}(\tau, x)} e^{-r \tau} K e^{-\frac{1}{2} y^{2}} d y \\
\\
&=\frac{x}{\sqrt{2 \pi}} \int_{-\infty}^{d_{-}(\tau, x)} \exp \left\{-\frac{1}{2}(y+\sigma \sqrt{\tau})^{2}\right\} d y-e^{-r \tau} K N\left(d_{-}(\tau, x)\right) \\
\\
&=\frac{x}{\sqrt{2 \pi}} \int_{-\infty}^{d_{-}(\tau, x)+\sigma \sqrt{\tau}} \exp \left\{-\frac{z^{2}}{2}\right\} d z-e^{-r \tau} K N\left(d_{-}(\tau, x)\right) \\
\\
&=x N\left(d_{+}(\tau, x)\right)-e^{-r \tau} K N\left(d_{-}(\tau, x)\right)
\end{aligned} \tag{4.9}
$$

其中

$$
N(x) = \frac{1}{\sqrt{2 \pi}} \int_{-\infty}^{x} e^ {-\frac{t^2}{2}} d t \tag{5.4}
$$

$$
d_+(\tau, x) = \frac{1}{\sigma \sqrt{\tau}}\left[\log \frac{x}{K}+\left(r+\frac{1}{2} \sigma^{2}\right) \tau\right] \tag{4.10}
$$

$$
d_-(\tau, x) = d_+(\tau, x) - \sigma \sqrt{\tau} \tag{4.11}
$$

由此我们得到了欧式看涨期权定价公式

$$
\text{BSM}(\tau, x, K, r, \sigma) = x N\left(d_{+}(\tau, x)\right) - K e^{-r \tau} N\left(d_{-}(\tau, x)\right) \tag{4.12}
$$


## 5. 期权定价的C++实现 

### 5.1 实现$N(x)$函数 - *From Scratch to Boost Library*

定价公式$(1)$-$(5)$中涉及$N(x)$、$\exp(x)$和$\log(x)$等3个函数，其中$\exp(x)$和$\log(x)$已在标准库`<cmath>`中实现，可以直接使用。因此只剩$N(x)$需要我们在标准库外自己实现或寻求其他库的支持。我们按照如下三种方式分别进行实现：

- 多项式逼近法
- 数值积分法
- 调用`Boost`库

#### 5.1.1 多项式逼近法

如$x > 0$，定义$k = 1/(1 + 0.2316419x)$，则$N(x)$可用如下关于$k$多项式进行逼近
$$
\small 1 − \frac{1}{\sqrt{2 \pi}} \exp(-\frac{x^2} {2}) k(0.319381530 \\
+ k(−0.356563782 + k(1.781477937 + k(−1.821255978 + 1.330274429k)))) \tag{5.1}
$$

借助$N(x) + N(-x) = 1$可以求得$N(-x)$即$x < 0$时$N(x)$的值。

C++实现如下：

```cpp
const double Pi = 3.141592653589793;
const double root2Pi = sqrt(2.0 * Pi);

double normcdf_poly(double x)
{
    if (x < 0)
        return 1 - normcdf_ploy(-x);
    double k = 1 / (1 + 0.2316419 * x);
    double poly = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 
                  + k * (-1.821255978 + 1.330274429 * k))));
    double approx = 1.0 - 1.0 / root2Pi * exp(-0.5 * x * x) * poly;
    return approx;
}
```

#### 5.1.2 数值积分法

对于实数域上$\mathbb R \to \mathbb R$的黎曼积分

$$
F(b) - F(a) = \int_{a}^{b} f(x) d x \tag{5.2}
$$

令$\Delta x = \frac{(b-a)}{N}$，根据黎曼积分定义，我们可以用以下矩形的面积和逼近$\int_{a}^{b} f(x) d x$：

$$
\lim_{N \to \infty} \sum_{i=0}^{N-1} f\left(a+(i+\frac{1}{2}) \Delta x\right) \Delta x \tag{5.3}
$$

只要$N$取值足够大即可很好的逼近积分值，可以通过定义$f(x) = x ^ 2 + 1$进行简单验证，C++实现如下：

```cpp
#include <iostream>
#include <cmath>
using namespace std;

double f_x(double x)
{
    return pow(x, 2) + 1;
}

double integrate_fx(double a, double b, int N)
{
    double delta = (b - a) / N;
    double integrate = 0.0;
    int i = 0;
    while (i < N)
    {
        double x = a + (i + 0.5) * delta;
        integrate += f_x(x) * delta;
        i++;
    }
    return integrate;
}

int main()
{
    double a = 1;
    double b = 3;
    int N = 1000;
    double integrate_ = 0.0;
    
    integrate_ = integrate_fx(a, b, N);
    cout << integrate_ << endl; //返回10+2/3(约等于10.666667)
    
    return 0;
}
```

下面我们考虑如何用数值积分法逼近

$$
N(x)=\frac{1}{\sqrt{2 \pi}} \int_{-\infty}^{x} \exp \left(-\frac{t^{2}}{2}\right) d t
$$

由于积分下限为$- \infty$，我们需要考虑进行变量替换，定义$t = x + 1 - \frac{1} {s}$，则$N(x)$等价于如下积分：

$$
\frac{1} {\sqrt {2 \pi}} \int_{0}^{1} \frac{1}{s^{2}} \exp \left(-\frac{\left(x+1-\frac{1}{s}\right)^{2}}{2}\right) \mathrm{d} s \tag{5.4}
$$

C++实现如下：

```cpp
double normcdf_integrate(double x)
{
    int a = 0;
    int b = 1;
    double N = 1000;
    double delta = 1 / N;
    double integrate = 0.0;
    for (int i = 0; i < N; i++)
    {
        double s = a + (i + 0.5) * delta;
        double t = x + 1 - pow(s, -1);
        double f = pow(2 * PI, -0.5) * pow(s, -2) * exp(-0.5 * t * t);
        integrate += f * delta;
    }
    return integrate;
}
```

#### 5.1.3 调用`boost`库

我们可以通过如下方式调用`boost`库中的累积分布函数（**C**umulative **D**istribution **F**unction，简称$\mathrm{cdf}$)的实现

```cpp
#include "boost\math\distributions\normal.hpp"
#include <random>

double normcdf_boost(double x)
{
	double cdf = 0.0;
	
	boost::math::normal_distribution<> norm(0, 1); //生成标准正态分布
	cdf = boost::math::cdf(norm, x); //计算N(x)值
	return cdf;
}
```

### 5.2 期权费计算的C++实现

实现$N(x)$后，计算公式$(1)$中欧式看涨期权$c(S_0, K, T, r, \sigma)$的障碍就清除了。下面我们分别定义`blackScholesCallPrice`和`blackScholesPutPrice`函数完成计算$c(S_0, K, T, r, \sigma)$和$p(S_0, K, T, r, \sigma)$。

```cpp
double blackScholesCallPrice(double K, double T, double S_t, double vol, double r)
{
    double d1 = pow((vol * sqrt(T)), -1) * (log(S_t / K) + (r + 0.5 * pow(vol, 2)) * sqrt(T));
    double d2 = d1 - vol * sqrt(T);

    //以下normcdf(x)函数3选1，默认选normcdf_poly(x)，使用normcdf_boost(x)时请自行配置boost库
    double callPrice = normcdf_poly(d1) * S_t - normcdf_poly(d2) * K * exp(-r * T); 
    //double callPrice = normcdf_integrate(d1) * S_t - normcdf_integrate(d2) * K * exp(-r * T);
    //double callPrice = normcdf_boost(d1) * S_t - normcdf_boost(d2) * K * exp(-r * T);
    return callPrice;
}
```

```cpp
double blackScholesPutPrice(double K, double T, double S_t, double vol, double r)
{
    double d1 = pow((vol * sqrt(T)), -1) * (log(S_t / K) + (r + 0.5 * vol * vol) * sqrt(T));
    double d2 = d1 - vol * sqrt(T);

    //以下normcdf(x)函数3选1，默认选normcdf_poly(x)，使用normcdf_boost(x)时请自行配置boost库
    double putPrice = normcdf_poly(-d2) * K * exp(-r * T) - normcdf_poly(-d1) * S_t;
    //double putPrice = normcdf_integrate(-d2) * K * exp(-r * T) - normcdf_integrate(-d1) * S_t;
    //double putPrice = normcdf_boost(-d2) * K * exp(-r * T) - normcdf_boost(-d1) * S_t;
    return putPrice;
}
```

### 5.3 完整工程文件及测试

头文件*stdafx.h*

```cpp
#pragma once
#include <iostream>
#include <cmath>
//请在Visual Studio中自行配置好boost库后取消注释
//#include "boost\math\distributions\normal.hpp"
//#include <random>
```

头文件*OptionPricing.h*

```cpp
#pragma once

static const double PI = 3.14159265358979; //定义Pi值

//S_t: spotPrice
//T: maturity
//K: strikePrice
//vol: volatility
//r: riskFreeInterestRate

//声明看涨期权定价公式
double blackScholesCallPrice(double K, double T, double S_t, double vol, double r);

//声明看跌期权定价公式
double blackScholesPutPrice(double K, double T, double S_t, double vol, double r);
```

源文件*OptionPricing.cpp*

```cpp

#include "OptionPricing.h"
#include "stdafx.h"

static const double root2Pi = sqrt(2.0 * PI);

static inline double normcdf_poly(double x)
{
    if (x < 0)
        return 1 - normcdf_poly(-x);
    double k = 1 / (1 + 0.2316419 * x);
    double poly = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937
        + k * (-1.821255978 + 1.330274429 * k))));
    double approx = 1.0 - 1.0 / root2Pi * exp(-0.5 * x * x) * poly;
    return approx;
}

static inline double normcdf_integrate(double x)
{
    int a = 0;
    int b = 1;
    double N = 1000;
    double delta = 1 / N;
    double integrate = 0.0;
    for (int i = 0; i < N; i++)
    {
        double s = a + (i + 0.5) * delta;
        double t = x + 1 - pow(s, -1);
        double f = pow(2 * PI, -0.5) * pow(s, -2) * exp(-0.5 * t * t);
        integrate += f * delta;
    }
    return integrate;
}

//请在Visual Studio中自行配置好boost库后取消注释
/*
static inline double normcdf_boost(double x)
{
    double cdf = 0.0;

    boost::math::normal_distribution<> norm(0, 1); //生成标准正态分布
    cdf = boost::math::cdf(norm, x); //计算N(x)值
    return cdf;
}
*/

//S_t: spotPrice
//T: maturity
//K: strikePrice
//vol: volatility
//r: riskFreeInterestRate

double blackScholesCallPrice(double K, double T, double S_t, double vol, double r)
{
    double d1 = pow((vol * sqrt(T)), -1) * (log(S_t / K) + (r + 0.5 * pow(vol, 2)) * sqrt(T));
    double d2 = d1 - vol * sqrt(T);

    //以下normcdf(x)函数3选1，默认选normcdf_poly(x)，使用normcdf_boost(x)时请自行配置boost库
    double callPrice = normcdf_poly(d1) * S_t - normcdf_poly(d2) * K * exp(-r * T); 
    //double callPrice = normcdf_integrate(d1) * S_t - normcdf_integrate(d2) * K * exp(-r * T);
    //double callPrice = normcdf_boost(d1) * S_t - normcdf_boost(d2) * K * exp(-r * T);
    return callPrice;
}

double blackScholesPutPrice(double K, double T, double S_t, double vol, double r)
{
    double d1 = pow((vol * sqrt(T)), -1) * (log(S_t / K) + (r + 0.5 * vol * vol) * sqrt(T));
    double d2 = d1 - vol * sqrt(T);

    //以下normcdf(x)函数3选1，默认选normcdf_poly(x)，使用normcdf_boost(x)时请自行配置boost库
    double putPrice = normcdf_poly(-d2) * K * exp(-r * T) - normcdf_poly(-d1) * S_t;
    //double putPrice = normcdf_integrate(-d2) * K * exp(-r * T) - normcdf_integrate(-d1) * S_t;
    //double putPrice = normcdf_boost(-d2) * K * exp(-r * T) - normcdf_boost(-d1) * S_t;
    return putPrice;
}
```

测试文件*main.cpp*

```cpp
#include <iostream>
#include "OptionPricing.h"
using namespace std;

int main()
{
    double K = 100.0;
    double S_t = 110.0;
    double vol = 0.1;
    double r = 0.03;
    double T = 0.5;
    double call;
    double put;
    double putCallParity;

    call = blackScholesCallPrice(K, T, S_t, vol, r);
    cout << "The call option price is: " << call << endl; //11.6725
    put = blackScholesPutPrice(K, T, S_t, vol, r);
    cout << "The put option price is: " << put << endl; //0.183688
    cout << endl;

    cout << "<Put-Call-Parity = Call Price - Put Price> Testing" << endl;
    cout << "<Call Price - Put Price> is: " << call - put << endl;
    putCallParity = S_t - exp(-r * T) * K;
    cout << "<Put-Call-Parity> is: " << putCallParity << endl;
    cout << "Wonderful! It's correct." << endl;

    return 0;
}
```

其中看涨-看跌期权平价公式为

$$
S_0 - \exp(-rT)K = c(S_0, K, T, r, \sigma) - p(S_0, K, T, r, \sigma) \tag{5.5}
$$

## 6. 期权定价的Python实现

```python
import math
import numpy as np
from scipy.stats import norm

//定义以下简写变量
//K: strikePrice
//T: maturity
//S_t: spotPrice
//vol: volatility
//r: riskFreeRate

//定义d1和d2
def d1f(K, T, S_t, vol, r) :
    return math.pow(vol * np.sqrt(T), -1) * \
    (np.log(S_t / K) + (r + 0.5 * math.pow(vol, 2) * np.sqrt(T)))

def d2f(K, T, S_t, vol, r) :
    return d1f(K, T, S_t, vol, r) - vol * np.sqrt(T)

//定义看涨-看跌期权计算公式
def blackScholesCallPrice(K, T, S_t, vol, r) :
    d1 = d1f(K, T, S_t, vol, r)
    d2 = d2f(K, T, S_t, vol, r)

    callPrice = norm.cdf(d1) * S_t - norm.cdf(d2) * K * np.exp(-r * T)
    return callPrice

def blackScholesPutPrice(K, T, S_t, vol, r) :
    d1 = d1f(K, T, S_t, vol, r)
    d2 = d2f(K, T, S_t, vol, r)

    putPrice = norm.cdf(-d2) * K * np.exp(-r * T) - norm.cdf(-d1) * S_t
    return putPrice

//定义看涨-看跌期权平价测试公式
def callPutParity(K, T, S_t, r) :
    return S_t - np.exp(-r * T) * K

//对给定变量进行测试
K = 100.0
S_t = 110.0
vol = 0.1
r = 0.03
T = 0.5

call = blackScholesCallPrice(K, T, S_t, vol, r)
put = blackScholesPutPrice(K, T, S_t, vol, r)
callPutParity_ = callPutParity(K, T, S_t, r)

print("The call option price is: {0}".format(call))
print("The put option price is: {0}".format(put))
print("Call price - put price is: {0}".format(call - put))
print("The Call-Put-Parity is: {0}".format(callPutParity_))
```

完结











