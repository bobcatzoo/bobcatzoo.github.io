---
layout:     post
title:      "Black-Scholes-Merton期权定价公式推导"
subtitle:   "BSM偏微分方程求解"
date:       2023-09-29
author:     "Alex"
catalog: true
tags:
  - 金融数学
  - 期权定价
---

## Black-Scholes-Merton期权定价公式推导

10几年前上大学期间，因为机缘巧合结识了期权，对期权的定价无疑充满了金融学和数学的集体智慧，初学时囿于知识所限，只知定价公式，无法从头开始进行推导，后来学习了《Stochastic Calculus for Finance》，对其大致过程有了初步了解，书中对期权定价公式导出的基本思路是基于无风险套利原则在货币市场和股票市场进行复制对冲未来行权风险，涉及具体方法主要分两种方法（i）导出BSM偏微分方程，然后直接求解析解（p153-159），比较可惜的是，导出偏微分方程后，未给出方程求解过程而是一笔带过直接给出了解析解，让人意犹未尽，（ii）借助风险中性测度进行推导（p210-220），相比于求解偏微分方程所涉及的技巧，该种方法相对简单直观。

本文主要针对方法（i），先后参考了《期权定价的数学模型和方法》（第二版）第五章，给出了BSM偏微分方程的求解过程，但美中不足是对热传导方程基本解一笔带过，少了傅里叶变换求解过程（p81）。于是又参考了《数学物理方程》（第二版）季孝达等编第四章（p137-138）对热传导方程的求解方法，虽然对其中的傅里叶变换还理解不深，但完整的思路算是齐了。以下就将散落在这三本书里的知识汇总在一起，对方法（i）做一个相对完整的注记。

### 一、投资组合价值演化

假设投资人A在初始时刻仅靠出售欧式看涨期权，获得期权费然后在货币市场和股票市场进行组合投资以对冲未来期权行权风险，即确保该组合价值在行权日等于$(S_T - K)^+$.

假设投资人A在时刻$t$的投资组合价值为$X_t$，该组合投资于固定利率为$r$的货币市场（比如银行存款等）和股票市场（简化为一只股票$S_t$），该股票由几何布朗运动驱动即：
$$
d S_t = \alpha S_t dt + \sigma S_t d W_t \tag{1.1}
$$
或
$$
\frac{d S_t}{S_t} = \alpha dt + \sigma d W_t \tag{1.2}
$$

> 公式$\text{(1.1)}$由正收益项$\alpha S_t dt$和随机波动项$\sigma S_t dW_t$构成。可以理解为瞬时股票收益$\frac{d S_t}{S_t}$隐含一个正的收益项$\alpha dt$，同时包含一个波动项$\sigma d W_t$，两者叠加后的效果就是股票价格围绕一个向上的趋势上下波动。而且由于股票市场风险要高于货币市场，因此满足$\alpha > r$。

假设在时刻$t$，投资人A持有$\Delta_t$份额的股票，$\Delta_t$可以是随机的(random)，但由布朗运动$W_t, t \geq 0$的的域流（信息流）决定，该投资组合剩余金额$X_t - \Delta_t S_t$ 投资于货币市场。

该投资组合随时间的变动$d X_t$由两部分构成：

（i）股票市场资本所得$\Delta_t (S_t + d S_t) - \Delta_t S_t = \Delta_t dS_t$ ,

（ii）货币市场利息收益$r (X_t - \Delta_t S_t) dt$，

即：
$$
\begin{equation}
\begin{aligned}
d X_t &= \Delta_t dS_t + r (X_t - \Delta_t S_t) dt \\
&= \Delta_t (\alpha S_t dt + \sigma S_t dW_t) + r (X_t - \Delta_t S_t) dt \\
&= r X_t dt + \Delta_t (\alpha - r) S_t dt + \Delta_t \sigma S_t d W_t
\end{aligned}
\tag{1.3}
\end{equation}
$$

> $d X_t$由以下三部分构成：
>
> 1. 投资组合$X_t$基于利率$r$的隐含固定收益$r X_t dt$；
> 2. 股票市场的风险溢价(risk premium $\alpha - r$)收益$\Delta_t (\alpha - r) S_t dt$；
> 3. $\Delta_t$份额股票的波动收益$\Delta_t \sigma S_t d W_t$；

下面我们考虑下按照连续复利折现后的股价$e^{-rt} S_t$和投资组合价值$e^{-rt} X_t$的微分。

> 假设$f(t, x) = e^{-rt}x$，则$f_t(t, x) = -r e^{-rt} x$，$f_x(t, x) = e^{-rt}$，$f_{xx}(t, x) = 0$，

根据伊藤公式
$$
d f(t, X_t) = f_t(t, X_t) dt + f_x(t, X_t) dX_t + \frac{1}{2} f_{xx}(t, X_t) dX_t dX_t \tag{1.4}
$$


可得：
$$
\begin{equation}
\begin{aligned}
d (e^{-rt} S_t) &= d f(t, S_t) \\
&= f_t(t, S_t) dt + f_x(t, S_t) dS_t + \frac{1}{2}f_{xx}(t, S_t) dS_t dS_t \\
&= -r e^{-rt} S_t dt + e^{-rt} dS_t \\
&= -r e^{-rt} S_t dt + e^{-rt} \alpha S_t dt + e^{-rt} \sigma S_t dW_t \\
&= e^{-rt} (\alpha - r) S_t dt + e^{-rt} \sigma S_t dW_t
\end{aligned}
\tag{1.5}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
d (e^{-rt} X_t) &= d f(t, X_t) \\
&= f_t(t, X_t) dt + f_x(t, X_t) dX_t + \frac{1}{2}f_{xx}(t, X_t) dX_t dX_t \\
&= -r e^{-rt} X_t dt + e^{-rt} dX_t \\
&= \left[-r e^{-rt} X_t dt + e^{-rt} r X_t dt\right] + e^{-rt} \Delta_t (\alpha - r) S_t dt + e^{-rt} \Delta_t \sigma S_t d W_t \\
&= \Delta_t e^{-rt} (\alpha - r) S_t dt + \Delta_t e^{-rt} \sigma S_t d W_t \\
&= \Delta_t d(e^{-rt} S_t)
\end{aligned}
\tag{1.6}
\end{equation}
$$

由于货币市场收益也是按照连续复利进行计算，所以对其进行折现后为初始值常数，不会对折现投资组合价值变动（微分）产生影响，此即公式$\text{(1.5)}$所反映的：折现投资组合价值变动$d (e^{-rt} X_t)$仅取决于折现股票价格变动$d (e^{-rt} S_t)$。

### 二、期权价值（价格）演化

我们考虑在行权日期$T$支付为$(S_T - K)^+$的欧式看涨期权，假设在时刻$t$其价格由函数$c(t, S_t)$决定，下面我们分别对对$c(t, S_t)$和其折现过程$e^{-rt} c(t, S_t)$求微分：
$$
\begin{equation}
\begin{aligned}
d c(t, S_t) &= c_t(t, S_t) dt + c_x(t, S_t) dS_t + \frac{1}{2} c_{xx}(t, S_t) dS_t dS_t \\
&= c_t(t, S_t) dt + c_x(t, S_t) [\alpha S_t dt + \sigma S_t dW_t] + \frac{1}{2} c_{xx}(t, S_t) \sigma^2 S_t^2 dt \\
&= \big[c_t(t, S_t) + \alpha S_t c_x(t, S_t) + \frac{1}{2} \sigma^2 S_t^2 c_{xx}(t, S_t)\big]dt + \sigma S_t c_x(t, S_t) dW_t
\end{aligned}
\tag{2.1}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
d e^{-rt} c(t, S_t) &= df\big(t, c(t, S_t)\big) \\
&= f_t\big(t, c(t, S_t)\big) dt \\
& + f_x\big(t, c(t, S_t)\big) dc(t, S_t) + \frac{1}{2} f_{xx}\big(t, c(t, S_t)\big) dc(t, S_t) dc(t, S_t) \\
&= -r e^{-rt} c(t, S_t) dt + e^{-rt} d c(t, S_t) \\
&= -r e^{-rt} c(t, S_t) dt \\
& + e^{-rt} \left\{ \big[c_t(t, S_t) + \alpha S_t c_x(t, S_t) + \frac{1}{2} \sigma^2 S_t^2 c_{xx}(t, S_t)\big]dt + \sigma S_t c_x(t, S_t) dW_t \right\} \\
&= e^{-rt} \big[ -rc(t, S_t) + c_t(t, S_t) + \alpha S_t c_x(t, S_t) + \frac{1}{2} \sigma^2 S_t^2 c_{xx}(t, S_t) \big]dt \\
& + e^{-rt} \sigma S_t c_x(t, S_t) dW_t
\end{aligned}
\tag{2.2}
\end{equation}
$$

### 三、投资组合价值与期权价值的等价

我们的投资始于在时刻$t=0$时卖出期权（short option）获得$c(0 ,S_0)$，然后投资于货币市场和股票市场形成投资组合$X_0$，满足及时立刻行权也无套利$c(0, S_0) = X_0 = (S_0 - K)^+$，同时为了保证我们的期权买卖交易在任意时刻$t, t \in (0, T]$无套利，还需要满足：

$$
\begin{aligned}
& c(t, S_t) = X_t \\\
& e^{-rt}c(t, S_t) = e^{-rt} X_t
\end{aligned}  \tag{3.1}
$$
> 例如我们可以假设存在一个时空穿越者B，其时间总比现实中的投资者A快$t$，也就是他总能提前准确获知未来的信息，但其投资行为也必须满足无套利原则，则他会告诉你上式必须满足，否则市场就会存在套利。比如$c(t, S_t) \geq X_t$，那么投资者B就可以告诉投资者A在时刻$t=0$卖出看涨期权，获得$e^{-rt}c(t, S_t)$，同时抽走$e^{-rt}c(t, S_t) - e^{-rt} X_t$，仅用$\left( e^{-rt}c(t, S_t) - \left[ e^{-rt}c(t, S_t) - e^{-rt} X_t \right]\right) == e^{-rt} X_t = X_0$投资于货币市场和股票市场进行对冲卖出看涨期权的交易行为，即投资者A在时刻$t=0$获得了无风险收益$e^{-rt}c(t, S_t) - e^{-rt} X_t$。

如满足：
$$
d\left(e^{-rt} X_t\right) = d \left(e^{-rt}c(t, S_t)\right) \quad \text{ for all } t \in [0, T] \tag{3.2}
$$
则上式必然成立。展开公式$\text{(3.2)}$并约减掉两边的$e^{-rt}$得到如下公式：
$$
\begin{equation}
\begin{aligned}
& \Delta_t (\alpha - r) S_t dt + \Delta_t \sigma S_t d W_t \\
& = \big[ -rc(t, S_t) + c_t(t, S_t) + \alpha S_t c_x(t, S_t) + \frac{1}{2} \sigma^2 S_t^2 c_{xx}(t, S_t) \big]dt \\
& + \sigma S_t c_x(t, S_t) dW_t
\end{aligned}
\tag{3.3}
\end{equation}
$$

首先令$dW_t$项系数相等，得到德尔塔对冲法则（***delta -hedging rule***）：
$$
\Delta_t = c_x(t, S_t) \quad \text{for all } t \in [0, T]. \tag{3.4}
$$
即在任意时刻$t$，持有股票份额$\Delta_t$等于期权价格对股价的偏导数$c_x(t, S_t)$。

然后再令$dt$项系数相等，并将公式$\text{(3.4)}$带入其中，得到：
$$
\begin{equation}
\begin{aligned}
& c_x(t,S_t) (\alpha - r) S_t = -rc(t, S_t) + c_t(t, S_t) + \alpha S_t c_x(t, S_t) \\
& + \frac{1}{2} \sigma^2 S_t^2 c_{xx}(t, S_t) \quad \text{for all } t \in [0, T]. 
\end{aligned}
\tag{3.5}
\end{equation}
$$

$$
-r c_x(t,S_t) S_t = -rc(t, S_t) + c_t(t, S_t) + \frac{1}{2} \sigma^2 S_t^2 c_{xx}(t, S_t) \tag{3.6}
$$

$$
rc(t, S_t) = c_t(t, S_t) + r S_t c_x(t,S_t) +  \frac{1}{2} \sigma^2 S_t^2 c_{xx}(t, S_t)  \quad \text{for all } t \in [0, T]. \tag{3.7}
$$

将$S_t$用哑变量$x$代替，便得到了BSM偏微分方程;


$$
\left\{
\begin{aligned}
& rc(t, x)  = c_t(t, x) + r x c_x(t,x) +  \frac{1}{2} \sigma^2 x^2 c_{xx}(t, x)   \quad \text{for all } t \in [0, T], x \geq 0.\\\
& c(T, x)  = (x - K)^+.
\end{aligned}  \tag{3.8}
\right.
$$

下面我们通过适当变换将公式$\text{(3.8)}$写成常见的热传导方程形式并进行求解：
$$
\left\{
\begin{aligned}
& \frac{\partial u}{\partial t} = a^2 \frac{\partial^2 u}{\partial x^2},   \quad t > 0,  -\infty < x < \infty.\\\
& u\vert_{t = 0} = \varphi(x).
\end{aligned}  \tag{3.9}
\right.
$$

### 四、BSM偏微分方程求解

我们先将BSM偏微分方程转换为常规的热传导方程形式，然后借助傅里叶变换对其进行求解。

#### BSM方程转换

为方便起见，将$\text{(3.8)}$作如下替换：$c(t, x)$简写成$V(t,S_t)$，并用偏导数$\frac{\partial V}{\partial S}$形式取代$c_x(t, x)$得到：
$$
\left\{
\begin{aligned}
& \frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + r S \frac{\partial V}{\partial S} - r V = 0 \quad 0 \leq S \leq \infty,  0 \leq t \leq T.\\\ 
& V\vert_{t =T} = (S_T - K)^+.
\end{aligned}  \tag{4.1}
\right.
$$
对$V(t,S_t)$做如下变量替换
$$
\begin{aligned}
& x = \ln S \\
& \tau = T-t
\end{aligned} \tag{4.2}
$$
计算$V(\tau, x)$的各阶偏导数如下：
$$
\frac{\partial V}{\partial t} = \frac{\partial V}{\partial \tau} \frac{\partial \tau}{\partial t} = -\frac{\partial V}{\partial \tau} \tag{4.3}\\
\frac{\partial V}{\partial S} = \frac{\partial V}{\partial x} \frac{\partial x}{\partial S} = \frac{1}{S} \frac{\partial V}{\partial x} \\
\frac{\partial^2 V}{\partial S^2} = \frac{\partial \frac{\partial V}{\partial S}}{\partial S} = -\frac{1}{S^2} \frac{\partial V}{\partial x} + \frac{1}{S} \frac{\partial^2 V}{\partial x^2} \frac{\partial x}{\partial S} = \frac{1}{S^2} \frac{\partial^2 V}{\partial x^2} - \frac{1}{S^2} \frac{\partial V}{\partial x}
$$
将上述公式带入$\text{(4.1)}$得到：
$$
\left\{
\begin{aligned}
& \frac{\partial V}{\partial \tau} - \frac{\sigma^2}{2} \frac{\partial^2 V}{\partial x^2} - (r - \frac{\sigma^2}{2}) \frac{\partial V}{\partial x} + r V = 0 \quad 0 \leq S \leq \infty,  0 \leq t \leq T.\\\ 
& V\vert_{\tau = 0} = (e^x - K)^+.
\end{aligned}  \tag{4.4}
\right.
$$
然后再做如下函数变化：
$$
V(\tau, x) = V = u e^{\alpha \tau + \beta x} = u(\tau, x) e^{\alpha \tau + \beta x} \tag{4.5}
$$

$$
\begin{aligned}
& \frac{\partial V}{\partial \tau} = e^{\alpha \tau + \beta x} \frac{\partial u}{\partial \tau} + e^{\alpha \tau + \beta x} \alpha u = e^{\alpha \tau + \beta x} (\frac{\partial u}{\partial \tau} + \alpha u) \\
& \frac{\partial V}{\partial x} = e^{\alpha \tau + \beta x} (\frac{\partial u}{\partial x} + \beta u) \\
& \frac{\partial^2 V}{\partial x^2} = e^{\alpha \tau + \beta x} (\frac{\partial^2 u}{\partial x^2} + 2 \beta \frac{\partial u}{\partial x} +\beta^2 u)
\end{aligned}  \tag{4.6}
$$

将以上公式带入$\text{(4.4)}$并消去$e^{\alpha \tau + \beta x}$得到：
$$
e^{\alpha \tau + \beta x} (\frac{\partial u}{\partial \tau} + \alpha u) - \frac{\sigma^2}{2} e^{\alpha \tau + \beta x} (\frac{\partial^2 u}{\partial x^2} + 2 \beta \frac{\partial u}{\partial x} +\beta^2 u) \\ - (r - \frac{\sigma^2}{2}) e^{\alpha \tau + \beta x} (\frac{\partial u}{\partial x} + \beta u) + r e^{\alpha \tau + \beta x} u = 0 \tag{4.7}
$$

$$
\frac{\partial u}{\partial \tau} - \frac{\sigma^2}{2} \frac{\partial^2 u}{\partial x^2} - (\beta \sigma^2 + r - \frac{\sigma^2}{2}) \frac{\partial u}{\partial x} + \left(r - \beta(r - \frac{\sigma^2}{2}) - \frac{\sigma^2}{2} \beta^2 + \alpha \right) u = 0
\tag{4.8}
$$



通过选取适当的$\alpha$和$\beta$的值，消去$\frac{\partial u}{\partial x}$和$u$项，即：
$$
(\beta \sigma^2 + r - \frac{\sigma^2}{2}) = 0 \quad \rightarrow \quad \beta = \frac{1}{2} - \frac{r}{\sigma^2} \tag{4.9}\\
\left(r - \beta(r - \frac{\sigma^2}{2}) - \frac{\sigma^2}{2} \beta^2 + \alpha \right) = 0 \quad \rightarrow \quad \alpha = -r + \beta(r - \frac{\sigma^2}{2}) + \frac{\sigma^2}{2} \beta^2 = -r - \frac{1}{2 \sigma^2}(r - \frac{\sigma^2}{2})^2
$$
则$\text{(4.8)}$变换为$\text{(4.10)}$形式的的常规热传导方程：
$$
\left\{
\begin{aligned}
& \frac{\partial u}{\partial \tau} = \frac{\sigma^2}{2} \frac{\partial^2 u}{\partial x^2},   \quad 0 \leq \tau \leq T,  0 < x < \infty.\\\
& u\vert_{\tau = 0} = e ^{-\beta x} V \vert_{\tau = 0} = e ^{-\beta x} (e ^x - K)
 ^+ .
 \end{aligned}  \tag{4.10}
\right.
$$
接下来便是求解该方程得到$u(\tau, x)$的解析表达式，然后利用$V(\tau, x) = u(\tau, x) e^{\alpha \tau + \beta x}$求得$V(\tau, x)$。

#### BSM方程求解

在偏微分方程（或者叫数学物理方程）理论中，有成熟的傅里叶变换方法可以求解上述热传导方程，为了简化傅里叶变换公式的呈现，我们将边界条件进行简写如下：
$$
\begin{equation}
\begin{aligned}
u(\tau, x)\vert_{\tau = 0} &= e ^{-\beta x} V \vert_{\tau = 0} \\
&= e ^{-\beta x} (e ^x - K)^+ \\
&= \varphi(x).
\end{aligned} \tag{4.11}
\end{equation}
$$
记$\mathcal{F}$为傅里叶变换算符，$\mathcal{F^\text{-1}}$为傅里叶逆变换算符，$u(\tau, x)$的傅里叶变换为：
$$
\widehat{u}(\tau, \lambda) = \mathcal{F}[u(\tau, x)] = \int^{+\infty}_{-\infty} u(\tau, x) e^{-i \lambda x} dx \tag{4.12}
$$
对$\text{(4.10)}$方程两边进行傅里叶变换得到**像函数**$\widehat u(\tau, \lambda)$满足的常微分方程初值问题：
$$
\left\{
\begin{aligned}
& \frac{\text{d} \widehat{u}}{\text{d}\tau} = \frac{\sigma^2}{2} (i \lambda)^2 \widehat{u} = -\frac{\sigma^2}{2} \lambda^2 \widehat{u},   \quad 0 \leq \tau \leq T,  0 < x < \infty.\\\
& \widehat{u}\vert_{\tau = 0} = \widehat \varphi(\lambda).
\end{aligned}  \tag{4.13}
\right.
$$
上述傅里叶变换将偏微分方程变成关于$\widehat u(\tau, \lambda)$的常微分方程。求解此方程可得像函数：
$$
\widehat u(\tau, \lambda) = \widehat \varphi (\lambda) e ^{-\frac{\sigma^2}{2} \lambda^2 \tau} \tag{4.14}
$$
对上述像函数进行傅里叶逆变换，并进行卷积运算可得：
$$
\begin{equation}
\begin{aligned}
u(\tau, x) &= \mathcal{F^\text{-1}} [\widehat \varphi (\lambda) e^{-\frac{\sigma^2}{2} \lambda^2 \tau}] \\
&= \mathcal{F^\text{-1}} [\widehat \varphi (\lambda)] * \mathcal{F^\text{-1}} [e^{-\frac{\sigma^2}{2} \lambda^2 \tau}]  \\
&= \varphi(x) * \frac{1}{\sigma \sqrt{2 \pi \tau}} \exp(-\frac{x^2}{2 \sigma^2 \tau}) \\
&= \frac{1}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{-\infty} \varphi(\xi) \exp(-\frac{(x-\xi)^2}{2 \sigma^2 \tau}) d\xi
\end{aligned}
\tag{4.15}
\end{equation}
$$
其中：
$$
\begin{equation}
\begin{aligned}
\mathcal{F^\text{-1}} [e^{-\frac{\sigma^2}{2} \lambda^2 \tau}] &= \frac{1}{2 \pi} \int^{\infty}_{-\infty} e^{-\frac{\sigma^2}{2} \lambda^2 \tau} e^{i \lambda x} d\lambda \\
&= \frac{1}{\pi} \int^{+\infty}_{0} e^{-\frac{\sigma^2}{2} \lambda^2 \tau} \cos(\lambda x) d\lambda \\
&= \frac{1}{\sigma \sqrt{2 \pi \tau}} \exp(-\frac{x^2}{2 \sigma^2 \tau})
\end{aligned} \tag{4.16}
\end{equation}
$$
将初值$\varphi(\xi)= e ^{-\beta \xi} (e ^\xi - K)^+$带入$\text{(4.15)}$得到：
$$
\begin{equation}
\begin{aligned}
u(\tau, x) &= \frac{1}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{-\infty} \exp(-\frac{(x-\xi)^2}{2 \sigma^2 \tau}) \left[e ^{-\beta \xi} (e ^\xi - K)^+ \right] d\xi \\
&= \frac{1}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{\ln K} \exp(-\frac{(x-\xi)^2}{2 \sigma^2 \tau}) \left[e ^{(1-\beta) \xi} - K e^{-\beta \xi} \right] d\xi
\end{aligned}
\tag{4.17}
\end{equation}
$$


然后将$\alpha = -r - \frac{1}{2 \sigma^2}(r - \frac{\sigma^2}{2})^2$和$\beta = \frac{1}{2} - \frac{r}{\sigma^2}$带入$V(\tau, x) = u(\tau, x) e^{\alpha \tau + \beta x}$得到：
$$
\begin{equation}
\begin{aligned}
V(\tau, x) &= \exp(-r \tau - \frac{1}{2 \sigma^2}(r - \frac{\sigma^2}{2})^2 \tau - \frac{1}{\sigma^2}(r - \frac{\sigma^2}{2}) x) u(\tau, x) \\
&= I_1 + I_2
\end{aligned}
\tag{4.18}
\end{equation}
$$
其中：
$$
\begin{equation}
\begin{aligned}
I_1 &= \frac{1}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{\ln K} \exp(-r \tau - \frac{1}{2 \sigma^2}(r - \frac{\sigma^2}{2})^2 \tau - \frac{1}{\sigma^2}(r - \frac{\sigma^2}{2}) x) \exp(-\frac{(x-\xi)^2}{2 \sigma^2 \tau}) \exp \bigg( {(1-\beta) \xi} \bigg) d\xi \\
&= \frac{1}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{\ln K} \exp(-r \tau - \frac{1}{2 \sigma^2}(r - \frac{\sigma^2}{2})^2 \tau - \frac{1}{\sigma^2}(r - \frac{\sigma^2}{2}) x) \exp(-\frac{(x-\xi)^2}{2 \sigma^2 \tau}) \exp \bigg( \left[ 1+\frac{1}{\sigma^2}(r - \frac{\sigma^2}{2}) \right] \xi \bigg) d\xi \\
&= \frac{e^{-r \tau}}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{\ln K} \exp( - \frac{1}{2 \sigma^2}(r - \frac{\sigma^2}{2})^2 \tau - \frac{1}{\sigma^2}(r - \frac{\sigma^2}{2}) x -\frac{(x-\xi)^2}{2 \sigma^2 \tau} + \left[ \xi +\frac{1}{\sigma^2}(r - \frac{\sigma^2}{2}) \xi \right] )d\xi \\
&= \frac{e^{-r \tau}}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{\ln K} \exp(-\frac{(x-\xi)^2}{2 \sigma^2 \tau} - \frac{1}{\sigma^2}(r - \frac{\sigma^2}{2})(x - \xi) - \frac{1}{2 \sigma^2}(r - \frac{\sigma^2}{2})^2 \tau + \xi) d\xi \\
&= \frac{e^{-r \tau}}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{\ln K} \exp(-\frac{1}{2 \sigma^2 \tau} \left[(x-\xi)^2 + 2 (r - \frac{\sigma^2}{2}) \tau (x - \xi) + (r - \frac{\sigma^2}{2})^2 \tau^2 \right] + \xi) d\xi \\
&= \frac{e^{-r \tau}}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{\ln K} \exp(-\frac{1}{2 \sigma^2 \tau} \left[(x-\xi) + (r - \frac{\sigma^2}{2}) \tau \right] + \xi) d\xi \\
\end{aligned}
\tag{4.19}
\end{equation}
$$
令$\eta = x-\xi + (r - \frac{\sigma^2}{2}) \tau$并带入$\text{(4.19)}$：
$$
\begin{equation}
\begin{aligned}
I_1 &= \frac{e^{-r \tau}}{\sigma \sqrt{2 \pi \tau}} \int^{-\infty}_{x - \ln K + (r - \frac{\sigma^2}{2}) \tau} - \exp(-\frac{\eta^2}{2 \sigma^2 \tau} + (x - \eta + (r - \frac{\sigma^2}{2}) \tau)) d\eta \\
&= \frac{1}{\sigma \sqrt{2 \pi \tau}} \int^{x - \ln K + (r - \frac{\sigma^2}{2}) \tau}_{-\infty} e^x \exp(-\frac{\eta^2}{2 \sigma^2 \tau} - \eta - \frac{\sigma^2}{2} \tau) d\eta \\
&= \frac{1}{\sigma \sqrt{2 \pi \tau}} \int^{x - \ln K + (r - \frac{\sigma^2}{2}) \tau}_{-\infty} e^x \exp(-\frac{\eta^2 + 2 \sigma^2 \tau \eta + (\sigma^2 \tau)^2}{2 \sigma^2 \tau}) d\eta \\
&= \frac{e^x}{\sigma \sqrt{2 \pi \tau}} \int^{x - \ln K + (r - \frac{\sigma^2}{2}) \tau}_{-\infty} \exp(-\frac{(\eta + \sigma^2 \eta)^2}{2 \sigma^2 \tau}) d\eta
\end{aligned}
\tag{4.20}
\end{equation}
$$
令
$$
\omega = \frac{\eta + \sigma^2 \tau}{\sigma \sqrt{\tau}} \tag{4.21}\\
N(x) = \frac{1}{\sqrt{2\pi}} \int^{x}_{-\infty} \exp(-\frac{\omega^2}{2}) d\omega
$$
带入$\text{(4.20)}$，求得：
$$
\begin{equation}
\begin{aligned}
I_1 &= \frac{e^x}{\sigma \sqrt{2 \pi \tau}} \int^{\frac{x - \ln K + (r + \frac{\sigma^2}{2}) \tau}{\sigma \sqrt{\tau}}}_{-\infty} \sigma \sqrt{\tau} \exp(-\frac{\omega^2}{2}) d\omega \\
&= e^x N(\frac{x - \ln K + (r + \frac{\sigma^2}{2}) \tau}{\sigma \sqrt{\tau}})
\end{aligned}
\tag{4.22}
\end{equation}
$$
同理
$$
\begin{equation}
\begin{aligned}
I_2 &= \frac{1}{\sigma \sqrt{2 \pi \tau}} \int^{+\infty}_{\ln K} - \exp(-r \tau - \frac{1}{2 \sigma^2}(r - \frac{\sigma^2}{2})^2 \tau - \frac{1}{\sigma^2}(r - \frac{\sigma^2}{2}) x) \exp(-\frac{(x-\xi)^2}{2 \sigma^2 \tau}) K e^{-\beta \xi} d\xi \\
&= -\frac{e^{-r\tau} K}{\sigma \sqrt{2 \pi \tau}} \int^{x - \ln K + (r - \frac{\sigma^2}{2}) \tau}_{-\infty} \exp(-\frac{\eta^2}{2 \sigma^2 \tau}) d\eta \\
&= -K e^{-r \tau} N(\frac{x - \ln K + (r - \frac{\sigma^2}{2}) \tau}{\sigma \sqrt{\tau}})
\end{aligned}
\tag{4.23}
\end{equation}
$$
根据$x = \ln S, \tau = T-t$，最终得到$V(t, S_t)$：
$$
V(t, S_t) = S_t N(\frac{\ln \frac{S_t}{K} + (r + \frac{\sigma^2}{2}) (T-t)}{\sigma \sqrt{T-t}}) - K e^{-r (T-t)} N(\frac{\ln \frac{S_t}{K} + (r - \frac{\sigma^2}{2}) (T-t)}{\sigma \sqrt{(T-t)}})
\tag{4.24}
$$
或者简写为：
$$
\begin{equation}
\begin{aligned}
V(t, S_t) &= S_t N(d_1) - K e^{-r (T-t)} N(d_2) \\
\\
& d_1 = \frac{\ln \frac{S_t}{K} + (r + \frac{\sigma^2}{2}) (T-t)}{\sigma \sqrt{T-t}} \\
& d_2 = d_1 - \sigma \sqrt{T-t} = \frac{\ln \frac{S_t}{K} + (r - \frac{\sigma^2}{2}) (T-t)}{\sigma \sqrt{(T-t)}}
\end{aligned}
\tag{4.25}
\end{equation}
$$
由此，我们便求得了BSM偏微分方程的解析解，亦即欧式看涨期权定价公式。

### 注记

以上便整理完了对期权定价公式的偏微分方程直接求解方法，仅做形式化推导（有些地方确实考验对常规数学分析知识的掌握及细心程度例如公式$\text{4.19-4.20}$）的推导）。对偏微分方程求解过程中所涉及的傅里叶变换的适用性不做讨论，傅里叶变换目前除了印象中可以将非病态函数分解成三角函数的叠加以外所知甚少，再者公式$\text(4.16)$第一行到第二行的转化，第二行含参变量积分的计算（曾大言不惭想直接借助分部积分推导下，后来发现没那么容易），还有就是《Stochastic Calculus for Finance》关于对伊藤积分采用黎曼分割的合理性问题等等，总之想要彻底搞清楚其中的细节，还需不忘初心，继续精进。

