---
layout: post
title: "关于风险中性测度下Black-Scholes-Merton Model推导的总结"
tags: 金融数学
date: 2021-05-27 22:00:00 +0545
categories: Finance
comments: 1
---

2008年一次考证经历，第一次知道了期权并看到期权定价的Black-Scholes-Merton公式开始，被其深深迷住，断断续续翻阅了许多书籍及文献资料，包括Hull那本及基础的数学分析、实分析及基于测度论的概率论等等，以其能踩到看懂其背后原理的垫脚石。但由于大学专业和数学不沾边，看起来很吃力，直到工作后于2015年左右遇见世图影印出版的《Stochastic Calculus for Finance》才看到完整的讲解，似懂非懂中看到第5章开头就因种种原因放弃了，时间一直来到了2021年，趁COVID-19疫情期间工作之余继续怼，从头开始又来到了第5章，这次感觉仿佛看明白了许多。。。

因此做个总结，希望若干年以后还能不忘初心，顺着笔记还能摸到瓜。🍉


#### 1. 风险中性测度下的股价

假设 $W(t), 0 \leq t \leq T $是概率空间$ (\Omega, \mathcal{F}, \mathbb{P})$上的布朗运动，$\mathcal{F}(t), 0\leq t \leq T$是该布朗运动的域流，考虑股价$S(t)$，其微分如下：

$$
d S(t)=\alpha(t) S(t) d t+\sigma(t) S(t) d W(t), \quad 0 \leq t \leq T.  \tag{1.1}
$$

其中平均回报率$\alpha(t)$和股价波动率$\sigma(t)$为适应性过程，则$S(t)$满足以下等式：

$$
S(t)=S(0) \exp \left\{\int_{0}^{t} \sigma(s) d W(s)+\int_{0}^{t}\left(\alpha(s)-\frac{1}{2} \sigma^{2}(s)\right) d s\right\} \tag{1.2}
$$

假设我们有适应性利率过程$R(t)$，定义折现过程

$$
D(t)=e^{-\int_{0}^{t} R(s) d s} \tag{1.3}
$$

则

$$
d D(t)=-R(t) D(t) d t \tag{1.4}
$$

$D(t)S(t)$及其微分$d (D(t)S(t))$分别为：

$$
D(t) S(t)=S(0) \exp \left\{\int_{0}^{t} \sigma(s) d W(s)+\int_{0}^{t}\left(\alpha(s)-R(s)-\frac{1}{2} \sigma^{2}(s)\right) d s\right\} \tag{1.5}
$$

$$
\begin{aligned}
d(D(t) S(t)) &=(\alpha(t)-R(t)) D(t) S(t) d t+\sigma(t) D(t) S(t) d W(t) \\
&=\sigma(t) D(t) S(t)[\Theta(t) d t+d W(t)]
\end{aligned} \tag{1.6}
$$

其中定义风险的市场价格$\Theta(t)=\frac{\alpha(t)-R(t)}{\sigma(t)} $。

根据Girsanov's Theorem，在概率测度$\widetilde{\mathbb{P}}$下，$d \widetilde W(t) = \Theta(t) d t+d W(t) $，因此公式1.6也可以写为

$$
d (D(t)S(t)) = \sigma(t) D(t) S(t) d \widetilde W(t) \tag{1.7}
$$

两边同时积分

$$
D(t) S(t)=S(0)+\int_{0}^{t} \sigma(u) D(u) S(u) d \widetilde{W}(u) \tag{1.8}
$$

由于在概率测度$ \widetilde{\mathbb{P}}$下，$\int_{0}^{t} \sigma(u) D(u) S(u) d  \widetilde{W}(u)$是伊藤过程，因此是一个鞅。由此我们称Girsanov's Theorem下的概率测度$\widetilde{\mathbb{P}}$为风险中性测度（risk-neutral measure）。

将$d \widetilde W(t) = \Theta(t) d t+d W(t) $带入公式1.1，可以得到在概率测度$\widetilde{\mathbb{P}}$下，公式1.1和1.2分别可以改写成公式1.9和公式1.10的形式

$$
d S(t)=\alpha(t) S(t) d t+\sigma(t) S(t) d \widetilde {W}(t) \tag{1.9}
$$

$$
S(t)=S(0) \exp \left\{\int_{0}^{t} \sigma(s) d \widetilde{W}(s)+\int_{0}^{t}\left(R(s)-\frac{1}{2} \sigma^{2}(s)\right) d s\right\} \tag{1.10}
$$

#### 2. 未完待续...



感谢妈咪叔LaTeX助力 [在线LaTeX](https://www.latexlive.com)

参考:
				[金融随机分析](https://book.douban.com/subject/3369848/)
				[Stochastic Calculus for Finance II - Continuous-Time Models](https://book.douban.com/subject/2187921/)











