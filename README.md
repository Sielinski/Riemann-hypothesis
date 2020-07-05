# Riemann-hypothesis
An implementation of Riemann's hypothesis in R, exploring the relationship of Riemann's zeta function: 

<img src="https://latex.codecogs.com/gif.latex?\zeta(s)%20=%20\sum_n%20n^{-s}"/>

to prime numbers through his prime counting function:

<img src="https://latex.codecogs.com/gif.latex?\pi(x)%20=%20\sum_n%20\frac{\mu(n)}{n}%20J(\sqrt[n]{x})"/>

where μ(n) is the Möbius function, and J(x) incorporates the zeroes of the zeta function (ρ): 

<img src="https://latex.codecogs.com/gif.latex?J(x)%20=%20\operatorname{li}(x)%20-%20\sum_{\rho}^{\infty}\operatorname{li}(x^{\rho})%20-%20\operatorname{log}2%20+%20\int_{x}^{\infty}\frac{1}{t(t^{2}-1)\operatorname{log}t}dt"/>

[comment]: <> ($$\pi(x) = \operatorname{R}(x) - \sum_{\rho}\operatorname{R}(x^{\rho}) \tag{3}$$)

