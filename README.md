# Riemann-hypothesis
An implementation of Riemann's hypothesis in R, exploring the relationship of Riemann's zeta function: 

$$\zeta(s) = \sum_n n^{-s} \tag{1}$$

to prime numbers through his prime counting function:

$$\pi(x) = \sum_n \frac{\mu(n)}{n} J(\sqrt[n]{x}) \tag{2}$$

where $\mu(n)$ is the Möbius function, and $J(x)$ incorporates the zeroes of the zeta function ($\rho$): 

$$J(x) = \operatorname{li}(x) - \sum_{\rho}^{\infty}\operatorname{li}(x^{\rho}) - \operatorname{log}2 + \int_{x}^{\infty}\frac{1}{t(t^{2}-1)\operatorname{log}t}dt \tag{3}$$

[comment]: <> ($$\pi(x) = \operatorname{R}(x) - \sum_{\rho}\operatorname{R}(x^{\rho}) \tag{3}$$)

