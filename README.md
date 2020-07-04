# Riemann-hypothesis
An implementation of Riemann's hypothesis in R, exploring the relationship of Riemann's zeta function: 

$$\zeta(s) = \sum_n n^{-s} \tag{1}$$

to prime numbers through his prime counting function:

$$R(x) = \sum_n \mu(n) \frac{J(\sqrt[n]{x})}{n} \tag{2}$$

where $\mu(n)$ is the Möbius function. And, specifically, how the zeroes of the zeta function ($\rho$) improve the accuracy of Riemann's prime counting function: 

$$\pi(x) = \operatorname{R}(x) - \sum_{\rho}\operatorname{R}(x^{\rho}) \tag{3}$$