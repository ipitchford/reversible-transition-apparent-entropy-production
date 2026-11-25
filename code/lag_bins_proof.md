# Finite-speed certificates from waiting-time bins

Take the unit-period triangle described in the manuscript, with fast rate $R q$, $q=1+(16/17)\cos\theta$, $h=4/5$, $a=2\pi(1+h\sin2\theta)$, $b=2\pi(1-h\sin2\theta)$, and $p=1/2+h(\sin2\theta-\cos2\theta)/4$. The reversed experiment sends $h$ to $-h$. Angular brackets denote uniform phase average. The reporter event phase has density $2pq$ relative to that measure, and the reporter event rate is $R/2$.

## A common witness for either observation scheme

Write $P$ for the normalised forward law of a consecutive ordered pair of marks and its intervening waiting time, and $Q$ for the reverse-experiment law with the two marks reversed and exchanged. In this example every directed mark and its reversed-experiment counterpart have the same mean occurrence rate. The two marginals of the pair law agree in each experiment by stationarity under a one-event shift. Consequently the estimator in Eq. (14) of Harunari, Fiore and Barato equals $K D(P\Vert Q)$, where $K$ is the activity of all observed jumps. This holds for the single reporter edge and for all edges.

Select positive bin weights $w_i$, and assign $w_i$ to a pair of consecutive reporter marks when the scaled lag $u=R\tau$ falls in bin $i$. Assign weight one to every pair containing a slow mark. The pair test is invariant under reversed exchange of the marks. The entropy variational inequality and $\log z\le z-1$ imply
$$
K D(P\Vert Q)\ge K\{E_P\log w-E_Q(w-1)\}.
$$
Only pairs starting at a reporter event contribute. Thus we may average both right-hand terms conditional on starting with a reporter and multiply by $R/2$. Under full observation, a slow event occurring first sets the score to zero. The test already incorporates the exact distribution of directions through the equality of mark activities.

## Elementary certificate for a fast reporter

The script `verify_lag_bins_exact.py` uses the eight finite thresholds $1/4,1/2,1,2,3,5,8,12$, with nine weights $(90,94,101,114,132,132,103,77,56)/100$.

For the comparison experiment, keep the reporter enabled after the initial event, independently of the slow rates. Given phase $s$ in unit-period time, this reporter has survival
$$
\exp\!\left[-R\int_0^{u/R}q(s+v)\,dv\right].
$$
Its difference in log survival from the frozen reporter with survival $e^{-u q(s)}$ is at most $\epsilon=\pi(16/17)u^2/R$. Multiplying lower and upper enclosures of the frozen mean survival by $1-\epsilon$ and $1/(1-\epsilon)$ therefore bounds its true mean survival whenever $\epsilon<1$.

Write $L(s)$ for the unscaled mean waiting time to this always-enabled reporter given initial phase $s$. The periodic backward equation $L'=RqL-1$ gives $\langle qL\rangle=1/R$. The probability that a slow exit precedes the comparison reporter, averaged over the phase of a typical initial reporter event, is at most
$$
2 b_{\max}\langle pqL\rangle\le 2 b_{\max}(4/5)/R.
$$
The bound $p\le4/5$ follows from $\sqrt2<3/2$. Couple the actual observation to this comparison reporter up to its first slow exit. The expected score then differs by at most this exit probability times the range of the score. For the all-edge observation, include the neutral score zero in that range. For the chosen weights zero already lies inside both score ranges. This argument applies separately to the logs in the forward arm and the linear weights in the reverse arm.

The remaining finite enclosures use rational arithmetic. For each bin boundary the frozen survival is $\langle2pq e^{-uq}\rangle$. The degree-111 and degree-112 Taylor polynomials sandwich the exponential on the nonnegative ray. Their phase averages follow from $\langle\cos^{2j}\theta\rangle=\binom{2j}{j}/4^j$, vanishing odd moments, and $\cos2\theta=2\cos^2\theta-1$. For a step score with values $v_i$, its expectation is $v_0+\sum_i (v_{i+1}-v_i)\Pr[U>u_i]$. The code takes the appropriate endpoint of each survival interval according to the sign of its coefficient. The logarithms are bounded below by truncating the power series for $2\operatorname{atanh}[(w-1)/(w+1)]$ and subtracting a geometric bound on the absolute tail. It uses $\pi<22/7$ throughout.

At $R=100000$, the certified common bound is greater than 366.44498 per unit time. The ratio exceeds 145.74516 even with 88/35 in the denominator as an upper bound on the actual dissipation.

## Directly validating the dynamics at moderate speed

The second script `verify_moderate_speed.py` directly computes enclosures of finite-time event probabilities. It uses finite bin thresholds $1/4,1/2,3/4,1,3/2,2,3,4,6,10,16$, with respective weights $(90,94,99,104,111,121,132,136,121,90,67)/100$ for the bins ending at these thresholds. Unresolved waits beyond 16 receive unit weight (zero score).

At scaled time $u=R\tau$, the phase evolves with angular speed $\delta=2\pi/R$. Let $D$ denote the active doublet and $C$ the inactive singleton. The killed backward generator acts on periodic pairs of functions as
$$
 A f = \delta\partial_\theta f+
 \begin{pmatrix}-q-b/R & b/R\\ a/R&-a/R\end{pmatrix} f.
$$
The reporter jump kills the process. For the all-edge statistic a slow exit from the doublet also kills the test. In that case use the scalar generator $A=\delta\partial_\theta-q-b/R$ on the doublet. Write $g=(q,0)^T$, or simply $q$ in the scalar case. Given a starting doublet state at phase $\theta$, the probability of a reporter event by $u$, with any required adjacency, is the doublet component of
$$
 F(u)=\int_0^u e^{vA} g\,dv.
$$
Here $e^{vA}$ is the stopped-process semigroup and contracts the supremum norm. For every initial state, the expectation of a bounded function up to survival is at most its sup norm in absolute value.

Put $U=16$, and build the trigonometric polynomials
$$
 V_n=\frac{U^{n+1}}{(n+1)!}A^n g,\qquad
 V_{n+1}=\frac{U}{n+2}A V_n.
$$
Taylor's formula with integral remainder and contraction give, uniformly over $0\le u\le U$,
$$
 \left\|F(u)-\sum_{n=0}^{N-1}(u/U)^{n+1}V_n\right\|_\infty
 \le (u/U)^{N+1}\|V_N\|_\infty.
$$
Each $A^ng$ is a trigonometric polynomial in the domain of every further power of the generator. The calculation therefore requires only finite differentiations. The script encloses the sup norm by the sum of the absolute real and imaginary parts of the Fourier coefficients and stops when the entire norm enclosure is below $10^{-12}$.

The recurrence differentiates Fourier mode $k$ by multiplication with $ik$ and multiplies trigonometric polynomials by finite convolution. The code retains every generated mode, and projects $V_n$ against $2pq$ before evaluating the Taylor polynomial. Real and imaginary coefficient intervals use signed integers with fixed denominator $2^{224}$. Every arithmetic operation rounds both endpoints outward. Machin's identity for $π$, with adjacent alternating-series enclosures of the two arctangents, encloses the phase drift. Bin mass intervals follow by subtracting cumulative probabilities. Rational bounds on the scores and the entropy variational inequality finish the certificate. No ordinary differential equation solver or floating-point operation enters the proof.

At $R=750$, the certified one-edge lower bound exceeds 2.746942 and the all-edge bound exceeds 2.755619. At $R=1000$, they exceed 3.654098 and 3.661074, respectively. The exact physical dissipation is $4\pi/5<88/35$. Each comparison therefore gives a finite-speed violation.
