# IRT_Intervention
Bayesian Item Response Model on Intervention Data

## IRT Model
This is a Probit Regression model. For $i = 1,\dots,N, \ k = i,\dots K$
\begin{align*}
Z_{ik} \sim N(\lambda_k^TM_k\theta_i - b_k, 1) \\
Y_{ik} &=
\left\{
\begin{array}{cc}
1, & Z_{ik} > 0 \\
0, & Z_{ik} \leq 0
\end{array}
\right.
\end{align*}
where $\lambda_k$ and $\theta_i$ are $3 \times 1$ vectors and $b_k$ is scalar, and $M_k$ is a diagonal matrix with binary elements \par

The priors on parameters are
\begin{align*}
\theta_i  \overset{\text{i.i.d.}}{\sim} & N_3 (0, I_3) \\  
\lambda_k \overset{\text{i.i.d.}}{\sim} & N_3 (0, I_3) \\  
b_k      \overset{\text{i.i.d.}}{\sim} & N (0, 1)
\end{align*}

## Intervention Data


## Model Diagnostic

