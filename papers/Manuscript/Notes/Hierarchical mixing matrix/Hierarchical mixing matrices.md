
we assume that each region shares the proportion of infections that occur in the home region, $\delta_{H}$, and those that occur from outside the region, $\delta^{A} = 1 - \delta^{H}$, to produce the matrix to produce the single level mixing matrix

Similarly, $\delta^A$, the proportion of infections occurring outside an individuals home region can be divided into those infections which occur inside the home region at the next highest scale, $\delta^{H_2}$, and those that occur at outside this higher level home region $\delta^{A_2}$, to create the two level mixing matrix:

$$ M = \begin{pmatrix}  
        \delta^{H} & \delta^{A}\cdot L_{1,2}^{2} & \cdots & \delta^{A}\cdot L_{1,j}^{2} \\          
        \delta^{A}\cdot L_{2,1}^{2} & \delta^{H} & \cdots & \delta^{A}\cdot L_{2,j}^{2} \\         \vdots & \vdots & \ddots & \vdots \\         \delta^{A}\cdot L_{j,1}^{2} & \delta^{A}\cdot L_{j,2}^{2} & \cdots & \delta^{H} \end{pmatrix} $$

where

$$ L_{i,j}^{k} = \begin{cases} 
 \frac{\delta^{H}}{N_{j,k}} & \mbox{if } j \in R^{k}_{i} \\ 
\frac{\delta^{A}}{N_{j,k}} & \mbox{if } j \notin R^{k}_{i} \mbox{ \& } k = n \\
\delta^{A}L_{i,j}^{k+1} & otherwise  \\
\end{cases} $$

Where $R_{i}^{k}$, is the set of patches in the same $k$ level region as $i$, and $N_{j,k}$ is the number of patches in the same level $k$ region as $j$, and $n$ is the total number of levels. 

The sum of each row is given by

$$
\sum_{j=1}^{N} M_{i,j} = \delta^{H} + \delta^{A}\delta^{H} + \cdots + \delta^{A^{n-1}}\delta^{H} + \delta^{A^{n}}
$$
Which factorises to 

$$
\sum_{j=1}^{N} M_{i,j} = \delta^{H} + \delta^{A}( \delta^{H} + \delta^{A}\delta^{H} +\cdots + \delta^{A^{n-2}}\delta^{H} + \delta^{A^{n-1}})
$$
iteratively until 
$$
\sum_{j=1}^{N} M_{i,j} = \delta^{H} + \delta^{A}( \delta^{H} + \delta^{A}(\delta^{H} + \delta^{A}(\cdots ( \delta^{H} + \delta^{A})))
$$

since $\delta^{H}+ \delta^{A} =1$, then 
$$
\sum_{j=1}^{N} M_{i,j} = \delta^{H} + \delta^{A} = 1
$$

