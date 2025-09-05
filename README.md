# wwCRP.R
$w$-way coupled representative points based on Hybrid energy distance criterion
 w=1 +++++ Marginally coupled representative points (default)
 w=2 ----- Doubly coupled representative points
 w=q ***** Sliced representative points


# uniformcriteriaformixtures.R
uniform criteria: 
mean ($\Delta_{\mu}(\mathcal{P})=\|\boldsymbol{\hat{\mu}}(\mathcal{P})-\boldsymbol{\mu}\|_2 $)  standard deviation ($ \Delta_{\sigma}(\mathcal{P})=\|\boldsymbol{\hat{\sigma}}(\mathcal{P})-\boldsymbol{\sigma}\|_2 $) 
root mean squared distance ($\mathrm{RMSD}(\mathcal{P})=\sqrt{\frac{1}{N}\sum_{t=1}^N\underset{\boldsymbol{x}\in\mathcal{P}}{\min}\|\boldsymbol{u}_t
 -\boldsymbol{x}\|_2^2}$); 
 maximum packing distance ($\mathrm{MaD}(\mathcal{P})=\max_{t\in\{1,2,\dots,N\}}\underset{\boldsymbol{x}\in\mathcal{P}}{\min}\|\boldsymbol{u}_t-\boldsymbol{x}\|_2$); minimum separation distance ($\mathrm{MiD}(\mathcal{P})=\underset{\boldsymbol{x}_i,\boldsymbol{x}_j\in\mathcal{P}\atop \boldsymbol{x}_i\neq\boldsymbol{x}_j}{\min}\|\boldsymbol{x}_i-\boldsymbol{x}_j\|_2$)

# clcv.R    
necessary function in wwCRP.R, among each of the w variables selection, rows that have same level combination as row i

# classify_rows.R

# findindex.R
 for row i, find_identical_rows_with_submatrix_include_self 
