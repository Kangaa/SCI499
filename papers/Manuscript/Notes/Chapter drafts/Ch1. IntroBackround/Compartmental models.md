* A [[compartmental model]] treat a individuals as belonging to finite number of discrete states (Diekmann et al. 2009)
* For example, the SIR compartmental model is a common epidemiological model that describes individuals as belonging to one of three potential states: Susceptible (S),  Infected (I), or Recovered (R). 
* Individuals in this model can only move between compartments in a specified way. In particular, susceptible individuals can become infected, and infected individuals can recover.  Thus the model can be represented by  the flow diagram (fig.)
```mermaid
flowchart LR
    S --> I
    I --> R
```

* We can specify parameters that govern the probability that an individual will move between compartments
| $\beta$| probability of infection|
