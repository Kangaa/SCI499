Abstract
1. Intro and Background
	1. [[Disease spread]]
	2. [[Compartmental models]]
		1. [[SIR]]
		2. Structure (compartments)
		3. Parameters
			1. $\beta$ Transmission coefficient 
			2. $\gamma$ recovery coefficient
		3. Important properties/results
			1. $R_{0}$
			2. Final Size
			3. Duration
	3. Deterministic ([[Ordinary Differential Equation|ODE]]) implimentation
		1. Equations
		2. limitations
	4. [[Stochastic Process]] ([[Markov Chain]]) formulation
		1. [[Stochastic process]] definition
		2. [[Markov chain]]
		3. [[Discrete time Markov chain]] (DTMC)
			1. Transition matrix
		4. [[ Branching Process]] ?
		5. [[Continuous time Markov chain]] (CTMC)
			1. Transition rates
			2. [[Continuous Time Markov Chain#Infitetestimal Generator (Q)]]
			4. $P(t) = e^{Qt}$
		8. Stochastic simulation
			1. [[Gillespie algorithm]] for homogenous mixing model
	5. Example: Melbourne region compartmental SIR model
2.  Spatially structured Metapopulation models
	1. [[Metapopulation models]]
		1. Motivation
		2. Mixing matrix
	3. Spatial structure in mixing matrix
		1. e.g. transport info
	4. Case study: Melbourne city (replication if Moss et. al )
3. Hierarchically Structured metapopulation model
	2. Hierarchical mixing matrix
		1.  [[SA classification]]
			1. Hierarchical structure 
			2. Tree representation
		2. Population based models
			1. Pop-prop mixing matrices
			2. SA2-SA4 and homogenous comparison
		3. HMM construction
			2. 4. [[Hierarchical mixing matrices]]
			3. Population scaling
			4. [[Revised Spatial Mixing Matrix definition]]
			5. Case Study: Inner east area ()
	3. Hierarchical spatial models
		1. Specification & implimentation
		2. Results
			1. intra-level mixing parameter $\xi$
			2. effect of mm-blending parameter ?
			3. comparisons between models with differen spatial resolution
	