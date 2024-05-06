# Shepherding/herding control problem

The code in this folder reproduces the shepherding dynamics as discussed in [LAMA2023]. All the codes are written in MATLAB.
The shepherding (herding) control problem consists of a group of agents, the *herders*, that need to coordinate so as to steer the collective dynamics of a second group of agents, the *targets*.
Specifically, in [LAMA2023] it is discussed a 2d case where the herders need to steer to a prescribed goal region and contain therin a group of targets. The herders achieve this result by i) selecting one target to chase according to their distance from the goal region and ii) steering the selected target towards the goal region exploiting a local repulsion that the herders can exert on close enough targets.

The role of all the variables in the codes is commented in detail therin, and the nomenclature is as coherent as possible with that used in [LAMA2023] (some variables may have the name of built-in matlab functions).

- main.m function launches a simulation of the shepherding problem, where herders can have (A) finite or (B) infinite sensing capabilities. 
- (A) Shepherding_finXi.m   simulates a shepherding dynamics where herders have finite sensing radius (passed as a parameter)
- (B) Shepherding_infXi.m   simulates a shepherding dynamics where herders have infinite sensing radius 

Both "Shepherding_finXi.m" and "Shepherding_infXi.m" produce a video analogous to the ones reported in the supplementary material of [LAMA2023], and also save some output quantities that are relevant to carry out the analysis discussed in [LAMA23].

The remaining functions in the folder are "auxiliary" functions
  - initial_pos_circle.m   generates the initial condition where the agents are randomly and uniformly distributed in a circle around the origin
  - cartesian_to_pol.m     converts an array of 2d cartesian coordinates into an array of 2d polar coordinates
  - chi.m computes the     fraction of target agents that are within the goal region
  - saturation.m limits    herders' maximum speed to a prefixed value

References:
- [LAMA2023](https://arxiv-org.translate.goog/abs/2307.16797?_x_tr_sl=en&_x_tr_tl=it&_x_tr_hl=it&_x_tr_pto=sc) Lama, A., Di Bernardo, M. (2023) "Shepherding control and herdability in complex multiagent systems", arXiv:2307.16797

    
    
    
