%%%%% Starts ONE shepherding  simulation %%%%%
% The outputs are 
% i) a simulation video (not saved) analogous to the Supplementary Video 1
% and 2 of the paper, and
% ii) a file  "R_t_M_N_xi.mat" containing the values of the parameters used for the simulation and
% the values in time of the fraction of targets in the goal region and the
% average distance and relative dispersion of the herders from the origin.

% All the variables are named coherently with those presented in the paper,
% unless the name coincides with a built-in MATLAB function

%%%% INPUT PARAMETERS MEANING %%%%%%%
%  M is the number of targets and N is the number of herders.
% At t=0 the agents are randomly and uniformly distributed in a circle of
% radius R.
% The total duration (in arbitrary units) of the simulation is t
% xi is the sensing radius of the herders. If the input xi is a number you
% simulate the finite sensing case with the size of herders' sensing radius being
% the one specified. If the input xi is whatever string ("trial",
% "infinite", etc...), then you will simulate the case of herders with
% infinite sensing radius.


close all

tStart=cputime;     %%% Initialize CPU time



%%% Parameters of the dynamics of the targets
M=50;               % Number of targets
D=1;                % Regulates noise strength
beta=3;             % Regulates the strenght of the repulsion exerted onto a target by one herder
lambda=2.5;         % Maximum distance at which a target is repelled by a herder

%%% Parameters of the dynamics of the herders
N=10;               % Number of herders
v_H=7.5;            % Maximum speed achieved by one herder
a=5;                % $\alpha$ in the paper. Regulates the attraction strenght one herder feels with respect to an observed target
g=10;               % $\gamma$ in the paper. It is used to reproduce the selection rule via a weighted average.
delta=lambda*.5;    % Chasing distance at which each herder tries to place itself behind the selected target.
xi=10;              % Sensing radius of the herders
% xi="infi"          % if xi is a string then the simulation will consider herders endowed with unlimited sensing

%%% 
R=50;               % Size of the circle where the agents are distributed at t=0
t=100;              % Duration of the simulation in arbitrary units
dt=.01;                     % Time step of the integration scheme
t_save=.1;                  % Sampling time (in absolute units) of the snapshots plotted in the output shepherding animation video              
rg=R/5;                     % Size of the goal region centered arounf the origin where you want to collect the targets                     


if class(xi)=="string"
    Shepherding_infXi(N,M,R,t,dt,v_H,a,g,delta,lambda,beta,D,t_save,rg);
else
    Shepherding_finXi(N,M,R,t,dt,v_H,a,g,delta,xi,lambda,beta,D,t_save,rg)
end

TComp=cputime-tStart %%% Prints elapsed CPU time

