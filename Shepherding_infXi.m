function Shepherding_infXi(N,M,R,t,dt,v_H,a,g,delta,lambda,beta,D,t_save,rg)


t_start=cputime;

%%%% UNCOMMENT HERE AND AT LINE 191 TO SAVE THE VIDEO ( .avi format with
%%%% name file_name) %%%%%%%%%
% axis square
% axis tight manual
% file_name=sprintf("shepherding");
% v = VideoWriter(file_name);
% v.FrameRate=10;
% open(v);

%%%%% Check on the time
step_save=floor(t_save/dt);

if mod(t,t_save)==0
    number_save=t/t_save+1;
else
    fprintf("t must be divisible by t_check")
    return

end
%%% Name of the .mat file  the  output data is saved
mkdir("Data_infinite_xi")
filename=sprintf("Data_infinite_xi/%d_%d_%d_%d",round(R),round(t),M,N);

%%% Arrays where the output data is saved
T_final=zeros(M,2);
H_final=zeros(N,2);

active_herders  =zeros(number_save,1);  % Number of active herders that are chasing at least one target
chi_r           =zeros(number_save,1);  % Fraction of targets in the goal region (circle of radius rg)
chi_rDr         =zeros(number_save,1);  % Fraction of targets in the goal region accounting for some buffer distance (circle of radius rgT)
avg_rhoH        =zeros(number_save,1);  % Avg distance of the herders from the origin
std_rhoH        =zeros(number_save,1);  % Std of the distance of the herders from the origin
time            =zeros(number_save,1);  % Array of sample times
success_indicator=0;                    % return 1 if shepherding is successfull within a time t, 0 otherwise

drg=rg*.5;
rgT=rg+drg;
dR=2*sqrt(R*D/(v_H));

time_steps=round(t/dt);                 %%% Number of time integrations

%%% For plotting
th          =   0:2*pi/100:50;
goal_color  =   [0 0.4470 0.7410];



active_herder_counter=0;
active_herder_idx=false(N,1);

%%% Generating the initial conditions of N herders and M targets randomly
%%% and uniformly distributed  in a circle of radius R centered around the
%%% origin

H=initial_pos_circle(N,R);
T=initial_pos_circle(M,R);

% Initialization of the variable that counts the amout of sample times
% saved
save_counter=1;

% Initialize the array storing the output data

chi_r(save_counter)=chi(T,rg);
chi_rDr(save_counter)=chi(T,rgT);
active_herders(save_counter)=N;
rho=vecnorm(H,2,2);
avg_rhoH(save_counter)=mean(rho);
std_rhoH(save_counter)=std(rho);
time(save_counter)=0;


for i=1:time_steps %for over time integration

    right_to_herd=zeros(M,1);   % Specifies which herder is tasked to chase each target

    i_T=zeros(M,2);             % Array containing the total repulsion exerted by all the herders on each target

    i_H=zeros(N,2);             % Array containing the total attraction towards observed targets for each herder. $I$ in Eq.2 of the paper.
    i_H_num=zeros(N,2);         % Array containing the numerator of $I$  computed via the weighted average according to Eq. 5 of the paper.
    i_H_den=zeros(N,1);         % Array containing the denominator of $I$  computed via the weighted average according to Eq. 5 of the paper.
    f_H=zeros(N,2);             % Array containing the term $F$ in Eq. 2 of the paper, describing the dynamics of herders that are not chasing any target.

    active_herder_counter=0;        % Counts how many herders are chasing at least one target
    active_herder_idx=false(N,1);   % Logical array where the active herders are labeled with 1

    diff_HT_IH=zeros(M,2);          % Array containing the distances between each target and the closest herder

    for it=1:M          % for over the M targets

        diff_HT_IT=H-T(it,:);                   % Distances between target it and all the herders
        dHT_IT=vecnorm(diff_HT_IT,2,2);

        [~,right]=min(dHT_IT);                  % For each target, the herder tasked to chase it is the closest one in the infinite sensing case
        right_to_herd(it)=right;

        diff_HT_IH(it,:)=diff_HT_IT(right,:);   % Save the distance between target it and the closes herder, used later to compute herders' feedback

        idx=dHT_IT<=lambda;                     % Which herders are close enough to target it to repell it

        i_T(it,1)=-beta*dot((lambda-dHT_IT(idx)),(diff_HT_IT(idx,1)./dHT_IT(idx)));
        i_T(it,2)=-beta*dot((lambda-dHT_IT(idx)),(diff_HT_IT(idx,2)./dHT_IT(idx)));

    end



    for ih=1:N           % for over the N herders

        dHerdGoal=H(ih,:);                          % Distance of herder ih from the origin
        distHerdGoal=vecnorm(dHerdGoal,2,2);

        idx=right_to_herd==ih;                      % The targets that herder ih has to chase are the ones for which ih is the closest herder

        if any(idx)

            T_pol=T(idx,:);
            T_pol=cartesian_to_polar(T_pol);

            i_H_num(ih,:)= -sum(exp(g*( T_pol(:,1) - distHerdGoal)).* ( diff_HT_IH(idx,:)- delta.*[cos(T_pol(:,2)),sin(T_pol(:,2))] )   ,1);

            i_H_den(ih)= sum(exp(g*( T_pol(:,1) - distHerdGoal)),'all');

            i_H(ih,:)=saturation(a,i_H_num(ih,:)./i_H_den(ih),v_H);     % Limits herders maximum speed to v_H

            active_herder_counter=active_herder_counter+1;              % Herder ih is chasing at least one target
            active_herder_idx(ih)=1;

        elseif  distHerdGoal>=rg                                        % If ih is not the closest herder to any target

            f_H(ih,1)=-v_H*(H(ih,1)/distHerdGoal);
            f_H(ih,2)=-v_H*(H(ih,2)/distHerdGoal);

        end
    end

    %%% UPDATES HERDERS AND TARGETS STATES %%%%%

    H_new(:,1)=H(:,1)+(i_H(:,1)+f_H(:,1))*dt;
    H_new(:,2)=H(:,2)+(i_H(:,2)+f_H(:,2))*dt;

    T_new(:,1)=T(:,1)+(i_T(:,1))*dt+randn(M,1)*sqrt(2*D*dt);
    T_new(:,2)=T(:,2)+(i_T(:,2))*dt+randn(M,1)*sqrt(2*D*dt);

    H=H_new;
    T=T_new;

    if mod(i,step_save)==0    %%% If this is a time step to plot/ where to compute the quantities that will be saved

        figure(1)

        O0=fill(R*cos(th),R*sin(th),'y',EdgeColor='none');      % Plots initial region of size R as yellow shaded circle
        alpha(O0,.3)

        hold on
        goal=plot(rg*cos(th),rg*sin(th),color=goal_color,LineWidth=5); % Plots the goal region of size rg as a blue circle
        tar=scatter(T(:,1),T(:,2),Marker="o",SizeData=15, LineWidth=.05,MarkerFaceColor="magenta",MarkerEdgeColor="k",DisplayName="M Targets");
        her=scatter(H(:,1),H(:,2),Marker="diamond",SizeData=22, LineWidth=2.2,MarkerFaceColor="blue",MarkerEdgeColor="blue",DisplayName="N Herders");
        title( ["Fraction of captured targets $\chi=$"+ sprintf(" %.2f",chi(T,rg)) ; "$t=$"+sprintf("%.2f",i*dt)], FontSize=18,Interpreter="latex" )


        %%%% Uncomment if you want the sensing regions plotted around each herder %%%%

        %         for herd=1:N
        %             sens=fill(H(herd,1)+xi*cos(th),H(herd,2)+xi*sin(th),"magenta",EdgeColor='none');
        %             alpha(sens,.15);
        %         end

        hold off

        axis square
        xlim([-1.15*R,1.15*R])
        ylim([-1.15*R,1.15*R])
        xline(-1.15*R,Color=[.9,.9,.9],LineWidth=.000001)
        yline(-1.15*R,Color=[.9,.9,.9],LineWidth=.000001)
        xticks([])
        yticks([])
        ax=gca;
        xticks([])
        yticks([])
        ax.XColor='w';
        ax.YColor='w';
        set(gcf,'color','w');
        axis square

        %%%%TT UNCOMMENT TO SAVE VIDEO %%%%%%%%%%%
        %         frame=getframe(gcf);
        %         writeVideo(v,frame);

        save_counter=save_counter+1;

        % COMPUTES AT THE CURRENT TIME THE OUTPUT QUANTITIES THAT WILL BE SAVED
        chi_r(save_counter)=chi(T,rg);
        chi_rDr(save_counter)=chi(T,rgT);
        active_herders(save_counter)=active_herder_counter;
        rho=vecnorm(H(active_herder_idx,:),2,2);
        avg_rhoH(save_counter)=mean(rho);
        std_rhoH(save_counter)=std(rho);
        time(save_counter)=i*dt;

        %%% CHECK ON EXIT CONDITIONS 
        if max(rho)<rg               %%% Herders' full success: all herders arrive in the goal region (some targets may or may not be lost)

            success_indicator=1;

            H_final=H;
            T_final=T;

            %%% Since an early exit condition is met, cut empty elemets of
            %%% the saved output arrays, and saves the arrays

            time=time(1:save_counter);
            active_herders=active_herders(1:save_counter);
            avg_rhoH=avg_rhoH(1:save_counter);
            std_rhoH=std_rhoH(1:save_counter);
            chi_r=chi_r(1:save_counter);
            chi_rDr=chi_rDr(1:save_counter);

            cputime_elapsed=cputime- t_start;
            time_elapsed=i*dt;
            fprintf("All herdders in Og, simulation ended at t=%.1f. cputime=%.1f\n", i*dt,cputime_elapsed)

            save(filename,"time","active_herders","avg_rhoH","std_rhoH","chi_r","chi_rDr","H_final","T_final","time_elapsed","N","M","R","dt","v_H","a","g","delta","lambda","beta","D","t_save","rg","drg","dR","cputime_elapsed","success_indicator")

            return


        end

        if min(rho)>R+dR            %%% Herders' failure: herders cannot balance the diffusion of the targets and fail to contain them in a circular region smaller than the initial region 
            success_indicator=-1;

            H_final=H;
            T_final=T;

            %%% Since an early exit condition is met, cut empty elemets of
            %%% the saved output arrays, and saves the arrays

            time=time(1:save_counter);
            active_herders=active_herders(1:save_counter);
            avg_rhoH=avg_rhoH(1:save_counter);
            std_rhoH=std_rhoH(1:save_counter);
            chi_r=chi_r(1:save_counter);
            chi_rDr=chi_rDr(1:save_counter);

            cputime_elapsed=cputime- t_start;
            time_elapsed=i*dt;

            fprintf("Herders failed (too far), simulation ended at t=%.1f. Cpu time =%.1f\n", i*dt,cputime_elapsed)

            save(filename,"time","active_herders","avg_rhoH","std_rhoH","chi_r","chi_rDr","H_final","T_final","time_elapsed","N","M","R","dt","v_H","a","g","delta","lambda","beta","D","t_save","rg","drg","dR","cputime_elapsed","success_indicator")

            return


        end
    end


end

%%%% END OF THE SIMULATION. SAVES OUTPUT ARRAYS
time_elapsed=i*dt;
cputime_elapsed=cputime- t_start;
H_final=H;
T_final=T;

fprintf("Simulation ended. cputime=%.1f\n",cputime_elapsed)


save(filename,"time","active_herders","avg_rhoH","std_rhoH","chi_r","chi_rDr","H_final","T_final","time_elapsed","N","M","R","dt","v_H","a","g","delta","lambda","beta","D","t_save","rg","drg","dR","cputime_elapsed","success_indicator")


end