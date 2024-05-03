function fraction=chi(T,rg)
% T is an Nx2 array containing the cartesian coordinates of N agents
% Returns the fraction of the N agents in the plane that are within a circle of
% radius rg centered aroun the origin

    captured_T=0;
        
        for i=1:size(T,1)
            T_pol=cartesian_to_polar(T(i,:));
            if T_pol(1)<=rg
                captured_T=captured_T+1;
            end
        end

    fraction=captured_T/size(T,1);

end