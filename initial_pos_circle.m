function pos=initial_pos_circle(N,R)
% Randomly and uniformlu distributes N agents in a circle of radius R around the origin 
% Returns the cartesian 2d coordinates

    pos=zeros(N,2);

    r=zeros(N);         % Radial coordinate
    theta=zeros(N);     % Angular coordinate

    
    for i=1:N
        r(i)=sqrt(R*R*rand());
        theta(i)=rand()*2*pi;
    end
    pos=zeros(N,2);
                                % Polar to caresian transformation
    for i=1:N
        pos(i,1)=r(i)*cos(theta(i));
        pos(i,2)=r(i)*sin(theta(i));
    end



end