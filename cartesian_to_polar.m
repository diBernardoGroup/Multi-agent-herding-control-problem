function X_polar=cartesian_to_polar(X)

% Converts an array of 2d cartesian coordinates to an array of 2d polar
% coordinates

    X_polar=X;

    for k=1:size(X,1)
        
        X_polar(k,1)=((X(k,1)^2+ X(k,2)^2 )^.5);

        if X(k,2)<0
            X_polar(k,2)=2*pi-acos(X(k,1)/X_polar(k,1));
        else 
            X_polar(k,2)=acos(X(k,1)/X_polar(k,1));
        end

    end


end

    
