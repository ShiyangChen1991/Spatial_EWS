function [ V_c, D_c ] = analytical( model, bifurcation_parameter, spatial_parameter, steady_state, Diffusion)
% Calculate the analytical solution of the covariance matrix using the
% linearized Fokker-Planck equation
% Parameters:
% model: model_type
% steady_state: steady state solution of the model
% diffusion: diffusion matrix for the Fokker-Planck equation

    % model factory
    if model == "harvest"
        fun = @harvest;
    elseif model == "eutrophication"
        fun = @eutrophication;
    elseif model == "veg_turb"
        fun = @veg_turb;
    end
    c = bifurcation_parameter;
    r = spatial_parameter;
    
    % reshape the steady state solution and calculate the Jocabian matrix
    shape = size(steady_state);
    size_x =shape(1); size_y = shape(2);
    x0 = reshape(steady_state,[size_x*size_y,1]);

    % calculate the Jacobian matrix
    options = optimoptions('fsolve','Display','iter');
    [~,~,~,~,jacobian]  = fsolve(fun,x0,options);

    % calculate the eigenvalues and eigenvectors of the matrix
    [V,D] = eig(jacobian);
    
    % coordinate transformation
    V_tilde = inv(V) * Diffusion * inv(V');

    % calculate the covariance matrix
    for i = 1:size_x*size_y
        for j = 1:size_x*size_y
            Sigma_tilde(i,j) = -V_tilde(i,j)*(1+1/(D(i,i)+D(j,j))*(D(i,i)-D(j,j)))/2/D(i,i);
        end
    end
    Cov_Matrix = V*Sigma_tilde*V';

    % calculate the eigenvalues and eigenvectors of the covariance matrix
    [V_c,D_c] = eig(Cov_Matrix);
    
    function dpop = harvest(pop)
    % harvesting model
        K = 10;
        R = 0.2;

        % reshape
        data = zeros(size_x, size_y);
        for k = 1:size_x
            data(k,:) = pop((k-1)*size_x+1:k*size_x);
        end

        % define dpop with reflective boundary condition
        dpop = zeros(size_x*size_y,1);
        for h=1:size_x
            for k=1:size_y
                X = data(h,k);
                dpop((h-1)*size_x+k) = r(h,k)*X*(1-X/K)-c*X^2/(X^2+1)+R*diffusion(size_x, data, h, k);    
            end
        end    
    end

    function dpop = eutrophication(pop)
    % Eutrophication model
        a = 0.5;
        d_r = 0.2;

        % reshape
        data = zeros(size_x, size_y);
        for k = 1:size_x
            data(k,:) = pop((k-1)*size_x+1:k*size_x);
        end

        % define dpop with reflective boundary condition
        dpop = zeros(size_x*size_y,1);
        for h=1:size_x
            for k=1:size_y
                X = data(h,k);
                dpop((h-1)*size_x+k) = a - r(h,k)*X+c*X^8/(X^8+1)+d_r*diffusion(size_x, data, h, k);    
            end
        end    
    end
    function dpop = veg_turb(pop)
    % Vegetation?turbidity model
        d_r = 0.2;
        h_v = 0.2;

        % reshape
        data = zeros(size_x, size_y);
        for k = 1:size_x
            data(k,:) = pop((k-1)*size_x+1:k*size_x);
        end

        % define dpop with reflective boundary condition
        dpop = zeros(size_x*size_y,1);
        for h=1:size_x
            for k=1:size_y
                X = data(h,k);
                E = c*h_v/(h_v+X);
                dpop((h-1)*size_x+k) = 0.5*X*(1-X*(r(h,k)^4+E^4)/r(h,k)^4)+d_r*diffusion(size_x, data, h, k);    
            end
        end    
    end
end

function diff = diffusion(size, data, j, k)
    % calculate diffusion from neighboring cells
    X = data(j,k);
    if j+1>size
        X1 = X;
    else 
        X1 = data(j+1,k);
    end

    if j-1<1
        X2 = X;
    else 
        X2 = data(j-1,k);
    end

    if k-1<1
        X3 = X;
    else 
        X3 = data(j,k-1);
    end

    if k+1>size
        X4 = X;
    else 
        X4 = data(j,k+1);
    end
    diff = X1+X2+X3+X4-4*X;
end
