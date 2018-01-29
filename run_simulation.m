function [ S, T ] = run_simulation( model, size, bifurcation_parameter, sigma, spatial_parameter, nPeriods, dt)
    % This function runs the simulation of a spatially distributed ecological system.  
    % The system is modelled by a stochastic differential equation
    % Parameters
    % model: model_type
    % size: size*size is the total number of cells in the model
    % bifurcation_parameter: [c_0, delta_c], c_0 is the initial value of the paprameter and delta_c is the varying rate
    % sigma: magnitude of stochastic noise
    % nPeriods and dt specify the time frame of the simulation

    % bifurcation parameters
    c0 = bifurcation_parameter(1);
    delta_c = bifurcation_parameter(2);

    % spatial parameter
    r = spatial_parameter;

    % initial condition
    x_ini = 10*rand(size*size, 1);

    % stochastic simulation
    % drift
    if model == 'harvest'
        F = @ harvest;
    end
    % diffusion
    G0 = eye(size*size, size*size);
    G = @(t,X) sigma*G0;

    SDE = sde(F,G,'StartState', x_ini);
    [S,T] = simulate(SDE, nPeriods, 'DeltaTime', dt);

        function dpop = harvest(t, pop)
        % harvesting model
            K = 10;
            D = 0.2;
            c = c0 + delta_c*t;

            % reshape
            data = zeros(size, size);
            for j = 1:size
                data(j,:) = pop((j-1)*size+1:j*size);
            end

            % define dpop with reflective boundary condition
            dpop = zeros(size*size,1);
                for j=1:size
                    for k=1:size
                        X = data(j,k);
                        dpop((j-1)*size+k) = r(j,k)*X*(1-X/K)-c*X^2/(X^2+1)+D*diffusion(size, data, j, k);    
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
