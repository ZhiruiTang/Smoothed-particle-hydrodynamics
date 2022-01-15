function SPH3D

%% The purpose of this code
%  Simulate the evolution of N body, according to this paper:
%  https://pmocz.github.io/manuscripts/pmocz_sph.pdf

%%
	% Simulation parameters
	N         = 200;    % Number of particles
	tEnd      = 60;     % time at which simulation ends
	dt        = 0.04;   % timestep
	M         = 2;      % star mass
	R         = 0.75;   % star radius
	h         = 0.1;    % smoothing length
	k         = 0.1;    % equation of state constant
	n         = 1;      % polytropic index
	nu        = 1;      % damping
	
    
    %initialization
    pos     = randn(N, 3)*2;
    vel     = zeros(size(pos));
    v_mhalf = zeros(size(pos));
    m       = M/n; %single particle mass
    Nt      = ceil(tEnd/dt); %number of timestep
    lambda  = 2*k*(1+n)*pi^(-3/(2*n)) * (M*gamma(5/2+n)/R^3/gamma(1+n))^(1/n)/R^2;
    
    %initialize accelerations
    rho = Cal_Density(pos, m, h);
    P   = k * rho.^(1+1/n);
    acc = Cal_Acc(pos, vel, m, rho, P, nu, lambda, h);
    
    %For a better quality of animation, we store all the position and
    %density info. to make animation later.
    POS = zeros(N, 3, Nt+1);
    RHO = zeros(N, Nt+1);
    
    %Store the initial values
    POS(:, :, 1) = pos;
    RHO(:, 1)    = rho;
    
    %Draw images on time

    ax = plot3(pos(:,1), pos(:,2), pos(:,3), 'MarkerFaceColor','#D9FFFF');
    ax.Marker = 'o';
    ax.Color = 'b';
    ax.LineStyle = 'none';
    axis([-30 30 -30 30 -30 30]);
    hold on 
    
    %Main loop: the leap-frog scheme
    for i = 1 : Nt
        v_phalf = v_mhalf + acc * dt;
        pos     = pos + v_phalf * dt;
        vel     = 0.5*(v_mhalf + v_phalf);
        v_mhalf = v_phalf;
        
        %update densities, pressures, accelerations
        rho = Cal_Density(pos, m, h);
        P   = k * rho.^(1+1/n);
        acc = Cal_Acc(pos, vel, m, rho, P, nu, lambda, h);
        
        %store info.
        POS(:, :, i+1) = pos;
        RHO(:, i+1) = rho;
        
        %scatter3(pos(:,1), pos(:,2), pos(:,3), [], rho)
         ax.XData = pos(:,1);
         ax.YData = pos(:,2);
         ax.ZData = pos(:,3);
        pause(0.01)
        
    end
    
    %Save data
    save('pos_info.mat', 'pos');
    save('rho_info.mat', 'rho')

end

%% Gaussian kernel and its gradient
function w= Kernel(rij, h)
    w = 1.0/(h*pi^(3/2)) * exp(-norm(rij)^2/h^2);
end

function dw = Grad_Kernel(rij, h)
    dw = -2.0/(h^5*pi^(3/2)) * exp(-norm(rij)^2/h^2) * rij;
end

%% Calculate density 
function rho = Cal_Density(pos, m, h)
    N = size(pos, 1);
    rho = zeros(N, 1);
    for i = 1 : N
        % initialize density with i = j contribution 
        rho(i) = m * Kernel(0, h);
        for j = i + 1 : N
        % calculate vector between two particles
            uij = pos(i, :) - pos(j, :);
            rho_ij = m * Kernel(uij, h);
            % add contribution to density
            rho(i) = rho(i) + rho_ij;
            rho(j) = rho(j) + rho_ij; 
        end
    end
end

%% Calculate acceleration
function a = Cal_Acc(pos, vel, m, rho, P, nu, lambda, h)
% initialize accelerations
    N = size(pos,1);
    a = zeros(N, 3); % add damping and gravity

    a = a - nu * vel + - lambda*pos; % repmat([0 0 -9.8*m], num_p, 1) ;

    % add pressure
    for i =1: N
        for j = i+1 : N
            % calculate vector between two particles
            uij = pos(i, :) - pos(j, :);

            % calculate acceleration due to pressure
            p_a = -m * (P(i)/rho(i)^2 + P(j)/rho(j)^2) * Grad_Kernel(uij, h);
            a(i, :) = a(i, :) + p_a;
            a(j, :) = a(j, :) - p_a;
        end
    end
end



