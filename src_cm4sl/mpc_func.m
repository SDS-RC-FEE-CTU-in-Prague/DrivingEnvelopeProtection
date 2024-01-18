function u = mpc_func(input, params, controller_on)

Ts = 0.005; % used sampling time
mu = 1; % friction estimation

% vehicle parameters
lf   = params(1); % GC-front distance
lr   = params(2); % GC-rear distance
m    = params(3); % mass
I    = params(4); % inertia
ha   = params(5); % half of the axle 
p    = params(6); % wheel radius
caf  = params(7); % cornering stiffness front
car  = params(8); % cornering stiffness rear

% input parsing
x1 = [input(1); input(2)];
v = input(3);
wf_now = input(4);
r = input(6:7);
u_prev = input(8:9);

thr_pedal = input(10);
brk_pedal = input(11);

if(controller_on ~= 1 || v == 0) % check v_min and start_time boundary
    u = r;
else

    %% define the system

    a11 = -mu*(caf+car)/m/v;
    a12 = mu*(car*lr-caf*lf)/m/v/v-1;
    b11 = mu*caf/m/v;

    a21 = (mu*car*lr-caf*lf)/I;
    a22 = -mu*(lf^2*caf+car*lr^2)/I/v;
    b21 = mu*caf*lf/I;

    A = [a11, a12;
         a21, a22];
    B = [b11, 0;
         b21, 0];

    Ad = eye(2) + Ts*A;
    Bd = Ts*B;

    % mpc stuff
    nx = 2; % number of states
    nu = 2; % number of inputs
    
    % Weights
    steering_wheel_track = 20;
    omega_track          = 0.008; % 0.005
    
    Qs   = 1e3; % slack weight
    Qenv = 1e4; % envelope weight
    Qenv_ar = 1e6;
    
    Qdu = [50, 0.45]; % 0.1
    
    % MPC data
    R = [steering_wheel_track,omega_track];
    N = 3;
    
    % changable envelope restrictions
    aFmax = 0.4;
    aRmax = 0.4;
    lFLmax = 0.3;
    lFRmax = 0.3;
    deltaFmax = 0.65;
    deltaFmin = -0.65;
    
    % slew restrictions
    delta_slew = deg2rad(120)*Ts; %rad/s
    omega_slew = 1000*Ts; %400
    slew = [delta_slew; omega_slew];

    
    % u11,u12,...,u1N
    % u21,u22,...,u2N
    % x12,...,x1N+1
    % x22,...,x2N+1
    % sAR2,sAR3,...sARN+1 % in step 1 AR cannot be optimized
    % sD1,sD2,...,sDN
    % sW1,sW2,...,sWN
    % sSFL1,sSFL2,...sSFLN+1
    % sSFR1,sSFR2,...sSFRN+1
    % c11,c12,...,c1N
    % c21,c22,...,c2N
    
    abs_num = 3; % count of abs in objective function
    var_number = (nx+nu)*N + 3*N + 2*(N+1) + abs_num*nu; %(states+inputs)*N+slacks+abs_num*inputs
    Q = zeros(var_number);
    c = zeros(var_number, 1);
    A = zeros(2*nu*N+3*N+3*nu*N+2*5*(N+1)+abs_num*nu*2,var_number); % inequality constraints (min/max, aR and its slacks, du and its slacks, sigmaFL and sigmaFR and its slacks abs) 
    b = zeros(2*nu*N+3*N+3*nu*N+2*5*(N+1)+abs_num*nu*2,1);
    Aeq = zeros(nx*N, var_number); %equality constraints (dynamics)
    beq = zeros(nx*N,1);
    lb = zeros(var_number, 1); % boundaries for variables (seems to be unused in yalmip)
    ub = zeros(var_number, 1);
    
    cons_count = 0;
    
    %% min/max inputs 
    for j = 1:N
        for i = 1:nu
            A(2*nu*(j-1)+nu*(i-1)+(1:nu), (j-1)*nu+i) = [-1;1];
            if(i==1)
                b(2*nu*(j-1)+nu*(i-1)+(1:nu)) = [-deltaFmin; deltaFmax];
            end
            if(i==2)
                b(2*nu*(j-1)+nu*(i-1)+(1:nu)) = [-(wf_now-j*omega_slew); wf_now+j*omega_slew];
            end
        end
    end
    cons_count = cons_count + 2*N*nu;
    %% aR envelope and slack constraints
    for j = 2:N+1
        % ar
        A(cons_count+1+(j-2)*3+(0:2),2*N+(j-2)*nx+(1:nx)) = [-1, lr/v;
                                                1, -lr/v;
                                                0, 0];
        A(cons_count+1+(j-2)*3+(0:2), nu*N+nx*N+j-1) = -ones(3,1);
        b(cons_count+1+(j-2)*3+(0:2)) = [aRmax;
                             aRmax;
                             0];
        Q(nu*N+nx*N+j-1, nu*N+nx*N+j-1) = 2*Qenv_ar;
    end
    cons_count = cons_count + 3*N;
    
    %% slew protection
    for j = 1:N
        if j == 1
            A(cons_count+(1:nu),1:nu) = -eye(2);
            A(cons_count+(1:nu),nu*N+nx*N+N+(1:nu)) = -eye(2);
            b(cons_count+(1:nu)) = -u_prev + slew;
            A(cons_count+nu+(1:nu),1:nu) = eye(2);
            A(cons_count+nu+(1:nu),nu*N+nx*N+N+(1:nu)) = -eye(2);
            b(cons_count+nu+(1:nu)) = u_prev+slew;     
            A(cons_count+nu+nu+(1:nu),nu*N+nx*N+N+(1:nu)) = -eye(2);
        else
            A(cons_count+(j-1)*nu*3+(1:nu),(j-1)*nu-nu+(1:nu)) = eye(2);
            A(cons_count+(j-1)*nu*3+(1:nu),(j-1)*nu+(1:nu)) = -eye(2);
            A(cons_count+(j-1)*nu*3+(1:nu),nu*N+nx*N+N+(j-1)*nu+(1:nu)) = -eye(2);
            b(cons_count+(j-1)*nu*3+(1:nu)) = slew;
            A(cons_count+nu+(j-1)*nu*3+(1:nu),(j-1)*nu-nu+(1:nu)) = -eye(2);
            A(cons_count+nu+(j-1)*nu*3+(1:nu),(j-1)*nu+(1:nu)) = +eye(2);
            A(cons_count+nu+(j-1)*nu*3+(1:nu),nu*N+nx*N+N+(j-1)*nu+(1:nu)) = -eye(2);
            b(cons_count+nu+(j-1)*nu*3+(1:nu)) = slew;
            A(cons_count+nu+(j-1)*nu*3+nu+(1:nu),nu*N+nx*N+N+(j-1)*nu+(1:nu)) = -eye(2);                                 
        end
        Q(nu*N+nx*N+N+(j-1)*nu+(1:nu),nu*N+nx*N+N+(j-1)*nu+(1:nu)) = 2*Qs*eye(nu);
    end
    cons_count = cons_count + 3*nu*N;
    
    %% sigmaFL envelope and slack constrains
%     mu = 1;
    X = (1-lFLmax)/lFLmax;
    Y = (1-lFLmax)/tan(aFmax);
    inputsFL = [-Y, -X*p/v;
                 Y,  X*p/v;
                 Y, -X*p/v;
                -Y,  X*p/v;];
    statesFL = [ Y, (X*ha+Y*lf)/v;
                -Y, (-X*ha-Y*lf)/v;
                -Y, (X*ha-Y*lf)/v;
                 Y, (-X*ha+Y*lf)/v;];
             
    for i = 1:N
        if i == 1
            A(cons_count+(1:5),(i-1)*nu+(1:nu)) = [inputsFL; 0, 0];
            A(cons_count+(1:5),nu*N+nx*N+3*N+i) = -ones(5,1);
            b(cons_count+(1:5)) = mu*[ones(4,1);0] + 1*[-X; X; -X; X; 0] - ...
                                  [statesFL; 0, 0] * x1;
        else
            A(cons_count+(1:5),(i-1)*nu+(1:nu)) = [inputsFL; 0, 0];
            A(cons_count+(1:5),nu*N+(i-2)*nx+(1:nx)) = [statesFL; 0, 0];
            A(cons_count+(1:5),nu*N+nx*N+3*N+i) = -ones(5,1);
            b(cons_count+(1:5)) = mu*[ones(4,1);0] + 1*[-X; X; -X; X; 0];
        end
        cons_count = cons_count + 5;
        Q(nu*N+nx*N+3*N+i,nu*N+nx*N+3*N+i) = 2*Qenv;
    end
    % final state
            A(cons_count+(1:5),(N-1)*nu+(1:nu)) = [inputsFL; 0, 0];
            A(cons_count+(1:5),nu*N+(N-1)*nx+(1:nx)) = [statesFL; 0, 0];
            A(cons_count+(1:5),nu*N+nx*N+3*N+N+1) = -ones(5,1);
            b(cons_count+(1:5)) = mu*[ones(4,1);0] + 1*[-X; X; -X; X; 0];   
            cons_count = cons_count + 5;
            Q(nu*N+nx*N+3*N+i,nu*N+nx*N+3*N+i) = 2*Qenv;
    
    %% sigmaFR envelope and slack constrains
%     mu = .4;
    X = (1-lFRmax)/lFRmax;
    Y = (1-lFRmax)/tan(aFmax);
    inputsFR = [-Y, -X*p/v;
                 Y,  X*p/v;
                 Y, -X*p/v;
                -Y,  X*p/v;];
    statesFR = [ Y, (-X*ha+Y*lf)/v;
                -Y, (X*ha-Y*lf)/v;
                -Y, (-X*ha-Y*lf)/v;
                 Y, (X*ha+Y*lf)/v;];
    for i = 1:N
        if i == 1
            A(cons_count+(1:5),(i-1)*nu+(1:nu)) = [inputsFR; 0, 0];
            A(cons_count+(1:5),nu*N+nx*N+4*N+1+i) = -ones(5,1);
            b(cons_count+(1:5)) = mu*[ones(4,1);0] + 1*[-X; X; -X; X; 0] - ...
                                  [statesFR; 0, 0] * x1;
        else
            A(cons_count+(1:5),(i-1)*nu+(1:nu)) = [inputsFR; 0, 0];
            A(cons_count+(1:5),nu*N+(i-2)*nx+(1:nx)) = [statesFR; 0, 0];
            A(cons_count+(1:5),nu*N+nx*N+4*N+1+i) = -ones(5,1);
            b(cons_count+(1:5)) = mu*[ones(4,1);0] + 1*[-X; X; -X; X; 0];
        end
        cons_count = cons_count + 5;
        Q(nu*N+nx*N+4*N+i,nu*N+nx*N+4*N+i) = 2*Qenv;
    end
    % final state
            A(cons_count+(1:5),(N-1)*nu+(1:nu)) = [inputsFR; 0, 0];
            A(cons_count+(1:5),nu*N+(N-1)*nx+(1:nx)) = [statesFR; 0, 0];
            A(cons_count+(1:5),nu*N+nx*N+4*N+N+1+1) = -ones(5,1);
            b(cons_count+(1:5)) = mu*[ones(4,1);0] + 1*[-X; X; -X; X; 0];   
            cons_count = cons_count + 5;
            Q(nu*N+nx*N+4*N+1+i,nu*N+nx*N+4*N+1+i) = 2*Qenv;
    
    
    %% tracking objective (for the first steps only)
%         obj = obj + R*abs(r - u{1}) + R*(r - u{1}).^2;
%         obj = obj + R*abs(r - u{2}) + R*(r - u{2}).^2;
    
    Q(1:nu*abs_num,1:nu*abs_num) = 2*diag(repmat(R,1,abs_num));
    
    c(1:nu*abs_num) = repmat(-2*R.*r',1,abs_num);
    c(end-(nu*abs_num-1):end) = repmat(R,1,abs_num);
    
    lb(:) = -Inf;
    ub(:) = Inf;
    
    % abs in obj 
    for j = 1:abs_num
        for i = 1:nu
            A(cons_count+2*nu*(j-1)+nu*(i-1)+(1:2), (j-1)*nu+i) = [-1;1];
            A(cons_count+2*nu*(j-1)+nu*(i-1)+(1:2), (nx+nu)*N + 3*N + 2*(N+1)+(j-1)*nu+i) = -ones(nu,1);
            b(cons_count+2*nu*(j-1)+nu*(i-1)+(1:2)) = [-r(i);r(i)];
        end
    end
    
    %% general dynamics
    Aeq(1:nx, 1:nu) = -Bd;
    Aeq(1:nx, nu*N+(1:nx)) = eye(nx);
    beq(1:nx) = Ad*x1;
    for k = 2:N
        Aeq((k-1)*nx+(1:nx), (k-1)*nu+(1:nu)) = -Bd;
        Aeq((k-1)*nx+(1:nx), (k-2)*nx+nu*N+(1:nx)) = -Ad;
        Aeq((k-1)*nx+(1:nx), nu*N+(k-1)*nx+(1:nx)) = eye(nx);
    end

    %% du objective
    c(1:nu) = c(1:nu) - (2*Qdu.*u_prev')';
    Q_tmp = diag(repmat(Qdu*4,1,N))+diag(repmat(-2*Qdu,1,N-1),2)+diag(repmat(-2*Qdu,1,N-1),-2);
    for i = 0:nu-1
        Q_tmp(end-i,end-i) = Q_tmp(end-i,end-i)/2;
    end
    Q(1:nu*N,1:nu*N) = Q(1:nu*N,1:nu*N)+Q_tmp;
    
    %%
    Qs = sparse(Q);  
    lbA = [beq; -inf(length(b),1)];
    ubA = [beq; b];
    Asolver = sparse([Aeq; A]);
    
    options = qpOASES_options('default', 'printLevel', 3);
    [x,~,exitflag,~,~] = qpOASES(Qs, c, Asolver, lb, ub, lbA, ubA, {options, qpOASES_auxInput('x0',r)});

    %% obtain a solution
    if (exitflag~=0)
        u = u_prev;
    else
        u = x(1:2);
    end
end
end
   