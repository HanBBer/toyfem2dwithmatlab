% This demo analyzes the convergence of NS example

process = 4;
%% Produce data
if process == 1
    errdata = cell(6, 1);
    
    dt = 1e-3;
    N = 2;
    for i = 1:3
        tic;
        N = N*2;
        errdata{i} = CalNSp(dt, N, 10);
        fprintf('%d/6 costs %.4f s\n', i, toc)
    end
    
    N = 4;
    dt = 0.2;
    for i = 4:6
        tic
        dt = dt/2;
        errdata{i} = CalNSp(dt, N, 10);
        fprintf('%d/6 costs %.4f s\n', i, toc)
    end
    
    save errdata.mat errdata
end

%% Visualization  for part.1.
if process == 2
    load errdata.mat
    T = (0.1:0.1:10)';
    type = {'k-','b-.','r.'};
    
    u_errH1 = zeros(100,1);
    u_errL2 = zeros(100,1);
    p_L2 = zeros(100,1);
    figure(1)
    hold on
    for i = 1:3
        data = errdata{i};
        Nu = data.Nu;
        Np = data.Np;
        Mu = data.Mu;
        Ku = data.Ku;
        Mp = data.Mp;
        err = data.err;
        
        uH1 = [Mu+Ku, sparse(Nu, Nu);sparse(Nu, Nu), Mu+Ku];
        uL2 = [Mu, sparse(Nu, Nu);sparse(Nu, Nu), Mu];
        for j = 1:100
            u_errH1(j) = sqrt(err{j}(1:2*Nu)'*uH1*err{j}(1:2*Nu));
            u_errL2(j) = sqrt(err{j}(1:2*Nu)'*uL2*err{j}(1:2*Nu));
            p_L2(j) = sqrt(err{j}(2*Nu+1:end)'*Mp*err{j}(2*Nu+1:end));
        end
        %plot(T, u_errH1, type{i})
        %plot(T, u_errL2, type{i})
        plot(T, p_L2, type{i})
    end
    legend({'h = 0.25', 'h = 0.125', 'h = 0.0625'}, 'location', 'northeast');
    hold off
    
    figure(2)
    hold on
    for i = 4:6
        data = errdata{i};
        Nu = data.Nu;
        Np = data.Np;
        Mu = data.Mu;
        Ku = data.Ku;
        Mp = data.Mp;
        err = data.err;
        
        uH1 = [Mu+Ku, sparse(Nu, Nu);sparse(Nu, Nu), Mu+Ku];
        uL2 = [Mu, sparse(Nu, Nu);sparse(Nu, Nu), Mu];
        for j = 1:100
            u_errH1(j) = sqrt(err{j}(1:2*Nu)'*uH1*err{j}(1:2*Nu));
            u_errL2(j) = sqrt(err{j}(1:2*Nu)'*uL2*err{j}(1:2*Nu));
            p_L2(j) = sqrt(err{j}(2*Nu+1:end)'*Mp*err{j}(2*Nu+1:end));
        end
        %plot(T, u_errH1, type{i-3})
        %plot(T, u_errL2, type{i-3})
        plot(T, p_L2, type{i-3})
    end
    legend({'\Deltat = 0.1', '\Deltat = 0.05', '\Deltat = 0.025'}, 'location', 'northeast');
    hold off
end

%% ALE part
if process == 3
    errdata2 = cell(3, 1);
    
    dt = 0.1;
    N = 2;
    for i = 1:3
        tic;
        N = N*2;
        errdata2{i} = CalNSALEp(dt, N, 10);
        fprintf('%d/3 costs %.4f s\n', i, toc)
    end
    
    save errdata2.mat errdata2
end
%% Visualization for part.3.
if process == 4
    load errdata2.mat
    T = (0.1:0.1:10)';
    type = {'k-','b-.','r.'};
    
    u_errH1 = zeros(100,1);
    p_L2 = zeros(100,1);
    figure(1)
    hold on
    figure(2)
    hold on
    for i = 1:3
        data = errdata2{i};
        Nu = data.Nu;
        Np = data.Np;
        
        err = data.err; 
        
        for j = 1:100
            Mu = data.Mu{j};
            Ku = data.Ku{j};
            Mp = data.Mp{j};
            uH1 = [Mu+Ku, sparse(Nu, Nu);sparse(Nu, Nu), Mu+Ku];
            u_errH1(j) = sqrt(err{j}(1:2*Nu)'*uH1*err{j}(1:2*Nu));
            p_L2(j) = sqrt(err{j}(2*Nu+1:end)'*Mp*err{j}(2*Nu+1:end));
        end
        figure(1)
        plot(T, u_errH1, type{i})
        figure(2)
        plot(T, p_L2, type{i})
    end
    legend({'h = 0.25', 'h = 0.125', 'h = 0.0625'}, 'location', 'northeast');
    hold off
    figure(1)
    legend({'h = 0.25', 'h = 0.125', 'h = 0.0625'}, 'location', 'northeast');
    hold off
end

