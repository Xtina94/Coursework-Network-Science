%% %%%%%%%%%%%%%%%%% Ranking %%%%%%%%%%%%%%%%%%%%
close all
clear all

y = 1;
switch y
    case 1
        G = importdata('wiki-Vote.txt', '\t', 4);
        N_red = max(max(G.data));
        A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N_red,N_red);
        Au = 1*(A+A'>0); % undirected network
        W = A;
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
    case 2
        G = importdata('bioCel.txt');
        N_red = max(max(G));
        A = sparse(G(:,2),G(:,1),ones(size(G,1),1),N_red,N_red);
        Au = 1*(A+A'>0); % undirected network
        W = A;
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
    case 3
        G = importdata('ca_sandi_auths.txt');
        N_red = max(max(G));
        A = sparse(G(:,2),G(:,1),ones(size(G,1),1),N_red,N_red);
        Au = 1*(A+A'>0); % undirected network
        W = A;
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
    case 4
        G = importdata('collaboration.edgelist.txt', '\t', 1);
        G.data = G.data + 1;
        N_red = max(max(G.data));
        A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N_red,N_red);
        Au = 1*(A+A'>0); % undirected network
        W = A;
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
end

% First, remove the nodes that are isolated
pos = find(sum(A)~=0);
A_red = A(pos,pos);
N_red = size(A_red,1);

% Create the equivalent undirected network, only if the network is directed
if directed
    pos = find(sum(Au)~=0);
    Au = Au(pos,pos);
    A_red = A(pos,pos);
    N_red = size(A_red,1);
    
    if y == 1
        % Dead ends removal
        exit = false;
        while (~exit)
            pos = find(sum(A_red)~=0);
            A_red = A_red(pos,pos);
            Au = Au(pos,pos);
            N_red = size(A_red,1);
            exit = isempty(find(sum(A_red)==0, 1));
        end
    end
else
    Au = A_red;
end


if directed % find the largest connected component for the directed graph
    e1 = [1;zeros(N_red-1,1)];
    exit = false;
    while(~exit)
        e1_old = e1;
        e1 = 1*(Au*e1>0);
        exit = (sum(e1-e1_old)==0);
    end
    idx = find(e1);
    A_red = A_red(idx,idx);
    N_red = size(A_red,1);
% else % find the largest connected component for the undirected graph
%     e1 = [1;zeros(N_red-1,1)];
%     exit = false;
%     while(~exit)
%         e1_old = e1;
%         e1 = 1*(A_red*e1>0);
%         exit = (sum(e1-e1_old)==0);
%     end
%     idx = find(e1);
%     A_red = A_red(idx,idx);
%     N_red = size(A_red,1);
end

%% PageRank equation evaluation 

c = 0.85; %damping factor
q = ones(N_red,1)/N_red; %teleportation vector
d_out = (A_red)'*ones(size(A_red,1),1); %output degree vector
d_out_2 = sum(A_red,2);
diff = d_out - d_out_2;
d_out = sparse(ones(length(d_out),1)./d_out);
d_out_2 = sparse(ones(length(d_out_2),1)./d_out_2);
AdjM = A_red*sparse(diag(d_out)); %Weighted adjacency matrix
M_11 = c*AdjM;
M_12 = (1-c)*q;

% Through linear system solution
tic
pr_ls = (sparse(eye(size(AdjM,1))) - M_11)\M_12;
pr_ls = pr_ls/sum(pr_ls);
toc

% Through power iteration
t = 25; % number of iterations
tic
pr = ones(N_red,1)/N_red; % initial guess on p_0(i)
s = zeros(t,1);
for i = 1:t
    pr_old = pr;
    % iterative step
    pr = M_11*pr + M_12;
    pr = pr/sum(pr);
    s(i) = norm(pr - pr_ls)/sqrt(N_red);
end
toc
distance = s;

%%%%%%%%%%%% Checking for convergence %%%%%%%%%%%%%

% Extract eigenvalues
disp('6 biggest eigenvalues extraction')
tic
% dbstop if naninf
[V, lambdas Flag]= eigs(AdjM);
toc
if Flag
    disp('Eigenvalues not found, eigs does not converge!')
end

lambda_2 = lambdas(2,2); %take the second eigenvalue, since the first one is 1 and stands for the connected component

thr = (c*abs(lambda_2)).^(1:t);
convergence = (distance' <= thr); %Checking if there is convergence in the solution

figure('Name','Convergence of the power iteration method for PageRank equation')
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
semilogy((1:t),distance)
grid
xlabel('k [iteration \#]')
ylabel('$\|r_k - r_\infty\|$')
title('PageRank convergence')

lambdas = eigs(AdjM,size(AdjM,1));

figure('Name','Eigenvalues representation')
plot(real(lambdas(2:end)),imag(lambdas(2:end)),'o')
hold on
viscircles([0 0],lambdas(1),'Color','g');
hold off
grid

%% HITS equation evaluation - authorities scores

AdjM = sparse(A_red*A_red');

% Through linear system solution
tic
[hits_ls_vect, hits_ls] = eigs(AdjM,2);
toc
pr_dir_aut = -hits_ls_vect(:,1)/norm(hits_ls_vect(:,1)); % Extract the most importand dimension and project the values wrt to that (SVD)

% Through power iteration
tic
hits = ones(N_red,1)/N_red; % initial guess
s = zeros(t,1);
for i = 1:t
    temp = hits;
    hits = AdjM*hits; % iterative step
    hits = hits/norm(hits); % normalization step
    s(i) = norm(hits - temp)/sqrt(N_red);
end
toc

thr_hubs = (abs(hits_ls(2,2)/hits_ls(1,1))).^(1:t);
norm_thr = thr_hubs/thr_hubs(end); % normalized threshold
convergence_hits = (s <= thr); %Checking if there is convergence in the solution

figure('Name','Convergence of the power iteration method for Hits authorities')
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
semilogy((1:t),s)
hold on
semilogy((1:t),s(end)*norm_thr);
hold off
grid
xlabel('k [iteration \#]')
ylabel('$\|hits_k - hits_\infty\|$')
title('Hits convergence')

%% HITS equation evaluation - hubs scores

% Through linear system solution
tic
M_2 = c*AdjM;
[hits_ls_vect, hits_ls] = eigs(M_2,2);
pr_dir_hub = -hits_ls_vect(:,1)/norm(hits_ls_vect(:,1)); % Extract the most importand dimension and project the values wrt to that (SVD)
toc

% Through power iteration
tic
hits = ones(N_red,1)/N_red; % initial guess
s = zeros(t,1);
for i = 1:t
    temp = hits;
    hits = M_2*hits; % iterative step
    hits = hits/norm(hits); % normalization step
    s(i) = norm(hits - temp)/sqrt(N_red);
end
toc

thr_hubs = (abs(hits_ls(2,2)/hits_ls(1,1))).^(1:t);
norm_thr = thr_hubs/thr_hubs(end);
convergence_hits = (s <= thr); %Checking if there is convergence in the solution

figure('Name','Convergence of the power iteration method for Hits hubs')
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
semilogy((1:t),s)
hold on
semilogy((1:t),s(end)*norm_thr);
hold off
grid
xlabel('k [iteration \#]')
ylabel('$\|hits_k - hits_\infty\|$')
title('Hits convergence')

%% %%%%%%%%%%% Comparing PageRank and Hits %%%%%%%%%%%%%%

figure('Name','Comparison between PageRank and Hits authorities scores')
plot(pr_dir_aut/sum(pr_dir_aut),pr_ls/sum(pr_ls),'o')
grid
title('PageRank vs HITS ')

figure('Name','Comparison between PageRank and hits authorities scores')
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
plot(1:N_red,pr_ls/sum(pr_ls))
grid
hold on
plot(1:N_red,-pr_dir_aut/sum(pr_dir_aut))
hold off
legend('PageRank','Hits')
title('PageRank vs HITS ')

% % % figure('Name','Comparison between PageRank and Hits hubs scores')
% % % set(0,'defaultTextInterpreter','latex') % to use LaTeX format
% % % % plot(hits,pr,'o')
% % % plot(pr_dir_hub/sum(pr_dir_hub),pr_ls/sum(pr_ls),'o')
% % % grid
% % % xlabel('Hits scores')
% % % ylabel('PageRank scores')
% % % title('Hits convergence')
% % % 
% % % figure('Name','Comparison between PageRank and hits hubs scores')
% % % set(0,'defaultTextInterpreter','latex') % to use LaTeX format
% % % plot(1:N,pr_ls/sum(pr_ls))
% % % grid
% % % hold on
% % % plot(1:N,-pr_dir_hub/sum(pr_dir_hub))
% % % hold off
% % % legend('PageRank','Hits')
% % % title('PageRank vs HITS ')
