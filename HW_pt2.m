%% %%%%%%%%%%%%%%%%%% ASSORTATIVITY ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
addpath('C:\Users\cryga\Documents\GitHub\HomeworkNS\Datasets');

y = 5;
switch y
    case 1
        G = importdata('wiki-Vote.txt', '\t', 4);
        N = max(max(G.data));
        A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
        W = A; %Save adjacency matrix for later work on communities
        clear G;
        directed = 1;
    case 2
        G = importdata('bioCel.txt');
        N = max(max(G));
        A = sparse(G(:,2),G(:,1),ones(size(G,1),1),N,N);
        W = A;
        clear G;
        directed = 0;
    case 3
        G = importdata('ca_sandi_auths.txt');
        N = max(max(G));
        A = sparse(G(:,2),G(:,1),ones(size(G,1),1),N,N);
        W = A;
        clear G;
        directed = 0;
    case 4
        G = importdata('collaboration.edgelist.txt', '\t', 1);
        G.data = G.data + 1;
        N = max(max(G.data));
        A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
        A = 1*(A+A'>0); % build undirected network
        clear G;
        directed = 1;
    case 5
        A = importdata('occupyWs.txt');
        N = max(max(A));
        A = sparse(A(:,2),A(:,1),ones(size(A,1),1),N,N);
        A = 1*(A+A'>0); % build undirected network
        clear G;
        directed = 1;
end

%% Avg neigh degree k_nn (it is given by k_nn = mean(q_m))

% First, remove the nodes that are isolated
pos = find(sum(A)~=0);
A_red = A(pos,pos);
N_red = size(A_red,1);

% Then find the largest connected component in the graph
%%%%%%%%%%%%%%%% Riattivare una volta capito qual è l'errore %%%%%%%%%%%%%%
if directed
    e1 = [1;zeros(N_red-1,1)];
    exit = false;
    while(~exit)
        e1_old = e1;
        e1 = 1*(A_red*e1+e1>0);
        exit = (sum(e1-e1_old)==0);
    end
    idx = find(e1);
    A_red = A_red(idx,idx);
    N_red = size(A_red,1);
end

% Estimation of avg k_nn = mean on degree

% degrees
d_in = sum(A_red,2);
d_out = sum(A_red)';

% averages of neighbours
k_tmp_oo = (A_red'*d_out)./d_out;
k_tmp_oi = (A_red'*d_in)./d_out;
k_tmp_ii = (A_red*d_in)./d_in;
k_tmp_io = (A_red*d_out)./d_in;

% extract averages for each value of k
u_out = unique(d_out);
for i = 1:length(u_out)
    k_nn_oo(i) = mean(k_tmp_oo(d_out==u_out(i)));
    k_nn_oi(i) = mean(k_tmp_oi(d_out==u_out(i))); % this is the mean of the neighbouring degree 
                                                      % of all the nodes having the same degree
end

u_in = unique(d_in);
for i = 1:length(u_in)
    k_nn_io(i) = mean(k_tmp_io(d_in == u_in(i)));
    k_nn_ii(i) = mean(k_tmp_ii(d_in == u_in(i)));
end

% do the linear fittings
p_oo = polyfit(log(u_out(2:end)'),log(k_nn_oo(2:end)),1);
disp(['Assortativity factor out-out ' num2str(p_oo(1))])
p_io = polyfit(log(u_in(2:end)'),log(k_nn_io(2:end)),1);
disp(['Assortativity factor in-out ' num2str(p_io(1))])
p_oi = polyfit(log(u_out(2:end)'),log(k_nn_oi(2:end)),1);
disp(['Assortativity factor out-in ' num2str(p_oi(1))])
p_ii = polyfit(log(u_in(2:end)'),log(k_nn_ii(2:end)),1);
disp(['Assortativity factor in-in ' num2str(p_ii(1))])

%% Plot the results
figure('Name','Avg neighbouring degree')
subplot(2,2,1)
loglog(d_out,k_tmp_oo,'b.');
hold on
loglog(u_out,exp(p_oo(2)+log(u_out)*p_oo(1)),'r-');
loglog(u_out,k_nn_oo,'g.');
grid
xlabel('k_{out}')
ylabel('k_{nn,out}')
title('Network Assortativity. out --> out')

subplot(2,2,2)
loglog(d_out,k_tmp_oi,'b.');
hold on
loglog(u_out,exp(p_oi(2)+log(u_out)*p_oi(1)),'r-');
loglog(u_out,k_nn_oi,'g.');
grid
xlabel('k_{out}')
ylabel('k_{nn,out}')
title('Network Assortativity. out --> in')

subplot(2,2,3)
loglog(d_in,k_tmp_io,'b.');
hold on
loglog(u_in,exp(p_io(2)+log(u_in)*p_io(1)),'r-');
loglog(u_in,k_nn_io,'g.');
grid
xlabel('k_{in}')
ylabel('k_{nn,in}')
title('Network Assortativity. in --> out')

subplot(2,2,4)
loglog(d_in,k_tmp_ii,'b.');
hold on
loglog(u_in,exp(p_ii(2)+log(u_in)*p_ii(1)),'r-');
loglog(u_in,k_nn_ii,'g.');
grid
xlabel('k_{in}')
ylabel('k_{nn,in}')
title('Network Assortativity. in --> in')