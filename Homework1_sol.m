close all
clear all

y = 3;
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
end


% adjacency matrix
% NOTA PER ME: in G.data non sono elencati tutti i nodi che sono isolati,
% per questo motivo il numero massimo di nodi presente è maggiore del
% numero di elementi effettivamente presenti in G.data, perchè i nodi non
% elencati sono isolati



%% Extract the distribution

% distribution. degree is a 2 columns vector with node# --> deg
d = full(sum(A,1));
% Remove the first elements of k since 0 is not an acceptable degree
d = d(d > 0);
k = unique(d); % degree samples
k_min = min(k);
M = length(k);

% Linear PDF plot
pdf = histc(d,k)';
n_links = sum(d); % the total number of links in the network
pdf = pdf/n_links; % normalize to 1

% cumulative distribution
Pdf = cumsum(pdf,'reverse');

% log binning
klog = 10.^(0:0.1:ceil(log10(max(k))));
pklog = histc(d,klog)'; % counts occurrences
pklog = pklog/sum(pklog); % normalize to 1

%% Plot the results

figure('Name','Pdf, logaritmic Pdf and CCDF plots')
subplot(2,2,1)
plot(k,pdf,'.')
grid
xlabel('k')
ylabel('PDF')
title('linear PDF plot')
subplot(2,2,2)
loglog(k,pdf,'.')
grid
xlabel('k')
ylabel('PDF')
title('logarithmic PDF plot')
subplot(2,2,3)
loglog(klog,pklog,'.')
grid
xlabel('k')
ylabel('PDF')
title('logarithmic PDF plot (log bins)')
subplot(2,2,4)
loglog(k,Pdf,'.')
grid
xlabel('k')
ylabel('CCDF')
title('logarithmic CCDF plot')


%% PDF Logarithmic CCDM plot

%With fixed gamma = 3.5
gamma = 3.5;
c = (gamma-1)*k_min^(gamma-1);

%With ML estimation of gamma
k_min = floor(max(d)/2);
d2 = d(d>=k_min); % restrict range

j = d2/k_min;
den = sum(log(j));
n_nodes = sum(length(d2));
gamma_ML = 1+(n_nodes/den);

p_k_35 = k.^(-gamma);
p_k = k.^(-gamma_ML);

ccdm_pdf_35 = k.^(1-gamma);
ccdm_pdf = k.^(1-gamma_ML);

% Plot the results
interval = [ceil(length(k)/3*2):length(k)];
figure('Name','Logarithmic pdf and CCDF plot through ML estimation')
subplot(1,2,1)
thr = datasample(interval,1);
loglog(k,pdf,'o')
hold on
loglog(k,p_k_35/p_k_35(thr)*pdf(thr),'--');
loglog(k,p_k/p_k(thr)*pdf(thr),'--');
hold off
axis([xlim min(pdf/2) 2*max(pdf)])
grid
subplot(1,2,2)
thr = datasample(interval,1);
loglog(k,Pdf,'o')
hold on
loglog(k,ccdm_pdf_35/ccdm_pdf_35(thr)*Pdf(thr),'--');
loglog(k,ccdm_pdf/ccdm_pdf(thr)*Pdf(thr),'--');
hold off
axis([xlim min(pklog/2) 2*max(Pdf)])
grid

%% ML fitting with saturation
% See it later
% % % 
% % % for ks = 1:max(k)
% % %     k_min = min(d);
% % %     tmp = mean(log((d+ks)/(k_min+ks)));
% % %     ga2(ks) = 1+1/tmp;
% % %     de(ks) = log(ga2(ks)-1)-log(k_min+ks)-ga2(ks)*tmp;
% % % end
% % % [~,ks] = max(de);
% % % disp(['k_sat ML sat = ' num2str(ks)])
% % % disp(['gamma ML sat = ' num2str(ga2(ks))])
% % % 
% % % 
% % % %% Plot the results
% % % 
% % % figure(2)
% % % semilogy(de)
% % % grid
% % % xlabel('k_{sat}')
% % % ylabel('ML target function')
% % % title('best k_{sat} value')
% % % 
% % % figure(3)
% % % % data
% % % loglog(k,Pk,'.')
% % % hold on
% % % % ML fitting (we make sure that the plot follows the data)
% % % s1 = k.^(1-gamma); % build the CCDF signal
% % % loglog(k,s1/s1(150)*Pk(150));
% % % % ML fitting with saturation
% % % s1 = ((k+ks)/(k_min+ks)).^(1-ga2(ks));
% % % loglog(k,s1)
% % % hold off
% % % axis([xlim min(Pk/2) 2])
% % % grid
% % % xlabel('k')
% % % ylabel('CCDF')
% % % title('ML fittings')
% % % legend('data','ML','ML with sat.')

%% Probability of connecting to node i

p_conn = d/sum(d);

% % % %% Degree distribution for Barabàsi-Albert (BA) model
% % % 
% % % m = 1:5; % number of trials for connections
% % % pdf_BA = [pdf' zeros(M,length(m)-1)];
% % % pdf_BA(:,2:size(pdf_BA,2)) = (2*m.^2)./(pdf(:,1).^3);
% % % 
% % % % Plotting the BA distribution
% % % figure('Name','BA model log distribution')
% % % loglog(pdf_BA(:,1),pdf_BA(:,2),'o');
% % % hold on
% % % for i = 3:size(pdf_BA,2)
% % %     loglog(pdf_BA(:,1),pdf_BA(:,i),'o');
% % % end
% % % hold off
% % % 
% % % % DOMANDA: perchè mi vengono probabilità > 1?????
% % % 
% % % % Expected gamma exponent A/m
% % % gamma_BA = size(A,1)./m;


% % % % generating gen_m random values according to the distribution pdf
% % % gen_m = n_nodes;
% % % x = rand(1,gen_m);
% % % edges = cumsum([0; pdf_gen(:,2)]);
% % % y = histc(x,edges);
% % % 
% % % figure()
% % % plot(y,'o')

%% %%%%%%%%%%%%%%%%%% ASSORTATIVITY ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%

%% Avg neigh degree k_nn (it is given by k_nn = mean(q_m))

% First, remove the nodes that are isolated
pos = find(sum(A)~=0);
A_red = A(pos,pos);
N = size(A_red,1);

% Then find the largest directed component in the graph
%%%%%%%%%%%%%%%% Riattivare una volta capito qual è l'errore %%%%%%%%%%%%%%
if directed
    e1 = [1;zeros(N-1,1)];
    exit = false;
    while(~exit)
        e1_old = e1;
        e1 = 1*(A_red*e1+e1>0);
        exit = (sum(e1-e1_old)==0);
    end
    idx = find(e1);
    A_red = A_red(idx,idx);
    N = size(A_red,1);
end

% Estimation of avg k_nn = mean on degree

% degrees
d_red_in = sum(A_red,2);
d_red_out = sum(A_red)';

% averages of neighbours
k_tmp_oo = (A_red'*d_red_out)./d_red_out;
k_tmp_oi = (A_red'*d_red_in)./d_red_out;
k_tmp_ii = (A_red*d_red_in)./d_red_in;
k_tmp_io = (A_red*d_red_out)./d_red_in;

% extract averages for each value of k
u_out = unique(d_red_out);
for i = 1:length(u_out)
    k_nn_oo(i) = mean(k_tmp_oo(d_red_out==u_out(i)));
    k_nn_oi(i) = mean(k_tmp_oi(d_red_out==u_out(i))); % this is the mean of the neighbouring degree 
                                                      % of all the nodes having the same degree
end

u_in = unique(d_red_in);
for i = 1:length(u_in)
    k_nn_io(i) = mean(k_tmp_io(d_red_in == u_in(i)));
    k_nn_ii(i) = mean(k_tmp_ii(d_red_in == u_in(i)));
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
loglog(d_red_out,k_tmp_oo,'b.');
hold on
loglog(u_out,exp(p_oo(2)+log(u_out)*p_oo(1)),'r-');
loglog(u_out,k_nn_oo,'g.');
grid
xlabel('k_{out}')
ylabel('k_{nn,out}')
title('Network Assortativity. out --> out')

subplot(2,2,2)
loglog(d_red_out,k_tmp_oi,'b.');
hold on
loglog(u_out,exp(p_oi(2)+log(u_out)*p_oi(1)),'r-');
loglog(u_out,k_nn_oi,'g.');
grid
xlabel('k_{out}')
ylabel('k_{nn,out}')
title('Network Assortativity. out --> in')

subplot(2,2,3)
loglog(d_red_in,k_tmp_io,'b.');
hold on
loglog(u_in,exp(p_io(2)+log(u_in)*p_io(1)),'r-');
loglog(u_in,k_nn_io,'g.');
grid
xlabel('k_{in}')
ylabel('k_{nn,in}')
title('Network Assortativity. in --> out')

subplot(2,2,4)
loglog(d_red_in,k_tmp_ii,'b.');
hold on
loglog(u_in,exp(p_ii(2)+log(u_in)*p_ii(1)),'r-');
loglog(u_in,k_nn_ii,'g.');
grid
xlabel('k_{in}')
ylabel('k_{nn,in}')
title('Network Assortativity. in --> in')


%% %%%%%%%%%%%%%%%%% Ranking %%%%%%%%%%%%%%%%%%%%

% Create the equivalent undirected network, only if the network is directed
if directed
    Au = 1*(A+A'>0); % undirected network
    idx = find(sum(Au)~=0);
    A_red = A_red(idx,idx);
    Au = Au(idx,idx);
    
    % remove dead ends
    exit = false;
    while (~exit)
        idx = find(sum(A)~=0);
        A_red = A_red(idx,idx);
        Au = Au(idx,idx);
        N = size(A_red,1);
        exit = isempty(find(sum(A_red)==0, 1));
    end
else
    Au = A_red;
end

% find the largest connected component for the undirected graph
if directed
    e1 = [1;zeros(N-1,1)];
    exit = false;
    while(~exit)
        e1_old = e1;
        e1 = 1*(Au*e1>0);
        exit = (sum(e1-e1_old)==0);
    end
    idx = find(e1);
    A_red = A_red(idx,idx);
    N = size(A_red,1);
end

%% PageRank equation evaluation 

c = 0.85; %damping factor
q = ones(N,1)/N; %teleportation vector
d_out = A_red'*ones(size(A_red,1),1); %output degree vector
AdjM = A_red*sparse(diag(d_out.^(-1))); %Weighted adjacency matrix
M_11 = c*AdjM;
M_12 = (1-c)*q;

% Through linear system solution
tic
pr_ls = (sparse(eye(size(AdjM,1))) - M_11)\((1-c)*q);
pr_ls = pr_ls/sum(pr_ls);
toc

% Through power iteration
t = 40; % number of iterations
tic
pr = ones(N,1)/N; % initial guess on p_0(i)
s = zeros(t,1);
for i = 1:t
    pr_old = pr;
    % iterative step
    pr = M_11*pr + M_12;
    pr = pr/sum(pr);
    s(i) = norm(pr - pr_ls)/sqrt(N);
end
toc
distance = s;

%%%%%%%%%%%% Checking for convergence %%%%%%%%%%%%%

% Extract eigenvalues
disp('Extracting the eigenvalues')
tic
lambdas = eigs(AdjM,size(AdjM,1));
toc

lambda_2 = lambdas(2); %take the second eigenvalue, since the first one is 1 and stands for the connected component

thr = (c*abs(lambda_2)).^(1:t);
convergence = (distance' <= thr); %Checking if there is convergence in the solution

figure('Name','Convergence of the power iteration method for PageRank equation')
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
semilogy((1:t),distance)
grid
xlabel('k [iteration \#]')
ylabel('$\|r_k - r_\infty\|$')
title('PageRank convergence')

figure('Name','Eigenvalues representaion')
plot(real(lambdas(2:end)),imag(lambdas(2:end)),'o')
hold on
viscircles([0 0],lambdas(1),'Color','g');
hold off
grid

%% HITS equation evaluation - authorities scores

AdjM = sparse(A*A');

% Through linear system solution
tic
[hits_ls_vect, hits_ls] = eigs(AdjM,2);
pr_dir_aut = -hits_ls_vect(:,1)/norm(hits_ls_vect(:,1)); % Extract the most importand dimension and project the values wrt to that (SVD)
toc

% Through power iteration
tic
N = size(AdjM,1);
hits = ones(N,1)/N; % initial guess
s = zeros(t,1);
for i = 1:t
    temp = hits;
    hits = AdjM*hits; % iterative step
    hits = hits/norm(hits); % normalization step
    s(i) = norm(hits - temp)/sqrt(N);
end
toc

thr_hits = (abs(hits_ls(2)/hits_ls(1))).^(1:t);
convergence_hits = (s <= thr); %Checking if there is convergence in the solution

figure('Name','Convergence of the power iteration method for Hits equation')
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
semilogy((1:t),s)
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
N = size(AdjM,1);
hits = ones(N,1)/N; % initial guess
s = zeros(t,1);
for i = 1:t
    temp = hits;
    hits = M_2*hits; % iterative step
    hits = hits/norm(hits); % normalization step
    s(i) = norm(hits - temp)/sqrt(N);
end
toc

thr_hits = (abs(hits_ls(2)/hits_ls(1))).^(1:t);
convergence_hits = (s <= thr); %Checking if there is convergence in the solution

%% %%%%%%%%%%% Comparing PageRank and Hits %%%%%%%%%%%%%%

figure('Name','Comparison between PageRank and Hits authorities scores')
plot(pr_dir_aut/sum(pr_dir_aut),pr_ls/sum(pr_ls),'o')
grid
title('PageRank vs HITS ')

figure('Name','Comparison between PageRank and hits authorities scores')
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
plot(1:N,pr_ls/sum(pr_ls))
grid
hold on
plot(1:N,-pr_dir_aut/sum(pr_dir_aut))
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

%% %%%%%%%%%%%%%%%%%%% SPECTRAL CLUSTERING %%%%%%%%%%%%%%%%%%%%%

% food web -> 0.19
% https://www.nceas.ucsb.edu/interactionweb/html/thomps_towns.html
% G = importdata('./Narrowdaletxt.txt', '\t');
% W = sparse(1*(G>0));
% clear G;

%% Preprocessing the network

% Remove non-connected components
pos = find(sum(W,1));
W = W(pos,pos);

% Equivalent undirected network, if it is the case
if directed
    Wu = 1*(W+W'>0);
    Wu = Wu - diag(diag(Wu)); % clear diagonal
    pos = find(sum(Wu,1));
    Wu = Wu(pos,pos);
else
    Wu = W;
end
% Network size
N = size(Wu,1);

% Spectral domain info
d = full(sum(Wu)); % the degree vector
D = diag(d); % degree diagonal matrix
% spidiags(ones(N,1),0,N,N)
L = eye(N) - D^(-0.5)*Wu*D^(-0.5); % the normalized Laplacian matrix

if N < 2e3
    [eigv, lambdas] = eig(L); %extracting eigenvalues and eigenvectors
    
    figure('Name','Eigenvalues of the Laplacian')
    plot(diag(lambdas),'o');
    grid
    title('Eigenvalues of the normalized Laplacian')
else
    [eigv, lambdas] = eigs(L,3);
end

% Normalize the eigenvectors
Neigv = D*eigv;
% Find the Fiedler vector x_N-1 and normalize it
fv1 = Neigv(:,2)/norm(Neigv(:,2));
% Extract the following vector
fv2 = Neigv(:,3)/norm(Neigv(:,3));

% Identifying communities wrt to fv1
comm1 = find(fv1>0);
comm2 = find(fv1<0);

%% Reordering the nodes wrt fiedler's vector

% % % B = [fv1 (1:size(Wu,1))'];
% % % bsort(:,1) = sort(B(:,1));
% % % for i = 1:length(bsort)
% % %     temp = find(B(:,1) == bsort(i,1));
% % %     bsort(i,2) = B(temp,2); % bsort has the ordered values in the first column
% % %                             % and the indexes in the Wu matrix in the
% % %                             % second one
% % %     Wu_ord(i,:) = Wu(bsort(i,2),:); % Adjacency matrix whose rows have been 
% % %                                     % reordered following the fiedler's
% % %                                     % vector
% % % end

% reorder the adjacency matrix
[B,pos] = sort(fv1);
Wu_ord = Wu(pos,pos);

% Finding the conductance measure
a = sum(triu(Wu_ord)); % the degree of the columns
b = sum(tril(Wu_ord)); % the degree of the rows
d = a + b; % the total degree
newD = sum(d);

% association update
assoc = cumsum(d);
den = min(assoc,newD - assoc);

% cut update
cut = cumsum(b - a);

% the conductance
cond = cut./den;
cond = cond(1:end-1);

% plotting the conductance measure
figure('Name','Conductance measure')
plot(cond,'g--')
grid
title('Conductance')

% minimum conductance
disp(['the minimum conductance is: ', num2str(min(cond))]);

% The cheeger's bound
cb = sqrt(2*lambdas(2,2));

disp(['The cheeger''s bound is : ', num2str(cb)]);
if (min(cond) < cb)
    disp('Congrats, you found a good community!');
else
    disp('Are you sure it is the right way?');
end

%% 

% show network with partition
[mincond,mpos] = min(cond);
threshold = mean(B(mpos:mpos+1));

figure('Name','Network communities')
plot(fv1,fv2,'.')
grid
hold on
plot(threshold*[1,1],ylim,'r-')
hold off
title('communities')

% show network with partition and links
if newD < 1e4 % only if edges are not too many!!!
    figure(4)
    plot(fv1,fv2,'.')
    grid
    [I,J,~] = find(Wu);
    hold on
    plot([fv1(I),fv1(J)]',[fv2(I),fv2(J)]')
    plot(threshold*[1,1],ylim,'k-')
    hold off
    title('communities (with links)')
end

% save network of largest community 
if (mpos>=N-mpos)
    A = W(pos(1:mpos),pos(1:mpos));
else
    A = W(pos(mpos+1:end),pos(mpos+1:end));
end
save('previous_community','W')

toc







