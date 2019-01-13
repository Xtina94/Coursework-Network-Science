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
        Au = 1*(A+A'>0); % undirected network
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
    case 2
        G = importdata('bioCel.txt');
        N = max(max(G));
        A = sparse(G(:,2),G(:,1),ones(size(G,1),1),N,N);
        W = A;
        Au = 1*(A+A'>0); % undirected network
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
    case 3
        G = importdata('ca_sandi_auths.txt');
        N = max(max(G));
        A = sparse(G(:,2),G(:,1),ones(size(G,1),1),N,N);
        W = A;
        Au = 1*(A+A'>0); % undirected network
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
    case 4
        G = importdata('collaboration.edgelist.txt', '\t', 1);
        G.data = G.data + 1;
        N = max(max(G.data));
        A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
        A = 1*(A+A'>0); % build undirected network
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
    case 5
        A = importdata('occupyWs.txt');
        N = max(max(A));
        A = sparse(A(:,2),A(:,1),ones(size(A,1),1),N,N);
        Au = 1*(A+A'>0); % undirected network
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
end

%% %%%%%%%%%%%%%%%% DISTRIBUTION EXTRACTION AND EXPONENTIAL LAW %%%%%%%%%%%%%%%

% distribution. degree is a 1 column vector with node degrees
d = full(sum(A,1));
% Remove the first elements of k since 0 is not an acceptable degree
d = d(d > 0);
k = unique(d); % degree samples
k_min = min(k);
k_max = max(k);
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

saveas(gcf,'pdfPlots.png')

%% PDF Logarithmic CCDM plot

%With fixed gamma = 3.5
gamma = 3.5;
c = (gamma-1)*k_min^(gamma-1);

%With ML estimation of gamma
k_min_temp = floor(max(d)/2);
d2 = d(d>=k_min_temp); % restrict range

j = d2/k_min_temp;
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
xlabel('k');
ylabel('log(pdf)');
subplot(1,2,2)
thr = datasample(interval,1);
loglog(k,Pdf,'o')
hold on
loglog(k,ccdm_pdf_35/ccdm_pdf_35(thr)*Pdf(thr),'--');
loglog(k,ccdm_pdf/ccdm_pdf(thr)*Pdf(thr),'--');
hold off
axis([xlim min(pklog/2) 2*max(Pdf)])
grid
xlabel('k');
ylabel('log(CCDF)');

saveas(gcf,'MLestimate.png')

%% ML fitting with saturation

for i = 1:max(k)
    tmp = mean(log((d+i)/(k_min+i)));
    ga2(i) = 1+1/tmp;
    de(i) = log(ga2(i)-1)-log(k_min+i)-ga2(i)*tmp;
end

[~,i] = max(de);
disp(['k_sat ML sat = ' num2str(i)])
disp(['gamma ML sat = ' num2str(ga2(i))])

%% Plot the results

figure('Name','ML target function')
semilogy(de)
grid
xlabel('k_{sat}')
ylabel('ML target function')
title('best k_{sat} value')

saveas(gcf,'MLtargetFunction,png');

figure('Name','cumulative Pdf')
% data
loglog(k,Pdf,'.')
hold on
% ML fitting (we make sure that the plot follows the data)
s1 = k.^(1-gamma_ML); % build the CCDF signal
loglog(k,s1/s1(thr)*Pdf(thr));
% ML fitting with saturation
s1 = ((k+i)/(k_min+i)).^(1-ga2(i));
loglog(k,s1)
hold off
axis([xlim min(Pdf/2) 2])
grid
xlabel('k')
ylabel('CCDF')
title('ML fittings')
legend('data','ML','ML with sat.')

 saveas(gcf,'MLfitting.png')

%% Estimation of other parameters

%%%% Number of nodes %%%%
n_nodes = N;
disp(['The number of nodes is: ', num2str(n_nodes)]);

%%%% Number of links %%%%
disp(['The total number of links is: ', num2str(n_links)]);

%%%% Moments of the degree distribution %%%%

% average degree <k>
avg_k = mean(k);
% variance of the degree <k^2> (the spread)
sigma_k = var(k);
% skewness (for symmetry)
skew_k = pdf'*(k.^3)';

%%%% Average distance %%%%
if (2 < gamma_ML) && (gamma_ML < 3)
    avg_dist = log(log(N));
elseif gamma_ML == 3
    avg_dist = log(N)/log(log(N));
elseif gamma_ML > 3
    avg_dist = log(N);
end
        
disp(['Moments of the degree dist (average, variance, skewness): ', num2str(avg_k), ' ', num2str(sigma_k), ' ' num2str(skew_k)]);
disp(['Average distance: ',num2str(avg_dist)]);

%%%%% Natural cutoff estimation %%%%%
k_nat = k_min*N^(1/(gamma_ML-1));

%%%%% Inhomogeneity ratio %%%%%
inom = sigma_k/avg_k;

%%%% Breaking point estimation %%%%%
if gamma_ML > 3
    f_c = 1 - 1/((gamma_ML-2)/(gamma_ML - 3)*k_min - 1);
else
    f_c = 1 - 1/((gamma_ML-2)/(3-gamma_ML)*k_min^(gamma_ML-2)*k_max^(3-gamma_ML) - 1);
end

save('Outputs.mat','n_nodes','n_links','avg_k','sigma_k','skew_k','inom');

%% %%%%%%%%%%%%%%%%%% ASSORTATIVITY ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%

%% Avg neigh degree k_nn (it is given by k_nn = mean(q_m))

% First, remove the nodes that are isolated
pos = find(sum(A)~=0);
A_red = A(pos,pos);
A_isol = A_red;
N_red = size(A_red,1);

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

% Take the average
assort_avg = mean([p_oo(1), p_io(1), p_oi(1), p_ii(1)]);

%% Plot the results
figure('Name','Assortativity behaviour')
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
% legend('Average degree of neighbours','Average neighbouring degree', 'Assortativity interpolation');
title('Network Assortativity. in --> in')

saveas(gcf,'Assortativity.png')

%% %%%%%%%%%%%%%%%%%%%%% SPECTRAL CLUSTERING %%%%%%%%%%%%%%%%%%%%%%%

% Equivalent undirected network, if it is the case
if directed
    pos = find(sum(Au)~=0);
    Au = Au(pos,pos);
    A_red = A(pos,pos);
else
    pos = find(sum(A)~=0);
    A_red = A(pos,pos);
    Au = A_red;
end
% Network size
N_red = size(Au,1);

% Spectral domain info
d = full(sum(Au)); % the degree vector
d = ones(N_red,1)./sqrt(d);
D = spdiags(d,0,N_red,N_red); % degree diagonal matrix
L = spdiags(ones(N_red,1),0,N_red,N_red) - D*Au*D; % the normalized Laplacian matrix

if N_red < 2e3
    [eigv, lambdas] = eigs(L); %extracting eigenvalues and eigenvectors
    
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

% reorder the adjacency matrix
[B,pos] = sort(fv1);
Wu_ord = Au(pos,pos);

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
ylabel('Cond');

saveas(gcf,'Conductance.png');

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

%% Network partitioning

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

saveas(gcf,'Communities.png');

% show network with partition and links
if newD < 1e4 % only if edges are not too many!!!
    figure('Name','Communities connected')
    plot(fv1,fv2,'.')
    grid
    [I,J,~] = find(Au);
    hold on
    plot([fv1(I),fv1(J)]',[fv2(I),fv2(J)]')
    plot(threshold*[1,1],ylim,'k-')
    hold off
    title('communities (with links)')
    
    saveas(gcf,'NetworkCommunities.png');
end

% save network of largest community 
if (mpos>=N_red-mpos)
    A_red = A_red(pos(1:mpos),pos(1:mpos));
else
    A_red = A_red(pos(mpos+1:end),pos(mpos+1:end));
end


%% %%%%%%%%%%%%%%% Clustering coefficient %%%%%%%%%%%%%%%%%
% % % N_red = 5; % decomment only in case of algorithm check
% % % Au = [0 1 1 0 0;
% % %       1 0 1 0 0;
% % %       1 1 0 1 1;
% % %       0 0 1 0 1;
% % %       0 0 1 1 0];

% Counting the triangles which each node belongs to
tri = zeros(N_red,1);
upper_tri = triu(Au);
for i = 1:N_red
    pos = find(Au(i,:));   
    comp = Au(i,:);
    for j = 1:length(pos)
        pos_2 = upper_tri(pos(j),:);
        temp = comp;
        temp(pos(j)) = 0;
        bridge = find(temp + pos_2 == 2);
        tri(i) = tri(i) + length(bridge);
    end
end

% total number of triangles
n_tri = sum(tri);

%% Clustering coefficient
d = (full(sum(Au,1)))';
C = 2*tri./(d.*(d-1));
rm = find(isnan(C));
C(rm) = 0;
C_u = unique(C);

% Average clustering coefficient
C_avg = mean(C);

% plot the distribution
pdf_C = histc(C,C_u);
figure('Name','pdf of Clustering coeff')
plot(C_u,pdf_C,'.')
grid
xlabel('C');
ylabel('pdf(C)');

saveas(gcf,'CLusteringCpdf.png');

save('OutputClustering','cond','A','n_tri','C_avg');

%% %%%%%%%%%%%%%% Clustering coefficient for the whole network %%%%%%%%%%%%%%%%%

% Counting the triangles which each node belongs to
tri = zeros(N,1);
if directed
    Au = 1*(A+A'>0); % undirected full network
else 
    Au = A;
end
upper_tri = triu(Au);
for i = 1:N
    pos = find(Au(i,:));   
    comp = Au(i,:);
    for j = 1:length(pos)
        pos_2 = upper_tri(pos(j),:);
        temp = comp;
        temp(pos(j)) = 0;
        bridge = find(temp + pos_2 == 2);
        tri(i) = tri(i) + length(bridge);
    end
end

% total number of triangles
n_tri = sum(tri);

%% Clustering coefficient
d = (full(sum(Au,1)))';
C = 2*tri./(d.*(d-1));
rm = find(isnan(C));
C(rm) = 0;
C_u = unique(C);

% Average clustering coefficient
C_avg = mean(C);

% plot the distribution
pdf_C = histc(C,C_u);
figure('Name','Pdf of Clustering coeff for the full network')
plot(C_u,pdf_C,'.')
grid
saveas(gcf,'CLusteringCpdfFull.png');
xlabel('C');
ylabel('pdf(C)');

save('OutputClusteringull','cond','A','n_tri','C_avg');