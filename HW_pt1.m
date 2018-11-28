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
% See it later

for ks = 1:max(k)
    k_min = min(d);
    tmp = mean(log((d+ks)/(k_min+ks)));
    ga2(ks) = 1+1/tmp;
    de(ks) = log(ga2(ks)-1)-log(k_min+ks)-ga2(ks)*tmp;
end
[~,ks] = max(de);
disp(['k_sat ML sat = ' num2str(ks)])
disp(['gamma ML sat = ' num2str(ga2(ks))])


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
s1 = ((k+ks)/(k_min+ks)).^(1-ga2(ks));
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

% Number of nodes
n_nodes = N;
disp(['The number of nodes is: ', num2str(n_nodes)]);

% Number of links
disp(['The total number of links is: ', num2str(n_links)]);

%%%%%%%%% Moments of the degree distribution %%%%%%%%%%%
% average degree <k>
avg_k = mean(k);
% variance of the degree <k^2> (the spread)
sigma_k = var(k);
% skewness (for symmetry)
skew_k = pdf'*(k.^3)';

disp(['Moments of the degree dist (average, variance, skewness): ', num2str(avg_k), ' ', num2str(sigma_k), ' ' num2str(skew_k)]);

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
