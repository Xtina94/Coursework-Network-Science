clear all;
close all;

% Load the dataset in example
load bioCel.txt
% % % G = importdata('wiki-Vote.txt', '\t', 4);
% % % bioCel = G.data;

%% Extract the number of nodes and links
nodes_firstCol = sort(bioCel(:,1));
nodes_secondCol = sort(bioCel(:,2));
nodes = [nodes_firstCol, nodes_secondCol];
nodes = unique(nodes);
n_nodes = length(nodes); % number of nodes

degree = zeros(n_nodes,2); % a n_nodes x 2 matrix with each row having nodeName --> nodeDegree
for i = 1:n_nodes
    degree(i,:) = [nodes(i), histc(bioCel(:,1),nodes(i))];
end

n_links = sum(degree(:,2));

%% Extract the PDF and the CCDM measure
A = sparse(bioCel(:,2),bioCel(:,1),ones(size(bioCel,1),1),n_nodes,n_nodes);

figure('Name','Sparse representation of the adjacency matrix')
spy(A)

% Linear PDF plot
k = unique(degree(:,2)); % all the degree types
k_gen = k; %saved value of degrees for later estimation on the BA model network generation

% Remove the first element of k since 0 is not an acceptable degree
k = k(2:length(k));
k_min = min(k);
M = length(k);
pdf_gen = zeros(M+1,2); % a n_nodes x 2 matrix with each row having degree --> #Nodes with that degree/tot nodes

for i = 1:M+1
    pdf_gen(i,:) = [k_gen(i), histc(degree(:,2),k_gen(i))];
end
pdf_gen(:,2) = pdf_gen(:,2)./n_links;
pdf = pdf_gen(2:size(pdf_gen,1),:); %remove the first value which corresponds to a degree = 0

figure('Name','Linear pdf plot')
scatter(pdf(:,1),pdf(:,2));

%PDF logarithmic plot
figure('Name','Logarithmic pdf plot')
loglog(pdf(:,1),pdf(:,2),'oblue');

%PDF log-binning plot

% % BinNum = 25; %Like professor's chart
% % [midpoints,Freq,eFreq] = lnbin(degree, BinNum);

BinNum = 8;
n = M/BinNum;
k_bins = zeros(BinNum,2);
for i = 1:n:M
    k_bins(ceil(i/n),1) = (k(i)+k(i+n-1))/2; %vector containing the average degree in the bin observed
    k_bins(ceil(i/n),2) = (pdf(i)+pdf(i+n-1))/2; %vector containing the average pdf in the bin observed
end

figure('Name','Log-binning plot')
loglog(k_bins(:,1),k_bins(:,2),'o');
% axis([0 250 0 0.05]);

%% PDF Logarithmic CCDM plot

%With fixed gamma = 3.5
gamma = 3.5;
c = (gamma-1)*k_min^(gamma-1);


%With ML estimation of gamma
j = degree(:,2)/k_min;
for i = 1:length(degree) 
    if j(i) == 0
        j(i) = 1;
    end
end
den = sum(log(j));
gamma_ML = 1+(n_nodes/den);


ccdm_pdf_35 = zeros(M,2);
ccdm_pdf_35(:,1) = pdf(:,1);
ccdm_pdf = zeros(M,2);
ccdm_pdf(:,1) = pdf(:,1);

p_k_35 = k.^(-gamma-1);
p_k = k.^(-gamma_ML-1);

for i = 1:M
    ccdm_pdf_35(i,2) = sum(p_k_35(i:M));
    ccdm_pdf(i,2) = sum(p_k(i:M));
end

% % % ccdm_pdf_35(:,2) = c*k.^(1-gamma);
% % % ccdm_pdf(:,2) = c*k.^(1-gamma_ML);

figure('Name','Logarithmic CCDM_35')
loglog(ccdm_pdf_35(:,1),ccdm_pdf_35(:,2),'o');

figure('Name','Logarithmic CCDM')
loglog(ccdm_pdf(:,1),ccdm_pdf(:,2),'o');

%% Probability of connecting to node i

p_conn = degree(:,2)/sum(degree(:,2));

%% Degree distribution for Barab�si-Albert (BA) model

m = 1:5; % number of trials for connections
pdf_BA = [pdf zeros(M,length(m)-1)];
pdf_BA(:,2:size(pdf_BA,2)) = (2*m.^2)./(pdf(:,1).^3);

% Plotting the BA distribution
figure('Name','BA model log distribution')
loglog(pdf_BA(:,1),pdf_BA(:,2),'o');
hold on
for i = 3:size(pdf_BA,2)
    loglog(pdf_BA(:,1),pdf_BA(:,i),'o');
end
hold off

% DOMANDA: perch� mi vengono probabilit� > 1?????

% Expected gamma exponent A/m
gamma_BA = size(A,1)./m;


% generating gen_m random values according to the distribution pdf
gen_m = n_nodes;
x = rand(1,gen_m);
edges = cumsum([0; pdf_gen(:,2)]);
y = histc(x,edges);

figure()
plot(y,'o')


