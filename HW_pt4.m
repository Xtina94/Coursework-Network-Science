%% %%%%%%%%%%%%%%%%%%% SPECTRAL CLUSTERING %%%%%%%%%%%%%%%%%%%%%
close all
clear all

y = 1;
switch y
    case 1
        G = importdata('wiki-Vote.txt', '\t', 4);
        N_red = max(max(G.data));
        A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N_red,N_red);
        clear G;
        Au = 1*(A+A'>0); % undirected network
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
        clear G;
        Au = 1*(A+A'>0); % undirected network
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
        clear G;
        Au = 1*(A+A'>0); % undirected network
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
    case 4
        G = importdata('./Narrowdaletxt.txt', '\t');
        A = sparse(1*(G>0));
        Au = 1*(A+A'>0); % undirected network
        clear G;
        if (isequal(A,Au))
            directed = 0;
        else
            directed = 1;
        end
end

%% Preprocessing the network

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
    [I,J,~] = find(Au);
    hold on
    plot([fv1(I),fv1(J)]',[fv2(I),fv2(J)]')
    plot(threshold*[1,1],ylim,'k-')
    hold off
    title('communities (with links)')
end

% save network of largest community 
if (mpos>=N_red-mpos)
    A_red = A_red(pos(1:mpos),pos(1:mpos));
else
    A_red = A_red(pos(mpos+1:end),pos(mpos+1:end));
end
save('previous_community','A')







