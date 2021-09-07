% fixed graph that discretizes (compact) space 
% sample n nodes
% connect to k nearest neighbors 
% parameters: n, k

n = 1000; % number of nodes 
node_x = rand(n,1); 
node_y = rand(n,1);
node = horzcat(node_x, node_y); % 1000 pts from uniform distribution (0,1) 

V = zeros(n,2); % velocity field evaluated at node locations 
% for i = 1:n
% %     V(i,1) = node_x(i)/(2*pi*(node_x(i)^2+node_y(i)^2));
% %     V(i,2) = node_y(i)/(2*pi*(node_x(i)^2+node_y(i)^2));
%     V(i,1) = 1/(4*pi)*exp(-((node_x(i)-0.5)^2+(node_y(i)-0.5)^2)/(4))*(-2*(node_x(i)-0.5)/4);
%     V(i,2) = 1/(4*pi)*exp(-((node_x(i)-0.5)^2+(node_y(i)-0.5)^2)/(4))*(-2*(node_y(i)-0.5)/4);
% end
dx = 1/(4*pi)*exp(-((node_x -0.5).^2+(node_y -0.5).^2)./4).*(-2.*(node_x -0.5)./4);
dy = 1/(4*pi)*exp(-((node_x -0.5).^2+(node_y -0.5).^2)./4).*(-2.*(node_y -0.5)/4);
V(:,1) = dx; 
V(:,2) = dy; 

% % different vector fields to play with 
% [x,y] = meshgrid(0:.1:1,0:.1:1);

% % 1. for vector field toward top right corner
% figure; 
% dfx = 1/(4*pi)*exp(-(x.^2+y.^2)./(4))*((x./2)) + 0.5;
% dfy = 1/(4*pi)*exp(-(x.^2+y.^2)./(4))*((y./2)) + 0.5;
% quiver(x, y, dfx, dfy); 
% xlim([0 1]) 
% ylim([0 1])

% 2. movement toward center 
% figure; 
% dfx = 1/(4*pi)*exp(-((x-0.5).^2+(y-0.5).^2)./(4))*(-2.*(x-0.5)./4);
% dfy = 1/(4*pi)*exp(-((x-0.5).^2+(y-0.5).^2)./(4))*(-2.*(y-0.5)./4);
% quiver(x, y, dfx, dfy); 
% xlim([0 1]) 
% ylim([0 1])

k = 10; % k-1 closest neighbors; (k - 1) because trivially the closest neighbor is itself
dist = zeros(n,n); % symmetric matrix of distances between two nodes, dist(i,j) = distance between node i and node j
ind = zeros(n,k); % index matrix where ind(i,j) = the jth (out of k) closest neighbor of node i
node_nb_x = zeros(n,k); % x-coord of neighbors of node i=1:n
node_nb_y = zeros(n,k); % y-coord of neighbors of node i=1:n

for i = 1:n
for j = 1:n 
    dist(i,j) = norm(node(i,:) - node(j,:)); % distance between i and j
end
end

w = zeros(n,n); % holds velocities projected onto edges where w(i,j) is the velocity at node i projected onto the edge from i to j
A = zeros(n,n); 

figure; % for graph
for i = 1:n
[M, I] = mink(dist(i,:),k); % M contains the k nearest neighbors of node i and I contains their indices 
ind(i,:) = I;  % store the indices of the k nearest neighbors of node i in the ith row of matrix ind
node_nb_i = [node_x(ind(i,:)), node_y(ind(i,:))]; % holds the locations of the k nearest neighbors of node i
for j = 1:k
    e_ik = node_nb_i(j,:) - node(i,:); % edge vector is going from node i to the jth closest neighbor of node i 
    w(i,ind(i,j)) = max(dot(e_ik, V(i,:)),0)/norm(e_ik); 
    
    A(ind(i,j),i) = max(dot(e_ik, V(i,:)),0)/norm(e_ik);
end 
scatter(node_nb_i(:,1), node_nb_i(:,2)) % comment/uncomment to plot/turn-off plot graph
line(node_nb_i(:,1), node_nb_i(:,2)) 

hold on
end
% to visualize velocity vectors evaluated at nodes on the graph plot 
hold on 
quiver(node_x,node_y,dx,dy)
xlim([0 1]) 
ylim([0 1])

TF = isnan(A);
A(TF) = 0; 
TF = isnan(w);
w(TF) = 0; 

d_test = zeros(n,1); % test using velocity matrix w for the diagonal 
for i = 1:n
    d_test(i) = sum(w(i,:));
end

d_test2 = zeros(n,1); % test using A for the diagonal 
for i = 1:n
d_test2(i) = sum(A(:,i)); 
end

if d_test == d_test2
    D = diag(d_test);
    D_inv = ones(n,n)./D;
    TF = isinf(D_inv); 
    D_inv(TF) = 0;
end

L = D - A; % L is the advection Laplacian 
P = A * D_inv; % P is a left stochastic matrix

for i = 1:n
    if sum(P(:,i)) == 0
        P(i,i) = 1;
    end
end

sanity_check = zeros(n,1); % to check if all columns sum to one 
for i = 1:n
sanity_check(i) = sum(P(:,i));
end 

%% 
mc = dtmc(P'); % need to do P transpose b/c P is left stochastic and dtmc takes in right stochastic matrices 
numSteps = 1000; % timesteps 
X0 = ones(n,1)./n; % uniform distribution initial condition 
% X0 = zeros(mc.NumStates,1); % delta distribution initial condition 
% X0(25) = 100; % 100 random walks starting from state 50 only

S = redistribute(mc,numSteps,'X0',X0); % Rows correspond to time steps, and columns correspond to states.
% figure; % for histogram plot of markov chain 
% dp = distplot(mc,S,'Type','histogram','FrameRate',0.5); % comment/uncomment to turn on/off visualization 

[r_eigenv, evalues, l_eigenv] = eig(P); % [right eigenvectors, eigenvalues, left eigenvectors]
evalues = diag(evalues); 
[evalues, ev_ind] = sort(evalues, 'descend'); 

lc = lazy(mc); % to force laziness i.e. self loops so that with prob 0.5, chain feels lazy and stays in same state 
st = asymptotics(lc);

% If mc is ergodic, and numSteps is sufficiently large, X(end,:) approximates x = asymptotics(mc)

% TO-DO: try with more nodes (probability of having very long edge will go
% diminishingly smaller) 
% savefig instead of visualizing plots 
% change definition of inverse from regular inverse to pseudoinverse 
% cross-check with making matrix P by defining each value separately 