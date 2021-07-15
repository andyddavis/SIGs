clear

% number of nodes in the line
n = 9;

% initialise the probabilities
p = 1/4;         % probability of moving right
k = 1/4;        % probability of staying stationary
p = (1-k) / 2;  % uniform dist p

% initialise the transition matrix
tMatrix = zeros(n,n);
for i = 1 : n
    tMatrix(i,i) = k;
end
for i = 2 : n-1
    tMatrix(i,i-1) = 1-p-k;
    tMatrix(i,i+1) = p;
end
tMatrix(1,1) = 1-p;
tMatrix(1,2) = p;
tMatrix(n,n-1)=k;
tMatrix(n,n) = 1-k;

% initialise the agents on the far left
A = 100;
agents = ones(1,A);

% simulation data
time = 100;
dt = 1;
timesteps = time/dt;

% initial prob. of being at node i 
pi = rand(1, n);
pi = pi/sum(pi);

% live plotting
for timestep = 1 : timesteps
    
    pi = pi*tMatrix;
    
    for agent = 1 : A
        counter = 1;
        rVal = rand();
        cVal = tMatrix(agents(agent), counter);
        while (cVal < rVal)
            counter = counter + 1;
            cVal = cVal + tMatrix(agents(agent), counter);
        end
        agents(agent) = counter;
    end
    
    histogram(agents,n,'BinWidth', 1);
    xlim([1,n+1]);
    ylim([0,A+1]);
    pause(0.1);
    
    %histogram(agents,n, 'Normalization', 'pdf');%,'BinWidth', 1);
    %xlim([1,n+1]);
    %ylim([0,1]);
    %pause(0.5);
    %drawnow;
    
end