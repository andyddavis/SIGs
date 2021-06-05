clear
% number of nodes in the line
n = 9;

% initialise the transition matrix
tMatrix = zeros(n,n);
for i = 1 : n
    tMatrix(i,i) = 1/3;
end
tMatrix(1,2) = 2/3;
for i = 2 : n-1
    tMatrix(i,i-1) = 1/6;
    tMatrix(i,i+1) = 1/2;
end
tMatrix(n,n-1)=1/3;
tMatrix(n,n) = 2/3;

% initialise the agents on the far left
A = 20;
agents = ones(1,A);

% simulation data
time = 50;
dt = 1;
timesteps = time/dt;

% live plotting
for timestep = 1 : timesteps
    
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
    pause(0.5);
    
end
