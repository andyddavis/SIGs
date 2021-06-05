n = 9;

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

A = 9;
agents = ones(1,A);

time = 20;
dt = 1;
timesteps = time/dt;

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
    %plot(agents);
    xlim([1,n+1]);
    ylim([1,A+1]);
    pause(1);
    
end

histogram(agents,n,'BinWidth', 1);
  xlim([1,n+1]);
  ylim([1,A+1]);
  
