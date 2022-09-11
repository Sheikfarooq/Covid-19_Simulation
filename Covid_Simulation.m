%% Student names

% Student: Sheik Farooq (419396) 


%% Command to clean the screen
clc;
clear;
close All;

%% Task 1

% Initialization
a = 0.01; % minimum x and y values
b = 1.99; % maximum x value (Actually 2, but better representation in graph, 1.99 was used)
c = 0.99; % maximum y value (Actually 1, but better representation in graph, 0.99 was used)
n = 100;  % number of individuals
dt = 0.001;
vel_min = -1;
vel_max = 1;
infected_initial = 10; % in percentage

n_infec = infected_initial*n/100;    % no. of infected people
n_hel = (100-infected_initial)*n/100;  % no. of healthy people
n_recov = n-n_infec-n_hel;           % no. of recovered people

Time = 0:dt:1;  % Time vector

% initial X and Y positions of the points
pos_x_infec = (b-a).*rand(n_infec,1) + a;  % Healthy
pos_x_infec = [pos_x_infec, dt*ones(size(pos_x_infec))];
pos_y_infec = (c-a).*rand(n_infec,1) + a;

pos_x_hel = (b-a).*rand(n_hel,1) + a;  % Infected
pos_y_hel = (c-a).*rand(n_hel,1) + a;

pos_x_rec = [];
pos_y_rec = [];
vel_x_rec = [];
vel_y_rec = [];

% initial point positions of healthy and infected points shown in graph
figure(1);
hold on;
plot(pos_x_hel,pos_y_hel, 'o','MarkerFaceColor', 'g','MarkerEdgeColor', 'g')
plot(pos_x_infec(:,1),pos_y_infec, 'o','MarkerFaceColor', 'r','MarkerEdgeColor', 'r')

% Velocities
vel_x_hel = (vel_max-vel_min).*rand(n_hel,1) + vel_min;  % Velocities of Healthy
vel_y_hel = (vel_max-vel_min).*rand(n_hel,1) + vel_min;

vel_x_infec = (vel_max-vel_min).*rand(n_infec,1) + vel_min;  % Velocities of infected
vel_y_infec = (vel_max-vel_min).*rand(n_infec,1) + vel_min;


v=VideoWriter('Pandamic_1.avi');    % saves dynamic graph as a video file
open(v);


for k=1:length(Time)   % Main "for" loop (1 loop (dt)== 0.001 sec 
  
    % new positions of the persons based on velocity and old position
    
    pos_x_new_hel = pos_x_hel + vel_x_hel.*dt;   % Healthy
    pos_y_new_hel = pos_y_hel + vel_y_hel.*dt;
    pos_x_hel = pos_x_new_hel;
    pos_y_hel = pos_y_new_hel;
    
    
    pos_x_new_infec = pos_x_infec(:,1) + vel_x_infec.*dt;   % Unhealthy
    pos_y_new_infec = pos_y_infec + vel_y_infec.*dt;
    pos_x_infec = [pos_x_new_infec pos_x_infec(:,2)+dt];
    pos_y_infec = pos_y_new_infec;
    
    pos_x_new_rec = pos_x_rec + vel_x_rec.*dt;   % Healthy
    pos_y_new_rec = pos_y_rec + vel_y_rec.*dt;
    pos_x_rec = pos_x_new_rec;
    pos_y_rec = pos_y_new_rec;
    
    % Recovery of the infected points
    
    for j=1:length(pos_x_infec(:,1))
        if pos_x_infec(j,2)>=0.35
            pos_x_rec = [pos_x_rec; pos_x_infec(j,1)];
            pos_y_rec = [pos_y_rec; pos_y_infec(j,1)];
            vel_x_rec = [vel_x_rec; vel_x_infec(j,1)];
            vel_y_rec = [vel_y_rec; vel_y_infec(j,1)];
            pos_x_infec(j,:) = [];
            % pos_x_infec(j,2) = [];
            pos_y_infec(j) = [];
            vel_x_infec(j) = [];
            vel_y_infec(j) = [];
            pos_x_infec = [pos_x_infec;[0,0]];
            pos_y_infec = [pos_y_infec;0];
            vel_x_infec = [vel_x_infec;0];
            vel_y_infec = [vel_y_infec;0];
            j=j-1;
        end
    end
    
    % deleting unnecessory zero elements in infected arrays 
    % This deletes the unnecessory points at the origin
    
    z=find(pos_x_infec(:,1));           % provides the positions of non zero elements
    pos_x_infec = pos_x_infec(1:length(z),:); % deleting the elements, which are after non-zero elements (i.e removing zero elements)
    pos_y_infec(size(z)+1:end)=[]; % removing zero elements in the array
    vel_x_infec(size(z)+1:end)=[]; % removing zero elements in the array
    vel_y_infec(size(z)+1:end)=[]; % removing zero elements in the array
    
    
    % Checking social distancing ( code to changing the green elements to red colour)
    
    for i=1:length(pos_y_infec)
        for j=1:n_hel
            if (sqrt((pos_x_hel(j)-pos_x_infec(i)).^2+(pos_y_hel(j)-pos_y_infec(i)).^2))<=0.055
                pos_x_infec = [pos_x_infec; pos_x_hel(j), dt];
                pos_y_infec = [pos_y_infec; pos_y_hel(j)];
                vel_x_infec = [vel_x_infec; vel_x_hel(j)];
                vel_y_infec = [vel_y_infec; vel_y_hel(j)];
                pos_x_hel(j) = [];
                pos_y_hel(j) = [];
                vel_x_hel(j) = [];
                vel_y_hel(j) = [];
                pos_x_hel = [pos_x_hel;0];
                pos_y_hel = [pos_y_hel;0];
                vel_x_hel = [vel_x_hel;0];
                vel_y_hel = [vel_y_hel;0];
                n_infec = n_infec+1;
                i=i+1;
                n_hel = n_hel-1;
                j=j-1;
            end
        end
    end
    
    % deleting unnecessory zero elements in health arrays 
    % This deletes the unnecessory points at the origin
    z=find(pos_x_hel);           % provides the positions of non zero elements
    pos_x_hel(size(z)+1:end)=[]; % deleting the elements, which are after non-zero elements (i.e removing zero elements)
    pos_y_hel(size(z)+1:end)=[]; % removing zero elements in the array
    vel_x_hel(size(z)+1:end)=[]; % removing zero elements in the array
    vel_y_hel(size(z)+1:end)=[]; % removing zero elements in the array
    
    
    
    % Code for points to bounce back at the boundaries
    % checking, whether the new positions are with in range (X=0 to 2) (Y=0 to 1)
    for j=1:n_hel         % Healthy
        if pos_x_hel(j)>=b
            vel_x_hel(j)=-1*vel_x_hel(j);
        end
        if pos_x_hel(j)<=a
            vel_x_hel(j)=-1*vel_x_hel(j);
        end
        if pos_y_hel(j)>=c
            vel_y_hel(j)=-1*vel_y_hel(j);
        end
        if pos_y_hel(j)<=a
            vel_y_hel(j)=-1*vel_y_hel(j);
        end  
    end
    
    for j=1:length(pos_y_infec)   % infected
        if pos_x_infec(j)>=b
            vel_x_infec(j)=-1*vel_x_infec(j);
        end
        if pos_x_infec(j)<=a
            vel_x_infec(j)=-1*vel_x_infec(j);
        end
        if pos_y_infec(j)>=c
            vel_y_infec(j)=-1*vel_y_infec(j);
        end
        if pos_y_infec(j)<=a
            vel_y_infec(j)=-1*vel_y_infec(j);
        end  
    end
    
    for j=1:length(pos_x_rec)         % Recovered
        if pos_x_rec(j)>=b
            vel_x_rec(j)=-1*vel_x_rec(j);
        end
        if pos_x_rec(j)<=a
            vel_x_rec(j)=-1*vel_x_rec(j);
        end
        if pos_y_rec(j)>=c
            vel_y_rec(j)=-1*vel_y_rec(j);
        end
        if pos_y_rec(j)<=a
            vel_y_rec(j)=-1*vel_y_rec(j);
        end  
    end
    
    
    % Grapf for the moving points
    
    figure(2);
    axis([0,2,0,1])
    plot(pos_x_hel,pos_y_hel, 'o','MarkerFaceColor', 'g','MarkerEdgeColor', 'g')
    hold on
    plot(pos_x_infec(:,1),pos_y_infec, 'o','MarkerFaceColor', 'r','MarkerEdgeColor', 'r')
    plot(pos_x_rec,pos_y_rec, 'o','MarkerFaceColor', 'b','MarkerEdgeColor', 'b')
    frame = getframe(gcf);
    writeVideo(v,frame);
    hold off
    
    % Data of no. of healthy, infected, recovered people along with time period is stored in the below matrix
    person(k,:) = [dt*k length(pos_x_hel) length(pos_y_infec) length(pos_x_rec)];
end


figure(2);
hold on;
xlabel('X direction Length')
ylabel('y direction Length')
legend('healthy', 'infected', 'recovered')
close(v);
figure(3);
h=area(person(:,1),person(:,2:4));  % Area is used for the area representation in graph
h(1).FaceColor=[0,1,0];    % RGB triplet for green
h(2).FaceColor=[1,0,0];   % RGB triplet for red
h(3).FaceColor=[0,0,1];    % RGB triplet for blue
axis([0,1,0,n]);         
xlabel('Time')
ylabel('Individuals')
legend('healthy', 'infected', 'recovered')
saveas(gcf,'pandamic_1.png')

