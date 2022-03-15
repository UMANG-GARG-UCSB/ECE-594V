clear all;
clc;


%% Problem 1
  
% Ntrials = 10000;
% Nruns = 1;
% delta = 0.25;
% 
% for k =1:1:length(delta)
%     
%     x = zeros(1, Ntrials);
%     y = zeros(1, Ntrials);
% 
%     %   x(1) = -1 + 2*rand();
%     %   y(1) = -1 + 2*rand();
% 
%     for j=1:1:Nruns
% 
%           x(1) = 1;
%           y(1) = 1;
%           Nhit(1) = 0;
%           accept_moves(j) = 0;
%           reject_moves(j) = 0;
% 
%         for i=1:1:Ntrials
% 
%             delx = -delta(k) + 2*delta(k)*rand();
%             dely = -delta(k) + 2*delta(k)*rand();
%             N(i) = i;
% 
%                  x(i+1) = x(i);
%                  y(i+1) = y(i);
% 
%             if  (abs(x(i)+ delx) < 1) && (abs(y(i)+ dely) < 1) 
%                  x(i+1) = x(i)+ delx;
%                  y(i+1) = y(i)+ dely;
%                  accept_moves(j) = accept_moves(j) + 1;
% 
%             else     
%                  reject_moves(j) = reject_moves(j) + 1;    
%             end
% 
%             if ((x(i+1))^2 + (y(i+1))^2 ) < 1
%                 Nhit(i) = Nhit(i) + 1 ;
%             end
%                 Nhit(i+1) = Nhit(i);
% 
%                 temp(i) = 4*Nhit(i)/N(i);
% 
%         end
% 
%     Run(j) = Nhit(Ntrials);
%     AR(j) = accept_moves(j)/ N ;
% 
%     end    
% 
% avg_AR(k) = mean(AR)
% end
% 
% 
% M2= zeros(1,Ntrials);
% M2(:)= pi;
% 
% Nhit = Nhit(1:Ntrials);
% figure(1);% Will give trials for only the last run
% p1= plot(N, temp);
% ylabel('Nhits / Ntrials');
% xlabel(' Trials ' );
% p1. LineWidth = 3;
% hold on;
% p2= plot(N, M2,'-');
% p2. LineWidth = 3;
% 
% figure(2);
% p3 = plot(delta, avg_AR);
% ylabel('Average acceptance ratio (AR)');
% xlabel(' \delta ' );
% p3. LineWidth = 3;
% 
% var(avg_AR);
% sigma = sqrt(var(avg_AR));
% 
% 


%% Problem 1 part e

% Ntrials = 10000;
% Nruns = 100;
% delta = [0.1, 0.5, 1, 1.5, 2]
% 
% for k =1:1:length(delta)
%     
%     x = zeros(1, Ntrials);
%     y = zeros(1, Ntrials);
% 
%     for j=1:1:Nruns
% 
%           x(1) = 1;
%           y(1) = 1;
%           Nhit(1) = 0;
%           accept_moves(j) = 0;
%           reject_moves(j) = 0;
% 
%         for i=1:1:Ntrials
% 
%             delx = -delta(k) + 2*delta(k)*rand();
%             dely = -delta(k) + 2*delta(k)*rand();
%             N(i) = i;
% 
%                  x(i+1) = x(i);
%                  y(i+1) = y(i);
% 
%             if  (abs(x(i)+ delx) < 1) && (abs(y(i)+ dely) < 1) 
%                  x(i+1) = x(i)+ delx;
%                  y(i+1) = y(i)+ dely;
%                  accept_moves(j) = accept_moves(j) + 1;
%                 
%             else     
%                  reject_moves(j) = reject_moves(j) + 1;    
%             end
% 
%             if ((x(i+1))^2 + (y(i+1))^2 ) < 1
%                 Nhit(i) = Nhit(i) + 1 ;
%             end
%                 Nhit(i+1) = Nhit(i);
% 
%                 temp(i) = 4*Nhit(i)/N(i);
% 
%         end
%         
%     AR(j) = accept_moves(j)/ Ntrials ;      % Acceptance ratio is accept/ Total
%     err_in_run(j) = (4*Nhit(Ntrials) - floor(Ntrials * pi)); 
%     Run(j) = Nhit(Ntrials);
%     
%     end    
% 
% avg_AR(k) = mean(AR)            
% sigma(k) = sqrt(var(err_in_run))
% end
% 
% 
% M2= zeros(1,Ntrials);
% M2(:)= pi;
% 
% Nhit = Nhit(1:Ntrials);
% figure(1);% Will give trials for only the last run
% p1= plot(N, temp);
% ylabel('Nhits / Ntrials');
% xlabel(' Trials ' );
% p1. LineWidth = 3;
% hold on;
% p2= plot(N, M2,'-');
% p2. LineWidth = 3;
% 
% figure(2);
% p3 = plot(delta, avg_AR*100);
% ylabel('Average acceptance ratio (AR)');
% xlabel(' \delta ' );
% hold on;
% 
% p3. LineWidth = 3;
% 
% figure(3);
% p4 = plot(delta, sigma);
% ylabel('Standard deviation in error');
% xlabel(' \delta ' );
% hold on;
% 
% p4. LineWidth = 3;


%% Problem 2

% finding Cdf for normal distribution -- taken from Wikipedia

% x = -10:0.01:10 ;
% y = (1 + erf(x/sqrt(2)) )/2 ; 
% 
% figure(5);
% p5 = plot (x,y);
% ylabel('Cdf function');
% xlabel('x' ); 
% p5.LineWidth = 3;

% Nsamples = 100000;
% 
% for i=1:1:Nsamples
% 
%  y(i)= rand();
%  x(i) = (sqrt(2)) * erfinv(2*y(i) -1) ;
%  
% end
% 
% h = histfit(x);
% avg = mean(x)
% sigma = sqrt(var(x))

%% Problem 3 

% N = 100000;
% 
% x = -1:2/N:1;
% f = (1/(N+1))*ones(1,N+1);
% avg = 0;
% expectation_x_squared  =0;
% for i =1:1:N+1
%     
%     avg = avg + (x(i)*f(i)) ;
%     expectation_x_squared = expectation_x_squared + f(i)*(x(i))^2 ;
% end
% 
% % Using Var(X) = E(X2) – [E(X)]2  
% 
% 
% variance = expectation_x_squared - avg^2 ;
% SD = sqrt(variance);


%% Problem 3 part b

clc; 
clear all;

N = 100;
a=-1;b =1;

x = -1:2/N:1;
f = (1/(N+1))*ones(1,N+1);
avg = 0;
expectation_x_squared  = 0;
% for i =1:1:N+1
%     
%     avg = avg + (x(i)*f(i)) ;
%     expectation_x_squared = expectation_x_squared + f(i)*(x(i))^2 ;
% end
% 
% % Using Var(X) = E(X2) – [E(X)]2  
% 
% 
% variance = expectation_x_squared - avg^2 ;
% SD = sqrt(variance);

% N_distributions = 100;
% N_samples = 1e5;
% 
% for j= 1:1: N_distributions    
%     
%         spectrum(j,:) = a + (b-a)*rand(1,N_samples) ;   
% end  
% 
% for i = 1:N_samples
%     spectrum_final(i) = sum(spectrum(:,i));
% end
% 
% % Normalize the spectrum
%     spectrum_final = spectrum_final/N_distributions;
% 
%     % Theoritical and Actual comparison
%     
%     variance = var(spectrum(1,:));
%     sd = sqrt(var(spectrum(1,:)));
%     theory_sd = sd*(1/sqrt(N_distributions));
%     theory_var = (1/N_distributions)*variance;
%     
%     sd_Actual = sqrt(var(spectrum_final));
%     var_Actual = var(spectrum_final);
%     
%     figure();
%     histogram(spectrum_final,'Normalization','probability');
%     xlabel('Sample Sum'); 
%     ylabel('Probability pdf');
%   
%     
