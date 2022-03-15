clear all;

%% Take Home Part e
% 
 J =1;
 L =22;
 %      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
Adj = [ 0 1 1 1 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0                %1 
        1 0 0 0 1 1 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %2 
        1 0 0 0 0 0 1 1 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %3 
        1 0 0 0 0 0 0 0 1 1  0  0  0  0  0  0  0  0  0  0  0  0                 %4 
        0 1 0 0 0 0 0 0 0 0  1  1  0  0  0  0  0  0  0  0  0  0                 %5 
        0 1 0 0 0 0 0 0 0 0  0  0  1  1  0  0  0  0  0  0  0  0                 %6 
        0 0 1 0 0 0 0 0 0 0  0  0  0  0  1  1  0  0  0  0  0  0                 %7 
        0 0 1 0 0 0 0 0 0 0  0  0  0  0  0  0  1  1  0  0  0  0                 %8
        0 0 0 1 0 0 0 0 0 0  0  0  0  0  0  0  0  0  1  1  0  0                 %9 
        0 0 0 1 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  1  1                 %10 
        0 0 0 0 1 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %11
        0 0 0 0 1 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %12
        0 0 0 0 0 1 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %13
        0 0 0 0 0 1 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %14
        0 0 0 0 0 0 1 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %15
        0 0 0 0 0 0 1 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %16
        0 0 0 0 0 0 0 1 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %17
        0 0 0 0 0 0 0 1 0 0  0  0  0  0  0  0  0  0  0  0  0  0                 %18
        0 0 0 0 0 0 0 0 1 0  0  0  0  0  0  0  0  0  0  0  0  0                 %19
        0 0 0 0 0 0 0 0 1 0  0  0  0  0  0  0  0  0  0  0  0  0                 %20
        0 0 0 0 0 0 0 0 0 1  0  0  0  0  0  0  0  0  0  0  0  0                 %21
        0 0 0 0 0 0 0 0 0 1  0  0  0  0  0  0  0  0  0  0  0  0 ];              %22     
        
        
figure(1); imagesc(Adj)
title('Adjacency Matrix');

figure(2);
plot(graph(Adj));

allowed_spins = [-1 1];
spin_matrix = zeros(L,1);

% Random spin matrix initialization
for i = 1:1:L
        
        spin = randsample(allowed_spins,1,true);
        spin_matrix(i) = spin;

end

disp('Initial spin matrix is: ')
disp(spin_matrix);

figure(3);
imagesc(spin_matrix);

% Calculate initial energy
E =0;

NT =1e5;
N = 1:1:L ;
i = 0;

for beta=0:0.05:2   
    
    i= i +1;
      % inverse temp
    for t= 1:1:NT


        % Select a random spin to flip
        %select_spin_number = 10;
        select_spin_number = randsample(N,1,true);
          h = 0;
          for k = 1:1:L       % this loop search for neighbors in adjacency matrix
                    if ( Adj(k,select_spin_number) ==1 )  
                        h = h + spin_matrix(k,1);
                    end        
          end

          del_E = 2*h*spin_matrix(select_spin_number) ;        % Change in energy expected on flip of selected spin
          gamma = exp(-1*beta*del_E);

          random_number = rand;
          if (random_number < gamma)
              spin_matrix(select_spin_number)= -1*spin_matrix(select_spin_number) ;
              %Net_energy = Net_energy + del_E ;
          end

%           if (t ==1)
%               figure();
%               imagesc(spin_matrix);
%           end
% 
% 
%           if (t ==1000)
%               figure();
%               imagesc(spin_matrix);
%           end
% 
% 
%           if (t ==1e5)
%               figure();
%               imagesc(spin_matrix);
%           end

    end
    mean_spins(i) = mean(spin_matrix);
end

beta_plot=0:0.05:2  ;
figure(4);
plot(beta_plot, mean_spins);
xlabel('Beta');
ylabel('Mean of spins');


%% Take Home Part f

% J =1;
% L =16;
%  %      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 
% Adj = [ 0 1 1 1 0 0 0 0 0 0  0  0  0  0  0  0                   %1 
%         1 0 0 0 1 1 0 0 0 0  0  0  0  0  0  0                   %2 
%         1 0 0 0 0 0 1 1 0 0  0  0  0  0  0  0                   %3 
%         1 0 0 0 0 0 0 0 1 1  0  0  0  0  0  0                   %4 
%         0 1 0 0 0 0 0 0 0 0  1  0  0  0  0  0                   %5 
%         0 1 0 0 0 0 0 0 0 0  0  1  0  0  0  0                   %6 
%         0 0 1 0 0 0 0 0 0 0  0  0  1  0  0  0                   %7 
%         0 0 1 0 0 0 0 0 0 0  0  0  0  1  0  0                   %8
%         0 0 0 1 0 0 0 0 0 0  0  0  0  0  1  0                   %9 
%         0 0 0 1 0 0 0 0 0 0  0  0  0  0  0  1                   %10 
%         0 0 0 0 1 0 0 0 0 0  0  0  0  0  0  0                   %11
%         0 0 0 0 0 1 0 0 0 0  0  0  0  0  0  0                   %12
%         0 0 0 0 0 0 1 0 0 0  0  0  0  0  0  0                   %13
%         0 0 0 0 0 0 0 1 0 0  0  0  0  0  0  0                   %14
%         0 0 0 0 0 0 0 0 1 0  0  0  0  0  0  0                   %15
%         0 0 0 0 0 0 0 0 0 1  0  0  0  0  0  0       ];          %16
%         
%         
% figure(1); imagesc(Adj)
% title('Adjacency Matrix');
% 
% figure(2);
% plot(graph(Adj));
% 
% allowed_spins = [-1 1];
% spin_matrix = zeros(L,1);
% 
% % Random spin matrix initialization
% for i = 1:1:L
%         
%         spin = randsample(allowed_spins,1,true);
%         spin_matrix(i) = spin;
% 
% end
% 
% spin_matrix(11) = 1;
% spin_matrix(12) = 1;
% spin_matrix(13) = 1;
% spin_matrix(14) = 1;
% spin_matrix(15) = 1;
% spin_matrix(16) = 1;
% 
% disp('Initial spin matrix is: ')
% disp(spin_matrix);
% 
% figure(3);
% imagesc(spin_matrix);
% 
% % Calculate initial energy
% E =0;
% NT =1e5;
% N = 1:1:L ;
% i = 0;
% 
% for beta=2   
%     i= i + 1;
%       % inverse temp
%     for t= 1:1:NT
%         % Select a random spin to flip
%         %select_spin_number = 10;
%         select_spin_number = randsample(N,1,true);
%           h = 0;
%           if(select_spin_number <= 16)
%               
%               for k = 1:1:L       % this loop search for neighbors in adjacency matrix
%                         if ( Adj(k,select_spin_number) ==1 )  
%                             h = h + spin_matrix(k,1);
%                         end        
%               end
% 
%               del_E = 2*h*spin_matrix(select_spin_number) ;        % Change in energy expected on flip of selected spin
%               gamma = exp(-1*beta*del_E);
% 
%               random_number = rand;
%               if (random_number < gamma)
%                   spin_matrix(select_spin_number)= -1*spin_matrix(select_spin_number) ;
%                   % Net_energy = Net_energy + del_E ;
%               end
%   
%           end
%            
%               if (t ==1)
%                   figure();
%                   imagesc(spin_matrix);
%               end
% 
% 
%               if (t ==1000)
%                   figure();
%                   imagesc(spin_matrix);
%               end
% 
% 
%               if (t ==1e5)
%                   figure();
%                   imagesc(spin_matrix);
%               end
%          
%     end
%     mean_spins(i) = mean(spin_matrix);
% end
% 
% % beta_plot=0:0.05:2  ;
% % figure(4);
% % plot(beta_plot, mean_spins);
% % xlabel('Beta');
% % ylabel('Mean of spins');
% % title('Leaf Nodes clamped to -1');
% 

%% Convert mm to system state S

% Boltzmann Law  
beta = 1;
Z0 = 2;
allowed_spins = [-1 +1];
count  = 0
E = zeros(2^4,1);
for i = [-1 1]
    for j = [-1 1]
        for k = [-1 1]
            for l =[-1 1]
                count = count +1 ;
                E(count) = -1*i*(j+k+l) ;
            end
        end
    end
end
                

P=exp(-E*beta);

Formulated_value = Z0*(2*cosh(beta))^(3) ;

%  Z = ... Given all probabilities in P , calculate the partition function

Z = sum (P,'all');
