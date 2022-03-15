clear all;

%% Problem 

% mu = 0;
% sigma = 1;
% Ii = -5:0.1:5;
% total_time = 10000;
% time = 1:1:total_time ;
% 
% for i = 1:1:length(Ii)
%     for t = 1:1:total_time
% 
%         r = normrnd(mu, sigma);
%         temp = Ii(i) + r;
% 
%         if (temp<0)
%             s(i,t) = 0;
%         else
%             s(i,t) = 1;
%         end
%     end
%     avg_s(i) = sum(s(i,:))/total_time ;
% end
% 
% % figure(1);
% % p1 = plot(time,s);
% % xlabel('Time');
% % ylabel('s_{i}');
% % title('Output for I_{i} = -1');
% 
% figure(2);
% p2 = plot(Ii, avg_s);
% xlabel('I_{i}');
% ylabel('<s_{i}>');
% hold on; 
% plot(Ii, (1+tanh(Ii))/2);
% legend('<s_{i}>', 'Scaled Tanh(I_{i})');
% 
% %% Problem 1 part c
% 
% x = [-10:.1:10];
% gaussian_distribution = normpdf(x,0,1);
% figure();
% plot(x,gaussian_distribution);
% 
% for i = 1:1:length(Ii)
% 
%     r = normrnd(mu, sigma);
%     Arg(i) = r + Ii(i) ;
% end
%  
% figure();
% plot(Ii, Arg);
% xlabel('I_{i}');
% ylabel('Argument');
% 
% 

%% Problem 2

L = 30;
x = L;
y = L;
adj = zeros(x*y);

for i=1:x
    for j=1:y
        k = sub2ind([x y],i,j);
        if i>1
            ii=i-1; jj=j;
            adj(k,sub2ind([x y],ii,jj)) = 1; 
        end
        if i<x
            ii=i+1; jj=j;
            adj(k,sub2ind([x y],ii,jj)) = 1; 
        end
        if j>1
            ii=i; jj=j-1;
            adj(k,sub2ind([x y],ii,jj)) = 1; 
        end
        if j<y
            ii=i; jj=j+1;
            adj(k,sub2ind([x y],ii,jj)) = 1; 
        end
    end
end

figure(1); imagesc(adj)
title('Adjacency Matrix');

%% Problem 2 part b 

% allowed_spins = [-1 1];
% spin_matrix = zeros(L);
% beta = 1;
% NT = 1;
% J = 1;
% 
% % Random spin matrix initialization
% for i = 1:1:L
%     for j =1:1:L
%         
%         spin = randsample(allowed_spins,1,true);
%         spin_matrix(i,j) = spin;
%  
%     end
% end
% 
% disp('Initial spin matrix is: ')
% disp(spin_matrix);
% 
% figure();
% imagesc(spin_matrix);
% 
% % Calculate initial energy
% E =0;
% 
% for i = 1:1:L               %% First two loops select one spin in the matrix
%     for j = 1:1:L           %%
%         
%         index = sub2ind([x y],i,j);
%         nbr_sum = 0;
%         
%         for k = 1:1:L       % these 2 loops search for neighbors in adjacency matrix
%             for l = 1:1:L
%                 if ( adj(index, sub2ind([x y],k,l)) ==1 )  
%                     nbr_sum = nbr_sum + spin_matrix(k,l);
%                 end     
%             end
%         end
% 
%             E = E + (-1*J*spin_matrix(i,j)*nbr_sum);
%         
%     end
% end
% 
% Net_energy = E/2;           % Initial energy
% disp('Initial Total energy')
% disp(Net_energy);
% 
% NT =1e5;
% N = 1:1:100 ;    
% for t= 1:1:NT
%     
%     beta = 1*beta;
%     % Select a random spin to flip
%     %select_spin_number = 10;
%     select_spin_number = randsample(N,1,true);
%     [row,col] = ind2sub([x y],select_spin_number);
%       h = 0;
%       for k = 1:1:L       % these 2 loops search for neighbors in adjacency matrix
%             for l = 1:1:L
%                 if ( adj(sub2ind([x y],k,l),select_spin_number) ==1 )  
%                     h = h + spin_matrix(k,l);
%                 end     
%             end
%       end
%         
%       del_E = 2*h*spin_matrix(row,col) ;        % Change in energy expected on flip of selected spin
%       gamma = exp(-1*beta*del_E);
%       
%       random_number = rand;
%       if (random_number < gamma)
%           spin_matrix(row,col)= -1*spin_matrix(row,col) ;
%           Net_energy = Net_energy + del_E ;
%       end
%       
%       if (t ==1)
%           figure();
%           imagesc(spin_matrix);
%       end
%       
%             
%       if (t ==1000)
%           figure();
%           imagesc(spin_matrix);
%       end
%       
%             
%       if (t ==1e5)
%           figure();
%           imagesc(spin_matrix);
%       end
%    
% end
% 
% disp('Timesteps'); 
% disp(NT);
% disp('Updated Total energy')
% disp(Net_energy);
% disp('Updated spin matrix is: ')
% disp(spin_matrix);

%% Problem 2 part c

allowed_spins = [-1 1];
spin_matrix = zeros(L);
beta = 0.01;
J = 1;

% Random spin matrix initialization
for i = 1:1:L
    for j =1:1:L
        
        spin = randsample(allowed_spins,1,true);
        spin_matrix(i,j) = spin;
 
    end
end

disp('Initial spin matrix is: ')
disp(spin_matrix);

figure();
imagesc(spin_matrix);

% Calculate initial energy
E =0;

for i = 1:1:L               %% First two loops select one spin in the matrix
    for j = 1:1:L           %%
        
        index = sub2ind([x y],i,j);
        nbr_sum = 0;
        
        for k = 1:1:L       % these 2 loops search for neighbors in adjacency matrix
            for l = 1:1:L
                if ( adj(index, sub2ind([x y],k,l)) ==1 )  
                    nbr_sum = nbr_sum + spin_matrix(k,l);
                end     
            end
        end

            E = E + (-1*J*spin_matrix(i,j)*nbr_sum);
        
    end
end

Net_energy = E/2;           % Initial energy
disp('Initial Total energy')
disp(Net_energy);

NT =1e5;
N = 1:1:L*L ;    
for t= 1:1:NT
    
    
    % Select a random spin to flip
    %select_spin_number = 10;
    select_spin_number = randsample(N,1,true);
    [row,col] = ind2sub([x y],select_spin_number);
      h = 0;
      for k = 1:1:L       % these 2 loops search for neighbors in adjacency matrix
            for l = 1:1:L
                if ( adj(sub2ind([x y],k,l),select_spin_number) ==1 )  
                    h = h + spin_matrix(k,l);
                end     
            end
      end
        
      del_E = 2*h*spin_matrix(row,col) ;        % Change in energy expected on flip of selected spin
      gamma = exp(-1*beta*del_E);
      
      random_number = rand;
      if (random_number < gamma)
          spin_matrix(row,col)= -1*spin_matrix(row,col) ;
          Net_energy = Net_energy + del_E ;
      end
      
      if (t ==1)
          figure();
          imagesc(spin_matrix);
      end
      
            
      if (t ==1000)
          figure();
          imagesc(spin_matrix);
      end
      
            
      if (t ==1e5)
          figure();
          imagesc(spin_matrix);
      end
   
end

disp('Timesteps'); 
disp(NT);
disp('Updated Total energy')
disp(Net_energy);
disp('Updated spin matrix is: ')
disp(spin_matrix);



