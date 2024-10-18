%calculates CD in BA testing inoculation treatment 
%DHT procedure

function[time,tau_v,Q_beliefs]=inoc_CD_BA(MG_A,MG_B,MG_C,f1,f2,f3,N,tt)
%for 2 mega-nodes
%%%Define 
N=100; %nodes, 
T=150; %time periods
M=4; %number of hypotheses *Must be even integer*
D=5; %dim of parameters
s=5; %parameter for generating well-defined semi-def pos var mat
MG_A=0/N; %one meganode
MG_B=0/N;
MG_C=0/N;
f1=5/N;%conspirators
f2=0/N;%debunkers
f3=5/N;%prebunkers
tt=20;
V=100;
links=50;
eta=10^-5;


Q_R=zeros(V,T,M);
Q_v=zeros(T,M);
tau_v=zeros(T,1);



Q_beliefs=zeros(V,T);
B_beliefs=zeros(V,T);
CD=zeros(V,T);


for v=1:V

       %Weight matrix for BA 
         % Initialising adjacency matrix of initial core (fully connected)
        W = ones(links) - eye(links);
        W = sparse(W);
        
        % Loop on number of nodes to be added 
        for i = 1:N-links
        
            k = full(sum(W)); % Degree sequence
            p = k/sum(k); % Probability of linking to a node
        
            %%% Our strategy to simulate preferential attachement will be
            %%% to convert the vector p of probabilities into a vector of
            %%% cumulative probabilities (i.e., a vector that sums up to 1,
            %%% which is akin to the partition of the unit interval into k
            %%% sub-intervals, whose length is proportional to the
            %%% probability of connecting to nodes).
            
            p = cumsum(p); 
            r = rand(links,1); % m random numbers in [0,1]
        
            %%% Adding new nodes and links
            ind = [];
        
            %%% Loop on new links formed by new node: the links are formed
            %%% with the preexisting nodes whose corresponding
            %%% sub-interval in the vector p contains the random numbers in
            %%% r
            for d = 1:links
                        
                auxx = p - r(d);
                auxx = find(auxx > 0);
                ind = [ind; auxx(1)];
        
            end
        
            ind = unique(ind);
        
            %%% Creating new rows and columns in adjacency matrix
            W = [W; zeros(1,size(W,2))];
            W = [W zeros(size(W,1),1)];
        
            %%% Adding new links
            W(end,ind) = 1;
            W(ind,end) = 1;
        
        end 

        k = full(sum(W))'; % Degree sequence
        auxe = (k+1+eta).*ones(1,N);
        auxe = auxe.^-1; % the . here forces matlab to do operations element-wise
        W = W.*auxe;




          %%%calculate centrality of each node
         G=graph(W,'upper');
        %G=graph(W,'lower');
        C=centrality(G,'eigenvector'); %compute eig vec cent for each node
        [BB,ind_h]= maxk(C, round((f1+f3)*N)); %index of most influential nodes 

            
              %%%assign conditions to nodes
             ind_h1=ind_h(1:round((MG_A*N)));
         


             ind_h2=ind_h(1:round((MG_B*N)));

           
             ind_h3=randsample(ind_h,round((f1*N)));

             ind_h4=ind_h(1:round((f2*N)));

             
             ind_h5=setdiff(ind_h,ind_h3);




             ind_h6=ind_h(1:round((MG_C*N)));


             
            hackers1=[0.88,.0001,0.001,0.1189];
            hackers_m1=repmat(hackers1,round((MG_A)*N),1);
            hackers_m11=repmat(hackers1,round((f1)*N),1);
            
            hackers2=[0.1189,0.001,0.001,0.88];
            hackers_m2=repmat(hackers2,round((MG_B)*N),1);
            hackers_m22=repmat(hackers2,round((f2)*N),1);
           
            %mildly misinforming 
            rho=0.40;
            hr=rand(1,M-1);
           
            hr=hr/sum(hr)*(1-rho);
            hackers3=[rho,hr];
        
            hackers_m33=repmat(hackers3,round((f3)*N),1); % dist nodes
            hackers_m44=repmat(hackers3,round((MG_C)*N),1); %mega-node



         
      
        
        %initialize storage --an NxM matrix for each time period 
        %define as cell
        B_M=cell(T,1);
        Q_M=cell(T,1);
        
        %sigma, the covariance, is always the same, positive semi-def matrix
        sig=randn(D,s*D); % generate semi-def pos matrix
        sig=(sig*sig');
        
        %mu's
        mu_mat=cell(N,M/2);
        %initial beliefs
        b=zeros(N,M); %public beliefs
        q = rand(N,M); %private beliefs
        %normalize
        for i=1:N
          q(i,:)=q(i,:)/sum(q(i,:));
        end 
        q(ind_h1,:)=hackers_m1; %assign beliefs to meganode(s) 
        b(ind_h1,:)=hackers_m1;

        q(ind_h2,:)=hackers_m2; %assign beliefs to meganode(s) 
        b(ind_h2,:)=hackers_m2;

        q(ind_h3,:)=hackers_m11; %assign beliefs to conspirators
        b(ind_h3,:)=hackers_m11;

        q(ind_h4,:)=hackers_m22; %assign beliefs to debunkers
        b(ind_h4,:)=hackers_m22;

        q(ind_h5,:)=hackers_m33; %assign beliefs to debunkers
        b(ind_h5,:)=hackers_m33;

        q(ind_h6,:)=hackers_m44; %assign beliefs to cent prebunkers
        b(ind_h6,:)=hackers_m44;


      A=zeros(N,M);
         aux = [ones(1,M/2) 2*ones(1,M/2)];
         for i=1:N
          A(i,:) = aux(randperm(M));
         end
        
        
        
        
        
        for i=1:N
              for j=1:(M/2)
                  mu=randn(1,D);
                  mu_mat{i,j}=mu;
        
              end
        end  







        %%%%%%Inoculation Phase%%%%%




prob_t=zeros(N,1);
N_p=1:N;
mgs_t=zeros(N,M);
    
    
    
    
      for i=1:N
        array=[(MG_C*N)*hackers2;(MG_C*N)*hackers2];
        r_index=randi([1,size(array,1)]);
       
        mgs_t(i,:)=array(r_index,:);
      end 


      mgs_t(ind_h6,:)=hackers_m44;




  for t=1:tt





     for i=1:N
      prob_t(i)=dot(q(i,:),mgs_t(i,:)); %calculate probability of a node listening to mega-node 

     end 

     prob_t(ind_h6,1)=1;

      vu_Np=N_p(find(rand<prob_t)); %index of nodes that would listen to prebunking meganode in this iteration
      vu_Np=setdiff(vu_Np,ind_h6);



    obs=zeros(N,D); %observations
    pdf=zeros(N,(M/2)); %pdf
    
      for i=1:N
              %generate and save observation for each agent at each time 
              %each agent has M/2 possibilities
               
               obs(i,:)=mvnrnd(mu_mat{i,A(i,M)},sig);
        
               b_den = 0;% Denominator for public beliefs (Eq. 2 of the paper)
               aux = zeros(1,M); % Auxiliary vector for the numerators in Eq. 2
        
                for j=1:(M/2)
        
        
                    pdf(i,j)=mvnpdf(obs(i,:),mu_mat{i,j},sig);
        
                    ind = find(A(i,:) == j); % Finding indices associated to the hypotheses in the j-th subset
        
                    b_den = b_den + pdf(i,j)*sum(q(i,ind)); % Updating the denominator in Eq. 2 of the paper
                    
                    aux(ind) = pdf(i,j);
        
                   

                end
            
                  
                           
              %%%public beliefs%%%
        
              aux = aux.*q(i,:);


        
                    if any(i==ind_h6)
            
                        %aub=b; %hackers
                        b(i,:)=hackers3;


                    elseif any(i==vu_Np)
                    b(i,:)=hackers3; 


                   elseif any(i==ind_h5)
                    b(i,:)=hackers3;   
                  
                    else    
                
                        b(i,:) = aux/b_den; % M dimensional vector for agent i's public beliefs
                   
                    end

    
                
        
              
              
        
     end
            
           
        
            B_M{t}=b; %%% HERE you can compute the product matrix W*log(B), which you can use to compute the denominators in the q's
        
            
            q_den=(exp(W*(log(b))));
        
            %%%%private beliefs%%%%
        
            for i=1:N
                 for m=1:M
                    

                     if any(i==ind_h6)
            
                        q(i,:)=hackers3;  


                     elseif any(i==vu_Np)
                        q(i,:)=hackers3;  


                     elseif any(i==ind_h5)
                        q(i,:)=hackers3;
     
    
                     else  
                       q(i,m)=(exp(sum(((W(i,:)*log(b(:,m)))))))/(sum(q_den(i,:)));
                     end  
        %          
                   
                 end
        
            end
        
            
            
           
            Q_M{t}=q;




  end    







        %%%%%%%%%%%%%%%%%

        %prob=zeros(,1); %for each simulation v 
        N_L=1:N; %index of agents...subtract mega-nodes ***
        ind_m=[ind_h1,ind_h2];
        N_L=setdiff(N_L,ind_m);
      
        mgs=zeros(N,M);
        
        prob=zeros(N,1); %for each simulation v 


          for i=1:N
            array=[(MG_A*N)*hackers1;(MG_B*N)*hackers2];
            r_index=randi([1,size(array,1)]);
            
            mgs(i,:)=array(r_index,:);
          end 

          mgs(ind_h1,:)=hackers_m1;
          mgs(ind_h2,:)=hackers_m2; 

         


         A=zeros(N,M);
         aux = [ones(1,M/2) 2*ones(1,M/2)];
         for i=1:N
          A(i,:) = aux(randperm(M));
         end
        
        
        
        
        
        for i=1:N
              for j=1:(M/2)
                  mu=randn(1,D);
                  mu_mat{i,j}=mu;
        
              end
        end   
   
%        
         
        for t=1:T


            %initialize storage
            
            obs=zeros(N,D); %observations
            pdf=zeros(N,(M/2)); %pdf
            nodes=1:N;

% 
              for i=1:N
                prob(i,1)=dot(q(i,:),mgs(i,:)); %calculate probability of a node listening to mega-node 

              end 
                prob(ind_h1,1)=1;
                prob(ind_h2,1)=1;
                

              vu_N=nodes(find(rand<prob)); %index of nodes that would listen to one of the meganodes in this iteration

              vu_N=setdiff(vu_N,ind_h);

%              


           
         
            for i=1:N
              %generate and save observation for each agent at each time 
              %each agent has M/2 possibilities
               
               obs(i,:)=mvnrnd(mu_mat{i,A(i,M)},sig);
        
               b_den = 0;% Denominator for public beliefs (Eq. 2 of the paper)
               aux = zeros(1,M); % Auxiliary vector for the numerators in Eq. 2
        
                for j=1:(M/2)
        
        
                    pdf(i,j)=mvnpdf(obs(i,:),mu_mat{i,j},sig);
        
                    ind = find(A(i,:) == j); % Finding indices associated to the hypotheses in the j-th subset
        
                    b_den = b_den + pdf(i,j)*sum(q(i,ind)); % Updating the denominator in Eq. 2 of the paper
                    
                    aux(ind) = pdf(i,j);
        
                   

                end
            
                  
                           
              %%%public beliefs%%%
        
              aux = aux.*q(i,:);

%             
        
                    if any(i==ind_h1)
            
                       
                        b(i,:)=hackers1;
    
                    elseif any(i==vu_N) && isequal(mgs(i,:),hackers1)
                        b(i,:)=hackers1;
                    elseif any(i==vu_N) && isequal(mgs(i,:),hackers2)      
                        b(i,:)=hackers2;    
    
                    else    
                
                        b(i,:) = aux/b_den; % M dimensional vector for agent i's public beliefs
                   
                    end

%                  end
% 
%               end     
                
        
              
              
        
            end
            
           
        
            B_M{t}=b; %%% HERE you can compute the product matrix W*log(B), which you can use to compute the denominators in the q's
        %     
            
            q_den=(exp(W*(log(b))));
        
            %%%%private beliefs%%%%
        
            for i=1:N
                 for m=1:M
                     if any(i==ind_h1)
        
                       q(i,:)=hackers1;

                     elseif any(i==vu_N) && isequal(mgs(i,:),hackers1)
                        q(i,:)=hackers1;

                     elseif any(i==vu_N) && isequal(mgs(i,:),hackers2)
                        q(i,:)=hackers2;    
                     else  
                       q(i,m)=(exp(sum(((W(i,:)*log(b(:,m)))))))/(sum(q_den(i,:)));
                     end  
        %          
                   
                 end
        
            end
        
            
           
           
            Q_M{t}=q;
       
        
        end
        



          time=1:T;
          Q_E=zeros(T,N);
          B_E=zeros(T,N);
          tau=zeros(T,1); %storage for 'truthfulness' at each time for all players 
          tau_b=zeros(T,1); %storage for 'truthfulness' at each time for all players 
          aa=find(~cellfun(@isempty,Q_M));
          Q_M=Q_M(aa);
          B_M=B_M(aa);

     for t=1:T
            for i=1:N
     
            Q_E(t,i)=Q_M{t}(i,M); %belief in true hyp at each time for each N
            B_E(t,i)=B_M{t}(i,M);
            
            end    
           
           
    end
   
      length_h=length(ind_h);
      %%%%%% %%%
      au = Q_E;
      au(:,ind_h) = [];
      aub=B_E;
      aub(:,ind_h)=[];


      %%%
        for t=1:T
            
             tau(t,1)=sum(au(t,:))/(N-length_h);
             tau_b(t,1)=sum(aub(t,:))/(N-length_h);
            
        end    

     
        Q_beliefs(v,:)=tau'; 
        B_beliefs(v,:)=tau_b';




%        
end

   

    
    
    
      for t=1:T
       tau_v(t,1)=sum(abs((Q_beliefs(:,t)-B_beliefs(:,t))))/V; %CD
      end