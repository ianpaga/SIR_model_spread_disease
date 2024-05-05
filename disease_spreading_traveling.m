
clear all;      
format long;    

% --- Input parameters ---
r0=4.5;                                                    %Epidemiological parameters
beta=1;
alpha=beta/r0;
%gamma = 0;
gamma= 1-1/r0;
h = 0.1;                                                   %Step size
Ntot = [88788,278523,124188];                              %Population of Lund, Malmö and Helsingborg     
S0 =[0.95*Ntot(1), 0.95*Ntot(2), 0.95*Ntot(3)];            %Initial conditions     
I0 =[0.05*Ntot(1), 0.05*Ntot(2), 0.05*Ntot(3)];
R0 =[0, 0, 0];                                            

N = 200;                                                   %Number of steps

w = [ 0,     0.3,   0.1;  %Probability matrix
      0.3,     0,   0.1;
      0.28, 0.18,   0];


% --- Initialize vectors ---  

tvec = zeros(N,1);   
S = zeros(N,3);    
I = zeros(N,3);
R = zeros(N,3);

S2h = zeros(N,3);     
I2h = zeros(N,3);
R2h = zeros(N,3);

Euler_errorS = zeros(N,3);
Euler_errorI = zeros(N,3);
Euler_errorR = zeros(N,3);

           
for p=1:1:3
    
    tvec(1)=0; 
    S(1,p)=S0(p);
    I(1,p)=I0(p);
    R(1,p)=R0(p);
    
    S2h(1,p)=S0(p);              %cuidado aquí. iguales. problema?
    I2h(1,p)=I0(p);
    R2h(1,p)=R0(p);
    
end

% --- The Euler algorithm ---

for n=1:1:N-1          % loop over different times 
   
    tvec(n+1) = tvec(n) + h;
    
    for p=1:1:3        % loop over different cities 
                       
        
       S(n+1,p) =  S(n,p) + (-beta*(I(n,p)/Ntot(p))-gamma)*S(n,p)*h + S(n,:)*w(p,:)'*h - sum(S(n,p)*w(:,p))*h;
       I(n+1,p) =  I(n,p) + (beta*(S(n,p)/Ntot(p))-alpha)*I(n,p)*h  + I(n,:)*w(p,:)'*h - sum(I(n,p)*w(:,p))*h ;
       R(n+1,p) =  R(n,p) + (gamma*S(n,p)+alpha*I(n,p))*h           + R(n,:)*w(p,:)'*h - sum(R(n,p)*w(:,p))*h;
       
       S2h(n+1,p) =  S(n,p) + (-beta*(I(n,p)/Ntot(p))-gamma)*S(n,p)*2*h + S(n,:)*w(p,:)'*2*h - sum(S(n,p)*w(:,p))*2*h;
       I2h(n+1,p) =  I(n,p) + (beta*(S(n,p)/Ntot(p))-alpha)*I(n,p)*2*h  + I(n,:)*w(p,:)'*2*h - sum(I(n,p)*w(:,p))*2*h;
       R2h(n+1,p) =  R(n,p) + (gamma*S(n,p)+alpha*I(n,p))*2*h           + R(n,:)*w(p,:)'*2*h - sum(R(n,p)*w(:,p))*2*h;
  
       Euler_errorS(n,p) = abs( S(n,p) - S2h(n,p) ); 
       Euler_errorI(n,p) = abs( I(n,p) - I2h(n,p) ); 
       Euler_errorR(n,p) = abs( R(n,p) - R2h(n,p) ); 
       
    end
end   

% --- Plots ---

figure(1)
title('Amount of susceptible S in Lund, Malmö and Helsingborg');
set(gca,'Fontsize',15); 
plot(tvec,S(:,1),'bo',tvec,S(:,2),'b-',tvec,S(:,3),'b+'); hold on;
plot(tvec,S(:,1),'b-',tvec,S(:,2),'b-',tvec,S(:,3),'b-'); hold on; 
legend('Lund','Malmö', 'Helsingborg');
xlabel('time (t)','Fontsize',15);                                  
ylabel('Amount of Suceptibles','Fontsize',15);
title(' Amount of Susceptible people in Lund, Malmö and Helsingborg');

figure(2)
title('Amount of Infected I in Lund, Malmö and Helsingborg');
set(gca,'Fontsize',15); 
plot(tvec,I(:,1),'ro',tvec,I(:,2),'r-',tvec,I(:,3),'r+'); hold on;
plot(tvec,I(:,1),'r-',tvec,I(:,2),'r-',tvec,I(:,3),'r-'); hold on; 
legend('Lund','Malmö', 'Helsingborg');
xlabel('time (t)','Fontsize',15);                                  
ylabel('Amount of Infected','Fontsize',15); 
title(' Amount of Infected people in Lund, Malmö and Helsingborg');

figure(3)
set(gca,'Fontsize',15); 
plot(tvec,R(:,1),'go',tvec,R(:,2),'g-',tvec,R(:,3),'g+'); hold on;
plot(tvec,R(:,1),'g-',tvec,R(:,2),'g-',tvec,R(:,3),'g-'); hold on; 
legend('Lund','Malmö', 'Helsingborg');
xlabel('time (t)','Fontsize',15);                                  
ylabel('Amount of Recovered','Fontsize',15);
title('Amount of Recovered people in Lund, Malmö and Helsingborg');


 figure(4)
 errorbar(tvec, S(:,1), Euler_errorS(:,1),'-s','MarkerSize',5, 'MarkerFaceColor', 'black','Color', 'b');
 hold on;
 errorbar(tvec, I(:,1), Euler_errorI(:,1),'-s','MarkerSize',5, 'MarkerFaceColor', 'black','Color', 'r');
 hold on;
 errorbar(tvec, R(:,1), Euler_errorR(:,1),'-s','MarkerSize',5, 'MarkerFaceColor', 'black','Color', 'g');
 title('Amount of susceptible, infected and recovered people in Lund');
 legend('S(t) susceptible','I(t) infected', 'R(t) recovered');
 xlabel('time (t)','Fontsize',15);               
 ylabel('Number of persons','Fontsize',15);
 hold off;
 
 figure(5)
 errorbar(tvec, S(:,2), Euler_errorS(:,1),'-s','MarkerSize',5, 'MarkerFaceColor', 'black','Color', 'b');
 hold on;
 errorbar(tvec, I(:,2), Euler_errorI(:,1),'-s','MarkerSize',5, 'MarkerFaceColor', 'black','Color', 'r');
 hold on;
 errorbar(tvec, R(:,2), Euler_errorR(:,1),'-s','MarkerSize',5, 'MarkerFaceColor', 'black','Color', 'g');
 title('Amount of susceptible, infected and recovered people in Malmö');
 legend('S(t) susceptible','I(t) infected', 'R(t) recovered');
 xlabel('time (t)','Fontsize',15);               
 ylabel('Number of persons','Fontsize',15);
 hold off;
 
 figure(6)
 errorbar(tvec, S(:,3), Euler_errorS(:,1),'-s','MarkerSize',5, 'MarkerFaceColor', 'black', 'Color', 'b');
 hold on;
 errorbar(tvec, I(:,3), Euler_errorI(:,1),'-s','MarkerSize',5, 'MarkerFaceColor', 'black', 'Color', 'r');
 hold on;
 errorbar(tvec, R(:,3), Euler_errorR(:,1),'-s','MarkerSize',5, 'MarkerFaceColor', 'black', 'Color', 'g');
 title('Amount of susceptible, infected and recovered people in Helsingborg');
 legend('S(t) susceptible','I(t) infected', 'R(t) recovered');
 xlabel('time (t)','Fontsize',15);               
 ylabel('Number of persons','Fontsize',15);
 hold off;
 
 




   
