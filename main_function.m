 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Построение векторной авторегрессии с семью переменными (6 + 1 экзо).
 % последняя итоговая версия        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d1=xlsread('r.xlsx');
d2=xlsread('cpi.xlsx');
d3=xlsread('y_min.xlsx');
d4=xlsread('y.xlsx');
d5=xlsread('d.xlsx'); %dlog total
d6=xlsread('d1.xlsx'); %dlog old
d7=xlsread('param.xlsx'); %dlog old

data=[d1,d2,d3,d4,d5,d6];
exo=xlsread('p.xlsx'); 

N=size(data,2); %Число эндогенных переменных
m=2; % Число экзогенных переменных (константа + p_oil)
L=12;   %number of lags in the VAR
Y=data;
k=N*L+m;   % Число оцениваемых коэффциентов в одном уравнении
y_init=data';
y_init=y_init(:, 1:L);
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
%X=[X ones(size(Y,1),1)];
X=[X ones(size(Y,1),1) exo];
Y=Y(L+1:end,:);
X=X(L+1:end,:);
T=rows(X);
irf_length=T+L;

%% Compute standard deviation of each series residual via an ols regression to be used in setting the prior
%first variable
y=Y(:,1);
x=X(:,[N*L+1 1]); % Из X берём первую переменную (то есть лаг первой переменной) и константу, которая по номеру идёт N*L+1
b0=inv(x'*x)*(x'*y);
s1=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));  %std of residual standard error
%second variable
y=Y(:,2);
x=X(:,[N*L+1 2]); 
b0=inv(x'*x)*(x'*y);
s2=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));
%third variable
y=Y(:,3);
x=X(:,[N*L+1 3]); 
b0=inv(x'*x)*(x'*y);
s3=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));
%fourth variable
y=Y(:,4);
x=X(:,[N*L+1 4]); 
b0=inv(x'*x)*(x'*y);
s4=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));
%fifth variable
y=Y(:,5);
x=X(:,[N*L+1 5]); 
b0=inv(x'*x)*(x'*y);
s5=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));
%sixth variable
y=Y(:,6);
x=X(:,[N*L+1 6]); 
b0=inv(x'*x)*(x'*y);
s6=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));

%% Parameters to control the prior
lamda1=0.1;  %tightness prior on the AR coefficients
lamda2=0.5;
lamda3=2;   %tightness of prior on higher lags 
lamda4=100;  %tightness of prior on the constant term

%specify the prior mean of the coefficients of the Two equations of the VAR

B0=zeros((N*L+m),N);
for i=1:N
    B0(i,i)=1;
end
% Пока для экзогенных переменных и для константы предпологаем нулевое
% среднее
B0(N*L+1,1)=0; 
B0(N*L+1,2)=0;
B0=vec(B0);

%% Определим априорную ковариоционную матрицу для параметров BVAR
 % Код для ковариоционной матрицы распредления Миннесоты взят из кодов
 % BEAR Toolbox, разработанным ЕЦБ. Номера уравнений из Technical guide
 % пакетов
H=zeros(k*N,k*N);
arvar=[s1, s2, s3, s4, s5, s6]';
% set the variance on coefficients trelated to own lags, using (1.3.5)
for ii=1:N
   for jj=1:L
   H((ii-1)*k+(jj-1)*N+ii,(ii-1)*k+(jj-1)*N+ii)=(lamda1/jj^lamda3)^2;
   end
end


%  set variance for coefficients on cross lags, using (1.3.6)
for ii=1:N
   for jj=1:L
      for kk=1:N
      if kk==ii
      else
      H((ii-1)*k+(jj-1)*N+kk,(ii-1)*k+(jj-1)*N+kk)=(arvar(ii,1)/arvar(kk,1))^2*(((lamda1*lamda2)/(jj^lamda3))^2);
      end
      end
   end
end

% finally set the variance for exogenous variables, using (1.3.7)
for ii=1:N 
   for jj=1:m
   H(ii*k-m+jj,ii*k-m+jj)=(arvar(ii,1))^2*((lamda1*lamda4)^2);
   end
end

%prior scale matrix for sigma the VAR covariance
S=eye(N);
%prior degrees of freedom
alpha=N+1;

%starting values for the Gibbs sampling algorithm
Sigma=eye(N);
betaols=vec(inv(X'*X)*(X'*Y));
REPS=2000;
BURN=1000;

out1=zeros(REPS-BURN,irf_length,N); % Для irf на шок дкп
out01=zeros(1001,irf_length,1); % Для irf на шок дкп
%out2=zeros(1001,irf_length,1); % Для irf на шок дкп (old debt)
%out3=zeros(1001,irf_length,1); % Для irf на шок дкп (new credit)
%out4=zeros(1001,irf_length,1); % Для irf на шок дкп (total credit)

out8=zeros(REPS-BURN,1,N*k); % Для хранения сэмплов коэффициентов beta;
out9=zeros(REPS-BURN,N,N); % Для хранения сэмплов коэффициентов дисперсий;

out10=zeros(1001,irf_length,1); % Для irf долговой нагрузки момент
out010=zeros(1001,irf_length,1); % Для irf долговой нагрузки накоп

fevd_record1=zeros(REPS-BURN,144,N); % Для FEVD вклад шока дкп
hd_record1=zeros(REPS-BURN,irf_length-L,N); % Для HD вклад шока дкп

likelihood=[];

%% Гиббс сэмплинг
jj=1;
for i=1:REPS
M=inv(inv(H)+kron(inv(Sigma),X'*X))*(inv(H)*B0+kron(inv(Sigma),X'*X)*betaols);
V=inv(inv(H)+kron(inv(Sigma),X'*X));
%check for stability of the VAR
check=-1;
while check<0 
beta=M+(randn(1,N*(N*L+m))*chol(V))';
CH=stability(beta,N,L,m);
if CH==0
    check=10;
end
end
%draw sigma from the IW distribution
e=Y-X*reshape(beta,N*L+m,N);
%scale matrix
scale=e'*e+S;
Sigma=IWPQ(T+alpha,inv(scale));
    
    if i>=BURN
        %impose sign restrictions
        chck=-1;
        while chck<0
            K=randn(N,N);
            Q=getqr(K);
            A0hat=chol(cov(e));
            A0hat1=(A0hat*Q);  %candidate draw
            %check signs
            e1=A0hat1(1,1)>0;  %отклик r на шок дкп
            e2=A0hat1(1,2)<0;  %отклик cpi на шок дкп
            e3=A0hat1(1,3)<0;  %отклик y_min на шок дкп
            e4=A0hat1(1,4)<0;  %отклик y на шок дкп
            e5=A0hat1(1,5)<0;  %отклик d_min на шок дкп
            %e6=A0hat1(1,6)>=0;  %отклик d_min1 на шок дкп

            %e7=A0hat1(2,2)<0;  %отклик r на шок предложения
            e8=A0hat1(2,2)<0;  %отклик cpi на шок предложения
            %e9=A0hat1(2,3)>0;  %отклик y_min на шок предложения
            e10=A0hat1(2,4)>0;  %отклик y на шок предложения
            %e11=A0hat1(2,5)?0;  %отклик d_min на шок предложения

            e13=A0hat1(3,1)>0;  %отклик r на шок спроса1
            %e14=A0hat1(3,2)>0;  %отклик cpi на шок спроса1
            e15=A0hat1(3,3)>0;  %отклик y_min на шок спроса1
            %e16=A0hat1(3,4)>0;  %отклик y на шок спроса1
            %e17=A0hat1(3,5)>0;  %отклик d_min на шок спроса1

            e19=A0hat1(4,1)>0;  %отклик r на шок спроса
            %e20=A0hat1(4,2)>0;  %отклик cpi на шок спроса
            %e21=A0hat1(4,3)>0;  %отклик y_min на шок спроса
            e22=A0hat1(4,4)>0;  %отклик y на шок спроса
            %e23=A0hat1(4,5)>0;  %отклик d_min на шок спроса

            e25=A0hat1(5,1)>0;  %отклик r на шок кредитования
            %e26=A0hat1(5,2)=0;  %отклик cpi на шок кредитования
            %e27=A0hat1(5,3)?0;  %отклик y_min на шок кредитования
            %e28=A0hat1(5,4)>0;  %отклик y на шок кредитования
            e29=A0hat1(5,5)>0;  %отклик d_min на шок кредитования

            e31=A0hat1(6,1)>0;  %отклик r на шок кредитования1
            e32=A0hat1(6,2)>0;  %отклик cpi на шок кредитования1
            %e33=A0hat1(5,3)?0;  %отклик y_min на шок кредитования1
            %e34=A0hat1(5,4)>0;  %отклик y на шок кредитования1
            %e36=A0hat1(6,6)>0;  %отклик d_min на шок кредитования1

            if e1+e2+e3+e4+e5+e8+e10+e13+e15+e19+e22+e25+e29+e31+e32==15
                chck=10;
            else
                e1=-A0hat1(1,1)>0;  
                e3=-A0hat1(1,3)<0; 
                e2=-A0hat1(1,2)<0;
                e4=-A0hat1(1,4)<0; 
                e5=-A0hat1(1,5)<0;
                %e6=-A0hat1(1,6)>=0;
          
                e8=-A0hat1(2,2)<0;
                e10=-A0hat1(2,4)>0;  
            
                e13=-A0hat1(3,1)>0;  
                %e12=-A0hat1(3,2)>0;  
                e15=-A0hat1(3,3)>0; 
                %e14=-A0hat1(3,4)>0;  
                %e15=-A0hat1(3,5)>0; 

                e19=-A0hat1(4,1)>0; 
                e22=-A0hat1(4,4)>0;
            
                e25=-A0hat1(5,1)>=0; 
                e29=-A0hat1(5,5)>0;

                e31=-A0hat1(6,1)>=0; 
                e32=-A0hat1(6,2)>0;
                
                if e1+e2+e3+e4+e5+e8+e10+e13+e15+e19+e22+e25+e29+e31+e32==15
                    A0hat1(1,1)=-A0hat1(1,1);
                    A0hat1(1,2)=-A0hat1(1,2);
                    A0hat1(1,3)=-A0hat1(1,3); 
                    A0hat1(1,4)=-A0hat1(1,4);
                    A0hat1(1,5)=-A0hat1(1,5);
                    %A0hat1(1,6)=-A0hat1(1,6);
          
                    %A0hat1(1,2)=-A0hat1(1,2); 
                    A0hat1(2,2)=-A0hat1(2,2);  
                    %A0hat1(2,3)=-A0hat1(2,3);
                    A0hat1(2,4)=-A0hat1(2,4);
                   
                    A0hat1(3,1)=-A0hat1(3,1);  
                    %A0hat1(3,2)=-A0hat1(3,2);  
                    A0hat1(3,3)=-A0hat1(3,3); 
                    %A0hat1(3,4)=-A0hat1(3,4); 
                    %A0hat1(3,5)=-A0hat1(3,5);  

                    A0hat1(4,1)=-A0hat1(4,1);  
                    %A0hat1(4,2)=-A0hat1(4,2);  
                    A0hat1(4,4)=-A0hat1(4,4);  
            
                    A0hat1(5,1)=-A0hat1(5,1);  
                    A0hat1(5,5)=-A0hat1(5,5);

                    A0hat1(6,1)=-A0hat1(6,1);  
                    A0hat1(6,2)=-A0hat1(6,2);

                    chck=10;
                end
           end
        end
        
 % Далее считаем имульсные отклики       
yhat1=zeros(irf_length,N); 
vhat1=zeros(irf_length,N); 
vhat1(L+1,1)=1; %Шок дкп
for j=L+1:irf_length
 yhat1(j,:)=[reshape(wrev(yhat1(j-L:j-1,:))', 1, L*N) 0 0]*reshape(beta,N*L+m,N)+vhat1(j,:)*A0hat1;
end


 out1(jj,:,:)=yhat1; % Отклики на шок дкп, 
 %1-й индекс в out1 отвечает за итерацию сэмплирования Гибса, второй за
 %горизонт прогноза отклика, 3 за номер переменной

 out8(jj,:,:)=beta;
 out9(jj,:,:)=Sigma;

 ETA=(inv(A0hat1)*e')'; %Матрица структурных шоков T на N

 l=loglik(reshape(beta,N*L+m,N),Sigma,Y,X);
 likelihood=[likelihood;l];

  %hd1=zeros(irf_length-L, N); %Для каждой переменной считаем накопленный вклад шока дкп
 %for t=L+1:irf_length
    %hd1(t-L, :)=[((a+c)/(a+b+c))*wrev(out1(jj, L+1:t, 1))', (-1)*wrev(out1(jj, L+1:t, 3))', wrev(out1(jj, L+1:t, 5))', wrev(out1(jj, L+1:t, 6))'];
 %end
 %hd_record1(jj, :, :)=hd1;

 jj=jj+1; 
    end 
end 

betam=squeeze(mean(out8,1));
sigmam=squeeze(mean(out9,1));
lik=loglik(reshape(betam,N*L+m,N),sigmam,Y,X);
BIC=log(N)*(N*L+m)-2*lik;

%% Долговая нагрузка
t=d7(1,1); %total
o=d7(2,1); %old
n=d7(3,1); %new

a=0.45*t*(1+mean(d1)/100); %взвешивание тотал
b=0.55*o*(1+mean(d1)/100); %взвешивание old
c=0.55*n*(1+mean(d1)/100); %взвешиние new

e=1*o*(1+mean(d1)/100); %взвешивание old
f=1*n*(1+mean(d1)/100); %взвешиние new

u=t/n; %взвешивание для new
v=o/n; %взвешивание для new

%for i=1:irf_length
    %if i <= 12
        %out01(:,i,1)=out1(:,i,5);
        %out01(:,i,2)=out1(:,i,6);
    %else
       %out01(:,i,1)=out1(:,i,5)-out1(:,i-1,5); %old в первых разностях
       %out01(:,i,2)=out1(:,i,6)-out1(:,i-1,6); %new в первых разностях
    %end
%end

%out4=out1(:,:,5)*0.829587462+out1(:,:,6)*0.170286756; %прирост total
%out5=out1(:,:,5)*0.829587462+out1(:,:,6)*0.170286756; %прирост new

for i=1:irf_length
    if i <= 12
       out01(:,i,1)=out1(:,i,1);
    else
       out01(:,i,1)=out1(:,i,1)-out1(:,i-1,1); %ставка в первых разностях
    end
end

kasha=(a/(a+b+c))*(out1(:,:,5)+out01(:,:,1))+(b/(a+b+c))*out1(:,:,6)+(c/(a+b+c))*(u*out1(:,:,5)-v*out1(:,:,6)+out01(:,:,1));

out10=kasha-out1(:,:,3); %моментальный отклик долговой нагрузки на шок дкп
out010=cumsum(out10,2); %накопленный отклик долговой нагрузки на шок дкп

%% Графический вывод импульсных откликов для базового сценария
figure(1);
 subplot(4,3,1)
  temp01=out1(:,:,1);
  temp1=squeeze(prctile(temp01,[16 84],1))';
  temp11=squeeze(prctile(temp01,50,1))';
  plot(temp1(L+1:L+25,:), '--', 'Color',"#7E2F8E", 'LineWidth',1);
  hold on
  plot(temp11(L+1:L+25,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('r, mp shock');
 axis tight

 subplot(4,3,2)
  temp02=out1(:,:,2);
  temp2=squeeze(prctile(temp02,[16 84],1))';
  temp22=squeeze(prctile(temp02,50,1))';
  plot(temp2(L+1:L+20,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp22(L+1:L+20,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('cpi, mp shock');
 axis tight

 subplot(4,3,3)
  temp03=out1(:,:,3);
  temp3=squeeze(prctile(temp03,[16 84],1))';
  temp33=squeeze(prctile(temp03,50,1))';
  plot(temp3(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp33(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('y ind, mp shock');
 axis tight

 subplot(4,3,4)
  temp04=out1(:,:,4);
  temp4=squeeze(prctile(temp04,[16 84],1))';
  temp44=squeeze(prctile(temp04,50,1))';
  plot(temp4(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp44(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('y, mp shock');
 axis tight

 subplot(4,3,5)
  temp05=out1(:,:,6);
  %temp05=out01(:,:,1);
  temp5=squeeze(prctile(temp05,[16 84],1))';
  temp55=squeeze(prctile(temp05,50,1))';
  plot(temp5(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp55(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('old debt, mp shock');
 axis tight

 subplot(4,3,6)
  temp06=u*out1(:,:,5)-v*out1(:,:,6);
  %temp06=out01(:,:,2);
  temp6=squeeze(prctile(temp06,[16 84],1))';
  temp66=squeeze(prctile(temp06,50,1))';
  plot(temp6(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp66(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('new credit, mp shock');
 axis tight 

 subplot(4,3,7)
  temp012=out1(:,:,5);
  temp12=squeeze(prctile(temp012,[16 84],1))';
  temp1212=squeeze(prctile(temp012,50,1))';
  plot(temp12(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp1212(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('total credit, mp shock');
 axis tight

 subplot(4,3,8)
  temp07=(a/(a+b+c))*(out1(:,:,5)+out01(:,:,1))+(b/(a+b+c))*out1(:,:,6)+(c/(a+b+c))*(u*out1(:,:,5)-v*out1(:,:,6)+out01(:,:,1))-out1(:,:,3);
  temp7=squeeze(prctile(temp07,[16 84],1))';
  temp77=squeeze(prctile(temp07,50,1))';
  plot(temp7(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp77(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('cd, mp shock');
 axis tight  

 subplot(4,3,9)
  temp08=cumsum((a/(a+b+c))*(out1(:,:,5)+out01(:,:,1))+(b/(a+b+c))*out1(:,:,6)+(c/(a+b+c))*(u*out1(:,:,5)-v*out1(:,:,6)+out01(:,:,1))-out1(:,:,3),2);
  temp8=squeeze(prctile(temp08,[16 84],1))';
  temp88=squeeze(prctile(temp08,50,1))';
  plot(temp8(L+1:L+37,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp88(L+1:L+37,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('cd, cum, mp shock');
 axis tight

 subplot(4,3,10)
  temp08=cumsum((a/(a+b+c))*((prctile(out1(:,:,5),50,1))+(prctile(out01(:,:,1),50,1)))+(b/(a+b+c))* ...
      (prctile(out1(:,:,6),50,1))+(c/(a+b+c))*(u*(prctile(out1(:,:,5),50,1))-v*(prctile(out1(:,:,6),50,1))+(prctile(out01(:,:,1),50,1)))-(prctile(out1(:,:,3),50,1)),2)';
  plot(temp08(L+1:L+37,:), 'Color',"#7E2F8E", 'LineWidth',1);
  yline(0,':');
  title('cd1, cum, mp shock');
 axis tight

 %% Графический вывод для жесткого сценария
figure(2);
 subplot(4,3,1)
  temp01=out1(:,:,1);
  temp1=squeeze(prctile(temp01,[16 84],1))';
  temp11=squeeze(prctile(temp01,50,1))';
  plot(temp1(L+1:L+25,:), '--', 'Color',"#7E2F8E", 'LineWidth',1);
  hold on
  plot(temp11(L+1:L+25,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('r, mp shock');
 axis tight

 subplot(4,3,2)
  temp02=out1(:,:,2);
  temp2=squeeze(prctile(temp02,[16 84],1))';
  temp22=squeeze(prctile(temp02,50,1))';
  plot(temp2(L+1:L+20,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp22(L+1:L+20,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('cpi, mp shock');
 axis tight

 subplot(4,3,3)
  temp03=out1(:,:,3);
  temp3=squeeze(prctile(temp03,[16 84],1))';
  temp33=squeeze(prctile(temp03,50,1))';
  plot(temp3(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp33(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('y ind, mp shock');
 axis tight

 subplot(4,3,4)
  temp04=out1(:,:,4);
  temp4=squeeze(prctile(temp04,[16 84],1))';
  temp44=squeeze(prctile(temp04,50,1))';
  plot(temp4(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp44(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('y, mp shock');
 axis tight

 subplot(4,3,5)
  temp05=out1(:,:,5);
  temp5=squeeze(prctile(temp05,[16 84],1))';
  temp55=squeeze(prctile(temp05,50,1))';
  plot(temp5(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp55(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('d ind, mp shock');
 axis tight

 subplot(4,3,6)
  temp06=out01(:,:,1)+out1(:,:,5)-out1(:,:,3);
  temp6=squeeze(prctile(temp06,[16 84],1))';
  temp66=squeeze(prctile(temp06,50,1))';
  plot(temp6(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp66(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('cd, mp shock');
 axis tight  

 subplot(4,3,7)
  temp07=cumsum(out01(:,:,1)+out1(:,:,5)-out1(:,:,3),2);
  temp7=squeeze(prctile(temp07,[16 84],1))';
  temp77=squeeze(prctile(temp07,50,1))';
  plot(temp7(L+1:L+37,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp77(L+1:L+37,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('cd, cum, mp shock');
 axis tight

 %% Графический вывод для мягкого сценария
figure(3);
 subplot(4,3,1)
  temp01=out1(:,:,1);
  temp1=squeeze(prctile(temp01,[16 84],1))';
  temp11=squeeze(prctile(temp01,50,1))';
  plot(temp1(L+1:L+25,:), '--', 'Color',"#7E2F8E", 'LineWidth',1);
  hold on
  plot(temp11(L+1:L+25,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('r, mp shock');
 axis tight

 subplot(4,3,2)
  temp02=out1(:,:,2);
  temp2=squeeze(prctile(temp02,[16 84],1))';
  temp22=squeeze(prctile(temp02,50,1))';
  plot(temp2(L+1:L+20,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp22(L+1:L+20,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('cpi, mp shock');
 axis tight

 subplot(4,3,3)
  temp03=out1(:,:,3);
  temp3=squeeze(prctile(temp03,[16 84],1))';
  temp33=squeeze(prctile(temp03,50,1))';
  plot(temp3(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp33(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('y ind, mp shock');
 axis tight

 subplot(4,3,4)
  temp04=out1(:,:,4);
  temp4=squeeze(prctile(temp04,[16 84],1))';
  temp44=squeeze(prctile(temp04,50,1))';
  plot(temp4(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp44(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('y, mp shock');
 axis tight

 subplot(4,3,5)
  temp05=out1(:,:,6);
  %temp05=out01(:,:,1);
  temp5=squeeze(prctile(temp05,[16 84],1))';
  temp55=squeeze(prctile(temp05,50,1))';
  plot(temp5(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp55(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('old debt, mp shock');
 axis tight

 subplot(4,3,6)
  temp06=u*out1(:,:,5)-v*out1(:,:,6);
  %temp06=out01(:,:,2);
  temp6=squeeze(prctile(temp06,[16 84],1))';
  temp66=squeeze(prctile(temp06,50,1))';
  plot(temp6(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp66(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('new credit, mp shock');
 axis tight 

 subplot(4,3,7)
  temp07=(e/(e+f))*out1(:,:,6)+(f/(e+f))*(u*out1(:,:,5)-v*out1(:,:,6)+out01(:,:,1))-out1(:,:,3);
  temp7=squeeze(prctile(temp07,[16 84],1))';
  temp77=squeeze(prctile(temp07,50,1))';
  plot(temp7(L+1:L+15,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp77(L+1:L+15,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('cd, mp shock');
 axis tight  

 subplot(4,3,8)
  temp08=cumsum((e/(e+f))*out1(:,:,6)+(f/(e+f))*(u*out1(:,:,5)-v*out1(:,:,6)+out01(:,:,1))-out1(:,:,3),2);
  temp8=squeeze(prctile(temp08,[16 84],1))';
  temp88=squeeze(prctile(temp08,50,1))';
  plot(temp8(L+1:L+37,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
  hold on
  plot(temp88(L+1:L+37,:), 'Color',"#7E2F8E", 'LineWidth',1);
  hold off
  yline(0,':');
  title('cd, cum, mp shock');
 axis tight

 %% Графический вывод импульсных откликов
figure(4);
 subplot(4,3,1)
  temp08=cumsum((a/(a+b+c))*(out1(:,:,5)+out01(:,:,1))+(b/(a+b+c))*out1(:,:,6)+(c/(a+b+c))*(u*out1(:,:,5)-v*out1(:,:,6)+out01(:,:,1))-out1(:,:,3),2);
  temp88=cumsum((e/(e+f))*out1(:,:,6)+(f/(e+f))*(u*out1(:,:,5)-v*out1(:,:,6)+out01(:,:,1))-out1(:,:,3),2);
  temp0=cumsum(out01(:,:,1)+out1(:,:,5)-out1(:,:,3),2);
  temp088=squeeze(prctile(temp08,50,1))';
  temp0888=squeeze(prctile(temp88,50,1))';
  temp00=squeeze(prctile(temp0,50,1))';
  plot(temp0888(L+1:L+37,:), '--', 'Color',"#000000", 'LineWidth',1 );
  hold on
  plot(temp088(L+1:L+37,:), 'Color',"#0072BD", 'LineWidth',1);
  hold on
  plot(temp00(L+1:L+37,:), '--', 'Color',"#A2142F", 'LineWidth',1 );
  hold off
  yline(0,':');
  title('cd, cum, mp shock');
 axis tight

  %% Таблица с итоговыми коэффициентами
  temp11=squeeze(prctile(out1(:,:,1),50,1))'; %накоп первой разности ставки
  temp08=cumsum((a/(a+b+c))*(out1(:,:,5)+out01(:,:,1))+(b/(a+b+c))*out1(:,:,6)+(c/(a+b+c))*(u*out1(:,:,5)-v*out1(:,:,6)+out01(:,:,1))-out1(:,:,3),2);
  temp88=cumsum((e/(e+f))*out1(:,:,6)+(f/(e+f))*(u*out1(:,:,5)-v*out1(:,:,6)+out01(:,:,1))-out1(:,:,3),2);
  temp0=cumsum(out01(:,:,1)+out1(:,:,5)-out1(:,:,3),2);
  temp088=squeeze(prctile(temp08,50,1))'; %базовый сценарий
  temp0888=squeeze(prctile(temp88,50,1))'; %мягкий сценарий
  temp00=squeeze(prctile(temp0,50,1))'; %жесткий сценарий

  soft=d7(4,1);
  base=d7(5,1);
  hard=d7(6,1);

  temp_5pp0=[(1+5/temp11(L+1, 1)*temp0888(L+1, 1)/100)*soft (1+5/temp11(L+1, 1)*temp088(L+1, 1)/100)*base (1+5/temp11(L+1, 1)*temp00(L+1, 1)/100)*hard];
  temp_10pp0=[(1+10/temp11(L+1, 1)*temp0888(L+1, 1)/100)*soft (1+10/temp11(L+1, 1)*temp088(L+1, 1)/100)*base (1+10/temp11(L+1, 1)*temp00(L+1, 1)/100)*hard];
  temp_15pp0=[(1+15/temp11(L+1, 1)*temp0888(L+1, 1)/100)*soft (1+15/temp11(L+1, 1)*temp088(L+1, 1)/100)*base (1+15/temp11(L+1, 1)*temp00(L+1, 1)/100)*hard];
  temp_all0=[temp_5pp0 temp_10pp0 temp_15pp0]';
  %мягкий/базовый/жесткий сценарий моментально

  temp_5pp6=[(1+5/temp11(L+1, 1)*temp0888(L+7, 1)/100)*soft (1+5/temp11(L+1, 1)*temp088(L+7, 1)/100)*base (1+5/temp11(L+1, 1)*temp00(L+7, 1)/100)*hard];
  temp_10pp6=[(1+10/temp11(L+1, 1)*temp0888(L+7, 1)/100)*soft (1+10/temp11(L+1, 1)*temp088(L+7, 1)/100)*base (1+10/temp11(L+1, 1)*temp00(L+7, 1)/100)*hard];
  temp_15pp6=[(1+15/temp11(L+1, 1)*temp0888(L+7, 1)/100)*soft (1+15/temp11(L+1, 1)*temp088(L+7, 1)/100)*base (1+15/temp11(L+1, 1)*temp00(L+7, 1)/100)*hard];
  temp_all6=[temp_5pp6 temp_10pp6 temp_15pp6]';
  %мягкий/базовый/жесткий сценарий через полгода

  temp_5pp12=[(1+5/temp11(L+1, 1)*temp0888(L+13, 1)/100)*soft (1+5/temp11(L+1, 1)*temp088(L+13, 1)/100)*base (1+5/temp11(L+1, 1)*temp00(L+13, 1)/100)*hard];
  temp_10pp12=[(1+10/temp11(L+1, 1)*temp0888(L+13, 1)/100)*soft (1+10/temp11(L+1, 1)*temp088(L+13, 1)/100)*base (1+10/temp11(L+1, 1)*temp00(L+13, 1)/100)*hard];
  temp_15pp12=[(1+15/temp11(L+1, 1)*temp0888(L+13, 1)/100)*soft (1+15/temp11(L+1, 1)*temp088(L+13, 1)/100)*base (1+15/temp11(L+1, 1)*temp00(L+13, 1)/100)*hard];
  temp_all12=[temp_5pp12 temp_10pp12 temp_15pp12]';
  %мягкий/базовый/жесткий сценарий через год

  temp_5pp24=[(1+5/temp11(L+1, 1)*temp0888(L+25, 1)/100)*soft (1+5/temp11(L+1, 1)*temp088(L+25, 1)/100)*base (1+5/temp11(L+1, 1)*temp00(L+25, 1)/100)*hard];
  temp_10pp24=[(1+10/temp11(L+1, 1)*temp0888(L+25, 1)/100)*soft (1+10/temp11(L+1, 1)*temp088(L+25, 1)/100)*base (1+10/temp11(L+1, 1)*temp00(L+25, 1)/100)*hard];
  temp_15pp24=[(1+15/temp11(L+1, 1)*temp0888(L+25, 1)/100)*soft (1+15/temp11(L+1, 1)*temp088(L+25, 1)/100)*base (1+15/temp11(L+1, 1)*temp00(L+25, 1)/100)*hard];
  temp_all24=[temp_5pp24 temp_10pp24 temp_15pp24]';
  %мягкий/базовый/жесткий сценарий через два года

  temp_5pp36=[(1+5/temp11(L+1, 1)*temp0888(L+37, 1)/100)*soft (1+5/temp11(L+1, 1)*temp088(L+37, 1)/100)*base (1+5/temp11(L+1, 1)*temp00(L+37, 1)/100)*hard];
  temp_10pp36=[(1+10/temp11(L+1, 1)*temp0888(L+37, 1)/100)*soft (1+10/temp11(L+1, 1)*temp088(L+37, 1)/100)*base (1+10/temp11(L+1, 1)*temp00(L+37, 1)/100)*hard];
  temp_15pp36=[(1+15/temp11(L+1, 1)*temp0888(L+37, 1)/100)*soft (1+15/temp11(L+1, 1)*temp088(L+37, 1)/100)*base (1+15/temp11(L+1, 1)*temp00(L+37, 1)/100)*hard];
  temp_all36=[temp_5pp36 temp_10pp36 temp_15pp36]';
  %мягкий/базовый/жесткий сценарий на третий год
  
  temp=[temp_all0 temp_all6 temp_all12 temp_all24 temp_all36]; %уровни долговой нагрузки при изменении ставки

  del_temp_5pp0=[(5/temp11(L+1, 1)*temp0888(L+1, 1)/100)*soft (5/temp11(L+1, 1)*temp088(L+1, 1)/100)*base (5/temp11(L+1, 1)*temp00(L+1, 1)/100)*hard];
  del_temp_10pp0=[(10/temp11(L+1, 1)*temp0888(L+1, 1)/100)*soft (10/temp11(L+1, 1)*temp088(L+1, 1)/100)*base (10/temp11(L+1, 1)*temp00(L+1, 1)/100)*hard];
  del_temp_15pp0=[(15/temp11(L+1, 1)*temp0888(L+1, 1)/100)*soft (15/temp11(L+1, 1)*temp088(L+1, 1)/100)*base (15/temp11(L+1, 1)*temp00(L+1, 1)/100)*hard];
  del_temp_all0=[del_temp_5pp0 del_temp_10pp0 del_temp_15pp0]';

  %мягкий/базовый/жесткий сценарий моментально

  del_temp_5pp6=[(5/temp11(L+1, 1)*temp0888(L+7, 1)/100)*soft (5/temp11(L+1, 1)*temp088(L+7, 1)/100)*base (5/temp11(L+1, 1)*temp00(L+7, 1)/100)*hard];
  del_temp_10pp6=[(10/temp11(L+1, 1)*temp0888(L+7, 1)/100)*soft (10/temp11(L+1, 1)*temp088(L+7, 1)/100)*base (10/temp11(L+1, 1)*temp00(L+7, 1)/100)*hard];
  del_temp_15pp6=[(15/temp11(L+1, 1)*temp0888(L+7, 1)/100)*soft (15/temp11(L+1, 1)*temp088(L+7, 1)/100)*base (15/temp11(L+1, 1)*temp00(L+7, 1)/100)*hard];
  del_temp_all6=[del_temp_5pp6 del_temp_10pp6 del_temp_15pp6]';
  %мягкий/базовый/жесткий сце нарий через полгода

  del_temp_5pp12=[(5/temp11(L+1, 1)*temp0888(L+13, 1)/100)*soft (5/temp11(L+1, 1)*temp088(L+13, 1)/100)*base (5/temp11(L+1, 1)*temp00(L+13, 1)/100)*hard];
  del_temp_10pp12=[(10/temp11(L+1, 1)*temp0888(L+13, 1)/100)*soft (10/temp11(L+1, 1)*temp088(L+13, 1)/100)*base (10/temp11(L+1, 1)*temp00(L+13, 1)/100)*hard];
  del_temp_15pp12=[(15/temp11(L+1, 1)*temp0888(L+13, 1)/100)*soft (15/temp11(L+1, 1)*temp088(L+13, 1)/100)*base (15/temp11(L+1, 1)*temp00(L+13, 1)/100)*hard];
  del_temp_all12=[del_temp_5pp12 del_temp_10pp12 del_temp_15pp12]';
  %мягкий/базовый/жесткий сценарий через год

  del_temp_5pp24=[(5/temp11(L+1, 1)*temp0888(L+25, 1)/100)*soft (5/temp11(L+1, 1)*temp088(L+25, 1)/100)*base (5/temp11(L+1, 1)*temp00(L+25, 1)/100)*hard];
  del_temp_10pp24=[(10/temp11(L+1, 1)*temp0888(L+25, 1)/100)*soft (10/temp11(L+1, 1)*temp088(L+25, 1)/100)*base (10/temp11(L+1, 1)*temp00(L+25, 1)/100)*hard];
  del_temp_15pp24=[(15/temp11(L+1, 1)*temp0888(L+25, 1)/100)*soft (15/temp11(L+1, 1)*temp088(L+25, 1)/100)*base (15/temp11(L+1, 1)*temp00(L+25, 1)/100)*hard];
  del_temp_all24=[del_temp_5pp24 del_temp_10pp24 del_temp_15pp24]';
  %мягкий/базовый/жесткий сценарий через два года

  del_temp_5pp36=[(5/temp11(L+1, 1)*temp0888(L+37, 1)/100)*soft (5/temp11(L+1, 1)*temp088(L+37, 1)/100)*base (5/temp11(L+1, 1)*temp00(L+37, 1)/100)*hard];
  del_temp_10pp36=[(10/temp11(L+1, 1)*temp0888(L+37, 1)/100)*soft (10/temp11(L+1, 1)*temp088(L+37, 1)/100)*base (10/temp11(L+1, 1)*temp00(L+37, 1)/100)*hard];
  del_temp_15pp36=[(15/temp11(L+1, 1)*temp0888(L+37, 1)/100)*soft (15/temp11(L+1, 1)*temp088(L+37, 1)/100)*base (15/temp11(L+1, 1)*temp00(L+37, 1)/100)*hard];
  del_temp_all36=[del_temp_5pp36 del_temp_10pp36 del_temp_15pp36]';
  %мягкий/базовый/жесткий сценарий на третий год

  del_temp=[del_temp_all0 del_temp_all6 del_temp_all12 del_temp_all24 del_temp_all36]; %изменение уровней долговой нагурзки при изменении ставки

 %% Графический вывод для разложения долговой нагрузки
% три фактора объясняют долговую нагрузку -> нужна декомпозиция отклика на три составляющих: 
% эффект повышения ставки, эффект снижения долга, эффект сокращения выпуска

figure(5);
A = datetime(2010,L+8,01);
B = datetime(2014,8,01);
C=A:calmonths(1):B;
  temp1=cumsum((prctile(out01(:,:,1),50,1)).*((a+c)/(a+b+c)),2); %ставка
  temp2=cumsum((prctile(out1(:,:,6),50,1)).*((a/(a+b+c))*(o/t)+(b/(a+b+c))),2); %старый долг
  temp3=cumsum((u*(prctile(out1(:,:,5),50,1))-v*(prctile(out1(:,:,6),50,1))).*((a/(a+b+c))*(n/t)+(c/(a+b+c))),2); %новый долг
  temp4=cumsum((prctile(out1(:,:,3),50,1)).*(-1),2); %выпуск

  temp01=out1(:,:,1);
  temp11=squeeze(prctile(temp01,50,1));
  coef=5/temp11(1,13); %необходимо для соизмеримости откликов

  temp11=(temp1(:,L+1:L+37)*coef)';
  temp21=(temp2(:,L+1:L+37)*coef)';
  temp31=(temp3(:,L+1:L+37)*coef)';
  temp41=(temp4(:,L+1:L+37)*coef)';
  bar(C', [temp11, temp21, temp41, temp31], "stacked");
  hold on
  l=temp11+temp21+temp31+temp41;
  plot(C, l,'Color',"#000000", 'LineWidth', 2.5);
  hold on
  title('Debt burden');
  ylim([-20 20]);
  legend('r','old debt', "y", "new debt");
  hold off;

  temp_all=[temp11 temp21 temp31 temp41]; %значения по hd по отрасли (r, old, new, y)
