r = 0.06 #interest rate
#volatility = 0.2 #volatility
s_0 = 40 #price of the underlying at time 0
#strike = 40 #strike
numberofpaths = 1000 #number of paths
N = 50 #number of exercise times
#T = 1 #Echéance

#Data<- data.frame("s","sigma","T","BSEuro","estimated_american","diff_american_BSeuro");
Vs=c();
Vsigma=c();
VT=c();
VBSEuro=c();
Vestimated_american=c();
Vdiff_american_BSeuro=c();
Price=c();
for(strike in c(38,40,42,44)){
  for(volatility in c(0.2,0.4)){
    for(T in c(1,2)){
    

#preallocation of the below vectors and matrices
cashflow = matrix(0L,numberofpaths,1);
W = matrix(0L,numberofpaths,N);  #matrice des sous-jacent

s_t = matrix(0L,numberofpaths,N)

M=matrix(0L,numberofpaths,N);#Matrice des cash flow actualisé 


set.seed(0) #to use seed 0
t=1:N; #N exercise events
time=t/N; #normalising the time frame to 1 year

s_t[,1]=s_0*matrix(1,numberofpaths,1)
dt=T/N


#Simulation des sous-jacents
for (i in 2:N){
  increments = as.vector(rnorm(numberofpaths,mean=0,sd=1)) 
  #the stock price following a geometric Brownian motion
  s_t[,i]=  s_t[,i-1]* exp((r - 0.5 * volatility^2) * dt + volatility*sqrt(dt)*increments) 
  W[,i] = s_t[,i]
}
W[,1]=s_t[,1]

#Etape 1   
#calcul des cashflow à l'instant N s'i est in the money
for( i in 1:numberofpaths){
  if( W[i,N]<strike ){
    cashflow[i] = (strike - W[i,N]) 
  }
  else{
    cashflow[i] = 0
  }
}
estimated_american = sum(cashflow[1:numberofpaths])/numberofpaths # european_value = cashflow[N]/numberofpaths

R0=function(x){rep(1,length(x))}
R1=function(x){(1-x)}
R2=function(x){(1/2*(x^2-4*x+2))}

for(K in 1:(N-2)){
  print(K);
  if(K>1){
    cashflow = matrix(0L,numberofpaths,1);
  for( i in 1:numberofpaths){
    if( W[i,N-K+1]<strike ){
      cashflow[i] = M[i,N-K+1] ;
   
    }
    else{
      cashflow[i] = 0
    }
  }
  }
X= c()  ;#Vecteur des sous-jacent
Y=c();#vecteur des cash-flow actualisé
index= list();
for(i in 1:numberofpaths){
  if(W[i,N-K]<strike){
    index=c(index,i)
    X=c(X,W[i,N-K]);
    Y=c(Y,cashflow[i]*exp(-r));
  }
}

Espcond= matrix(0L,length(index),1);

Phi=cbind(R0(X),R1(X),R2(X)) 
A=t(Phi)%*%Phi 
B=t(Phi)%*%Y
alpha=solve(A,B) 
Espcond=Phi%*%alpha 
Espcond
alpha
### etape 5 calcule valeur d'exercice K-X
exercice= matrix(0L,length(index),1);
for(i in 1:length(index)){
  exercice[i]=strike-X[i]
}
exercice
###Etape 6 
Cn_1=c();
for(i in 1:length(index)){
if(Espcond[i]>exercice[i]){
  Cn_1[i]=0;
}else{
  Cn_1[i]=exercice[i];
}
 ##Remplissage de la matrice M 
}
for(i in 1:numberofpaths){
  if(i %in% index){
    ind= which(index==i)
    M[i,N-K]= Cn_1[ind]
  }else{
    M[i,N-K]=0
  }
}}
for( i in 1:numberofpaths){
  if( W[i,N]<strike ){
    M[i,N] = (strike - W[i,N]) 
  }
  else{
    M[i,N] = 0
  }
}
M

M1=M

maximum=0
vect_maximum=c()
ligne = c(length = N)
for(i in 1:numberofpaths){
  ligne = M1[i,]
  maximum = max (ligne)
  vect_maximum=c(vect_maximum,maximum)
  for(j in 1:N) {
    if (M1[i,j] < maximum){
      M1[i,j] = 0
    }
  }
}
M1
valeur_americaine = sum(vect_maximum[1:numberofpaths])/numberofpaths


###Black & Scholes
value = 0
euro_vect = c()

for( j in (1:numberofpaths)){
  d1= (log(W[j,N]/strike)+(r + (volatility^2)/2)*T ) / (volatility*sqrt(T));
  d2 = d1 - (volatility*sqrt(T));
  value = -(W[j,N]*pnorm(-d1) - strike*exp(-r*T)*pnorm(-d2));
  euro_vect=c(euro_vect,value);
}


print(estimated_american);
print(valeur_americaine);

BSeuro = sum(euro_vect)/numberofpaths
print(BSeuro);
V=estimated_american-BSeuro;
print(strike)
print(V)
print(T)
Vs=c(Vs,strike);
Vsigma=c(Vsigma,volatility);
VT=c(VT,T);
VBSEuro=c(VBSEuro,BSeuro);
Vestimated_american=c(Vestimated_american,estimated_american);
Vdiff_american_BSeuro=c(Vdiff_american_BSeuro,V);
Price=c(Price,valeur_americaine)


# Data[nrow(Data)+1,]<- c(S=strike,sigma=volatility,T=T,BSEuro=BSeuro,estimated_american=estimated_american,estimated_american.BSeuro=V)
# Data
    }
  }
}

df <- data.frame(Vs,Vsigma,VT,VBSEuro,Vestimated_american,Vdiff_american_BSeuro,Price)
names(df) <- c('Strike', 'Sigma', 'T', 'Closed from European','Simulated American','difference','Price')

mean(Price)

