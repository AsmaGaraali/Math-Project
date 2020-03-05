r = 0.06                  #Taux d'intérêt
s_0 = 40                  #Prix du sous-jacent initial
numberofpaths = 100000    #Nombre des trajectoires simulées
N = 50                    #Nombre desinstants d'exercice possibles

#Création du dataset des résultats
df = data.frame("s","sigma","T","BSEuro","estimated_american","estimated_american-BSeuro")

strikes = c(38,40,42,44)  #Initialisation du vecteur des Prix d'exercice
volatilities = c(0.2,0.4) #Initialisation du vecteur des Volatilités
maturities = c(1,2)       #Initialisation du vecteur des Prix d'exercice


for( strike in 1:length(strikes)){
  for (volatility in 1:length(volatilities)){
    for(T in 1:length(maturities)){
      
      #Preallocation des vecteurs et des matrices et leurs initialisation
      cashflow = matrix(0L,numberofpaths,1);    #Vecteur des flux de trésorerie
      W = matrix(0L,numberofpaths,N);           #Matrice des sous-jacents
      s_t = matrix(0L,numberofpaths,N);         #Matrice des sous-jacents
      M=matrix(0L,numberofpaths,N);             #Matrice des cash flow actualisé 
      
      set.seed(0);          #to use seed 0
      t=1:N;                #N évènements d'exercice
      dt=maturities[T]/N;   #Discrétisation de l'intervalle de temps
      
      #Initialisation de la colonne du premier instant par la valeur de S0
      s_t[,1]=s_0*matrix(1,numberofpaths,1);   
      
      #Simulation des sous-jacents pour le reste des instants
      for (i in 2:N){
        increments = as.vector(rnorm(numberofpaths,mean=0,sd=1));
        s_t[,i]=  s_t[,i-1]* exp((r - 0.5 * (volatilities[volatility]^2)) * dt + (volatilities[volatility])*sqrt(dt)*increments) ;
        W[,i] = s_t[,i];
      }
      W[,1]=s_t[,1];
      
      
      #Calcul des flux de trésorerie à l'instant N si l'option est dans la monnaie
      for(i in 1:numberofpaths){
        if(W[i,N] < strikes[strike]){
          cashflow[i] = (strikes[strike] - W[i,N]) 
        }
        else{
          cashflow[i] = 0
        }
      }
      
      #Définition/Réinitialisation des polynômes de Laguerre (Pour l'approximation moindres carrés)
      R0=function(x){rep(1,length(x))}
      R1=function(x){(1-x)}
      R2=function(x){(1/2*(x^2-4*x+2))}
      
      #Rétrograde
      for(K in 1:(N-2)){
        #Calcul du flux de trésorerie à l'instant N-K+1
        if(K>1){
          cashflow = matrix(0L,numberofpaths,1);
          for( i in 1:numberofpaths){
            if( W[i,N-K+1]<strikes[strike] ){
              cashflow[i] = (strikes[strike]- W[i,N-K+1]) ;
              W[i,N-K+1];
            }
            else{
              cashflow[i] = 0
            }
          }
        }
        
        X=c();          #Vecteur des sous-jacent dans la monnaie pour l'instant N-K+1
        Y=c();          #vecteur des cash-flow actualisés des options dans la monnaie pour l'instant N-K+1
        index= list();  #Liste des indices des chemins dans la monnaie
        
        #Sauvegarde des sous-jacents et des cash-flow actualisés lorsque l'option est dans la monnaie
        for(i in 1:numberofpaths){
          if(W[i,N-K]<strikes[strike]){
            index=c(index,i);
            X=c(X,W[i,N-K]);
            Y=c(Y,cashflow[i]*exp(-r));
          }
        }
        
        #Allocation de l'espace memoire et Initialisation du vecteur de l'esperance conditionnelle
        Espcond= matrix(0L,length(index),1);
        
        #Résolution du système matriciel pour déterminer les alpha_i
        Phi=cbind(R0(X),R1(X),R2(X)) ;
        A=t(Phi)%*%Phi ;
        B=t(Phi)%*%Y;
        alpha=solve(A,B); 
        Espcond=Phi%*%alpha; 

        #Détermination du cash-flow par rapport au sous-jacent pour les chemins dans la monnaie
        exercice= matrix(0L,length(index),1);
        for(i in 1:length(index)){
          exercice[i]=strikes[strike]-X[i]
        }
        
        #Vecteur de la valeur de l'option pour chaque chemin dans la monnaie
        Cn_1=c();         
        
        #Détermination des cash-flows de chaque chemin dans la monnaie 
        #(Possibilité d'exercice ou non)
        for(i in 1:length(index)){
          
          # Comparaison de la valeur espérée avec celle calculée
          if(Espcond[i]>exercice[i]){
            Cn_1[i]=0;            #On n'exerce pas l'option à cet instant pour ce chemin
          }
          else{
            Cn_1[i]=exercice[i];  #L'option peut être exercée à cet instant pour ce chemin
          }
        }
        
        #Sauvegarde du vecteur obtenu dans une matrice M dans les bons indices des chemins
        #Et remplissage du reste des cases par des 0
        for(i in 1:numberofpaths){
          if(i %in% index){
            ind= which(index==i);
            M[i,N-K]= Cn_1[ind];
            cashflow[i] = Cn_1[ind];
          }
          else{
            M[i,N-K]=0;
            cashflow[i]=0;
          }
        }
      }
      #Fin de la boucle Rétrograde
      
      #Remplissage de la dernière colonne de la matrice M
      for( i in 1:numberofpaths){
        if( W[i,N]<strikes[strike] ){
          M[i,N] = (strikes[strike]- W[i,N]) ;
        }
        else{
          M[i,N] = 0
        }
      }
      M1=M;
      
      maximum=0;
      vect_maximum=c();
      ligne = c(length = N);
      
      #On garde seulement la valeur maximale d'exercice pour chaque chemin si elle existe
      for(i in 1:numberofpaths){
        ligne = M1[i,];
        maximum = max (ligne);
        vect_maximum=c(vect_maximum,maximum);
        for(j in 1:N) {
          if (M1[i,j] < maximum){
            M1[i,j] = 0
          }
        }
      }
      
      
      #Calul de la valeur finale estimée par LSM à l'echéance
      p =0;
      for(b in 1:numberofpaths) {
        if (M1[b,N] > 0){
          p=p+1
        }
      }
      estimated_american = sum(M1[,N])/p
      
      
      #Calcul de la valeur finale par Black & Scholes à l'echéance (Cas de l'option européenne)
      value = 0;
      euro_vect = c();
      
      for( j in (1:numberofpaths)){
        d1 =( log(W[j,N]/strikes[strike])+(r + (volatilities[volatility]^2)/2)*maturities[T] ) / (volatilities[volatility]*sqrt(maturities[T]));
        d2 = d1 - (volatilities[volatility]*sqrt(maturities[T]));
        value = -(W[j,N]*pnorm(-d1) - strikes[strike]*exp(-r*maturities[T])*pnorm(-d2));
        euro_vect=c(euro_vect,value);
      }
      BSeuro = sum(euro_vect)/numberofpaths;
      
      
      #Ajout du résultat à la dataset des résultats finaux
      vv= data.frame(toString(strikes[strike]),toString(volatilities[volatility]),toString(maturities[T]),toString(BSeuro),toString(estimated_american),toString(estimated_american-BSeuro));
      names(vv)= names(df)
      df = rbind(df,vv)
    }
  }
}
df

