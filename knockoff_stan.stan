// The dummy variable should include the baseline

data {
  int<lower = 1> N;         // num obs
  int<lower = 1> lengthZ;   // total number of covariates after we convert categorical covariates into dummy variables, after removing the baseline
  int<lower = 0> contNum;   // number of continuous covariates
  int<lower = 1> cateNum;   // number of categorical covariates
  int<lower = 1> cateDummyNum; // The number of categorical covariates after we convert them into binary dummy variables 

  int cumCateDummyNum[cateNum+1]; // cumulative number of the factors
  int cumCateDummyNumInter[cateNum+1]; //  cumulative number of the factors, after removing the baseline
  vector[lengthZ] mu0;                
  matrix[lengthZ,lengthZ] tau0;            
  matrix[N, contNum] contX; // The data matrix that contains only the continuous covariates
  int<lower=1> cateX[N, cateNum]; 
}

parameters {
  vector<lower=-6.0, upper=6.0>[lengthZ] Z; 
  cov_matrix[lengthZ] Sigma; 
  real<lower=0> tau;
}

transformed parameters {
	real<lower=0> sigma;
	real theta[cateDummyNum];
	
    for (i in 1:cateNum) {
        theta[1+cumCateDummyNum[i]] = 0;
        for (j in 2:(cumCateDummyNum[i+1]-cumCateDummyNum[i])) {
            theta[j+cumCateDummyNum[i]] = Z[j-1+cumCateDummyNumInter[i]+contNum];
        }
	}
	sigma=inv_sqrt(tau);
}

model {
    for(n in 1:N){
		if(contNum>0){
			for(s in 1:contNum){
				contX[n,s] ~ normal(Z[s],sigma);
			}
		}
		for (i in 1:cateNum) {
			cateX[n,i] ~ categorical_logit(to_vector(theta[(1+cumCateDummyNum[i]):cumCateDummyNum[i+1]]));
		}
	}
	Z ~  multi_normal_prec(mu0 , Sigma);
	Sigma ~ wishart(lengthZ,tau0);
    tau ~ gamma(0.001, 0.001);
}