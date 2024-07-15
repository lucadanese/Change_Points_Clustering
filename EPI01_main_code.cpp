#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// -----------------------------------------------------------------------------
// UTILS
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
double my_min(double a, 
              double b){
  if(a <= b){
    return a;
  } else {
    return b;
  }
}

// [[Rcpp::export]]
int rint(vec freqs){
  vec probs = cumsum(freqs / sum(freqs));
  double u = randu();
  for(uword i = 0; i < freqs.n_elem; i++){
    if(u <= probs(i)) {
      return i;
    }
  }
  return -1;
}

// [[Rcpp::export]]
int rbinom(int n, 
           double p){
  int out = 0;
  for(uword i = 0; i < n; i++){
    if(randu() < p){
      out += 1;
    }
  }
  return out;
}

// [[Rcpp::export]]
vec rmultin(int n, 
            vec probs){
  vec out(probs);
  out.fill(0);
  int temp;
  for(uword i = 0; i < n; i++){
    temp = rint(probs);
    out(temp) += 1;
  }
  return out;
}

// [[Rcpp::export]]
double log_sum_exp(vec log_vals){
  double M = max(log_vals);
  double lse = M + log(sum(exp(log_vals - M)));
  return lse;
}

// [[Rcpp::export]]
int vec_equal(mat v1, 
			  mat v2){
	int a = 0;
	if(std::equal(v1.begin(), v1.end(), v2.begin())){
		a = 1;
	}
	return a;
}


// -----------------------------------------------------------------------------
// LIKELIHOOD
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
vec DSA_curve(double dt, 
              int T,  
              double rhoval, 
              double gamma, 
              vec Betat, 
              double S0 = 1, 
              double R0 = 0){

  vec out(T);
  double t = 0;
  int j = 0;

  double St = S0 + dt * (-Betat(j) * rhoval);
  double It = rhoval + dt * (Betat(j) * rhoval - gamma * rhoval);
  double Rt = R0 + dt * (gamma * rhoval);
  
  for(uword i = 0; i < T/dt; i++){

    t = t + dt;
    St = St + dt * (-Betat(j) * St * It);
    It = It + dt * (Betat(j) * St * It - gamma*It);
    Rt = Rt + dt * (gamma * It);

    if(floor(t) == j + 1){
      out(j) = log(Betat(j)) + log(St) + log(It) ;
      j = j + 1;
    }
  }
  return(out - log(1 - St));
}

// [[Rcpp::export]]
Rcpp::List SIR_curve(double dt, 
              int T,  
			  vec order,
              double rhoval, 
              double gamma, 
			  double a0,
			  double b0,
              double S0 = 1, 
              double R0 = 0){

  vec out(T);
  
  
  vec t_beta(T); 
  
  t_beta(0) = randg(1, distr_param(a0, 1.0 / b0) )[0];
    for(uword i = 1; i < T; i++){
      if(order(i) == order(i - 1)){
        t_beta(i) = t_beta(i - 1);
      } else {
        t_beta(i) = randg(1, distr_param(a0, 1.0 / b0) )[0];
      }
    }
	
  double t = 0;
  int j = 0;
  
  vec Betat = t_beta;

  double St = S0 + dt * (-Betat(j) * rhoval);
  double It = rhoval + dt * (Betat(j) * rhoval - gamma * rhoval);
  double Rt = R0 + dt * (gamma * rhoval);
  
  for(uword i = 0; i < T/dt; i++){

    t = t + dt;
    St = St + dt * (-Betat(j) * St * It);
    It = It + dt * (Betat(j) * St * It - gamma*It);
    Rt = Rt + dt * (gamma * It);

    if(floor(t) == j + 1){
      out(j) = St;
      j = j + 1;
    }
  }
  
  
    return Rcpp::List::create(
    Rcpp::Named("st") = out,
    Rcpp::Named("beta") = Betat
  );
  
  //return(out);
}


// [[Rcpp::export]]
mat integrated_curves_mat(double dt, 
                          vec order, 
                          double a0, 
                          double b0, 
                          double gamma, 
                          double rhoval,
                          int M,
                          double S0 = 1, 
                          double R0 = 0){
  int T = order.n_elem, L = max(order) + 1;
  mat out(M, T + 1);
  vec t_beta(T);
  vec t_Imp(T);    
  double rho_temp, gamma_temp;
  
  for(uword m = 0; m < M; m++){
    
    t_beta(0) = randg(1, distr_param(a0, 1.0 / b0))[0];
	t_Imp(0) = R::dgamma(t_beta(0), a0, 1/b0, true);
	
	for(uword i = 1; i < T; i++){
      if(order(i) == order(i - 1)){
        t_beta(i) = t_beta(i - 1);
      } else {
		t_beta(i) = randg(1, distr_param(a0, 1.0 / b0) )[0];
		t_Imp.resize(t_Imp.n_elem + 1);
		t_Imp(t_Imp.n_elem - 1) = R::dgamma(t_beta(i), a0, 1/ b0, true);
      }
    }
	
    out.row(m).cols(0,T-1) = DSA_curve(dt, T, rhoval, gamma, t_beta, S0, R0).t();
	out.row(m).col(T) = sum(t_Imp);
	
	t_Imp.resize(1);
	
  }
  return out;
}

// [[Rcpp::export]]
vec generate_random_order(int T,
                          double p){
					   
					   
  vec freqs, temp_probs, new_order, cfreq;
  int k;
  double ord_lprob;
	
  k = 1 + rbinom(T, p);
  temp_probs.resize(k);
    
  for(uword l = 0; l < k; l++){
    temp_probs(l) = randg(1, distr_param(1.0, 0.2) )[0];
  }
  freqs = rmultin(T, temp_probs);
  new_order.resize(T);
  cfreq = cumsum(freqs);
    
  for(uword i = 0; i < cfreq(0); i++){
      new_order(i) = 0;
  }
  
  for(uword j = 1; j < cfreq.n_elem; j++){
      for(uword i = cfreq(j-1); i < cfreq(j); i++){
        new_order(i) = j;
      }
    }
    
  while(new_order(0) > 0){
     new_order -= 1;
  }
	
  return new_order;

}

// [[Rcpp::export]]
vec norm_constants(mat data,
                   double dt,
                   double a0,
                   double b0,
                   double gamma, 
                   vec rho,
                   int M,
                   int R,
                   double S0 = 1,
                   double R0 = 0, 
                   double p = 0.03){

  vec temp_llik_vec(data.n_rows), freqs, temp_probs, new_order, cfreq;
  mat temp_llik_mat(R,data.n_rows), curve_mat;
  int k, T = data.n_cols;
  double ord_lprob;

  //loop
  int res_index = 0;
  int start_s = clock();
  int current_s;
  int nupd = round(R / 10);
  for(uword r = 0; r < R; r++){

    k = 1 + rbinom(T, p);
    temp_probs.resize(k);
    for(uword l = 0; l < k; l++){
      temp_probs(l) = randg(1, distr_param(1.0, 0.2) )[0];
    }
    freqs = rmultin(T, temp_probs);
    ord_lprob = lgamma(T + 1) - lgamma(k + 1) - lgamma(T - k + 1) + 
      k * log(p) + (T - k) * log(1 - p) + lgamma(T + 1) - sum(lgamma(freqs + 1)) + 
      sum(freqs.t() * temp_probs / sum(temp_probs));
    new_order.resize(T);
    cfreq = cumsum(freqs);
    for(uword i = 0; i < cfreq(0); i++){
      new_order(i) = 0;
    }
    for(uword j = 1; j < cfreq.n_elem; j++){
      for(uword i = cfreq(j-1); i < cfreq(j); i++){
        new_order(i) = j;
      }
    }
    while(new_order(0) > 0){
      new_order -= 1;
    }
	
    //curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho(1), M, S0, R0);
    for(uword i = 0; i < data.n_rows; i++){
	  curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho(i), M, S0, R0);
      temp_llik_mat(r,i) = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M)  - ord_lprob;
    }
    
    // print time
    if((r + 1) % nupd == 0){
      current_s = clock();
      Rcpp::Rcout << "Normalization constant - completed:\t" << (r + 1) << "/" << R << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  }
  
  for(uword i = 0; i < data.n_rows; i++){
    temp_llik_vec(i) = log_sum_exp(temp_llik_mat.col(i));
  }

  return temp_llik_vec;
}

// -----------------------------------------------------------------------------
// UPDATE RHO
// -----------------------------------------------------------------------------


void update_rho(mat data, 
                vec &rho, 
                double a0, 
                double b0, 
                double c0, 
                double d0, 
                double MH_var,
                double gamma,
                double dt,
                int M,
                double S0, 
                double R0,
                vec &llik, 
                vec clust,
                mat orders){
					
  for(uword i = 0; i < data.n_rows; i++){
	  
	  double rho_temp = rho(i); 
	  double llik_temp = llik(i);
	  double new_llik;
	  double log_rho = log(rho_temp), log_rho_new, rho_new, acc_rate;
	  int T = data.n_cols;
	  mat curve_mat;
	  
	  log_rho_new = log_rho + randn() * sqrt(MH_var);
	  rho_new = exp(log_rho_new);
	  
	  
	  int clust_obs = clust(i);
	  curve_mat = integrated_curves_mat(dt, orders.row(clust_obs).t(), a0, b0, gamma, rho_new, M, S0, R0);
	  
	  new_llik = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);
	  
	  acc_rate = my_min(0, new_llik - llik_temp + (b0 - 1) * log(rho_new) - d0 * rho_new - (b0 - 1) * log(rho_temp) + d0 * rho_temp + rho_new - rho_temp);
	
	  if(log(randu()) < acc_rate){
			rho(i) = rho_new;
			llik(i) = new_llik;
		}	 
	} 
}




// -----------------------------------------------------------------------------
// UPDATE ORDERS
// -----------------------------------------------------------------------------


void update_single_order(mat data,
                         vec clust,
                         int clust_id,
                         mat &orders,
                         vec &llik,
                         double q,
                         double dt,
                         double a0,
                         double b0,
                         double gamma, 
                         vec rho,
                         int M,
                         double S0 = 1,
                         double R0 = 0,
						 double coarsening = 1){

  int k = max(orders.row(clust_id)) + 1, temp_id, temp_obs,
    T = orders.n_cols, bound = 0, temp_count;
  vec temp_llik = llik, temp_prob, new_order, freq_temp;
  mat curve_mat;
  bool check;
  double u = randu(), acc_rate;

  if(k == 1 || (u < q && k < T)){
    check = true;
  } else {
    check = false;
  }

  if(check == true){

    freq_temp.resize(k);
    for(uword l = 0; l < freq_temp.n_elem; l++){
	  if(std::count(orders.row(clust_id).begin(), orders.row(clust_id).end(), l) > 1){
        freq_temp(l) = 1;
      } else {
        freq_temp(l) = 0;
      }
    }

    temp_id = rint(freq_temp);
    temp_prob.resize(T);
    temp_prob.fill(0.0);
    for(uword i = 0; i < T; i++){
      if(orders(clust_id,i) == temp_id){
        temp_prob(i) += 1.0;
        bound = i;
      }
    }
    temp_prob(bound) = 0.0;
    temp_obs = rint(temp_prob);
    new_order = orders.row(clust_id).t();
    for(uword i = temp_obs + 1; i < T; i++){
      new_order(i) += 1;
    }

    //curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho, M, S0 = 1, R0 = 0);
    for(uword i = 0; i < clust.n_elem; i++){
      if(clust(i) == clust_id){
		curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho(i), M, S0 = 1, R0 = 0);
        temp_llik(i) = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);
      }
    }

    if(k == 1){
      acc_rate = my_min(0, log(1 - q) - log(q) + (coarsening * sum(temp_llik)) - (coarsening * sum(llik)) +
		log(sum(freq_temp)) + log(std::count(orders.row(clust_id).begin(), orders.row(clust_id).end(), temp_id)) - log(k));
    } else {
      acc_rate = my_min(0, log(1 - q) + log(T - 1) + (coarsening * sum(temp_llik)) - (coarsening * sum(llik)));
    }

    u = randu();
    if(u <= exp(acc_rate)){
      orders.row(clust_id) = new_order.t();
      llik = temp_llik;
    }

  } else {

    freq_temp.resize(k - 1);
    freq_temp.fill(1.0);
    temp_id = rint(freq_temp);
    new_order = orders.row(clust_id).t();
    for(uword i = 0; i < T; i++){
      if(orders(clust_id, i) > temp_id){
        new_order(i) -= 1;
      }
    }

    //curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho, M, S0 = 1, R0 = 0);
    for(uword i = 0; i < clust.n_elem; i++){
      if(clust(i) == clust_id){
		curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho(i), M, S0 = 1, R0 = 0);
        temp_llik(i) = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);
		
      }
    }

    temp_count = 0;
    for(uword i = 0; i < k - 1; i++){
	  if(std::count(orders.row(clust_id).begin(),orders.row(clust_id).end(), i) > 1){	  
        temp_count += 1;
      }
    }

    if(k == 1){
      acc_rate = my_min(0, log(q) - log(1 - q) + (coarsening * sum(temp_llik)) - (coarsening * sum(llik)) +
		log(k - 1) - log(temp_count) - log(std::count(new_order.begin(), new_order.end(), temp_id) - 1));
    } else {
      acc_rate = my_min(0, log(q) + log(T - 1) + (coarsening * sum(temp_llik)) - (coarsening * sum(llik)));
    }

    u = randu();
    if(u <= exp(acc_rate)){
      orders.row(clust_id) = new_order.t();
      llik = temp_llik;
    }
  }

  k = max(orders.row(clust_id)) + 1;
  if(k > 1){
    freq_temp.resize(k - 1);
    freq_temp.fill(1.0);
    temp_id = rint(freq_temp);

    temp_prob.resize(T);
    temp_prob.fill(0.0);
    for(uword i = 0; i < T; i++){
      if(orders(clust_id,i) == temp_id || orders(clust_id,i) == temp_id + 1){
        temp_prob(i) += 1.0;
        bound = i;
      }
    }
    temp_prob(bound) = 0.0;
    temp_obs = rint(temp_prob);
    new_order = orders.row(clust_id).t();
    for(uword i = temp_obs + 1; i < T; i++){
      if(new_order(i) == temp_id || new_order(i) == temp_id + 1){
        new_order(i) = temp_id + 1;
      }
    }

    //curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho, M, S0, R0);
    for(uword i = 0; i < clust.n_elem; i++){
      if(clust(i) == clust_id){
		curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho(i), M, S0, R0);  
        temp_llik(i) = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);
      }
    }

    acc_rate = my_min(0, (coarsening * sum(temp_llik)) - (coarsening * sum(llik)));

    u = randu();
    if(u <= exp(acc_rate)){
      orders.row(clust_id) = new_order.t();
      llik = temp_llik;
    }
  }
}

// -----------------------------------------------------------------------------
// UPDATE PARTITION
// -----------------------------------------------------------------------------

void update_partition(mat data,
                      vec &clust,
                      mat &orders,
                      vec &llik,
                      vec norm_const,
                      double alpha,
                      double p,
                      double q,
                      double dt,
                      double a0,
                      double b0,
                      double gamma, 
                      vec rho,
                      int M,
                      int L,
                      double S0 = 1,
                      double R0 = 0,
					  double coarsening = 1){
  vec temp_llik = llik, temp_clust = clust, freq_temp, prob_temp(clust.n_elem), 
    temp_vec1(clust.n_elem), temp_vec2(clust.n_elem), temp_vec3(clust.n_elem);
  int id1, id2, id3, id4, k, n, u_bound;
  double acc_rate;
  mat temp_order(2, data.n_cols), curve_mat1, curve_mat2, curve_mat3;
  temp_order.fill(0);

  freq_temp.resize(clust.n_elem);
  freq_temp.fill(1.0);
  id1 = rint(freq_temp);
  freq_temp(id1) = 0.0;
  id2 = rint(freq_temp);
  k = max(clust) + 1;
  n = clust.n_elem;

  if(clust(id1) != clust(id2)){
    
    prob_temp.fill(0.0);
    for(uword i = 0; i < clust.n_elem; i++){
      if(clust(i) == clust(id1) || clust(i) == clust(id2)){
        prob_temp(i) = 1.0;
      }
    }
		
    id3 = rint(prob_temp);
	
    temp_order.row(0) = orders.row(clust(id3));
    temp_clust.fill(1);
	
    temp_clust(id3) = 0;
    
    for(uword l = 0; l < L; l++){
      update_single_order(data, temp_clust, 0, temp_order, temp_llik, q, dt, 
                          a0, b0, gamma, rho, M, S0, R0);
    }
        	
	for(uword i = 0; i < clust.n_elem; i++){
	    curve_mat1 = integrated_curves_mat(dt, temp_order.row(0).t(), a0, b0, gamma, rho(i), M, S0, R0);
		curve_mat2 = integrated_curves_mat(dt, orders.row(clust(id1)).t(), a0, b0, gamma, rho(i), M, S0, R0);
		curve_mat3 = integrated_curves_mat(dt, orders.row(clust(id2)).t(), a0, b0, gamma, rho(i), M, S0, R0);	
      if(clust(i) == clust(id1) || clust(i) == clust(id2)){
		int T = data.n_cols;
        temp_llik(i) = log_sum_exp(curve_mat1.cols(0,T-1) * data.row(i).t() - curve_mat1.col(T)) - log(M);
        temp_vec1(i) = temp_llik(i) - norm_const(i);
		
		if(clust(i) == clust(id1)){
			llik(i) = log_sum_exp(curve_mat2.cols(0,T-1) * data.row(i).t() - curve_mat2.col(T)) - log(M);
			
			temp_vec2(i) = log_sum_exp(curve_mat2.cols(0,T-1) * data.row(i).t() - curve_mat2.col(T)) - log(M) - norm_const(i);
		} else if (clust(i) == clust(id2)){
			llik(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M);
			temp_vec2(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M) - norm_const(i);
		}
		
      } else {
		int T = data.n_cols;
        temp_vec1(i) = log_sum_exp(curve_mat1.cols(0,T-1) * data.row(i).t() - curve_mat1.col(T)) - log(M) - norm_const(i);
        temp_vec2(i) = log_sum_exp(curve_mat1.cols(0,T-1) * data.row(i).t() - curve_mat1.col(T)) - log(M) - norm_const(i);
      }
    }
    
	
	vec vec1 = regspace(1,k);
	vec vec2 = regspace(0,k-1);
	vec vec3 = regspace(1,k+1);
	vec vec4 = regspace(0,k);
	
	vec vecNpar1(vec1.n_elem);
	vecNpar1.fill(std::pow(2,(data.n_cols-1)));
	
	vec vecNpar2(vec2.n_elem);
	vecNpar2.fill(std::pow(2,(data.n_cols-1)));
	
	vec vecNpar3(vec3.n_elem);
	vecNpar3.fill(std::pow(2,(data.n_cols-1)));
	
	vec vecNpar4(vec4.n_elem);
	vecNpar4.fill(std::pow(2,(data.n_cols-1)));
	
	
	vecNpar1 = vecNpar1 - vec1;
	vecNpar2 = vecNpar2 - vec2;
	vecNpar3 = vecNpar3 - vec3;
	vecNpar4 = vecNpar4 - vec4;

	vecNpar1 = log(vecNpar1);
	vecNpar2 = log(vecNpar2);
	vecNpar3 = log(vecNpar3);
	vecNpar4 = log(vecNpar4);
	
	vec diffvec1 = vecNpar1 - vecNpar2;
	vec diffvec2 = vecNpar3 - vecNpar4;
	vec diffvec3 = vecNpar3 - vecNpar4;
		
	int sum1 = sum(diffvec1);
	int sum2 = sum(diffvec2);
	int sum3 = sum(diffvec3);
	
	vec sum1vec(2);
	sum1vec(0) = log(1); 
	sum1vec(1) = -sum1; 
	
	vec sum2vec(2);
	sum2vec(0) = log(1); 
	sum2vec(1) = -sum2; 
	
	vec sum3vec(2);
	sum3vec(0) = log(1); 
	sum3vec(1) = -sum3; 
	
	acc_rate = my_min(0, ((std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1)) +  std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2)) - 2) * log(0.5) + 
	  lgamma(alpha) + 
	  log_sum_exp(sum1vec) -
	  log_sum_exp(sum2vec) +
      lgamma(alpha + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1)) +  std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2))) - 
      lgamma(alpha + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1))) - lgamma(alpha + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2))) + 
      (coarsening * sum(temp_llik)) - (coarsening * sum(llik)) + 
	  log_sum_exp(temp_vec1) - log_sum_exp(temp_vec2)));
	  
	  
    if(log(randu()) < acc_rate){
      clust.elem(find(clust == clust(id2))).fill(clust(id1));
      orders.row(clust(id1)) = temp_order.row(0);
      llik = temp_llik;
    }
          
  } else {
    
	
    prob_temp.fill(0.0);
    for(uword i = 0; i < clust.n_elem; i++){
      if(clust(i) == clust(id1)){
        prob_temp(i) = 1.0;
      }
    }
	
	id3 = rint(prob_temp);
	prob_temp(id3) = 0;
	id4 = rint(prob_temp);
	
	prob_temp(id3) = 1;
	
	prob_temp(id1) = 0;
	prob_temp(id2) = 0;
	
	temp_clust(id2) = k;
	
	for(uword i = 0; i < clust.n_elem; i++){
      if(prob_temp(i) == 1){
		  
		if(randu() < 0.5){
			temp_clust(i) = temp_clust(id1);
		} else {
			temp_clust(i) = temp_clust(id2);
		}
      }
    }
	
	
    temp_order.row(0)= orders.row(temp_clust(id1)); 
	temp_order.row(1) = orders.row(temp_clust(id1));
	
	vec temp_clust_update(clust.n_elem);
	temp_clust_update.fill(2);
	
	temp_clust_update(id3) = 0; 
	temp_clust_update(id4) = 1;
	
    
    for(uword l = 0; l < L; l++){
						  
	  update_single_order(data, temp_clust_update, 0, temp_order, temp_llik, q, dt, 
                          a0, b0, gamma, rho, M, S0, R0);
						  
						  
						  
	  update_single_order(data, temp_clust_update, 1, temp_order, temp_llik, q, dt, 
                          a0, b0, gamma, rho, M, S0, R0);
    }
    
    for(uword i = 0; i < temp_clust.n_elem; i++){
		
		if(temp_clust(i) == clust(id1) && i != id1 && i != id2){

			if(randu() < 0.5){
				temp_clust(i) = clust(id1);
			  } else {
				temp_clust(i) = k;
			  }
			  
		}		
		
      
    }
	
	for(uword i = 0; i < clust.n_elem; i++){
		curve_mat1 = integrated_curves_mat(dt, temp_order.row(0).t(), a0, b0, gamma, rho(i), M, S0, R0);
		curve_mat2 = integrated_curves_mat(dt, temp_order.row(1).t(), a0, b0, gamma, rho(i), M, S0, R0);
		curve_mat3 = integrated_curves_mat(dt, orders.row(clust(id1)).t(), a0, b0, gamma, rho(i), M, S0, R0);
      if(temp_clust(i) == clust(id1)){
		int T = data.n_cols;
        temp_llik(i) = log_sum_exp(curve_mat1.cols(0,T-1) * data.row(i).t() - curve_mat1.col(T)) - log(M);
					
		llik(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M);
        temp_vec1(i) = temp_llik(i) - norm_const(i);
        temp_vec3(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M) - norm_const(i);
		
      } else if(temp_clust(i) == k){
		  
		int T = data.n_cols;
        temp_llik(i) = log_sum_exp(curve_mat2.cols(0,T-1) * data.row(i).t() - curve_mat2.col(T)) - log(M);
		
		llik(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M);
		
        temp_vec1(i) = temp_llik(i) - norm_const(i);;
        temp_vec3(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M) - norm_const(i);
		
		
      } else {
        temp_vec1(i) = temp_llik(i) - norm_const(i);
        temp_vec3(i) = temp_llik(i) - norm_const(i);
      }
    }
    
	vec vec1 = regspace(1,k+1);
	vec vec2 = regspace(0,k);
	vec vec3 = regspace(1,k);	
	vec vec4 = regspace(0,k-1); 
	
	vec vecNpar1(vec1.n_elem);
	vecNpar1.fill(std::pow(2,(data.n_cols-1)));
	
	vec vecNpar2(vec2.n_elem);
	vecNpar2.fill(std::pow(2,(data.n_cols-1)));
	
	vec vecNpar3(vec3.n_elem);
	vecNpar3.fill(std::pow(2,(data.n_cols-1)));
	
	vec vecNpar4(vec4.n_elem);
	vecNpar4.fill(std::pow(2,(data.n_cols-1)));
	
	vecNpar1 = vecNpar1 - vec1;
	vecNpar2 = vecNpar2 - vec2;
	vecNpar3 = vecNpar3 - vec3;
	vecNpar4 = vecNpar4 - vec4;
	
	vecNpar1 = log(vecNpar1);
	vecNpar2 = log(vecNpar2);
	vecNpar3 = log(vecNpar3);
	vecNpar4 = log(vecNpar4);
	
	vec diffvec1 = vecNpar1 - vecNpar2;
	vec diffvec2 = vecNpar3 - vecNpar4;
	vec diffvec3 = vecNpar3 - vecNpar4;

	int sum1 = sum(diffvec1);
	int sum2 = sum(diffvec2);
	int sum3 = sum(diffvec3);
	
	vec sum1vec(2);
	sum1vec(0) = log(1); 
	sum1vec(1) = -sum1; 
	
	vec sum2vec(2);
	sum2vec(0) = log(1); 
	sum2vec(1) = -sum2; 
	
	vec sum3vec(2);
	sum3vec(0) = log(1); 
	sum3vec(1) = -sum3; 
	
	acc_rate = my_min(0, (- (std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1)) +  std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2)) - 2) * log(0.5) + 
		- lgamma(alpha + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1)) + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2))) +
		log_sum_exp(sum1vec) -
	    log_sum_exp(sum3vec) +
		(coarsening * sum(temp_llik)) - (coarsening * sum(llik)) + log_sum_exp(temp_vec1) - 
		log_sum_exp(temp_vec3)));

    if(log(randu()) < acc_rate){
      clust = temp_clust;
      orders.resize(k + 1, orders.n_cols);
	  orders.row(clust(id1)) = temp_order.row(0);
	  
      orders.row(k) = temp_order.row(1);
      llik = temp_llik;
    }
	
	
  }
  
  k = orders.n_rows;
  for(uword i = 0; i < k; i++){
    
	if((int) std::count(clust.begin(), clust.end(), i) == 0){	
      for(uword j = k; j > i; j--){
        if((int) std::count(clust.begin(), clust.end(), j) != 0){  
          clust(find(clust == j)).fill(i);
          orders.swap_rows(i,j);
          break;
        }
      }
    }
  }
 
  u_bound = 1;
  for(uword i = 1; i < k + 1; i++){
	if(std::count(clust.begin(), clust.end(),i) > 0){	
      u_bound += 1;
    }
  }
  orders.resize(u_bound, orders.n_cols);

}

// -----------------------------------------------------------------------------
// FIX POSSIBLE EQUAL ORDERS AFTER ACC STEP
// -----------------------------------------------------------------------------


void fix_equal_order(vec &clust,
                           mat &orders){
							   
	int k = orders.n_rows;
	
	for(uword i = 0; i < k; i++){

		for(uword j = k-1; j > i; j--){
			
			if(vec_equal(orders.row(j),orders.row(i)) == 1.0){
				clust(find(clust == j)).fill(i);
				
			}	

		}
		
	}
	
	k = orders.n_rows;
	for(uword i = 0; i < k; i++){
    
		if((int) std::count(clust.begin(), clust.end(), i) == 0){	
		  for(uword j = k; j > i; j--){
			if((int) std::count(clust.begin(), clust.end(), j) != 0){  
			  clust(find(clust == j)).fill(i);
			  orders.swap_rows(i,j);
			  break;
			}
		  }
		}
	  }
	 
	  int u_bound = 1;
	  for(uword i = 1; i < k + 1; i++){
		if(std::count(clust.begin(), clust.end(),i) > 0){	
		  u_bound += 1;
		}
	  }
	  orders.resize(u_bound, orders.n_cols);
	}


// -----------------------------------------------------------------------------
// ESTIMATE MARGINAL CPs
// -----------------------------------------------------------------------------


// [[Rcpp::export]]
Rcpp::List marginal_CP(mat data,
                       int niter,
                       int nburn,
                       double alpha,
                       double q,
                       double dt,
                       double a0,
                       double b0,
                       double c0, 
                       double d0,
                       double gamma, 
                       double MH_var,
                       int M,
                       int R,
                       int L,
                       double S0 = 1,
                       double R0 = 0,
                       double p = 0.003, 
                       int nupd = 0,
					   double coarsening = 1){


  
	int n = data.n_rows, T = data.n_cols;  
	vec rho(n);
    rho.fill(0.001);
    vec clust(n), llik(n);
    clust.fill(0);
    mat orders(n, T);
    //mat orders(1, T);
    orders.fill(0);
	
	for(uword i = 0; i < orders.n_rows; i++){
	  
	  orders.row(i) = generate_random_order(T, p).t();
	  
	}
   
	for(uword i = 0; i < clust.n_elem; i++){
	  
		mat curve_mat = integrated_curves_mat(dt, orders.row(0).t(), a0, b0, gamma, rho(i), M, S0, R0);
		llik(i) = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);
	
	}
  
  
	mat res_clust(niter - nburn, n);
	cube res_orders(n, T, niter - nburn);
	mat res_llik(niter - nburn, n);
	mat res_rho(niter - nburn, n);
  
	//loop
	int res_index = 0;
	int start_s = clock();
	int current_s;
	if(nupd == 0){
		int nupd = round(niter / 10);
	}

	Rcpp::Rcout << "\n------ MAIN LOOP ------\n\n";
	// start
	for(uword iter = 0; iter < niter; iter++){
	  
		update_rho(data, rho, a0, b0, c0, d0, MH_var, gamma, dt, M,
		           S0, R0, llik, clust, orders);
    
		for(uword j = 0; j < orders.n_rows; j++){
		update_single_order(data, clust, j, orders, llik,
		                    q, dt, a0, b0, gamma, rho, M, S0, R0, coarsening);
		}
	
		if(iter >= nburn){
		  res_clust.row(iter-nburn) = clust.t();
		  res_orders.slice(iter-nburn) = orders;
		  res_llik.row(iter-nburn).cols(0,n-1) = llik.t();
		  res_rho.row(iter-nburn).cols(0,n-1) = rho.t();
		}
		// print time
		if((iter + 1) % nupd == 0){
		  current_s = clock();
		  Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << " - in " <<
			double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
		}
		Rcpp::checkUserInterrupt();
	}
	double time = double(current_s-start_s)/CLOCKS_PER_SEC;

	Rcpp::List results;
	results["clust"] = res_clust;
	results["orders"] = res_orders;
	results["llik"] = res_llik;
	results["rho"] = res_rho;
	results["time"] = time;
  return results;
}


// -----------------------------------------------------------------------------
// MAIN FUNCTION
// -----------------------------------------------------------------------------


// [[Rcpp::export]]
Rcpp::List main_function(mat data,
                         int niter,
                         int nburn,
                         double alpha,
                         double q,
                         double dt,
                         double a0,
                         double b0,
                         double c0, 
                         double d0,
                         double gamma, 
                         double MH_var,
                         int M,
                         int B,
                         int L,
                         double S0 = 1,
                         double R0 = 0,
                         double p = 0.003, 
                         int nupd = 0,
						 double coarsening = 1){


  // preprocessing - constant
  int n = data.n_rows, T = data.n_cols;
  
  vec rho(n);
  rho.fill(0.0045);
  
  vec clust(n), llik(n);
  clust = regspace(0, n-1);
  mat orders(n, T);
  orders.fill(0);
  
  vec norm_vec = norm_constants(data, dt, a0, b0, gamma, rho, M, B, S0, R0, p);
  
    
  for(uword i = 0; i < orders.n_rows; i++){
	  
	  orders.row(i) = generate_random_order(T, p).t();
	  
  }
  
  
  //mat curve_mat = integrated_curves_mat(dt, orders.row(0).t(), a0, b0, gamma, rho, M, S0, R0);
  for(uword i = 0; i < clust.n_elem; i++){
	mat curve_mat = integrated_curves_mat(dt, orders.row(0).t(), a0, b0, gamma, rho(i), M, S0, R0);
    llik(i) = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);
  }
  for(uword l = 0; l < 1; l++){
    for(uword j = 0; j < orders.n_rows; j++){
     update_single_order(data, clust, j, orders, llik,
                        q, dt, a0, b0, gamma, rho, M, S0, R0, coarsening);
    }
  }
  mat res_clust(niter - nburn, n);
  cube res_orders(n, T, niter - nburn);
  mat res_llik(niter - nburn, n);
  //vec res_rho(niter - nburn);
  mat res_rho(niter - nburn, n);
  
  //loop
  int res_index = 0;
  int start_s = clock();
  int current_s;
  if(nupd == 0){
    int nupd = round(niter / 10);
  }

  Rcpp::Rcout << "\n------ MAIN LOOP ------\n\n";
  // start
  for(uword iter = 0; iter < niter; iter++){
	
	
    update_rho(data, rho, a0, b0, c0, d0, MH_var, gamma, dt, M,
               S0, R0, llik, clust, orders);
    
	
    update_partition(data, clust, orders, llik, norm_vec, alpha, p, q,
                     dt, a0, b0, gamma, rho, M, L, S0, R0,coarsening);
    
	
    for(uword j = 0; j < orders.n_rows; j++){
      update_single_order(data, clust, j, orders, llik,
                          q, dt, a0, b0, gamma, rho, M, S0, R0, coarsening);
    }

    if(iter >= nburn){
      res_clust.row(iter-nburn) = clust.t();
      res_llik.row(iter-nburn).cols(0,n-1) = llik.t();
      res_rho.row(iter-nburn).cols(0,n-1) = rho.t();
	  res_orders.slice(iter-nburn).rows(0, orders.n_rows-1) = orders;
    }
    // print time
    if((iter + 1) % nupd == 0){
      current_s = clock();
      Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  
  }
  double time = double(current_s-start_s)/CLOCKS_PER_SEC;
  
  Rcpp::List results;
  results["clust"] = res_clust;
  results["orders"] = res_orders;
  results["llik"] = res_llik;
  results["rho"] = res_rho;
  results["time"] = time;
  return results;
}

// -----------------------------------------------------------------------------
// PARTITION ESTIMATE
// -----------------------------------------------------------------------------

//[[Rcpp::export]]
mat psm(mat M){
  // initialize results
  mat result(M.n_cols, M.n_cols, fill::zeros);
  
  for(uword i = 0; i < M.n_cols; i++){
    for(uword j = 0; j <= i; j++){
      result(i,j) = sum(M.col(i) == M.col(j));
	  
      result(j,i) = result(i,j);
    }
    Rcpp::checkUserInterrupt();
  }
  return(result / M.n_rows);
}

//[[Rcpp::export]]
vec VI_LB(mat C_mat, mat psm_mat){
  
  vec result(C_mat.n_rows);
  double f = 0.0;
  int n = psm_mat.n_cols;
  vec tvec(n);
  
  for(uword j = 0; j < C_mat.n_rows; j++){
    f = 0.0;
    for(uword i = 0; i < n; i++){
      tvec = psm_mat.col(i);
      f += (log2(sum(C_mat.row(j) == C_mat(j,i))) +
        log2(sum(tvec)) -
        2 * log2(sum(tvec.elem(find(C_mat.row(j).t() == C_mat(j,i))))))/n;
    }
    result(j) = f;
    Rcpp::checkUserInterrupt();
  }
  return(result);
}

