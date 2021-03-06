functions {
  vector SEIR_matrix_sym(real time, vector y, real[] params) {
    vector[30] dydt;
    real InR1;
    real RR1;
    real InR2;
    real RR2;
    real InR3;
    real RR3;
    real InR4;
    real RR4;
    real k21;
    real k31;
    real k32;
    real k41;
    real k42;
    real k43;
    real RI1;
    real RI2;
    real RI3;
    real RI4;
    real lambda1;
    real InR5;
    real RR5;
    real RI5;
    real lambda2;
    real lambda3;
    real lambda4;
    real k54;
    real k53;
    real k51;
    real k52;
    real InR6;
    real RR6;
    real RI6;
    real k65;
    real k64;
    real k63;
    real k61;
    real k62;
    real IR1;
    real IR2;
    real IR3;
    real IR4;
    real lambda5;
    real lambda6;
    real IR5;
    real IR6;
    InR1 = y[2]/1;
    RR1 = y[3]/2;
    InR2 = y[6]/1;
    RR2 = y[7]/2;
    InR3 = y[10]/1;
    RR3 = y[11]/2;
    InR4 = y[14]/1;
    RR4 = y[15]/2;
    k21 = params[2];
    k31 = params[3];
    k32 = params[8];
    k41 = params[4];
    k42 = params[9];
    k43 = params[13];
    RI1 = InR1;
    RI2 = InR2;
    RI3 = InR3;
    RI4 = InR4;
    lambda1 = params[1]*y[3]+params[2]*y[7]+params[3]*y[11]+params[4]*y[15]+params[5]*y[23]+params[6]*y[28];
    InR5 = y[22]/1;
    RR5 = y[23]/2;
    RI5 = InR5;
    lambda2 = k21*y[3]+params[7]*y[7]+params[8]*y[11]+params[9]*y[15]+params[10]*y[23]+params[11]*y[28];
    lambda3 = k31*y[3]+k32*y[7]+params[12]*y[11]+params[13]*y[15]+params[14]*y[23]+params[15]*y[28];
    lambda4 = k41*y[3]+k42*y[7]+k43*y[11]+params[16]*y[15]+params[17]*y[23]+params[18]*y[28];
    k54 = params[17];
    k53 = params[14];
    k51 = params[5];
    k52 = params[10];
    InR6 = y[27]/1;
    RR6 = y[28]/2;
    RI6 = InR6;
    k65 = params[20];
    k64 = params[18];
    k63 = params[15];
    k61 = params[6];
    k62 = params[11];
    IR1 = y[1]*lambda1/1e+05;
    IR2 = y[5]*lambda2/1e+05;
    IR3 = y[9]*lambda3/1e+05;
    IR4 = y[13]*lambda4/1e+05;
    lambda5 = k51*y[3]+k52*y[7]+k53*y[11]+k54*y[15]+params[19]*y[23]+params[20]*y[28];
    lambda6 = k61*y[3]+k62*y[7]+k63*y[11]+k64*y[15]+k65*y[23]+params[21]*y[28];
    IR5 = y[21]*lambda5/1e+05;
    IR6 = y[26]*lambda6/1e+05;
    dydt[1] = -IR1;
    dydt[2] = IR1-InR1;
    dydt[3] = InR1-RR1;
    dydt[4] = RR1;
    dydt[5] = -IR2;
    dydt[6] = IR2-InR2;
    dydt[7] = InR2-RR2;
    dydt[8] = RR2;
    dydt[9] = -IR3;
    dydt[10] = IR3-InR3;
    dydt[11] = InR3-RR3;
    dydt[12] = RR3;
    dydt[13] = -IR4;
    dydt[14] = IR4-InR4;
    dydt[15] = InR4-RR4;
    dydt[16] = RR4;
    dydt[17] = RI1;
    dydt[18] = RI2;
    dydt[19] = RI3;
    dydt[20] = RI4;
    dydt[21] = -IR5;
    dydt[22] = IR5-InR5;
    dydt[23] = InR5-RR5;
    dydt[24] = RR5;
    dydt[25] = RI5;
    dydt[26] = -IR6;
    dydt[27] = IR6-InR6;
    dydt[28] = InR6-RR6;
    dydt[29] = RR6;
    dydt[30] = RI6;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs; // Number of days sampled
  int<lower = 1> n_params; // Number of model parameters
  int<lower = 1> n_difeq; // Number of differential equations in the system
  int y1[n_obs];
  int y2[n_obs];
  int y3[n_obs];
  int y4[n_obs];
  int y5[n_obs];
  int y6[n_obs];
  real t0; // Initial time point (zero)
  real ts[n_obs]; // Time points that were sampled
  vector[n_difeq] y0;
}
parameters {
  real<lower = 0> params[n_params]; // Model parameters
  real<lower = 0, upper = 1> rho;
}
transformed parameters{
  vector[n_difeq] y_hat[n_obs]; // Output from the ODE solver
  real incidence1[n_obs];
  real incidence2[n_obs];
  real incidence3[n_obs];
  real incidence4[n_obs];
  real incidence5[n_obs];
  real incidence6[n_obs];
  y_hat = ode_rk45(SEIR_matrix_sym, y0, t0, ts, params);
  incidence1[1] =  rho * (y_hat[1, 17]  - y0[17]);
  incidence2[1] =  rho * (y_hat[1, 18]  - y0[18]);
  incidence3[1] =  rho * (y_hat[1, 19]  - y0[19]);
  incidence4[1] =  rho * (y_hat[1, 20]  - y0[20]);
  incidence5[1] =  rho * (y_hat[1, 25]  - y0[25]);
  incidence6[1] =  rho * (y_hat[1, 30]  - y0[30]);
  for (i in 1:n_obs-1) {
    incidence1[i + 1] = rho * (y_hat[i + 1, 17] - y_hat[i, 17] + 1e-5);
    incidence2[i + 1] = rho * (y_hat[i + 1, 18] - y_hat[i, 18] + 1e-5);
    incidence3[i + 1] = rho * (y_hat[i + 1, 19] - y_hat[i, 19] + 1e-5);
    incidence4[i + 1] = rho * (y_hat[i + 1, 20] - y_hat[i, 20] + 1e-5);
    incidence5[i + 1] = rho * (y_hat[i + 1, 25] - y_hat[i, 25] + 1e-5);
    incidence6[i + 1] = rho * (y_hat[i + 1, 30] - y_hat[i, 30] + 1e-5);
   }
}
model {
    params ~ normal(0, 10);
    rho    ~ normal(0.5, 0.5);
    y1     ~ poisson(incidence1);
    y2     ~ poisson(incidence2);
    y3     ~ poisson(incidence3);
    y4     ~ poisson(incidence4);
    y5     ~ poisson(incidence5);
    y6     ~ poisson(incidence6);
}
generated quantities {
  real log_lik;
  log_lik = poisson_lpmf(y1 | incidence1) + poisson_lpmf(y2 | incidence2) + poisson_lpmf(y3 | incidence3) + poisson_lpmf(y4 | incidence4) + poisson_lpmf(y5 | incidence5) + poisson_lpmf(y6 | incidence6);
}
