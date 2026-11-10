#include <TMB.hpp>
#include "positive_taylor_forward.h"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_MATRIX(y);
  DATA_IVECTOR(trajectory);
  DATA_VECTOR(times);
  DATA_IVECTOR(steps);
  DATA_MATRIX(eta_basis);
  DATA_MATRIX(design);
  DATA_VECTOR(offset);
  DATA_SPARSE_MATRIX(precision);
  DATA_SCALAR(logdet);
  DATA_VECTOR(coef_sd);
  DATA_IVECTOR(trans_from);
  DATA_IVECTOR(trans_to);
  DATA_VECTOR(trans_weight);
  PARAMETER_MATRIX(zeta);
  PARAMETER_VECTOR(u);
  PARAMETER_VECTOR(coefficients);
  const int n=u.size(), r=zeta.cols(), m=y.cols();
  Type nll=0;
  vector<Type> f_latent=offset+u;
  for(int i=0;i<n;i++) for(int j=0;j<coefficients.size();j++)
    f_latent(i)+=design(i,j)*coefficients(j);
  vector<Type> f_report(n);
  for(int i=0;i<n;i++) f_report(i)=exp(f_latent(i));
  matrix<Type> pi0(n,r), predicted(n,m);
  for(int k=0;k<r;k++) {
    vector<Type> eta(n);eta.setZero();
    for(int i=0;i<n;i++) for(int j=0;j<zeta.rows();j++) eta(i)+=eta_basis(i,j)*zeta(j,k);
    Type mx=eta.maxCoeff(),total=0;
    for(int i=0;i<n;i++){pi0(i,k)=exp(eta(i)-mx);total+=pi0(i,k);}
    for(int i=0;i<n;i++) pi0(i,k)/=total;
  }
  vector<Type> log_growth(m);log_growth.setZero();
  for(int o=0;o<m;o++) {
    vector<Type> p0=pi0.col(trajectory(o));
    if(asDouble(times(o))==0) predicted.col(o)=p0;
    else {
      DivisionForwardResult<Type> out=division_forward(p0,f_report,trans_from,trans_to,trans_weight,times(o),steps(o));
      predicted.col(o)=out.frequency;log_growth(o)=out.log_total;
    }
    for(int i=0;i<n;i++) nll-=y(i,o)*log(predicted(i,o)+Type(1e-14));
  }
  vector<Type> Qu=precision*u;
  nll+=Type(.5)*(u*Qu).sum()-Type(.5)*logdet+Type(.5*n)*log(Type(2*M_PI));
  for(int j=0;j<coefficients.size();j++) nll-=dnorm(coefficients(j),Type(0),coef_sd(j),true);
  for(int j=0;j<zeta.size();j++) nll-=dnorm(zeta(j),Type(0),Type(1),true);
  ADREPORT(f_latent);ADREPORT(f_report);
  REPORT(pi0);REPORT(predicted);REPORT(log_growth);REPORT(f_latent);REPORT(f_report);
  return nll;
}
