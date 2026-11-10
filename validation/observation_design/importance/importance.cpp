#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

static std::vector<double> square_product(const std::vector<double>& A,const std::vector<double>& B,int n) {
  std::vector<double> C(n*n,0);
  for(int i=0;i<n;i++)for(int k=0;k<n;k++)if(A[i*n+k]!=0)
    for(int j=0;j<n;j++)C[i*n+j]+=A[i*n+k]*B[k*n+j];
  return C;
}

static std::vector<double> mean_action(const NumericMatrix& P,const std::vector<double>& birth,
    const std::vector<double>& initial,double dt,bool beta_zero,double& growth) {
  int n=birth.size();std::vector<double> state(initial),B(n*n),term(n),acc(n),next(n);
  growth=0;if(dt==0)return state;
  if(beta_zero) {
    double mx=R_NegInf;
    for(int i=0;i<n;i++)if(initial[i]>0)mx=std::max(mx,std::log(initial[i])+dt*birth[i]);
    double total=0;
    for(int i=0;i<n;i++){state[i]=initial[i]>0?std::exp(std::log(initial[i])+dt*birth[i]-mx):0;total+=state[i];}
    if(!std::isfinite(total)||!(total>0))stop("Nonfinite analytic mean action");
    growth=mx+std::log(total);for(double& v:state)v/=total;return state;
  }
  double shift=R_PosInf,norm=0;
  for(int i=0;i<n;i++)shift=std::min(shift,(2*P(i,i)-1)*birth[i]);
  for(int from=0;from<n;from++) {
    double column=0;
    for(int to=0;to<n;to++) {
      double x=2*P(from,to)*birth[from];if(from==to)x-=birth[from]+shift;
      if(x< -1e-10*std::max(1.0,birth[from]))stop("Mean generator is not Metzler");
      if(x<0)x=0; // roundoff of an analytically zero shifted diagonal only
      B[to*n+from]=x;column+=x;
    }
    norm=std::max(norm,column);
  }
  double extent=dt*norm;
  if(!std::isfinite(extent))stop("Nonfinite finite-proposal mean generator");
  if(extent<=8) {
    int blocks=std::max(1,(int)std::ceil(extent));double h=dt/blocks;
    for(int s=0;s<blocks;s++) {
      term=acc=state;bool converged=false;
      for(int order=1;order<=64;order++) {
        std::fill(next.begin(),next.end(),0);double maximum=0,scale=0;
        for(int to=0;to<n;to++) {
          for(int from=0;from<n;from++)next[to]+=B[to*n+from]*term[from]*h/order;
          acc[to]+=next[to];maximum=std::max(maximum,std::abs(next[to]));scale=std::max(scale,std::abs(acc[to]));
        }
        term=next;
        if(maximum<=1e-15*std::max(1.0,scale)){converged=true;break;}
      }
      if(!converged)stop("Independent mean series did not converge");
      double total=0;for(double x:acc)total+=x;
      if(!std::isfinite(total)||!(total>0))stop("Nonfinite mean action");
      growth+=h*shift+std::log(total);for(int i=0;i<n;i++)state[i]=acc[i]/total;
    }
    return state;
  }
  // Positive scaling/squaring avoids an exponential number of blocks for
  // heavy-tailed proposal points; every finite proposal remains evaluated.
  int scaling=std::max(0,(int)std::ceil(std::log2(extent)));
  double h=std::ldexp(dt,-scaling);std::vector<double> H=B,M(n*n,0),T(n*n,0);
  for(double& x:H)x*=h;for(int i=0;i<n;i++)M[i*n+i]=T[i*n+i]=1;
  bool converged=false;
  for(int order=1;order<=64;order++) {
    T=square_product(H,T,n);double maximum=0,scale=0;
    for(int j=0;j<n*n;j++){T[j]/=order;M[j]+=T[j];maximum=std::max(maximum,std::abs(T[j]));scale=std::max(scale,std::abs(M[j]));}
    if(maximum<=1e-15*std::max(1.0,scale)){converged=true;break;}
  }
  if(!converged)stop("Scaled mean matrix series did not converge");
  double normalizer=*std::max_element(M.begin(),M.end()),log_scale=std::log(normalizer);
  for(double& x:M)x/=normalizer;
  for(int j=0;j<scaling;j++) {
    M=square_product(M,M,n);normalizer=*std::max_element(M.begin(),M.end());
    if(!std::isfinite(normalizer)||!(normalizer>0))stop("Scaled matrix squaring failed");
    for(double& x:M)x/=normalizer;log_scale=2*log_scale+std::log(normalizer);
  }
  double total=0;std::fill(state.begin(),state.end(),0);
  for(int to=0;to<n;to++){for(int from=0;from<n;from++)state[to]+=M[to*n+from]*initial[from];total+=state[to];}
  if(!std::isfinite(total)||!(total>0))stop("Scaled mean action failed");
  for(double& x:state)x/=total;growth=dt*shift+log_scale+std::log(total);return state;
}

// [[Rcpp::export]]
List is_forward_cpp(NumericVector pi0,NumericVector birth,NumericMatrix probability,double dt,bool beta_zero=false) {
  std::vector<double> p(pi0.begin(),pi0.end()),b(birth.begin(),birth.end());double growth=0;
  auto q=mean_action(probability,b,p,dt,beta_zero,growth);
  return List::create(_["pi1"]=wrap(q),_["log_growth"]=growth);
}

// [[Rcpp::export]]
NumericVector is_loglikelihood_cpp(NumericMatrix draws,List specification) {
  int n=as<int>(specification["n"]),q=as<int>(specification["q"]),N=draws.nrow();
  if(draws.ncol()!=q)stop("Importance coordinates do not match target");
  NumericVector offset=specification["offset"],y0=specification["y0"],y1=specification["y1"],known=specification["known_pi0"];
  NumericMatrix L=specification["lower"],E=specification["eta_basis"],P=specification["probability"];
  double dt=as<double>(specification["dt"]);bool zero=as<bool>(specification["beta_zero"]);
  NumericVector result(N);std::vector<double> birth(n),p(n),eta(n);
  for(int k=0;k<N;k++) {
    if(k%1000==0)checkUserInterrupt();
    for(int i=0;i<n;i++) {
      double g=offset[i];for(int j=0;j<n;j++)g+=L(i,j)*draws(k,j);
      birth[i]=std::exp(g);if(!std::isfinite(birth[i]))stop("Nonfinite proposal birth: no sample omitted");
    }
    if(known.size()==n)for(int i=0;i<n;i++)p[i]=known[i];
    else {
      for(int i=0;i<n;i++){eta[i]=0;for(int j=0;j<q-n;j++)eta[i]+=E(i,j)*draws(k,n+j);}
      double mx=*std::max_element(eta.begin(),eta.end()),total=0;
      for(int i=0;i<n;i++)total+=(p[i]=std::exp(eta[i]-mx));
      for(double& x:p)x/=total;
    }
    double growth=0;auto end=mean_action(P,birth,p,dt,zero,growth);double ll=0;
    for(int i=0;i<n;i++)ll+=y0[i]*std::log(p[i]+1e-14)+y1[i]*std::log(end[i]+1e-14);
    if(!std::isfinite(ll))stop("Nonfinite proposal likelihood: no sample omitted");
    result[k]=ll;
  }
  return result;
}
