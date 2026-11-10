#ifndef OBSERVATION_POSITIVE_TAYLOR_FORWARD_H
#define OBSERVATION_POSITIVE_TAYLOR_FORWARD_H

// Positive Taylor action for the complete Metzler generator. Subtracting its
// minimum diagonal leaves a nonnegative matrix while factoring common growth
// out exactly. Order eight reduces the step refinements (and nested-AD tape)
// required by the former second-order splitting at large common birth scales.
template<class Type>
struct DivisionForwardResult {
  vector<Type> frequency;
  Type log_total;
};

template<class Type>
DivisionForwardResult<Type> division_forward(vector<Type> initial,vector<Type> birth,
    vector<int> from,vector<int> to,vector<Type> probability,Type dt,int steps) {
  const int n=initial.size();
  vector<Type> diagonal(n);
  matrix<Type> off(n,n);off.setZero();
  for(int i=0;i<n;i++)diagonal(i)=-birth(i);
  for(int e=0;e<probability.size();e++) {
    if(from(e)==to(e))diagonal(from(e))+=Type(2)*probability(e)*birth(from(e));
    else off(to(e),from(e))+=Type(2)*probability(e);
  }
  Type shift=diagonal(0);
  for(int i=1;i<n;i++)shift=CppAD::CondExpLt(diagonal(i),shift,diagonal(i),shift);
  Type h=dt/Type(steps),log_total=0;
  vector<Type> state=initial;
  for(int s=0;s<steps;s++) {
    vector<Type> term=state,accum=state;
    for(int order=1;order<=8;order++) {
      matrix<Type> rhs(n,1);
      for(int i=0;i<n;i++)rhs(i,0)=birth(i)*term(i);
      matrix<Type> product=atomic::matmul(off,rhs);
      vector<Type> next(n);
      for(int i=0;i<n;i++)next(i)=(product(i,0)+(diagonal(i)-shift)*term(i))*h/Type(order);
      accum+=next;term=next;
    }
    state=accum;log_total+=h*shift;
    if((s+1)%16==0||s+1==steps) {
      Type total=state.sum();log_total+=log(total);state/=total;
    }
  }
  DivisionForwardResult<Type> answer;
  answer.frequency=state;answer.log_total=log_total;return answer;
}
#endif
