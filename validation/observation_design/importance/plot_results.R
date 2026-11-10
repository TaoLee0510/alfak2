#!/usr/bin/env Rscript
out<-normalizePath(Sys.getenv("ALFAK2_IS_OUTPUT"));m<-read.csv(file.path(out,"mechanism/overall_metrics.csv"))
m<-m[m$distance>=0,];arms<-c("A_finite_unknown","B_mean_unknown","C_finite_oracle","D_mean_oracle")
labels<-c("Finite / unknown pi0","Mean / unknown pi0","Finite / oracle pi0","Mean / oracle pi0")
colors<-c("#0072B2","#D55E00","#009E73","#CC79A7");types<-c(1,1,2,2)
draw<-function() {
 par(mfrow=c(2,3),mar=c(4,4.3,2,1),oma=c(4.5,1,4,1))
 keys<-c("rmse","shape_rmse","edge_rmse","score95","coverage95","width95")
 titles<-c("Relative fitness RMSE","Globally centered Shape RMSE","Edge RMSE","95% interval score","95% interval coverage","95% interval width")
 for(k in seq_along(keys)) {
  key<-keys[k];yl<-range(c(0,m[[key]],if(key=="coverage95")1 else numeric()))
  plot(range(m$distance),yl,type="n",xlab="Original-control distance",ylab=titles[k],xaxt="n")
  axis(1,sort(unique(m$distance)))
  if(key=="coverage95")abline(h=.95,col="grey50",lty=3)
  for(j in seq_along(arms)) {a<-m[m$design==arms[j],];a<-a[order(a$distance),];lines(a$distance,a[[key]],type="b",col=colors[j],lty=types[j],pch=14+j,lwd=1.8)}
 }
 mtext("Fixed-layout development mechanism study: six landscapes; family-balanced means",side=1,outer=TRUE,line=2,cex=.85)
 mtext("Oracle arms are diagnostic and are excluded from deployable-method rankings",side=1,outer=TRUE,line=3.3,cex=.8)
 par(oma=c(0,0,0,0),fig=c(0,1,.94,1),new=TRUE,mar=c(0,0,0,0));plot.new()
 legend("center",legend=labels,col=colors,lty=types,pch=15:18,lwd=1.8,horiz=TRUE,bty="n",cex=.82)
}
pdf(file.path(out,"mechanism/mechanism_summary.pdf"),width=12,height=8);draw();dev.off()
png(file.path(out,"mechanism/mechanism_summary.png"),width=2400,height=1600,res=200);draw();dev.off()
