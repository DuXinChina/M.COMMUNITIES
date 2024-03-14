#' Title 基于点格局分析的方法，判定以hc聚类划分小群落时的合理聚类尺度
#' @description 基于点格局分析的方法，判定以hc聚类划分小群落时的合理聚类尺度
#' @usage hc_m_community_scale_select(entity,minx,miny,maxx,maxy,window,two_nested_clu=FALSE)
#'
#' @param entity 植物个体的分布坐标与物种
#' @param minx 研究范围的最小x
#' @param miny 研究范围的最小y
#' @param maxx 研究范围的最大x
#' @param maxy 研究范围的最大y
#' @param window 点格局分析时的窗口大小
#' @param two_nested_clu 逻辑变量，可为T或F，判定是否模拟双尺度镶嵌托马斯过程
#'
#' @return 聚类尺度
#' @export
#'
#' @examples
#' entity=M.COMMUNITIES::entity
#' hc_m_community_scale_select(entity,minx=0,miny=0,maxx=50,maxy=50,window=35,two_nested_clu=T)
hc_m_community_scale_select=function(entity,minx,miny,maxx,maxy,window,two_nested_clu=FALSE)
{
  library(minpack.lm)
  library(spatstat)
  library(zoo)
  max_dist=max(dist(entity[,1:2]))*0.5
  window_=window
  if(window>max_dist)
  {
    window_= max_dist
  }
  point=ppp(entity[,1],entity[,2],window=owin(c(minx,maxx),c(miny,maxy)))
  G_r=envelope(point,pcf,nsim=99,r=seq(0,max_dist,length.out=200),correction="Ripley")
  L_r=envelope(point,Lest,nsim=99,r=seq(0,max_dist,length.out=200),correction="Ripley")

  par(mfrow = c(2, 2))
  plot(L_r,.-theo~r,legend = FALSE,ylab=expression(italic("L(r)")),xlim=c(0,window_))
  plot(G_r, legend = FALSE,xlim=c(0,window_))
  density=density(dist(entity[,1:2]),kernel ="triangular")
  x=runif(nrow(entity),minx,maxx)
  y=runif(nrow(entity),miny,maxy)
  point=cbind(x,y)
  dist=dist(point)
  density_ran=density(dist,kernel="triangular")

  density1=density$y/(density$x)^2
  plot(density,main="Pair point dist density",xlab=expression(italic("r")),ylab=expression(italic("Density")),xlim=c(0,window_),ylim=c(0,1.1*max(density$y,density_ran$y)))
  for (i in 1:99)
  {
    x=runif(nrow(entity),minx,maxx)
    y=runif(nrow(entity),miny,maxy)
    point=cbind(x,y)
    dist=dist(point)
    density_i=density(dist,kernel="triangular")
    lines(density_i, col = "grey")
  }
  lines(density, col = "black")
  if(missing(two_nested_clu))
  {
    two_nested_clu=FALSE
  } else {
    two_nested_clu=TRUE
  }
  if (two_nested_clu==TRUE) {
    point=ppp(entity[,1],entity[,2],window=owin(c(minx,maxx),c(miny,maxy)))
    window=window

    pcf_=pcf(point,r=seq(0,window,length.out=200),correction="Ripley")
    r=pcf_$r
    pcf_=pcf_$iso
    pcf_=cbind(pcf_,r)
    pcf_=pcf_[-1,]
    r=r[-1]
    x=pcf_[,2]
    y=pcf_[,1]
    fun=function(x,sigma_s,sigma_l,ps,pl){1+exp(-0.5*x^2/(2*sigma_s^2))/(4*pi*sigma_s^2*ps)+exp(-0.5*x^2/(2*sigma_s^2+2*sigma_l^2))/(2*pi*(2*sigma_s^2+2*sigma_l^2)*pl)}
    nls.out <- nlsLM(y~fun(x,sigma_s,sigma_l,ps,pl),start = list(sigma_s =0.1,sigma_l=1,ps=0.1,pl=0.01))
    result=summary(nls.out)

    #result=summary(nls)
    sigma_s=result$parameters[1,1]
    sigma_l=result$parameters[2,1]
    pl=result$parameters[4,1]
    ps=result$parameters[3,1]

    print("OK")
    fit_y=1+exp(-0.5*pcf_[,2]^2/(2*sigma_s^2))/(4*pi*sigma_s^2*ps)+exp(-0.5*pcf_[,2]^2/(2*sigma_s^2+2*sigma_l^2))/(2*pi*(2*sigma_s^2+2*sigma_l^2)*pl)
    plot(pcf_[,2],pcf_[,1],typ="l",xlab=expression(italic("r")),ylab=expression(italic("g(r)")),main="Two nested scales of clustering")
    abline(h=1,col="red",lty=2)
    lines(r,fit_y,lty=2,col="blue")
    square_R=1-sum((fit_y-pcf_[,1])^2)/sum((pcf_[,1]-mean(pcf_[,1]))^2)
    sigma_s=round(sigma_s,6)
    sigma_l=round(sigma_l,6)
    pl=round(pl,6)
    ps=round(ps,6)
    scale_s=1.96*2*sigma_s
    scale_l=1.96*2*sigma_l
    square_R=round(square_R,3)
    text=paste0("Sigma_s:  ",sigma_s,"\n",
                "Sigma_l:  ",sigma_l,"\n",
                "ps:  ",ps,"\n",
                "pl:  ",pl,"\n",
                "scale_s:  ",scale_s,"\n",
                "scale_l:  ",scale_l,"\n",
                "square_R:  ",square_R)
    text(0.8*mean(r),0.7*max(pcf_[,1]),text,adj=c(0))
  }
  par(mfrow = c(1, 1))
}
