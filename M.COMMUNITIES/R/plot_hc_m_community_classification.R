#' 绘制基于hc聚类划分的小群落的分布图
#' @description 绘制基于hc聚类划分的小群落的分布图
#' @usage plot_hc_m_community_classification(entity,minx,miny,maxx,maxy,scale,edge)
#' @param entity 植物个体的分布坐标与物种
#' @param minx 核心区边界，最小x
#' @param miny 核心区边界，最小y
#' @param maxx 核心区边界，最大x
#' @param maxy 核心区边界，最大y
#' @param scale hc聚类的scale参数
#' @param edge 基于边界筛选小群落的方法，可选"inside"与"centroid"
#'
#' @return 小群落空间分布图
#' @export
#'
#' @examples
#' entity=M.COMMUNITIES::entity
#' plot_hc_m_community_classification(entity,minx=5,miny=5,maxx=45,maxy=45,scale=5,edge="inside")
plot_hc_m_community_classification=function(entity,minx,miny,maxx,maxy,scale,edge)
{
library(ggplot2)
entity=data.frame(entity)

hc_m_community_classification=function(entity,minx,miny,maxx,maxy,scale,edge)
{
  library(stats)
  area=(max(entity[,1])-min(entity[,1]))*(max(entity[,2])-min(entity[,2]))
  n_population=nrow(entity)
  D=n_population/area
  if (missing(scale)) {
    scale= sqrt(D)
    minPts=2
  }
  if (missing(minx))
  {
    minx=min(entity[,1])
  }
  if (missing(maxx))
  {
    maxx=max(entity[,1])
  }

  if (missing(miny))
  {
    miny=min(entity[,2])
  }
  if (missing(maxy))
  {
    maxy=max(entity[,2])
  }
  dist_mat=dist(entity[,1:2])

  hc=hclust(dist_mat,method = "complete")
  m_community_class=cutree(hc, h = scale)
  result=cbind(entity,m_community_class)
  result=result[order(result[,4]),]
  result

  class=c(NA)

  for(i in 1:max(result$m_community_class))
  {
    result_i=result[which(result$m_community_class==i),]
    if(edge=="inside")
    {    if(
      length(which(result_i[,1]>minx & result_i[,1]<maxx & result_i[,2]>miny & result_i[,2]<maxy))==0
    )
    {class=c(class,i)}
    }
    if(edge=="centroid")
    {
      centroid=apply(result_i[,1:2],2,mean)
      if(
        length(which(centroid[1]>minx & centroid[1]<maxx & centroid[2]>miny & centroid[2]<maxy))==0
      )
      {class=c(class,i)}
    }

  }

  class=class[-1]
  result=result[! result$m_community_class %in% class,]
  result$m_community_class=match(result$m_community_class, sort(unique(result$m_community_class)))
  result
}


point=hc_m_community_classification(entity,minx,miny,maxx,maxy,scale,edge)
p=ggplot(point,aes(x=x,y=y))+
  geom_label(aes(fill=factor(m_community_class),label=m_community_class),color="white",fontface="bold",size=2,nudge_x = 0.5,nudge_y=-1)+
  geom_point(data=entity,aes(x=x,y=y),color="grey",size=0.75)+
  geom_point(aes(color=factor(m_community_class)),alpha=0.5,size=0.75)+
  theme_bw()+
  theme(legend.position="none")+
  scale_x_continuous(expand= c(0, 0))+
  scale_y_continuous(expand= c(0, 0))+
  geom_vline(xintercept =c(minx,maxx),linetype=2)+
  geom_hline(yintercept =c(miny,maxy),linetype=2)
p
}
