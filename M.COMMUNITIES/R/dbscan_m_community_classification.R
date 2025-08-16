#' 基于dbscan聚类对小群落进行划分 / Divide mirror communities based on DBSCAN clustering
#' @description 基于dbscan聚类对小群落进行划分
#' @usage dbscan_m_community_classification(entity,minx,miny,maxx,maxy,eps,minPts,edge)
#' @param entity 植物个体的分布坐标与物种
#' @param minx 核心区边界，最小x
#' @param miny 核心区边界，最小y
#' @param maxx 核心区边界，最大x
#' @param maxy 核心区边界，最大y
#' @param eps dbscan聚类的eps参数
#' @param minPts dbscan聚类的minPts参数
#' @param edge 基于边界筛选小群落的方法，可选"inside"与"centroid"
#'
#' @return 小群落的划分
#' @export
#'
#' @examples
#' entity=M.COMMUNITIES::entity
#' dbscan_m_community_classification(entity,minx=5,miny=5,maxx=45,maxy=45,eps=2,minPts=1,edge="centroid")
dbscan_m_community_classification=function(entity,minx,miny,maxx,maxy,eps,minPts,edge)
{
  area=(max(entity[,1])-min(entity[,1]))*(max(entity[,2])-min(entity[,2]))
  n_population=nrow(entity)
  D=n_population/area
  if (missing(eps) & missing(minPts)) {
    eps = 1/2*sqrt(D)
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
  library(dbscan)
  m_community= dbscan(entity[,1:2], eps = eps, minPts =minPts)
  result=cbind(entity,m_community$cluster)
  result=result[order(result[,4]),]
  colnames(result)=c("x","y","species","m_community_class")
  result
  result_1=subset(result,result$m_community_class==0)
  if (!nrow(result_1)==0){
    result_1$m_community_class=1:nrow(result_1)
    result_2=subset(result,!result$m_community_class==0)
    result_2$m_community_class=result_2$m_community_class+max(result_1$m_community_class)
    result=rbind(result_1,result_2)}
  ####
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

