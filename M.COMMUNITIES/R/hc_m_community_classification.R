#' 基于hc聚类对小群落进行划分 / Divide mirror communities based on HC clustering
#' @description 基于hc聚类对小群落进行划分 / Divide mirror communities based on HC clustering
#' @usage hc_m_community_classification(entity,minx,miny,maxx,maxy,scale,edge)
#' @param entity 植物个体的分布坐标与物种 / Distribution coordinates and species of plant individuals
#' @param minx 核心区边界，最小x / Core area boundary, minimum x
#' @param miny 核心区边界，最小y / Core area boundary, minimum y
#' @param maxx 核心区边界，最大x / Core area boundary, maximum x
#' @param maxy 核心区边界，最大y / Core area boundary, maximum y
#' @param scale hc聚类的scale参数 / scale parameter of HC clustering
#' @param edge 基于边界筛选小群落的方法，可选"inside"与"centroid" / Method for screening mirror communities based on boundaries, with optional "inside" and "centroid"
#'
#' @return 小群落的划分 / Division of mirror communities
#' @export
#'
#' @examples
#' entity=M.COMMUNITIES::entity
#' hc_m_community_classification(entity,minx=5,miny=5,maxx=45,maxy=45,scale=5,edge="centroid")
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


