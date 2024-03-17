#' 计算由hc聚类所产生小群落的质心
#' @description 计算由hc聚类所产生小群落的质心
#' @usage centroid_hc_m_community(entity,minx,miny,maxx,maxy,scale,edge)
#' @param entity 植物个体的分布坐标与物种
#' @param minx 核心区边界，最小x
#' @param miny 核心区边界，最小y
#' @param maxx 核心区边界，最大x
#' @param maxy 核心区边界，最大y
#' @param scale hc聚类的尺度
#' @param edge 基于边界筛选小群落的方法，可选"inside"与"centroid"
#'
#' @return 各个小群落的质心坐标
#' @export
#'
#' @examples
#' entity=M.COMMUNITIES::entity
#' centroid_hc_m_community(entity,minx=5,miny=5,maxx=45,maxy=45,scale=5,edge="inside")
centroid_hc_m_community=function(entity,minx,miny,maxx,maxy,scale,edge)
{
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

  ###########################
  ###########################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###########################
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
  ########################
  ########################~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ########################
  abundance_hc_m_communities=function(entity,minx,miny,maxx,maxy,scale,edge)
  {
    library(FD)
    #####################

    #####################
    area=(max(entity[,1])-min(entity[,1]))*(max(entity[,2])-min(entity[,2]))
    n_population=nrow(entity)
    D=n_population/area
    if (missing(scale)) {
      scale= sqrt(D)
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

    ################################
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
    ##################
    class=hc_m_community_classification(entity,minx,miny,maxx,maxy,scale,edge)
    ###################
    spe=levels(as.factor(entity$species))
    abundance=rep(NA,length(spe)+1)
    abundance_i=numeric()
    n_species=length(spe)
    ##################
    n=max(class$m_community_class)
    for (i in 1:n)
    {
      m_community_i=subset(class,class$m_community_class==i)

      for (j in 1:n_species)
      {
        m_community_i_species=subset(m_community_i,m_community_i$species==spe[j])
        abundance_i[j] =nrow(m_community_i_species)
        abundance_i[j+1]=i

      }
      abundance=rbind(abundance,abundance_i)
      colnames(abundance)=c(spe,"m_community_ID")
    }

    abundance=abundance[-1,]
    rownames(abundance)=paste("m_community",1:nrow(abundance))
    abundance
  }
  ########################
  ########################~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ########################
  abundance=abundance_hc_m_communities(entity,minx,miny,maxx,maxy,scale,edge)
  m_community_class=hc_m_community_classification(entity,minx,miny,maxx,maxy,scale,edge)
  centroid_x=numeric()
  centroid_y=numeric()
  for (i in 1:max(m_community_class$m_community_class))
  {
    m_community_i=subset(m_community_class,m_community_class==i)
    centroid_x[i]=mean(m_community_i$x)
    centroid_y[i]=mean(m_community_i$y)
  }
  centroid=cbind(centroid_x,centroid_y)
  result=cbind(abundance,centroid)
  result
}
