#' Title 计算由hc聚类所产生小群落间的平均配对功能性状距离(SES.贝塔MPD)
#' @description 计算由hc聚类所产生小群落间的平均配对功能性状距离(SES.贝塔MPD)
#' @usage hc_m_community_Beta_MPD_pair(entity,traits,minx,miny,maxx,maxy,scale,replicate,edge)
#'
#' @param entity 植物个体的分布坐标与物种
#' @param traits 物种的功能性状
#' @param minx 核心区边界，最小x
#' @param miny 核心区边界，最小y
#' @param maxx 核心区边界，最大x
#' @param maxy 核心区边界，最大y
#' @param scale hc聚类的scale参数
#' @param replicate 零模型重复次数
#' @param edge 基于边界筛选小群落的方法，可选"inside"与"centroid"
#'
#' @return SES.贝塔MPD
#' @export
#'
#' @examples
#' entity=M.COMMUNITIES::entity
#' traits=M.COMMUNITIES::traits
#' hc_m_community_Beta_MPD_pair(entity,traits,minx=5,miny=5,maxx=45,maxy=45,scale=5,replicate=99,edge="inside")
hc_m_community_Beta_MPD_pair=
function(entity,traits,minx,miny,maxx,maxy,scale,replicate,edge)
{
  library(FD)
  library(picante)
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

  ######################
  ######################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ######################
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

  ######################
  ######################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ######################
  abundance=abundance_hc_m_communities(entity,minx,miny,maxx,maxy,scale,edge)
  abundance=abundance[,-ncol(abundance)]
  traits=traits[which(rownames(traits) %in% colnames(abundance)),]
  trait_dist=gowdis(traits)
  trait_dist=as.matrix(trait_dist)
  Beta_MPD=comdist(abundance, trait_dist, abundance.weighted=TRUE)
  Beta_MPD=as.matrix(Beta_MPD)
  Beta_MPD=round(Beta_MPD,3)
  diag(Beta_MPD)=NA
  Beta_MPD_obv=Beta_MPD
  prob=colSums(abundance)
  prob=as.numeric(prob)
  size=rowSums(abundance)
  size=as.numeric(size)

  null_list=list()
  for (i in 1:replicate)
  {
    null_comm = matrix(0, nrow = length(prob), ncol = length(size))
    for (j in 1:length(size)) {

      null_comm[, j] = rmultinom(1, size = size[j], prob = prob)
      rownames(null_comm)=colnames(abundance)
      colnames(null_comm)=rownames(abundance)
    }
    null_name = paste("m_com", i)
    null_list[[null_name]]=null_comm
    null_name = paste("m_com", i)
    null_comm=t(null_comm)
    null_comm=comdist(null_comm, trait_dist, abundance.weighted=TRUE)
    null_comm=as.matrix(null_comm)
    diag(null_comm)=NA
    null_list[[null_name]]=null_comm
  }
  null_Beta_MPD_mean=apply(simplify2array(null_list), MARGIN = c(1, 2), FUN = mean)
  null_Beta_MPD_mean=round(null_Beta_MPD_mean,3)
  null_Beta_MPD_std=apply(simplify2array(null_list), MARGIN = c(1, 2), FUN = sd)
  null_Beta_MPD_std=round(null_Beta_MPD_std,3)
  SES_Beta_MPD=(Beta_MPD_obv-null_Beta_MPD_mean)/null_Beta_MPD_std
  SES_Beta_MPD=round(SES_Beta_MPD,3)
  result=list(SES_Beta_MPD=SES_Beta_MPD,Beta_MPD_obv=Beta_MPD_obv,null_Beta_MPD_mean=null_Beta_MPD_mean,null_Beta_MPD_std=null_Beta_MPD_std)
  result
}
