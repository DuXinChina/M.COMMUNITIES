#' 统计有hc聚类形成小群落中的物种丰度排序信息
#' @description 统计有hc聚类形成小群落中的物种丰度排序信息
#' @usage stat_hc_m_community_abundance_composition(entity,minx,miny,maxx,maxy,scale,edge)
#' @param entity 植物个体的分布坐标与物种
#' @param minx 核心区边界，最小x
#' @param miny 核心区边界，最小y
#' @param maxx 核心区边界，最大x
#' @param maxy 核心区边界，最大y
#' @param scale hc聚类的scale参数
#' @param edge 基于边界筛选小群落的方法，可选"inside"与"centroid"
#'
#' @return 各类物种丰度相关信息
#' @export
#'
#' @examples
#' entity=M.COMMUNITIES::entity
#' stat_hc_m_community_abundance_composition(entity,minx=5,miny=5,maxx=45,maxy=45,scale=5,edge="centroid")
stat_hc_m_community_abundance_composition=function(entity,minx,miny,maxx,maxy,scale,edge)
{
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


  C11=abundance_hc_m_communities(entity,minx,miny,maxx,maxy,scale,edge)
  C11=C11[,-ncol(C11)]
  sp_n=NA
  for (i in 1:nrow(C11))
  {
    Ci=C11[i,]
    sp_ni=subset(Ci,Ci>0)
    sp_ni=length(sp_ni)
    sp_n=c(sp_n,sp_ni)
  }
  sp_n=sp_n[-1]
  max_sp_n=max(sp_n)
  abun=rep(NA,max_sp_n)
  for (i in 1:nrow(C11))
  {
    Ci=C11[i,]
    sp_ni=subset(Ci,Ci>0)
    abun_i=sp_ni
    abun_i=abun_i[order(abun_i,decreasing=T)]
    sp_ni=length(sp_ni)
    abun_i=c(abun_i,rep(0,(max_sp_n-sp_ni)))
    abun=rbind(abun,abun_i)
  }
  colnames(abun)=1:max_sp_n
  abun=abun[-1,]
  total_abun=apply(abun,1,sum)
  mean_total_abun=mean(total_abun)
  sd_total_abun=sd(total_abun)
  mean_abun=colMeans(abun)
  sd_abun=apply(abun,2,sd)
  mean_sp_n=mean(sp_n)
  sd_sp_n=sd(sp_n)
  library(ggplot2)
  mean_abun_=cbind(1:max_sp_n,mean_abun)
  mean_abun_=data.frame(mean_abun_)
  colnames(mean_abun_)=c("x","mean_abun")
  p=ggplot(mean_abun_, aes(x=x, y=mean_abun)) +
    geom_point()+
    geom_errorbar(aes(ymin=mean_abun-sd_abun, ymax=mean_abun+sd_abun), width=.1) +
    geom_line() +
    theme_classic()+
    ylab("Species abundance")+
    xlab("Species rank in abundance")+
    scale_x_continuous(breaks = 1:max_sp_n)
  print(p)
  abun=as.data.frame(abun)
  abun$m_community_ID=1:nrow(abun)
  list(sp_n=sp_n,abun=abun,total_abun=total_abun,mean_total_abun=mean_total_abun,sd_total_abun=sd_total_abun,mean_abun=mean_abun,sd_abun=sd_abun,mean_sp_n=mean_sp_n,sd_sp_n=sd_sp_n)
}
