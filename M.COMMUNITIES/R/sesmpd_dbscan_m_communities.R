#' 计算由dbscan聚类所划分小群落的标准化平均配对性状距离(SES.MPD) / Calculate the Standardized Effect Size of Mean Pairwise Trait Distance (SES.MPD) for the small clusters partitioned by DBSCAN clustering
#' @description 计算由dbscan聚类所划分小群落的标准化平均配对性状距离(SES.MPD) / Calculate the Standardized Effect Size of Mean Pairwise Trait Distance (SES.MPD) for the small clusters partitioned by DBSCAN clustering
#' @usage sesmpd_dbscan_m_communities(entity,traits,minx,miny,maxx,maxy,eps,minPts,replicate,edge)
#'
#' @param entity 植物个体的分布坐标与物种 / Distribution coordinates and species of plant individuals
#' @param traits 物种的功能性状 / Functional traits of species
#' @param minx 核心区边界，最小x / Core area boundary, minimum x
#' @param miny 核心区边界，最小y / Core area boundary, minimum y
#' @param maxx 核心区边界，最大x / Core area boundary, maximum x
#' @param maxy 核心区边界，最大y / Core area boundary, maximum y
#' @param eps dbscan聚类的eps参数 / eps parameter of DBSCAN clustering
#' @param minPts dbscan聚类的minPts参数 / minPts parameter of DBSCAN clustering
#' @param replicate 零模型重复次数 / Number of null model repetitions 
#' @param edge 基于边界筛选小群落的方法，可选"inside"与"centroid" / Method for screening small clusters based on boundaries, with optional "inside" and "centroid"
#'
#' @return 标准化平均配对性状距离 / Standardized Mean Pairwise Trait Distance
#' @export
#'
#' @examples
#' entity=M.COMMUNITIES::entity
#' traits=M.COMMUNITIES::traits
#' sesmpd_dbscan_m_communities(entity,traits,minx=5,miny=5,maxx=45,maxy=45,eps=2,minPts=1,replicate=99,edge="inside")
sesmpd_dbscan_m_communities=function(entity,traits,minx,miny,maxx,maxy,eps,minPts,replicate,edge)
{
  library(FD)
  library(picante)
  #####################
  #####################
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
  ################################
  ################################
  abundance_dbscan_m_communities=function(entity,minx,miny,maxx,maxy,eps,minPts,edge)
  {
    library(FD)
    #####################

    #####################
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

    ################################
    dbscan_m_community_classification=function(entity,minx,miny,maxx,maxy,eps,minPts)

    {
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
      ####?ų?ȫ??λ?ڱ߽??ڵ?СȺ??
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
    class=dbscan_m_community_classification(entity,minx,miny,maxx,maxy,eps,minPts)
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

  abundance=abundance_dbscan_m_communities(entity,minx,miny,maxx,maxy,eps,minPts,edge)
  abundance=abundance[,-ncol(abundance)]
  species=levels(as.factor(entity$species))
  Tot_species_abun=colSums(abundance)
  Tot_species_abun=as.data.frame(Tot_species_abun)
  Tot_species_abun=cbind(Tot_species_abun,species)
  m_community_population_num=rowSums(abundance)
  replicate_m_community=rep(NA,replicate)
  traits_dist=gowdis(traits)
  traits_dist=as.matrix(traits_dist)
  colnames(traits_dist)=rownames(traits)
  rownames(traits_dist)=rownames(traits)
  mpd_obv=mpd(abundance, traits_dist,abundance.weighted=T)
  mpd_obv[is.na(mpd_obv)]=0

  for (i in 1:nrow(abundance))
  {
    replicate_i=numeric()
    for (j in 1:replicate )
    {
      random_abu=t(as.matrix(rep(0,length(species))))
      random_m_community= sample(Tot_species_abun$species,size=m_community_population_num[i],replace = T,prob =Tot_species_abun$Tot_species_abun )
      random_m_community=table(random_m_community)
      random_m_community= t(as.matrix(random_m_community))
      replicate_i[j]=mpd(random_m_community, traits_dist,abundance.weighted=T)
      replicate_i[is.na(replicate_i)]=0
    }
    replicate_m_community=rbind(replicate_m_community,replicate_i)
  }
  replicate_m_community=replicate_m_community[-1,]
  null_mpd_mean=apply(replicate_m_community,1,mean)
  null_mpd_std=apply(replicate_m_community,1,sd)
  ses_MPD=-1*(mpd_obv-null_mpd_mean)/null_mpd_std
  result=cbind(abundance,mpd_obv,null_mpd_mean,null_mpd_std,ses_MPD)
  result[is.nan(result)]=0
  result
}

