#' 以基于物种存在的SES.MPD及SES.E对基于物种丰度的SES.MPD进行分解
#' @description 以基于物种存在的SES.MPD及SES.E对基于物种丰度的SES.MPD进行分解
#' @usage VPA_priority_enviroment_hc(entity,traits,minx,miny,maxx,maxy,scale,replicate,edge="centroid")
#' @param entity 植物个体的分布坐标与物种
#' @param traits 物种的功能性状
#' @param minx 核心区边界，最小x
#' @param miny 核心区边界，最小y
#' @param maxx 核心区边界，最大x
#' @param maxy 核心区边界，最大y
#' @param scale hc聚类的尺度
#' @param replicate 零模型重复次数
#' @param edge 基于边界筛选小群落的方法，可选"inside"与"centroid"
#'
#' @return 基于物种丰度(存在)的SES.MPD、SES.E
#'
#' @examples entity=M.COMMUNITIES::entity
#' @examples traits=M.COMMUNITIES::traits
#' @examples VPA_priority_enviroment_hc(entity,traits,minx=5,miny=5,maxx=45,maxy=45,scale=5,replicate=99,edge="centroid")
VPA_priority_enviroment_hc=function(entity,traits,minx,miny,maxx,maxy,scale,replicate,edge)
{
library(glmm.hp)
library(eulerr)
############################
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
################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################
################
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
################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################
################
sesmpd_hc_m_communities=function(entity,traits,minx,miny,maxx,maxy,scale,replicate,edge)
{
  library(FD)
  library(picante)
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

  abundance=abundance_hc_m_communities(entity,minx,miny,maxx,maxy,scale,edge)
  abundance=abundance[,-ncol(abundance)]
  species=levels(as.factor(entity$species))
  Tot_species_abun=colSums(abundance)
  Tot_species_abun=as.data.frame(Tot_species_abun)
  Tot_species_abun=cbind(Tot_species_abun,species)
  m_community_population_num=rowSums(abundance)
  replicate_m_community=rep(NA,replicate)
  traits_dist=gowdis(traits)
  traits_dist=as.matrix(traits_dist)

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
################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################
################
sdiversity_hc_m_communities=function(entity,minx,miny,maxx,maxy,scale,edge)
{
  library(vegan)
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

  abundance=abundance_hc_m_communities(entity,minx,miny,maxx,maxy,scale,edge)
  abundance=abundance[,-ncol(abundance)]
  Shannon=diversity(abundance)
  Simpson=diversity(abundance, "simpson")
  Pielou_E=Shannon/log(ncol(abundance))
  Max_relative_richness=apply(abundance,1,max)/rowSums(abundance)
  diversity=cbind(abundance,Shannon,Simpson,Pielou_E,Max_relative_richness)
  diversity
}
#########################////////////////##############\\\\\\\\\\\\\\\\######
#########################
#########################
rich_mpd=sesmpd_hc_m_communities(entity,traits,minx,miny,maxx,maxy,scale,replicate,edge)
ses_rich_mpd=rich_mpd[,"ses_MPD"]
#########################////////////////##############\\\\\\\\\\\\\\\\######
#########################
#########################
abundance=abundance_hc_m_communities(entity,minx,miny,maxx,maxy,scale,edge)
abundance=abundance[,-ncol(abundance)]
Pielou_E=sdiversity_hc_m_communities(entity,minx,miny,maxx,maxy,scale,edge)
Pielou_E_obv=Pielou_E[,"Pielou_E"]
species=levels(as.factor(entity$species))
Tot_species_abun=colSums(abundance)
Tot_species_abun=as.data.frame(Tot_species_abun)
Tot_species_abun=cbind(Tot_species_abun,species)
m_community_population_num=rowSums(abundance)
replicate_m_community=rep(NA,replicate)
for (i in 1:nrow(abundance))
{
  replicate_i=numeric()
  for (j in 1:replicate )
  {
    random_abu=t(as.matrix(rep(0,length(species))))
    random_m_community= sample(Tot_species_abun$species,size=m_community_population_num[i],replace = T,prob =Tot_species_abun$Tot_species_abun )
    random_m_community=table(random_m_community)
    random_m_community= t(as.matrix(random_m_community))
    replicate_i[j]=diversity(random_m_community)/log(length(species))
    replicate_i[is.na(replicate_i)]=0
  }
  replicate_m_community=rbind(replicate_m_community,replicate_i)
}
replicate_m_community=replicate_m_community[-1,]
null_Pielou_E_mean=apply(replicate_m_community,1,mean)
null_Pielou_E_std=apply(replicate_m_community,1,sd)
cbind(Pielou_E_obv,null_Pielou_E_mean,null_Pielou_E_std)
ses_Pielou_E=-1*(Pielou_E_obv-null_Pielou_E_mean)/null_Pielou_E_std
ses_Pielou_E[is.nan(ses_Pielou_E)]=0

#########################////////////////##############\\\\\\\\\\\\\\\\######
#########################
#########################
abundance=abundance_hc_m_communities(entity,minx,miny,maxx,maxy,scale,edge)
abundance=abundance[,-ncol(abundance)]
species=levels(as.factor(entity$species))
Tot_species_abun=colSums(abundance)
Tot_species_abun=as.data.frame(Tot_species_abun)
Tot_species_abun=cbind(Tot_species_abun,species)
m_community_population_num=rowSums(abundance)
replicate_m_community=rep(NA,replicate)
traits_dist=gowdis(traits)
traits_dist=as.matrix(traits_dist)
abundance1=abundance
abundance1[which(abundance1>1)]=1
occur_mpd_obv=mpd(abundance1, traits_dist,abundance.weighted=F)
occur_mpd_obv[is.na(occur_mpd_obv)]=0

for (i in 1:nrow(abundance))
{
  replicate_i=numeric()
  for (j in 1:replicate )
  {
    random_abu=t(as.matrix(rep(0,length(species))))
    random_m_community= sample(Tot_species_abun$species,size=m_community_population_num[i],replace = T,prob =Tot_species_abun$Tot_species_abun )
    random_m_community=table(random_m_community)
    random_m_community= t(as.matrix(random_m_community))
    random_m_community[which(random_m_community>1)]=1
    replicate_i[j]=mpd(random_m_community, traits_dist,abundance.weighted=F)
    replicate_i[is.na(replicate_i)]=0
  }
  replicate_m_community=rbind(replicate_m_community,replicate_i)
}
replicate_m_community=replicate_m_community[-1,]
null_mpd_mean=apply(replicate_m_community,1,mean)
null_mpd_std=apply(replicate_m_community,1,sd)
ses_occur_mpd=-1*(occur_mpd_obv-null_mpd_mean)/null_mpd_std
ses_occur_mpd[is.nan(ses_occur_mpd)]=0

data=cbind(ses_rich_mpd,ses_occur_mpd,ses_Pielou_E)
data=as.data.frame(data)
Data<<-subset(data,data[,3]>0)
data=cbind(abundance,data)
data
}
