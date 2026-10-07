PhyloTax=function(tax,tree,reftax,bootstrap_treshold=0.75,check_consistency=T,unassigned_level='Phylum',useNCBI=T){
  
  require(tidytree)
  require(ape)
  require(treeio)
  if(useNCBI==T){
    require(taxize)
  }
  
  # in the tax table change '' to 'Unassigned'
  tax$Kingdom[tax$Kingdom=='']='Unassigned'
  tax$Phylum[tax$Phylum=='']='Unassigned'
  tax$Class[tax$Class=='']='Unassigned'
  tax$Order[tax$Order=='']='Unassigned'
  tax$Family[tax$Family=='']='Unassigned'
  tax$Genus[tax$Genus=='']='Unassigned'
  tax$Species[tax$Species=='']='Unassigned'
  
  # in the reference tax table change '' to 'Unclassified'
  reftax[reftax$Kingdom=='','Kingdom']='Unclassified'
  reftax[reftax$Phylum=='','Phylum']='Unclassified'
  reftax[reftax$Class=='','Class']='Unclassified'
  reftax[reftax$Order=='','Order']='Unclassified'
  reftax[reftax$Family=='','Family']='Unclassified'
  reftax[reftax$Genus=='','Genus']='Unclassified'
  reftax[reftax$Species=='','Species']='Unclassified'
  
  #in the reference tax table change 'Unclassified' to the first known name
  for (i in 1:nrow(reftax)){
    if('Unclassified'%in%as.vector(reftax[i,2:8])){
      unclassified=grep('Unclassified',reftax[i,2:8])+1
      Lowest_classified=reftax[i,(min(unclassified)-1)]
      reftax[i,(grep('Unclassified',reftax[i,2:8])+1)]=paste('Unclassified',Lowest_classified)
    }
  }
  
  #first, create a duplicate of the taxonomy file
  phylo_tax=tax
  #use "Unassigned" as standardized name
  tax[] <- lapply(phylo_tax, function(x) as.character(x))
  tax[is.na(tax)]<-"Unassigned"
  tax[tax=='unassigned']<-"Unassigned"
  tax[tax=='unknown']<-"Unassigned"
  tax[tax=='Unknown']<-"Unassigned"
  tax[tax=='']<-"Unassigned"
  tax[tax==' ']<-"Unassigned"
  tax[tax=='uncultured']<-"Unassigned"
  tax[tax=='Uncultured']<-"Unassigned"
  
  
  while(nrow(tax[grep('Unassigned',tax[,unassigned_level]),])>0){
    
    #find rows containing "Unassigned"
    Unassigned_ASVs=tax[grep('Unassigned',tax[,unassigned_level]),]
    
    print(paste(nrow(Unassigned_ASVs),'Unassigned ASVs left to search'))
    
    #find the name of the Unassigned ASVs to be used as query
    Unassigned=rownames(Unassigned_ASVs)[1]
    print(paste('Searching tree using',Unassigned,'as query'))
    #find the parent node containing this ASV
    parent <- tree$edge[tree$edge[,2] == which(tree$tip.label == Unassigned), 1]
    
    #Climb up the tree until a known relative is found
    LEVELBACK=0
    known_relative=F
    LEVEL=NULL
    while(known_relative==F){
      tree.sub=tree_subset(tree,node =parent,levels_back = LEVELBACK)
      Tax.relatives=tax[rownames(tax)%in%tree.sub$tip.label,]
      if(length(Tax.relatives$Phylum[!Tax.relatives$Phylum%in%'Unassigned'])>0){
        #Once known relatives have been found, find the level at which they all have the same name
        if(length(unique(Tax.relatives$Species)[!unique(Tax.relatives$Species)%in%'Unassigned'])==1){
          LEVEL=7}else if(length(unique(Tax.relatives$Genus)[!unique(Tax.relatives$Genus)%in%'Unassigned'])==1){
            LEVEL=6}else if(length(unique(Tax.relatives$Family)[!unique(Tax.relatives$Family)%in%'Unassigned'])==1){
              LEVEL=5}else if(length(unique(Tax.relatives$Order)[!unique(Tax.relatives$Order)%in%'Unassigned'])==1){
                LEVEL=4}else if(length(unique(Tax.relatives$Class)[!unique(Tax.relatives$Class)%in%'Unassigned'])==1){
                  LEVEL=3}else if(length(unique(Tax.relatives$Phylum)[!unique(Tax.relatives$Phylum)%in%'Unassigned'])==1){
                    LEVEL=2}else if(length(unique(Tax.relatives$Kingdom)[!unique(Tax.relatives$Kingdom)%in%'Unassigned'])==1){
                      LEVEL=1}
        known_relative=T
      }else{LEVELBACK=LEVELBACK+1}
    }
    
    #find the bootstrap of the node that contains all unknown and the known relatives
    bs_tibble <- data.frame(tibble(node=1:Nnode(tree.sub) + Ntip(tree.sub),bootstrap = ifelse(tree.sub$node.label < 0, "", tree.sub$node.label)))
    
    if(nrow(bs_tibble)==1&as.character(bs_tibble[1,2])==''){bs_tibble[1,2]=1}
    
    #check that the boostratp value supports the clustering
    BOOT=bs_tibble[bs_tibble$node==min(as.numeric(bs_tibble$node)),'bootstrap']
    
    if(BOOT=='root'|BOOT==''|is.null(LEVEL)){BOOT=0}else{BOOT=as.numeric(BOOT)}
    
    if(BOOT>bootstrap_treshold){
      print(paste('known relatives of ASV',Unassigned,'found at',paste(colnames(Tax.relatives)[LEVEL]),'level'))
      relatives.unknown=rownames(Tax.relatives[Tax.relatives$Phylum=='Unassigned',])
      
      #check taxonomic consistency at each phylogeny level (l)
      Tax.relatives.known=Tax.relatives[!Tax.relatives[,LEVEL]%in%'Unassigned',]
      if(check_consistency){
        for (l in rev(1:(LEVEL-1)) ){
          inconsistencies=ifelse(length(unique(Tax.relatives.known[,l]))>1,T,F)
          
          if(inconsistencies){
            print(paste('inconsistent taxonomy detected at',colnames(Tax.relatives.known)[l],'level'))
            #if inconsistencies are found, replace the name at level l by the name in the reference database by searching the level l-1
            ## first, check that the taxonomy can be found in the reference database
            if(identical(unique(reftax[reftax[,colnames(Tax.relatives.known)[l+1]]%in%unique(Tax.relatives.known[,l+1]),l+1]),character(0))){stop('not found in ref tax')}
            ## check that the taxonomy is consistent at level l-1 in the reference taxonomy
            if(length(unique(reftax[reftax[,colnames(Tax.relatives.known)[l+1]]%in%unique(Tax.relatives.known[,l+1]),l+1]))==1){
              ### if it is consistent, use it to re-assign the taxonomy to the known relatives
              print(paste('Reference taxonomy used to re-assign taxonomy at',colnames(Tax.relatives.known)[l],'level'))
              Tax.relatives.known[,l]=unique(reftax[reftax[,colnames(Tax.relatives.known)[l+1]]%in%unique(Tax.relatives.known[,l+1]),l+1])
            }else{
              ### if it is not consistent, print a warnoing message
              print(paste('Inconsistencies in reference database for',colnames(Tax.relatives.known)[l+1],unique(Tax.relatives.known[,l+1])))
              if(useNCBI==T){
                #### if useNCBI is enabled, search the taxonomy against the NCBI reference database 
                print(paste('Searching',colnames(Tax.relatives.known)[l+1],unique(Tax.relatives.known[,l+1]),'against NCBI taxonomy'))
                search=classification(unique(Tax.relatives.known[,l+1]), db = 'ncbi')
                result=as.data.frame(search[1])
                ##### if a result is found, use the NCBI taxonomy to re-assign the taxonomy to the known relatives
                if(nrow(result)>1){
                  if('domain'%in%result[,2]){DO=result[result[,2]=='domain',1]}else{DO=NA}
                  if('phylum'%in%result[,2]){PH=result[result[,2]=='phylum',1]}else{PH=NA}
                  if('class'%in%result[,2]){CL=result[result[,2]=='class',1]}else{CL=NA}
                  if('order'%in%result[,2]){OR=result[result[,2]=='order',1]}else{OR=NA}
                  if('family'%in%result[,2]){FA=result[result[,2]=='family',1]}else{FA=NA}
                  if('genus'%in%result[,2]){GE=result[result[,2]=='genus',1]}else{GE=NA}
                  if('species'%in%result[,2]){SP=result[result[,2]=='genus',1]}else{SP=NA}
                  
                  result.df=data.frame(Species=SP,Genus=GE,Family=FA,Order=OR,Class=CL,Phylum=PH,Kingdom=DO)
                  Tax.relatives.known[,l]=result.df[,grep(colnames(Tax.relatives.known)[l+1],colnames(result.df))]
                }
              }else{
                #### if useNCBI is not enabled, stop the function and prompt to correct he reference taxonomy
                stop('Check consistency in reference database')
              }
            }
            
            Tax.relatives=Tax.relatives[!rownames(Tax.relatives)%in%rownames(Tax.relatives.known),];Tax.relatives=rbind(Tax.relatives,Tax.relatives.known)
          }
        }
      }else{
        Tax.relatives.known=Tax.relatives.known[1,,drop=F]
      }
      
      #Finally, re-assign the taxonomy of the Unassigned AVSs (in the phylo_tax tax file)
      print(paste('re-assigning taxonomy to',length(relatives.unknown)-1,'other ASVs'))
      for(j in 1:7){
        if(j<LEVEL){
          phylo_tax[relatives.unknown,j]=unique(Tax.relatives.known[,j])
        }else if(j==LEVEL){
          phylo_tax[relatives.unknown,j]=paste(unique(Tax.relatives.known[,j]),'like',sep='_')
        }else{
          phylo_tax[relatives.unknown,j]=paste('Unassigned',unique(Tax.relatives.known[,LEVEL]),'like',sep='_')
        }
      }
      #remove the ASVs that were previously unassigned from the original tax file so we don't look for them again
      tax=tax[!rownames(tax)%in%relatives.unknown,]
    }else{
      print(paste('Bootstrap value supporting AVS',Unassigned,'phylogeny is lower than',bootstrap_treshold,'; taxonomy remains Unassigned'))
      #remove the ASV that remains unassigned from the original tax file  so we don't look for it again
      tax=tax[!rownames(tax)%in%Unassigned,]
    }
  }
  
  if (nrow(tax[grep('Unassigned',tax[,unassigned_level]),])==0){write.csv(phylo_tax,paste('phylo_tax_bootstrapTreshold_',bootstrap_treshold,'.csv',sep=''),row.names = T)}
  
}
