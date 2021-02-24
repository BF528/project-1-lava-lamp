#!/usr/bin/env/Rscript

#Working path
#---------------
setwd('/projectnb/bf528/users/lava_lamp/project_1/esaake_pr1')

#packages
#----------
#installing R-package
#install.packages("ggplot2")
#install.packages("ggpubr")
# install gplots package
#install.packages("gplots")
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)
library(ggarrange)
library(ggpubr)
suppressPackageStartupMessages(library(dplyr))
library("gplots")

#csv file path
#-----------------
respath<-"esaake_pr1/newresultsf/"
lavalampdatapath<- "/projectnb2/bf528/users/lava_lamp/project_1/expression_data.csv"
mydata<-read.csv(lavalampdatapath) #reading csv data into mydata variable


#setting proper rowname
#-----------------------------
#print (mydata[,1])#first column
rownames(mydata)=mydata[,1];
mydata=mydata[,-1]  #dropping the duplicate column used to list_index the rows
gexp_data<-mydata; #Just a copy

#Dimension Reduction Function Definition
#---------------------------------------

#Test 1: Log Fold change Test
  #Function Execute_Logfoldtest()
    #Input: Dataframe or matrix
    #Output: List of row indices 


execute_logfoldtest <- function(givendata){

  #variable initialization
    numofrows<-nrow(givendata)
    numofcolumns<-ncol(givendata)
    logthreshold<-log2(15)
    topass_percent<- (numofcolumns)*0.2
    test1passers<-c();
    count=0;
    list_index=1;
    
  #row traversing loop
    for (i in 1:numofrows)
    {
        count=0 #counter reset after each row
        
        #column traversing loop
        for (j in 1:numofcolumns){
            if(givendata[i,j]>logthreshold){
              count=count+1
            }
        }
        if(count>(topass_percent)){
            test1passers[list_index]<-i;
            list_index<-list_index+1;
        }
    }
  return(test1passers);
}


#Test 2: Chisquare test for variance
#------------------------------------
  #Function execute_chisquarevariancetest()
    #Input: Dataframe or matrix
    #Output: List of row-indices

execute_chisquarevariancetest<- function(givendata) {
  
  #Variable Initialization
    numofrows<-nrow(givendata)
    numofcolumns<-ncol(givendata)
    count<-0;
    list_index<-1;
    degree_of_freedom<-numofcolumns-1;
    alpha<-0.01;
    test2passers<-c(); #empty list to hold list of passers
    qchisq_values_lower<-qchisq(alpha/2, degree_of_freedom);
    qchisq_values_upper<-qchisq((1-alpha)/2, degree_of_freedom);
  
  #from https://www.itl.nist.gov/div898/handbook/eda/section3/eda358.htm:
    #T=degree of freedom * square(sample_stddev/targetstddev)^2
  
  #sample standard deviation
    sample_stddev<-sd(as.matrix(givendata))
    
  for (i in 1:numofrows){
    
    test_statistic_per_gene<-(numofcolumns-1)*((sample_stddev/sd(givendata[i,])**2))
    
    #if test_stat is not within critical region
    
    if ((test_statistic_per_gene< qchisq_values_lower)||(test_statistic_per_gene>qchisq_values_upper)){
      test2passers[list_index]<-i;
      list_index=list_index+1 #increment of indexing variable
    }
  }
  return(test2passers)
}


#Test 2: Functin to test Coefficient of variation > 0.186
#------------------------------------
  #Function execute_CoVtest
    #Input: Dataframe or matrix
    #Output: List of row-indices


execute_CoVtest<-function(givendata){
  #Variable Initialization
    numofrows<-nrow(givendata)
    numofcolumns<-ncol(givendata)
    test3passers<-c();
    count<-0;
    list_index<-1;
    CoVthreshold<-0.186
    
  #loop to traverse rows
  for (i in 1:numofrows){
    gexpmean<-mean(as.numeric(givendata[i,]))
    std_dev<-sd(givendata[i,])
    CoV<-std_dev/gexpmean
    #test for condition of cov>CoVthreshold
    if (CoV>CoVthreshold){
      test3passers[list_index]<-i
      list_index=list_index+1
    }
  }
  return(test3passers)
}

#Function calls:
#----------------

  #Test1 execution
  logfoldtestpassers <-execute_logfoldtest(gexp_data);
  
  #Test2 execution
  chisquaretestpassers<-execute_chisquarevariancetest(gexp_data);

  #Test3 execution
  CoVtestpassers<-execute_CoVtest(gexp_data)

#Final Results and summary
#--------------------------
  #Combining passers of all three tests
        #Test1 execution
          #logfoldtestpassers #chisquaretestpassers #CoVtestpassers
    
          #using intersect we find intersection of all 3 sets of values
          #intersect(a,b) take two args #intersect all three matrices: intersect() fails with dataframe: conv to vector [[]]/[,1]
    ind_passed_all_threetests<-intersect(intersect(logfoldtestpassers, CoVtestpassers),chisquaretestpassers)
          
    #converting the list to dataframe so I can use indexes to truncate data out of original
    ind_passed_all_threetests<-data.frame(ind_passed_all_threetests);
    
    #Filtering out gnes that satisfied all three set i.e. the row from the original data
    dimreducedgene_expdata<- gexp_data[ind_passed_all_threetests[[1]], ]

    #data from 4.2 Only test2 passers
    dimreducedata4_2<- gexp_data[chisquaretestpassers, ]
    
#Deliverables for part4:
#----------------------------
    #write instruction
    #write.csv(dimreducedgene_expdata,"../gene_expression_sec4.csv")
    
    
    #Number that passed all three tests
    #Q4.2
    passalltests=count(ind_passed_all_threetests);
    print (paste("Genes that passed all three tests =",passalltests , sep=" ")) #1323
    
    #summary
    #logfoldtestpassers #chisquaretestpassers #CoVtestpassers
        countorig<-as.numeric(nrow(mydata)) #54675
        countpasslog<-as.numeric(length(logfoldtestpassers)) #39661
        countpasschisq<-as.numeric(length(chisquaretestpassers)) #54416
        countpassCoV<-as.numeric(length(CoVtestpassers)) #1699
        passall3<-as.numeric(passalltests) #1323
        datatab=c(countorig, countpasslog, countpasschisq, countpassCoV, passall3) #creating datafram for writing to csv
    
    theheadernames<-c("Num_of_total_genes_from_preprocessed", "Num_of_genes_passed_log_test", "Num_of_genes_outside_critical_region", "num_of_genes_passed_CoV_greater_0.186", "Num_of_genes_passed_all_three_tests")
    summary_result<-data.frame(count=unlist(datatab), row.names = unlist(theheadernames))
    ##write.csv(summary_result, "newresultsf/summarytable_p4.csv")

#==================================  
#ANALYST: PART 5
#==================================
    
    #Setting seed for reproducibility
      set.seed(35)
    #Variable Initialization
    #Transposing data to allow for clustering users
    transposeddata<-t(dimreducedgene_expdata);
    p5_data<-data.frame(transposeddata);
    
   #testing for null values
    any(is.na(p5_data))#To check for missing or null values
    
  #Scaling values
    data_scaled<-as.data.frame(scale(p5_data))
    
    #Finding distance between data points by using Euclidean
    distanced_data<-dist(data_scaled, method='euclidean')
    
  #computing hierarchical cluster average
    cluster_average<-hclust(distanced_data, method='average')
    plot(cluster_average)
    
    
 #Metadata #needed to get the cit.coloncancermolecularsubtype
      metadata<-read.csv("/project/bf528/project_1/doc/proj_metadata.csv");
      metadata<-data.frame(metadata);
      ##write.csv(metadata, "newresultsf/metacopy.csv")
    
    
    #nleaves(cluster_average)
 #Cutting dendo into two clusters
    new_trim<-cutree(cluster_average, k=2)
    #storing the cluster return as dataframe
    new_trim_df<-data.frame(new_trim)
    
    
#Lets redraw the cluster and include two cut lines
    plot(cluster_average)
    rect.hclust(cluster_average, k=2, border=2:4)
    abline(h=NULL, col='red');#
    
#Draw the dendogram in different colors
    suppressPackageStartupMessages(library(dendextend)); #required package
    dendo_obj_avg<-as.dendrogram(cluster_average) #converting cluster average to dendogram object
    colored_dendo_of_dendo_obj_avg<-color_branches(dendo_obj_avg, k=2) #using two diff colors to draw the dendogram
    
#opening an image object for saving
    ##jpeg('newresultsf/clusterplot11.jpg', width = 600, height = 600)
    plot(colored_dendo_of_dendo_obj_avg)
    ##dev.off()
    
    
#Create a new dataframe storing results
    p5_data_cluster<-mutate(as.data.frame(p5_data), cluster=new_trim)
    table(p5_data_cluster$cluster)
    labelz<-metadata[ ,'geo_accession'];#get sample shortname frm metadata 
    clustertable<-table(p5_data_cluster$cluster,labelz)
    #write to csv
    ##write.csv(clustertable, "newresultsf/clustertable.csv")
    
    
#Plotting results
    suppressPackageStartupMessages(library(ggplot2))
    
    
    #must be cluster #solution to Error: 'x' must be a numeric matrix #heatmap(as.matrix(dataset[, -1]))
    p5_data_cluster[, -1]
    
#setting up color scheme
    colsidecolfunc<-function(somedata) if(somedata['cit.coloncancermolecularsubtype']=='C3')"red" else "blue"
    colside_colors<-apply(metadata,1,colsidecolfunc) 

#setting color for colsidecolors #ColSideColors only takes vectors
    rowheaders<-row.names(p5_data);
    columnames<-names(p5_data);
    
#saving heatmap to jpeg
      ##jpeg('newresultsf/heatmapplotlabel.jpg',  width = 900, height = 600 )
      #heatmap(as.matrix(t(p5_data_cluster[,-1])), ColSideColors = colside_colors, label=columnames)
    heatmap.2(as.matrix(t(p5_data_cluster[,-1])), ColSideColors=colside_colors , main="Gene Expression Heatmap", ylab="Genes",  xlab="Samples")
      ##dev.off()
    
    
 #Welch Test
#=============
      new_trim<-data.frame(cluster=new_trim); # Renaming Data Frame columc cluster 
    
    #welchdata<-truncated_data_original_form;

    #setting function for p.values 
      mytest_pval<- function(welchdata) t.test(welchdata[new_trim$cluster==1], welchdata[new_trim$cluster==2], var.equal=T)$p.value
    
    #setting function for t statistics 
     mytest_stat<- function(welchdata) t.test(welchdata[new_trim$cluster==1], welchdata[new_trim$cluster==2], var.equal=T)$statistic
        #using apply to apply the t.test on data from part 4 #apply(X, M, FUN) #-M1`: the manipulation is on rows #-M2`: the manipulation is on columns
           pvals<- apply(dimreducedgene_expdata, 1, mytest_pval)
          t_statz<-apply(dimreducedgene_expdata, 1, mytest_stat)

    
#p-values less than 0.05
#---------------------------
      sum(pvals< 0.05)
      hist(pvals);#draw histogram with pvals
    
  #computing p adjusted values
  #------------------------------
      p_adj<-p.adjust(pvals, method="fdr")
      p_adj_less005_count<-sum( p.adjust(pvals, method="fdr") < 0.05 );# counting those that are less than 0.05
      p_adj_less005_data<-p.adjust(pvals, method="fdr") < 0.05 ;#for getting the list of best rep later

    
  #Storing all results in a dataframe
      sec5_4_df<-data.frame(t_statistics=t_statz , p_values=pvals, adjusted_p=p_adj)
      sec5_4_df<- sec5_4_df %>% arrange(desc(sec5_4_df$t_statistics))
      #write.csv(sec5_4_df,"../differential_expression_s5.csv")

      
    
  #Sec 5 Answers
  #=================
    #computing number of samples in each cluster
    #q1
    print(paste("Number of samples in cluster 1 =", count(p5_data_cluster[p5_data_cluster$cluster==1,]), sep=" ")) #59
    print(paste("Number of samples in cluster 2 =", count(p5_data_cluster[p5_data_cluster$cluster==2,]), sep=" ")) #75

    #q3
    print (paste("The number of differentially expressed genes at p<0.05 =", p_adj_less005_count, sep=" ")) #1064
    
    #q4
   
    p_adj_less005_data<-p_adj_less005_data[p_adj_less005_data==TRUE]
    p_val_less005<-pvals<0.005;
    p_val_less005<-p_val_less006[p_val_less005==TRUE];
    
    best_define_clust<-p_adj_less005_data;
    print("genes that best represent clusters are:")
    print(as.data.frame(best_define_clust))
    write.csv(best_define_clust, "newresultsf/p_adjust_less005_bestclust.csv")
    write.csv(p_val_less005, "newresultsf/p_val_less005.csv")
    
    #q5
    cat("This Genes have their p_adjust values less than 0.05. They are below their own threshhold[(0.05/m) * 14]")
    
    
    #summary
    distance_method<-'euclidean' ; cluster_method<-'average'; num_of_clusters<-2; cluster_type<-'hierarchical'
    pval_method<-'fdr';clust1<-as.numeric(ngroup1); clust2<-as.numeric(ngroup2); 
    num_of_c3_colon_cancer<-sum(subtypecolors=="red")
    num_of_other_subtype_of_coloncancer<-sum(subtypecolors=="blue")
    subtypecolors<-data.frame(subtype=subtypecolors, row.names= NULL)
    data_par=c(distance_method, cluster_method, cluster_type, pval_method, num_of_clusters, 
               clust1, clust2, num_of_c3_colon_cancer, num_of_other_subtype_of_coloncancer)
    parnames<-c('distance_method', 'cluster_method', 'cluster_type', 'pval_method',
                'number_of_clusters', 'number_of_samples_in_Cluster1',
                'number_of_samples_in_Cluster2', 'num_of_c3_colon_cancer', 'num_of_not-C3_coloncancer')
    num_of_genes_with_pval_less_0.05<-as.numeric(p_adj_less005_count)
    
    
    
    summary_result_p5<-data.frame(count=data_par, row.names = unlist(parnames))
    ##write.csv(summary_result_p5, "newresultsf/summarytable_p5.csv")
    
  #section 5.6
    dimreducedata4_2
    #setting function for p.values 
    mytest42_pval<- function(wdata) t.test(wdata[new_trim$cluster==1], wdata[new_trim$cluster==2], var.equal=T)$p.value
    
    #setting function for t statistics 
    mytest42_stat<- function(wdata) t.test(wdata[new_trim$cluster==1], wdata[new_trim$cluster==2], var.equal=T)$statistic
    
    #using apply to apply the t.test on data from part 4 #apply(X, M, FUN) #-M1`: the manipulation is on rows #-M2`: the manipulation is on columns
    pvals42<- apply(dimreducedata4_2, 1,  mytest42_pval)
    t_statz42<-apply(dimreducedata4_2, 1,  mytest42_stat)
    
    #p_adjust
    p_adj42<-p.adjust(pvals42, method="fdr")
    sec5_6_df<-data.frame(t_statistics=t_statz42 , p_values=pvals42, adjusted_p=p_adj42)
    sec5_6_df<- sec5_6_df %>% arrange(desc(sec5_6_df$t_statistics))
    #write.csv(sec5_6_df,"../diffexp_s5_6.csv")
    
#gene generation 
#---------------------
  library(tidyverse)
  library(tidyr)

  # get affymetrix probe id names
  #probe.ids <- row.names(eData)
  #----------------------------
    eData<-dimreducedgene_expdata
    dat<-data.frame(p_adj_less005_data)
    probe.ids <- rownames( dat)
    
  # IDs to Chromosome location, options: columns(hgu133plus2.db)
    genedata <- tibble(probes=probe.ids,
                   symbols=AnnotationDbi::mapIds(hgu133plus2.db,keys=probe.ids,
                                                  column='SYMBOL',keytype='PROBEID',multiVals='first'),
                   map=AnnotationDbi::mapIds(hgu133plus2.db,keys=probe.ids,
                                                 column='MAP',keytype='PROBEID',multiVals='first'))

     
     #Deletion of  NA  
     genedata <- genedata %>% filter(symbols!='') %>% dplyr::select(-map)
      
     # Get Chromosome location annotation data (a little more work)
     xchroloc <- hgu133plus2CHRLOC
     
     # Get the probe identifiers that are mapped to chromosome locations
     mapped.probes <- mappedkeys(xchroloc)
     
     # Convert to a list
     xxchrloclist <- as.list(xchroloc[mapped.probes])
     
     match.probes <- which(names(xxchrloclist) %in% genedata$probes)
     
     genedata <- genedata %>% filter(probes %in% names(xxchrloclist[match.probes]))
     
     genedata$loc <- lapply(xxchrloclist[match.probes],'[[',1) %>% as.numeric %>% abs
      
      # Drop any duplicated values
      genedata <- genedata[-which(duplicated(genedata[,-1])),]
      
      # Lineing up expression data
      eData <- eData[which(probe.ids %in% genedata$probes),]
      
      ##write.csv(genedata,"newresultsf/finalgeneoutput.csv")
    

      