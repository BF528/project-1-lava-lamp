#!/usr/bin/env Rscript


#setting working directory in R
setwd('/projectnb/bf528/users/lava_lamp/project_1')


#installing R-package
install.packages("ggplot2")
#install.packages("ggpubr")
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
#library(ggplot2)
#library(ggarrange)
library(ggpubr)
suppressPackageStartupMessages(library(ggplot2));
suppressPackageStartupMessages(library(dplyr))

#Microarray analysis "large p, small n" problem
# number of features (p, e.g. genes) exceeds the number of samples (n, e.g. patients).
# clustering, will yield little or no information when performed in scenarios where p >> n due to the low signal-to-noise ratio. 
#setting path to variables:
samplepath<-"/project/bf528/project_1/data/example_intensity_data.csv"
lavalampdatapath<- "/projectnb2/bf528/users/lava_lamp/project_1/expression_data.csv"

#reading csv data into mydata variable
mydata<-read.csv(lavalampdatapath)
mydata_rowdone<-mydata;

#Setting row names
#print (mydata[,1])#first column
rownames(mydata_rowdone)=mydata_rowdone[,1];

#This deletes the first column that was duplicated above as an index column
mydata_rowdone=mydata_rowdone[,-1] #for statistical computation


#Transpose data so I can use column functions like sd(), var(), mean(); 
transposed_data=t(mydata_rowdone); #NB: t() returns matrix 


#co-efficient of variation = (standard deviation)/Mean

genes_that_pass_log15_test<-function(somedata){
  
  #calculating number of columns
  number_of_column = ncol(somedata)
  
  #Getting sample size; number of rows for each gene
  N=nrow(somedata);
  
  #Computing our degree of freedom
  df= N-1;
  
  gene_expr_lowerbound= log2(15);
  numb_percentage= (N)*0.2;
  count_greater_than_log=0;
  
  #lists to store the column index of genes that specify conditions
  #Empty list to store the index numbers of genes that pass the log(15) test
  columns_that_pass_log15<-c();
  
  #index to store the column numbers of genes that pass the log2(15) test
  control_v=1;#R index starts from 1
  
  for (i in 1:number_of_column){
    #reset for each gene
    count_greater_than_log=0;
    
    for (x in 1:N){
      if (somedata[x,i]> gene_expr_lowerbound){
        count_greater_than_log= count_greater_than_log + 1;
      }
      
    }
    
    if (count_greater_than_log > numb_percentage){
      columns_that_pass_log15[control_v]<-i;
      #increasing the index to the next;
      control_v=control_v+1;
      
    }
    
  }#end of second for loop
  
  return (columns_that_pass_log15);
} #end of first for loop



##Calc Cofficient of var/ each gene

compute_coefficient_of_variation_each_gene<-function(somedata){
  co_of_var = c();
  
  #for iter
  Nc<-ncol(somedata);
  
  
  #loop to run through each gene/column to computer coefficient of variation
  for (i in 1:Nc){
    cv <- sd(somedata[,i]) / mean(somedata[,i])
    co_of_var[i]<-cv;
    
  }
  
  #returning the list or df of coefficient of variatin for each gene
  return(co_of_var);
}


##TEST STATISTICS

compute_T<-function(somedata){
  #whatever somedata we passed is the transposed version of the original dataset
  #to allow us use column functions
  #hence rows are columns and columns are rows
  
  #Number of genes (each orginal row -> but now rows are columns in the transposed data)
  Lr= ncol(somedata);
  
  #initalizing empty list to store T values
  T_values=c();
  
  #Number of samples (i.e. each column)
  #but now rows are columns so to find number of samples or columns we use nrows()
  N=nrow(somedata)
  Ndf=N-1;
  
  for(i in 1:Lr){
    t_val<-(Ndf)*((sd(somedata[,i])/var(somedata[,i]))**2) 
    T_values[i]<-t_val
  }
  
  return (T_values);
} #END OF FUNCTION



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##################FUNCTION TO COMPUT QCHISQ()#############
#Takes transposed data

#Threshold P < 0.01
#Test statistics T= (N-1)*(s/@)^2
#Observed Variance
#Standard deviation
#Degree of freedom = N-1
#N- sample size
#key aspect is ratio of sample standard deviation to the
#...target standard deviation (s/@)
# The more the ratio above deviates from 1, the more reject-likely 
#...the null hypo
#significance level - alpha


compute_QCHISQ<-function(somedata){
  #data is transposed: use column functions
  
  #Number of genes (each orginal row -> but now rows are columns in the transposed data)
  Lr= ncol(somedata);
  
  #initalizing empty list to hold qchisq values both upper and lower
  #list to store list of upper and lower T qchisq values and both in qchisq_values
  qchisq_values_upper=c();
  qchisq_values_lower=c();
  #Number of samples (i.e. each column)
  #but now rows are columns so to find number of samples or columns we use nrows()
  N=nrow(somedata)
  degree_of_freedom=N-1;
  
  #level of significance
  alpha=0.01;
  
  for(i in 1:Lr){
    qchisq_lower<-qchisq(alpha/2, degree_of_freedom);
    qchisq_upper<-qchisq((1-alpha)/2, degree_of_freedom);
    qchisq_values_lower[i]<-c(qchisq_lower);
    qchisq_values_upper[i]<-c(qchisq_upper);
  }
  qchisq_value=data.frame(qchisq_values_lower=qchisq_values_lower,qchisq_values_upper=qchisq_values_upper)
  return (qchisq_value);
}#End of TCHISQ


genes_outside_the_critical_region<-function(qchisq_est, T_Statistics){
  nrows<-nrow(T_Statistics); #Both qchisq_est and T_Staistics have same length
  
  indices_genes_outside_critical_region<-c();
  control=1;# R index starts from 1
  for(i in 1:nrows){
    
    if((T_Statistics[i,1] < qchisq_est[i,1]) || (T_Statistics[i,1]> qchisq_est[i,2])){
      indices_genes_outside_critical_region[control]<-i;
      control=control+1;
    }
  }
  return(indices_genes_outside_critical_region);
} #end of function


#function to return the indices after testing the coefficient of var output
Indices_of_Genes_that_passed_CoV_threshold<-function(co_of_var_result){
  #initializing variables
  CoV_threshold<-0.186; #threshold for coefficient of variation
  thelength=nrow(co_of_var_result); #upper bound for looping
  indexer=1; #initializing indexer for my list; starts from 1
  indices_of_genes_that_passed<-c(); #An empty list
  
  for (i in 1:thelength){
    
    if(co_of_var_result[i,1] > CoV_threshold){
      indices_of_genes_that_passed[indexer]<-i
      indexer=indexer+1;
    }
    
  }
  return(indices_of_genes_that_passed);
}#end of function


##FUNCTIONS RUNS

#%%%%%%%TEST 1
#calling functions genes_that_pass_log15_test
#returns index of columns that pass the test
columns_that_pass_log15<-genes_that_pass_log15_test(transposed_data);
columns_that_pass_log15<-data.frame(columns_that_pass_log15);


#%%%%%%TEST 2

#calling function compute_TCHISQ for  each gene
qchisq_est <-compute_QCHISQ(transposed_data);

#converting tchisq to dataframe 
#qschiq_df<-data.frame(qschiq_vals);

#calling function compute_T to compute T_values for each gene
Test_Statistics <-compute_T(transposed_data);

#converting T_values to dataframe 
T_Statistics<-data.frame(Test_Statistics);


indices_of_genes_outside_CR<-genes_outside_the_critical_region(qchisq_est, T_Statistics);
indices_of_genes_outside_CR<-data.frame(indices_of_genes_outside_CR);

#%%%%%%TEST 3
#Calling function for coefficient of variation
#returns matrix of each genes co-efficient of variation
co_of_var_result<-compute_coefficient_of_variation_each_gene(transposed_data);
co_of_var_result<-data.frame(co_of_var_result);

#calling function to return index of genes that passed CoV test
Genes_that_passed_CoV<-Indices_of_Genes_that_passed_CoV_threshold(co_of_var_result);
Genes_that_passed_CoV<-as.data.frame(Genes_that_passed_CoV);

#COMBINING PASSERS OF ALL THREE TESTS
#All genes that passed all test are the intersection of all three list of passers
#columns_that_pass_log15
#Genes_that_passed_CoV
#indices_of_genes_outside_CR

#using intersect we find intersection of all 3 sets of values
#intersect(a,b) can only take two arguments at a time
#intersect -> !fails with dataframe:to vector using [[]] or [,1]
ind_passed_all_threetests<-intersect(intersect(columns_that_pass_log15[[1]],Genes_that_passed_CoV[[1]] ),indices_of_genes_outside_CR[[1]])

#converting the list to dataframe so I can use indexes to truncate data out of original
ind_passed_all_threetests<-data.frame(ind_passed_all_threetests);

#truncating original data set
truncated_data<- transposed_data[ , ind_passed_all_threetests[[1]]]

#transposing the result to its original form
truncated_data_original_form<-t(truncated_data)

#DELIVERABLES PART 4
#write instruction
write.csv(truncated_data_original_form,"/projectnb/bf528/users/lava_lamp/project_1/gene_expression_sec4.csv")

#Number that passed all three tests
#q4.2
passalltests=count(ind_passed_all_threetests);
print (paste("Genes that passed all three tests =",passalltests , sep=" ")) #1193



###PART 5
#@@@@@@@@@

#Detailing informatin about dataset
p5_data<-data.frame(truncated_data);

cat(names(p5_data))
#testing for null values
any(is.na(p5_data))#To check for missing or null values

#Scaling values
data_scaled<-as.data.frame(scale(p5_data))
#summary(data_scaled)


#Finding distance between data points by using Euclidean
distanced_data<-dist(data_scaled, method='euclidean')

#computing hierarchical cluster average
cluster_average<-hclust(distanced_data, method='average')
plot(cluster_average)

#nleaves(cluster_average)
#Cutting into two cluster of dendograms using cutree()
#In clustering with hclust, you cluster columns according to rows
new_trim<-cutree(cluster_average, k=2)
#storing the cluster return as dataframe
new_trim_df<-data.frame(new_trim)


#Lets redraw the cluster and include two cut lines
plot(cluster_average)
rect.hclust(cluster_average, k=2, border=2:4)
abline(h=NULL, col='red');#

#Draw the dendogram in different colors
#install.packages('package_name', dependencies =True);
suppressPackageStartupMessages(library(dendextend));

#setting seed for reproducibility
set.seed(1)

#converting cluster average to dendogram object
dendo_obj_avg<-as.dendrogram(cluster_average)

#using two different colors to draw the dendogram
colored_dendo_of_dendo_obj_avg<-color_branches(dendo_obj_avg, k=2)
plot(colored_dendo_of_dendo_obj_avg)

#saving cluster data
suppressPackageStartupMessages(library(dplyr))

#Appending the new cluster to the initial data and storing in p5_data_cluster
#Create a new dataframe storing results
p5_data_cluster<-mutate(as.data.frame(p5_data), cluster=new_trim)



#Plotting results
suppressPackageStartupMessages(library(ggplot2))

#setting seed for reproducing output
set.seed(110)

#must be cluster
#solution to Error: 'x' must be a numeric matrix
#heatmap(as.matrix(dataset[, -1]))
p5_data_cluster[, -1]


##function to set the subtype colors
setcolsidecolors<-function(typedata){
  colside_colors=c();
  
  #nrow gave errors due to data type of form vector
  u_limit=length(typedata);
  for (i in 1:u_limit)
  {
    if (typedata[i]=="C3"){
      colside_colors[i]<-"red"
    }else{
      colside_colors[i]<-"blue"
    }
  }
  
  return(colside_colors)
}


#Metadata 
metadata<-read.csv("/project/bf528/project_1/doc/proj_metadata.csv");
metadata<-data.frame(metadata);
cancer_mol_subtype<-metadata[ ,"cit.coloncancermolecularsubtype"]

#setting color for colsidecolors with function secolsidecolors
#ColSideColors only takes vectors
rowheaders<-row.names(p5_data);
columnames<-names(p5_data);
subtypecolors<-setcolsidecolors(cancer_mol_subtype)


heatmap(as.matrix(t(p5_data_cluster[,-1])), ColSideColors = subtypecolors, label=columnames)

#Welch Test

# Renaming Data Frame  
new_trim<-data.frame(cluster=new_trim);

welchdata<-truncated_data_original_form;

set.seed(150)

#t.test produces stnd output
#result var:(these are column names you can use)
#$p.val and $statistic names from t.test output

#setting function for p.values
mytest_pval<- function(welchdata) t.test(welchdata[new_trim$cluster==1], welchdata[new_trim$cluster==2], var.equal=T)$p.value

#setting function for t statistrics only
mytest_stat<- function(welchdata) t.test(welchdata[new_trim$cluster==1], welchdata[new_trim$cluster==2], var.equal=T)$statistic

#using apply to apply the t.test on a given data
#apply(X, M, FUN) #-M1`: the manipulation is on rows #-M2`: the manipulation is on columns
pvals<- apply(truncated_data_original_form, 1, mytest_pval)
t_statz<-apply(truncated_data_original_form, 1, mytest_stat)


#p-values less than 0.05
sum(pvals< 0.05)
hist(pvals);#draw histogram with pvals

#computing p adjusted values
p_adj<-p.adjust(pvals, method="fdr") 
p_adj_less005_count<-sum( p.adjust(pvals, method="fdr") < 0.05 );# counting those that are less than 0.05
p_adj_less005_data<-p.adjust(pvals, method="fdr") < 0.05 ;#for getting the list of best rep later


#Storing all results in a dataframe
sec5_4_df<-data.frame(t_statistics=t_statz , p_values=pvals, adjusted_p=p_adj)
sec5_4<- sec5_4_df %>% arrange(desc(sec5_4_df$t_statistics))
#sec5_4_df_ext<-data.frame(t_statistics=t_statz , p_values=pvals, adjusted_p=p_adj, less_than_0.05=p_adj_less_005)

write.csv(sec5_4_df,"/projectnb/bf528/users/lava_lamp/project_1/differential_expression_s5.csv")

#function to return name that have significant p-values
best_rep_cluster<-function(data){
  genes_best_rep_clusters<-c();
  t=1;
  for (i in 1:length(data)){
    if (data[i] == TRUE){
      genes_best_rep_clusters[t]<-names(data[i]);
      t=t+1;
    }
  }
  return (genes_best_rep_clusters)
}



#Answers
#computing number of samples in each cluster
ngroup1=count(p5_data_cluster[p5_data_cluster$cluster==1,])
ngroup2=count(p5_data_cluster[p5_data_cluster$cluster==2,])
print(count(p5_data_cluster))

#printing out results
#q1
print(paste("Number of samples in cluster 1 =", ngroup1, sep=" ")) #57
print(paste("Number of samples in cluster 2 =", ngroup2, sep=" ")) #77

#q2-Go above

#q3
print (paste("The number of differentially expressed genes at p<0.05 =", p_adj_less005_count, sep=" ")) #962

#q4
genes_best_define_clust <- best_rep_cluster(p_adj_less005_data);
cat("genes that best represent clusters are:")
print(as.data.frame(genes_best_define_clust))

#q5
cat("This Genes have their p_adjust values less than 0.05. They are below their own threshhold[(0.05/m) * 14]")
