

get_ds_transcript_vect<-function(in_col_vect1, target_depth){
  trans_vect<-get_sim_transcripts(in_col_vect1, target_depth)
  new_vect<-rep(0,length(in_col_vect1))
  for (i in trans_vect)  {
    new_vect[i]<-new_vect[i]+1
  }
  return(new_vect)
}


get_sim_transcripts<-function(in_col_vect, target_depth){
  nz_rows<-which(in_col_vect!=0)
  transcript_vect<-c()
  for (i in 1:length(nz_rows)){
    temp_id<-nz_rows[i]
    transcript_vect<-c(transcript_vect,rep(temp_id,in_col_vect[nz_rows[i]]))
  }
  return(sample(transcript_vect,target_depth))
}


#' Downsample columns of a matrix
#'
#' This function takes in a matrix of counts & will downsample 
#' each column such that their colSums are equivalent. Note that you can't just
#' randomly subtract 1 from cells until they're equivalent. To recapitulate the
#' actual biases from Poisson sampling, we actually create the "bag of marbles"
#' and sample each until we get the the target depth. You can either pass in 
#' the target depth, or if not, it will be automatically set to the minimum of 
#' the colSums, or if you are doing an integration & want to downsample 
#' different datasets, you should first FILTER OUT the columns whose colSums 
#' are less than the cutoff, then set this value to the cutoff.
#'
#' @param in_mat An input matrix
#' @param target_depth This can either be nothing or a specific target to 
#' which we will normalize to. Note that this should be the minimum of colsum 
#' across all datasets.
#' @return A matrix of the normalized matrix
#' @examples
#'    in_mat<-matrix(rnbinom(1000*20000,.1,c(0.5,0.5,0.5)),ncol=1000,nrow=20000)
#'    out_mat<-downsample_mat(in_mat)
#' @importFrom parallel detectCores makeCluster parApply stopCluster
#' @name downsample_mat
#' @export
downsample_mat<-function(in_mat, target_depth=NULL, rm_less=TRUE, parallel=TRUE, cores=NULL, quiet=FALSE){
  if (is.null(target_depth)){
    target_depth<-min(colSums(in_mat))
  }
  if (parallel){
    if (is.null(cores)){
        cores<- detectCores()
    }
    cl <- makeCluster(cores)
    clusterExport(cl, c("get_ds_transcript_vect", "get_sim_transcripts"),
              envir=environment())
    if (!quiet){    
      message("        running with ",cores," cores")
    }
    out_mat <- parApply(cl, X=in_mat, MARGIN=2, FUN=get_ds_transcript_vect, target_depth=target_depth)
    stopCluster(cl)
  } else {
    out_mat <- apply(in_mat, 2, function(x) get_ds_transcript_vect(x, target_depth))
  }
  rownames(out_mat)<-rownames(in_mat)
  colnames(out_mat)<-colnames(in_mat)
  return(out_mat)
}

##############################################

get_min_count_norm_per_col<-function(in_col_vect1, target_depth){
  loading_factor<-sum(in_col_vect1)/target_depth
  return(in_col_vect1/loading_factor)
}


#' Normalize columns to the min-count of a matrix
#'
#' This function takes in a matrix of counts & will 
#'
#' @param in_mat An input matrix
#' @param target_depth This can either be nothing or a lower cutoff
#' @param return_log If we should log2(n+1) transform after normalizing
#' @param rm_less Will remove cells less than the cutoff rather than throwing an error.
#' @return A matrix of the downsampled in_mat
#' @examples
#'    in_mat<-matrix(rnbinom(1000*20000,.1,c(0.5,0.5,0.5)),ncol=1000,nrow=20000)
#'    out_mat<-min_count_norm_mat(in_mat)
#' @name downsample_mat
#' @export
min_count_norm_mat<-function(in_mat, target_depth=NULL, return_log=FALSE, rm_less=TRUE){
  if (is.null(target_depth)){
    target_depth<-min(colSums(in_mat))
  }
  out_mat <- apply(in_mat, 2, function(x) get_min_count_norm_per_col(x, target_depth))
  if (return_log){
    out_mat<-log2(out_mat+1)
  }
  return(out_mat)
}




