

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
#' @param target_depth This can either be nothing or a 
#' @return A matrix of the downsampled in_mat
#' @examples
#'    in_mat<-matrix(rnbinom(1000*20000,.1,c(0.5,0.5,0.5)),ncol=1000,nrow=20000)
#'    out_mat<-downsample_mat(in_mat)
#' @name downsample_mat
#' @export
downsample_mat<-function(in_mat, target_depth=NULL){
  if (is.null(target_depth)){
    target_depth<-min(colSums(in_mat))
  }
  out_mat <- apply(in_mat, 2, function(x) get_ds_transcript_vect(x, target_depth))
  return(out_mat)
}




