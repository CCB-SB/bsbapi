library(reticulate)
library(R.utils)

# a small hack to find a conda environment that contains the package (works only on linux)
bsbapi_find_command="conda search bsbapi --envs | tail -n1  | awk '{print $NF}'"
bsb_env<-system(bsbapi_find_command,intern=T)
use_condaenv(bsb_env)
bsbapi<-import("bsbapi")$bsbapi

fa_file<-"example.fa.gz"

#create a job setting various parameters and analysis
embedding<-list()
embedding$embedding_method<-"umap"
embedding_params<-list()
embedding_params$N_Neighbors<-15
embedding_params$Init<-"spectral"
embedding_params$Min_dist<-0.15
embedding_params$Metric<-"euclidian"
embedding$embedding_params<- embedding_params
custom_annot_file<-"example_annot.tsv"

job_id<-bsbapi$start_job(fa_file=fa_file, taxonomic=T, functional=T, plasmids=T, binquality=T,embedding=reticulate::dict(embedding),custom_annot_file=custom_annot_file)

#wait for computation to finish before retrieving results
repeat{
	states<-bsbapi$job_state(job_id)
	if(all(states$State=="SUCCESS"| states$State=="FAILURE" )){
	      break
	}
	#wait 5 minutes if results are not ready and test again
	print(states)
	Sys.sleep(5*60)
}

#try to retrieve the results of interest
binning_result=bsbapi$result_bin(job_id)
functional_result=bsbapi$result_functional(job_id)
plasmid_result=bsbapi$result_plasmid(job_id)
binquality_result=bsbapi$result_binquality(job_id)

#continue procssing as desired from here
