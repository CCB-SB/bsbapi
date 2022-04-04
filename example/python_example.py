from bsbapi import bsbapi
import time
def main():
    fa_file="example.fa.gz"
    #create a job setting various parameters and analysis
    embedding={"embedding_method":"umap", "embedding_params": {
            "N_Neighbors": 15,
            "Init": "spectral",
            "Min_dist":0.15,
            "Metric":"euclidian"
            }}
    custom_annot_file="example_annot.tsv"
    job_id=bsbapi.start_job(fa_file=fa_file, taxonomic=True, functional=True, plasmids=True, binquality=True,embedding=embedding,custom_annot_file=custom_annot_file)
    #wait for computation to finish before retrieving results
    while True:
        states=bsbapi.job_state(job_id)
        ready=True
        for i in states["State"]:
            if i != "SUCCESS" and i != "FAILURE":
                ready=False
                break
        if ready:
            break
        #communicate progress and wait 5 minutes if results are not ready and test again
        print(states)
        time.sleep(5*60) 

    #try to retrieve the results of interest
    binning_result=bsbapi.result_bin(job_id)
    functional_result=bsbapi.result_functional(job_id)
    plasmid_result=bsbapi.result_plasmid(job_id)
    binquality_result=bsbapi.result_binquality(job_id)
    
    #continue procssing as desired from here
    


if __name__ == "__main__":
    main()
