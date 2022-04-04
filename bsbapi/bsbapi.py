import requests
import pandas
import logging
import re
import os
import json

SERVER_URL = 'https://ccb-microbe.cs.uni-saarland.de/busybee-update/api/'
#MB in bytes
MB_200 = 209715200


###########################
### Helper functions
###########################

def response_validity_checker(response):
    if response.status_code==102:
        Exception("Result not ready for retrieval. Please wait.")
    if response.status_code!=200:
        exception_text="Server returned status code:"+str(response.status_code)+"."
        json_response=response.json()
        if "data" in json_response and "error" in json_response["data"]:
            exception_text+= " Error message: "+json_response["data"]["error"]        
        raise Exception(exception_text)
    json_response=response.json()
    return(json_response)

###########################
### API requests
###########################

def result_bin(job_id):
    assert(isinstance(job_id,str))
    response = requests.get(url=SERVER_URL+'/bin', params={"job_id":job_id})
    json_response=response_validity_checker(response)
    if "data_points" in json_response:
        result=pandas.json_normalize(json_response["data_points"])
        return (result)
    raise Exception("Incorrect content in server response")

def result_taxonomy(job_id):
    assert(isinstance(job_id,str))
    response = requests.get(url=SERVER_URL+'/taxonomy', params={"job_id":job_id})
    json_response=response_validity_checker(response)
    if "taxonomy" in json_response:
        result=pandas.DataFrame.from_dict(json_response["taxonomy"],orient="index")
        result.reset_index(inplace=True)
        result=result.rename({"index":"sequence"},axis="columns")
        return (result)
    raise Exception("Incorrect content in server response")


def result_functional(job_id):
    assert(isinstance(job_id,str))
    response = requests.get(url=SERVER_URL+'/functional', params={"job_id":job_id})
    json_response=response_validity_checker(response)
    if "functional" in json_response:
        result={}
        for key,values in json_response["functional"].items():
            result[key]= ", ".join(values)
        result=pandas.DataFrame.from_dict(result,orient="index")
        result.reset_index(inplace=True)
        result.columns=["sequence","annotation(s)"]
        return (result)
    raise Exception("Incorrect content in server response")

def result_plasmid(job_id):
    assert(isinstance(job_id,str))
    response = requests.get(url=SERVER_URL+'/plasmid', params={"job_id":job_id})
    json_response=response_validity_checker(response)
    if "plasmids" in json_response:
        result=pandas.json_normalize(json_response["plasmids"])        
        return (result.drop(columns=["index"]))
    raise Exception("Incorrect content in server response")
    

def result_binquality(job_id):
    assert(isinstance(job_id,str))
    response = requests.get(url=SERVER_URL+'/binquality', params={"job_id":job_id})
    json_response=response_validity_checker(response)
    if "binquality" in json_response:
        result=pandas.json_normalize(json_response["binquality"])        
        return (result)
    raise Exception("Incorrect content in server response")

def job_state(job_id):
    assert(isinstance(job_id,str))
    response = requests.get(url=SERVER_URL+'/job_state', params={"job_id":job_id})
    json_response=response.json()
    if response.status_code==200:
        if "job_status" and "jobs" in json_response:
            result=pandas.DataFrame.from_dict(json_response["jobs"],orient='index')
            result.reset_index(inplace=True)
            result.columns=["Job","State"]
            return(result)
    else:
        exception_text="Server returned status code:"+str(response.status_code)+"."
        if "data" in json_response and "error" in json_response["data"]:
            exception_text+= " Error message: "+json_response["data"]["error"]
        raise Exception(exception_text)
    raise Exception("Unexpected Error")
    

def start_job (
        fa_file,
        custom_annot_file="",
        name="",
        rename_fasta_headers=True,
        taxonomic=False,
        functional=False,
        plasmids=False,
        binquality=False,
        min_contig=500,
        border_min_contig=1000,
        cluster_min_contig=2000,
        kmer_length=5,
        probability=0,
        transformation="clr",
        compression=0,
        seed=0,
        embedding={"embedding_method":"tsne", "embedding_params": {
            "Perplexity": 30,
            "Iterations": 1000,
            "Theta":0.5}},
        clustering={
            "clustering_method": "hdbscan",
            "clustering_params": {
                "minPts": 30
                  }
            }):

    """
    Submits job to the server.
    :param fa_file : Input fasta file
    :param custom_annot_file : Optional tab seperated file containing annotations for each contig
    :param rename_fasta_headers : Boolean value indicating whether fasta header names should be adjusted
    :param name : Optional job name
    :param taxonomic : Boolean value indicating whether taxonomic analysis should be performed
    :param functional : Boolean value indicating whether functional analysis should be performed
    :param plasmids : Boolean value indicating whether plasmid analysis should be performed
    :param binquality : Boolean value indicating whether bin quality analysis should be performed
    :param min_contig : Minimum length of a contig. Must be >=500bp.
    :param border_min_contig : Minimum length of a contig to be considered a border contig
    :param cluster_min_contig : Minimum length of a contig to be considered a cluster contig
    :param kmer_length :  k-mer length. Either 4 or 5
    :param probability : Minimum bin prediction probability
    :param transformation : Normalization method applied to k-mer counts 
    :param compression : Factor by which points should be aggregated
    :param seed : Seed for random number generation
    :param embedding : Dictionary specifying embedding method and necessary parameters
    :param clustering : Dictionary specifying clustering method an necessary paramaters
    :return: A job_id from the server in case of successfull submission
    """
    import time
    from tempfile import NamedTemporaryFile

    if not os.path.isfile(fa_file):
        raise Exception("Input fasta file does not exist")
    if os.path.getsize(fa_file)>MB_200:
        raise Exception("Input fasta file must be below 200MB")
    
    if custom_annot_file!="" and not os.path.isfile(custom_annot_file):
        raise Exception("Input custom annotation file does not exist")
    assert(isinstance(rename_fasta_headers,bool)) 
    assert(isinstance(taxonomic,bool)) 
    assert(isinstance(functional,bool)) 
    assert(isinstance(plasmids,bool))
    assert(isinstance(binquality,bool))
    assert(500<=min_contig)
    assert(min_contig<= border_min_contig)
    assert(border_min_contig<=cluster_min_contig)
    assert(kmer_length<6 and kmer_length>3)
    assert(probability<=1 and probability>=0)
    assert(transformation in ["clr","standard"])
    assert(0<=compression)
    assert(isinstance(seed,int))
    

    PARAMS={            
    "rename_fasta_headers" : rename_fasta_headers,
    "name" : name ,
    "taxonomic" : taxonomic ,
    "functional" : functional ,
    "plasmids" : plasmids ,
    "binquality" : binquality ,
    "min_contig" : min_contig ,
    "border_min_contig" : border_min_contig ,
    "cluster_min_contig" : cluster_min_contig ,
    "kmer_length" : kmer_length ,
    "probability" : probability ,
    "transformation" : transformation ,
    "compression" : compression ,
    "seed" : seed ,
    "embedding" : json.dumps(embedding) ,
    "clustering" : json.dumps(clustering) }



    # search with ifile
    FILES = {'fa_file': open(fa_file, 'rb')}
    if custom_annot_file!="":
        FILES["custom_annot_file"]=custom_annot_file
    response = requests.post(url=SERVER_URL+'/start_job', files=FILES, data=PARAMS)
    response_json=response.json()
    if response.status_code==200:
        if "job_id" in response_json:
            return(response_json["job_id"])
        raise Exception("Incorrect content in server response")
    exception_text="Server returned status code:"+str(response.status_code)+"."
    if "data" in json_response and "error" in json_response["data"]:
        exception_text+= " Error message: "+json_response["data"]["error"]        
    raise Exception(exception_text)
