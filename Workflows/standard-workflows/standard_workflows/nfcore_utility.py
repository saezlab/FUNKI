#from re import S
import sys, yaml, dill, re
#sys.path.insert(1, '/Users/hanna/Documents/projects/Workflows/Python/scUtilities/v4')
#from scUtilities import analysis_params
#import sc_analysis_loops as scl
from copy import deepcopy
from pathlib import Path
from os.path import exists
from os import path, makedirs
import collections
import scanpy as sc, numpy as np, decoupler as dc, matplotlib.pyplot as plt, seaborn as sns, matplotlib as mpl, pandas as pd
from .sc_analysis_baseclass import AnalysisI
from .sc_analysis_baseclass import Baseanalysis


class NfCore(AnalysisI):
    """ This class prepares all needed files for running a nf-core pipeline. """
    def __init__(self):
        super().__init__()


    def prepare_run(self, pipeline_name:str):
        """_summary_

        Args:
            pipeline_name (_type_): name of the nf-core pipeline
        """
        pipeline = f"nfcore_{pipeline_name}"
        #pipeline_version = self.analysis_params['nfcore']['pipeline'][pipeline_name]['version']           # used for naming the folder
        pipeline_version_name = self.analysis_params['nfcore']['pipeline'][pipeline_name]['version_name'] # used in the nf-core config file and must match the version given by nf-core
        self.paths['nfcore']['nfcore_path'] = path.join(pipeline, pipeline_version_name)
        self.paths['nfcore']['exec_env_nfcore_path'] = path.join(self.paths['exec_env_data_path'], self.paths['nfcore']['nfcore_path'])
        self.paths['full_local_nfcore_path'] = path.join(self.paths['datapath_tmp'], self.paths['nfcore']['nfcore_path'])
        self.paths['full_raw_path'] = path.join(self.paths['datapath_tmp'], self.paths['rawpath'])

    def create_sample_sheet(self, pipeline_name:str):
        """Create sample sheet for nf-core pipeline from fastq files. 
            The files must be in the following format: 
            rawpath must lead to folders named by sample. 
            The folders must contain forward and reverse reads in the format *_1.fq.gz.

        Args:
            pipeline_name (str): name of nf-core pipeline
        """
        import os, pathlib
        import pandas as pd

        if pipeline_name == 'rnaseq':
            samplesheet = pd.DataFrame()
            samplesheet["sample"] = [p.name for p in pathlib.Path(self.paths['full_raw_path']).absolute().glob('*') if os.path.isdir(p)]
            samplesheet["fastq_1"] = [file for file in pathlib.Path(self.paths['full_raw_path']).absolute().glob('*/*_1.fq.gz')]#[os.path.abspath(name) for name in os.listdir(path) if os.path.isdir(name)]
            samplesheet["fastq_2"] = [file for file in pathlib.Path(self.paths['full_raw_path']).absolute().glob('*/*_2.fq.gz')]
            samplesheet["strandedness"] = "auto"
            samplesheet["fastq_1"] = path.join(self.paths['exec_env_data_path'], self.paths['rawpath'],'/'.join(samplesheet["fastq_1"][0].parts[-2:]))
            samplesheet["fastq_2"] = path.join(self.paths['exec_env_data_path'], self.paths['rawpath'],'/'.join(samplesheet["fastq_2"][0].parts[-2:]))
            #samplesheet["fastq_1"] = self.paths['exec_env_data_path'] + samplesheet["fastq_1"].astype(str)
            self.paths["nfcore"]["samplesheet_path"] = path.join(self.paths['metadatapath'], f"{self.paths['nfcore']['samplesheet_name']}.csv")
            samplesheet.to_csv(self.paths["nfcore"]["samplesheet_path"], index=False)


    def get_mouse_ref(self):
        """Get reference genome paths

        Returns:x
            string, string, string: fasta, gtf, blacklist paths
        """
        fasta = f"{self.paths['references_path']}/{self.analysis_params['organism']}/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
        gtf = f"{self.paths['references_path']}/{self.analysis_params['organism']}/Mus_musculus.GRCm39.109.gtf.gz"
        blacklist = f"{self.paths['references_path']}/{self.analysis_params['organism']}/mm39.excluderanges_lexicographical_nochr202122_nochr.bed"
        return fasta, gtf, blacklist

    def get_reference(self):
        if self.analysis_params['organism'] == "mouse":
            fasta, gtf, blacklist = self.get_mouse_ref()
            # this is the same for human so move this outside of ifelse when human is available, too. And use eval(get_{organism}_ref())
            self.paths['fasta'] = fasta
            self.paths['gtf'] = gtf
            self.paths['blacklist'] = blacklist
        elif self.analysis_params['organism'] == "human":
            print("Please provide a reference genome for this organism.")
        else:
            print("Please provide a reference genome for this organism.")


    def get_nfcore_params_rnaseq(self, foldername, pipeline)->dict:
        """
        nf-core params that are not used: 
        - BBSplit params
        - "genome": "mm10" # name of iGenomes reference
        - "additional_fasta" # containing spike-in sequences
        - "transcript_fasta" # only needed when salmon pseudo-alignment shall not use the provided genome fasta and gtf files.
        - "salmon_index" # see "transcript_fasta"
        - "--skip_alignment" # only needed when salmon shall be run in isolation
        - "extra_star_align_args" 
        - "extra_trimgalore_args"
        - "extra_fastp_args"

        Args:
            foldername (_type_): _description_
            pipeline (_type_): _description_

        Returns:
            dict: nf-core pipeline parameters 
        """
        folderpath = f"{self.paths['nfcore']['exec_env_nfcore_path']}/{foldername}/"
        samplesheet_path_ext = path.join(self.paths["exec_env_data_path"], "metadata", f"{self.paths['nfcore']['samplesheet_name']}.csv")

        nfcore_params = {
            "input": samplesheet_path_ext,
            "outdir": f"{folderpath}/",
            "email": self.analysis_params['cluster']['email'],
            "multiqc_title": f"{self.analysis_params['proj_id']} {pipeline}",
            # reference genome
            "fasta": self.paths['fasta'],
            "gtf": self.paths['gtf'],
            "gencode": False,            # if GTF annotation is in GENCODE format
            "gtf_extra_attributes": "gene_name",
            "gtf_group_features": "gene_id",
            "featurecounts_group_type": "gene_biotype",
            "featurecounts_feature_type": "exon",
            "save_reference": True,         # STAR index
            # adapter trimming
            "trimmer": "trimgalore",
            "min_trimmed_reads": 10000,
            "skip_trimming": False,
            "save_trimmed": True,
            # alignment
            "aligner": "star_rsem",
            "pseudo_aligner": "salmon",
            "bam_csi_index": False,         # required for genomes with larger chromosome sizes like wheat
            "star_ignore_sjdbgtf": False,
            "min_mapped_reads": 5,
            "extra_salmon_quant_args": "--seqBias --gcBias",
            "stringtie_ignore_gtf": False,   # StringTie is not that popular? Using this option expands the results but not necessarily with true data.
            "rseqc_modules": "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication,tin"
            }
        if self.analysis_params['hasUmi'] == 'true':
            nfcore_params.update({
                "with_umi": true,
                "umitools_extract_method": "regex",
                "skip_umi_extract": true,
                "umitools_bc_pattern": "NNNNNN",
                "umitools_bc_pattern2": "NNNNNN",
                "umi_discard_read": 1,
                "umitools_umi_separator": ";",
                "umitools_grouping_method": "cluster",
                "umitools_dedup_stats": true,
                "save_umi_intermeds": true
            })
        return nfcore_params

    def get_nfcore_params_atacseq(self, foldername, pipeline)->dict:
        folderpath = f"{self.paths['nfcore']['exec_env_nfcore_path']}/{foldername}/"
        samplesheet_path_ext = path.join(self.paths["exec_env_data_path"], "metadata", f"{self.paths['nfcore']['samplesheet_name']}.csv")

        return {"input": samplesheet_path_ext,
            "outdir": f"{folderpath}/",
            "multiqc_title": f"{self.analysis_params['proj_id']} {pipeline}",
            "email": self.analysis_params['cluster']['email'],
            # sequencing
            "with_control": False,
            "fragment_size": 300,
            "read_length": 50,
            # reference genome
            "fasta": self.paths['fasta'],
            "gtf": self.paths['gtf'],
            "blacklist": self.paths['blacklist'],
            "save_reference": True,
            # mitochondria
            "mito_name": "chrM",
            "keep_mito": False,
            # trimming
            "trim_nextseq": 20,
            "skip_trimming": False,
            "save_trimmed": False,
            # alignment
            "aligner": "bwa",
            "keep_dups": False,
            "keep_multi_map": False,
            "skip_merge_replicates": False,
            "save_align_intermeds": False,
            "save_unaligned": False,
            # peak calling
            "narrow_peak": False,
            "broad_cutoff": 0.1,
            "min_reps_consensus": 1,
            "save_macs_pileup": True,
            "skip_peak_qc": False
            }

    def get_config(self, pipeline_name)->str:
        """Get content of config file

        Returns:
            str: filecontent
        """
        conf = self.analysis_params['nfcore']['pipeline'][pipeline_name]
        return f"params {{\n\
max_cpus   = {conf['max_cpus']}\n\
max_memory = {conf['max_memory']}\n\
max_time   = {conf['max_time']}\n\
}}\n\n\
process{{\n\
{conf['process_config']}\n\
executor = '{conf['executor']}'\n\
}}"

    def get_run_sh(self, foldername:str, pipeline_name:str)->str:
        """Create content for run.sh file.

        Args:
            foldername (str): name of destination folder

        Returns:
            str: run.sh content
        """
        folderpath = f"{self.paths['nfcore']['exec_env_nfcore_path']}/{foldername}/"

        pbs = f"#PBS -N {pipeline_name}_{foldername}\n\
#PBS -l walltime=06:00:00\n\
#PBS -l nodes=1:ppn=2\n\
# Ensures that the Linux environment for the job is the same as the one we're working in:\n\
#PBS -V\n\
# stderr redirection\n\
#PBS -e {folderpath}/errors.err\n\
# stdout redirection\n\
#PBS -o {folderpath}/log.log\n\
"

        slurm = f"#SBATCH --job-name={pipeline_name}_{foldername}\n\
#SBATCH --partition=single\n\
#SBATCH --time=06:00:00\n\
#SBATCH --ntasks=1\n\
#SBATCH --ntasks-per-node=2\n\
#SBATCH --cpus-per-task=1 \n\
#SBATCH --nodes=1\n\
#SBATCH --output={folderpath}/log.log\n\
"
        
        return f"#!/bin/sh\n\
{eval(self.analysis_params['nfcore']['pipeline'][pipeline_name]['executor'])}\n\
\n\
PATH=$PATH:{self.paths['exec_env_path']}\n\
export HOME={self.paths['exec_env_path']}\n\
export NXF_SINGULARITY_CACHEDIR=$HOME/singularityImages/\n\
export NXF_WORK=$HOME/nxf/\n\
module load {self.analysis_params['nfcore']['pipeline'][pipeline_name]['module_java']}\n\
{self.analysis_params['nfcore']['pipeline'][pipeline_name]['load_modules']}\n\
exec_env_nfcore_path={folderpath}\n\
{self.paths['nfcore']['nextflow_executable']} run nf-core/{pipeline_name} -r {self.analysis_params['nfcore']['pipeline'][pipeline_name]['version']} -profile {self.analysis_params['nfcore']['pipeline'][pipeline_name]['profile']} -c $exec_env_nfcore_path/config.config -resume -params-file $exec_env_nfcore_path/params.json\n\
"


    def init_nfcore_run(self, update_content:dict, foldername:str):
        """Creates the three needed files for running a nf-core pipeline. 

        Args:
            update_content (dict): 'analysis_params' get updated with 'update_content'
            foldername (str): foldername for the three files and the results of the nf-core pipeline
        """
        import json

        match self.seq_type:
            case 'bulkRNA':
                pipeline_name = 'rnaseq'
            case 'bulkHic':
                pipeline_name = 'hic'
            case 'bulkAtac':
                pipeline_name = 'atacseq'

        run_sh = self.get_run_sh(foldername, pipeline_name)
        config = self.get_config(pipeline_name)
        analysis_params = eval(f"self.get_nfcore_params_{pipeline_name}('{foldername}', '{self.paths['nfcore']['samplesheet_name']}')")
        analysis_params.update(update_content)

        # prepare folder
        path_to_folder = f"{self.paths['full_local_nfcore_path']}/{foldername}/"
        if not path.exists(path_to_folder):
            makedirs(path_to_folder)
            
        # params
        params = json.dumps(analysis_params, indent=4)
        with open(f'{path_to_folder}/params.json', 'w+') as file: 
            file.write(params)

        # config
        with open(f'{path_to_folder}/config.config', 'w+') as file: 
            file.write(config)

        # run.sh
        with open(f'{path_to_folder}/run.sh', 'w+') as file:
            file.write(run_sh)
