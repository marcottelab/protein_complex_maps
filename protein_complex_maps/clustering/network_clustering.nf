#!/usr/bin/env nextflow

//kdrew: original commandline for network_clustering.py
//python ~/scripts/protein_complex_maps/protein_complex_maps/clustering/network_clustering.py --input_network ../orig9k_bioplex2_hygeo_bioid_hygeo_boldt_apms_hygeo_treiber_hygeo_wgt2_youn_hygeo_trimCols_groupbyMean.train_labeled.libsvm.scale_c512.0_g0.001953125.noTR.pairsWprob --random_seed 1235 --procs 30 --ppi_threshold_score 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 --twostep_combination clusterone mcl --output_file orig9k_bioplex2_bioid_boldt_treiber_youn_c512_g0.001953125.noTR.cluster_train.psweep.startii1701_20190711.txt --temp_dir /dev/shm/tmp/ --mcl_inflation 1.2 2 3 4 5 7 9 11 15 --clusterone_density 0.1 0.2  --clusterone_max_overlap 0.6  --clusterone_jar /home/kdrew/hopper/programs/clusterone/cluster_one-1.0.jar  --trim2threshold --starting_id_num 1701 &> orig9k_bioplex2_bioid_boldt_treiber_youn_c512_g0.001953125.noTR.cluster_train.psweep.startii1701_20190711.out


//kdrew: parameter definitions
params.input_network = 'pairsWprob'
params.random_seed = 1235

params.ppi_threshold_score = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
params.mcl_inflation = [1.2, 2, 3, 4, 5, 7, 9, 11, 15]
params.clusterone_density = [0.1, 0.2, 0.3, 0.4]
params.clusterone_max_overlap = [0.6, 0.7]
params.twostep_combination = 'clusterone mcl'

params.additional_args = '--trim2threshold'

params.clusterone_jar = '/home/kdrew/hopper/programs/clusterone/cluster_one-1.0.jar'
params.temp_dir = './tmp'

params.starting_id_num = 0
params.protein_complex_maps = '/protein_complex_maps/protein_complex_maps/'

params.results_path = "./results"
params.output_file = 'clustering_psweep.txt'

//kdrew: setup channels for each input parameter (currently only a subset of all possible but these are the main ones)
ppi_threshold_score_channel = Channel.from(params.ppi_threshold_score)
mcl_inflation_channel = Channel.from(params.mcl_inflation)
clusterone_density_channel = Channel.from(params.clusterone_density)
clusterone_max_overlap_channel = Channel.from(params.clusterone_max_overlap)
//kdrew: create all possible combinations of input parameters for full parameter sweep
parameters = ppi_threshold_score_channel.combine(mcl_inflation_channel).combine(clusterone_density_channel).combine(clusterone_max_overlap_channel)

parameters.into { parameters_channel; parameters_channel2}

//kdrew: determine how many total parameter combinations there are
//kdrew: for whatever reason it was difficult to query parameter channel to determine how many total elements were in it, probably a better way out there
parameters_len = params.ppi_threshold_score.size()
parameters_len = parameters_len * params.mcl_inflation.size()
parameters_len = parameters_len * params.clusterone_density.size()
parameters_len = parameters_len * params.clusterone_max_overlap.size()


//kdrew: calculate the starting id and ending id
start_num = params.starting_id_num
end_num = parameters_len + params.starting_id_num

//kdrew: generate all ids within range
parameter_ids = Channel.from( start_num .. end_num )

//kdrew: link ids to parameter channel
parameter_tuples = parameter_ids.merge( parameters_channel)

process setupParameterCombinations {

    tag "${parameter_tuple[0]}"
    publishDir "${params.results_path}"

    input:
    val parameter_tuple from  parameter_tuples

    output:
    file '*.txt' into cluster_files
    file '*.params' into parameter_files
    stdout result

    """
    echo ${parameter_tuple[0]}
    python ${params.protein_complex_maps}/clustering/network_clustering.py --input_network ${params.input_network} --output_file ${params.output_file} --random_seed ${params.random_seed} --temp_dir ${params.temp_dir} --procs 1 --ppi_threshold_score ${parameter_tuple[1]} --twostep_combination ${params.twostep_combination} --mcl_inflation ${parameter_tuple[2]} --clusterone_density ${parameter_tuple[3]}  --clusterone_max_overlap ${parameter_tuple[4]}  --clusterone_jar ${params.clusterone_jar}  ${params.additional_args} --starting_id_num ${parameter_tuple[0]} 

    ofile=${params.output_file}
    python ${params.protein_complex_maps}/preprocessing_util/complexes/complex_merge.py --cluster_filename \${ofile/txt/ii${parameter_tuple[0]}.txt} --output_filename \${ofile/txt/ii${parameter_tuple[0]}.nr.txt} --merge_threshold 1.0
    mv \${ofile/txt/params} \${ofile/txt/ii${parameter_tuple[0]}.params} 
    """
}

process combineParams {

    publishDir "${params.results_path}"

    input:
    file param_file_list from parameter_files.collect()

    output:
    file "${params.output_file}.params" into param_file
    stdout result2

    """
    echo $param_file_list
    head -n 1 ${param_file_list[0]} > ${params.output_file}.params
    for i in $param_file_list; do tail -n 1 \$i >> ${params.output_file}.params; done
    """
}


result.subscribe {
    println it
}

result2.subscribe {
    println it
}

