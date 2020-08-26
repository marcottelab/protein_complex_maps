
import os
import argparse

def main():

    parser = argparse.ArgumentParser(description="Creates bash script run libsvm for different parameters")
    parser.add_argument("--input_train_files", action="store", dest="train_files", nargs='+', required=True, 
                                    help="Filenames of training files")
    parser.add_argument("--input_leaveout_files", action="store", dest="leaveout_files", nargs='+', required=True, 
                                    help="Filenames of leaveout files, keep in order of matching input_train_files")
    parser.add_argument("--output_script_name", action="store", dest="output_script_name", required=True, 
                                    help="Filename of output script")
    parser.add_argument("--scale_file", action="store", dest="scale_file", required=False, default=None,
                                    help="Filename of libsvm scale parameter file")
    parser.add_argument("--scripts_dir", action="store", dest="scripts_dir", required=False, default="~/scripts/protein_complex_maps/",
                                    help="Directory of protein_complex_maps scripts")
    parser.add_argument("--libsvm_dir", action="store", dest="libsvm_dir", required=False, default="",
                                    help="Directory of svm-scale svm-train and svm-predict")
    parser.add_argument("--c_values", action="store", dest="c_values", nargs='+', type=float, required=False, default=[0.00390625, 0.0078125, 2, 32, 128],
                                    help="List of C parameter values to evaluate, default = 0.00390625 0.0078125 2 32 128")
    parser.add_argument("--gamma_values", action="store", dest="gamma_values", nargs='+', type=float, required=False, default=[0.015625, 0.03125, 0.0625, 0.125, 0.5],
                                    help="List of gamma parameter values to evaluate, default = 0.015625 0.03125 0.0625 0.125 0.5")
    parser.add_argument("--kernel_values", action="store", dest="kernel_values", nargs='+', type=int, required=False, default=[2],
                                    help="List of kernel values 0=linear, 1=poly, 2=rbf (see svm-train docs for more), default = 2 (rbf)")
    parser.add_argument("--id_columns", action="store", dest="id_columns", nargs='+', required=True, 
                                    help="Names of columns that store protein ids")
    parser.add_argument("--split_train_predict_scripts", action="store_true", dest="split_train_predict_scripts", required=False,  default=False,
                                    help="Split train/predict portion into separate scripts for every parameter pair, default = False")
    parser.add_argument("--group_train_predict_scripts", action="store_true", dest="group_train_predict_scripts", required=False,  default=False,
                                    help="Group train commands and predict commands into their own scripts, default = False ")
    parser.add_argument("--features", action="store", dest="features_str", nargs='+', required=False, 
                                    help="Names of features")
    parser.add_argument("--feature_file", action="store", dest="feature_file", required=False, 
                                    help="File containing one feature name per line")
    parser.add_argument("--feature_identifier", action="store", dest="feature_id", required=False, default="sel", 
                                    help="String to ID choice of features in output models and predictions")
    parser.add_argument("--ignore_asserts", action="store_true", dest="ignore_asserts", required=False, default=False,
                                    help="Overrides code to check files are present, this flag is useful when running script generation on a different machine/file system than the compute machine")
 
    args = parser.parse_args()

    if not args.features_str and not args.feature_file:
        print("Provide either a space separated string of features using --features or a file of features 1 per line using --feature_file")
        return
    elif args.features_str and args.feature_file:
        print("Provide either space separated string of features using --features or a file of features 1 per line using --feature_file")
        return
    elif args.features_str:
        feats = args.features_str
        print(feats)
    elif args.feature_file:
        featfile = open(args.feature_file, "r")
        feats = featfile.read().split("\n")
        featfile.close()
        print(feats)


    if args.group_train_predict_scripts:
        feature2libsvm_group_outfile = open("%s.%s.sh" % ('.'.join(args.output_script_name.split('.')[:-1]), "feature2libsvm_group"), "wb")
        scale_group_outfile = open("%s.%s.sh" % ('.'.join(args.output_script_name.split('.')[:-1]), "scale_group"), "wb")
        train_group_outfile = open("%s.%s.sh" % ('.'.join(args.output_script_name.split('.')[:-1]), "train_group"), "wb")
        predict_group_outfile = open("%s.%s.sh" % ('.'.join(args.output_script_name.split('.')[:-1]), "predict_group"), "wb")
        results2pairs_group_outfile = open("%s.%s.sh" % ('.'.join(args.output_script_name.split('.')[:-1]), "results2pairs_group"), "wb")
    else:
        outfile = open(args.output_script_name,"wb")

    for i in xrange(len(args.train_files)):
        train_file = args.train_files[i]
        leaveout_file = args.leaveout_files[i]

        base_train_file = train_file.replace(".lfeatmat", "")
        base_train_file = base_train_file + "." + args.feature_id
        base_leaveout_file = leaveout_file.replace(".lfeatmat", "")
        base_leaveout_file = base_leaveout_file + "." +  args.feature_id

        #kdrew: convert to libsvm format
        if args.group_train_predict_scripts:
            feature2libsvm_outfile = feature2libsvm_group_outfile
        else:
            feature2libsvm_outfile = outfile

        feature2libsvm_outfile.write("#create_cv_commands: converting to libsvm format\n")
        feature2libsvm_outfile.write("python %s/protein_complex_maps/model_fitting/svm_utils/feature2libsvm.py --input_feature_matrix %s --output_filename %s.libsvm --features %s --label_column label --sep ,\n" % (args.scripts_dir, train_file, base_train_file, ' '.join(feats)))
        feature2libsvm_outfile.write("python %s/protein_complex_maps/model_fitting/svm_utils/feature2libsvm.py --input_feature_matrix %s --output_filename %s.libsvm --features %s --label_column label --sep , \n" % (args.scripts_dir, leaveout_file, base_leaveout_file, ' '.join(feats)))

        #kdrew: if no libsvm_dir is passed in then use the one in path, if it is passed in, add a '/' so it is a directory
        libsvm_dir = args.libsvm_dir
        if libsvm_dir != "":
            libsvm_dir = libsvm_dir + "/"

            if not args.ignore_asserts:
                assert os.path.isfile("%ssvm-scale" % libsvm_dir) 
                assert os.path.isfile("%ssvm-train" % libsvm_dir) 
                assert os.path.isfile("%ssvm-predict" % libsvm_dir) 

        #kdrew: scale
        if args.group_train_predict_scripts:
            scale_outfile = scale_group_outfile
        else:
            scale_outfile = outfile
        scale_outfile.write("#create_cv_commands: scaling train and leaveout files\n")
        #kdrew: keep commands on same line separated by ';' so that they are run together in sequential order. The second command relies on the first
        if args.scale_file == None:
            scale_outfile.write("%ssvm-scale -s %s.scale_parameters %s.libsvm > %s_scale.libsvm;" % (libsvm_dir, base_train_file, base_train_file, base_train_file))
            scale_outfile.write("%ssvm-scale -r %s.scale_parameters %s.libsvm > %s_scale.libsvm\n" % (libsvm_dir, base_train_file, base_leaveout_file, base_leaveout_file))
        else:
            scale_outfile.write("%ssvm-scale -r %s.scale_parameters %s.libsvm > %s_scale.libsvm;" % (libsvm_dir, args.scale_file, base_train_file, base_train_file))
            scale_outfile.write("%ssvm-scale -r %s.scale_parameters %s.libsvm > %s_scale.libsvm\n" % (libsvm_dir, args.scale_file, base_leaveout_file, base_leaveout_file))

        for cvalue in args.c_values:
            for gvalue in args.gamma_values:
                for kvalue in args.kernel_values:

                    #kdrew: output plumbing, point output to their proper files
                    if args.split_train_predict_scripts:
                        #kdrew: create separate scripts for each C and gamma parameter pair
                        parameter_outfile = open("%s.kernel%s.C%s.g%s.%s.sh" % ('.'.join(args.output_script_name.split('.')[:-1]), kvalue, cvalue, gvalue, i), "wb")
                        parameter_outfile.write("#create_cv_commands: commands for kernel=%s, C=%s and gamma=%s\n" % (kvalue, cvalue, gvalue))
                        train_file = parameter_outfile
                        predict_file = parameter_outfile
                        results2pairs_file = parameter_outfile
                    elif args.group_train_predict_scripts:
                        train_file = train_group_outfile
                        predict_file = predict_group_outfile
                        results2pairs_file = results2pairs_group_outfile
                    else:
                        outfile.write("#create_cv_commands: commands for kernel=%s, C=%s and gamma=%s\n" % (kvalue, cvalue, gvalue))
                        train_file = outfile
                        predict_file = outfile
                        results2pairs_file = outfile

                    #kdrew: train
                    train_file.write("%ssvm-train -t %s -b 1 -c %s -g %s %s_scale.libsvm %s.libsvm.scale.model_k%s_c%s_g%s\n" % (libsvm_dir, kvalue, cvalue, gvalue, base_train_file, base_train_file, kvalue, cvalue, gvalue))
                    #kdrew: predict
                    predict_file.write("%ssvm-predict -b 1 %s_scale.libsvm %s.libsvm.scale.model_k%s_c%s_g%s %s.libsvm.scale.k%s_c%s_g%s.resultsWprob\n" % (libsvm_dir,base_leaveout_file, base_train_file, kvalue, cvalue, gvalue, base_leaveout_file, kvalue, cvalue, gvalue))
                    #kdrew: pairwise formatting
                    #parameter_outfile.write("#create_cv_commands: reformat results into pairs\n")
                    results2pairs_file.write("python %s/protein_complex_maps/model_fitting/svm_utils/svm_results2pairs.py --input_feature_matrix %s --input_results %s.libsvm.scale.k%s_c%s_g%s.resultsWprob --output_file %s.libsvm.scale.k%s_c%s_g%s.resultsWprob_pairs_noself_nodups.txt --sep , --id_columns %s --label_not0 --add_prob\n" % (args.scripts_dir, leaveout_file, base_leaveout_file, kvalue, cvalue, gvalue, base_leaveout_file, kvalue, cvalue, gvalue, ' '.join(args.id_columns)))

                    if args.split_train_predict_scripts:
                        parameter_outfile.close()

                #kdrew: end of kernel_values for loop
            #kdrew: end of gamma_values for loop
        #kdrew: end of cvalues for loop

        if args.group_train_predict_scripts:
            train_group_outfile.write("\n")
            predict_group_outfile.write("\n")
            results2pairs_group_outfile.write("\n")

        outfile.write("\n")

    #kdrew: end of train_files for loop

    if args.group_train_predict_scripts:
        train_group_outfile.close()
        predict_group_outfile.close()
        results2pairs_group_outfile.close()

    outfile.close()


if __name__ == "__main__":
    main()

