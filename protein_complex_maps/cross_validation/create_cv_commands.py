
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
    parser.add_argument("--id_columns", action="store", dest="id_columns", nargs='+', required=True, 
                                    help="Names of columns that store protein ids")
    parser.add_argument("--split_train_predict_scripts", action="store_true", dest="split_train_predict_scripts", required=False,  default=False,
                                    help="Split train/predict portion into separate scripts for every parameter pair")
    parser.add_argument("--features", action="store", dest="features_str", nargs='+', required=False, 
                                    help="Names of features")
    parser.add_argument("--feature_file", action="store", dest="feature_file", required=False, 
                                    help="File containing one feature name per line")
 
    parser.add_argument("--feature_identifier", action="store", dest="feature_id", required=False, default="sel", 
                                    help="String to ID choice of features in output models and predictions")
 
    args = parser.parse_args()

    outfile = open(args.output_script_name,"wb")

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



    for i in xrange(len(args.train_files)):
        train_file = args.train_files[i]
        leaveout_file = args.leaveout_files[i]

        base_train_file = train_file.replace(".lfeatmat", "")
        base_train_file = base_train_file + "." + args.feature_id
        base_leaveout_file = leaveout_file.replace(".lfeatmat", "")
        base_leaveout_file = base_leaveout_file + "." +  args.feature_id
        #kdrew: convert to libsvm format
        outfile.write("#create_cv_commands: converting to libsvm format\n")
        outfile.write("python %s/protein_complex_maps/svm_utils/feature2libsvm.py --input_feature_matrix %s --output_filename %s.%s.libsvm --features %s --label_column label --sep ,\n" % (args.scripts_dir, train_file, base_train_file, ' '.join(feats)))
        outfile.write("python %s/protein_complex_maps/svm_utils/feature2libsvm.py --input_feature_matrix %s --output_filename %s.%s.libsvm --features %s --label_column label --sep , \n" % (args.scripts_dir, leaveout_file, base_leaveout_file, ' '.join(feats)))

        #kdrew: scale
        outfile.write("#create_cv_commands: scaling train and leaveout files\n")
        if args.scale_file == None:
            outfile.write("%s/svm-scale -s %s.scale_parameters %s.libsvm > %s_scale.libsvm\n" % (args.libsvm_dir, base_train_file, base_train_file, base_train_file))
            outfile.write("%s/svm-scale -r %s.scale_parameters %s.libsvm > %s_scale.libsvm\n" % (args.libsvm_dir, base_train_file, base_leaveout_file, base_leaveout_file))
        else:
            outfile.write("%s/svm-scale -r %s.scale_parameters %s.libsvm > %s_scale.libsvm\n" % (args.libsvm_dir, args.scale_file, base_train_file, base_train_file))
            outfile.write("%s/svm-scale -r %s.scale_parameters %s.libsvm > %s_scale.libsvm\n" % (args.libsvm_dir, args.scale_file, base_leaveout_file, base_leaveout_file))

        for cvalue in args.c_values:
            for gvalue in args.gamma_values:
                #kdrew: create separate scripts for each C and gamma parameter pair
                if args.split_train_predict_scripts:
                    parameter_outfile = open("%s.C%s.g%s.%s.sh" % ('.'.join(args.output_script_name.split('.')[:-1]), cvalue, gvalue, i), "wb")
                else:
                    parameter_outfile = outfile

                parameter_outfile.write("#create_cv_commands: commands for C=%s and gamma=%s\n" % (cvalue, gvalue))

                #kdrew: train
                parameter_outfile.write("%s/svm-train -b 1 -c %s -g %s %s_scale.libsvm %s.libsvm.scale.model_c%s_g%s\n" % (args.libsvm_dir, cvalue, gvalue, base_train_file, base_train_file, cvalue, gvalue))
                #kdrew: predict
                parameter_outfile.write("%s/svm-predict -b 1 %s_scale.libsvm %s.libsvm.scale.model_c%s_g%s %s.libsvm.scale.c%s_g%s.resultsWprob\n" % (args.libsvm_dir,base_leaveout_file, base_train_file, cvalue, gvalue, base_leaveout_file, cvalue, gvalue))
                #kdrew: pairwise fle
                parameter_outfile.write("#create_cv_commands: reformat results into pairs\n")
                parameter_outfile.write("python %s/protein_complex_maps/svm_utils/svm_results2pairs.py --input_feature_matrix %s --input_results %s.libsvm.scale.c%s_g%s.resultsWprob --output_file %s.libsvm.scale.c%s_g%s.resultsWprob_pairs_noself_nodups.txt --sep , --id_columns %s --label_not0 --add_prob\n" % (args.scripts_dir, leaveout_file, base_leaveout_file, cvalue, gvalue, base_leaveout_file, cvalue, gvalue, ' '.join(args.id_columns)))

                if args.split_train_predict_scripts:
                    parameter_outfile.close()


        outfile.write("\n")

if __name__ == "__main__":
    main()

