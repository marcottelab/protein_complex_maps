
import os
import argparse

def main():

    parser = argparse.ArgumentParser(description="Creates bash script to run (scale, train, predict) full model")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of labeled feature matrix, ext = .featmat")
    parser.add_argument("--output_script_name", action="store", dest="output_script_name", required=True, 
                                    help="Filename of output script")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default=',',
                                    help="Field separator, default=,")
    parser.add_argument("--label_column", action="store", dest="label_column", required=False, default='label',
                                    help="Label column, default=label")
    parser.add_argument("--features", action="store", dest="features_str", nargs='+', required=False, 
                                    help="Names of features")
    parser.add_argument("--feature_file", action="store", dest="feature_file", required=False, 
                                    help="File containing one feature name per line")
    parser.add_argument("--scripts_dir", action="store", dest="scripts_dir", required=False, default="~/scripts/protein_complex_maps/",
                                    help="Directory of protein_complex_maps scripts")
    parser.add_argument("--libsvm_dir", action="store", dest="libsvm_dir", required=False, default="",
                                    help="Directory of svm-scale svm-train and svm-predict")
    parser.add_argument("--c_value", action="store", dest="c_value", type=float, required=True, 
                                    help="C parameter value to train svm")
    parser.add_argument("--gamma_value", action="store", dest="gamma_value", type=float, required=True,
                                    help="Gamma parameter value to train svm")
    parser.add_argument("--kernel_value", action="store", dest="kernel_value", type=int, required=True,
                                    help="Kernel values 0=linear, 1=poly, 2=rbf (see svm-train docs for more)")
    parser.add_argument("--numOfThreads", action="store", dest="numOfThreads", type=int, required=False, default=1,
                                    help="Number of threads to use for svm-train")
    parser.add_argument("--id_columns", action="store", dest="id_columns", nargs='+', required=True, 
                                    help="Names of columns that store protein ids")
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

    outfile = open(args.output_script_name,"wb")
    
    outfile.write("#kdrew: trim feature matrix to just labeled entries\n")
    train_only_featmat_filename =  "%s.only.featmat" % os.path.basename(args.feature_matrix).replace(".featmat","")
    trim_command = "python %s/protein_complex_maps/model_fitting/svm_utils/trim_unlabeled_featmat.py --input_feature_matrix %s --output_filename %s --label_column %s --sep %s\n" % (args.scripts_dir, args.feature_matrix, train_only_featmat_filename, args.label_column, args.sep)
    outfile.write(trim_command)
    outfile.write("\n")

    outfile.write("#kdrew: converting labeled only feature matrix to libsvm format\n")
    train_only_libsvm_filename = train_only_featmat_filename.replace(".featmat",".libsvm.txt")
    only_libsvm_fmt_command = "python %s/protein_complex_maps/model_fitting/svm_utils/feature2libsvm.py --input_feature_matrix %s --output_filename %s --features %s --label_column %s --sep %s\n" % (args.scripts_dir, train_only_featmat_filename, train_only_libsvm_filename, ' '.join(feats), args.label_column, args.sep)
    outfile.write(only_libsvm_fmt_command)
    outfile.write("\n")

    outfile.write("#kdrew: converting full feature matrix to libsvm format\n")
    train_libsvm_filename = "%s.libsvm.txt" % os.path.basename(args.feature_matrix).replace(".featmat","")
    libsvm_fmt_command = "python %s/protein_complex_maps/model_fitting/svm_utils/feature2libsvm.py --input_feature_matrix %s --output_filename %s --features %s --label_column %s --sep %s\n" % (args.scripts_dir, args.feature_matrix, train_libsvm_filename, ' '.join(feats), args.label_column, args.sep)
    outfile.write(libsvm_fmt_command)
    outfile.write("\n")

    #kdrew: if no libsvm_dir is passed in then use the one in path, if it is passed in, add a '/' so it is a directory
    libsvm_dir = args.libsvm_dir
    if libsvm_dir != "":
        libsvm_dir = libsvm_dir + "/"

        if not args.ignore_asserts:
            assert os.path.isfile("%ssvm-scale" % libsvm_dir) 
            assert os.path.isfile("%ssvm-train" % libsvm_dir) 
            assert os.path.isfile("%ssvm-predict" % libsvm_dir) 

    outfile.write("#kdrew: set number of threads for libsvm\n")
    outfile.write("export OMP_NUM_THREADS=%s\n" % args.numOfThreads)
    outfile.write("\n")

    outfile.write("#kdrew: scale labeled only feature matrix \n")
    train_only_scale_parameter_filename = train_only_libsvm_filename.replace(".libsvm.txt",".scale_parameters")
    train_only_scale_filename = train_only_libsvm_filename.replace(".libsvm.txt",".libsvm.scale.txt")
    train_only_scale_command = "%ssvm-scale -s %s %s > %s\n" % (libsvm_dir, train_only_scale_parameter_filename, train_only_libsvm_filename, train_only_scale_filename)
    outfile.write(train_only_scale_command)
    outfile.write("\n")

    outfile.write("#kdrew: scale full feature matrix by scale parameters used for training\n")
    train_scale_filename = train_libsvm_filename.replace(".libsvm.txt",".libsvm.scale.txt")
    train_scale_command = "%ssvm-scale -r %s %s > %s\n" % (libsvm_dir, train_only_scale_parameter_filename, train_libsvm_filename, train_scale_filename)
    outfile.write(train_scale_command)
    outfile.write("\n")

    outfile.write("#kdrew: commands for training model given specific C and gamma \n")
    train_only_model_filename = "%s_c%s_g%s" % (train_only_scale_filename.replace(".scale.txt",".scale.model"), args.c_value, args.gamma_value)
    train_model_command = "%ssvm-train -b 1 -t %s -c %s -g %s %s %s\n" % (libsvm_dir, args.kernel_value, args.c_value, args.gamma_value, train_only_scale_filename, train_only_model_filename)
    outfile.write(train_model_command)
    outfile.write("\n")

    outfile.write("#kdrew: predict entire feature matrix\n")
    resultsWprob_filename = "%s_c%s_g%s.resultsWprob" % (train_scale_filename.replace(".txt",""), args.c_value, args.gamma_value)
    predict_command = "%ssvm-predict -b 1 %s %s %s\n" % (libsvm_dir, train_scale_filename, train_only_model_filename, resultsWprob_filename)
    outfile.write(predict_command)
    outfile.write("\n")
    
    outfile.write("#kdrew: convert svm results to pairs with probability file\n")
    pairsWprob_filename = resultsWprob_filename.replace("resultsWprob","pairsWprob")
    results2pairs_command = "python %s/protein_complex_maps/model_fitting/svm_utils/svm_results2pairs.py --input_feature_matrix %s --input_results %s --output_filename %s --sep %s --id_columns %s --add_prob --nofilter_label\n" % (args.scripts_dir, args.feature_matrix, resultsWprob_filename, pairsWprob_filename, args.sep, ' '.join(args.id_columns))
    outfile.write(results2pairs_command)
    outfile.write("\n")

    outfile.close()


if __name__ == "__main__":
    main()

