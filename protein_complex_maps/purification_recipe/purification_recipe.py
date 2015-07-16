
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpld3

import logging
import numpy as np
import itertools as it

import argparse
import pickle
import pandas as pd

logging.basicConfig(level = logging.INFO,format='%(asctime)s %(levelname)s %(message)s')

def main():

    parser = argparse.ArgumentParser(description="Create purification recipe")
    parser.add_argument("--input_msds_pickles", action="store", nargs='+', dest="msds_filenames", required=True,
                                help="Filenames of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
    parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True,
                                help="Protein ids in which to anaylze")
    parser.add_argument("--protein_percent_threshold", action="store", dest="protein_percent_threshold", type=float, required=False, default=0.9,
                                help="Percentage of proteins that need to be present in input fractions (expressed as decimal), default 0.9 (ie 90%)")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                                help="Filename of plot")
    parser.add_argument("--plotly", action="store_true", dest="plotly", required=False, default=False,
                                help="Plot to plot.ly")
    parser.add_argument("--fractionation_type", action="store", dest="fractionation_type_file", required=False, default="./fractionation_type.txt",
                                help="File that describes the type of column used for a fractionation experiment")

    args = parser.parse_args()

    pr_obj = PurificationRecipe( args.msds_filenames, args.proteins, args.protein_percent_threshold, args.plot_filename, args.fractionation_type_file)
    pr_obj.create_recipes()
    pr_obj.show_results()
    if args.plot_filename != None:
        if args.plotly:
            pr_obj.plot_results_to_plotly()
        else:
            pr_obj.plot_results()


class PurificationRecipe(object):

    def __init__(self, msds_filenames, proteins, protein_percent_threshold, plot_filename=None, plot_title=None, fractionation_type_file=None, normalize_by_fraction_flag=True, normalize_by_protein_flag=False):
        self.msds_filenames = msds_filenames
        self.input_proteins = proteins
        self.proteins = []
        self.protein_map = dict()
        self.protein_percent_threshold = protein_percent_threshold
        self.normalize_by_fraction_flag = normalize_by_fraction_flag
        self.normalize_by_protein_flag = normalize_by_protein_flag
        self.plot_filename = plot_filename
        self.plot_title = plot_title
        self.fractionation_type_file = fractionation_type_file
        self.plot_colors = ['rbg(240,163,255)','rbg(0,117,220)','rbg(153,63,0)','rbg(76,0,92)','rbg(25,25,25)','rbg(0,92,49)','rbg(43,206,72)','rbg(255,204,153)','rbg(128,128,128)','rbg(148,255,181)','rbg(143,124,0)','rbg(157,204,0)','rbg(194,0,136)','rbg(0,51,128)','rbg(255,164,5)','rbg(255,168,187)','rbg(66,102,0)','rbg(255,0,16)','rbg(94,241,242)','rbg(0,153,143)','rbg(224,255,102)','rbg(116,10,255)','rbg(153,0,0)','rbg(255,255,128)','rbg(255,255,0)','rbg(255,80,5)']
        self.plot_hexcolors = ['#f0a3ff', '#0075dc', '#993f00', '#4c005c', '#191919', '#005c31', '#2bce48', '#ffcc99', '#808080', '#94ffb5', '#8f7c00', '#9dcc00', '#c20088', '#003380', '#ffa405', '#ffa8bb', '#426600', '#ff0010', '#5ef1f2', '#00998f', '#e0ff66', '#740aff', '#990000', '#ffff80', '#ffff00', '#ff5005']

        self.results_list = []
        self.df_dict = dict()
        self.df_index_dict = dict()
        self.frac_type_dict = dict()

        if self.fractionation_type_file != None:
            self.map_fractionation_type()

        #kdrew: initialize
        for prot_id in self.input_proteins:
            self.protein_map[prot_id] = prot_id

        self.create_data_frames()


    def map_fractionation_type(self,):

        ft_file = open(self.fractionation_type_file, "rb")
        for line in ft_file.readlines():
            ls = line.split('\t')

            fname = None
            for msds_filename in self.msds_filenames:
                if ls[0] in msds_filename:
                    fname = msds_filename
                    break

            if fname != None:
                if len(ls) == 2:
                    self.frac_type_dict[fname] = ls[1].strip()
                else:
                    self.frac_type_dict[fname] = None

        print self.frac_type_dict


    def create_data_frames(self,):

        msds_dict= dict()
        for msds_filename in self.msds_filenames:
            msds = pickle.load( open( msds_filename, "rb" ) )
            msds_dict[msds_filename] = msds

            #kdrew: map input protein ids to protein ids used for indexing
            #kdrew: id_dict is the mapping from protein_id -> index
            id_dict = msds.get_id_dict()
            #kdrew: name list will map id -> protein_id used for indexing
            name_list = msds.get_name_list()
            #kdrew: this will fill protein_map with the id used in msds
            for prot_id in self.input_proteins:
                try:
                    self.protein_map[prot_id] = name_list[id_dict[prot_id]]
                except KeyError:
                    continue

        self.proteins = self.protein_map.values()
        print "self.proteins"
        print self.proteins

        #kdrew: store data in data frames
        for msds_filename in msds_dict.keys():
            msds = msds_dict[msds_filename]
            df = msds.get_data_frame()

            #print "filename: %s" % msds_filename

            missing_proteins = [ self.protein_map[i] for i in self.input_proteins if self.protein_map[i] not in df.columns ]

            print "final mapped proteins: %s" % self.proteins
            print "final missing proteins: %s" % missing_proteins

            #kdrew: add in 0.0 rows for missing entries that are still unmapped
            for missing_prot in missing_proteins:
                df[missing_prot] = np.repeat( 0.0, len(df.index) )


            ##kdrew: normalize by fractions (max of all fraction vectors = 1.0)
            #df = df.div( df.max(axis=1), axis=0 )

            if self.normalize_by_fraction_flag:
                #kdrew: normalize by fractions (sum of all fraction vectors = 1.0, i.e. probability)
                df = df.div( df.sum(axis=1), axis=0 )

            if self.normalize_by_protein_flag:
                #kdrew: from a probability perspective, it makes sense to normalize by protein
                #kdrew: unless we want to preserve the total amount of spectral counts as a proxy for yield
                df = df.div( df.sum(axis=0), axis=1 )

            #kdrew: threshold out fraction vectors that do not have enough of the input proteins
            #print self.proteins
            #print df
            bool_vector = df[self.proteins] > 0.0
            sum_vector = bool_vector.sum(axis=1)
            percent_vector = 1.0*(sum_vector)/len(self.proteins)
            #percent_vector = 1.0*(df[self.proteins] > 0.0).sum(axis=1)/len(self.proteins)
            #df = df[percent_vector > self.protein_percent_threshold]


            #kdrew: store list of fractions that pass threshold
            self.df_index_dict[msds_filename] = list(df[percent_vector > self.protein_percent_threshold].index)

            #kdrew: store dataframe
            self.df_dict[ msds_filename ] = df


    def create_recipes(self,):

        #kdrew: for different possible numbers of sequential experiments, focusing on 1 exp or combinations of 2 or more ...
        for r in range(1, len(self.df_dict)+1):
            for exp_list in it.combinations( self.df_dict.keys(), r ):
                self.create_recipe_helper(exp_list)



    def create_recipe_helper(self, exp_list):

        if self.fractionation_type_file != None:
            #kdrew: filter out experiment sets that have multiple of the same fractionation type
            fractionation_types = set()
            for exp in exp_list:
                if self.frac_type_dict[exp] not in fractionation_types and self.frac_type_dict[exp] != None:
                    fractionation_types.add(self.frac_type_dict[exp])
                else:
                    #kdrew: if the fractionation type is already present than do not calculate
                    return

        df_list = [ self.df_dict[exp] for exp in exp_list ]
        df_index_list = [ self.df_index_dict[exp] for exp in exp_list ]

        #kdrew: for different sets of fractions
        for fractions in it.product(*df_index_list):

            #kdrew: initialize vector to be the first fraction
            y_df = df_list[0]
            y_df = y_df.div(y_df.sum(axis=0))
            y_vector = y_df.loc[fractions[0]]
            #kdrew: multiply additional fraction vectors
            for i in range(1,len(fractions)):
                y_df2 = df_list[i]
                y_df2 = y_df2.div(y_df2.sum(axis=0))
                y_vector2 = y_df2.loc[fractions[i]]
                y_vector = y_vector*y_vector2

            yield_mean = y_vector[self.proteins].mean()
            #print "average yield: %s" % (yield_mean)


            #kdrew: combine vectors as distributions across fractions
            df = df_list[0]

            #kdrew: do not do normalization across fractions on the first one because we want the initial concentration to be present
            #df = df.div(df.sum(axis=0),axis=1)
            vector = df.loc[fractions[0]]

            for i in range(1,len(fractions)):
                df = df_list[i]
                df = df.div(df.sum(axis=0),axis=1)
                vector = vector.mul( df.loc[fractions[i]], fill_value = 0.0 )


            score_dict = self.score_fraction( vector )

            pr = PurificationRecipeResult()
            pr.proteins = score_dict['proteins_present']
            #kdrew: protein_percent, purity_percent and yield_percent are not really precents but fractionals.
            #kdrew: The term "fraction" is already used for column separation so keep with percent to avoid confusion.
            pr.protein_percent = score_dict['protein_percent']
            pr.purity_percent = score_dict['purity_percent']
            pr.experiments=exp_list
            pr.fractions=fractions
            pr.yield_percent = yield_mean

            if self.fractionation_type_file != None:
                pr.fraction_types = fractionation_types

            self.results_list.append(pr)




    def show_results(self,):
        sorted_results_list = sorted(self.results_list, key = lambda result : result.purity_percent)
        for result in sorted_results_list:
            print result

    def has_results(self,):
        return len(self.results_list) > 0


    def plot_results(self,):
        css = """
        table
        {
          border-collapse: collapse;
        }
        th
        {
          color: #ffffff;
          background-color: #000000;
        }
        td
        {
          background-color: #FFFFFF;
            padding: 5px;
        }
        table, th, td
        {
          font-family:Arial, Helvetica, sans-serif;
          border: 1px solid grey;
          text-align: center;
        }
        body
        {
          font-family:Arial, Helvetica, sans-serif;
        }

        """

        fig, ax = plt.subplots(subplot_kw=dict(axisbg='#FFFFFF'))
        ax.grid(color='#EEEEEE', linestyle='solid')
        ax.set_title(self.plot_title, size=15)

        ax.set_xlabel("Purity %", size=15)
        ax.set_ylabel("Yield %", size=15)

        traces = []
        fraction_types_set = set( [ tuple(result.fraction_types) for result in self.results_list ] )
        for i, frac_types_tup in enumerate(fraction_types_set):

            #kdrew: make "percents" really percents
            purity_percent = [result.purity_percent*100 for result in self.results_list if tuple(result.fraction_types) == frac_types_tup ]
            yield_percent = [result.yield_percent*100 for result in self.results_list if tuple(result.fraction_types) == frac_types_tup ]
            fraction_text = ["<table> <tr> <td>Fraction Types</td><td> %s </td></tr><tr><td>Fractions</td><td> %s </td></tr><tr><td>Protein Percent</td><td> %s </td></tr></table>" % ("<br>&nbsp&nbsp".join(result.fraction_types),"<br>&nbsp&nbsp".join(result.fractions), result.protein_percent*100) for result in self.results_list if tuple(result.fraction_types) == frac_types_tup ]

            scatter = ax.scatter( purity_percent, yield_percent, c=self.plot_hexcolors[i], s=70, edgecolor=self.plot_hexcolors[i], label=" ".join(frac_types_tup))

            #tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=fraction_text)
            tooltip = mpld3.plugins.PointHTMLTooltip(scatter, labels=fraction_text, css=css)
            mpld3.plugins.connect(fig, tooltip)



        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #legend = ax.legend(loc=2, bbox_to_anchor=(.85, 1),borderaxespad=0.)
        mpld3.save_html(fig, self.plot_filename)

        #fig.savefig(self.plot_filename)


    def plot_results_to_plotly(self,):

        import plotly.plotly as py
        import plotly.graph_objs as go

        traces = []
        fraction_types_set = set( [ tuple(result.fraction_types) for result in self.results_list ] )
        for i, frac_types_tup in enumerate(fraction_types_set):

            #kdrew: make "percents" really percents
            purity_percent = [result.purity_percent*100 for result in self.results_list if tuple(result.fraction_types) == frac_types_tup ]
            yield_percent = [result.yield_percent*100 for result in self.results_list if tuple(result.fraction_types) == frac_types_tup ]
            fraction_text = ["Fraction Types: %s <br> Fractions: %s <br> Protein Percent: %s" % ("<br>".join(result.fraction_types),"<br>".join(result.fractions), result.protein_percent*100) for result in self.results_list if tuple(result.fraction_types) == frac_types_tup ]

            trace1 = go.Scatter(
                x=purity_percent,
                y=yield_percent,
                mode='markers',
                name=", ".join(frac_types_tup),
                text=fraction_text,
                marker=go.Marker(
                    color=self.plot_colors[i],
                    size=12,
                    line=go.Line(
                        color='white',
                        width=0.5
                    )
                )
            )

            traces.append(trace1)

        data = go.Data(traces)
        layout = go.Layout(
            title='Purification Recipe',
            xaxis=go.XAxis(
                title='Purity %',
                showgrid=True,
                zeroline=False
            ),
            yaxis=go.YAxis(
                title='Yield %',
                showline=False
            )
        )
        fig = go.Figure(data=data, layout=layout)
        plot_url = py.plot(fig, filename=self.plot_filename, auto_open=False, world_readable=False)
        print plot_url


    def score_fraction(self, vector ):

        ##kdrew: this scoring function gives a 1*abundance to proteins within in the complex and a -1*abundance to all other proteins
        ##kdrew: the sum is the total score
        #score_mask_1neg1 = [ 1 if i in proteins else -1 for i in vector.index]
        #vector_scored = vector * score_mask_1neg1
        #print "score %s" % (vector_scored.sum())


        #kdrew: count how many original proteins are still present (counts > 0.0)
        protein_count = (vector[self.proteins] > 0.0).sum()
        proteins_present = vector.index[(vector[self.proteins] > 0.0)].tolist()
        protein_percent = (1.0*protein_count)/len(self.proteins)
        #print "protein_count %s / %s = %s" % (protein_count, len(proteins), protein_percent)

        #kdrew: this scoring function gives what percent the final solution will have of proteins in the complex
        score_mask_1_0 = [ 1 if i in self.proteins else 0 for i in vector.index]
        vector_scored = vector * score_mask_1_0
        purity_percent = vector_scored.sum()/vector.sum()
        #print vector_scored
        #print vector
        #print "percent score %s" % (purity_percent)

        #pr = PurificationRecipeResult(proteins=proteins_present, protein_percent=protein_percent, purity_percent=purity_percent )
        score_dict = dict()
        score_dict['proteins_present'] = proteins_present
        score_dict['protein_percent'] = protein_percent
        score_dict['purity_percent'] = purity_percent

        return score_dict


class PurificationRecipeResult(object):

    def __init__(self, experiments=[], fractions=[], fraction_types=[], proteins=[], protein_percent=None, purity_percent=None, yield_percent=None):
        self.experiments = experiments
        self.fractions = fractions
        self.fraction_types = fraction_types
        self.proteins = proteins
        self.protein_percent = protein_percent
        self.purity_percent = purity_percent
        self.yield_percent = yield_percent

    def __str__(self,):
        return "experiments: %s, fractions: %s, fraction_types: %s, protein_percent: %s, purity_percent: %s, yield_percent: %s" % (self.experiments, self.fractions, self.fraction_types, self.protein_percent*100, self.purity_percent*100, self.yield_percent*100)

if __name__ == "__main__":
    main()
