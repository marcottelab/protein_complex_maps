
This is out of date


def make_annotation_table(annot_filename):
        annotation_table = pd.DataFrame(pd.read_table("/project/cmcwhite/protein_complex_maps/protein_complex_maps/db_data/geneid_uniprot_mapping.txt", sep="\t"))
        a= annotation_table[['gene_id', 'genename']]
        b = pd.DataFrame(a.gene_id.str.split(',').tolist(), index=a.genename).stack()
        b = b.reset_index()[[0, 'genename']] # var1 variable is currently labeled 0
        b.columns = ['gene_id', 'genename']
        print(b)
        #d = c.groupby('gene_id').head(1) 
        annotation_table = b.groupby('gene_id').head(1)

        print(annotation_table)



        annotation_table['gene_id'] = annotation_table['gene_id'].astype(int)

        print(annotation_table)

        return annotation_table


if __name__ == "__main__":
   make_annotation_table()
