
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.protein_util as pu

import protein_complex_maps.complex_map_website.complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Outputs table of protein complexes from database")
    parser.add_argument("--outputfile", action="store", dest="outputfile", required=False, default=None,
                                    help="Filename to store output, default=stdout")
    parser.add_argument("--field_delimiter", action="store", dest="field_delimiter", required=False, default=',',
                                    help="delimiter between fields, default = ,")
    parser.add_argument("--infield_delimiter", action="store", dest="infield_delimiter", required=False, default='\t',
                                    help="delimiter within fields, default = \\t")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    #kdrew: from SysCilia_1stNeighbors_complexmap.cys
    ciliary_complex_keys = set(['2960', '2647', '2844', '2615', '1786', '219', '4478', '4229', '2093', '3646', '491', '2458', '1832', '1774', '1771', '3713', '4383', '2328', '3254', '1549', '1370', '2019', '402', '2400', '933', '280', '2408', '2243', '2465', '2460', '3765', '2451', '4335', '4483', '1490', '1954', '2317', '3409', '3424', '3741', '3567', '4101', '58', '55', '1824', '2992', '3037', '2495', '3033', '928', '1559', '2144', '2140', '920', '1411', '358', '292', '2794', '319', '4646', '2772', '4021', '4325', '4582', '2777', '2823', '3029', '1714', '2939', '2789', '2381', '1257', '4113', '1520', '522', '1402', '1650', '2179', '529', '426', '2445', '1335', '3589', '3951', '4195', '2053', '1135', '2287', '245', '2352', '3262', '385', '243', '370', '103', '1997', '2277', '4429', '904', '3152', '2253', '3612', '4458', '1846', '2570', '2671', '437', '435', '2162', '3611', '2596', '3614', '1222', '2611', '2861', '335', '336', '331', '4223', '2453', '3695', '60', '4629', '1906', '175', '1945', '856', '1981', '1428', '1733', '4530', '859', '3788', '656', '1507', '2368', '2110', '3084', '3890', '3894', '4272', '750', '4278', '2427', '2244', '229', '932', '90', '2188', '95', '161', '3207', '3529', '1991', '1990', '2217', '1188', '3378', '3446', '2891', '4522', '3209', '1106', '3905', '1748', '4501', '1515', '2695', '2693', '2357', '3617', '882', '3420', '2868', '749', '2430', '608', '1723', '3061', '2621', '3066', '3064', '1896', '2236', '4613', '3231', '2133', '4161', '49', '950', '2348', '4168', '3738', '1112', '3974', '324', '1440', '2749', '4589', '3916', '1175', '776', '3912', '4095', '206', '1270', '3101', '4099', '3622', '4255', '4254', '3627', '1178', '3349', '209', '4466', '4461', '4612', '2331', '4151', '4419', '4154', '949', '947', '2233', '3114', '3425', '941', '3196', '689', '4396', '479', '3981', '682'])
    ciliary_protein_keys = set(['132320', '11190', '5727', '27285', '10210', '117177', '55112', '79140', '56890', '347733', '79023', '5902', '9688', '5901', '55764', '54919', '150483', '163786', '9371', '6993', '5991', '51684', '28981', '27077', '1499', '8021', '57728', '11020', '4735', '129401', '26005', '378', '3842', '129880', '26160', '55582', '55081', '989', '80199', '4647', '25886', '55835', '22919', '55722', '64792', '79809', '22994', '51199', '79960', '5311', '1855', '55212', '83894', '200894', '50807', '585', '583', '80184', '65062', '2932', '51175', '4218', '4117', '11127', '23090', '3312', '352909', '10806', '7272', '22981', '7277', '26146', '4750', '4751', '79583', '4291', '27095', '79598', '4952', '4957', '54536', '7428', '9731', '29922', '57560', '9738', '50855', '124602', '56912', '55329', '56288', '27241', '8195', '10426', '92104', '80776', '90410', '84197', '7419', '8100', '53340', '6809', '123811', '10640', '23432', '5108', '6616', '9742', '11336', '23233', '7846', '80173', '5347', '10844', '221421', '27031', '25911', '95681', '7280', '79659', '3064', '150737', '83696', '122664', '27185', '23636', '83853', '23224', '51715', '123016', '374654', '755', '91147', '6102', '55172', '22897', '27229', '57539', '23059', '84790', '112752', '9696', '10464', '5048', '7109', '261734', '54585', '51256', '6787', '85378', '60412', '1808', '8481', '8766', '27152', '51098', '4998', '9662'])

    complexes = db.session.query(cdb.Complex).filter(cdb.Complex.complex_id.in_(ciliary_complex_keys)).all()

    output_str = "complex_id%sgene_ids%sgenenames%sSysCilia_gene_ids%sSysCilia_genenames\n" % (args.field_delimiter,args.field_delimiter,args.field_delimiter,args.field_delimiter)
    for comp in complexes:
        print comp.complex_id
        geneid_str = args.infield_delimiter.join([p.gene_id for p in comp.proteins])
        genename_str = args.infield_delimiter.join([p.genename() for p in comp.proteins])
        ciliary_geneid_str = args.infield_delimiter.join([p.gene_id for p in comp.proteins if p.gene_id in ciliary_protein_keys])
        ciliary_genename_str = args.infield_delimiter.join([p.genename() for p in comp.proteins if p.gene_id in ciliary_protein_keys])
        output_str = output_str + "%s%s%s%s%s%s%s%s%s%s\n" % (comp.complex_id, args.field_delimiter,
                                                    geneid_str,args.field_delimiter,
                                                    genename_str,args.field_delimiter,
                                                    ciliary_geneid_str,args.field_delimiter,
                                                    ciliary_genename_str,args.field_delimiter
                                                    )

    if args.outputfile != None:
        f = open(args.outputfile, "wb")
        f.write(output_str)
        f.close()
    else:
        print output_str


if __name__ == "__main__":
    main()


