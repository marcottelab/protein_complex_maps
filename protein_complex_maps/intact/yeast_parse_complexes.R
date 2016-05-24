
library("RpsiXML")

complex_xml_files = list.files("/home/kdrew/data/intact/20160408/Saccharomyces_cerevisiae/", full.names=TRUE)


#intactComplexSet <- parsePsimi25Complex("/home/kdrew/data/intact/20160408/Saccharomyces_cerevisiae/exocyst_yeast.xml", INTACT.PSIMI25)
#names(interactors(intactComplexSet))

complex_list = lapply(complex_xml_files, parsePsimi25Complex, INTACT.PSIMI25)
interactor_list = lapply(complex_list, interactors)
name_list = lapply(interactor_list, names)


for(comp in name_list)
{
    cat(comp)
    cat("\n")
}

