import pandas as pd
from bokeh.charts import Line, show, save, output_file, vplot, hplot
from bokeh.charts import defaults
from bokeh.layouts import gridplot
from bokeh.plotting import figure
from bokeh.models import HoverTool, Span, Label 
from bokeh.embed import components


import feather 
defaults.width = 500
defaults.height = 100

#The other plots are linked to this first one for rescaling and zooming and such

#Adding a span
#http://stackoverflow.com/questions/38353986/adding-a-label-to-a-span-in-bokeh-0-12






def add_spans(p1):
    major = [339]
    minor =[84, 171, 248, 391]


    for loc in minor:
         p1.add_layout(Span(location=loc,
                              dimension='height', line_color='green',
                              line_dash='dashed', line_width=1)) 
    for loc in major:
         p1.add_layout(Span(location=loc,
                              dimension='height', line_color='green',
                               line_width=1)) 
    return p1

def plot_attributes(p1, df, columnID):
    p1.min_border = 0
    p1.outline_line_color = None
    p1.xgrid.grid_line_color = None
    p1.ygrid.grid_line_color = None
    p1.yaxis.axis_label = columnID
    p1.axis.visible = False
    p1.background_fill_alpha = 0.9
    xmax = len(df.index)
    txt = Label(x=0.01, y=0.01, text=columnID, text_color="black", text_font_size="8pt")
    p1.add_layout(txt)
 


    hover = p1.select_one(HoverTool)
    #hover.point_policy = "snap_to_data"
    hover.line_policy="nearest"
    hover.mode = "vline"
    #print df.head
    hover.tooltips = [
    ("bla", "$index"),("Experiment", "@experiment"), ("Species", "@spec")]
    p1 = add_spans(p1)
    return p1

def master_sparkline(df, tools, columnID, color):
    #print df
    # build the line plots
    print columnID
    #print df[columnID]
    p1=figure(width=800, height=60, tools=tools)
    p1.line('index', columnID, source=df, line_color=color)

    p1 = plot_attributes(p1, df, columnID)
        #p1.left[0].visible = False
    #p1.below[0].visible = False
    return p1

def header_line(masterplot, tools, species_longform):
   
    header = figure(width=800,height=50, tools=tools,  x_range=masterplot.x_range, y_range=masterplot.y_range)
    header.min_border = 0
    header.outline_line_color = None
    header.xgrid.grid_line_color = None
    header.ygrid.grid_line_color = None
    header.axis.visible = False
    txt = Label(x=0, y=1, text="InDark IEX", text_color="black", text_font_size="8pt")
    header.add_layout(txt)
    txt = Label(x=200, y=0, text="InDark IEX", text_color="black", text_font_size="10pt")
    header.add_layout(txt)
    txt = Label(x=0, y=0.05, text="InDark IEX", text_color="black", text_font_size="10pt")
    header.add_layout(txt)
    txt = Label(x=0, y=0.4, text=species_longform, text_color="black", text_font_size="14pt")
    header.add_layout(txt)
    return header

def spacer(masterplot, tools):
   
    spacer = figure(width=800,height=10, tools=tools,  x_range=masterplot.x_range, y_range=masterplot.y_range)
    spacer.min_border = 0
    spacer.outline_line_color = None
    spacer.xgrid.grid_line_color = None
    spacer.ygrid.grid_line_color = None
    spacer.axis.visible = False
    #background coloring not working
    #spacer.background_fill_color = "beige"
    txt = Label(x=0, y=1, text="InDark IEX", text_color="black", text_font_size="8pt")
    spacer.add_layout(txt)
    return spacer
 
 
def sparkline(df, plots, masterplot, tools, columnID, color):

    # build the line plots
    p1=figure(width=800, height=60, tools=tools,  x_range=masterplot.x_range, y_range=masterplot.y_range)
    p1.line('index', columnID, source=df, line_color=color)
    p1 = plot_attributes(p1, df, columnID)

    plots.append(p1)

    return plots






def make_protein_sparklines(groups, spec, spec_longform, conversion_tbl):

    print "Start make_protein_sparklines"
    print groups
    print spec
    groups = ["fraction", "fraction_order", "experiment", "spec"] + groups
    #groups = ["ENOG410IDXB" , "ENOG410IDXK", "ENOG410IDXN"]
    plots = []
   
   
    groupdf  = feather.read_dataframe("formatted_plants_euNOG_concat.csv.feather")
    #print groupdf.head

    if spec != "All plants":
        groupdf = groupdf[groupdf['spec'] == spec]

        prot_elution_name = "formatted_" + spec + "_proteins_concat.csv.feather"      
        protdf = feather.read_dataframe(prot_elution_name)


    #Remove a column if no observations of the group in a species   
    groupdf = groupdf[groupdf.columns[(groupdf != 0).any()]]
    #print groupdf.head 

    #remove queries that aren't a column name
    groups = [x for x in groups if x in groupdf.columns.tolist()]
    #for group in groups:
    #   if group not in groupdf.columns.tolist():
     #       groups.remove(group)
    print groups 
 
    #Filter the profile down to query groups
    groupdf = groupdf[groups]
    groupdf = groupdf.reset_index(drop=True)

    if spec != "All plants":
        prots_conv_tbl = conversion_tbl[conversion_tbl['ID'].isin(groups)]
        prots_conv_tbl = prots_conv_tbl[prots_conv_tbl['Species']==spec]
        prots = prots_conv_tbl['ProteinID'].tolist()
        prots = [x for x in prots if x in protdf.columns.tolist()]
        prots = ["fraction", "fraction_order", "experiment", "spec"] + prots
        protdf = protdf[prots]
        print protdf.columns

        first_prots = prots_conv_tbl[prots_conv_tbl['ID']==groups[4]]['ProteinID'].tolist()

        first_prots = [x for x in first_prots if x in protdf.columns.tolist()]
    
        
        first_prots = ["fraction", "fraction_order", "experiment", "spec"] + first_prots
    
        print first_prots
        first_prot_df = protdf[first_prots]
    
       # axis_len = len(protdf.index)
    #print groupdf.head
    tools = 'xwheel_zoom,reset, hover'
    palette = ["#0072B2","#E69F00","#009E24","#FF0000", "#979797","#5530AA"] * 100    
    #print palette

    #new_index = pd.Index(range(0,axis_len + 100 ,1), name="index")
    #print(new_index)
    ##protdf.reset_index()
    #protdf  = protdf.reindex(new_index)
    #groupdf.reset_index()
    #groupdf= groupdf.reindex(new_index)
    #print groupdf 
    #The rest of the plots' movement will be tied to the first plot
    #Insert a header
    print groups
    masterplot = master_sparkline(groupdf, tools, groups[4], color="black")
    header = header_line(masterplot,tools, spec_longform)
    plots.append(header)
 
    plots.append(masterplot)
    try:
        for num in range(4, len(first_prots)):
             plots = sparkline(first_prot_df, plots, masterplot, tools, first_prots[num], color = palette[num-5]) 
    except Exception as e:
        print e

    #Add a spacer
    plots.append(spacer(masterplot, tools))

    #now append its component proteins...
    

    for groupnum in range(5, len(groups)):
       if groups[groupnum] not in groupdf.columns.tolist():
           continue

       else: 
           #Insert a label and a break
           plots = sparkline(groupdf ,plots, masterplot, tools, groups[groupnum], color="black")
           #Also append protein sparklines
           if spec != "All plants":
               first_prots = prots_conv_tbl[prots_conv_tbl['ID']==groups[groupnum]]['ProteinID'].tolist()
    
               first_prots = [x for x in first_prots if x in protdf.columns.tolist()]
               first_prots = ["fraction", "fraction_order", "experiment", "spec"] + first_prots
    
               print first_prots
    
               first_prot_df = protdf[first_prots]
               for protnum in range(4, len(first_prots)):
    
                  try:
                      plots = sparkline(first_prot_df, plots, masterplot, tools, first_prots[protnum], color=palette[protnum-5]) 
                  except Exception as e:
                      print e
           plots.append(spacer(masterplot, tools))
    #output_file("lines.html", title="line.py example")
      
    final_plot = gridplot(plots, ncols=1)
    #save(final_plot)
    return components(final_plot)





if __name__ == "__main__":
    make_sparklines("ENOG410IDXB ENOG410IDXK ENOG410IDXN")
