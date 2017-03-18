     
#http://nbviewer.jupyter.org/urls/ep2016.europython.eu/media/conference/slides/networkx-visualization-powered-by-bokeh.ipynb


import networkx



from bokeh.plotting import show, figure
from bokeh.io import output_notebook
from bokeh.models import HoverTool
from bokeh.models import ColumnDataSource

from math import sqrt
network = networkx.read_gml('ep2016.gml')

layout = networkx.spring_layout(network,
                                k=1.1/sqrt(network.number_of_nodes()),
                                iterations=100)


nodes, nodes_coordinates = zip(*sorted(layout.items()))
nodes_xs, nodes_ys = list(zip(*nodes_coordinates))
nodes_source = ColumnDataSource(dict(x=nodes_xs, y=nodes_ys,
                                     name=nodes))



hover = HoverTool(tooltips=[('name', '@name'), ('id', '$index')])
plot = figure(plot_width=800, plot_height=400,
              tools=['tap', hover, 'box_zoom', 'reset'])
r_circles = plot.circle('x', 'y', source=nodes_source, size=10,
                        color='blue', level = 'overlay')

show(plot)





