import numpy as np
import networkx as nx

def twoD_directed_cubegraph_nodedict(cube,rmsinput,rmslevel,rmslevelup,searchradius=100): 
    """give 2d cube array, rmslevel, return pixel array and a directed Graph. searchradius in unit of pixel""" 
    xmin,ymin = [0,0]
    xmax,ymax = cube.shape
    nodemap = np.copy(cube)
    nodemap[:,:] = -1 # this stores the node number
    nodedict = {} # this stores the node number
    xcenter,ycenter = xmax/2,ymax/2
    rms = rmsinput
    # get ready for generating graph 
    print 'generating graph, adding nodes and edges' 
    node = 0
    G = nx.DiGraph()
    # generate the graph, add pixel as node if above rmslevel*rms
    for i in range(xmin,xmax):
        if i in range(xmin,xmax+1,(xmax-xmin)/5):
            print int(round(float(i+1)/float(xmax)*100./20.)*20),'%'
        for j in range(ymin,ymax):
            if cube[i,j] >= rmslevel*rms and cube[i,j] < rmslevelup*rms and math.isnan(cube[i,j]) == False and math.sqrt((i-xcenter)**2+(j-ycenter)**2) < searchradius: # ignore nan pixels, also ignore pixels outside searchradius
                G.add_node(node) # node --> pixel [i,j]
                G.node[node]['pix'] = [i,j]
                G.node[node]['intens'] = cube[i,j]
                nodemap[i,j] = node # pixel [i,j] --> node
                nodedict[str(i)+','+str(j)] = node # pixel [i,j] --> node
                if i > xmin:
                    if nodemap[i-1,j] != -1: # go back 1 in i 
                        if cube[i,j] > cube[i-1,j]:
                            G.add_edge(node,nodemap[i-1,j],grad=(cube[i,j]-cube[i-1,j])) 
                        else:
                            G.add_edge(nodemap[i-1,j],node,grad=(-cube[i,j]+cube[i-1,j])) 
                if j > ymin:
                    if nodemap[i,j-1] != -1: # go back 1 in j
                        if cube[i,j] > cube[i,j-1]:
                            G.add_edge(node,nodemap[i,j-1],grad=(cube[i,j]-cube[i,j-1])) 
                        else:
                            G.add_edge(nodemap[i,j-1],node,grad=(-cube[i,j]+cube[i,j-1])) 
                node = node + 1
    return nodedict,G

def showmasks(maskin,imagein,pdfname,(xmin,xmax),(ymin,ymax),cellsize,(xcenter,ycenter),(mincolor,maxcolor)):
    import os
    import aplpy
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    colors = ['red','blue','green','white']
    ypanels=1
    xpanels=1
    wid = (xmax - xmin) * cellsize / 3600.
    hei = (ymax - ymin) * cellsize / 3600.
    fig=plt.figure(figsize=(10.,10./wid*hei))
    ff = aplpy.FITSFigure(imagein,figure=fig)
    ff.recenter(xcenter,ycenter,width=wid,height=hei) 
    ff.set_theme('publication')
    ff.set_system_latex(True)
    ff.set_tick_labels_font(size='x-small')
    ff.set_axis_labels_font(size='small')
    ff.tick_labels.set_yformat('dd.d')
    ff.tick_labels.set_xformat('ddd.d')
    ff.show_colorscale(vmin=mincolor,vmax=maxcolor,cmap='gray',stretch='sqrt',interpolation='none')
    ff.add_colorbar() 
    for c in range(len(maskin)):
        ff.show_contour(maskin[c], levels=1, colors=colors[c], linewidths=0.4)
    fig.canvas.draw()
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    plt.close(fig)
    os.system('open '+pdfname)

