import numpy as np
import sys
from scipy.spatial import ConvexHull
from astropy.coordinates import SkyCoord
import pywcs
import networkx as nx
import pickle
import pyfits
import graphfunc as gf

imagein = 'gal_california_spire350_level3.fits'
hdulist = pyfits.open(imagein)
print hdulist[0].data.shape
data = hdulist[0].data[:,:]
header = hdulist[0].header
hdulist.close()

# image properties
contrms = 1 # 
cell='10arcsec' # cell size for imaging.
cellsize = float(cell[:2]) # arcsec
beamsize = 5.
searchradius = 10000
imsize = [2561,2414] # size of image in pixels. 
conpixels = 5.*beamsize # pixel number in a core

########################
########################
########################

makedirectedG = 0 #
connectedcomponents = 1 #

########################
########################
########################

if makedirectedG == 1:
    rmslevel = '10'
    nodedict,G = gf.twoD_directed_cubegraph_nodedict(data,contrms,float(rmslevel),1.e8,searchradius) 
    nx.write_gpickle(G, 'directedG_cutoff'+rmslevel+'.pickle')
    with open('nodedict_cutoff'+rmslevel+'.pickle', 'wb') as fp:
        pickle.dump(nodedict, fp)

directedgraphs = {
                  '10':{'gfile':'directedG_cutoff10.pickle','dictfile':'nodedict_cutoff10.pickle'},
          }

if connectedcomponents == 1:
    from astropy.io import fits
    directG = nx.read_gpickle(directedgraphs['10']['gfile'])
    G = directG.to_undirected()
    with open(directedgraphs['10']['dictfile'], 'rb') as fp:
        nodemap = pickle.load(fp)
    print 'finding connected components...'
    temp1 = sorted(nx.connected_components(G), key=len, reverse=True) # cores with largest number of pixels rank first, each item in temp has node numbers representing pixels
    print 'finish finding. total ',len(temp1),'now filter by size' 
    temp = [ii for ii in temp1 if len(ii) >= conpixels]
    print 'total ',len(temp),'after filtering' 
    hdu = fits.open(imagein)[0]
    mask = np.zeros(hdu.data.shape)
    xcenter=162.7101
    ycenter=-8.7114
    xmin,xmax = (0,2116)
    ymin,ymax = (849,1745)
    mincolor,maxcolor = (40,350)
    for i in range(3):
        corei = [tuple(G.node[int(j)]['pix']) for j in temp[i]] # i-th core as a list of pixels
        mask[[kk[0] for kk in corei],[kk[1] for kk in corei]] = 1 #
        mask_hdu = fits.PrimaryHDU(mask.astype('short'), hdu.header)
        mask_hdu.writeto('showmasks'+str(i)+'.fits',output_verify='exception',clobber=True,checksum=False)
        mask[:,:] = 0
    pdfname = "calicores.pdf"
    gf.showmasks(['showmasks0.fits','showmasks1.fits','showmasks2.fits'],imagein,pdfname,(xmin,xmax),(ymin,ymax),cellsize,(xcenter,ycenter),(mincolor,maxcolor))


