from numpy import loadtxt
from math import *
from ROOT import *
from PMTmap import pmt_map
from Plots import Plot4D, MakeGif

gStyle.SetOptStat('')

xmin, xmax = -220., 220.
ymin, ymax = -220., 220.
zmin, zmax = -300., 300.
nmin, nmax =    0., 50000/1e6
rmin, rmax =    0., 220.
pmin, pmax =   -pi, pi

xbin = int( (xmax-xmin) // 10 )
ybin = int( (ymax-ymin) // 10 )
zbin = int( (zmax-zmin) // 10 )
nbin = 500
rbin = 100
pbin = 100

xyz_map = { i : ( (i%1936)//44, (i%1936)%44 , -300. + 10*(i//1936) ) for i in range(118096) }

zs = [ -300 + 10*i for i in range(61) ]

S1table = loadtxt( '/Users/Gonzalo/github/S1parameterization/S1table.dat' )

def Make4Dplot( data = S1table ):
    x, y, z = zip( *map( xyz_map.get, zip(*data)[0] ) )
    n = [ sum( dat[1:] ) for dat in data ]
    return Plot4D( x, y, z, n )

def z_dependence( data = S1table ):
    rcor = 40.
    rcor2 = rcor**2
    ncor = int( (xmax/rcor)**2 )
    titles = [ '{0} to {1};z (mm);# photons'.format(r**0.5*rcor,(r+1)**0.5*rcor) for r in range(ncor) ]
    hs = [ TH2F( titles[r], titles[r], zbin, zmin, zmax, nbin, nmin, nmax ) for r in range(ncor) ]

    for dat in data:
        x, y, z = xyz_map[dat[0]]
        n = sum( dat[1:] )
        r2 = ( x**2 + y**2 )
        bin = int( r2 // rcor2 )
        hs[bin].Fill( z, n )

    ps = map( TH2.ProfileX, hs )
    canvas = TCanvas()
    canvas.Divide(7,7)
    for i in range(ncor):
        canvas.cd(i+1)
        hs[i].Draw()
        ps[i].SetLineColor(kRed)
        ps[i].SetLineWidth(2)
        ps[i].Draw('same')

    canvas.Modified()
    canvas.Update()
    return hs, ps, canvas

def r_dependence( data = S1table ):
    titles = { z : 'z = {0};r (mm);# photons'.format(z) for z in zs }
    hs = { z : TH2F( titles[z], titles[z], rbin, rmin, rmax, nbin, nmin, nmax ) for z in zs }
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        n = sum( dat[1:] )
        r = ( x**2 + y**2 )**0.5
        hs[z].Fill(r,n)

    ps = { z : hs[z].ProfileX() for z in zs }
    canvas = TCanvas()
    canvas.Divide(8,8)
    for i,z in enumerate(zs):
        canvas.cd(i+1)
        hs[z].Draw()
        ps[z].SetLineColor(kRed)
        ps[z].SetLineWidth(2)
        ps[z].Draw('same')

    canvas.Modified()
    canvas.Update()
    return hs, ps, canvas

def phi_dependence( data = S1table ):
    titles = { z : 'z = {0};phi (rad);# photons'.format(z) for z in zs }
    hs = { z : TH2F( titles[z], titles[z], pbin, pmin, pmax, nbin, nmin, nmax ) for z in zs }
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        n = sum( dat[1:] )
        phi = atan2( y, x )
        hs[z].Fill(phi,n)

    ps = { z : hs[z].ProfileX() for z in zs }
    canvas = TCanvas()
    canvas.Divide(8,8)
    for i,z in enumerate(zs):
        canvas.cd(i+1)
        hs[z].Draw()
        ps[z].SetLineColor(kRed)
        ps[z].SetLineWidth(2)
        ps[z].Draw('same')

    canvas.Modified()
    canvas.Update()
    return hs, ps, canvas

def xy_dependence( data = S1table, outdir = './forgif/' ):
    gROOT.SetBatch()
    titles = { z : 'z = {0};x (mm);y (mm);# photons'.format(z) for z in zs }
    names = { z : titles[z].split(';')[0] + '.pdf' for z in zs }
    hs = { z : TH2F( titles[z], titles[z], xbin, xmin, xmax, ybin, ymin, ymax ) for z in zs }
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        n = sum( dat[1:] )
        h = hs[z]
        for i in xrange(int(n)):
            h.Fill(x,y)

    for z in zs:
        c = TCanvas()
        hs[z].Draw('zcol')
        c.SaveAs( outdir + titles[z].split(';')[0] + '.png' )
        del c

    MakeGif( sorted(z.keys()), outdir )
