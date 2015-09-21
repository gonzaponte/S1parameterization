from numpy import loadtxt
from math import *
from ROOT import *
from PMTmap import pmt_map
from Plots import Plot4D, MakeGif
import scipy.integrate
from util import *

gStyle.SetOptStat('')

xmin, xmax = -220., 220.
ymin, ymax = -220., 220.
zmin, zmax = -300., 300.
nmin, nmax =    0., 30000/1e6
rmin, rmax =    0., 220.
pmin, pmax =   -pi, pi

xbin = int( (xmax-xmin) // 10 )
ybin = int( (ymax-ymin) // 10 )
zbin = int( (zmax-zmin) // 10 )
nbin = 500
rbin = 100
pbin = 100

xyz_map = { i : ( -215 + 10. * ( (i%1936)//44), -215. + 10 * ( (i%1936)%44 ), -300. + 10*(i//1936) ) for i in range(118096) }

zs = [ -300 + 10*i for i in range(61) ]

S1table = loadtxt( '/Users/Gonzalo/github/S1parameterization/S1table.dat' )
data = S1table


def Make4Dplot():
    x, y, z = zip( *map( xyz_map.get, zip(*data)[0] ) )
    n = [ sum( dat[1:] ) for dat in data ]
    return Plot4D( x, y, z, n )

def z_dependence():
    rcor = 40.
    rcor2 = rcor**2
    ncor = int( (xmax/rcor)**2 )
    titles = [ '{0} to {1};z (mm);# photons'.format(r**0.5*rcor,(r+1)**0.5*rcor) for r in range(ncor) ]
    hs = [ TH2F( titles[r], titles[r], zbin, zmin, zmax, nbin, nmin, nmax ) for r in range(ncor) ]
    for h in hs: h.GetXaxis().SetTitleSize(0.05)
    for h in hs: h.GetYaxis().SetTitleSize(0.05)
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        n = sum( dat[1:] )
        r2 = ( x**2 + y**2 )
        bin = int( r2 // rcor2 )
        if bin < ncor:
            hs[bin].Fill( z, n )

    ps = map( TH2.ProfileX, hs )
    canvas = TCanvas()
    canvas.Divide(6,5)
    for i in range(ncor):
        canvas.cd(i+1)
        hs[i].Draw()
        ps[i].SetLineColor(kRed)
        ps[i].SetLineWidth(2)
        ps[i].Draw('same')

    canvas.Modified()
    canvas.Update()

    return hs, ps, canvas

def z_dependence():
    title = '_'#s = [ '{0} to {1};z (mm);# photons'.format(r**0.5*rcor,(r+1)**0.5*rcor) for r in range(ncor) ]
    hs = TH2F( title, title, zbin, zmin, zmax, nbin, nmin, nmax )
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        r = (x**2 + y**2)**0.5
        n = sum( dat[1:] )
        if not n or r > 200.: continue
        hs.Fill( z, n )

    ps = hs.ProfileX()
    canvas = TCanvas()
    hs.Draw()
    ps.SetLineColor(kRed)
    ps.SetLineWidth(2)
    ps.Draw('same')

    return hs, ps, canvas


def r_dependence():
    titles = { z : 'z = {0};r (mm);# photons'.format(z) for z in zs }
    hs = { z : TH2F( titles[z], titles[z], rbin, rmin, rmax, nbin, nmin, nmax ) for z in zs }
    for h in hs.values(): h.GetXaxis().SetTitleSize(0.05)
    for h in hs.values(): h.GetYaxis().SetTitleSize(0.05)

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

def phi_dependence():
    titles = { z : 'z = {0};phi (rad);# photons'.format(z) for z in zs }
    hs = { z : TH2F( titles[z], titles[z], pbin, pmin, pmax, nbin, nmin, nmax ) for z in zs }
    for h in hs.values(): h.GetXaxis().SetTitleSize(0.05)
    for h in hs.values(): h.GetYaxis().SetTitleSize(0.05)

    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if (x**2 + y**2)**0.5 > 215: continue
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


def rphi( pmt = 3, zs = [-200., 0., 200.] ):
    titles = { z : 'z = {0};r (mm);phi (rad);# photons'.format(z) for z in zs }
    #hs     = { z : TH2F( titles[z], titles[z], rbin, rmin, rmax, pbin, pmin, pmax ) for z in zs }
    hs     = { z : TH2F( titles[z], titles[z], 2*rbin, rmin, 2*rmax, 100, 0, 2*pi ) for z in zs }

    x0, y0 = pmt_map[pmt]
    pmt_phi = atan3(y0,x0)
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if not z in zs: continue
        if x**2 + y**2 > 47306.0: continue
        #r   = (x**2 + y**2)**0.5
        r   = ( (x-x0)**2 + (y-y0)**2 )**0.5

        #phi = atan2( y, x )
        phi = PositiveAngle( atan3( y - y0, x - x0 ) - pmt_phi )
        n   = dat[1+pmt]
        bin = hs[z].Fill(r,phi)
        hs[z].SetBinContent( bin, n )

    canvas = TCanvas()
    canvas.Divide(2,2)
    for i,z in enumerate(zs):
        canvas.cd(i+1)
        hs[z].Draw('zcol')

    canvas.Modified()
    canvas.Update()
    return hs, canvas

def rphitrue( pmt = 0, zs = [-230., 0., 230.] ):
    x0, y0 = pmt_map[pmt]
    z0     = -250.#-382.
    xm, xp = x0 - 32., x0 + 32.

    titles = { z : 'true z = {0};r (mm);phi (rad);# photons'.format(z) for z in zs }
    hs     = { z : TH2F( titles[z], titles[z], rbin, rmin, rmax, pbin, pmin, pmax ) for z in zs }
    psf    = { z : lambda y,x: ((z-z0)/2*pi) / ( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )**1.5 for z in zs }

    integrate = scipy.integrate.dblquad
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if not z in zs: continue
        r   = (x**2 + y**2)**0.5
        phi = atan2( y, x )
        n   = integrate(psf[z],xm,xp, lambda x: y0 - sqrt(1024.0-(x-x0)**2), lambda x: y0 + sqrt(1024.0-(x-x0)**2) )[0]
        bin = hs[z].Fill(r,phi)
        hs[z].SetBinContent( bin, n )

    canvas = TCanvas()
    canvas.Divide(2,2)
    for i,z in enumerate(zs):
        canvas.cd(i+1)
        hs[z].Draw('zcol')

    canvas.Modified()
    canvas.Update()
    return hs, canvas


'''
def xy_dependence( outdir = './forgif/' ):
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
'''

def xy_dependence( outdir = './forgif/' ):
    titles = { z : 'z = {0};x (mm);y (mm);# photons'.format(z) for z in zs }
    hs = { z : TProfile2D( titles[z], titles[z], 10, xmin, xmax, 10, ymin, ymax ) for z in zs }
    for h in hs.values(): h.GetXaxis().SetTitleSize(0.05)
    for h in hs.values(): h.GetYaxis().SetTitleSize(0.05)


    for dat in data:
        x, y, z = xyz_map[dat[0]]
        for i,n in enumerate(map(float,dat[1:])):
            x, y = pmt_map[i]
            hs[z].Fill(x,y,n)

    canvas = TCanvas()
    canvas.Divide(8,8)
    for i,z in enumerate(zs):
        canvas.cd(i+1)
        hs[z].SetMinimum( nmin )
        hs[z].SetMaximum( -1111 )
        hs[z].Draw('zcol')

    canvas.Modified()
    canvas.Update()
    return canvas,hs

def zfit( rmax = 200. ):
    hz = TH2F( 'zfit', ';z (mm);# photons', zbin, zmin, zmax, nbin, nmin, nmax )
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        n = sum( dat[1:] )
        r = ( x**2 + y**2 )**0.5
        if r>rmax: continue
        hz.Fill( z, n )

    pz = hz.ProfileX()
    pz.SetLineWidth(2)
    pz.SetLineColor(kYellow)

    canvas = TCanvas()
    hz.Draw()
    pz.Draw('same')
    pz.Fit('pol9','','same',-240,280)

    expofit = pz.GetFunction('pol9')
    return hz, pz, expofit, canvas

def rfit( zmin = -240, zmax = 280 ):
    hr = TH2F( 'rfit', ';r (mm);# photons', rbin, rmin, rmax, nbin, nmin, nmax )
    z_fit = zfit(data,100.)
    raw_input()
    z_fit = z_fit[:-1]
    zfunc = z_fit[-1]
    nphmax = zfunc.Eval(zmin)
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if not (zmin <= z < zmax): continue
        #if not z == zmax-10: continue
        n = sum( dat[1:] ) * nphmax / zfunc.Eval(z)
        r = ( x**2 + y**2 )**0.5
        hr.Fill( r, n )

    pr = hr.ProfileX()
    pr.SetLineWidth(2)
    pr.SetLineColor(kYellow)

    canvas = TCanvas()
    hr.Draw()
    pr.Draw('same')
    #pr.Fit('pol9','','same',-240,280)
    expofit = 0
    #expofit = pz.GetFunction('pol9')
    return hr, pr, expofit, canvas

def r_profiles( zmin = -240, zmax = 280 ):
    N = int( (zmax-zmin)//10 )
    hs = [ TH2F( str(i), '', rbin, rmin, rmax, nbin, nmin, nmax ) for i in range(N) ]
    z_fit = zfit(data,50.)
    raw_input()
    z_fit = z_fit[:-1]
    zfunc = z_fit[-1]
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if not (zmin <= z < zmax ): continue
        n = sum( dat[1:] ) * zfunc.Eval(zmin) / zfunc.Eval(z)
        r = ( x**2 + y**2 )**0.5
        hs[int(z-zmin)//10].Fill(r,n)

    ps = [ h.ProfileX() for h in hs ]
    canvas = TCanvas()
    colors = [ kBlack, kRed, kBlue, kGreen, kOrange, kYellow, kMagenta ] * 10
    for i in range(N):
        ps[i].SetLineColor(colors[i])
        ps[i].SetLineWidth(2)
        ps[i].Draw('same')

    return hs, ps, canvas

def phi_fit( rmax = 200., zmin = -240, zmax = 280 ):
    hphi = TH2F( 'phifit', ';phi;# photons', pbin, pmin, pmax, nbin, nmin, nmax )
    z_fit = zfit(data,100.)
    #raw_input()
    z_fit = z_fit[:-1]
    zfunc = z_fit[-1]
    nphmax = zfunc.Eval(zmin)
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if not (zmin <= z < zmax): continue
        if ( x**2 + y**2 )**0.5 > rmax: continue
        n = sum( dat[1:] ) * nphmax / zfunc.Eval(z)
        phi = atan2(y,x)
        hphi.Fill( phi, n )

    pphi = hphi.ProfileX()
    pphi.SetLineWidth(2)
    pphi.SetLineColor(kYellow)

    canvas = TCanvas()
    hphi.Draw()
    pphi.Draw('same')
    #pr.Fit('pol9','','same',-240,280)
    expofit = 0
    #expofit = pz.GetFunction('pol9')
    return hphi, pphi, expofit, canvas

def xy_fit( zmin = -240., zmax = 280., rmax = 180. ):
    hphi = [ TH2F( 'PMT {0} fit'.format(ID), str(ID) + ';phi;# photons', pbin, pmin, pmax, nbin, nmin/10, nmax/10 ) for ID in sorted(pmt_map) ]
    z_fit = zfit(data,100.)
    z_fit = z_fit[:-1]
    zfunc = z_fit[-1]
    nphmax = zfunc.Eval(zmin)
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if not (zmin <= z < zmax): continue
        if ( x**2 + y**2 )**0.5 > rmax: continue
        for i in pmt_map:
            n = dat[1+i] * nphmax / zfunc.Eval(z)
            phi = atan2(y,x)
            hphi[i].Fill( phi, n )

    colors = [kAzure,kBlack,kGray+1,kBlue,kGreen,kOrange,kViolet,kRed,kYellow,kMagenta,kGray,kCyan]
    pphi = map(TH2.ProfileX,hphi)
    for ID in pmt_map:
        canvas = TCanvas()
        pphi[ID].SetLineWidth(2)
        pphi[ID].SetLineColor(colors[ID])
        pphi[ID].SetMinimum(0)
        hphi[ID].Draw()
        pphi[ID].Draw('same')
        raw_input()
        del canvas

    canvas = TCanvas()
    for p in pphi: p.Draw('same')
    canvas.BuildLegend()
    expofit = 0
    #expofit = pz.GetFunction('pol9')
    return hphi, pphi, expofit, canvas

def r3( pmt = 0, zmin = -200 ):
    h = TH2F( 'h', '', 1000, 0, 1000, 1000, nmin, nmax*0.12 )
    x0, y0 = pmt_map[pmt]
    z0     = -382.
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if z<zmin: continue
        r = ( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )**0.5
        n = dat[1+pmt]
        h.Fill( r, n )

    h.Draw()
    raw_input()

def rphi_z( pmt = 0, z = -200. ):
    global nmax
    nmax  *= 0.1
    z0     = z
    title  = 'z = {0};r (mm);phi (rad);# photons'.format(z)
    hrphi  = TH2F( title, title, rbin, rmin, rmax, pbin, pmin, pmax )
    titles = [ 'cor = {0};phi (rad);# photons'.format(i) for i in range(6) ]
    hr     = [ TH2F( titles[i], titles[i], pbin/4, pmin, pmax, nbin, nmin, nmax ) for i in range(6) ]
    titles = [ 'phi = ({0}#pm 1) pi/6;r (mm);# photons'.format(-5 + 2*i) for i in range(6) ]
    hphi   = [ TH2F( titles[i], titles[i], rbin, rmin, rmax, nbin, nmin, nmax ) for i in range(6) ]

    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if z != z0: continue
        r2  = x**2 + y**2
        r   = r2**0.5
        phi = atan2( y, x )
        n   = dat[1+pmt]
        hrphi.SetBinContent( hrphi.Fill(r,phi), n )
        if r2 < 15000: hr[int(r2)//2500].Fill(phi,n)
        hphi[int((phi+pi)//(0.4*pi))].Fill(r,n)


    fs = range(len(hr))
    for i in fs:
        f = TF1('f' + str(i),'[0]+[1]*sin([2]*(x-[3]))')
        f.SetParameters(0.002,0.001,1)
        f.SetParLimits(0,0,0.01)
        f.SetParLimits(1,0,0.01)
        f.SetParLimits(2,0,1.5)
        f.SetParLimits(3,-1,1)
        fs[i] = f

    graphs = [ TGraph() for j in range(4) ]

    hrp = map( TH2.ProfileX, hr )
    [h.SetMinimum(0.) for h in hrp]
    hrphi.Draw('zcol')
    canvas1 = TCanvas()
    canvas1.Divide(3,2)
    gcanvas = TCanvas()
    gcanvas.Divide(2,2)
    for i,h in enumerate(hrp):
        canvas1.cd(1+i)
        h.Draw()
        h.Fit(fs[i])
        for j in range(4):
            graphs[j].SetPoint(i,i,fs[i].GetParameter(j))

    canvas1.Modified()
    canvas1.Update()

    names = [ 'Amplitude offset', 'Amplitude', 'Frequency', 'Frequency offset' ]
    for j in range(4):
        gcanvas.cd(1+j)
        graphs[j].SetTitle(names[j])
        graphs[j].SetMarkerStyle(20)
        graphs[j].SetMarkerSize(1)
        graphs[j].Draw('ap')

    gcanvas.Modified()
    gcanvas.Update()


    raw_input('done')

    canvas2 = TCanvas()
    canvas2.Divide(3,2)
    for i,h in enumerate(hphi):
        canvas2.cd(1+i)
        h.Draw()
    canvas2.Modified()
    canvas2.Update()
    raw_input('done')
    return hs, canvas

def r_parameters( pmt = 0, z0 = -200., nsectors = 9 ):
    titles = [ 'phi = {0} -> {1} deg;r (mm);# photons'.format(i*360./nsectors,(i+1)*360./nsectors) for i in range(nsectors) ]
    hr     = [ TH2F( titles[i], titles[i], rbin, rmin, rmax, nbin, nmin, nmax/10 ) for i in range(nsectors) ]

    phi_sector = 2*pi/nsectors

    x0, y0 = pmt_map[pmt]
    pmt_phi = atan3(x0,y0)
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if z != z0: continue
        r   = ( (x-x0)**2 + (y-y0)**2 )**0.5
        phi = atan3( y, x ) - pmt_phi
        phi_bin = int( phi // phi_sector )
        n   = dat[1+pmt]
        hr[ phi_bin ].Fill( r, n )

    plot = Plot( hr, 3, 3 )
    raw_input()

    fs = range(len(hr))
    for i in fs:
        f = TF1('f' + str(i),'[0]+[1]*sin([2]*(x-[3]))')
        f.SetParameters(0.002,0.001,1)
        f.SetParLimits(0,0,0.01)
        f.SetParLimits(1,0,0.01)
        f.SetParLimits(2,0,1.5)
        f.SetParLimits(3,-1,1)
        fs[i] = f

    graphs = [ TGraph() for j in range(4) ]

    hrp = map( TH2.ProfileX, hr )
    [h.SetMinimum(0.) for h in hrp]
    hrphi.Draw('zcol')
    canvas1 = TCanvas()
    canvas1.Divide(3,2)
    gcanvas = TCanvas()
    gcanvas.Divide(2,2)
    for i,h in enumerate(hrp):
        canvas1.cd(1+i)
        h.Draw()
        h.Fit(fs[i])
        for j in range(4):
            graphs[j].SetPoint(i,i,fs[i].GetParameter(j))

    canvas1.Modified()
    canvas1.Update()

    names = [ 'Amplitude offset', 'Amplitude', 'Frequency', 'Frequency offset' ]
    for j in range(4):
        gcanvas.cd(1+j)
        graphs[j].SetTitle(names[j])
        graphs[j].SetMarkerStyle(20)
        graphs[j].SetMarkerSize(1)
        graphs[j].Draw('ap')

    gcanvas.Modified()
    gcanvas.Update()


    raw_input('done')

    canvas2 = TCanvas()
    canvas2.Divide(3,2)
    for i,h in enumerate(hphi):
        canvas2.cd(1+i)
        h.Draw()
    canvas2.Modified()
    canvas2.Update()
    raw_input('done')
    return hs, canvas

def rrel_parameters( z0 = -200., nsectors = 9 ):
    titles = [ 'phi = {0} -> {1} deg;r (mm);# photons'.format(i*360./nsectors,(i+1)*360./nsectors) for i in range(nsectors) ]
    hinner = [ TGraph() for i in range(nsectors) ]
    houter = [ TGraph() for i in range(nsectors) ]

    phi_sector = 2*pi/nsectors
    pinner = [ 0 ] * nsectors
    pouter = [ 0 ] * nsectors
    for dat in data:
        x, y, z = xyz_map[dat[0]]
        if z != z0: continue
        if x**2 + y**2 > 47306.0: continue
        for pmt in range(12):
            x0, y0 = pmt_map[pmt]
            pmt_phi = atan3(y0,x0)
            r   = ( (x-x0)**2 + (y-y0)**2 )**0.5
            phi = PositiveAngle( atan3( y-y0, x-x0 ) - pmt_phi )
            phi_bin = int( phi // phi_sector )
            n = dat[1+pmt]
            if not n: continue
            h = hinner if pmt < 3 else houter
            p = pinner if pmt < 3 else pouter
            h[ phi_bin ].SetPoint( p[phi_bin], r, n )
            p[phi_bin] += 1

    [ (g.SetMinimum(0), g.SetMarkerStyle(20), g.SetMarkerSize(0.5), g.SetTitle(titles[i]) ) for i,g in enumerate(hinner) ]
    [ (g.SetMinimum(0), g.SetMarkerStyle(20), g.SetMarkerSize(0.5), g.SetTitle(titles[i]) ) for i,g in enumerate(houter) ]

    plotinner = Plot( hinner, 3, 3, **{ str(i) : 'AP' for i in range(nsectors) } )
    plotouter = Plot( houter, 3, 3, **{ str(i) : 'AP' for i in range(nsectors) } )
    #raw_input()

    parinner = []
    errinner = []
    for h in hinner:
        h.Fit('pol2')
        f = h.GetFunction('pol2')
        parinner.append( [ f.GetParameter(i) for i in range(3) ] )
        errinner.append( [ f.GetParError(i) for i in range(3) ] )

    parouter = []
    errouter = []
    for h in houter:
        h.Fit('pol2')
        f = h.GetFunction('pol2')
        parouter.append( [ f.GetParameter(i) for i in range(3) ] )
        errouter.append( [ f.GetParError(i) for i in range(3) ] )

    finner = map( CreatePolynomial, parinner, errinner )
    fouter = map( CreatePolynomial, parouter, errouter )

    innerevol = ParamEvolution( finner, 3 )
    outerevol = ParamEvolution( fouter, 3 )

    plotinnerparam = Plot( innerevol, 2, 2, **{ str(i) : 'AP' for i in range(3) } )
    plotouterparam = Plot( outerevol, 2, 2, **{ str(i) : 'AP' for i in range(3) } )

    return finner, fouter

def rrel_param_z( nsectors = 9, zmin = -200., zmax = -100. ):
    z = zmin

    finner = []
    fouter = []
    while z <= zmax:
        fi, fo = rrel_parameters( z, nsectors )
        finner.append(fi[4])
        fouter.append(fo[4])
        z += 10.


    innerevol = ParamEvolution( finner, 3 )
    outerevol = ParamEvolution( fouter, 3 )
    plotinnerparam = Plot( innerevol, 2, 2, **{ str(i) : 'AP' for i in range(3) } )
    plotouterparam = Plot( outerevol, 2, 2, **{ str(i) : 'AP' for i in range(3) } )

    raw_input()


a = Make4Dplot
b = z_dependence
c = r_dependence
d = phi_dependence
e = xy_dependence
f = zfit
g = rfit
h = phi_fit
i = xy_fit
j = rphi
k = rphitrue
l = r3
m = rphi_z
n = r_parameters
o = rrel_parameters
p = rrel_param_z()
#x = b()
