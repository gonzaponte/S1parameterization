from numpy import loadtxt
from ROOT import *
from checkmap import xyz_map

a = loadtxt('log.py')[:,0]
nphmin = ( min( a ) // 100 ) * 100
nphmax = ( max( a ) // 100 ) * 101

xs,ys,zs = zip(*xyz_map.values())
xmin, xmax = min( xs ) - 5, max( xs ) + 5
ymin, ymax = min( ys ) - 5, max( ys ) + 5
zmin, zmax = min( zs ) - 5, max( zs ) + 5

xbin = int( ( xmax - xmin ) // 10 )
ybin = int( ( ymax - ymin ) // 10 )
zbin = int( ( zmax - zmin ) // 10 )

xs, ys, zs = map( list, map( set, [xs,ys,zs] ) )

hxy = { z : TH2F(str(z),'z = ' + str(z) +'; x (mm); y (mm)',xbin,xmin,xmax,ybin,ymin,ymax) for z in zs }
hr = TH2F('r','r dependence;r (mm);# photons',xbin*2,0,xmax,100,nphmin,nphmax)
hz = TH2F('z','z dependence;z (mm);# photons',zbin,zmin,zmax,100,nphmin,nphmax)

for i in range(len(a)):
    x,y,z = xyz_map[i]
    r = ( x**2 + y**2 )**0.5
    bin = hxy[z].Fill(x,y)
    hxy[z].SetBinContent(bin,a[i])
    hr.Fill( r, a[i] )
    hz.Fill( z, a[i] )

canvas = TCanvas()
canvas.Divide(3,2)
for i in range(len(zs)):
    z = zs[i]
    canvas.cd(i+1)
    hxy[z].SetMinimum(nphmin)
    hxy[z].Draw('colz')

canvas.cd(4)
hpr = hr.ProfileX()
hpr.SetLineColor(kRed)
hpr.SetLineWidth(2)
hr.Draw()
hpr.Draw('same')
canvas.cd(5)
hpz = hz.ProfileX()
hpz.SetLineColor(kRed)
hpz.SetLineWidth(2)
hz.Draw()
hpz.Draw('same')
raw_input()
