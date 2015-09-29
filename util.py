from ROOT import *
from math import *


_letters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q',
'R','S','T','U','V','W','X','Y','Z'] + ['a','b','c','d','e','f','g','h','i','j',
'k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

_nletters = len(_letters)
_random = TRandom3(0)

def PositiveAngle(a):
    '''
        Put angle in [0,2pi) range.
    '''
    return 2*pi+a if a<0 else a

def NegativeAngle(a):
    '''
        [0,2pi) to [-pi,pi) conversion.
    '''
    return a if a<pi else a - 2*pi

def atan3( y, x ):
    '''
        Return the phi angle [0,2pi) given y and x.
    '''
    a = atan2(y,x)
    return PositiveAngle(a)

def ParamEvolution( funs, npars, xvals = None ):
    '''
        Creates graphs to see the evolution of the parameters for each
        function.
    '''
    npoints = len(funs)
    if xvals is None:
        xvals = range(npoints)
    pars = [ [ fun.GetParameter(i) for i in range(npars) ] for fun in funs ]
    errs = [ [ fun.GetParError(i) for i in range(npars) ] for fun in funs ]
    graphs = [ TGraphErrors() for i in range(npars) ]
    [ (g.SetMarkerStyle(20),g.SetMarkerSize(1)) for g in graphs ]
    for i in range(npars):
        for j in range(npoints):
            graphs[i].SetPoint(j, xvals[j],pars[j][i])
            graphs[i].SetPointError(j,0,errs[j][i])

    return graphs

def Plot( hs, nh, nv, **options ):
    '''
        Draw hs in a single canvas with nh x nv plots with drawing options
        given by '0' = 'col', '1' = 'AP', etc.
    '''
    c = TCanvas()
    c.Divide(nh,nv)
    for i,h in enumerate(hs):
        c.cd(1+i)
        h.Draw( options.get(str(i),'') )
    return c

def CreatePolynomial( pars, errs = [], name = None ):
    formula = ' + '.join( '[{0}]*x^{0}'.format(i) for i in range(len(pars)) )
    #formula = ' + '.join( '{1}*x^{0}'.format(i,p) for i,p in enumerate(pars) )
    if name is None:
        name = ''.join( [ _letters[_random.Integer(_nletters)] for i in range(50) ] )
    f = TF1( name, formula )
    f.SetParameters(*pars)
    [ f.SetParError(i,err) for i,err in enumerate(errs) ]
    return f
