
ifilename = 'missing.txt'
ofilename = 'missing.dat'

ifile = open( ifilename )
ofile = open( ofilename, 'w' )

ostr = '#globID relID zindex zvalue\n'

zlimit = 200.
n = 0
for line in ifile:
    globalID    = int( line.split('_')[1].split('.')[0] )
    relativeID  = globalID % 1936
    zindex      = globalID // 1936
    zvalue      = -300. + zindex * 10.
    if zvalue > zlimit: continue
    n += 1
    ostr       += ' '.join( map(str,[globalID,relativeID,zindex,zvalue]) ) + '\n'

print n

ofile.write(ostr)
ofile.close()
    
