problemCases = []
p = 179
count = 0
while count < 1:
    R = E.base_ring().ring_of_integers()
    if R != ZZ:
        print R
    if R == ZZ:
        r = tates_algorithm( E , p , R[x] )
        t = r[0][0][0]
        if t.find('In*') != -1:
            if p == 2:
                t = t.replace('n' , '?')
            else:
                t = t.replace('n' , str( E.discriminant().ord(p) - 6 ) ) 
        elif t.find('In') != -1:
            t = t.replace('n', str( E.discriminant().ord(p) ) )
        s = str( E.kodaira_symbol(p) )
        if not t.endswith(s):
            print t , s , E
            problemCases.append(E)
    else:
        for f in R.number_field().factor(p):
            r = tates_algorithm( E , f[0] , R[x] )
            t = r[0][0][0]
            if t.find('In\*') != -1:
                if p == 2:
                    t = t.replace('n' , '?' )
                else:
                    t = t.replace('n' , str( E.discriminant().ord(p) - 6 ) ) 
            elif t.find('In') != -1:
                t = t.replace('n', str( E.discriminant().ord(f[0]) ) )
            s = str( E.kodaira_symbol(f[0]) )
            if not t.endswith(s):
                print t , s , E
                problemCases.append(E)
            print t , s , E
    count += 1
    if (10).divides(count):
        print count
