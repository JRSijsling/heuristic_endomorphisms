CC.<i> = ComplexField(300)

seqs = [ ]
#seq = [ -0.13595418147107261484344693819092277934956316175396 - 0.20244259713450634445443990213922119573791241591512*i, -0.27190836294214522968689387638184555869912632350792 ];
seq = [ 3.1380735450430275094334329000899057043678605681760, 1.5690367725215137547167164500449528521839302840880 + 1.0982040141252551721419720023716696564009347156950*i ];
seqs.append(seq)
seq = [ 1, seq[1]/seq[0] ]
seqs.append(seq)
seq = [ 2*seq[0], 2*seq[1] ]
seqs.append(seq)
seq = [ seq[1], seq[0] ]
seqs.append(seq)

for seq in seqs:
    seqm = magma(seq)
    seqmi = magma([ 1, seqm[2]/seqm[1] ])
    print "seq"
    print seqmi
    RR = magma.RealField(magma.Parent(seqm[1]))
    print "pers"
    print gp.elleisnum(seq, 4, flag=1)
    print 120 * (1/seqm[1])^4 * magma.ZetaFunction(RR, 4) * magma.Eisenstein(4, seqm);
    print 120 * (1/seqm[1])^4 * magma.ZetaFunction(RR, 4) * magma.Eisenstein(4, seqmi);
    print gp.elleisnum(seq, 6, flag=1)
    print 280 * (1/seqm[1])^6 * magma.ZetaFunction(RR, 6) * magma.Eisenstein(6, seqm);
    print 280 * (1/seqm[1])^6 * magma.ZetaFunction(RR, 6) * magma.Eisenstein(6, seqmi);
