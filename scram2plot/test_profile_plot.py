import scram2plot.profile_plot as pp
# import profile_plot as pp

def test_file_load():
    t1=pp.SingleAlignment(pp.DNA("TCGTTGTTACCTGATTGCTCT"),1,"-",1,[1.0,2.0,3.0])
    t2=pp.SingleAlignment(pp.DNA("AGAGCAATCAGGTAACAACGA"),1,"+",1,[2.0,4.0,6.0])
    t3=pp.SingleAlignment(pp.DNA("TCGTTGTTACCTGATTGCTCT"),1,"-",1,[1.0,2.0,3.0])
    t4=pp.SingleAlignment(pp.DNA("AGAGCAATCAGGTAACAACGA"),1,"+",1,[2.0,4.0,6.0])
    p1=pp.SingleRefProfile()
    p1.ref_len=3
    p1.all_alignments=[t1,t2]
    p2=pp.SingleRefProfile()
    p2.ref_len=3
    p2.all_alignments=[t3,t4]
    r = pp.RefProfiles("test")
    r.srna_len = 21
    r.single_ref_profiles["L"]=p1
    r.single_ref_profiles["S"]=p2
    a = pp.RefProfiles("test")

    a.load_single_ref_profiles("scram2plot/test_data/profile_test_data21.csv")
    assert a==r

def test_extract_from_profiles():
    a = pp.RefProfiles("test")
    a.load_single_ref_profiles("scram2plot/test_data/profile_test_data21.csv")
    b = pp.DataForPlot()
    b.extract_from_ref_profiles(a, 'L')
    assert b.fwd == [0,[2.0, 4.0, 6.0],0,0]