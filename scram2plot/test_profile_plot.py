import scram2plot.profile_plot as pp

# import profile_plot as pp
import numpy as np

t1 = pp.SingleAlignment(
    pp.DNA("TCGTTGTTACCTGATTGCTCT"), 1, "-", 1, np.array([1.0, 2.0, 3.0])
)
t2 = pp.SingleAlignment(
    pp.DNA("AGAGCAATCAGGTAACAACGA"), 1, "+", 1, np.array([2.0, 4.0, 6.0])
)
t3 = pp.SingleAlignment(
    pp.DNA("TCGTTGTTACCTGATTGCTCT"), 1, "-", 1, np.array([1.0, 2.0, 3.0])
)
t4 = pp.SingleAlignment(
    pp.DNA("AGAGCAATCAGGTAACAACGA"), 1, "+", 1, np.array([2.0, 4.0, 6.0])
)
p1 = pp.SingleRefProfile()
p1.ref_len = 3
p1.all_alignments = [t1, t2]
p2 = pp.SingleRefProfile()
p2.ref_len = 3
p2.all_alignments = [t3, t4]
r = pp.RefProfiles()
r.srna_len = 21
r.replicates = 3
r.single_ref_profiles["L"] = p1
r.single_ref_profiles["S"] = p2
s = pp.DataForPlot(r, "L")
s.extract_from_ref_profiles()


def test_file_load():
    a = pp.RefProfiles()
    a.load_single_ref_profiles("scram2plot/test_data/profile_test_data21.csv")
    assert a == r


def test_mean():
    assert t1.mean_alignments() == 2.0


def test_extract_from_profiles():
    a = pp.RefProfiles()
    a.load_single_ref_profiles("scram2plot/test_data/profile_test_data21.csv")
    b = pp.DataForPlot(a, "L")
    b.extract_from_ref_profiles()
    assert b == s

