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
p1.ref_len = 21
p1.all_alignments = [t1, t2]
p2 = pp.SingleRefProfile()
p2.ref_len = 21
p2.all_alignments = [t3, t4]
r = pp.RefProfiles()
r.srna_len = 21
r.replicates = 3
r.single_ref_profiles["L"] = p1
r.single_ref_profiles["S"] = p2
s = pp.DataForPlot(r, "L")


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

    assert b == s


def test_convert_to_coverage():
    t4 = pp.SingleAlignment(pp.DNA("AAA"), 1, "+", 1, np.array([1.0, 2.0, 3.0]))
    t5 = pp.SingleAlignment(pp.DNA("AAT"), 2, "+", 1, np.array([1.0, 2.0, 3.0]))
    p3 = pp.SingleRefProfile()
    p3.ref_len = 5
    p3.all_alignments = [t4, t5]
    r2 = pp.RefProfiles()
    r2.srna_len = 3
    r2.replicates = 3
    r2.single_ref_profiles["L"] = p3
    s2 = pp.DataForPlot(r2, "L")
    s2.convert_to_coverage()
    should_be = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 6.0],
            [2.0, 4.0, 6.0],
            [1.0, 2.0, 3.0],
            [0.0, 0.0, 0.0],
        ]
    )
    assert np.array_equal(s2.fwd, should_be)
