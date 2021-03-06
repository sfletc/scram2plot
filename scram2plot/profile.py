import numpy as np
import matplotlib.pyplot as plt  # @Reimport
import os.path
import csv


class DNA(object):
    """
    DNA class
    """

    dna_alphabet = set("AGCTN")

    def __init__(self, sequence):
        self.sequence = sequence.upper()

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __hash__(self):
        return hash(self.sequence)

    def __repr__(self):
        return self.sequence

    def __eq__(self, other):
        return self.sequence == other.sequence


class SingleAlignment(object):
    """
    Single sRNA read alignment class
    """

    def __init__(self, srna, position, strand, times_aligned, indv_alignments):
        self.srna = srna
        self.position = position
        self.strand = strand
        self.times_aligned = times_aligned
        self.indv_alignments = indv_alignments

    def srna_len(self):
        return len(self.srna)

    def standard_error(self):  # TODO: check - likely don't neeed here
        return np.std(self.indv_alignments, ddof=1) / np.sqrt(
            np.size(self.indv_alignments)
        )

    def mean_alignments(self):
        return np.mean(self.indv_alignments)

    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}".format(
            self.srna,
            self.position,
            self.strand,
            self.times_aligned,
            self.indv_alignments,
        )

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (
            self.srna == other.srna
            and self.position == other.position
            and self.strand == other.strand
            and self.times_aligned == other.times_aligned
            and np.array_equal(self.indv_alignments, other.indv_alignments)
        )


class SingleRefProfile(object):
    """
    Single reference sequence class
    """

    def __init__(self):
        self.ref_len = 0
        self.all_alignments = []

    def __str__(self):
        return "{0}\t{1}".format(self.ref_len, self.all_alignments)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (
            self.ref_len == other.ref_len
            and self.all_alignments == other.all_alignments
        )


class RefProfiles(object):
    """
    All references in a file class
    """

    def __init__(self):
        self.srna_len = 0
        self.replicates = 0
        self.single_ref_profiles = {}

    def __str__(self):
        return "{0}\t{1}\t{2}".format(
            self.srna_len, self.replicates, self.single_ref_profiles
        )

    def __eq__(self, other):
        return (
            self.replicates == other.replicates
            and self.single_ref_profiles == other.single_ref_profiles
            and self.srna_len == other.srna_len
        )

    def load_single_ref_profiles(self, in_file):
        """Loads a scram2 profile file for a single sRNA length

        Args:
            in_file (string): scram2 profile file path
        """
        with open(in_file, "r") as in_handle:
            reader = csv.reader(in_handle, delimiter=",")
            for row in reader:
                if row[0] == "Header":
                    continue
                header = row[0]
                ref_len = int(row[1])
                srna = DNA(row[2])
                position = int(row[3])
                strand = row[4]
                times_aligned = int(row[5])
                indv_alignments = np.array([float(x) for x in row[6:]])
                sa = SingleAlignment(
                    srna, position, strand, times_aligned, indv_alignments
                )
                if header not in self.single_ref_profiles:
                    self.single_ref_profiles[header] = SingleRefProfile()
                    self.single_ref_profiles[header].ref_len = ref_len
                    self.srna_len = len(srna)
                self.single_ref_profiles[header].all_alignments.append(sa)
            self.replicates = len(sa.indv_alignments)


class DataForPlot(object):
    """
    Data for plotting class
    """

    def __init__(self, ref_profiles, header):
        self.ref_profiles = ref_profiles
        self.header = header
        self.ref_len = ref_profiles.single_ref_profiles[header].ref_len
        self.srna_len = ref_profiles.srna_len
        self.replicates = ref_profiles.replicates
        self.fwd = np.zeros((self.ref_len + 1, self.replicates), dtype=float)
        self.rvs = np.zeros((self.ref_len + 1, self.replicates), dtype=float)
        self.x_axis = list(range(len(self.fwd)))
        self.y_flat = []
        self._extract_from_ref_profiles()

    def _extract_from_ref_profiles(self):
        """Extracts data from a ref profiles object

        Args:
            ref_profiles (RefProfiles): ref profiles object
            header (string): header of the reference sequence
        """
        for sa in self.ref_profiles.single_ref_profiles[self.header].all_alignments:
            if sa.strand == "+":
                self.fwd[sa.position] = sa.indv_alignments
            else:
                self.rvs[sa.position] = sa.indv_alignments

    def convert_to_coverage(self, abund=True):
        """_summary_

        Args:
            abund (bool, optional): _description_. Defaults to True.
        """
        self.fwd = self._coverage_per_strand(self.fwd, abund)
        self.rvs = self._coverage_per_strand(self.rvs, abund)

    def _coverage_per_strand(
        self, old_array, abund
    ): 
        """_summary_

        Args:
            old_array (_type_): _description_
            abund (_type_): _description_

        Returns:
            _type_: _description_
        """
        new_arr = np.zeros((old_array.shape), dtype=float)
        for i in range(len(new_arr) - self.srna_len + 1):
            for j in range(len(new_arr[i])):
                if not abund:
                    new_arr[i : i + self.srna_len, j] += old_array[i, j]
                else:
                    for k in range(self.srna_len):
                        if old_array[i, j] > new_arr[i + k, j]:
                            new_arr[i + k, j] = old_array[i, j]
        return new_arr

    def convert_to_error_bounds(self):  # TODO: document and test
        """_summary_
        """
        self.fwd = self._error_bounds(self.fwd)
        self.rvs = self._error_bounds(self.rvs)

    def _error_bounds(self, old_array):  # TODO: document and test
        """_summary_

        Args:
            old_array (_type_): _description_

        Returns:
            _type_: _description_
        """
        new_arr = np.zeros((self.ref_len + 1, 2), dtype=float)
        for i in range(len(new_arr)):
            new_arr[i, 0] = np.mean(old_array[i, :]) - (
                np.std(old_array[i, :]) / np.sqrt(self.replicates)
            )
            new_arr[i, 1] = np.mean(old_array[i, :]) + (
                np.std(old_array[i, :]) / np.sqrt(self.replicates)
            )
        return new_arr

    def flatten(self, d=1):
        if d == 1:
            for i in range(self.fwd.shape[1]):
                self.y_flat.append(self.fwd[:, [i]].flatten())
                self.y_flat.append(-self.rvs[:, [i]].flatten())
        else:
            for i in range(self.fwd.shape[1]):
                self.y_flat.append(self._smoothTriangle(self.fwd[:, [i]].flatten(), d))
                self.y_flat.append(self._smoothTriangle(-self.rvs[:, [i]].flatten(), d))

    @staticmethod
    def _smoothTriangle(data, degree):
        triangle = np.concatenate(
            (np.arange(degree + 1), np.arange(degree)[::-1])
        )  # up then down
        smoothed = []

        for i in range(degree, len(data) - degree * 2):
            point = data[i : i + len(triangle)] * triangle
            smoothed.append(np.sum(point) / np.sum(triangle))
        # Handle boundaries
        smoothed = [smoothed[0]] * int(
            degree + degree / 2
        ) + smoothed  # TODO: this can be better
        while len(smoothed) < len(data):
            smoothed.append(smoothed[-1])
        return smoothed

    def __str__(self):  # TODO: update
        return "{0}\t{1}\t{2}\t{3}\t{4}".format(
            self.header,
            self.ref_len,
            self.srna_len,
            self.fwd,
            self.rvs,
        )

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):  # TODO: update
        return (
            self.header == other.header
            and self.ref_len == other.ref_len
            and self.srna_len == other.srna_len
            and np.array_equal(self.fwd, other.fwd)
            and np.array_equal(self.rvs, other.rvs)
        )


def profile_plot(
    align_prefix, align_lens, header, smooth=1, cov=True, abund=True, se=True
):
    """Wrapper profile plot function

    Args:
        align_prefix (str): Output prefix from scram2 aligner
        align_lens (list[int]): list of alignment length integers  - e.g. [21,22,24]
        header (str): Reference sequence header to plot 
        smoothing_window (int, optional): Plot smoothing value. Defaults to 1.
        cov (bool, optional): Display coverage if true, otherwise show show count at 5' read position. Defaults to True.
        abund (bool, optional): Display count of the most abundant read at each position if true, otherwise standard coverage. Defaults to True.
        se (bool, optional): Show standard error bounds if true, otherwise show individual akignment files. Defaults to True.

    Raises:
        FileNotFoundError: if an alignment file is missing
    """
    file_paths = []
    set_up_plot()
    for i in align_lens:
        file_paths.append("{0}_{1}.csv".format(align_prefix, i))
    for i in file_paths:
        if not os.path.isfile(i):
            raise FileNotFoundError("File {0} not found".format(i))
        else:

            rp = RefProfiles()
            rp.load_single_ref_profiles(i)
            single_plot(rp, header, smooth, cov, abund, se)
    if se:
        plt.legend()
    plt.show()


def set_up_plot():
    """_summary_
    """
    plt.figure(figsize=(24, 6), dpi=80)
    plt.xlabel("Position")
    plt.ylabel("Abundance")
    plt.title("Abundance Profile")


def single_plot(ref_profiles, header, smoothing_window, cov, abund, se):
    """_summary_

    Args:
        ref_profiles (_type_): _description_
        header (_type_): _description_
        smoothing_window (_type_): _description_
        cov (_type_): _description_
        abund (_type_): _description_
        se (_type_): _description_
    """
    cols = {
        24: "darkgreen",
        21: "red",
        22: "blue",
        18: "#f781bf",
        19: "#a65628",
        20: "#984ea3",
        23: "#999999",
        25: "#e41a1c",
        26: "#dede00",
    }  

    spd = DataForPlot(ref_profiles, header)
    plt.plot(spd.x_axis, [0] * len(spd.x_axis), color="grey", linewidth=0.5)
    if cov:
        if abund:
            spd.convert_to_coverage(abund=True)
        else:
            spd.convert_to_coverage(abund=False)
    if se:
        spd.convert_to_error_bounds()
        spd.flatten(smoothing_window)
        plt.fill_between(
            spd.x_axis,
            spd.y_flat[0],
            spd.y_flat[2],
            color=cols[spd.srna_len],
            alpha=0.4,
            label=str(spd.srna_len) + " nt",
        )
        plt.fill_between(
            spd.x_axis,
            spd.y_flat[1],
            spd.y_flat[3],
            color=cols[spd.srna_len],
            alpha=0.4,
        )
    else:
        spd.flatten(smoothing_window)
        first = True
        for i in range(spd.replicates):
            if first:
                plt.plot(
                    spd.x_axis,
                    spd.y_flat[i],
                    spd.y_flat[i + spd.replicates],
                    color=cols[spd.srna_len],
                    alpha=0.8,
                    label=str(spd.srna_len) + " nt",
                )
                first = False
            else:
                plt.plot(
                    spd.x_axis,
                    spd.y_flat[i],
                    spd.y_flat[i + spd.replicates],
                    color=cols[spd.srna_len],
                    alpha=0.8,
                )
            plt.fill_between(
                spd.x_axis,
                spd.y_flat[i],
                spd.y_flat[i + spd.replicates],
                color=cols[spd.srna_len],
                alpha=0.05,
            )
