#!/usr/bin/env python3

from argparse import ArgumentParser

import compare_plot as cp
import profile_plot as pp

__all__ = []
__version__ = 0.1
__date__ = '2022-06-03'
__updated__ = '2018-06-03'


def main():
    """Command line options."""

    try:
        # Setup argument parser
        parser = ArgumentParser()

        subparsers = parser.add_subparsers(help="Select profile or scatter plot", dest="command")
        parser_profile = subparsers.add_parser("profile",
                                               help="Generates profile plot/s from a SCRAM profile alignment")

        # profile plot
        parser_profile.add_argument('-a', '--alignment',
                                    type=str,
                                    help="sRNA alignment file prefix used by SCRAM profile (i.e. exclude _21.csv, _22.csv, "
                                         "_24.csv)")
        parser_profile.add_argument('-cutoff', '--cutoff', type=int, default=1,
                                    help="Min. alignment RPMR from the most abundant profile (if multi) to generate plot")

        parser_profile.add_argument('-s', '--search', type=str, help="Full header or substring of header.  Without "
                                                                     "flag, all headers will be plotted", nargs='*',
                                    default="")

        parser_profile.add_argument('-l', '--length', type=str, help="Comma-separated list of sRNA lengths to plot.  "
                                                                     "SCRAM alignment files must be available for each "
                                                                     "sRNA "
                                                                     "length")

        parser_profile.add_argument('-ylim', '--ylim',
                                    type=float, help='+/- y axis limit',
                                    default=0)
        parser_profile.add_argument('-win', '--win', type=int, help='Smoothing window size (default=auto) - '
                                                                    'must be 1 (no smoothing) or an even integer above 5',
                                    default=0)

        parser_profile.add_argument('-pub', '--publish', action='store_true',
                                    default=False,
                                    help='Remove all labels from profiles for editing for publication')
        parser_profile.add_argument('-png', '--png', action='store_true',
                                    default=False,
                                    help='Export plot/s as 300 dpi .png file/s')

        parser_profile.add_argument('-bin_reads', '--bin_reads', action='store_true',
                                    default=False,
                                    help='For plotting large profiles (i.e. chromosomes).  Assigns reads to 10,000 bins prior '
                                         'to smoothing. X-axis shows bin, not reference position')

        # Compare plot
        parser_cdp = subparsers.add_parser("compare", help="Generates a scatter plot for a SCRAM compare alignment")
        parser_cdp.add_argument('-plot_type', '--plot_type', default="log_error",
                                help='Bokeh plot type to display (log, log_error or all)')
        parser_cdp.add_argument('-a', '--alignment',
                                type=str,
                                help="sRNA alignment file prefix used by SCRAM profile (i.e. exclude _21.csv, _22.csv, "
                                     "_24.csv)")
        parser_cdp.add_argument('-l', '--length', type=str, help="Comma-separated list of sRNA lengths to plot.  "
                                                                 "SCRAM alignment files must be available for each "
                                                                 "sRNA length. For an miRNA alignment file, "
                                                                 "use 'mir' instead of an integer")
        parser_cdp.add_argument('-xlab', '--x_label', default=["Treatment 1"],
                                help='x label - corresponds to -s1 treatment in SCRAM arguments', nargs='*')
        parser_cdp.add_argument('-ylab', '--y_label', default=["Treatment 2"],
                                help='y label - corresponds to -s2 treatment in SCRAM arguments', nargs='*')
        parser_cdp.add_argument('-html', '--html', default=False, action='store_true',
                                help='If not using Jupyter Notebook, output interactive plot to browser as save to .html')
        parser_cdp.add_argument('-pub', '--publish', action='store_true',
                                default=False,
                                help='Remove all labels from profiles for editing for publication')
        parser_cdp.add_argument('-png', '--png', action='store_true',
                                default=False,
                                help='Export plot/s as 300 dpi .png file/s')
        parser_cdp.add_argument('-xylim', '--xylim',
                                type=str, help="x and y max. axis limits", default="auto")
        parser_cdp.add_argument('-fig_size', '--fig_size',
                                type=int, help="Output plot dimensions", default=8)

        # Process arguments
        args = parser.parse_args()
        if args.alignment is None:
            print("Missing argument: alignment file prefix required (-a)")
        elif args.length is None:
            print("Missing argument: alignment lengths required (-l)")
        else:
            if args.command == "profile":
                search_term = args.search
                alignment_prefix = args.alignment
                cutoff = args.cutoff
                ylim = args.ylim
                pub = args.publish
                win = args.win
                nt_list = args.length.split(',')
                save_plot = args.png
                bin_reads = args.bin_reads
                pp.profile_plot(nt_list, search_term, alignment_prefix, cutoff, ylim, win, pub, save_plot, bin_reads)

            if args.command == "compare":
                alignment_prefix = args.alignment
                nt_list = args.length.split(',')
                xlab = " ".join(args.x_label)
                ylab = " ".join(args.y_label)
                plot_type = args.plot_type
                browser = args.html
                save_plot = args.png
                pub = args.publish
                fig_size = args.fig_size
                xylim = args.xylim
                cp.compare_plot(alignment_prefix, nt_list, xlab, ylab, plot_type, browser, save_plot, pub, fig_size,
                                xylim)

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0


if __name__ == "__main__":
    main()
