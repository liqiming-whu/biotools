#!/usr/bin/env python3
import sys
import argparse
import pandas as pd


def df_merge_row(df, by=None, columns=None, method="distinct", sep=",", count=False):
    if not columns:
        return df.drop_duplicates(subset=by, keep='first')

    def generate_content(x):
        data = {}
        for c, m in zip(columns, method):
            try:
                series = pd.to_numeric(x[c])
                if m == "mean":
                    data[c] = series.mean()
                elif m == 'medium':
                    data[c] = series.median()
            except ValueError:
                series = x[c]
            if m == 'distinct':
                data[c] = sep.join(map(str, series.unique()))
            elif m == 'collapse':
                data[c] = sep.join(map(str, series))
            elif m == 'max':
                data[c] = series.max()
            elif m == 'min':
                data[c] = series.min()
            elif m == 'sum':
                data[c] = series.sum()
            elif m == 'first':
                data[c] = series.to_list()[0]
            elif m == 'last':
                data[c] = series.to_list()[-1]
            elif m == 'count':
                data[c] = series.count()
            elif m == 'count_distinct':
                data[c] = pd.Series(series.unique()).count()
            else:
                if c not in data:
                    data[c] = sep.join(map(str, series.unique()))
        if count:
            data["count"] = len(x[columns[0]])
        return pd.Series(data)

    merge_df = df.groupby(by).apply(generate_content).reset_index()
    return merge_df


def main(inpath, outpath, header, by, columns=None, method="distinct", sep="\t", delim=",", count=False):
    if inpath.endswith(".xlsx"):
        df = pd.read_excel(inpath, header=header)
    else:
        df = pd.read_csv(inpath, header=header, delimiter=sep)
    by = by.split(",")
    method = method.split(",")
    if columns:
        by_name = [df.columns[int(i)] for i in by]
        columns = columns.split(",")
        if len(method) == 1:
            method = len(columns) * method
        assert len(columns) == len(method)
        try:
            column_name = [df.columns[int(i)] for i in columns if i not in by]
            output_name = [df.columns[int(i)] for i in sorted(by + columns)]
            if count:
                output_name.append("count")
        except Exception:
            column_name = None
            output_name = by_name
    else:
        by_name = [df.columns[int(i)] for i in by]
        column_name = [i for i in df.columns if i not in by_name]
        if len(method) == 1:
            method = len(column_name) * method
        assert len(column_name) == len(method)
        output_name = df.columns.to_list()
        if count:
            output_name.append("count")

    out_df = df_merge_row(df, by_name, column_name, method, delim, count)
    out_df = out_df[output_name]

    out_header = True if header is not None else False
    if outpath:
        if outpath.endswith(".csv"):
            out_df.to_csv(outpath, index=False, header=out_header, sep=",")
        elif outpath.endswith(".xlsx"):
            out_df.to_excel(outpath, index=False, header=out_header)
        else:
            out_df.to_csv(outpath, index=False, header=out_header, sep="\t")
    else:
        if out_header:
            print("\t".join(map(str, out_df.columns)))
        for _, row in out_df.iterrows():
            print("\t".join(map(str, dict(row).values())))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", dest="inpath",
                        type=str, required=True, help="Input table file. [required].")
    parser.add_argument("-o", dest="outpath", type=str,
                        help="Output file, discard this parameter to print output to stdout.")
    parser.add_argument("-header", dest="header", type=int,
                        default=None, help=("The selected row will be the header, "
                                            "discard this parameter if there is no header, 0 based."))
    parser.add_argument("-by", dest="by", type=str, required=True,
                        help="The same values in the selected columns will be merged, 0 based, use ',' to separate. [required].")
    parser.add_argument("-c", dest="columns", type=str, default=None,
                        help=("Other columns selected will be shown in the output, "
                              "discard this parameter to show all columns, 0 based, use ',' to separate. "
                              "If an illegal parameter is set, no other columns will be appended to the output table."))
    parser.add_argument("-m", dest="method", type=str, default="distinct",
                        help=("Specify the operation that should be applied to -c. "
                              "Valid operations: sum, min, max, mean, median, collapse, distinct, first, last, count, count_distinct. "
                              "use ',' to separate."))
    parser.add_argument("-sep", dest="sep", type=str, default="\t",
                        help="Specify a custom delimiter to separate the input file.")
    parser.add_argument("-delim", dest="delim", type=str, default=",",
                        help="Specify a custom delimiter for the collapse operations.")
    parser.add_argument("-count", dest="count", action="store_true",
                        help="Add a count column to the output table.")
    args = parser.parse_args()

    main(**args.__dict__)
