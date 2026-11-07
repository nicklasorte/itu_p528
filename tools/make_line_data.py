"""Regenerate src/p676_line_data_rev1.m from the NTIA reference C++ data files.

Usage:
    python3 make_line_data.py <NTIA p528 repo>/src/p676 ../src/p676_line_data_rev1.m
"""
import re
import sys


def parse(path, cls):
    text = open(path).read()
    out = {}
    pattern = r"const vector<double> " + cls + r"::(\w+)\s*=\s*\{(.*?)\};"
    for m in re.finditer(pattern, text, re.S):
        out[m.group(1)] = [v.strip() for v in m.group(2).replace("\n", " ").split(",") if v.strip()]
    return out


def mrow(vals, per=8):
    lines = ["    " + ", ".join(vals[i:i + per]) for i in range(0, len(vals), per)]
    return " ...\n".join(lines)


def main(src_dir, out_path):
    ox = parse(src_dir + "/OxygenData.cpp", "OxygenData")
    wv = parse(src_dir + "/WaterVapourData.cpp", "WaterVapourData")
    assert all(len(v) == 44 for v in ox.values()), "expected 44 oxygen lines"
    assert all(len(v) == 35 for v in wv.values()), "expected 35 water vapour lines"
    txt = open(out_path).read().split("persistent ox_c wv_c")[0] + "persistent ox_c wv_c\nif isempty(ox_c)\n"
    for k, v in ox.items():
        txt += f"    ox_c.{k} = [ ...\n{mrow(v)}];\n"
    for k, v in wv.items():
        txt += f"    wv_c.{k} = [ ...\n{mrow(v)}];\n"
    txt += "end\nox = ox_c;\nwv = wv_c;\nend\n"
    open(out_path, "w").write(txt)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
