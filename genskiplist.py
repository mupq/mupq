#!/usr/bin/env python3
"""This script is used to generate the skiplist.py files in pqm{4,3}. It
takes the output files from the mps2-an386 target, when run via the make
script, i.e., compile pqm{4,3} with the mps2-an38{6,5} plattform and the
run-stack-tests targets (remember to backup your benchmark results, if
you want to keep them!):

    cp -r benchmarks benchmarks.bak
    make clean # Because this will delete it
    make PLATFORM=mps2-an386 -j4 run-stack-tests
    find benchmarks -name frommake -exec python3 mupq/genskiplist.py {} + > skiplist.py

You can also simply update a single scheme:

    make PLATFORM=mps2-an386 benchmarks/stack/crypto_kem/mysuperscheme/opt/frommake
    python3 mupq/genskiplist.py benchmarks/stack/crypto_kem/mysuperscheme/opt/frommake
"""

import argparse
import pprint
import re
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate a skiplist.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "inputs",
        help="Stack benchmark output of mps2 platform",
        nargs="+",
        type=argparse.FileType("r"),
    )
    parser.add_argument(
        "-r",
        "--round",
        help="Round up the usage to N bytes",
        default=1024,
        type=int,
        metavar="N",
    )
    parser.add_argument(
        "-m",
        "--margin",
        help="Add N bytes of additional stack margin",
        default=32,
        type=int,
        metavar="M",
    )
    return parser.parse_args()


def parse_flashsize(contents):
    match = re.search(
        r"^\s+(?P<text>\d+)"
        r"\s+(?P<data>\d+)"
        r"\s+(?P<bss>\d+)"
        r"\s+(?P<dec>\d+)"
        r"\s+(?P<hex>[0-9A-Fa-f]+)"
        r"\s+(?P<filename>[a-zA-Z0-9_\-/\.]+)",
        contents,
        re.MULTILINE,
    )
    if match is None:
        raise Exception("Size output not found!")
    text = int(match.group("text"))
    data = int(match.group("data"))
    bss = int(match.group("bss"))
    return match.group("filename"), text + data, data + bss


def parse_stackusage(contents):
    match = re.search(
        r"^keypair stack usage:\s+(?P<usage>\d+)",
        contents,
        re.MULTILINE,
    )
    if match is None:
        raise Exception("keypair usage not found")
    keypair = int(match.group("usage"))
    match = re.search(
        r"^(sign|decaps) stack usage:\s+(?P<usage>\d+)",
        contents,
        re.MULTILINE,
    )
    if match is None:
        raise Exception("private-op usage not found")
    private = int(match.group("usage"))
    match = re.search(
        r"^(verify|decaps) stack usage:\s+(?P<usage>\d+)",
        contents,
        re.MULTILINE,
    )
    if match is None:
        raise Exception("public-op usage not found")
    public = int(match.group("usage"))
    return max(keypair, public, private)


def parse_filename(filename):
    match = re.match(
        r"elf/(?P<project>|mupq|mupq_pqclean)_?"
        r"crypto_(?P<type>kem|sign)_"
        r"(?P<scheme>[a-zA-Z0-9_\-]+)_"
        r"(?P<impl>[a-zA-Z0-9_\-]+)"
        r"_stack\.elf",
        filename,
    )
    if match is None:
        raise Exception("Can't parse filename!")
    project = match.group("project")
    scheme = match.group("scheme")
    impl = match.group("impl")
    return project, scheme, impl


def roundto(x, to):
    return x + (to - x % to)


def main():
    args = parse_arguments()
    schemes = []
    for f in args.inputs:
        contents = f.read()
        # Check for failed test
        try:
            filename, flashsize, ramsize = parse_flashsize(contents)
            project, scheme, impl = parse_filename(filename)
            match = re.search("HardFault", contents, re.MULTILINE)
            if match is None:
                stackusage = parse_stackusage(contents)
            else:
                stackusage = 4096 * 1024
        except Exception as e:
            if "usage not found" in str(e):
                stackusage = 4096 * 1024
            else:
                print(f"Error during parsing file {f.name}: {e}", file=sys.stderr)
                continue
        memoryusage = roundto(stackusage + ramsize + args.margin, args.round)
        print(
            f"Scheme: {scheme} Context: {project} Implementation: {impl} Flashsize: {flashsize} Memorysize: {memoryusage}",
            file=sys.stderr,
        )
        schemes.append(
            {
                "scheme": scheme,
                "implementation": impl,
                # 'project': (project + "/") if len(project) > 0 else "",
                "estmemory": memoryusage,
            }
        )
        # if project == "":
        #     del schemes[-1]['project']
    print("skip_list = [")
    for it in schemes:
        print(
            "    "
            + pprint.pformat(it, indent=4, compact=False, sort_dicts=False, width=110),
            end=",\n",
        )
    print("]")


if __name__ == "__main__":
    main()
